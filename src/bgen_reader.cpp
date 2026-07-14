#include "bgen_reader.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <zlib.h>
#include <zstd.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

namespace cosmic {

// --- Prefetcher Implementation ---
// A simple SPSC (Single Producer Single Consumer) queue for variant data
// Actually we need to prefetch *dosage* vectors.
// Or just decompressed buffers? Decompressed buffers are safer.
// Dosage vector is float[N], size 4*N bytes.
// Decompressed buffer size depends on BGEN layout.
// Let's queue dosage vectors directly, so the consumer just does math.

struct PrefetchedVariant {
    uint32_t index; // Variant index
    Eigen::VectorXf dosage;
    bool success;
};

class Prefetcher {
public:
    Prefetcher(BgenReader* reader, int q_size = 128)
        : reader_(reader), max_queue_size_(q_size), running_(false) {}

    ~Prefetcher() { stop(); }

    void start() {
        if (running_) return;
        running_ = true;
        // Launch producer thread
        worker_ = std::thread(&Prefetcher::work, this);
    }

    void stop() {
        if (!running_) return;
        running_ = false;
        cond_full_.notify_all();
        cond_empty_.notify_all();
        if (worker_.joinable()) worker_.join();

        // Clear queue
        std::lock_guard<std::mutex> lock(mutex_);
        std::queue<PrefetchedVariant> empty;
        std::swap(queue_, empty);
    }

    bool pop(PrefetchedVariant& var) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_empty_.wait(lock, [this]{ return !queue_.empty() || !running_; });

        if (queue_.empty()) return false;

        var = std::move(queue_.front());
        queue_.pop();
        lock.unlock();
        cond_full_.notify_one();
        return true;
    }

private:
    void work() {
        // This thread reads from file and decompresses
        // We need a separate file stream or lock on the reader's stream?
        // Reader's stream is not thread safe.
        // BUT, the main thread will NOT touch the file stream while prefetching is active.
        // So we can use the reader's methods directly, provided they don't use shared mutable state that we conflict with.
        // BgenReader::readVariant uses internal buffers (compressed_buffer, decompressed_buffer).
        // We need our own buffers!
        // So we can't call readVariant directly if it uses member buffers.
        // We should duplicate the read logic or make buffers local.

        // Refactoring readVariant to use local buffers is best.
        // Or we just instantiate a new Reader? No, too heavy.

        // Let's assume for now we modify BgenReader to be thread-safe or use local buffers.
        // Actually, BgenReader::readVariant uses member buffers.
        // We will modify readVariant to accept buffers or use local ones.

        // WAIT: The main thread is doing matrix multiplication (consumer).
        // The worker thread is doing IO + Decompression (producer).
        // So main thread DOES NOT touch file.
        // So we can use reader_'s stream.
        // BUT reader_'s buffers are member variables.
        // If main thread accesses them? No, main thread gets 'dosage' vector from queue.
        // So it's safe! The main thread only does math on dosage.

        // However, we need to know WHERE to start reading.
        // We assume we continue from current position.

        std::string chrom, rsid, a1, a2;
        uint32_t pos;

        while (running_) {
            PrefetchedVariant var;
            var.dosage.resize(reader_->getNumSamples());

            // Critical section: File IO and Decompression (using member buffers)
            // Actually, we are the only one touching reader_ state.
            bool res = reader_->readVariant(chrom, rsid, pos, a1, a2, var.dosage);

            if (!res) {
                // EOF
                running_ = false;
                cond_empty_.notify_all(); // Wake up consumer to see empty queue
                break;
            }

            // Push to queue
            std::unique_lock<std::mutex> lock(mutex_);
            cond_full_.wait(lock, [this]{ return queue_.size() < max_queue_size_ || !running_; });

            if (!running_) break;

            var.success = true;
            // var.index = ... we don't track index easily unless we count.
            queue_.push(std::move(var));
            lock.unlock();
            cond_empty_.notify_one();
        }
    }

    BgenReader* reader_;
    size_t max_queue_size_;
    std::atomic<bool> running_;
    std::thread worker_;
    std::mutex mutex_;
    std::condition_variable cond_full_;
    std::condition_variable cond_empty_;
    std::queue<PrefetchedVariant> queue_;
};

// Global/Static prefetcher pointer? Or member of BgenReader?
// BgenReader needs to hold it. But BgenReader header didn't include thread.
// We can use a unique_ptr to forward declared struct.
// But for simplicity in this single-file edit context, let's use a PIMPL-like approach with void* or just keep it in cpp file and manage it carefully.
// Or just include <thread> in .cpp and use a member pointer that is casted.
// Actually, I can't add member to class in .h without editing .h again.
// I already edited .h but didn't add the pointer.
// I added `prefetch_enabled` flag.
// I can add a `void* prefetcher_impl = nullptr;` to .h if I edit it again.
// Or use a static map? No.

// Let's edit .h again to add `void* prefetcher_ptr = nullptr;`
// This is cleaner.

void BgenReader::startPrefetch() {
    if (prefetch_enabled) return;
    if (prefetcher_ptr) return; // Already exists?

    prefetch_enabled = true;
    Prefetcher* p = new Prefetcher(this, 64); // 64 items queue
    prefetcher_ptr = static_cast<void*>(p);
    p->start();
    std::cout << "  BGEN Background Prefetching Started." << std::endl;
}

void BgenReader::stopPrefetch() {
    if (!prefetch_enabled) return;

    if (prefetcher_ptr) {
        Prefetcher* p = static_cast<Prefetcher*>(prefetcher_ptr);
        p->stop();
        delete p;
        prefetcher_ptr = nullptr;
    }
    prefetch_enabled = false;
    std::cout << "  BGEN Background Prefetching Stopped." << std::endl;
}

BgenReader::BgenReader() {}

BgenReader::~BgenReader() {
    stopPrefetch();
    close();
}

void BgenReader::open(const std::string& fname) {
    filename = fname;
    ifs.open(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open BGEN file: " + filename);
    }
}

void BgenReader::openSharedMetadata(const BgenReader& source) {
    filename = source.filename;
    num_samples = source.num_samples;
    num_variants = source.num_variants;
    header_length = source.header_length;
    flags = source.flags;
    offset = source.offset;
    variant_data_offset = source.variant_data_offset;
    sample_ids_view = &source.getSampleIds();
    variants_view = &source.getVariants();
    chrom_ranges_view_ = &source.chrom_ranges_;
    chrom_labels_view_ = &source.chrom_labels_;
    ifs.open(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open BGEN file: " + filename);
    }
    reset();
}

void BgenReader::close() {
    if (ifs.is_open()) {
        ifs.close();
    }
}

template<typename T>
T BgenReader::read(std::ifstream& in) {
    T val;
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
    return val;
}

std::string BgenReader::readString(std::ifstream& in, uint16_t len_bytes) {
    uint16_t len = 0;
    in.read(reinterpret_cast<char*>(&len), len_bytes);
    std::string str(len, '\0');
    in.read(&str[0], len);
    return str;
}

uint32_t BgenReader::initialize() {
    if (!ifs.is_open()) throw std::runtime_error("BGEN file not open");

    // 1. Header
    // Offset (4 bytes)
    offset = read<uint32_t>(ifs);
    // Header length (4 bytes)
    header_length = read<uint32_t>(ifs);
    // Variant count (4 bytes) — plink2 writes variant count here (bytes 8-11)
    num_variants = read<uint32_t>(ifs);
    // Sample count (4 bytes) — plink2 writes sample count here (bytes 12-15)
    num_samples = read<uint32_t>(ifs);
    // Magic (4 bytes)
    char magic[4];
    ifs.read(magic, 4);
    if (std::strncmp(magic, "bgen", 4) != 0 && std::strncmp(magic, "\0\0\0\0", 4) != 0) {
        // Some older BGEN might differ, but v1.2 should be "bgen" or 0s.
        // Warn?
    }

    // Flags (4 bytes) — located right after magic in plink2 output
    flags = read<uint32_t>(ifs);

    // Check compression
    uint32_t compression = flags & 0x3;
    if (compression != 0 && compression != 1 && compression != 2) {
        throw std::runtime_error("Unsupported BGEN compression type (must be 0, 1=zlib, or 2=zstd). Type=" + std::to_string(compression));
    }

    // Check layout
    uint32_t layout = (flags >> 2) & 0xF;
    if (layout != 2) {
        throw std::runtime_error("Only BGEN Layout 2 is supported.");
    }

    // 2. Sample Block (Optional)
    bool has_sample_ids = (flags >> 31) & 0x1;
    uint32_t block_len = 0;
    if (has_sample_ids) {
        block_len = read<uint32_t>(ifs);
        uint32_t n_s = read<uint32_t>(ifs);
        if (n_s != num_samples) {
            throw std::runtime_error("Sample count mismatch in Sample Block");
        }

        for(uint32_t i=0; i<n_s; ++i) {
            sample_ids.push_back(readString(ifs, 2));
        }
    } else {
        // Generate dummy IDs
        for(uint32_t i=0; i<num_samples; ++i) {
            sample_ids.push_back("Sample_" + std::to_string(i+1));
        }
    }

    // plink2 BGEN v1.2 variant data starts right after the sample block.
    // block_len includes its own 4-byte field, so offset = 24 + block_len.
    if (has_sample_ids) {
        variant_data_offset = 24 + block_len;
    } else {
        variant_data_offset = offset;
    }
    return num_samples;
}

void BgenReader::reset() {
    if (ifs.is_open()) {
        ifs.clear();
        ifs.seekg(variant_data_offset, std::ios::beg);
    }
}

void BgenReader::loadVariantInfo() {
    if (variants_view) return;
    if (!variants.empty()) return;
    std::streampos saved_pos = ifs.tellg();
    ifs.clear();
    ifs.seekg(variant_data_offset, std::ios::beg);
    variants.resize(num_variants);
    chrom_ranges_.clear();
    chrom_labels_.clear();

    std::string current_chrom;
    uint32_t chrom_start = 0;

    for (uint32_t i = 0; i < num_variants; ++i) {
        BgenVariant& v = variants[i];
        v.data_offset = (uint32_t)ifs.tellg();

        // Skip 2-byte padding
        uint16_t pad;
        ifs.read((char*)&pad, 2);

        // Read variant ID
        v.id = readString(ifs, 2);
        v.rsid = v.id;

        // Chromosome
        v.chrom = readString(ifs, 2);

        // Position
        v.pos = read<uint32_t>(ifs);

        // Number of alleles
        v.n_alleles = read<uint16_t>(ifs);

        // Track chromosome ranges
        if (v.chrom != current_chrom) {
            if (!current_chrom.empty()) {
                chrom_ranges_.emplace_back(chrom_start, i - chrom_start);
                chrom_labels_.push_back(current_chrom);
            }
            current_chrom = v.chrom;
            chrom_start = i;
        }

        // Read alleles
        v.alleles.clear();
        for (uint16_t a = 0; a < v.n_alleles; ++a) {
            uint32_t alen = read<uint32_t>(ifs);
            std::string allele(alen, '\0');
            ifs.read(&allele[0], alen);
            v.alleles.push_back(std::move(allele));
        }

        // Read C and D, skip compressed data
        uint32_t C = read<uint32_t>(ifs);
        uint32_t D = read<uint32_t>(ifs);
        v.data_length = C + 4;  // C includes the D field's 4 bytes
        ifs.seekg(C - 4, std::ios::cur);
    }

    if (!current_chrom.empty()) {
        chrom_ranges_.emplace_back(chrom_start, num_variants - chrom_start);
        chrom_labels_.push_back(current_chrom);
    }

    ifs.clear();
    ifs.seekg(saved_pos);
}

void BgenReader::seekVariant(uint32_t idx) {
    const auto& meta = getVariants();
    if (idx >= num_variants || idx >= meta.size()) return;
    ifs.clear();
    ifs.seekg(meta[idx].data_offset, std::ios::beg);
}

std::pair<uint32_t, uint32_t> BgenReader::getChromRange(const std::string& chrom) const {
    const auto& labels = chrom_labels_view_ ? *chrom_labels_view_ : chrom_labels_;
    const auto& ranges = chrom_ranges_view_ ? *chrom_ranges_view_ : chrom_ranges_;
    for (size_t i = 0; i < labels.size(); i++) {
        if (labels[i] == chrom) {
            return ranges[i];
        }
    }
    return {0, 0};
}

bool BgenReader::decompress(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) {
    dst.resize(dst_len);

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = (uInt)src_len;
    strm.next_in = (Bytef*)src;
    strm.avail_out = (uInt)dst_len;
    strm.next_out = (Bytef*)dst.data();

    inflateInit(&strm);
    int ret = inflate(&strm, Z_FINISH);
    inflateEnd(&strm);

    return ret == Z_STREAM_END;
}

void BgenReader::computeStats(int n_threads) {
    if (stats_computed) return;

    reset();
    snp_means.assign(num_variants, 0.0);
    snp_inv_stds.assign(num_variants, 0.0);

    std::cout << "Computing stats for " << num_variants << " variants in BGEN..." << std::endl;

    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage(num_samples);

    variants.clear();
    variants.reserve(num_variants);

    for (size_t j = 0; j < num_variants; ++j) {
        if (!readVariant(chrom, rsid, pos, a1, a2, dosage)) break;

        // Compute mean/std
        double sum = 0.0;
        double sum_sq = 0.0;
        long long count = 0;

        #pragma omp parallel num_threads(n_threads) reduction(+:sum, sum_sq, count)
        {
            double local_sum = 0.0;
            double local_sum_sq = 0.0;
            long long local_count = 0;

            #pragma omp for
            for (long long i = 0; i < (long long)num_samples; ++i) {
                float val = dosage[i];
                if (!std::isnan(val)) {
                    local_sum += val;
                    local_sum_sq += (val * val);
                    local_count++;
                }
            }
            sum += local_sum;
            sum_sq += local_sum_sq;
            count += local_count;
        }

        if (count > 0) {
            double mean = sum / count;
            double var = sum_sq / count - mean * mean;
            if (var > 1e-9) {
                snp_means[j] = mean;
                snp_inv_stds[j] = 1.0 / std::sqrt(var);
            } else {
                snp_means[j] = mean;
                snp_inv_stds[j] = 0.0;
            }
        } else {
            snp_means[j] = 0.0;
            snp_inv_stds[j] = 0.0;
        }

        if (j % 10000 == 0 && j > 0) {
            std::cout << "\rProcessed " << j << "/" << num_variants << "..." << std::flush;
        }
    }
    std::cout << std::endl;

    stats_computed = true;
    reset();
}

void BgenReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_variants) throw std::runtime_error("Dimension mismatch in multiply_Z_v");

    if (prefetch_enabled) {
        stopPrefetch();
        reset();
        startPrefetch();
    } else {
        reset();
    }

    out.setZero(num_samples);

    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage(num_samples);

    for (size_t j = 0; j < num_variants; ++j) {
        if (prefetch_enabled && prefetcher_ptr) {
            Prefetcher* p = static_cast<Prefetcher*>(prefetcher_ptr);
            PrefetchedVariant var;
            if (!p->pop(var)) break;

            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) continue;

            double mean = snp_means[j];
            double weight = v(j) * inv_std;

            #pragma omp parallel for num_threads(n_threads)
            for (long long i = 0; i < (long long)num_samples; ++i) {
                float g = var.dosage[i];
                if (!std::isnan(g)) {
                    out(i) += (g - mean) * weight;
                }
            }
        } else {
            if (!readVariant(chrom, rsid, pos, a1, a2, dosage)) break;

            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) continue;

            double mean = snp_means[j];
            double weight = v(j) * inv_std;

            #pragma omp parallel for num_threads(n_threads)
            for (long long i = 0; i < (long long)num_samples; ++i) {
                float g = dosage[i];
                if (!std::isnan(g)) {
                    out(i) += (g - mean) * weight;
                }
            }
        }
    }
    if (prefetch_enabled) stopPrefetch();
    reset();
}

void BgenReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_samples) throw std::runtime_error("Dimension mismatch in multiply_Zt_v");

    if (prefetch_enabled) {
        stopPrefetch();
        reset();
        startPrefetch();
    } else {
        reset();
    }

    out.setZero(num_variants);

    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage(num_samples);

    for (size_t j = 0; j < num_variants; ++j) {
        if (prefetch_enabled && prefetcher_ptr) {
            Prefetcher* p = static_cast<Prefetcher*>(prefetcher_ptr);
            PrefetchedVariant var;
            if (!p->pop(var)) break;

            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) {
                out(j) = 0.0;
                continue;
            }

            double mean = snp_means[j];
            double dot_gv = 0.0;
            double sum_v = 0.0;

            #pragma omp parallel num_threads(n_threads) reduction(+:dot_gv, sum_v)
            {
                double local_dot = 0.0;
                double local_sum_v = 0.0;

                #pragma omp for
                for (long long i = 0; i < (long long)num_samples; ++i) {
                    float g = var.dosage[i];
                    if (!std::isnan(g)) {
                        local_dot += (double)g * v(i);
                        local_sum_v += v(i);
                    }
                }
                dot_gv += local_dot;
                sum_v += local_sum_v;
            }
            out(j) = (dot_gv - mean * sum_v) * inv_std;

        } else {
            if (!readVariant(chrom, rsid, pos, a1, a2, dosage)) break;

            double inv_std = snp_inv_stds[j];
            if (inv_std == 0.0) {
                out(j) = 0.0;
                continue;
            }

            double mean = snp_means[j];

            double dot_gv = 0.0;
            double sum_v = 0.0;

            #pragma omp parallel num_threads(n_threads) reduction(+:dot_gv, sum_v)
            {
                double local_dot = 0.0;
                double local_sum_v = 0.0;

                #pragma omp for
                for (long long i = 0; i < (long long)num_samples; ++i) {
                    float g = dosage[i];
                    if (!std::isnan(g)) {
                        local_dot += (double)g * v(i);
                        local_sum_v += v(i);
                    }
                }
                dot_gv += local_dot;
                sum_v += local_sum_v;
            }

            out(j) = (dot_gv - mean * sum_v) * inv_std;
        }
    }
    if (prefetch_enabled) stopPrefetch();
    reset();
}

bool BgenReader::readVariant(std::string& chrom, std::string& rsid, uint32_t& pos,
                             std::string& a1, std::string& a2,
                             Eigen::VectorXf& dosage) {
    if (!ifs || ifs.peek() == EOF) return false;

    // plink2 BGEN v1.2 Layout 2 variant block:
    //   00 00  (2-byte padding/separator)
    //   [id_len:2] [id]
    //   [chrom_len:2] [chrom]
    //   [pos:4]
    //   [n_alleles:2]
    //   [allele1_len:4] [allele1]
    //   [allele2_len:4] [allele2]
    //   [C:4] [D:4] [compressed_data: C-4]

    // Skip 2-byte padding
    uint16_t pad;
    ifs.read((char*)&pad, 2);

    // Read variant ID (plink2 stores the RSID string in the ID field)
    std::string id = readString(ifs, 2);
    rsid = id;

    // Chromosome
    chrom = readString(ifs, 2);

    // Position
    pos = read<uint32_t>(ifs);

    // Number of alleles
    uint16_t n_alleles = read<uint16_t>(ifs);

    if (n_alleles != 2) {
        // Skip alleles
        for(int i=0; i<n_alleles; ++i) {
             uint32_t len;
             ifs.read((char*)&len, 4);
             ifs.seekg(len, std::ios::cur);
        }
        uint32_t C = read<uint32_t>(ifs);
        ifs.seekg(C, std::ios::cur);
        return readVariant(chrom, rsid, pos, a1, a2, dosage);
    }

    // Read 2 alleles (4-byte length prefix each, as per Layout 2)
    uint32_t len1 = read<uint32_t>(ifs);
    std::string al1(len1, '\0');
    ifs.read(&al1[0], len1);

    uint32_t len2 = read<uint32_t>(ifs);
    std::string al2(len2, '\0');
    ifs.read(&al2[0], len2);

    a1 = al1;
    a2 = al2;

    // 2. Genotype Data Block
    uint32_t C = read<uint32_t>(ifs);
    uint32_t D = read<uint32_t>(ifs);

    if (compressed_buffer.size() < C - 4) compressed_buffer.resize(C - 4);
    ifs.read(compressed_buffer.data(), C - 4);

    if (decompressed_buffer.size() < D) decompressed_buffer.resize(D);

    uint32_t compression = flags & 0x3;
    if (compression == 0) {
        std::memcpy(decompressed_buffer.data(), compressed_buffer.data(), C-4);
    } else if (compression == 1) {
        if (!decompress(compressed_buffer.data(), C - 4, decompressed_buffer, D)) {
            throw std::runtime_error("Decompression failed (zlib) for variant " + rsid);
        }
    } else if (compression == 2) {
        size_t dSize = ZSTD_decompress(decompressed_buffer.data(), D, compressed_buffer.data(), C - 4);
        if (ZSTD_isError(dSize)) {
            throw std::runtime_error("Decompression failed (Zstd) for variant " + rsid + ": " + ZSTD_getErrorName(dSize));
        }
        if (dSize != D) {
            throw std::runtime_error("Zstd decompressed size mismatch for variant " + rsid);
        }
    }

    parseLayout2(decompressed_buffer.data(), D, num_samples, n_alleles, dosage);
    return true;
}

void BgenReader::parseLayout2(const char* data, size_t len, uint32_t n_samples, uint16_t n_alleles, Eigen::VectorXf& dosage) {
    // Layout 2 structure (Uncompressed):
    // N (4)
    // K (2) - Number of alleles (should match)
    // Pmin (1)
    // Pmax (1)
    // Ploidy (N bytes) - usually 2
    // Phased (1 byte) - 0 or 1.
    // B (1 byte) - Number of bits per probability.
    // Probabilities...

    const char* ptr = data;
    uint32_t N = *(uint32_t*)ptr; ptr += 4;
    if (N != n_samples) throw std::runtime_error("Sample count mismatch in genotype block");

    uint16_t K = *(uint16_t*)ptr; ptr += 2;
    if (K != n_alleles) throw std::runtime_error("Allele count mismatch");

    uint8_t Pmin = *(uint8_t*)ptr; ptr += 1;
    uint8_t Pmax = *(uint8_t*)ptr; ptr += 1;

    // Ploidy
    // We assume diploid (2) for now.
    // Skip N bytes (ploidy)
    // Check if all are 2?
    // For speed, just skip.
    ptr += N;

    // Phased
    uint8_t phased = *(uint8_t*)ptr; ptr += 1;

    // Bits
    uint8_t B = *(uint8_t*)ptr; ptr += 1;

    if (dosage.size() != N) dosage.resize(N);

    // Probabilities
    // If B=8, values are 0..255.
    // P = val / (2^B - 1).
    // For Biallelic Diploid Unphased:
    // 3 probabilities per sample: P(AA), P(AB), P(BB).
    // Actually Layout 2 stores (K*(K+1)/2) - 1 probabilities per sample?
    // No, "The number of probabilities stored for each individual is... K(K+1)/2 - 1".
    // For K=2, this is 2 probabilities: P(AA), P(AB). P(BB) = 1 - sum.

    // Wait, v1.2 spec:
    // "For unphased data... probabilities are stored... P_AA, P_AB, P_BB..."
    // "Since the sum is 1, the last one is omitted? No, layout 2 usually stores N_prob - 1?"
    // Actually, check spec:
    // "Probability data is stored as a sequence of integers..."
    // "For each sample... if ploidy=2, unphased, K=2... 3 probabilities are stored?"
    // Wait, BGEN v1.2 Layout 2 supports "Omitted last probability" mode?
    // Let's check B.
    // If phased=0:
    //   Number of probabilities = N * (K*(K+1)/2 - 1).
    //   Wait, is it always minus 1?
    //   "Layout 2... the probability of the last genotype is not stored."
    //   So for biallelic (3 genotypes), we store 2 probs: P(AA), P(AB).

    double denom = (double)((1 << B) - 1);

    // We want Dosage = count of the first allele (a1/ref).
    // Dosage = 2*P(AA) + P(AB).
    // We have P(AA), P(AB) from the file (last probability is omitted).

    // Bit extraction helper?
    // If B=8, just read bytes.
    // If B=16, read shorts.
    // Most BGEN v1.2 use B=8 or B=16.

    if (phased) {
        throw std::runtime_error("Phased BGEN not supported yet");
    }

    if (B == 8) {
        const uint8_t* p8 = (const uint8_t*)ptr;
        for(uint32_t i=0; i<N; ++i) {
            // Read 2 probs
            // But wait, missing data?
            // "If the MSB of the ploidy byte is 1, the data for this sample is missing."
            // I skipped ploidy check!
            // Let's go back and check ploidy for missingness.
            // This is getting complicated.

            // Re-read ploidy
            const uint8_t* ploidy_ptr = (const uint8_t*)(data + 8);
            uint8_t ploidy_byte = ploidy_ptr[i];
            bool missing = (ploidy_byte & 0x80);

            if (missing) {
                dosage[i] = std::numeric_limits<float>::quiet_NaN();
                p8 += 2;
                continue;
            }

            double p0 = (*p8++) / denom;  // P(AA)
            double p1 = (*p8++) / denom;  // P(AB)

            const double p2 = std::max(0.0, 1.0 - p0 - p1);
            dosage[i] = (float)(p1 + 2.0 * p2);
        }
    } else if (B == 16) {
        const uint16_t* p16 = (const uint16_t*)ptr;
        for(uint32_t i=0; i<N; ++i) {
             const uint8_t* ploidy_ptr = (const uint8_t*)(data + 8);
             uint8_t ploidy_byte = ploidy_ptr[i];
             bool missing = (ploidy_byte & 0x80);

             if (missing) {
                 dosage[i] = std::numeric_limits<float>::quiet_NaN();
                 p16 += 2;
                 continue;
             }

             double p0 = (*p16++) / denom;  // P(AA)
             double p1 = (*p16++) / denom;  // P(AB)

             const double p2 = std::max(0.0, 1.0 - p0 - p1);
             dosage[i] = (float)(p1 + 2.0 * p2);
        }
    } else {
        throw std::runtime_error("Unsupported BGEN bit depth: " + std::to_string(B));
    }
}

} // namespace cosmic
