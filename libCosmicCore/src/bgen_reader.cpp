#include "bgen_reader.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <zlib.h>
#include <zstd.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <limits>
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

void BgenReader::stopPrefetcher() {
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
    stopPrefetcher();
    close();
}

void BgenReader::open(const std::string& fname) {
    filename = fname;
    ifs.open(filename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open BGEN file: " + filename);
    }
}

void BgenReader::close() {
    if (ifs.is_open()) {
        ifs.close();
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

template<typename T>
T BgenReader::read(std::ifstream& in) const {
    T val;
    in.read(reinterpret_cast<char*>(&val), sizeof(T));
    return val;
}

std::string BgenReader::readString(std::ifstream& in, uint16_t len_bytes) const {
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
    // Number of variants (4 bytes)
    num_variants = read<uint32_t>(ifs);
    // Number of samples (4 bytes)
    num_samples = read<uint32_t>(ifs);
    // Magic (4 bytes)
    char magic[4];
    ifs.read(magic, 4);
    if (std::strncmp(magic, "bgen", 4) != 0 && std::strncmp(magic, "\0\0\0\0", 4) != 0) {
        // Some older BGEN might differ, but v1.2 should be "bgen" or 0s.
        // Warn?
    }

    // Header length includes its own 4 bytes, followed by M, N, magic,
    // optional free data, and flags. We have already consumed the first
    // 16 bytes of that block (header_length, M, N, magic).
    if (header_length < 20) throw std::runtime_error("Invalid BGEN header length");
    uint32_t skip_len = header_length - 20;
    if (skip_len > 0) {
        ifs.seekg(skip_len, std::ios::cur);
    }
    
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
    if (has_sample_ids) {
        uint32_t block_len = read<uint32_t>(ifs);
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
    
    // Seek to first variant block
    // BGEN's leading offset is relative to byte 4 of the file.
    variant_data_offset = static_cast<size_t>(offset) + 4;
    ifs.seekg(static_cast<std::streamoff>(variant_data_offset), std::ios::beg);
    
    return num_samples;
}

void BgenReader::reset() const {
    if (ifs.is_open()) {
        ifs.clear();
        ifs.seekg(static_cast<std::streamoff>(variant_data_offset), std::ios::beg);
    }
}

bool BgenReader::decompress(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) const {
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
        
        // Store metadata
        BgenVariant v;
        v.chrom = chrom;
        v.rsid = rsid;
        v.pos = pos;
        v.alleles = {a1, a2}; // Assuming a1 is first allele (Ref?), a2 is second (Alt?)
        // BGEN usually stores Ref, Alt.
        variants.push_back(v);
        
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
        stopPrefetcher();
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
    if (prefetch_enabled) stopPrefetcher();
    reset();
}

void BgenReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_samples) throw std::runtime_error("Dimension mismatch in multiply_Zt_v");
    
    if (prefetch_enabled) {
        stopPrefetcher();
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
    if (prefetch_enabled) stopPrefetcher();
    reset();
}

long long BgenGenotypeProbs::observedCount() const {
    long long n = 0;
    for (int i = 0; i < size(); ++i) {
        if (missing[i] == 0.0f) ++n;
    }
    return n;
}

double BgenGenotypeProbs::meanDosage() const {
    double sum = 0.0;
    long long n = 0;
    for (int i = 0; i < size(); ++i) {
        if (missing[i] != 0.0f) continue;
        sum += static_cast<double>(dosage[i]);
        ++n;
    }
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    return sum / static_cast<double>(n);
}

double BgenGenotypeProbs::alleleFrequency() const {
    const long long n = observedCount();
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    return meanDosage() / 2.0;
}

bool BgenReader::readBlock(std::string& chrom, std::string& rsid, uint32_t& pos,
                           std::string& a1, std::string& a2, uint16_t& n_alleles,
                           bool skip_multiallelic, bool strict, uint32_t& d_out) const {
    if (!ifs || ifs.peek() == EOF) return false;

    // Variant identifying data (Layout 2).
    readString(ifs, 2);                    // id
    rsid = readString(ifs, 2);
    chrom = readString(ifs, 2);
    pos = read<uint32_t>(ifs);
    n_alleles = read<uint16_t>(ifs);

    if (n_alleles != 2) {
        if (!skip_multiallelic) {
            throw std::runtime_error(
                "BGEN genotype probability API does not support multiallelic records "
                "(n_alleles=" + std::to_string(n_alleles) + ", variant " + rsid + ")");
        }
        // Historical behavior: skip the record and recurse to the next one.
        for (int i = 0; i < n_alleles; ++i) {
            uint32_t len;
            ifs.read((char*)&len, 4);
            ifs.seekg(len, std::ios::cur);
        }
        uint32_t C = read<uint32_t>(ifs);
        ifs.seekg(C, std::ios::cur);
        return readBlock(chrom, rsid, pos, a1, a2, n_alleles, skip_multiallelic, strict, d_out);
    }

    uint32_t len1 = read<uint32_t>(ifs);
    std::string al1(len1, '\0');
    ifs.read(&al1[0], len1);

    uint32_t len2 = read<uint32_t>(ifs);
    std::string al2(len2, '\0');
    ifs.read(&al2[0], len2);

    // Scanner result paths report their second argument as the effect allele.
    // Match PLINK/PGEN A1 convention by returning BGEN allele 1 there.
    //
    // IMPORTANT: the returned dosage/probabilities below count copies of the
    // FIRST STORED BGEN allele, which is returned here as the SECOND out
    // parameter (a2). `G` in BgenGenotypeProbs is defined against a2, not a1.
    a1 = al2;
    a2 = al1;

    // Genotype data block: C = total length of the rest of the record,
    // D = uncompressed length, followed by C-4 bytes of payload.
    uint32_t C = read<uint32_t>(ifs);
    uint32_t D = read<uint32_t>(ifs);
    d_out = D;
    // The strict probability path validates the record framing. The legacy
    // dosage path deliberately keeps its historical (unvalidated) behavior.
    if (strict && C < 4) throw std::runtime_error("Malformed BGEN genotype block for variant " + rsid);
    if (strict && D == 0) throw std::runtime_error("Zero-length BGEN genotype block for variant " + rsid);

    if (compressed_buffer.size() < C - 4) compressed_buffer.resize(C - 4);
    ifs.read(compressed_buffer.data(), C - 4);
    if (strict && !ifs) throw std::runtime_error("Truncated BGEN genotype block for variant " + rsid);

    if (decompressed_buffer.size() < D) decompressed_buffer.resize(D);

    const uint32_t compression = flags & 0x3;
    if (compression == 0) {
        std::memcpy(decompressed_buffer.data(), compressed_buffer.data(), C - 4);
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
    return true;
}

bool BgenReader::readVariantProbs(std::string& chrom, std::string& rsid, uint32_t& pos,
                                  std::string& a1, std::string& a2,
                                  BgenGenotypeProbs& out,
                                  const BgenProbOptions& options) const {
    uint16_t n_alleles = 0;
    uint32_t D = 0;
    if (!readBlock(chrom, rsid, pos, a1, a2, n_alleles,
                   /*skip_multiallelic=*/!options.reject_multiallelic,
                   /*strict=*/true, D)) {
        return false;
    }
    parseLayout2Probs(decompressed_buffer.data(), D, num_samples,
                      n_alleles, out, options, rsid);
    return true;
}

void BgenReader::parseLayout2Probs(const char* data, size_t len, uint32_t n_samples,
                                   uint16_t n_alleles, BgenGenotypeProbs& out,
                                   const BgenProbOptions& options,
                                   const std::string& rsid) const {
    out.n_alleles = n_alleles;

    // Validate the entire fixed header, including sample ploidy bytes, before
    // dereferencing any field. The later probability-payload check cannot
    // protect a truncated header.
    const size_t minimum_header = 10u + static_cast<size_t>(n_samples);
    if (len < minimum_header) {
        throw std::runtime_error("Truncated BGEN genotype header for variant " + rsid);
    }
    const char* ptr = data;
    uint32_t N = *(uint32_t*)ptr; ptr += 4;
    if (N != n_samples) throw std::runtime_error("Sample count mismatch in genotype block");

    uint16_t K = *(uint16_t*)ptr; ptr += 2;
    if (K != n_alleles) throw std::runtime_error("Allele count mismatch");

    uint8_t Pmin = *(uint8_t*)ptr; ptr += 1;
    uint8_t Pmax = *(uint8_t*)ptr; ptr += 1;
    if (options.require_diploid && (Pmin != options.required_ploidy || Pmax != options.required_ploidy)) {
        throw std::runtime_error(
            "BGEN genotype probability API requires diploid records (Pmin=" +
            std::to_string(Pmin) + ", Pmax=" + std::to_string(Pmax) + ")");
    }

    const uint8_t* ploidy_ptr = reinterpret_cast<const uint8_t*>(ptr);
    ptr += N;

    uint8_t phased = *(uint8_t*)ptr; ptr += 1;
    if (phased) {
        throw std::runtime_error(
            "Phased BGEN records are not supported by the genotype probability API");
    }
    out.phased = false;

    uint8_t B = *(uint8_t*)ptr; ptr += 1;
    if (B != 8 && B != 16) {
        throw std::runtime_error("Unsupported BGEN bit depth: " + std::to_string(B));
    }
    out.prob_bits = B;
    out.ploidy = Pmax;

    // Record truncation guard: the payload must contain the fixed header, the
    // per-sample ploidy bytes, and the two stored probabilities per sample.
    // This is the check that actually protects the reads below, because the
    // legacy parser trusts the record's own sample count.
    const size_t header_bytes = 4 + 2 + 1 + 1;
    const size_t prob_bytes = (B == 8) ? 2u : 4u;
    const size_t need = header_bytes + static_cast<size_t>(N) + 1 + 1 +
                        static_cast<size_t>(N) * prob_bytes;
    if (len < need) {
        throw std::runtime_error(
            "Truncated BGEN genotype record for variant " + rsid +
            ": record declares " + std::to_string(len) + " bytes but requires at least " +
            std::to_string(need));
    }

    out.p0.setZero(N);
    out.p1.setZero(N);
    out.p2.setZero(N);
    out.missing.setZero(N);
    out.dosage.setConstant(N, std::numeric_limits<float>::quiet_NaN());
    out.missing_sample = false;

    const double denom = static_cast<double>((1 << B) - 1);

    if (B == 8) {
        const uint8_t* p8 = reinterpret_cast<const uint8_t*>(ptr);
        for (uint32_t i = 0; i < N; ++i) {
            const uint8_t ploidy_byte = ploidy_ptr[i];
            if (options.require_diploid && (ploidy_byte & 0x3F) != options.required_ploidy) {
                throw std::runtime_error(
                    "BGEN genotype probability API requires ploidy " +
                    std::to_string(options.required_ploidy) + " for every sample; found " +
                    std::to_string(static_cast<int>(ploidy_byte & 0x3F)) +
                    " at sample index " + std::to_string(i));
            }
            const bool missing = (ploidy_byte & 0x80) != 0;
            // Stored order is P(A0A0), P(A0A1) with P(A1A1) omitted.
            const double s0 = static_cast<double>(*p8++) / denom;
            const double s1 = static_cast<double>(*p8++) / denom;
            p8 += 0;
            if (missing) {
                out.missing[i] = 1.0f;
                out.missing_sample = true;
                continue;
            }
            const double p2 = s0;   // P(G=2) = P(counted/counted)
            const double p1 = s1;   // P(G=1) = P(het)
            const double p0 = 1.0 - s0 - s1;  // P(G=0) = P(non-counted/non-counted)
            if (options.validate_range) {
                if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(p2)) {
                    throw std::runtime_error("Non-finite BGEN probability at sample index " + std::to_string(i));
                }
                if (p0 < 0.0 || p1 < 0.0 || p2 < 0.0 || p0 > 1.0 || p1 > 1.0 || p2 > 1.0) {
                    throw std::runtime_error("BGEN probability outside [0,1] at sample index " +
                                             std::to_string(i));
                }
                if (std::fabs(p0 + p1 + p2 - 1.0) > options.sum_tolerance) {
                    throw std::runtime_error(
                        "BGEN probability sum deviates from 1 beyond tolerance at sample index " +
                        std::to_string(i) + " (sum=" + std::to_string(p0 + p1 + p2) + ")");
                }
            }
            out.p0[i] = static_cast<float>(p0);
            out.p1[i] = static_cast<float>(p1);
            out.p2[i] = static_cast<float>(p2);
            // Identical expression and operand order to the legacy dosage path.
            out.dosage[i] = static_cast<float>(2.0 * p2 + p1);
        }
    } else {
        const uint16_t* p16 = reinterpret_cast<const uint16_t*>(ptr);
        for (uint32_t i = 0; i < N; ++i) {
            const uint8_t ploidy_byte = ploidy_ptr[i];
            if (options.require_diploid && (ploidy_byte & 0x3F) != options.required_ploidy) {
                throw std::runtime_error(
                    "BGEN genotype probability API requires ploidy " +
                    std::to_string(options.required_ploidy) + " for every sample; found " +
                    std::to_string(static_cast<int>(ploidy_byte & 0x3F)) +
                    " at sample index " + std::to_string(i));
            }
            const bool missing = (ploidy_byte & 0x80) != 0;
            const double s0 = static_cast<double>(*p16++) / denom;
            const double s1 = static_cast<double>(*p16++) / denom;
            if (missing) {
                out.missing[i] = 1.0f;
                out.missing_sample = true;
                continue;
            }
            const double p2 = s0;
            const double p1 = s1;
            const double p0 = 1.0 - s0 - s1;
            if (options.validate_range) {
                if (!std::isfinite(p0) || !std::isfinite(p1) || !std::isfinite(p2)) {
                    throw std::runtime_error("Non-finite BGEN probability at sample index " + std::to_string(i));
                }
                if (p0 < 0.0 || p1 < 0.0 || p2 < 0.0 || p0 > 1.0 || p1 > 1.0 || p2 > 1.0) {
                    throw std::runtime_error("BGEN probability outside [0,1] at sample index " +
                                             std::to_string(i));
                }
                if (std::fabs(p0 + p1 + p2 - 1.0) > options.sum_tolerance) {
                    throw std::runtime_error(
                        "BGEN probability sum deviates from 1 beyond tolerance at sample index " +
                        std::to_string(i) + " (sum=" + std::to_string(p0 + p1 + p2) + ")");
                }
            }
            out.p0[i] = static_cast<float>(p0);
            out.p1[i] = static_cast<float>(p1);
            out.p2[i] = static_cast<float>(p2);
            out.dosage[i] = static_cast<float>(2.0 * p2 + p1);
        }
    }
}

bool BgenReader::readVariant(std::string& chrom, std::string& rsid, uint32_t& pos,
                             std::string& a1, std::string& a2,
                             Eigen::VectorXf& dosage) const {
    // Legacy dosage-only path. Behavior is intentionally unchanged, including
    // its historical lack of record-framing validation: strict=false.
    uint16_t n_alleles = 0;
    uint32_t D = 0;
    if (!readBlock(chrom, rsid, pos, a1, a2, n_alleles, /*skip_multiallelic=*/true,
                   /*strict=*/false, D)) {
        return false;
    }
    parseLayout2(decompressed_buffer.data(), decompressed_buffer.size(), num_samples,
                 n_alleles, dosage);
    return true;
}


void BgenReader::parseLayout2(const char* data, size_t len, uint32_t n_samples, uint16_t n_alleles, Eigen::VectorXf& dosage) const {
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
    
    // CosmicScanner uses dosage of allele 1: 2*P(AA) + P(AB).
    
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
                // Skip 2 bytes
                p8 += 2; 
                continue;
            }
            
            double p0 = (*p8++) / denom;
            double p1 = (*p8++) / denom;
            
            dosage[i] = static_cast<float>(2.0 * p0 + p1);
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
             
             double p0 = (*p16++) / denom;
             double p1 = (*p16++) / denom;
             
             dosage[i] = static_cast<float>(2.0 * p0 + p1);
        }
    } else {
        throw std::runtime_error("Unsupported BGEN bit depth: " + std::to_string(B));
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

        v.id = readString(ifs, 2);
        v.rsid = readString(ifs, 2);
        v.chrom = readString(ifs, 2);
        v.pos = read<uint32_t>(ifs);
        v.n_alleles = read<uint16_t>(ifs);

        if (v.chrom != current_chrom) {
            if (!current_chrom.empty()) {
                chrom_ranges_.emplace_back(chrom_start, i - chrom_start);
                chrom_labels_.push_back(current_chrom);
            }
            current_chrom = v.chrom;
            chrom_start = i;
        }

        v.alleles.clear();
        for (uint16_t a = 0; a < v.n_alleles; ++a) {
            uint32_t alen = read<uint32_t>(ifs);
            std::string allele(alen, '\0');
            ifs.read(&allele[0], alen);
            v.alleles.push_back(std::move(allele));
        }

        uint32_t C = read<uint32_t>(ifs);
        uint32_t D = read<uint32_t>(ifs);
        v.data_length = C + 4;
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

} // namespace cosmic

namespace cosmic {

void BgenReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    // Header
    geno_out << "VariantID";
    const auto& sample_ids = getSampleIds();
    for (const auto& sid : sample_ids) {
        geno_out << " " << sid;
    }
    geno_out << "\n";
    
    // Read all variants
    reset();
    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        geno_out << rsid;
        for (uint32_t i = 0; i < num_samples; ++i) {
            float val = dosage[i];
            if (std::isnan(val)) geno_out << " NA";
            else geno_out << " " << val;
        }
        geno_out << "\n";
    }
    
    // Variant info
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    reset();
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        var_out << rsid << " " << chrom << " " << pos << " " << a1 << " " << a2 << "\n";
    }
    
    // Sample info
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& sid : sample_ids) {
        size_t sep = sid.find(':');
        if (sep != std::string::npos) {
            sample_out << sid.substr(0, sep) << " " << sid.substr(sep + 1) << "\n";
        } else {
            sample_out << sid << " " << sid << "\n";
        }
    }
}

void BgenReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    // Collect all dosages
    std::vector<Eigen::VectorXf> all_dosages;
    std::vector<std::string> variant_ids;
    std::vector<std::string> chroms;
    std::vector<uint32_t> positions;
    std::vector<std::string> refs, alts;
    
    reset();
    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        all_dosages.push_back(dosage);
        variant_ids.push_back(rsid);
        chroms.push_back(chrom);
        positions.push_back(pos);
        refs.push_back(a1);
        alts.push_back(a2);
    }
    
    uint32_t n_variants = all_dosages.size();
    uint32_t n_samples = num_samples;
    
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2);
    
    // Write data
    for (uint32_t j = 0; j < n_variants; ++j) {
        for (uint32_t i = 0; i < n_samples; ++i) {
            float val = all_dosages[j][i];
            file.write(reinterpret_cast<const char*>(&val), 4);
        }
    }
    
    size_t data_end = file.tellp();
    
    std::ostringstream header;
    header << "{'descr': '<f4', 'fortran_order': False, 'shape': (" << n_variants << ", " << n_samples << "), }";
    std::string header_str = header.str();
    while ((header_str.size() + 10) % 64 != 0) header_str += " ";
    size_t header_len = header_str.size();
    
    file.seekp(header_len_pos);
    uint16_t header_len_le = static_cast<uint16_t>(header_len);
    file.write(reinterpret_cast<const char*>(&header_len_le), 2);
    file.write(header_str.c_str(), header_len);
    file.seekp(data_end);
    
    // Write var info
    std::ofstream var_out(var_info_file);
    var_out << "ID CHROM POS REF ALT\n";
    for (size_t j = 0; j < variant_ids.size(); ++j) {
        var_out << variant_ids[j] << " " << chroms[j] << " " << positions[j] << " " << refs[j] << " " << alts[j] << "\n";
    }
    
    // Sample info
    std::ofstream sample_out(sample_info_file);
    sample_out << "FID IID\n";
    for (const auto& sid : sample_ids) {
        size_t sep = sid.find(':');
        if (sep != std::string::npos) {
            sample_out << sid.substr(0, sep) << " " << sid.substr(sep + 1) << "\n";
        } else {
            sample_out << sid << " " << sid << "\n";
        }
    }
}

void BgenReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    // Export MAP
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    reset();
    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        map_out << chrom << " " << rsid << " 0 " << pos << " " << a1 << " " << a2 << "\n";
    }
    
    // Export PED
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    const auto& sample_ids = getSampleIds();
    reset();
    std::vector<std::vector<float>> all_dosages;
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        std::vector<float> var_dosage(num_samples);
        for (uint32_t i = 0; i < num_samples; ++i) var_dosage[i] = dosage[i];
        all_dosages.push_back(var_dosage);
    }
    
    for (uint32_t i = 0; i < num_samples; ++i) {
        ped_out << sample_ids[i] << " " << sample_ids[i] << " 0 0 1 -9";
        for (size_t j = 0; j < all_dosages.size(); ++j) {
            float val = all_dosages[j][i];
            if (std::isnan(val)) {
                ped_out << " 0 0";
            } else {
                // BGEN uses dosage, we need to convert to genotype calls
                // This is approximate
                if (val < 0.5) ped_out << " 0 0";
                else if (val < 1.5) ped_out << " 1 0";
                else ped_out << " 1 1";
            }
        }
        ped_out << "\n";
    }
}

void BgenReader::exportToVcf(const std::string& vcf_file) const {
    std::ofstream vcf_out(vcf_file);
    if (!vcf_out) throw std::runtime_error("Cannot open output VCF file: " + vcf_file);
    
    vcf_out << "##fileformat=VCFv4.2\n";
    vcf_out << "##source=CosmicScanner\n";
    vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    const auto& sample_ids = getSampleIds();
    for (const auto& sid : sample_ids) {
        vcf_out << "\t" << sid;
    }
    vcf_out << "\n";
    
    reset();
    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (readVariant(chrom, rsid, pos, a1, a2, dosage)) {
        vcf_out << chrom << "\t" << pos << "\t" << rsid << "\t" << a1 << "\t" << a2 << "\t.\tPASS\t.\tGT";
        
        for (uint32_t i = 0; i < num_samples; ++i) {
            float val = dosage[i];
            if (std::isnan(val)) {
                vcf_out << "\t./.";
            } else {
                if (val < 0.5) vcf_out << "\t0/0";
                else if (val < 1.5) vcf_out << "\t0/1";
                else vcf_out << "\t1/1";
            }
        }
        vcf_out << "\n";
    }
}

} // namespace cosmic
