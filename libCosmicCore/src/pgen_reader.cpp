#include "pgen_reader.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <zstd.h>
#include <cstring>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <omp.h>

namespace cosmic {

// PGEN Constants
static const uint16_t kPgenMagic = 0x1B6C; // 0x6C, 0x1B
static const uint8_t kPgenTypeUncompressed = 0;

PgenReader::PgenReader() {}

PgenReader::~PgenReader() {
    close();
}

void PgenReader::open(const std::string& pgen_file, const std::string& pvar_file, const std::string& psam_file) {
    auto has_suffix = [](const std::string& value, const std::string& suffix) {
        return value.size() >= suffix.size() &&
               value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
    };
    const std::string prefix = has_suffix(pgen_file, ".pgen")
        ? pgen_file.substr(0, pgen_file.size() - 5)
        : pgen_file;
    pgen_path = has_suffix(pgen_file, ".pgen") ? pgen_file : prefix + ".pgen";
    pvar_path = pvar_file.empty() ? prefix + ".pvar" : pvar_file;
    psam_path = psam_file.empty() ? prefix + ".psam" : psam_file;
    
    // Open MMap
    mmap_file.open(pgen_path);
    if (!mmap_file.is_open()) {
        throw std::runtime_error("Failed to open PGEN file: " + pgen_path);
    }
    
    readPgenHeader();
    readPvar();
    readPsam();
    
    // Validate
    if (variants.size() != n_variants_pgen) {
        std::cerr << "Warning: Number of variants in PGEN header (" << n_variants_pgen 
                  << ") does not match PVAR (" << variants.size() << "). Using PVAR count.\n";
        num_variants = (uint32_t)variants.size();
    } else {
        num_variants = n_variants_pgen;
    }
    
    if (sample_ids.size() != n_samples_pgen) {
        std::cerr << "Warning: Number of samples in PGEN header (" << n_samples_pgen 
                  << ") does not match PSAM (" << sample_ids.size() << "). Using PSAM count.\n";
        num_samples = (uint32_t)sample_ids.size();
    } else {
        num_samples = n_samples_pgen;
    }
    
    reset();
    stats_computed = false;
}

void PgenReader::close() {
    mmap_file.close();
    variants.clear();
    sample_ids.clear();
    snp_means.clear();
    snp_inv_stds.clear();
    stats_computed = false;
    resetMetadataViews();
}

void PgenReader::resetMetadataViews() {
    sample_ids_view_ = &sample_ids;
    variants_view_ = &variants;
    vrtypes_view_ = &vrtypes_;
    vlens_view_ = &vlens_;
    variant_offsets_view_ = &variant_offsets_;
    chrom_ranges_view_ = &chrom_ranges_;
    chrom_labels_view_ = &chrom_labels_;
}

void PgenReader::parseDirectory(const uint8_t* data, size_t size) {
    uint32_t fmt = (format_byte & 0x0F);
    bool vrtype_8bit = (fmt >= 4);
    uint32_t vrec_len_bytes = (fmt & 3) + 1;
    uint32_t allele_count_bytes = (format_byte >> 4) & 3;
    uint32_t provref_mode = (format_byte >> 6) & 3;

    uint32_t num_blocks = (n_variants_pgen + 65535) / 65536;
    size_t header_body_start = 12 + (size_t)num_blocks * 8;
    if (header_body_start >= size) {
        throw std::runtime_error("PGEN header body out of bounds");
    }

    std::vector<size_t> block_offsets(num_blocks);
    for (uint32_t b = 0; b < num_blocks; b++) {
        block_offsets[b] = (size_t)(*reinterpret_cast<const uint64_t*>(data + 12 + (size_t)b * 8));
        if (block_offsets[b] >= size) throw std::runtime_error("PGEN variant block offset out of bounds");
    }

    const uint8_t* ptr = data + header_body_start;
    const uint8_t* end = data + size;
    vrtypes_.clear();
    vrtypes_.reserve(n_variants_pgen);
    vlens_.assign(n_variants_pgen, 0);
    variant_offsets_.assign(n_variants_pgen, 0);

    for (uint32_t block = 0; block < num_blocks; block++) {
        uint32_t block_start = block * 65536u;
        uint32_t block_count = std::min<uint32_t>(65536u, n_variants_pgen - block_start);
        std::vector<uint8_t> block_vrtypes(block_count);

        if (vrtype_8bit) {
            if ((size_t)(end - ptr) < block_count) throw std::runtime_error("PGEN vrtype truncated");
            for (uint32_t i = 0; i < block_count; i++) block_vrtypes[i] = ptr[i];
            ptr += block_count;
        } else {
            uint32_t type_bytes = (block_count + 1) / 2;
            if ((size_t)(end - ptr) < type_bytes) throw std::runtime_error("PGEN vrtype truncated");
            for (uint32_t i = 0; i < block_count; i++) {
                uint8_t packed = ptr[i / 2];
                block_vrtypes[i] = (i % 2 == 0) ? (packed & 0x0F) : (packed >> 4);
            }
            ptr += type_bytes;
        }

        if ((size_t)(end - ptr) < (size_t)block_count * vrec_len_bytes) {
            throw std::runtime_error("PGEN vrec_len truncated");
        }
        size_t block_rel = block_offsets[block] - block_offsets[0];
        size_t cum_in_block = 0;
        for (uint32_t i = 0; i < block_count; i++) {
            uint32_t vlen = 0;
            for (uint32_t j = 0; j < vrec_len_bytes; j++) vlen |= (uint32_t)ptr[j] << (j * 8);
            ptr += vrec_len_bytes;
            uint32_t vidx = block_start + i;
            vrtypes_.push_back(block_vrtypes[i]);
            vlens_[vidx] = vlen;
            variant_offsets_[vidx] = block_rel + cum_in_block;
            cum_in_block += vlen;
        }

        if (allele_count_bytes != 0) {
            size_t bytes = (size_t)block_count * allele_count_bytes;
            if ((size_t)(end - ptr) < bytes) throw std::runtime_error("PGEN allele-count header truncated");
            ptr += bytes;
        }
        if (provref_mode == 3) {
            size_t bytes = (block_count + 7) / 8;
            if ((size_t)(end - ptr) < bytes) throw std::runtime_error("PGEN provisional-ref header truncated");
            ptr += bytes;
        }
    }
}

void PgenReader::reset() const {
    current_offset = variant_block_offset;
    current_variant_idx = 0;
    ldbase_genotypes_.clear();
}

// Helper to read PGEN length (1-5 bytes)
static uint32_t readPgenLength(const uint8_t*& ptr) {
    uint8_t b1 = *ptr++;
    if (b1 < 253) return b1;
    if (b1 == 253) {
        uint32_t len = *reinterpret_cast<const uint16_t*>(ptr);
        ptr += 2;
        return len;
    }
    if (b1 == 254) {
        uint32_t len = ptr[0] | (ptr[1] << 8) | (ptr[2] << 16);
        ptr += 3;
        return len;
    }
    // b1 == 255
    uint32_t len = *reinterpret_cast<const uint32_t*>(ptr);
    ptr += 4;
    return len;
}

// Helper to get raw variant data (decompressing if needed)
// Returns pointer to data and its length. Updates current_offset.
// type is updated to inner type if decompressed.
// decomp_buf is used for decompression storage.
static const uint8_t* getVariantData(const uint8_t* mmap_data, size_t file_size, 
                                     size_t& offset, uint8_t& type, size_t& len,
                                     std::vector<char>& decomp_buf, uint32_t num_samples) {
    if (offset >= file_size) return nullptr;
    
    type = mmap_data[offset++];
    
    if (type == kPgenTypeUncompressed) {
        // Uncompressed (Type 0)
        len = (num_samples + 3) / 4;
        if (offset + len > file_size) return nullptr;
        const uint8_t* ptr = mmap_data + offset;
        offset += len;
        return ptr;
    } else {
        // Compressed or Sparse (Type > 0)
        // Read length of following data
        const uint8_t* ptr_len = mmap_data + offset;
        uint32_t data_len = readPgenLength(ptr_len);
        offset = (ptr_len - mmap_data); // Update offset after length read
        
        if (offset + data_len > file_size) return nullptr;
        
        const uint8_t* data_ptr = mmap_data + offset;
        offset += data_len;
        
        // Check for Zstd Compression (Type >= 64)
        if (type >= 64) {
            uint8_t inner_type = type & 63;
            
            // Decompress
            unsigned long long const dSize = ZSTD_getFrameContentSize(data_ptr, data_len);
            if (dSize == ZSTD_CONTENTSIZE_ERROR || dSize == ZSTD_CONTENTSIZE_UNKNOWN) {
                // Fallback: Use expected size for Type 0 or grow?
                // For Type 0, expected size is (N+3)/4.
                // Let's assume we can resize to expected uncompressed size if Type 0.
                if (inner_type == kPgenTypeUncompressed) {
                    size_t expected = (num_samples + 3) / 4;
                    if (decomp_buf.size() < expected) decomp_buf.resize(expected);
                } else {
                    // Unknown size and unknown type -> unsafe.
                    // Try a reasonable buffer size?
                    if (decomp_buf.size() < 65536) decomp_buf.resize(65536);
                }
            } else {
                if (decomp_buf.size() < dSize) decomp_buf.resize(dSize);
            }
            
            size_t const decompressed_size = ZSTD_decompress(decomp_buf.data(), decomp_buf.size(), data_ptr, data_len);
            if (ZSTD_isError(decompressed_size)) {
                // std::cerr << "ZSTD Error: " << ZSTD_getErrorName(decompressed_size) << "\n";
                return nullptr;
            }
            
            type = inner_type;
            len = decompressed_size;
            return (const uint8_t*)decomp_buf.data();
        } else {
            // Not Zstd compressed (e.g. Type 1-63 stored directly)
            // Return raw data
            len = data_len;
            return data_ptr;
        }
    }
}

void PgenReader::computeStats(int n_threads) {
    if (stats_computed) return;

    (void)n_threads;
    reset();
    snp_means.assign(num_variants, 0.0);
    snp_inv_stds.assign(num_variants, 0.0);

    std::cout << "Computing stats for " << num_variants << " variants in PGEN..." << std::endl;
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    for (size_t j = 0; j < num_variants; ++j) {
        if (!readVariant(chrom, id, pos, ref, alt, dosage)) break;
        double sum = 0.0;
        double sum_sq = 0.0;
        size_t count = 0;
        for (Eigen::Index i = 0; i < dosage.size(); ++i) {
            if (!std::isfinite(dosage[i])) continue;
            sum += dosage[i];
            sum_sq += dosage[i] * dosage[i];
            ++count;
        }
        if (!count) continue;
        const double mean = sum / static_cast<double>(count);
        const double var = sum_sq / static_cast<double>(count) - mean * mean;
        snp_means[j] = mean;
        snp_inv_stds[j] = var > 1e-9 ? 1.0 / std::sqrt(var) : 0.0;
    }

    stats_computed = true;
    reset();
}

void PgenReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_variants) throw std::runtime_error("Dimension mismatch in multiply_Z_v");
    
    reset();
    out.setZero(num_samples);
    
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    for (size_t j = 0; j < num_variants; ++j) {
        if (!readVariant(chrom, id, pos, ref, alt, dosage)) break;
        const double inv_std = snp_inv_stds[j];
        if (inv_std == 0.0 || std::isnan(snp_means[j])) continue;
        const double weight = v(j) * inv_std;
        for (Eigen::Index i = 0; i < dosage.size(); ++i) {
            if (std::isfinite(dosage[i])) out[i] += (dosage[i] - snp_means[j]) * weight;
        }
    }
    reset();
}

void PgenReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_samples) throw std::runtime_error("Dimension mismatch in multiply_Zt_v");
    
    reset();
    out.setZero(num_variants);
    
    std::string chrom, id, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    for (size_t j = 0; j < num_variants; ++j) {
        if (!readVariant(chrom, id, pos, ref, alt, dosage)) break;
        const double inv_std = snp_inv_stds[j];
        if (inv_std == 0.0 || std::isnan(snp_means[j])) continue;
        double dot = 0.0;
        for (Eigen::Index i = 0; i < dosage.size(); ++i) {
            if (std::isfinite(dosage[i])) dot += (dosage[i] - snp_means[j]) * v[i];
        }
        out[j] = dot * inv_std;
    }
    reset();
}

void PgenReader::readPsam() {
    std::ifstream in(psam_path);
    if (!in) return;
    
    sample_ids.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string fid, iid;
        ss >> fid >> iid; // Assuming FID IID columns
        sample_ids.push_back(fid + ":" + iid);
    }
}

void PgenReader::readPvar() {
    std::ifstream in(pvar_path);
    if (!in) return;
    
    variants.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        PgenVariant v;
        ss >> v.chrom >> v.pos >> v.id >> v.ref >> v.alt;
        variants.push_back(v);
    }
}

void PgenReader::readPgenHeader() {
    const uint8_t* data = (const uint8_t*)mmap_file.get_data();
    size_t size = mmap_file.get_size();
    
    if (size < 12) throw std::runtime_error("PGEN file too small");
    
    size_t offset = 0;
    
    // Magic
    if (data[0] != 0x6C || data[1] != 0x1B) {
        throw std::runtime_error("Invalid PGEN magic number");
    }
    offset += 2;
    
    // Storage Mode
    storage_mode = data[offset++];
    
    // Num Variants
    n_variants_pgen = *reinterpret_cast<const uint32_t*>(data + offset);
    offset += 4;
    
    // Num Samples
    n_samples_pgen = *reinterpret_cast<const uint32_t*>(data + offset);
    offset += 4;
    
    if (storage_mode == 0x01) {
        // PLINK 1 BED mode
        variant_block_offset = 0;
        current_offset = 0;
        return;
    }

    // Header control/format byte immediately follows the sample count.
    format_byte = data[offset++];
    provisional_ref = ((format_byte >> 6) & 3) == 2;

    if (storage_mode >= 0x10) {
        // Mode 0x10/0x11: variable-width with offset table + directory
        uint32_t num_blocks = (n_variants_pgen + 65535) / 65536;
        if (offset + num_blocks * 8 > size) {
            throw std::runtime_error("PGEN offset table out of bounds");
        }
        variant_block_offset = *reinterpret_cast<const uint64_t*>(data + offset);
        current_offset = variant_block_offset;

        // Parse directory
        parseDirectory(data, size);
    } else {
        // Fixed-width mode (0x02): format_byte at offset 11, data starts at offset 12
        variant_block_offset = offset;
        current_offset = variant_block_offset;
        // All variants use the same encoding (from format_byte low nibble) and fixed length
        uint8_t fixed_vrtype = format_byte & 0x0F;
        uint32_t bit_width = (fixed_vrtype <= 1) ? (fixed_vrtype == 0 ? 2 : 1) : 2;
        uint32_t fixed_vlen = (n_samples_pgen * bit_width + 7) / 8;
        vrtypes_.assign(n_variants_pgen, fixed_vrtype);
        vlens_.assign(n_variants_pgen, fixed_vlen);
    }
}

void PgenReader::openSharedMetadata(const PgenReader& source) {
    close();
    pgen_path = source.pgen_path;
    pvar_path = source.pvar_path;
    psam_path = source.psam_path;
    mmap_file.open(pgen_path);
    if (!mmap_file.is_open()) {
        throw std::runtime_error("Failed to open PGEN file: " + pgen_path);
    }

    num_samples = source.num_samples;
    num_variants = source.num_variants;
    storage_mode = source.storage_mode;
    variant_block_offset = source.variant_block_offset;
    n_samples_pgen = source.n_samples_pgen;
    n_variants_pgen = source.n_variants_pgen;
    format_byte = source.format_byte;
    provisional_ref = source.provisional_ref;

    sample_ids_view_ = source.sample_ids_view_;
    variants_view_ = source.variants_view_;
    vrtypes_view_ = source.vrtypes_view_;
    vlens_view_ = source.vlens_view_;
    variant_offsets_view_ = source.variant_offsets_view_;
    chrom_ranges_view_ = source.chrom_ranges_view_;
    chrom_labels_view_ = source.chrom_labels_view_;

    reset();
    stats_computed = false;
    gt_buffer_.resize(num_samples);
}

void PgenReader::seekVariant(uint32_t idx) {
    if (idx >= num_variants) {
        current_variant_idx = num_variants;
        const auto& offs = *variant_offsets_view_;
        const auto& lens = *vlens_view_;
        current_offset = variant_block_offset + (offs.empty() ? 0 : offs.back() + lens.back());
        return;
    }
    if (storage_mode >= 0x10 && idx > 0) {
        reset();
        std::string chrom, id, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        while (current_variant_idx < idx && readVariant(chrom, id, pos, ref, alt, dosage)) {}
        return;
    }
    const auto& offs = *variant_offsets_view_;
    const auto& lens = *vlens_view_;
    if (!offs.empty() && idx < offs.size()) {
        current_offset = variant_block_offset + offs[idx];
    } else {
        current_offset = variant_block_offset;
        for (uint32_t i = 0; i < idx; i++) {
            if (i < lens.size()) current_offset += lens[i];
        }
    }
    current_variant_idx = idx;
}

void PgenReader::releaseVariantRange(uint32_t start_idx, uint32_t count) const {
    if (count == 0 || start_idx >= num_variants) return;
    uint32_t end_idx = std::min<uint32_t>(num_variants, start_idx + count);
    const auto& offs = *variant_offsets_view_;
    const auto& lens = *vlens_view_;
    if (offs.empty() || start_idx >= offs.size() || end_idx == 0) return;
    size_t start = variant_block_offset + offs[start_idx];
    size_t end = 0;
    uint32_t last = end_idx - 1;
    if (last < offs.size() && last < lens.size()) end = variant_block_offset + offs[last] + lens[last];
    if (end > start) mmap_file.adviseDontNeed(start, end - start);
}

std::pair<uint32_t, uint32_t> PgenReader::getChromRange(const std::string& chrom) const {
    const auto& labels = *chrom_labels_view_;
    const auto& ranges = *chrom_ranges_view_;
    for (size_t i = 0; i < labels.size(); i++) {
        if (labels[i] == chrom) return ranges[i];
    }
    return {0, 0};
}

std::string PgenRecordCapabilities::flag_names(unsigned int vrtype_byte) const {
    // Bit positions confirmed empirically against plink2 output.
    std::vector<std::string> names;
    if (vrtype_byte & 0x10u) names.push_back("hardcall phase");
    if (vrtype_byte & 0x20u) names.push_back("phase-set / multiallelic coding");
    if (vrtype_byte & 0x40u) names.push_back("dosage");
    if (vrtype_byte & 0x80u) names.push_back("genotype probabilities");
    const unsigned int other = vrtype_byte & 0x08u;
    if (other) names.push_back("reserved bit 3");
    if (names.empty()) names.push_back("non-hardcall record flags");
    std::string joined;
    for (size_t i = 0; i < names.size(); ++i) {
        if (i) joined += ", ";
        joined += names[i];
    }
    return joined;
}

std::string PgenRecordCapabilities::describe() const {
    if (hardcall_only) return "hardcall-only";
    std::string message = std::to_string(flagged_variants) + " variant record(s) carry " +
                          flag_names(first_flagged_vrtype);
    if (first_flagged_variant >= 0) {
        message += " (first: " + first_flagged_id + " at index " +
                   std::to_string(first_flagged_variant) +
                   ", vrtype=0x" + [this] {
                       char buf[8];
                       std::snprintf(buf, sizeof(buf), "%02X", first_flagged_vrtype & 0xFFu);
                       return std::string(buf);
                   }() + ")";
    }
    return message;
}

PgenRecordCapabilities PgenReader::getRecordCapabilities() const {
    PgenRecordCapabilities capabilities;
    const auto& types = *vrtypes_view_;
    if (types.empty()) return capabilities;
    const auto& variant_meta = *variants_view_;
    for (size_t j = 0; j < types.size(); ++j) {
        // Only the low 3 bits encode the hardcall storage type. Any higher bit
        // means the record carries data this reader does not decode.
        if ((types[j] & 0xF8u) == 0) continue;
        if (capabilities.hardcall_only) {
            capabilities.hardcall_only = false;
            capabilities.first_flagged_variant = static_cast<int>(j);
            capabilities.first_flagged_vrtype = types[j];
            if (j < variant_meta.size()) capabilities.first_flagged_id = variant_meta[j].id;
        }
        ++capabilities.flagged_variants;
    }
    return capabilities;
}

bool PgenReader::decompressZstd(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) const {
    size_t const dSize = ZSTD_decompress(dst.data(), dst_len, src, src_len);
    if (ZSTD_isError(dSize)) {
        std::cerr << "ZSTD Error: " << ZSTD_getErrorName(dSize) << std::endl;
        return false;
    }
    return dSize == dst_len;
}

bool PgenReader::readVariant(std::string& chrom, std::string& rsid, uint32_t& pos, 
                             std::string& ref, std::string& alt, 
                             Eigen::VectorXf& dosage) const {
    if (current_variant_idx >= num_variants) return false;
    if (current_variant_idx >= variants.size()) return false;
    
    const auto& v = variants[current_variant_idx];
    chrom = v.chrom;
    rsid = v.id;
    pos = v.pos;
    ref = v.ref;
    alt = v.alt;
    
    const uint8_t* mmap_data = (const uint8_t*)mmap_file.get_data();
    size_t file_size = mmap_file.get_size();
    
    dosage.resize(num_samples);
    
    if (storage_mode == 2) {
        // Fixed-width mode (format=2): all variants have the same type (vrtypes_[0]) and length (vlens_[0])
        uint8_t type = vrtypes_[0];
        size_t len = vlens_[0];
        size_t offset = variant_block_offset + current_variant_idx * len;

        if (offset + len > file_size) return false;
        const uint8_t* ptr = mmap_data + offset;

        if (type == 0) {
            // Type 0: 00=hom alt, 01=missing, 10=het, 11=hom ref
            for (size_t i = 0; i < num_samples; ++i) {
                size_t byte_idx = i / 4;
                size_t bit_idx = (i % 4) * 2;

                if (byte_idx < len) {
                    uint8_t geno = (ptr[byte_idx] >> bit_idx) & 0x03;
                    if (geno == 0) dosage[i] = 0.0f;
                    else if (geno == 1) dosage[i] = 1.0f;
                    else if (geno == 2) dosage[i] = 2.0f;
                    else dosage[i] = std::numeric_limits<float>::quiet_NaN();
                } else {
                    dosage[i] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        } else if (type == 6) {
            // Type 6 (format=2): 00=hom ref, 01=het, 10=hom alt, 11=missing
            for (size_t i = 0; i < num_samples; ++i) {
                size_t byte_idx = i / 4;
                size_t bit_idx = (i % 4) * 2;

                if (byte_idx < len) {
                    uint8_t geno = (ptr[byte_idx] >> bit_idx) & 0x03;
                    if (geno == 0) dosage[i] = 0.0f;      // hom ref
                    else if (geno == 1) dosage[i] = 1.0f; // het
                    else if (geno == 2) dosage[i] = 2.0f; // hom alt
                    else dosage[i] = std::numeric_limits<float>::quiet_NaN(); // missing
                } else {
                    dosage[i] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        } else {
            dosage.setConstant(std::numeric_limits<float>::quiet_NaN());
        }
    } else {
        const auto& types = *vrtypes_view_;
        const auto& lengths = *vlens_view_;
        const auto& offsets = *variant_offsets_view_;
        if (current_variant_idx >= types.size() || current_variant_idx >= lengths.size() ||
            current_variant_idx >= offsets.size()) return false;
        const uint8_t type = types[current_variant_idx] & 0x07;
        const size_t len = lengths[current_variant_idx];
        const size_t record_offset = variant_block_offset + offsets[current_variant_idx];
        if (record_offset + len > file_size) return false;
        const uint8_t* ptr = mmap_data + record_offset;

        auto apply_difflist = [&](std::vector<uint8_t>& decoded, const uint8_t* cursor) {
            const uint8_t* end = ptr + len;
            auto read_vint = [&]() -> uint32_t {
                uint32_t value = 0;
                int shift = 0;
                while (cursor < end && shift <= 28) {
                    const uint8_t byte = *cursor++;
                    value |= static_cast<uint32_t>(byte & 0x7f) << shift;
                    if (!(byte & 0x80)) return value;
                    shift += 7;
                }
                throw std::runtime_error("Malformed PGEN difference list");
            };
            const uint32_t difference_count = read_vint();
            const uint32_t group_count = (difference_count + 63) / 64;
            uint32_t sample_id_bytes = 1;
            while ((sample_id_bytes < 4) && (num_samples >> (8 * sample_id_bytes))) ++sample_id_bytes;
            const size_t group_info_bytes = difference_count ? group_count * (sample_id_bytes + 1) - 1 : 0;
            if (static_cast<size_t>(end - cursor) < group_info_bytes + (difference_count + 3) / 4) {
                throw std::runtime_error("Truncated PGEN difference list");
            }
            const uint8_t* group_info = cursor;
            cursor += group_info_bytes;
            const uint8_t* rare_genotypes = cursor;
            cursor += (difference_count + 3) / 4;
            uint32_t sample_index = 0;
            for (uint32_t k = 0; k < difference_count; ++k) {
                if ((k % 64) == 0) {
                    sample_index = 0;
                    const uint8_t* first = group_info + (k / 64) * sample_id_bytes;
                    for (uint32_t b = 0; b < sample_id_bytes; ++b) {
                        sample_index |= static_cast<uint32_t>(first[b]) << (8 * b);
                    }
                } else {
                    sample_index += read_vint();
                }
                if (sample_index >= num_samples) throw std::runtime_error("PGEN difference-list sample index out of range");
                decoded[sample_index] = (rare_genotypes[k / 4] >> ((k % 4) * 2)) & 0x03;
            }
        };

        if (type == kPgenTypeUncompressed) {
            ldbase_genotypes_.resize(num_samples);
            for (size_t i = 0; i < num_samples; ++i) {
                size_t byte_idx = i / 4;
                size_t bit_idx = (i % 4) * 2;

                if (byte_idx < len) {
                    uint8_t geno = (ptr[byte_idx] >> bit_idx) & 0x03;

                    ldbase_genotypes_[i] = geno;
                    if (geno <= 2) dosage[i] = static_cast<float>(geno);
                    else dosage[i] = std::numeric_limits<float>::quiet_NaN();
                } else {
                    dosage[i] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        } else if (type == 1) {
            const size_t bitmap_bytes = (num_samples + 7) / 8;
            if (len < 1 + bitmap_bytes) return false;
            const uint8_t common_codes = ptr[0];
            const uint8_t common0 = common_codes >> 2;
            const uint8_t common1 = common0 + (common_codes & 0x03);
            std::vector<uint8_t> decoded(num_samples);
            for (size_t i = 0; i < num_samples; ++i) {
                decoded[i] = (ptr[1 + i / 8] >> (i % 8)) & 1 ? common1 : common0;
            }
            apply_difflist(decoded, ptr + 1 + bitmap_bytes);
            ldbase_genotypes_ = decoded;
            for (size_t i = 0; i < num_samples; ++i) {
                dosage[i] = decoded[i] <= 2 ? static_cast<float>(decoded[i])
                                            : std::numeric_limits<float>::quiet_NaN();
            }
        } else if (type == 2 || type == 3) {
            if (ldbase_genotypes_.size() != num_samples) return false;
            std::vector<uint8_t> decoded = ldbase_genotypes_;
            apply_difflist(decoded, ptr);
            if (type == 3) {
                for (auto& genotype : decoded) {
                    if (genotype == 0) genotype = 2;
                    else if (genotype == 2) genotype = 0;
                }
            }
            for (size_t i = 0; i < num_samples; ++i) {
                dosage[i] = decoded[i] <= 2 ? static_cast<float>(decoded[i])
                                            : std::numeric_limits<float>::quiet_NaN();
            }
        } else if (type >= 4 && type <= 7) {
            const uint8_t base_genotype = type == 5 ? 0 : (type & 0x03);
            std::vector<uint8_t> decoded(num_samples, base_genotype);
            if (type != 5) apply_difflist(decoded, ptr);
            if (type != 5) ldbase_genotypes_ = decoded;
            for (size_t i = 0; i < num_samples; ++i) {
                dosage[i] = decoded[i] <= 2 ? static_cast<float>(decoded[i])
                                            : std::numeric_limits<float>::quiet_NaN();
            }
        } else {
            dosage.setConstant(std::numeric_limits<float>::quiet_NaN());
        }
    }
    
    current_variant_idx++;
    return true;
}

} // namespace cosmic

namespace cosmic {

void PgenReader::exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream geno_out(geno_file);
    if (!geno_out) throw std::runtime_error("Cannot open output geno file: " + geno_file);
    
    geno_out << "VariantID";
    const auto& sample_ids = getSampleIds();
    for (const auto& sid : sample_ids) {
        geno_out << " " << sid;
    }
    geno_out << "\n";
    
    // Need to use non-const reader to iterate
    const_cast<PgenReader*>(this)->reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;
    PgenReader* self = const_cast<PgenReader*>(this);
    while (self->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
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
    self->reset();
    while (self->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
        var_out << rsid << " " << chrom << " " << pos << " " << ref << " " << alt << "\n";
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

void PgenReader::exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const {
    std::ofstream file(npy_file, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open output NPY file: " + npy_file);
    
    std::vector<Eigen::VectorXf> all_dosages;
    std::vector<std::string> variant_ids;
    std::vector<std::string> chroms;
    std::vector<uint32_t> positions;
    std::vector<std::string> refs, alts;
    
    PgenReader* self = const_cast<PgenReader*>(this);
    self->reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (self->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
        all_dosages.push_back(dosage);
        variant_ids.push_back(rsid);
        chroms.push_back(chrom);
        positions.push_back(pos);
        refs.push_back(ref);
        alts.push_back(alt);
    }
    
    uint32_t n_variants = all_dosages.size();
    uint32_t n_samples = num_samples;
    
    file.write("\x93NUMPY", 6);
    file.put(1); file.put(0);
    
    size_t header_len_pos = file.tellp();
    file.write("\x00\x00", 2);
    
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

void PgenReader::exportToPedMap(const std::string& ped_file, const std::string& map_file) const {
    std::ofstream map_out(map_file);
    if (!map_out) throw std::runtime_error("Cannot open output MAP file: " + map_file);
    
    PgenReader* self = const_cast<PgenReader*>(this);
    self->reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;
    std::vector<std::vector<float>> all_dosages;
    while (self->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
        map_out << chrom << " " << rsid << " 0 " << pos << " " << ref << " " << alt << "\n";
        std::vector<float> var_dosage(num_samples);
        for (uint32_t i = 0; i < num_samples; ++i) var_dosage[i] = dosage[i];
        all_dosages.push_back(var_dosage);
    }
    
    std::ofstream ped_out(ped_file);
    if (!ped_out) throw std::runtime_error("Cannot open output PED file: " + ped_file);
    
    const auto& sample_ids = getSampleIds();
    for (uint32_t i = 0; i < num_samples; ++i) {
        ped_out << sample_ids[i] << " " << sample_ids[i] << " 0 0 1 -9";
        for (size_t j = 0; j < all_dosages.size(); ++j) {
            float val = all_dosages[j][i];
            if (std::isnan(val)) {
                ped_out << " 0 0";
            } else {
                if (val < 0.5) ped_out << " 0 0";
                else if (val < 1.5) ped_out << " 1 0";
                else ped_out << " 1 1";
            }
        }
        ped_out << "\n";
    }
}

void PgenReader::exportToVcf(const std::string& vcf_file) const {
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
    
    PgenReader* self = const_cast<PgenReader*>(this);
    self->reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;
    while (self->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
        vcf_out << chrom << "\t" << pos << "\t" << rsid << "\t" << ref << "\t" << alt << "\t.\tPASS\t.\tGT";
        
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
