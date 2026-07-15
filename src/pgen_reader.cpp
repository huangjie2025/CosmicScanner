#include "pgen_reader.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <zstd.h>
#include <cstring>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <omp.h>

namespace cosmic {

static const uint16_t kPgenMagic = 0x1B6C;

// Read a base-128 varint (protocol-buffers style)
uint32_t PgenReader::readVarint(const uint8_t*& ptr) {
    uint32_t result = 0;
    int shift = 0;
    while (true) {
        uint8_t b = *ptr++;
        result |= (uint32_t)(b & 0x7F) << shift;
        if (!(b & 0x80)) break;
        shift += 7;
    }
    return result;
}

// Parse a difflist from PGEN data.
// Returns bytes consumed, or 0 on error/failure.
// sample_ids and values are appended with up to num_samples entries.
uint32_t PgenReader::parseDifflistCounts(const uint8_t* data, size_t len,
                                          std::vector<uint32_t>& sample_ids,
                                          std::vector<uint8_t>& values,
                                          uint32_t num_samples) {
    const uint8_t* start = data;
    if (len < 1) return 0;

    const uint8_t* ptr = data;
    uint32_t L = readVarint(ptr);
    if (L == 0) return (uint32_t)(ptr - start);
    if (L > num_samples * 2) return 0;

    // Determine sample ID byte width
    uint32_t id_width;
    if (num_samples <= 256) id_width = 1;
    else if (num_samples <= 65536) id_width = 2;
    else if (num_samples <= 16777216) id_width = 3;
    else id_width = 4;

    uint32_t G = (L + 63) / 64;

    // Read group-start IDs (G entries, each id_width bytes)
    size_t needed = (size_t)G * id_width;
    if ((size_t)(ptr - start) + needed > len) return 0;

    std::vector<uint32_t> group_starts(G);
    for (uint32_t g = 0; g < G; g++) {
        uint32_t id = 0;
        for (uint32_t b = 0; b < id_width; b++) {
            id |= (uint32_t)ptr[b] << (b * 8);
        }
        ptr += id_width;
        group_starts[g] = id;
    }

    // Read group sizes (G-1 bytes)
    needed = (G > 0) ? (G - 1) : 0;
    if ((size_t)(ptr - start) + needed > len) return 0;

    std::vector<uint32_t> group_sizes(G);
    for (uint32_t g = 0; g < G - 1; g++) {
        group_sizes[g] = ptr[0] + 63;
        ptr++;
    }
    // Last group size is implicit (will be computed)

    // Check if values are stored (2-bit genotype category values)
    bool has_values = true; // For difflist types 1, 4, 6, 7, values ARE stored
    uint32_t val_bytes = (L + 3) / 4;

    size_t remaining_start = (size_t)(ptr - start);

    // Read 2-bit values if present
    const uint8_t* val_ptr = ptr;
    if (has_values) {
        ptr += val_bytes;
    }

    // Read delta varints for remaining entries
    for (uint32_t g = 0; g < G; g++) {
        uint32_t first_in_group = (g == 0) ? 0 : (g * 64);
        uint32_t last_in_group = std::min(first_in_group + 63, L - 1);
        uint32_t count_in_group = last_in_group - first_in_group + 1;

        if (count_in_group < 1) continue;

        // First ID in group is from group_starts[g]
        if (sample_ids.size() < L) {
            sample_ids.push_back(group_starts[g]);
        }

        // Read deltas for remaining entries in group
        for (uint32_t k = first_in_group + 1; k <= last_in_group; k++) {
            if ((size_t)(ptr - start) >= len) return 0;
            uint32_t delta = readVarint(ptr);
            uint32_t prev = sample_ids.empty() ? group_starts[g] : sample_ids.back();
            if (sample_ids.size() < L) {
                sample_ids.push_back(prev + delta);
            }
        }
    }

    // Extract values
    for (uint32_t i = 0; i < L && i < sample_ids.size(); i++) {
        uint8_t v = (val_ptr[i / 4] >> ((i % 4) * 2)) & 3;
        values.push_back(v);
    }

    return (uint32_t)(ptr - start);
}

void PgenReader::parseDirectory(const uint8_t* data, size_t size) {
    uint32_t fmt = (format_byte & 0x0F);
    if (fmt > 7) throw std::runtime_error("PGEN unsupported variant-record type/length format");
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

PgenReader::PgenReader() {}

PgenReader::~PgenReader() {
    close();
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

void PgenReader::open(const std::string& pgen_file, const std::string& pvar_file, const std::string& psam_file) {
    resetMetadataViews();
    pgen_path = pgen_file;
    pvar_path = pvar_file.empty() ? pgen_file.substr(0, pgen_file.find_last_of('.')) + ".pvar" : pvar_file;
    psam_path = psam_file.empty() ? pgen_file.substr(0, pgen_file.find_last_of('.')) + ".psam" : psam_file;

    mmap_file.open(pgen_path);
    if (!mmap_file.is_open()) {
        throw std::runtime_error("Failed to open PGEN file: " + pgen_path);
    }

    readPgenHeader();
    readPvar();
    readPsam();

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
    gt_buffer_.resize(num_samples);
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

void PgenReader::close() {
    mmap_file.close();
    variants.clear();
    sample_ids.clear();
    snp_means.clear();
    snp_inv_stds.clear();
    vrtypes_.clear();
    vlens_.clear();
    variant_offsets_.clear();
    chrom_ranges_.clear();
    chrom_labels_.clear();
    resetMetadataViews();
    stats_computed = false;
}

void PgenReader::reset() {
    current_offset = variant_block_offset;
    current_variant_idx = 0;
}

void PgenReader::seekVariant(uint32_t idx) {
    if (idx >= num_variants) {
        current_variant_idx = num_variants;
        const auto& offs = *variant_offsets_view_;
        const auto& lens = *vlens_view_;
        current_offset = variant_block_offset + (offs.empty() ? 0 : offs.back() + lens.back());
        return;
    }
    const auto& offs = *variant_offsets_view_;
    const auto& lens = *vlens_view_;
    if (!offs.empty() && idx < offs.size()) {
        current_offset = variant_block_offset + offs[idx];
    } else {
        // Fallback: compute from start (shouldn't happen if offsets are built)
        current_offset = variant_block_offset;
        for (uint32_t i = 0; i < idx; i++) {
            if (i < lens.size()) current_offset += lens[i];
        }
    }
    current_variant_idx = idx;
}

std::pair<uint32_t, uint32_t> PgenReader::getChromRange(const std::string& chrom) const {
    const auto& labels = *chrom_labels_view_;
    const auto& ranges = *chrom_ranges_view_;
    for (size_t i = 0; i < labels.size(); i++) {
        if (labels[i] == chrom) {
            return ranges[i];
        }
    }
    return {0, 0};
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

void PgenReader::readPgenHeader() {
    const uint8_t* data = (const uint8_t*)mmap_file.get_data();
    size_t size = mmap_file.get_size();

    if (size < 12) throw std::runtime_error("PGEN file too small");

    size_t offset = 0;
    if (data[0] != 0x6C || data[1] != 0x1B) {
        throw std::runtime_error("Invalid PGEN magic number");
    }
    offset += 2;

    storage_mode = data[offset++];

    n_variants_pgen = *reinterpret_cast<const uint32_t*>(data + offset);
    offset += 4;

    n_samples_pgen = *reinterpret_cast<const uint32_t*>(data + offset);
    offset += 4;

    if (storage_mode == 0x01) {
        // PLINK 1 BED mode
        variant_block_offset = 0;
        current_offset = 0;
        return;
    }

    // All modes after 0x01 have format_byte at current offset
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

// Return pointer to current variant's data, advance offset
// Sets len to the data length
const uint8_t* PgenReader::getVariantData(size_t& len) {
    const uint8_t* data = (const uint8_t*)mmap_file.get_data();
    size_t file_size = mmap_file.get_size();

    if (current_offset >= file_size) return nullptr;

    uint32_t vlen = (*vlens_view_)[current_variant_idx];
    if (current_offset + vlen > file_size) return nullptr;

    const uint8_t* ptr = data + current_offset;
    current_offset += vlen;
    current_variant_idx++;
    len = vlen;
    return ptr;
}

// Decode one variant's data into double genotype values (dosage 0..2)
uint32_t PgenReader::findLdBase(uint32_t vidx) const {
    // Scan backward for the most recent non-LD-compressed variant
    // Non-LD: (vrtype & 6) != 2
    int64_t base = (int64_t)vidx - 1;
    const auto& vrtypes = *vrtypes_view_;
    while (base >= 0) {
        if ((vrtypes[(uint32_t)base] & 6) != 2) return (uint32_t)base;
        base--;
    }
    // Should never reach here: first variant in each vblock is guaranteed non-LD
    return UINT32_MAX;
}

void PgenReader::loadLdBase(uint32_t vidx) {
    uint32_t base_vidx = findLdBase(vidx);
    if (base_vidx == UINT32_MAX) {
        throw std::runtime_error("PGEN: Cannot find LD base variant for variant " + std::to_string(vidx));
    }
    if (base_vidx == ldbase_vidx_ && !ldbase_genovec_.empty()) return; // Already cached

    // Seek to base variant and decode it
    size_t saved_offset = current_offset;
    uint32_t saved_variant_idx = current_variant_idx;
    seekVariant(base_vidx);
    size_t len;
    const uint8_t* ptr = getVariantData(len);
    if (!ptr) {
        current_offset = saved_offset;
        current_variant_idx = saved_variant_idx;
        throw std::runtime_error("PGEN: Failed to read LD base variant data");
    }

    ldbase_genovec_.resize(num_samples);
    uint8_t base_vrtype = (*vrtypes_view_)[base_vidx];
    decodeOneVariant(ptr, len, base_vrtype, num_samples, ldbase_genovec_.data());
    ldbase_vidx_ = base_vidx;
    current_offset = saved_offset;
    current_variant_idx = saved_variant_idx;
}

void PgenReader::applyDifflistToBuffer(const uint8_t* data_ptr, size_t len,
                                         uint32_t num_samps, double* dosage) {
    std::vector<uint32_t> diff_samples;
    std::vector<uint8_t> diff_vals;
    uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, num_samps);
    if (consumed > 0) {
        for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
            uint32_t sid = diff_samples[di];
            if (sid < num_samps) {
                uint8_t val = diff_vals[di];
                if (val == 0) dosage[sid] = 0.0;
                else if (val == 1) dosage[sid] = 1.0;
                else if (val == 2) dosage[sid] = 2.0;
                else dosage[sid] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}

void PgenReader::invertGenotypes(uint32_t num_samps, double* dosage) {
    for (uint32_t i = 0; i < num_samps; ++i) {
        if (dosage[i] == 0.0) dosage[i] = 2.0;
        else if (dosage[i] == 2.0) dosage[i] = 0.0;
        // 1.0 and NaN stay unchanged
    }
}

void PgenReader::invertGenotypesFloat(uint32_t num_samps, float* dosage) {
    for (uint32_t i = 0; i < num_samps; ++i) {
        if (dosage[i] == 0.0f) dosage[i] = 2.0f;
        else if (dosage[i] == 2.0f) dosage[i] = 0.0f;
    }
}

void PgenReader::decodeOneVariant(const uint8_t* data_ptr, size_t len, uint8_t vrtype,
                                   uint32_t num_samps, double* dosage) {
    uint32_t n = num_samps;
    uint32_t main_type = vrtype & 7;

    switch (main_type) {
        case 0: { // Simple 2-bit
            uint32_t nbytes = (n + 3) / 4;
            for (uint32_t i = 0; i < n; i++) {
                if (i / 4 < nbytes) {
                    uint8_t v = (data_ptr[i / 4] >> ((i % 4) * 2)) & 3;
                    if (v == 0) dosage[i] = 0.0;
                    else if (v == 1) dosage[i] = 1.0;
                    else if (v == 2) dosage[i] = 2.0;
                    else dosage[i] = std::numeric_limits<double>::quiet_NaN();
                } else {
                    dosage[i] = std::numeric_limits<double>::quiet_NaN();
                }
            }
            break;
        }
        case 1: { // 1-bit encoding
            if (len < 1) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
                break;
            }

            const uint8_t* ptr = data_ptr;
            uint8_t start = *ptr++;

            uint8_t cat_a, cat_b;
            switch (start) {
                case 1: cat_a = 0; cat_b = 1; break;
                case 2: cat_a = 0; cat_b = 2; break;
                case 3: cat_a = 0; cat_b = 3; break;
                case 5: cat_a = 1; cat_b = 2; break;
                case 6: cat_a = 1; cat_b = 3; break;
                case 9: cat_a = 2; cat_b = 3; break;
                default:
                    for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
                    return;
            }

            uint32_t bit_nbytes = (n + 7) / 8;
            if ((size_t)(ptr - data_ptr) + bit_nbytes > len) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
                break;
            }

            for (uint32_t i = 0; i < n; i++) {
                bool in_cat_b = (ptr[i / 8] >> (i % 8)) & 1;
                uint8_t val = in_cat_b ? cat_b : cat_a;
                dosage[i] = val <= 2 ? (double)val : std::numeric_limits<double>::quiet_NaN();
            }
            ptr += bit_nbytes;

            size_t remaining = len - (ptr - data_ptr);
            if (remaining > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(ptr, remaining, diff_samples, diff_vals, n);
                if (consumed > 0 && diff_samples.size() > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            dosage[sid] = val <= 2 ? (double)val : std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }

            break;
        }
        case 4: { // Difflist non-cat-0
            // Default category 0 for all samples
            for (uint32_t i = 0; i < n; i++) dosage[i] = 0.0;
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0;
                            else if (val == 1) dosage[sid] = 1.0;
                            else if (val == 2) dosage[sid] = 2.0;
                            else dosage[sid] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 6: { // Difflist non-cat-2
            for (uint32_t i = 0; i < n; i++) dosage[i] = 2.0;
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0;
                            else if (val == 1) dosage[sid] = 1.0;
                            else if (val == 2) dosage[sid] = 2.0;
                            else dosage[sid] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 7: { // Difflist non-cat-3 (missing)
            for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0;
                            else if (val == 1) dosage[sid] = 1.0;
                            else if (val == 2) dosage[sid] = 2.0;
                            else dosage[sid] = std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 5: { // All category 1 (het)
            for (uint32_t i = 0; i < n; i++) dosage[i] = 1.0;
            break;
        }
        case 2:
        case 3: {
            // LD compression: start with LD base, apply difflist, optionally invert
            if (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != n) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
                break;
            }
            std::copy(ldbase_genovec_.begin(), ldbase_genovec_.end(), dosage);
            if (len > 0) {
                applyDifflistToBuffer(data_ptr, len, n, dosage);
            }
            if (vrtype == 3) {
                invertGenotypes(n, dosage);
            }
            break;
        }
        default: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<double>::quiet_NaN();
            break;
        }
    }
}

void PgenReader::decodeOneVariantFloat(const uint8_t* data_ptr, size_t len, uint8_t vrtype,
                                        uint32_t num_samps, float* dosage) {
    uint32_t n = num_samps;
    uint32_t main_type = vrtype & 7;

    switch (main_type) {
        case 0: {
            uint32_t nbytes = (n + 3) / 4;
            for (uint32_t i = 0; i < n; i++) {
                if (i / 4 < nbytes) {
                    uint8_t v = (data_ptr[i / 4] >> ((i % 4) * 2)) & 3;
                    if (v == 0) dosage[i] = 0.0f;
                    else if (v == 1) dosage[i] = 1.0f;
                    else if (v == 2) dosage[i] = 2.0f;
                    else dosage[i] = std::numeric_limits<float>::quiet_NaN();
                } else {
                    dosage[i] = std::numeric_limits<float>::quiet_NaN();
                }
            }
            break;
        }
        case 1: {
            if (len < 1) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
                break;
            }
            const uint8_t* ptr = data_ptr;
            uint8_t start = *ptr++;
            uint8_t cat_a, cat_b;
            switch (start) {
                case 1: cat_a = 0; cat_b = 1; break;
                case 2: cat_a = 0; cat_b = 2; break;
                case 3: cat_a = 0; cat_b = 3; break;
                case 5: cat_a = 1; cat_b = 2; break;
                case 6: cat_a = 1; cat_b = 3; break;
                case 9: cat_a = 2; cat_b = 3; break;
                default:
                    for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
                    return;
            }
            uint32_t bit_nbytes = (n + 7) / 8;
            if ((size_t)(ptr - data_ptr) + bit_nbytes > len) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
                break;
            }
            for (uint32_t i = 0; i < n; i++) {
                bool in_cat_b = (ptr[i / 8] >> (i % 8)) & 1;
                uint8_t val = in_cat_b ? cat_b : cat_a;
                dosage[i] = val <= 2 ? (float)val : std::numeric_limits<float>::quiet_NaN();
            }
            ptr += bit_nbytes;
            size_t remaining = len - (ptr - data_ptr);
            if (remaining > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(ptr, remaining, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            dosage[sid] = val <= 2 ? (float)val : std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 4: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = 0.0f;
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0f;
                            else if (val == 1) dosage[sid] = 1.0f;
                            else if (val == 2) dosage[sid] = 2.0f;
                            else dosage[sid] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 6: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = 2.0f;
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0f;
                            else if (val == 1) dosage[sid] = 1.0f;
                            else if (val == 2) dosage[sid] = 2.0f;
                            else dosage[sid] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 7: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0f;
                            else if (val == 1) dosage[sid] = 1.0f;
                            else if (val == 2) dosage[sid] = 2.0f;
                            else dosage[sid] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }
            break;
        }
        case 5: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = 1.0f;
            break;
        }
        case 2:
        case 3: {
            // LD compression: requires base in float buffer
            if (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != n) {
                for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
                break;
            }
            for (uint32_t i = 0; i < n; i++) dosage[i] = (float)ldbase_genovec_[i];
            if (len > 0) {
                std::vector<uint32_t> diff_samples;
                std::vector<uint8_t> diff_vals;
                uint32_t consumed = parseDifflistCounts(data_ptr, len, diff_samples, diff_vals, n);
                if (consumed > 0) {
                    for (size_t di = 0; di < diff_samples.size() && di < diff_vals.size(); di++) {
                        uint32_t sid = diff_samples[di];
                        if (sid < n) {
                            uint8_t val = diff_vals[di];
                            if (val == 0) dosage[sid] = 0.0f;
                            else if (val == 1) dosage[sid] = 1.0f;
                            else if (val == 2) dosage[sid] = 2.0f;
                            else dosage[sid] = std::numeric_limits<float>::quiet_NaN();
                        }
                    }
                }
            }
            if (vrtype == 3) {
                invertGenotypesFloat(n, dosage);
            }
            break;
        }
        default: {
            for (uint32_t i = 0; i < n; i++) dosage[i] = std::numeric_limits<float>::quiet_NaN();
            break;
        }
    }
}

void PgenReader::computeStats(int n_threads) {
    if (stats_computed) return;
    if (num_variants == 0) return;

    // Pre-validate directory
    if (vlens_.size() != num_variants) {
        throw std::runtime_error("PGEN directory not initialized");
    }

    reset();
    snp_means.assign(num_variants, 0.0);
    snp_inv_stds.assign(num_variants, 0.0);

    std::cout << "Computing stats for " << num_variants << " variants in PGEN..." << std::endl;

    ldbase_vidx_ = UINT32_MAX;
    for (uint32_t j = 0; j < num_variants; ++j) {
        size_t len;
        const uint8_t* ptr = getVariantData(len);
        if (!ptr) break;

        uint8_t vrtype = vrtypes_[j];

        // Ensure LD base is loaded for type 2/3
        if ((vrtype & 6) == 2) {
            if (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != num_samples) {
                loadLdBase(j);
            }
        }

        decodeOneVariant(ptr, len, vrtype, num_samples, gt_buffer_.data());

        // Cache non-LD variant as potential LD base
        if ((vrtype & 6) != 2) {
            ldbase_vidx_ = j;
            ldbase_genovec_ = gt_buffer_;
        }

        double sum = 0.0;
        double sum_sq = 0.0;
        int64_t count = 0;

        for (uint32_t i = 0; i < num_samples; ++i) {
            double val = gt_buffer_[i];
            if (!std::isnan(val)) {
                sum += val;
                sum_sq += val * val;
                count++;
            }
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
    }

    stats_computed = true;
    reset();
}

void PgenReader::cacheAllDosages(int n_threads) {
    if (dosages_cached_) return;
    if (!stats_computed) computeStats(n_threads);

    std::cout << "  Caching all PGEN genotypes (" << num_samples << " x " << num_variants << ")..." << std::endl;
    dosages_cache_.resize(num_samples, num_variants);

    reset();
    ldbase_vidx_ = UINT32_MAX;
    for (uint32_t j = 0; j < num_variants; ++j) {
        size_t len;
        const uint8_t* ptr = getVariantData(len);
        if (!ptr) break;

        uint8_t vrtype = vrtypes_[j];

        if ((vrtype & 6) == 2 && (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != num_samples)) {
            loadLdBase(j);
        }

        decodeOneVariant(ptr, len, vrtype, num_samples, gt_buffer_.data());

        if ((vrtype & 6) != 2) {
            ldbase_vidx_ = j;
            ldbase_genovec_ = gt_buffer_;
        }

        for (uint32_t i = 0; i < num_samples; ++i) {
            dosages_cache_(i, j) = gt_buffer_[i];
        }
    }
    dosages_cached_ = true;
    reset();
}

void PgenReader::multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_variants) throw std::runtime_error("Dimension mismatch in multiply_Z_v");

    if (dosages_cached_) {
        // Use cache: out = (dosages_cache_ - mean * 1^T) * diag(inv_stds) * v
        Eigen::VectorXd w(num_samples); w.setConstant(0.0);
        #pragma omp parallel for num_threads(n_threads)
        for (long long i = 0; i < (long long)num_samples; ++i) {
            double sum = 0.0;
            for (uint32_t j = 0; j < num_variants; ++j) {
                double g = dosages_cache_(i, j);
                if (!std::isnan(g)) {
                    sum += (g - snp_means[j]) * snp_inv_stds[j] * v(j);
                }
            }
            w(i) = sum;
        }
        out = w;
        return;
    }

    reset();
    out.setZero(num_samples);
    ldbase_vidx_ = UINT32_MAX;

    for (uint32_t j = 0; j < num_variants; ++j) {
        size_t len;
        const uint8_t* ptr = getVariantData(len);
        if (!ptr) break;

        if (snp_inv_stds[j] != 0.0 && !std::isnan(snp_means[j])) {
            uint8_t vrtype = vrtypes_[j];

            if ((vrtype & 6) == 2 && (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != num_samples)) {
                loadLdBase(j);
            }

            decodeOneVariant(ptr, len, vrtype, num_samples, gt_buffer_.data());

            if ((vrtype & 6) != 2) {
                ldbase_vidx_ = j;
                ldbase_genovec_ = gt_buffer_;
            }

            double mean = snp_means[j];
            double weight = v(j) * snp_inv_stds[j];

            #pragma omp parallel num_threads(n_threads)
            {
                Eigen::VectorXd local_out = Eigen::VectorXd::Zero(num_samples);

                #pragma omp for
                for (long long i = 0; i < (long long)num_samples; ++i) {
                    double g = gt_buffer_[i];
                    if (!std::isnan(g)) {
                        local_out(i) += (g - mean) * weight;
                    }
                }

                #pragma omp critical
                {
                    out += local_out;
                }
            }
        }
    }
    reset();
}

void PgenReader::multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads) {
    if (!stats_computed) computeStats(n_threads);
    if (v.size() != (long long)num_samples) throw std::runtime_error("Dimension mismatch in multiply_Zt_v");

    if (dosages_cached_) {
        // Use cache: out(j) = (col_j(dosages_cache_) · v - mean_j * sum(v)) * inv_std_j
        double sum_v = v.sum();
        out.resize(num_variants);
        #pragma omp parallel for num_threads(n_threads)
        for (long long j = 0; j < (long long)num_variants; ++j) {
            if (snp_inv_stds[j] != 0.0 && !std::isnan(snp_means[j])) {
                double dot = 0.0;
                for (uint32_t i = 0; i < num_samples; ++i) {
                    double g = dosages_cache_(i, j);
                    if (!std::isnan(g)) {
                        dot += g * v(i);
                    }
                }
                out(j) = (dot - snp_means[j] * sum_v) * snp_inv_stds[j];
            } else {
                out(j) = 0.0;
            }
        }
        return;
    }

    reset();
    out.setZero(num_variants);
    ldbase_vidx_ = UINT32_MAX;

    for (uint32_t j = 0; j < num_variants; ++j) {
        size_t len;
        const uint8_t* ptr = getVariantData(len);
        if (!ptr) break;

        if (snp_inv_stds[j] != 0.0 && !std::isnan(snp_means[j])) {
            uint8_t vrtype = vrtypes_[j];

            if ((vrtype & 6) == 2 && (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != num_samples)) {
                loadLdBase(j);
            }

            decodeOneVariant(ptr, len, vrtype, num_samples, gt_buffer_.data());

            if ((vrtype & 6) != 2) {
                ldbase_vidx_ = j;
                ldbase_genovec_ = gt_buffer_;
            }

            double mean = snp_means[j];
            double inv_std = snp_inv_stds[j];
            double dot_gv = 0.0;
            double sum_v = 0.0;

            #pragma omp parallel num_threads(n_threads) reduction(+:dot_gv, sum_v)
            {
                double local_dot = 0.0;
                double local_sum = 0.0;

                #pragma omp for
                for (long long i = 0; i < (long long)num_samples; ++i) {
                    double g = gt_buffer_[i];
                    if (!std::isnan(g)) {
                        local_dot += g * v(i);
                        local_sum += v(i);
                    }
                }
                dot_gv += local_dot;
                sum_v += local_sum;
            }

            out(j) = (dot_gv - mean * sum_v) * inv_std;
        } else {
            out(j) = 0.0;
        }
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
        ss >> fid >> iid;
        sample_ids.push_back(fid + ":" + iid);
    }
}

void PgenReader::readPvar() {
    std::ifstream in(pvar_path);
    if (!in) return;

    variants.clear();
    chrom_ranges_.clear();
    chrom_labels_.clear();

    std::string line;
    std::string current_chrom;
    uint32_t chrom_start = 0;
    uint32_t chrom_count = 0;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        PgenVariant v;
        ss >> v.chrom >> v.pos >> v.id >> v.ref >> v.alt;

        if (v.chrom != current_chrom) {
            if (chrom_count > 0) {
                chrom_ranges_.emplace_back(chrom_start, chrom_count);
                chrom_labels_.push_back(current_chrom);
            }
            current_chrom = v.chrom;
            chrom_start = (uint32_t)variants.size();
            chrom_count = 0;
        }
        chrom_count++;
        variants.push_back(v);
    }
    // Last chromosome
    if (chrom_count > 0) {
        chrom_ranges_.emplace_back(chrom_start, chrom_count);
        chrom_labels_.push_back(current_chrom);
    }
}

bool PgenReader::decompressZstd(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) {
    size_t const dSize = ZSTD_decompress(dst.data(), dst_len, src, src_len);
    if (ZSTD_isError(dSize)) {
        std::cerr << "ZSTD Error: " << ZSTD_getErrorName(dSize) << std::endl;
        return false;
    }
    return dSize == dst_len;
}

bool PgenReader::readVariant(std::string& chrom, std::string& rsid, uint32_t& pos,
                             std::string& ref, std::string& alt,
                             Eigen::VectorXf& dosage) {
    if (current_variant_idx >= num_variants) return false;
    const auto& variants_ref = *variants_view_;
    const auto& vrtypes_ref = *vrtypes_view_;
    if (current_variant_idx >= variants_ref.size()) return false;

    uint32_t vidx = current_variant_idx;

    const auto& v = variants_ref[vidx];
    chrom = v.chrom;
    rsid = v.id;
    pos = v.pos;
    ref = v.ref;
    alt = v.alt;

    size_t len;
    const uint8_t* ptr = getVariantData(len);
    if (!ptr) return false;

    uint8_t vrtype = vrtypes_ref[vidx];
    dosage.resize(num_samples);

    if ((vrtype & 6) == 2) {
        if (ldbase_vidx_ == UINT32_MAX || ldbase_genovec_.size() != num_samples) {
            loadLdBase(vidx);
        }
    }

    decodeOneVariantFloat(ptr, len, vrtype, num_samples, dosage.data());

    if ((vrtype & 6) != 2) {
        ldbase_vidx_ = vidx;
        ldbase_genovec_.assign(dosage.data(), dosage.data() + num_samples);
    }

    return true;
}

} // namespace cosmic
