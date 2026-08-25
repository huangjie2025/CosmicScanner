#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "mmap_matrix.h"

namespace cosmic {

struct PgenVariant {
    std::string chrom;
    std::string id;
    uint32_t pos;
    std::string ref;
    std::string alt;
};

class PgenReader {
public:
    PgenReader();
    ~PgenReader();

    // Open .pgen, .pvar, .psam
    // If pvar/psam are empty, try to infer from pgen prefix
    void open(const std::string& pgen_file, const std::string& pvar_file = "", const std::string& psam_file = "");
    void openSharedMetadata(const PgenReader& source);
    void close();
    void reset();

    uint32_t getNumSamples() const { return num_samples; }
    uint32_t getNumVariants() const { return num_variants; }
    const std::vector<std::string>& getSampleIds() const { return *sample_ids_view_; }

    // Read next variant
    // Returns false if EOF
    // Fills variant metadata and dosage vector (size N)
    // Dosage = expected count of Alt allele (0..2)
    bool readVariant(std::string& chrom, std::string& rsid, uint32_t& pos,
                     std::string& ref, std::string& alt,
                     Eigen::VectorXf& dosage);

    // Compute allele frequencies and missingness stats
    // Required for standardized multiplication
    void computeStats(int n_threads = 1);

    // Matrix-vector multiplication: out = Z * v
    // Z is standardized genotype matrix (n x m)
    // v is vector of size m
    // out is vector of size n
    // Requires computeStats() called first
    void multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1);

    // Matrix-vector multiplication: out = Z' * v
    // Z is standardized genotype matrix (n x m)
    // v is vector of size n
    // out is vector of size m
    // Requires computeStats() called first
    void multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1);

    // Getters for stats
    const std::vector<double>& getMeans() const { return snp_means; }
    const std::vector<double>& getInvStds() const { return snp_inv_stds; }

    // Memory Caching
    // Use Memory Mapping for efficient access
    void cacheAllData() { /* No-op, mmap is already active */ }
    bool isCached() const { return mmap_file.is_open(); }

    // Decode and cache all genotypes as an n×m matrix (column-major)
    // Subsequent multiply_Z_v / multiply_Zt_v calls use the cache
    void cacheAllDosages(int n_threads = 1);
    bool hasDosageCache() const { return dosages_cached_; }
    const Eigen::MatrixXd& getDosageCache() const { return dosages_cache_; }

    // Get Variant Info (BIM-like)
    const std::vector<PgenVariant>& getVariants() const { return *variants_view_; }
    const std::string& getPgenPath() const { return pgen_path; }
    const std::string& getPvarPath() const { return pvar_path; }
    const std::string& getPsamPath() const { return psam_path; }
    void releaseVariantRange(uint32_t start_idx, uint32_t count) const;

    // Seek to a specific variant by index
    void seekVariant(uint32_t idx);

    // Get variant index range for a chromosome
    // Returns (start_idx, count); count=0 if chromosome not found
    std::pair<uint32_t, uint32_t> getChromRange(const std::string& chrom) const;

    // Get list of chromosome labels
    const std::vector<std::string>& getChromLabels() const { return *chrom_labels_view_; }

private:
    std::string pgen_path, pvar_path, psam_path;

    // Memory Mapped File
    MMapFile mmap_file;
    size_t current_offset = 0; // Current read position

    uint32_t num_samples = 0;
    uint32_t num_variants = 0;
    std::vector<std::string> sample_ids;
    std::vector<PgenVariant> variants;
    const std::vector<std::string>* sample_ids_view_ = &sample_ids;
    const std::vector<PgenVariant>* variants_view_ = &variants;

    // Stats
    bool stats_computed = false;
    std::vector<double> snp_means;
    std::vector<double> snp_inv_stds;

    uint32_t current_variant_idx = 0;

    // PGEN Header info
    uint8_t version = 0;
    uint8_t storage_mode = 0;
    size_t variant_block_offset = 0;
    uint32_t n_samples_pgen = 0;
    uint32_t n_variants_pgen = 0;
    uint8_t format_byte = 0;
    bool provisional_ref = false;

    // Directory from header body (mode 0x10+)
    std::vector<uint8_t> vrtypes_;
    std::vector<uint32_t> vlens_;
    std::vector<size_t> variant_offsets_; // cumulative offsets for seek
    const std::vector<uint8_t>* vrtypes_view_ = &vrtypes_;
    const std::vector<uint32_t>* vlens_view_ = &vlens_;
    const std::vector<size_t>* variant_offsets_view_ = &variant_offsets_;

    // Per-chromosome variant index ranges
    std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_; // (start_idx, count) per chrom
    std::vector<std::string> chrom_labels_; // chromosome labels
    const std::vector<std::pair<uint32_t, uint32_t>>* chrom_ranges_view_ = &chrom_ranges_;
    const std::vector<std::string>* chrom_labels_view_ = &chrom_labels_;

    // Reusable genotype buffer
    std::vector<double> gt_buffer_;

    // Dosage cache for fast GRM multiply operations
    bool dosages_cached_ = false;
    Eigen::MatrixXd dosages_cache_; // n × m, column-major: dosages_cache_(i, j) = genotype of sample i, variant j

    // LD compression cache (type 2/3)
    uint32_t ldbase_vidx_ = UINT32_MAX;
    std::vector<double> ldbase_genovec_; // cached genotypes of LD base variant

    // Helpers
    void resetMetadataViews();
    void readPsam();
    void readPvar();
    void readPgenHeader();

    // Parse directory for mode 0x10+
    void parseDirectory(const uint8_t* data, size_t size);

    // Get variant data pointer and advance offset
    const uint8_t* getVariantData(size_t& len);

    // Decode genotypes from raw PGEN data
    void decodeOneVariant(const uint8_t* data_ptr, size_t len, uint8_t vrtype,
                          uint32_t num_samps, double* dosage);
    void decodeOneVariantFloat(const uint8_t* data_ptr, size_t len, uint8_t vrtype,
                               uint32_t num_samps, float* dosage);

    // Decoder helpers
    static uint32_t readVarint(const uint8_t*& ptr);
    static uint32_t parseDifflistCounts(const uint8_t* data, size_t len,
                                         std::vector<uint32_t>& sample_ids,
                                         std::vector<uint8_t>& values,
                                         uint32_t num_samples);

    // LD compression helpers (type 2/3)
    // Find the LD base variant index for a given variant
    uint32_t findLdBase(uint32_t vidx) const;
    // Load LD base variant's genotypes into ldbase_genovec_, updating ldbase_vidx_
    void loadLdBase(uint32_t vidx);
    // Apply difflist on top of existing genotypes (used for type 2/3)
    void applyDifflistToBuffer(const uint8_t* data_ptr, size_t len,
                                uint32_t num_samps, double* dosage);
    // Invert genotypes in-place (0↔2, for type 3)
    static void invertGenotypes(uint32_t num_samps, double* dosage);
    static void invertGenotypesFloat(uint32_t num_samps, float* dosage);

    // Decompression buffer
    std::vector<char> compressed_buffer;
    std::vector<char> decompressed_buffer;

    // Helper to decompress Zstd
    // Returns true on success
    bool decompressZstd(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len);
};

} // namespace cosmic
