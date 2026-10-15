#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "mmap_matrix.h"
#include "matrix_adapter.h"

namespace cosmic {

struct PgenVariant {
    std::string chrom;
    std::string id;
    uint32_t pos;
    std::string ref;
    std::string alt;
};

// Per-file summary of what the PGEN genotype records actually carry.
//
// A PGEN variant-record type byte uses the low 3 bits for the hardcall
// storage type (0-7) and reserves every higher bit for optional per-variant
// payloads. Verified against plink2 output: BED-derived hardcall files only
// ever set the low 3 bits (0x00-0x07), a phased hardcall file sets 0x10, and a
// hardcall+dosage file sets 0x60. This reader decodes only the hardcall track,
// so a record with any higher bit set would otherwise be silently reduced to
// its hardcalls, discarding phase/dosage/probability information.
struct PgenRecordCapabilities {
    bool hardcall_only = true;
    unsigned int flagged_variants = 0;   // records with bits above the 3-bit storage type
    int first_flagged_variant = -1;      // 0-based variant index, -1 if none
    std::string first_flagged_id;        // variant ID at first_flagged_variant
    unsigned int first_flagged_vrtype = 0;  // raw record-type byte at first_flagged_variant

    // Human-readable flag names, derived from the observed bit positions.
    std::string flag_names(unsigned int vrtype_byte) const;

    std::string describe() const;
};

class PgenReader : public GenotypeMatrix {
public:
    PgenReader();
    ~PgenReader();

    // Open .pgen, .pvar, .psam
    // If pvar/psam are empty, try to infer from pgen prefix
    void open(const std::string& pgen_file, const std::string& pvar_file = "", const std::string& psam_file = "");
    void openSharedMetadata(const PgenReader& source);
    void close();
    void reset() const;

    uint32_t getNumSamples() const { return num_samples; }
    uint32_t getNumVariants() const { return num_variants; }
    const std::vector<std::string>& getSampleIds() const { return sample_ids_view_ ? *sample_ids_view_ : sample_ids; }
    int rows() const override { return static_cast<int>(num_samples); }
    int cols() const override { return static_cast<int>(num_variants); }

    // Read next variant
    // Returns false if EOF
    // Fills variant metadata and dosage vector (size N)
    // Dosage = expected count of Alt allele (0..2)
    bool readVariant(std::string& chrom, std::string& rsid, uint32_t& pos, 
                     std::string& ref, std::string& alt, 
                     Eigen::VectorXf& dosage) const;
    
    // Compute allele frequencies and missingness stats
    // Required for standardized multiplication
    void computeStats(int n_threads = 1);

    // Matrix-vector multiplication: out = Z * v
    // Z is standardized genotype matrix (n x m)
    // v is vector of size m
    // out is vector of size n
    // Requires computeStats() called first
    void multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;

    // Matrix-vector multiplication: out = Z' * v
    // Z is standardized genotype matrix (n x m)
    // v is vector of size n
    // out is vector of size m
    // Requires computeStats() called first
    void multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;

    // Getters for stats
    const std::vector<double>& getMeans() const { return snp_means; }
    const std::vector<double>& getInvStds() const { return snp_inv_stds; }
    
    // Memory Caching
    // Use Memory Mapping for efficient access
    void cacheAllData() { /* No-op, mmap is already active */ }
    bool isCached() const { return mmap_file.is_open(); }

    // Get Variant Info (BIM-like)
    const std::vector<PgenVariant>& getVariants() const { return variants_view_ ? *variants_view_ : variants; }
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
    const std::vector<std::string>& getChromLabels() const { return chrom_labels_view_ ? *chrom_labels_view_ : chrom_labels_; }

    // Report whether every genotype record is a plain hardcall record.
    // Callers that require exact hardcalls (e.g. dominance / additive-dominance
    // GWAS) must reject a file whose records carry phase, dosage, or genotype
    // probabilities rather than silently using only the hardcall track.
    PgenRecordCapabilities getRecordCapabilities() const;

private:
    std::string pgen_path, pvar_path, psam_path;
    
    // Memory Mapped File
    mutable MMapFile mmap_file;
    mutable size_t current_offset = 0; // Current read position
    
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
    
    mutable uint32_t current_variant_idx = 0;
    
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
    std::vector<size_t> variant_offsets_;
    const std::vector<uint8_t>* vrtypes_view_ = &vrtypes_;
    const std::vector<uint32_t>* vlens_view_ = &vlens_;
    const std::vector<size_t>* variant_offsets_view_ = &variant_offsets_;

    // Per-chromosome variant index ranges
    std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_;
    std::vector<std::string> chrom_labels_;
    const std::vector<std::pair<uint32_t, uint32_t>>* chrom_ranges_view_ = &chrom_ranges_;
    const std::vector<std::string>* chrom_labels_view_ = &chrom_labels_;

    // Reusable genotype buffer
    mutable std::vector<double> gt_buffer_;
    mutable std::vector<uint8_t> ldbase_genotypes_;
    
    // Helpers
    void resetMetadataViews();
    void readPsam();
    void readPvar();
    void readPgenHeader();
    void parseDirectory(const uint8_t* data, size_t size);
    
    // Decompression buffer
    mutable std::vector<char> compressed_buffer;
    mutable std::vector<char> decompressed_buffer;
    
    // Helper to decompress Zstd
    // Returns true on success
    bool decompressZstd(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) const;
    
public:
    // Export genotype data to various formats
    void exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToPedMap(const std::string& ped_file, const std::string& map_file) const;
    void exportToVcf(const std::string& vcf_file) const;
};

} // namespace cosmic
