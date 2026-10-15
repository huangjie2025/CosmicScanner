#ifndef COSMIC_BGEN_READER_H
#define COSMIC_BGEN_READER_H

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <Eigen/Dense>
#include "matrix_adapter.h"

namespace cosmic {

// Per-variant biallelic diploid genotype probability block.
//
// Orientation: G is the number of copies of the FIRST STORED BGEN allele.
// NOTE: readVariant()/readVariantProbs() return the stored alleles swapped
// (a1 = second stored allele, a2 = first stored allele), so G counts the allele
// in the SECOND out-parameter `a2`, not `a1`. `G` is defined here against the
// same allele that the existing production dosage path counts, which is exactly
// what makes the probability API dosage-identical to the legacy dosage API.
// Consumers that assume dosage counts `a1` are inverted; see blocker 6 in
// audit/BGEN_PROBABILITY_BACKEND_IMPLEMENTATION_REPORT.md.
//
//   p0 = P(G=0) = P(non-counted / non-counted)
//   p1 = P(G=1) = P(counted / non-counted)   <- heterozygosity probability
//   p2 = P(G=2) = P(counted / counted)
//
// Therefore E[G] = p1 + 2*p2, matching the frozen Phase E2.3 definition.
struct BgenGenotypeProbs {
    Eigen::VectorXf p0;
    Eigen::VectorXf p1;
    Eigen::VectorXf p2;
    // 1.0 for a hard-missing sample (BGEN ploidy-byte bit 7), else 0.0.
    // p0/p1/p2 are set to 0 for missing samples and must not be validated.
    Eigen::VectorXf missing;
    // E[G] = p1 + 2*p2, NaN for missing samples. Numerically identical to the
    // legacy readVariant() dosage for the same record.
    Eigen::VectorXf dosage;

    // Record-level facts, exposed so callers never have to re-parse.
    uint16_t n_alleles = 0;
    uint8_t ploidy = 0;        // validated ploidy (2 for supported records)
    uint8_t prob_bits = 0;     // BGEN probability bit depth (8 or 16)
    bool phased = false;      // always false for supported records
    bool missing_sample = false; // true if any sample is hard-missing

    int size() const { return static_cast<int>(p1.size()); }
    bool anyMissing() const { return missing_sample; }
    long long observedCount() const;
    // p = sum(E[G] over non-missing) / (2 * non-missing count).
    double alleleFrequency() const;
    // E[G] mean over non-missing samples; NaN if none observed.
    double meanDosage() const;
};

// Strict validation policy for the probability API.
struct BgenProbOptions {
    // Maximum allowed |p0+p1+p2-1| for a non-missing sample.
    // Layout 2 stores quantized probabilities, so this must be quantization
    // aware rather than a hard zero.
    double sum_tolerance = 1.0 / 255.0 + 1e-6;
    // Reject phased records (currently always true; phased is unsupported).
    bool reject_phased = true;
    // Reject any per-sample ploidy != required_ploidy.
    bool require_diploid = true;
    uint8_t required_ploidy = 2;
    // Reject n_alleles != 2 instead of silently skipping the record.
    bool reject_multiallelic = true;
    // Reject non-finite or out-of-range probabilities.
    bool validate_range = true;
};

// Structure to hold variant metadata from BGEN
struct BgenVariant {
    std::string id;
    std::string rsid;
    std::string chrom;
    uint32_t pos;
    uint16_t n_alleles;
    std::vector<std::string> alleles;
    uint32_t data_offset; // File offset to the start of genotype data block
    uint32_t data_length; // Length of compressed genotype data
};

class BgenReader : public GenotypeMatrix {
public:
    BgenReader();
    ~BgenReader();

    void open(const std::string& filename);
    void openSharedMetadata(const BgenReader& source);
    void close();

    // Read header and sample block
    // Returns number of samples
    uint32_t initialize();

    // Reset reader to the beginning of variant data
    void reset() const;

    // Iterate through variants
    // Returns false if end of file
    // Fills variant metadata and dosage vector (size N)
    // Dosage = P(AB) + 2*P(BB) (assuming Ref=AA, Alt=BB)
    // Actually standard is: Dosage of Alt allele.
    // If alleles are A1, A2. Dosage of A2 = 1*P(A1A2) + 2*P(A2A2).
    bool readVariant(std::string& chrom, std::string& rsid, uint32_t& pos, 
                     std::string& a1, std::string& a2, 
                     Eigen::VectorXf& dosage) const;

    // Strict, probability-preserving variant read.
    //
    // Unlike readVariant(), this NEVER silently skips or approximates:
    //   * multiallelic records throw instead of recursing to the next record;
    //   * phased records throw;
    //   * ploidy != 2 throws (per-sample ploidy bytes are validated);
    //   * non-finite or out-of-range probabilities throw;
    //   * |p0+p1+p2-1| > sum_tolerance throws.
    //
    // Variant/sample ordering is identical to readVariant(). Returns false only
    // at end of stream.
    bool readVariantProbs(std::string& chrom, std::string& rsid, uint32_t& pos,
                          std::string& a1, std::string& a2,
                          BgenGenotypeProbs& out,
                          const BgenProbOptions& options = BgenProbOptions()) const;
    
    // For random access (requires indexing, which we might skip for now)
    // We will stick to sequential reading for scanning.

    // Compute allele frequencies and missingness stats
    void computeStats(int n_threads = 1);

    // Matrix-vector multiplication: out = Z * v
    void multiply_Z_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;

    // Matrix-vector multiplication: out = Z' * v
    void multiply_Zt_v(const Eigen::VectorXd& v, Eigen::VectorXd& out, int n_threads = 1) override;

    // Seek to a specific variant by index (requires loadVariantInfo first)
    void seekVariant(uint32_t idx);

    // Get chrom range (start_idx, count) for a given chromosome label
    std::pair<uint32_t, uint32_t> getChromRange(const std::string& chrom) const;

    // Lazy-load variant info by scanning all variant headers
    void loadVariantInfo();

    const std::vector<BgenVariant>& getVariants() const { return variants_view ? *variants_view : variants; }
    uint32_t getNumSamples() const { return num_samples; }
    uint32_t getNumVariants() const { return num_variants; }
    uint32_t getCompressionType() const { return flags & 0x3u; }
    const std::vector<std::string>& getSampleIds() const { return sample_ids_view ? *sample_ids_view : sample_ids; }
    int rows() const override { return static_cast<int>(num_samples); }
    int cols() const override { return static_cast<int>(num_variants); }

// Start background prefetcher (optional)
    void startPrefetch();
    // Stop background prefetcher
    void stopPrefetcher();

public:
    // Export genotype data to various formats
    void exportToText(const std::string& geno_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToNpy(const std::string& npy_file, const std::string& var_info_file, const std::string& sample_info_file) const;
    void exportToPedMap(const std::string& ped_file, const std::string& map_file) const;
    void exportToVcf(const std::string& vcf_file) const;

private:
    std::string filename;
    mutable std::ifstream ifs;
    
    uint32_t num_samples = 0;
    uint32_t num_variants = 0;
    uint32_t header_length = 0;
    uint32_t flags = 0;
    uint32_t offset = 0; // Offset to first variant block
    size_t variant_data_offset = 0; // Actual offset after sample block
    
    std::vector<std::string> sample_ids;
    std::vector<BgenVariant> variants; // To store variant info if we parse them all
    const std::vector<std::string>* sample_ids_view = nullptr;
    
    // Chromosome range index (built during loadVariantInfo)
    std::vector<std::pair<uint32_t, uint32_t>> chrom_ranges_;
    std::vector<std::string> chrom_labels_;
    const std::vector<BgenVariant>* variants_view = nullptr;
    const std::vector<std::pair<uint32_t, uint32_t>>* chrom_ranges_view_ = nullptr;
    const std::vector<std::string>* chrom_labels_view_ = nullptr;
    
    // Stats
    bool stats_computed = false;
    std::vector<double> snp_means;
    std::vector<double> snp_inv_stds;
    
    // Internal buffers
    mutable std::vector<char> compressed_buffer;
    mutable std::vector<char> decompressed_buffer;
    mutable size_t current_variant_idx_ = 0;
    
    // Background Prefetcher
    bool prefetch_enabled = false;
    void* prefetcher_ptr = nullptr; // Opaque pointer to implementation
    
    // Helper to read simple types
    template<typename T> T read(std::ifstream& in) const;
    std::string readString(std::ifstream& in, uint16_t len_bytes = 2) const;
    
    // Decompress using zlib
    bool decompress(const char* src, size_t src_len, std::vector<char>& dst, size_t dst_len) const;
    
    // Parse Layout 2 genotype data (legacy dosage-only, behavior preserved)
    void parseLayout2(const char* data, size_t len, uint32_t n_samples, uint16_t n_alleles, Eigen::VectorXf& dosage) const;

    // Parse Layout 2 genotype data into full probability vectors.
    // Always strict: validates ploidy, phase, record length, probability range
    // and probability sums, and fails fast on any unsupported record.
    void parseLayout2Probs(const char* data, size_t len, uint32_t n_samples,
                           uint16_t n_alleles, BgenGenotypeProbs& out,
                           const BgenProbOptions& options,
                           const std::string& rsid) const;

    // Read variant identifying data + decompress the genotype block.
    // Returns false at EOF. `skip_multiallelic` reproduces the historical
    // recursive multiallelic skip; when false it throws instead. `strict`
    // enables record-framing validation; the legacy dosage path passes false so
    // its historical (unvalidated) behavior is preserved exactly. `d_out`
    // receives the record's uncompressed length D.
    bool readBlock(std::string& chrom, std::string& rsid, uint32_t& pos,
                   std::string& a1, std::string& a2, uint16_t& n_alleles,
                   bool skip_multiallelic, bool strict, uint32_t& d_out) const;
};

} // namespace cosmic

#endif // COSMIC_BGEN_READER_H
