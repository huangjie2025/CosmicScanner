#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include <fstream>

namespace cosmic::association::sumstats {

enum class SumstatsFormat {
    COSMIC,       // .cosmic.sumstats
    SSF,          // EBI GWAS-SSF
    METAL,        // METAL format
    PLINK         // PLINK .assoc.linear/.assoc.logistic
};

enum class AlleleMatchStatus {
    MATCH,           // Identical alleles
    SWAPPED,         // A1/A2 swapped, beta negated
    STRAND_FLIP,     // Strand complement (A<->T, C<->G)
    STRAND_FLIP_SWAPPED, // Strand complement with A1/A2 swapped
    AMBIGUOUS,       // A/T or C/G - cannot resolve strand
    MISMATCH,        // Alleles don't match
    MISSING_REF      // Not in reference
};

struct SumstatsRecord {
    std::string chrom;
    std::string snp_id;
    long long pos = 0;
    std::string a1;      // Effect allele (A1)
    std::string a2;      // Non-effect allele (A2)
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    long long n = 0;
    double freq = 0.0;   // Frequency of A1
    double log10_p = 0.0;
    std::string model;   // Model that generated this (LMM, LM, etc.)
    
    // Optional metadata
    double z_score = 0.0;
    double weight = 0.0;     // 1/SE^2
    double info = 1.0;       // Imputation info score
    std::string warning;     // Any warnings
};

struct SumstatsMetadata {
    std::string created_by = "CosmicScanner";
    std::string version = "1.0";
    std::string phenotype = "";
    std::string model = "";
    std::string reference_genome = "GRCh38";
    int n_samples = 0;
    int n_variants = 0;
    std::string date;
    std::string command_line;
};

struct AlleleHarmonizationResult {
    SumstatsRecord record;
    AlleleMatchStatus status = AlleleMatchStatus::MISMATCH;
    std::string reason;
    bool harmonized = false;
    double original_beta = 0.0;
    std::string original_a1;
    std::string original_a2;
};

class SumstatsWriter {
public:
    SumstatsWriter(const std::string& filepath, const SumstatsMetadata& meta = {});
    ~SumstatsWriter();
    
    void writeHeader();
    void writeRecord(const SumstatsRecord& record);
    void close();
    
private:
    std::ofstream file_;
    SumstatsMetadata meta_;
    bool header_written_ = false;
};

class SumstatsReader {
public:
    SumstatsReader(const std::string& filepath);
    ~SumstatsReader();
    
    bool readHeader();
    bool readRecord(SumstatsRecord& record);
    const SumstatsMetadata& getMetadata() const { return meta_; }
    bool hasColumn(const std::string& col) const;
    
private:
    std::ifstream file_;
    SumstatsMetadata meta_;
    std::vector<std::string> column_names_;
    std::unordered_map<std::string, int> col_idx_;
    bool header_read_ = false;
};

// Allele harmonization functions
std::string complementAllele(const std::string& allele);
bool isAmbiguousSNP(const std::string& a1, const std::string& a2);
AlleleMatchStatus checkAlleleMatch(const std::string& a1, const std::string& a2,
                                   const std::string& ref_a1, const std::string& ref_a2);

AlleleHarmonizationResult harmonizeAlleles(const SumstatsRecord& record,
                                           const std::string& ref_a1,
                                           const std::string& ref_a2);

std::vector<AlleleHarmonizationResult> harmonizeSumstats(
    const std::vector<SumstatsRecord>& records,
    const std::unordered_map<std::string, std::pair<std::string, std::string>>& ref_alleles);

// Convenience: write sumstats from GWAS results
void writeCosmicSumstats(const std::string& filepath,
                         const std::vector<SumstatsRecord>& records,
                         const SumstatsMetadata& meta = {});

} // namespace cosmic::association::sumstats
