#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"

namespace cosmic::association::cojo {

enum class CoJoMode {
    CONDITIONAL,
    STEPWISE,
    JOINT
};

struct CoJoOptions {
    CoJoMode mode = CoJoMode::CONDITIONAL;
    std::string ld_reference;
    std::string summary_stats;
    std::string condition_snps_file;
    double threshold = 5e-8;
    double collinearity = 0.9;
    int window_size = 10000000;
    int max_iter = 20;
    bool use_reference_ld = true;
    double maf_filter = 0.01;
    bool standardize = true;
    std::string region_chrom;
    long long region_start = -1;
    long long region_end = -1;
    bool canonical_summary_only = false;
};

struct SummaryRecord {
    std::string snp_id;
    std::string chrom;
    long long position = 0;
    std::string a1;
    std::string a2;
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    double freq = 0.0;
    int n = 0;
    double z_score = 0.0;
    double weight = 0.0;
};

struct LDReference {
    std::vector<std::string> snp_ids;
    std::vector<std::string> chroms;
    std::vector<long long> positions;
    std::vector<std::string> a1s;
    std::vector<std::string> a2s;
    Eigen::MatrixXd ld_matrix;
    std::unordered_map<std::string, int> snp_to_idx;
};

struct CoJoResult {
    std::string snp_id;
    std::string chrom;
    long long position = 0;
    std::string a1;
    std::string a2;
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    double freq = 0.0;
    int n = 0;
    int selected_order = 0;
    double conditional_beta = 0.0;
    double conditional_se = 0.0;
    double conditional_p = 1.0;
    bool selected = false;
    bool valid = false;
    std::string warning;
};

class CoJoEngine {
public:
    explicit CoJoEngine(const CoJoOptions& opts);

    void loadSummaryStats(const std::string& filepath);
    void loadLDReference(const std::string& filepath);
    void loadLDReferenceFromGenotype(const std::string& summary_stats_file,
                                     PlinkReader* plink_reader,
                                     BgenReader* bgen_reader,
                                     PgenReader* pgen_reader,
                                     VcfReader* vcf_reader,
                                     bool is_bgen, bool is_pgen, bool is_vcf);
    void loadConditionSnps(const std::string& filepath);
    std::vector<CoJoResult> run();

private:
    CoJoOptions options_;
    std::vector<SummaryRecord> summary_stats_;
    LDReference ld_ref_;
    std::vector<std::string> condition_snps_;
    std::unordered_map<std::string, int> summary_idx_;

    void alignAlleles();
    void validateLDReference() const;
    int ldIndexForSummary(int summary_index) const;
    std::vector<int> conditionIndices() const;
    Eigen::MatrixXd extractLDMatrix(const std::vector<int>& indices) const;
    Eigen::VectorXd computeConditionalEffects(const std::vector<int>& indices) const;
    Eigen::MatrixXd computeConditionalCovariance(const std::vector<int>& indices) const;
    std::vector<CoJoResult> runConditional() const;
    std::vector<CoJoResult> runStepwise() const;
    std::vector<CoJoResult> runJoint() const;
    std::vector<int> selectIndependentSignals() const;
    bool checkCollinearity(int idx1, int idx2) const;

public:
    // For testing
    void setLDReferenceForTesting(const LDReference& ref) { ld_ref_ = ref; }
    void setSummaryStatsForTesting(const std::vector<SummaryRecord>& stats) {
        summary_stats_ = stats;
        summary_idx_.clear();
        for (size_t i = 0; i < summary_stats_.size(); ++i)
            summary_idx_[summary_stats_[i].snp_id] = static_cast<int>(i);
    }
    void setConditionSnpsForTesting(const std::vector<std::string>& snps) { condition_snps_ = snps; }
    void alignAllelesForTesting() { alignAlleles(); }
    const std::vector<SummaryRecord>& getSummaryStatsForTesting() const { return summary_stats_; }
    const std::vector<SummaryRecord>& getSummaryStats() const { return summary_stats_; }
    const LDReference& getLDReference() const { return ld_ref_; }
};

void writeCoJoResults(
    const std::vector<CoJoResult>& results,
    const std::string& out_prefix,
    bool write_jma = true);

} // namespace cosmic::association::cojo
