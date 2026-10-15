#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

namespace cosmic::association::meta {

enum class MetaAnalysisMethod {
    FIXED_EFFECT,
    RANDOM_EFFECT_DL,
    RANDOM_EFFECT_REML
};

struct StudySummary {
    std::string study_id;
    std::string snp_id;
    std::string chrom;
    long long position = 0;
    std::string a1;
    std::string a2;
    double beta = 0.0;
    double se = 0.0;
    double p_value = 1.0;
    int n = 0;
    double freq = 0.0;
    double z_score = 0.0;
    double weight = 0.0;
};

struct MetaAnalysisOptions {
    MetaAnalysisMethod method = MetaAnalysisMethod::FIXED_EFFECT;
    double heterogeneity_threshold = 0.05;
    int min_studies = 2;
    bool use_gc_correction = false;
    double gc_lambda = 1.0;
    bool allow_sample_overlap = false;
    double sample_overlap_correction = 0.0;
};

struct HeterogeneityStats {
    double q_statistic = 0.0;
    double q_df = 0.0;
    double q_p_value = 1.0;
    double i2 = 0.0;
    double tau2 = 0.0;
    double h2 = 0.0;
    bool significant_heterogeneity = false;
};

struct MetaResult {
    std::string snp_id;
    std::string chrom;
    long long position = 0;
    std::string a1;
    std::string a2;
    double beta = 0.0;
    double se = 0.0;
    double z_score = 0.0;
    double p_value = 1.0;
    int n_studies = 0;
    int total_n = 0;
    double freq = 0.0;
    HeterogeneityStats heterogeneity;
    bool valid = false;
    std::string warning;
};

struct MetaAnalysisConfig {
    MetaAnalysisOptions options;
    std::vector<std::string> study_files;
    std::string output_prefix;
    std::vector<std::string> extract_snps;
    std::vector<std::string> exclude_snps;
    std::string chrom_filter;
    long long from_bp = -1;
    long long to_bp = -1;
};

class MetaAnalysisEngine {
public:
    explicit MetaAnalysisEngine(const MetaAnalysisOptions& opts);

    void addStudy(const std::string& study_id, const std::vector<StudySummary>& summaries);
    void addStudyFromFile(const std::string& study_id, const std::string& filepath);
    std::vector<MetaResult> run();

    const std::unordered_map<std::string, std::vector<StudySummary>>& getStudies() const;

private:
    MetaAnalysisOptions options_;
    std::unordered_map<std::string, std::vector<StudySummary>> studies_;
    std::unordered_map<std::string, int> study_indices_;
    std::vector<std::string> study_order_;

    MetaResult metaAnalyzeSnp(const std::vector<StudySummary>& snp_data) const;
    HeterogeneityStats computeHeterogeneity(
        const std::vector<double>& betas,
        const std::vector<double>& weights,
        const std::vector<int>& ns) const;
    double estimateTau2DL(
        const std::vector<double>& betas,
        const std::vector<double>& weights,
        const HeterogeneityStats& het) const;
    double estimateTau2REML(
        const std::vector<double>& betas,
        const std::vector<double>& weights,
        double initial_tau2) const;
    double computeQStatistic(
        const std::vector<double>& betas,
        const std::vector<double>& weights) const;
};

struct MetaWriterOptions {
    bool write_metal = true;
    bool write_gwama = false;
    bool write_simple = false;
};

void writeMetaResults(
    const std::vector<MetaResult>& results,
    const std::string& out_prefix,
    const MetaWriterOptions& opts = {});

} // namespace cosmic::association::meta
