#pragma once

#include "ad_engine.h" // Existing internal AdResult; no production schema changes.
#include "bgen_ad_null.h"
#include <optional>

namespace cosmic {

struct BgenAdPreparedMarker {
    Eigen::VectorXd additive;
    Eigen::VectorXd dominance;
    double frequency = 0.0;
    double missing_rate = 1.0;
    double info_a = 0.0;
    double info_d = 0.0;
    long long observed = 0;
};

// E3.2a independent internal engine. It is not wired to GwasScanner or CLI.
class BgenAdScanEngine {
public:
    struct Options {
        double min_maf = 0.01;
        double max_missing = 1.0;
        double min_info_a = 0.0;
        double min_info_d = 0.0;
        size_t block_size = 512;
        std::unordered_set<std::string> extract;
        std::unordered_set<std::string> exclude;
        std::string chromosome;
        int64_t from_bp = -1;
        int64_t to_bp = -1;
    };

    struct Statistics {
        long long tested = 0;
        long long qc_filtered = 0;
        long long target_filtered = 0;
        long long rank_deficient = 0;
    };

    BgenAdScanEngine();
    explicit BgenAdScanEngine(const Options& options);

    static BgenAdPreparedMarker encode(const BgenGenotypeProbs& probabilities,
                                       const std::vector<int>& sample_indices);

    void fitNull(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                 const Eigen::MatrixXf& covariates, const EigenVceAD::Options& options,
                 const std::string& excluded_chromosome = "");
    // Orchestration may fit once itself to persist a failure certificate, then
    // prepare this same selected null without a second fit.
    void prepareCertified(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                          const Eigen::MatrixXf& covariates, const BgenAdNullResult& fit,
                          const std::string& excluded_chromosome = "");

    // Explicit known/fixed-covariance path for internal reference validation.
    // The robust fitted path also rebuilds exactly this workspace using its
    // selected components, never the legacy warm-start components.
    void prepareFixed(const BgenAdKernelResult& kernels, const Eigen::VectorXd& y,
                      const Eigen::MatrixXf& covariates, double sigma_a2,
                      double sigma_d2, double sigma_e2,
                      const std::string& excluded_chromosome = "");
    void updateFixedPhenotype(const Eigen::VectorXd& y);

    void scan(BgenReader& reader, std::vector<AdResult>& results);
    const Statistics& statistics() const { return statistics_; }
    const std::optional<BgenAdNullResult>& nullResult() const { return null_result_; }
    long long nullFitCount() const { return null_fit_count_; }
    const Eigen::MatrixXd& rotation() const { return rotation_; }
    const Eigen::VectorXd& varianceEigenvalues() const { return variance_eigenvalues_; }

private:
    Options options_;
    Statistics statistics_;
    bool prepared_ = false;
    std::vector<int> sample_indices_;
    std::vector<std::string> sample_ids_;
    Eigen::MatrixXd rotation_;
    Eigen::VectorXd variance_eigenvalues_;
    Eigen::VectorXd inverse_variance_;
    Eigen::VectorXd weighted_y_;
    std::optional<BgenAdNullResult> null_result_;
    long long null_fit_count_ = 0;

    void evaluateBlock(const Eigen::MatrixXd& additive, const Eigen::MatrixXd& dominance,
                       const std::vector<int>& indices, std::vector<AdResult>& results);
};

} // namespace cosmic
