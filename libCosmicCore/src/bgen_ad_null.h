#pragma once

#include <map>
#include <limits>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include "bgen_reader.h"
#include "eigen_vce_ad.h"

namespace cosmic {

// E3.1 independent probability AD null backend. No association scan or CLI
// dispatch uses this API. Existing GRM/dominance-only paths are unaffected.
struct BgenAdKernelOptions {
    double min_maf = 0.01;
    double max_missing = 1.0;
    double min_info_a = 0.0;
    double min_info_d = 0.0;
    size_t block_size = 256;
    std::set<std::string> include_chromosomes;
    std::set<std::string> exclude_chromosomes;
    std::unordered_set<std::string> extract;
    std::unordered_set<std::string> exclude;
};

struct BgenAdRawPartition {
    Eigen::MatrixXd additive;
    Eigen::MatrixXd dominance;
    double denominator_a = 0.0;
    double denominator_d = 0.0;
    long long markers = 0;
};

struct BgenAdKernelPair {
    Eigen::MatrixXd additive;
    Eigen::MatrixXd dominance;
    double denominator_a = 0.0;
    double denominator_d = 0.0;
    long long markers = 0;
};

struct BgenAdKernelResult {
    BgenAdRawPartition global_raw;
    std::map<std::string, BgenAdRawPartition> chromosome_raw;
    std::vector<int> sample_indices;
    std::vector<std::string> sample_ids;
    std::vector<std::string> accepted_variant_ids;
    long long candidates = 0;
    long long rejected_missing = 0;
    long long rejected_maf = 0;
    long long rejected_info = 0;

    BgenAdKernelPair global() const;
    // Rejects unknown chromosomes and empty/degenerate outside partitions;
    // never substitutes the global pair for an unavailable LOCO pair.
    BgenAdKernelPair loco(const std::string& chromosome) const;
};

struct BgenAdRemlCandidate {
    std::string model;
    std::string start;
    std::string boundary_status;
    double sigma_a2 = 0.0;
    double sigma_d2 = 0.0;
    double sigma_e2 = 0.0;
    double log_likelihood = -std::numeric_limits<double>::infinity();
    double likelihood_delta = std::numeric_limits<double>::infinity();
    // Dimensionless average gradient in numerically scaled variance coordinates.
    // Free-variable projected gradient and full-domain KKT residual are separate.
    double gradient_norm = std::numeric_limits<double>::infinity();
    double kkt_residual = std::numeric_limits<double>::infinity();
    double raw_gradient_norm = std::numeric_limits<double>::infinity();
    bool face_converged = false;
    bool kkt_valid = false;
    int iterations = 0;
};

struct BgenAdNullResult {
    double sigma_a2 = 0.0;
    double sigma_d2 = 0.0;
    double sigma_e2 = 0.0;
    double log_likelihood = 0.0;
    bool converged = false;
    int iterations = 0;
    double likelihood_delta = std::numeric_limits<double>::infinity();
    double gradient_norm = std::numeric_limits<double>::infinity();
    double kkt_residual = std::numeric_limits<double>::infinity();
    double raw_gradient_norm = std::numeric_limits<double>::infinity();
    std::string boundary_status;
    std::string selected_start;
    double gradient_tolerance = 0.0;
    double likelihood_tolerance = 0.0;
    int information_rank = 0;
    double sigma_d2_local_se = std::numeric_limits<double>::infinity();
    std::vector<BgenAdRemlCandidate> candidates;
};

class BgenAdNullBackend {
public:
    // Empty sample_indices means all samples in reader order. A nonempty
    // selection is validated for uniqueness and bounds and retains its order.
    static BgenAdKernelResult build(BgenReader& reader,
                                   const BgenAdKernelOptions& options = {},
                                   const std::vector<int>& sample_indices = {});
    // BGEN-only boundary/multi-start refinement of the unchanged contrast-REML
    // likelihood. Legacy EigenVceAD and BED/PGEN dispatch remain untouched.
    // y and covariates must already be in the kernel's selected sample order.
    static BgenAdNullResult fit(const BgenAdKernelPair& kernels,
                               const Eigen::VectorXd& y,
                               const Eigen::MatrixXf& covariates,
                               const EigenVceAD::Options& options);
};

} // namespace cosmic
