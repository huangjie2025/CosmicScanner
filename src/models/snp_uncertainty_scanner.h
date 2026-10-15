#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "apy_inverse_operator.h"
#include "plink_reader.h"

namespace cosmic {

struct SnpUncertaintyScanInput {
    PlinkReader& reader;
    const ApyInverseOperator& core_inverse;
    const Eigen::MatrixXd& gw_cc;
    const Eigen::MatrixXd& pec_cc;
    const Eigen::MatrixXd* predictor_covariance = nullptr;
    const Eigen::VectorXd& effects;
    const std::vector<int>& snp_indices;
    const std::vector<int>& genotyped_keep;
    const std::vector<int>& core_genotyped_indices;
    double sigma2_a = 0.0;
    double sigma2_e = 0.0;
    double denominator = 0.0;
    int n_threads = 1;
    int block_size = 256;
};

struct SnpUncertaintyScanAudit {
    std::string mode;
    int block_size = 1;
    long long blocks = 0;
    long long variants_tested = 0;
    long long variants_decoded = 0;
    long long multiply_z_v_calls = 0;
    long long genotype_passes = 0;
    long long bed_bytes_read_estimate = 0;
    long long heap_allocations = 0;
    long long block_reallocations = 0;
    long long gcc_solve_calls = 0;
    long long gcc_rhs_columns = 0;
    long long gw_gemm_calls = 0;
    long long pec_gemm_calls = 0;
    long long missing_to_zero = 0;
    long long peak_transient_bytes = 0;
    double time_total = 0.0;
    double time_setup = 0.0;
    double time_onehot = 0.0;
    double time_genotype = 0.0;
    double time_core_gather = 0.0;
    double time_gcc_solve = 0.0;
    double time_gw_multiply = 0.0;
    double time_pec_multiply = 0.0;
    double time_quadratic = 0.0;
    double time_statistics = 0.0;
    double time_other = 0.0;
};

struct SnpUncertaintyScanResult {
    std::vector<double> variance;
    std::vector<double> se;
    std::vector<double> z;
    std::vector<double> p;
    std::vector<double> log10p;
    std::vector<std::string> status;
    int clamped = 0;
    int failures = 0;
    int zero_se = 0;
    SnpUncertaintyScanAudit audit;
};

class SnpUncertaintyScanner {
public:
    virtual ~SnpUncertaintyScanner() = default;
    virtual SnpUncertaintyScanResult scan(const SnpUncertaintyScanInput& input) const = 0;
};

class ScalarReferenceSnpUncertaintyScanner final : public SnpUncertaintyScanner {
public:
    SnpUncertaintyScanResult scan(const SnpUncertaintyScanInput& input) const override;
};

class BatchedOperatorSnpUncertaintyScanner final : public SnpUncertaintyScanner {
public:
    SnpUncertaintyScanResult scan(const SnpUncertaintyScanInput& input) const override;
};

// Exact selected covariance only. One inverse target must be supplied: a
// bounded full dense effect inverse, or the core-only Gcc inverse action.
// Never invoke genotype-backed full APY apply(MatrixXd) on SNP RHS columns.
struct ExactSnpUncertaintyScanInput {
    PlinkReader& reader;
    const Eigen::VectorXd& effects;
    const std::vector<int>& snp_indices;
    const std::vector<int>& sample_keep; // actual FAM/BED indices, in covariance order
    const Eigen::MatrixXd& evaluation_relationship;
    const Eigen::MatrixXd& predictor_covariance;
    const Eigen::MatrixXd* dense_effect_inverse = nullptr;
    const ApyInverseOperator* core_inverse = nullptr;
    double sigma2_a = 0.0;
    double denominator = 0.0;
    int block_size = 256;
};

class BatchedExactSnpUncertaintyScanner {
public:
    SnpUncertaintyScanResult scan(const ExactSnpUncertaintyScanInput& input) const;
};

} // namespace cosmic
