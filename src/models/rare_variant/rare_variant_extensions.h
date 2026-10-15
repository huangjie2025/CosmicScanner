#ifndef COSMIC_SCANNER_RARE_VARIANT_EXTENSIONS_H
#define COSMIC_SCANNER_RARE_VARIANT_EXTENSIONS_H

#include <Eigen/Dense>
#include <vector>

namespace cosmic {
namespace rare_variant {

inline constexpr const char* kSkatOGridMethod = "SKAT_O_EXPERIMENTAL_CONSERVATIVE_GRID";
inline constexpr const char* kSkatOCanonicalMethod = "SKAT_O_CANONICAL";
inline constexpr const char* kAcatVMethod = "ACAT_V";

struct MixtureTailResult {
    double p_value = 1.0;
    bool conservative_bound = false;
};

// Upper tail of sum(lambda_i * chi-square_1), lambda_i >= 0.
// Rank-one/equal spectra use the exact chi-square tail. Unequal spectra use a
// positive gamma-mixture series with an explicit remaining-mass upper bound,
// tightened by a Chernoff bound. No moment matching or unflagged approximation.
// Empty/all-zero spectra are a point mass at zero. Invalid input throws.
MixtureTailResult mixtureChiSquareUpperTail(double q, const Eigen::VectorXd& lambda);

struct SkatOGridResult {
    double p_value = 1.0;       // min(1, number of distinct rhos * min_p_value)
    double min_p_value = 1.0;   // Diagnostic ONLY, never the omnibus p-value
    double selected_rho = 0.0;
    double statistic = 0.0;     // Selected Q, with weights normalized to max(w)=1
    int n_variants = 0;        // Positive weight and positive score variance
    std::vector<double> rho_grid;
    std::vector<double> rho_statistics;
    std::vector<double> rho_p_values;
};

struct SkatOResult : SkatOGridResult {
    const char* p_method = "DEPENDENCE_ADJUSTED_GAUSSIAN_QMC";
    bool fallback_used = false;
    int null_draws = 0;
};

// Canonical fixed-grid SKAT-O statistic with a dependence-adjusted omnibus tail.
// The shared Gaussian score null is integrated by deterministic Halton QMC.
SkatOResult skatOCanonical(
    const Eigen::VectorXd& scores, const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights,
    const std::vector<double>& rho_grid = {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1.0},
    int null_draws = 65536);

// U ~ N(0, covariance); weights are nonnegative score multipliers, NOT squared
// kernel weights. Q_rho=(1-rho)*sum((w_i U_i)^2)+rho*(sum(w_i U_i))^2.
// The grid MUST be fixed independently of U. Duplicates are removed.
// This is an experimental conservative grid omnibus, NOT canonical SKAT-O.
// The Gaussian-score assumption is asymptotic for fitted LMM/GLM null models.
SkatOGridResult skatOConservativeGrid(
    const Eigen::VectorXd& scores, const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights,
    const std::vector<double>& rho_grid = {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1.0});

struct AcatVResult {
    double p_value = 1.0;
    long double statistic = 0.0L; // Weighted Cauchy statistic; may be infinite
    int n_variants = 0;
};

// Direct ACAT-V combination with nonnegative, pre-specified combination weights.
// Weights are normalized internally, NOT squared. Empty inputs return p=1.
// All-zero weights, nonfinite inputs, or p outside [0,1] throw invalid_argument.
// Active p=0 gives p=0; active p=1 gives p=1. Both endpoints together are
// undefined and throw. Zero-weight endpoints are ignored (but validated).
// Under arbitrary LD, the Cauchy reference is a tail approximation, not an
// exact finite-sample calibration. This does not pool ultra-rare variants.
AcatVResult acatV(const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights);

// Convenience adapter: p_i=erfc(|U_i|/sqrt(2*covariance_ii)). Zero-variance,
// zero-score variants are omitted. Weights here are ACAT combination weights.
// Both score APIs validate dimensions, finiteness, symmetry, PSD and nullspace
// consistency. Empty inputs return p=1; all-zero weights on nonempty input throw.
AcatVResult acatVFromScores(
    const Eigen::VectorXd& scores, const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights);

} // namespace rare_variant
} // namespace cosmic

#endif
