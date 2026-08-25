#pragma once
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "matrix_adapter.h"
#include "mme_builder.h"
#include "design.h"

namespace cosmic {

// ---------------------------------------------------------------------------
// Generic multi-random-effect LMM specification layer.
//
// All Gaussian "animal-model" GWAS models in CosmicScanner share the form:
//
//   y = Xb + x_j beta_j + sum_r Z_r u_r + e
//   u_r ~ N(0, K_r * sigma_r^2)
//   e  ~ N(0, I * sigma_e^2)
//
//   V = sum_r sigma_r^2 * Z_r K_r Z_r' + sigma_e^2 * I
//
// where each Z_r maps observations (records) to levels of random effect r,
// and K_r is the covariance kernel between levels:
//
//   direct additive genetic   : K = genomic relationship matrix (subject level)
//   permanent environment     : K = I over subjects      (Z is record->subject)
//   common environment/group  : K = I over groups        (Z is record->group)
//   maternal genetic          : K = GRM over dams        (optional cross-cov)
//   RRM coefficients          : Z rows scaled by Legendre/B-spline loadings,
//                               K = D (genetic) or P (permanent env)
//
// This header is intentionally a thin, testable specification + builder layer.
// It deliberately does NOT override the optimized solver paths of
// RepeatabilityModel / RrmModel / MultivariateModel: those keep their
// Woodbury / eigen / STCG / Kronecker implementations.
// ---------------------------------------------------------------------------

/// Covariance kernel between levels of a random effect.
enum class RandomEffectKernel {
    /// K_r = genomic relationship matrix built between subject levels.
    GRM,
    /// K_r = I (levels are exchangeable a priori).
    Identity,
    /// Custom dense kernel supplied by the caller (e.g. dam GRM, pedigree A).
    Custom
};

/// Specification of one random effect.
struct RandomEffectSpec {
    std::string name;           ///< e.g. "additive", "pe", "common", "maternal"

    /// Covariance kernel between levels.
    RandomEffectKernel kernel = RandomEffectKernel::Identity;

    /// Custom dense kernel (used when kernel == Custom). Row/col = levels.
    Eigen::MatrixXd custom_kernel;

    /// Observation -> level mapping. Must have one entry per record.
    /// Levels are dense-indexed 0..n_levels-1.
    std::vector<int> obs_to_level;

    /// Number of levels (n_levels). Must be > max(obs_to_level).
    int n_levels = 0;

    /// Optional per-record covariate loadings (RRM-style vector effects).
    /// If non-empty, record r contributes loading[r] to the level effect.
    std::vector<double> obs_loadings;

    /// Index of the variance component in the sigma2 vector (used when a
    /// single spec set is expanded to multiple sigma2 parameters, e.g.
    /// direct-maternal covariance pairing).
    int vc_index = -1;
};

/// Dense record-level null model for a FIXED covariance V.
///
/// V = sum_r sigma2_r * Z_r K_r Z_r' + sigma2_e * I  (records x records)
///
/// The factorization is done once; the P-projection pieces needed by the
/// per-SNP fixed-null scan (score / Wald / LRT) are cached:
///
///   P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1
///   score_j = x_j' P y
///   Var(score_j) = x_j' P x_j
class RecordLevelNullModel {
public:
    /// Build V and factorize it.
    ///   X_records : n_records x p fixed-effect design (double)
    ///   G_subject : subject-level GRM over the SAME level ordering used by
    ///               the "GRM" kernel specs (may be null if no GRM kernel).
    bool prepare(const Eigen::VectorXd& y,
                 const Eigen::MatrixXd& X_records,
                 const std::vector<RandomEffectSpec>& specs,
                 const std::vector<double>& sigma2,
                 double sigma2_e,
                 const Eigen::MatrixXd* G_subject);

    /// Aliased: build record covariance only (no factorization).
    static Eigen::MatrixXd buildRecordCovariance(
        const std::vector<RandomEffectSpec>& specs,
        const std::vector<double>& sigma2,
        double sigma2_e,
        const Eigen::MatrixXd* G_subject,
        int n_records);

    bool ok() const { return factorized_; }

    int n_records() const { return n_records_; }
    int n_covars() const { return n_covars_; }

    /// V^-1 y  (V fixed at null estimates)
    const Eigen::VectorXd& VinvY() const { return Vinv_y_; }
    /// V^-1 X
    const Eigen::MatrixXd& VinvX() const { return VinvX_; }
    /// (X' V^-1 X)^-1
    const Eigen::MatrixXd& XtVinvXInv() const { return XtVinvX_inv_; }
    /// P y = e / sigma2_e style projection residual
    const Eigen::VectorXd& Py() const { return Py_; }

    /// y' P y (null quadratic form)
    double nullQuad() const { return null_quad_; }
    double nEff() const { return n_eff_; }

    double sigma2E() const { return sigma2_e_; }
    const std::vector<double>& sigma2() const { return sigma2_; }
    /// Condition diagnostics: min eigenvalue of V and its jitter factor used.
    double minEig() const { return min_eig_; }
    double logDetV() const { return log_det_v_; }

    /// Apply V^-1 to a block of record-level vectors (e.g. genotype blocks).
    Eigen::MatrixXd applyVinv(const Eigen::MatrixXd& block) const;

    /// Per-SNP fixed-null statistics. x must be record-level and already
    /// mean-centered / missing-imputed by the caller.
    /// Returns beta, se, wald_stat, p (1-df chi2). out_var_score optional.
    bool testSnp(const Eigen::VectorXd& x, double& beta, double& se,
                 double& wald_stat, double& p,
                 double* out_var_score = nullptr, double* out_score = nullptr) const;

private:
    int n_records_ = 0;
    int n_covars_ = 0;
    Eigen::MatrixXd V_;
    Eigen::LDLT<Eigen::MatrixXd> V_ldlt_;   // V = LDL'
    bool factorized_ = false;
    double min_eig_ = 0.0;
    double log_det_v_ = 0.0;
    double n_eff_ = 1.0;
    double null_quad_ = 0.0;
    double sigma2_e_ = 1.0;
    std::vector<double> sigma2_;
    Eigen::VectorXd Vinv_y_;
    Eigen::MatrixXd VinvX_;
    Eigen::MatrixXd XtVinvX_inv_;
    Eigen::VectorXd Py_;
    Eigen::LDLT<Eigen::MatrixXd> XtVinvX_llt_; // reuse for per-SNP solves
};

/// Convert specs to libCosmicCore RandomComponent list (MME / AI_REML reuse).
/// Ginv: subject-level inverse GRM (required iff any spec uses GRM kernel).
/// group_I: sparse identity over group levels (built internally for Identity
///          kernels). The returned components are valid as long as Ginv and
///          the copied group identities stay alive.
std::vector<RandomComponent> toRandomComponents(
    const std::vector<RandomEffectSpec>& specs,
    const Eigen::MatrixXd* Ginv,
    std::vector<Eigen::SparseMatrix<double>>& scratch_group_I,
    std::vector<std::unique_ptr<AbstractMatrix>>& scratch_group_adapters);

} // namespace cosmic