#include "models/multi_random_effect_lmm.h"
#include "matrix_io.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmic {

// ---------------------------------------------------------------------------
// RecordLevelNullModel
// ---------------------------------------------------------------------------

Eigen::MatrixXd RecordLevelNullModel::buildRecordCovariance(
    const std::vector<RandomEffectSpec>& specs,
    const std::vector<double>& sigma2,
    double sigma2_e,
    const Eigen::MatrixXd* G_subject,
    int n_records) {
    if (n_records <= 0) return Eigen::MatrixXd();
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n_records, n_records);

    for (const auto& spec : specs) {
        if ((int)spec.obs_to_level.size() != n_records) {
            throw std::runtime_error(
                "RandomEffectSpec '" + spec.name +
                "' obs_to_level size does not match n_records");
        }
        int vc = (spec.vc_index >= 0) ? spec.vc_index
                                      : ((int)specs.size() > 1 ? (&spec - &specs[0]) : 0);
        if (vc < 0 || vc >= (int)sigma2.size()) {
            throw std::runtime_error("RandomEffectSpec '" + spec.name +
                                     "' references out-of-range variance component");
        }
        const double s2 = sigma2[vc];
        if (!(s2 >= 0.0) || !std::isfinite(s2)) continue;

        Eigen::MatrixXd K;
        switch (spec.kernel) {
            case RandomEffectKernel::GRM: {
                if (!G_subject) {
                    throw std::runtime_error(
                        "RandomEffectSpec '" + spec.name +
                        "' uses GRM kernel but no subject GRM was provided");
                }
                if (G_subject->rows() < spec.n_levels) {
                    throw std::runtime_error(
                        "RandomEffectSpec '" + spec.name +
                        "' GRM kernel has fewer levels than the GRM rows");
                }
                K = G_subject->topLeftCorner(spec.n_levels, spec.n_levels);
                break;
            }
            case RandomEffectKernel::Identity:
                K = Eigen::MatrixXd::Identity(spec.n_levels, spec.n_levels);
                break;
            case RandomEffectKernel::Custom:
                if (spec.custom_kernel.rows() < spec.n_levels) {
                    throw std::runtime_error(
                        "RandomEffectSpec '" + spec.name +
                        "' custom kernel is smaller than n_levels");
                }
                K = spec.custom_kernel.topLeftCorner(spec.n_levels, spec.n_levels);
                break;
        }

        // V += s2 * Z_r K_r Z_r'  (with optional per-record loadings)
        if (spec.obs_loadings.empty()) {
            for (int i = 0; i < n_records; ++i) {
                const int li = spec.obs_to_level[i];
                if (li < 0 || li >= spec.n_levels) continue;
                for (int j = 0; j < n_records; ++j) {
                    const int lj = spec.obs_to_level[j];
                    if (lj < 0 || lj >= spec.n_levels) continue;
                    V(i, j) += s2 * K(li, lj);
                }
            }
        } else {
            if ((int)spec.obs_loadings.size() != n_records) {
                throw std::runtime_error(
                    "RandomEffectSpec '" + spec.name +
                    "' obs_loadings size does not match n_records");
            }
            for (int i = 0; i < n_records; ++i) {
                const int li = spec.obs_to_level[i];
                if (li < 0 || li >= spec.n_levels) continue;
                const double wi = spec.obs_loadings[i];
                for (int j = 0; j < n_records; ++j) {
                    const int lj = spec.obs_to_level[j];
                    if (lj < 0 || lj >= spec.n_levels) continue;
                    V(i, j) += s2 * wi * K(li, lj) * spec.obs_loadings[j];
                }
            }
        }
    }

    if (sigma2_e > 0.0 && std::isfinite(sigma2_e)) {
        V.diagonal().array() += sigma2_e;
    }
    return 0.5 * (V + V.transpose()).eval();
}

bool RecordLevelNullModel::prepare(const Eigen::VectorXd& y,
                                   const Eigen::MatrixXd& X_records,
                                   const std::vector<RandomEffectSpec>& specs,
                                   const std::vector<double>& sigma2,
                                   double sigma2_e,
                                   const Eigen::MatrixXd* G_subject) {
    factorized_ = false;
    n_records_ = (int)y.size();
    n_covars_ = (int)X_records.cols();
    if (X_records.rows() != n_records_) {
        throw std::runtime_error("RecordLevelNullModel::prepare: X row mismatch");
    }
    sigma2_ = sigma2;
    sigma2_e_ = sigma2_e;

    V_ = buildRecordCovariance(specs, sigma2, sigma2_e, G_subject, n_records_);

    if (!factorize_with_jitter(V_, V_ldlt_)) {
        std::cerr << "RecordLevelNullModel: factorization of V failed\n";
        return false;
    }
    min_eig_ = V_ldlt_.vectorD().minCoeff();  // D of LDLT (after jitter)
    log_det_v_ = (V_ldlt_.vectorD().array().abs().cwiseMax(1e-300).log()).sum();

    Vinv_y_ = V_ldlt_.solve(y);
    VinvX_ = V_ldlt_.solve(X_records);
    Eigen::MatrixXd XtVinvX = X_records.transpose() * VinvX_;
    if (!factorize_with_jitter(XtVinvX, XtVinvX_llt_)) {
        std::cerr << "RecordLevelNullModel: factorization of X'V^(-1)X failed\n";
        return false;
    }
    XtVinvX_inv_ = XtVinvX_llt_.solve(Eigen::MatrixXd::Identity(n_covars_, n_covars_));

    Eigen::VectorXd bx = XtVinvX_llt_.solve(X_records.transpose() * Vinv_y_);
    Py_ = Vinv_y_ - VinvX_ * bx;

    null_quad_ = y.dot(Py_);
    n_eff_ = std::max(1.0, (double)n_records_ - (double)n_covars_);
    factorized_ = true;
    return true;
}

Eigen::MatrixXd RecordLevelNullModel::applyVinv(const Eigen::MatrixXd& block) const {
    if (!factorized_) {
        throw std::runtime_error("RecordLevelNullModel: not prepared");
    }
    return V_ldlt_.solve(block);
}

bool RecordLevelNullModel::testSnp(const Eigen::VectorXd& x, double& beta, double& se,
                                   double& wald_stat, double& p,
                                   double* out_var_score, double* out_score) const {
    if (!factorized_) return false;
    const double score = x.dot(Py_);

    Eigen::VectorXd Vinvx = V_ldlt_.solve(x);
    // x' V^-1 X -> p-vector
    Eigen::VectorXd xvX = VinvX_.transpose() * x;
    // x' P x = x'V^-1 x - x'V^-1 X (X'V^-1 X)^-1 X' V^-1 x
    Eigen::VectorXd proj = XtVinvX_inv_ * xvX;
    const double term2 = xvX.dot(proj);
    const double var_score = x.dot(Vinvx) - term2;
    if (out_score) *out_score = score;
    if (out_var_score) *out_var_score = var_score;
    if (!(var_score > 1e-12) || !std::isfinite(score) || !std::isfinite(var_score)) return false;

    beta = score / var_score;
    se = std::sqrt(1.0 / var_score);
    if (!(se > 0.0)) return false;
    wald_stat = beta * beta / (se * se);
    // 1-df chi2 upper tail
    p = std::erfc(std::sqrt(wald_stat) / std::sqrt(2.0));
    if (!std::isfinite(p)) p = 0.0;
    return true;
}

// ---------------------------------------------------------------------------
// toRandomComponents (MME / AI_REML reuse)
// ---------------------------------------------------------------------------

std::vector<RandomComponent> toRandomComponents(
    const std::vector<RandomEffectSpec>& specs,
    const Eigen::MatrixXd* Ginv,
    std::vector<Eigen::SparseMatrix<double>>& scratch_group_I,
    std::vector<std::unique_ptr<AbstractMatrix>>& scratch_group_adapters) {
    std::vector<RandomComponent> out;
    out.reserve(specs.size());
    scratch_group_I.clear();
    scratch_group_adapters.clear();

    for (size_t si = 0; si < specs.size(); ++si) {
        const auto& spec = specs[si];
        RandomComponent comp;
        comp.id_map = spec.obs_to_level;
        comp.covar_map = spec.obs_loadings;
        switch (spec.kernel) {
            case RandomEffectKernel::GRM: {
                if (!Ginv) {
                    throw std::runtime_error(
                        "toRandomComponents: GRM kernel requires Ginv matrix");
                }
                Eigen::MatrixXd Kinv = Ginv->topLeftCorner(spec.n_levels, spec.n_levels);
                scratch_group_adapters.push_back(
                    std::make_unique<DenseMatrixAdapter>(std::move(Kinv)));
                comp.Qinv = scratch_group_adapters.back().get();
                break;
            }
            case RandomEffectKernel::Identity: {
                Eigen::SparseMatrix<double> I_mat(spec.n_levels, spec.n_levels);
                I_mat.setIdentity();
                scratch_group_I.push_back(std::move(I_mat));
                scratch_group_adapters.push_back(
                    std::make_unique<SparseMatrixAdapter>(scratch_group_I.back()));
                comp.Qinv = scratch_group_adapters.back().get();
                break;
            }
            case RandomEffectKernel::Custom: {
                // Dense custom kernel -> precision via inverse.
                Eigen::MatrixXd Kinv_dense =
                    spec.custom_kernel.topLeftCorner(spec.n_levels, spec.n_levels);
                Kinv_dense.diagonal().array() += 0.001;
                Eigen::LDLT<Eigen::MatrixXd> kldlt(Kinv_dense);
                if (kldlt.info() != Eigen::Success) {
                    Kinv_dense = spec.custom_kernel.inverse();
                } else {
                    Kinv_dense = kldlt.solve(
                        Eigen::MatrixXd::Identity(spec.n_levels, spec.n_levels));
                }
                scratch_group_adapters.push_back(
                    std::make_unique<DenseMatrixAdapter>(Kinv_dense));
                comp.Qinv = scratch_group_adapters.back().get();
                break;
            }
        }
        out.push_back(std::move(comp));
    }
    return out;
}

} // namespace cosmic