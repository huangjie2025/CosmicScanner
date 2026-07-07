#include "gwas_scanner.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <omp.h>
#include "data_filter.h"
#include "design.h"
#include "mvlmm_engine.h"

namespace cosmic {

std::vector<std::string> sort_chrom_labels(const std::set<std::string>& chroms_set);

static double gamma_q_local(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7) {
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < 1e-7) {
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
            }
        }
    }
    return 0.0;
}

static double chisq_prob_local(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return gamma_q_local(df / 2.0, stat / 2.0);
}


struct RrmNullModelCache {
    Eigen::LDLT<Eigen::MatrixXd> solver;
    Eigen::MatrixXd VinvX;
    Eigen::MatrixXd XtVinvX_inv;
    Eigen::VectorXd P_y;
};

struct RrmJointScanCache {
    std::vector<Eigen::VectorXd> score_weights;
    std::vector<std::vector<Eigen::MatrixXd>> info_weights;
};


static Eigen::MatrixXd project_spd_gwas(const Eigen::MatrixXd& A, double eps = 1e-8) {
    Eigen::MatrixXd S = 0.5 * (A + A.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    if (es.info() != Eigen::Success) {
        Eigen::MatrixXd out = S;
        out.diagonal().array() += eps;
        return out;
    }
    Eigen::VectorXd eval = es.eigenvalues();
    for (int i = 0; i < eval.size(); ++i) eval(i) = std::max(eval(i), eps);
    return es.eigenvectors() * eval.asDiagonal() * es.eigenvectors().transpose();
}

static Eigen::MatrixXd sanitize_symmetric_finite(const Eigen::MatrixXd& A, double diag_floor = 1e-8) {
    Eigen::MatrixXd S = A;
    for (int i = 0; i < S.rows(); ++i) {
        for (int j = 0; j < S.cols(); ++j) {
            if (!std::isfinite(S(i, j))) S(i, j) = 0.0;
        }
    }
    S = 0.5 * (S + S.transpose());
    for (int i = 0; i < S.rows(); ++i) {
        if (!std::isfinite(S(i, i)) || S(i, i) < diag_floor) S(i, i) = diag_floor;
    }
    return S;
}

static RrmNullModelCache build_rrm_null_cache(const Eigen::MatrixXd& V,
                                              const Eigen::VectorXd& y,
                                              const Eigen::MatrixXf& X) {
    RrmNullModelCache cache;
    Eigen::MatrixXd V_work = sanitize_symmetric_finite(V, 1e-8);
    const double jitters[] = {1e-8, 1e-6, 1e-4};
    bool ok = false;
    for (double jitter : jitters) {
        Eigen::MatrixXd V_try = V_work;
        V_try.diagonal().array() += jitter;
        cache.solver.compute(V_try);
        if (cache.solver.info() == Eigen::Success) {
            ok = true;
            break;
        }
    }
    if (!ok) throw std::runtime_error("RRM GWAS null model failed to factorize covariance matrix.");

    Eigen::MatrixXd Xd = X.cast<double>();
    cache.VinvX = cache.solver.solve(Xd);
    Eigen::VectorXd Vinv_y = cache.solver.solve(y);
    Eigen::MatrixXd XtVinvX = sanitize_symmetric_finite(Xd.transpose() * cache.VinvX, 1e-8);
    Eigen::LDLT<Eigen::MatrixXd> x_solver;
    bool x_ok = false;
    for (double jitter : jitters) {
        Eigen::MatrixXd Xt_try = XtVinvX;
        Xt_try.diagonal().array() += jitter;
        x_solver.compute(Xt_try);
        if (x_solver.info() == Eigen::Success) {
            x_ok = true;
            break;
        }
    }
    if (x_ok) {
        cache.XtVinvX_inv = x_solver.solve(Eigen::MatrixXd::Identity(XtVinvX.rows(), XtVinvX.cols()));
    } else {
        cache.XtVinvX_inv = XtVinvX.completeOrthogonalDecomposition().pseudoInverse();
    }
    Eigen::VectorXd beta_null = cache.XtVinvX_inv * (Xd.transpose() * Vinv_y);
    cache.P_y = Vinv_y - cache.VinvX * beta_null;
    return cache;
}

static Eigen::MatrixXd build_rrm_observation_kernel(const Eigen::MatrixXd& G_subject,
                                                    const std::vector<int>& obs_to_subject,
                                                    const Eigen::MatrixXd& Phi_obs,
                                                    const Eigen::MatrixXd& K_coeff,
                                                    const Eigen::MatrixXd& P_coeff = Eigen::MatrixXd()) {
    int n_obs = (int)obs_to_subject.size();
    Eigen::MatrixXd time_cov = Phi_obs * K_coeff * Phi_obs.transpose();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(n_obs, n_obs);
    for (int i = 0; i < n_obs; ++i) {
        int si = obs_to_subject[i];
        K(i, i) = time_cov(i, i) * G_subject(si, si);
        for (int j = 0; j < i; ++j) {
            int sj = obs_to_subject[j];
            double val = time_cov(i, j) * G_subject(si, sj);
            K(i, j) = val;
            K(j, i) = val;
        }
    }
    if (P_coeff.size() > 0) {
        Eigen::MatrixXd pe_cov = Phi_obs * P_coeff * Phi_obs.transpose();
        for (int i = 0; i < n_obs; ++i) {
            int si = obs_to_subject[i];
            K(i, i) += pe_cov(i, i);
            for (int j = 0; j < i; ++j) {
                if (obs_to_subject[j] == si) {
                    K(i, j) += pe_cov(i, j);
                    K(j, i) += pe_cov(i, j);
                }
            }
        }
    }
    return K;
}

static RrmJointScanCache build_rrm_joint_scan_cache(const RrmNullModelCache& cache,
                                                    const std::vector<int>& obs_to_subject,
                                                    const Eigen::MatrixXd& Phi_obs,
                                                    int n_subjects) {
    int n_obs = (int)obs_to_subject.size();
    int n_basis = Phi_obs.cols();
    RrmJointScanCache out;
    out.score_weights.assign(n_basis, Eigen::VectorXd::Zero(n_subjects));
    out.info_weights.resize(n_basis);
    for (int a = 0; a < n_basis; ++a) {
        out.info_weights[a].resize(n_basis);
        for (int b = 0; b < n_basis; ++b) {
            out.info_weights[a][b] = Eigen::MatrixXd::Zero(n_subjects, n_subjects);
        }
    }

    Eigen::MatrixXd Vinv = cache.solver.solve(Eigen::MatrixXd::Identity(n_obs, n_obs));
    Eigen::MatrixXd P = Vinv - cache.VinvX * cache.XtVinvX_inv * cache.VinvX.transpose();
    P = 0.5 * (P + P.transpose());

    for (int i = 0; i < n_obs; ++i) {
        int si = obs_to_subject[i];
        for (int a = 0; a < n_basis; ++a) {
            out.score_weights[a](si) += Phi_obs(i, a) * cache.P_y(i);
        }
    }

    for (int i = 0; i < n_obs; ++i) {
        int si = obs_to_subject[i];
        for (int j = 0; j < n_obs; ++j) {
            int sj = obs_to_subject[j];
            double pij = P(i, j);
            if (pij == 0.0) continue;
            for (int a = 0; a < n_basis; ++a) {
                double pia = Phi_obs(i, a) * pij;
                for (int b = 0; b < n_basis; ++b) {
                    out.info_weights[a][b](si, sj) += pia * Phi_obs(j, b);
                }
            }
        }
    }
    for (int a = 0; a < n_basis; ++a) {
        out.info_weights[a][a] = 0.5 * (out.info_weights[a][a] + out.info_weights[a][a].transpose());
        for (int b = 0; b < a; ++b) {
            Eigen::MatrixXd sym = 0.5 * (out.info_weights[a][b] + out.info_weights[b][a].transpose());
            out.info_weights[a][b] = sym;
            out.info_weights[b][a] = sym.transpose();
        }
    }
    return out;
}

struct RrmNullVarianceFit {
    Eigen::MatrixXd K_coeff;
    Eigen::MatrixXd P_coeff;
    double sigma2_e = 1.0;
    double logL = -std::numeric_limits<double>::infinity();
};

static double calc_rrm_reml_logl(const Eigen::MatrixXd& V,
                                 const Eigen::VectorXd& y,
                                 const Eigen::MatrixXf& X) {
    Eigen::LDLT<Eigen::MatrixXd> v_solver(sanitize_symmetric_finite(V, 1e-8));
    if (v_solver.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
    Eigen::MatrixXd Xd = X.cast<double>();
    Eigen::VectorXd Vinv_y = v_solver.solve(y);
    Eigen::MatrixXd VinvX = v_solver.solve(Xd);
    Eigen::MatrixXd XtVinvX = sanitize_symmetric_finite(Xd.transpose() * VinvX, 1e-8);
    Eigen::LDLT<Eigen::MatrixXd> x_solver(XtVinvX);
    if (x_solver.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
    Eigen::VectorXd beta = x_solver.solve(Xd.transpose() * Vinv_y);
    Eigen::VectorXd resid = y - Xd * beta;
    double quad = resid.dot(v_solver.solve(resid));
    double log_det_V = v_solver.vectorD().array().log().sum();
    double log_det_X = x_solver.vectorD().array().log().sum();
    return -0.5 * (log_det_V + log_det_X + quad);
}

static std::vector<Eigen::MatrixXd> build_rrm_derivative_kernels(const Eigen::MatrixXd& G_subject,
                                                                 const std::vector<int>& obs_to_subject,
                                                                 const Eigen::MatrixXd& Phi_obs) {
    int n_obs = (int)obs_to_subject.size();
    int n_basis = Phi_obs.cols();
    std::vector<Eigen::MatrixXd> derivs;
    derivs.reserve(n_basis * (n_basis + 1) / 2);
    for (int a = 0; a < n_basis; ++a) {
        for (int c = 0; c <= a; ++c) {
            Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_obs, n_obs);
            for (int i = 0; i < n_obs; ++i) {
                int si = obs_to_subject[i];
                double phi_ia = Phi_obs(i, a);
                double phi_ic = Phi_obs(i, c);
                D(i, i) = G_subject(si, si) * phi_ia * phi_ic;
                if (a != c) D(i, i) += G_subject(si, si) * phi_ic * phi_ia;
                for (int j = 0; j < i; ++j) {
                    int sj = obs_to_subject[j];
                    double val = G_subject(si, sj) * phi_ia * Phi_obs(j, c);
                    if (a != c) val += G_subject(si, sj) * phi_ic * Phi_obs(j, a);
                    D(i, j) = val;
                    D(j, i) = val;
                }
            }
            derivs.push_back(std::move(D));
        }
    }
    return derivs;
}

static std::vector<Eigen::MatrixXd> build_rrm_pe_derivative_kernels(const std::vector<int>& obs_to_subject,
                                                                    const Eigen::MatrixXd& Phi_obs) {
    int n_obs = (int)obs_to_subject.size();
    int n_basis = Phi_obs.cols();
    std::vector<Eigen::MatrixXd> derivs;
    derivs.reserve(n_basis * (n_basis + 1) / 2);
    for (int a = 0; a < n_basis; ++a) {
        for (int c = 0; c <= a; ++c) {
            Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_obs, n_obs);
            for (int i = 0; i < n_obs; ++i) {
                int si = obs_to_subject[i];
                double phi_ia = Phi_obs(i, a);
                double phi_ic = Phi_obs(i, c);
                D(i, i) += phi_ia * phi_ic;
                if (a != c) D(i, i) += phi_ic * phi_ia;
                for (int j = 0; j < i; ++j) {
                    if (obs_to_subject[j] == si) {
                        double val = phi_ia * Phi_obs(j, c);
                        if (a != c) val += phi_ic * Phi_obs(j, a);
                        D(i, j) = val;
                        D(j, i) = val;
                    }
                }
            }
            derivs.push_back(std::move(D));
        }
    }
    return derivs;
}

static RrmNullVarianceFit estimate_rrm_null_variance(const Eigen::MatrixXd& G_subject,
                                                     const std::vector<int>& obs_to_subject,
                                                     const Eigen::MatrixXd& Phi_obs,
                                                     const Eigen::VectorXd& y,
                                                     const Eigen::MatrixXf& X,
                                                     const Eigen::MatrixXd& init_K,
                                                     double init_sigma2_e,
                                                     int max_ai_iter = 8,
                                                     const Eigen::MatrixXd& init_P = Eigen::MatrixXd()) {
    int n_basis = Phi_obs.cols();
    RrmNullVarianceFit fit;
    fit.K_coeff = project_spd_gwas(init_K, 1e-6);
    fit.P_coeff = (init_P.size() == n_basis * n_basis)
        ? project_spd_gwas(init_P, 1e-6)
        : Eigen::MatrixXd::Zero(n_basis, n_basis);
    fit.sigma2_e = (init_sigma2_e > 0.0 && std::isfinite(init_sigma2_e)) ? init_sigma2_e : 1.0;

    auto build_V = [&](const Eigen::MatrixXd& K, const Eigen::MatrixXd& P, double sigma2_e) {
        Eigen::MatrixXd V = build_rrm_observation_kernel(G_subject, obs_to_subject, Phi_obs, K, P) +
                            sigma2_e * Eigen::MatrixXd::Identity((int)obs_to_subject.size(), (int)obs_to_subject.size());
        return project_spd_gwas(sanitize_symmetric_finite(V, 1e-8), 1e-8);
    };

    Eigen::MatrixXd V0 = build_V(fit.K_coeff, fit.P_coeff, fit.sigma2_e);
    fit.logL = calc_rrm_reml_logl(V0, y, X);
    if (max_ai_iter <= 0) {
        return fit;
    }

    auto basis_matrix = [](int dim, int i, int j) {
        Eigen::MatrixXd E = Eigen::MatrixXd::Zero(dim, dim);
        E(i, j) = 1.0;
        if (i != j) E(j, i) = 1.0;
        return E;
    };

    int n_basis_loc = n_basis;
    int n_k = n_basis_loc * (n_basis_loc + 1) / 2;
    int n_pe = n_k;
    int n_params = n_k + n_pe + 1;

    for (int iter = 0; iter < max_ai_iter; ++iter) {
        auto k_derivs = build_rrm_derivative_kernels(G_subject, obs_to_subject, Phi_obs);
        auto pe_derivs = build_rrm_pe_derivative_kernels(obs_to_subject, Phi_obs);
        Eigen::MatrixXd V = build_V(fit.K_coeff, fit.P_coeff, fit.sigma2_e);
        RrmNullModelCache cache = build_rrm_null_cache(V, y, X);
        Eigen::MatrixXd Vinv = cache.solver.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));
        Eigen::MatrixXd P = Vinv - cache.VinvX * cache.XtVinvX_inv * cache.VinvX.transpose();
        Eigen::VectorXd Py = cache.P_y;

        Eigen::VectorXd grad = Eigen::VectorXd::Zero(n_params);
        Eigen::MatrixXd AI = Eigen::MatrixXd::Zero(n_params, n_params);
        std::vector<Eigen::MatrixXd> PD(n_k + n_pe);
        for (int i = 0; i < n_k; ++i) {
            PD[i] = P * k_derivs[i];
            grad(i) = -0.5 * PD[i].trace() + 0.5 * Py.dot(k_derivs[i] * Py);
        }
        for (int i = 0; i < n_pe; ++i) {
            PD[n_k + i] = P * pe_derivs[i];
            grad(n_k + i) = -0.5 * PD[n_k + i].trace() + 0.5 * Py.dot(pe_derivs[i] * Py);
        }
        for (int i = 0; i < n_k + n_pe; ++i) {
            for (int j = 0; j <= i; ++j) {
                double aij = 0.5 * (PD[i].cwiseProduct(PD[j].transpose())).sum();
                AI(i, j) = aij;
                AI(j, i) = aij;
            }
        }
        grad(n_params - 1) = -0.5 * P.trace() + 0.5 * Py.squaredNorm();
        for (int i = 0; i < n_k + n_pe; ++i) {
            double aie = 0.5 * (PD[i].cwiseProduct(P.transpose())).sum();
            AI(i, n_params - 1) = aie;
            AI(n_params - 1, i) = aie;
        }
        AI(n_params - 1, n_params - 1) = 0.5 * P.squaredNorm();
        AI = project_spd_gwas(AI, 1e-8);

        Eigen::LDLT<Eigen::MatrixXd> ai_solver(AI);
        Eigen::VectorXd delta;
        if (ai_solver.info() == Eigen::Success) delta = ai_solver.solve(grad);
        else delta = AI.inverse() * grad;
        if (!delta.allFinite()) break;

        Eigen::MatrixXd oldK = fit.K_coeff;
        Eigen::MatrixXd oldP = fit.P_coeff;
        double oldSigma = fit.sigma2_e;
        double oldLogL = fit.logL;
        bool accepted = false;
        for (int step = 0; step < 12; ++step) {
            Eigen::MatrixXd candK = oldK;
            int pidx = 0;
            for (int a = 0; a < n_basis; ++a)
                for (int c = 0; c <= a; ++c)
                    candK += delta(pidx++) * basis_matrix(n_basis, a, c);
            candK = project_spd_gwas(candK, 1e-6);

            Eigen::MatrixXd candP = oldP;
            for (int a = 0; a < n_basis; ++a)
                for (int c = 0; c <= a; ++c)
                    candP += delta(pidx++) * basis_matrix(n_basis, a, c);
            candP = project_spd_gwas(candP, 1e-6);

            double candSigma = oldSigma + delta(n_params - 1);
            if (!(candSigma > 0.0) || !std::isfinite(candSigma)) {
                delta *= 0.5;
                continue;
            }
            double candLogL = calc_rrm_reml_logl(build_V(candK, candP, candSigma), y, X);
            if (std::isfinite(candLogL) && candLogL >= oldLogL - 1e-8) {
                fit.K_coeff = candK;
                fit.P_coeff = candP;
                fit.sigma2_e = candSigma;
                fit.logL = candLogL;
                accepted = true;
                break;
            }
            delta *= 0.5;
        }
        if (!accepted) break;
    }
    return fit;
}

template <typename ResultAppender>
static void scan_plink_rrm_joint(const PlinkReader& reader,
                                 const std::string& chrom_filter,
                                 const std::vector<int>& keep_subject_indices,
                                 const RrmJointScanCache& scan_cache,
                                 ResultAppender append_result) {
    int n_basis = (int)scan_cache.score_weights.size();
    const auto& bim = reader.getBimInfo();
    Eigen::VectorXf g_raw;
    Eigen::VectorXd subject_g(keep_subject_indices.size());

    for (size_t j = 0; j < bim.size(); ++j) {
        if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
        reader.getSnpGenotypes(j, g_raw);
        double sum_g = 0.0;
        int valid_n = 0;
        for (size_t s = 0; s < keep_subject_indices.size(); ++s) {
            float val = g_raw[keep_subject_indices[s]];
            subject_g((int)s) = val;
            if (!std::isnan(val)) {
                sum_g += val;
                valid_n++;
            }
        }
        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        double mean = 2.0 * freq;
        for (int s = 0; s < subject_g.size(); ++s) {
            if (std::isnan(subject_g(s))) subject_g(s) = 0.0;
            else subject_g(s) -= mean;
        }

        Eigen::MatrixXd M(n_basis, n_basis);
        Eigen::VectorXd score(n_basis);
        for (int a = 0; a < n_basis; ++a) {
            score(a) = scan_cache.score_weights[a].dot(subject_g);
            for (int b = 0; b < n_basis; ++b) {
                M(a, b) = subject_g.dot(scan_cache.info_weights[a][b] * subject_g);
            }
        }
        M = project_spd_gwas(M, 1e-8);
        Eigen::LDLT<Eigen::MatrixXd> m_solver(M);
        if (m_solver.info() != Eigen::Success) continue;
        Eigen::VectorXd beta_joint = m_solver.solve(score);
        double stat = score.dot(beta_joint);
        Eigen::MatrixXd M_inv = m_solver.solve(Eigen::MatrixXd::Identity(n_basis, n_basis));

        GwasResult res;
        res.chrom = bim[j].chrom;
        res.snp_id = bim[j].id;
        res.pos = bim[j].bp_pos;
        res.a1 = bim[j].alt;
        res.a2 = bim[j].ref;
        res.freq = freq;
        res.n_obs = (int)keep_subject_indices.size();
        res.beta = beta_joint.size() > 0 ? beta_joint(0) : 0.0;
        res.se = (M_inv.rows() > 0 && M_inv(0, 0) > 0.0) ? std::sqrt(M_inv(0, 0)) : 0.0;
        res.p_value = chisq_prob_local(stat, n_basis);
        res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
        append_result(res, beta_joint, M_inv, stat);
    }
}

template <typename Reader, typename ReadFn, typename ResultAppender>
static void scan_seq_rrm_joint(Reader& reader,
                               size_t m,
                               ReadFn read_variant,
                               const std::string& chrom_filter,
                               const std::vector<int>& keep_subject_indices,
                               const RrmJointScanCache& scan_cache,
                               ResultAppender append_result) {
    int n_basis = (int)scan_cache.score_weights.size();
    Eigen::VectorXf dosage_f;
    Eigen::VectorXd subject_g(keep_subject_indices.size());

    reader.reset();
    for (size_t j = 0; j < m; ++j) {
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        if (!read_variant(chrom, rsid, pos, ref, alt, dosage_f)) break;
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;

        double sum_g = 0.0;
        int valid_n = 0;
        for (size_t s = 0; s < keep_subject_indices.size(); ++s) {
            float val = dosage_f[keep_subject_indices[s]];
            subject_g((int)s) = val;
            if (!std::isnan(val)) {
                sum_g += val;
                valid_n++;
            }
        }
        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        double mean = 2.0 * freq;
        for (int s = 0; s < subject_g.size(); ++s) {
            if (std::isnan(subject_g(s))) subject_g(s) = 0.0;
            else subject_g(s) -= mean;
        }

        Eigen::MatrixXd M(n_basis, n_basis);
        Eigen::VectorXd score(n_basis);
        for (int a = 0; a < n_basis; ++a) {
            score(a) = scan_cache.score_weights[a].dot(subject_g);
            for (int b = 0; b < n_basis; ++b) {
                M(a, b) = subject_g.dot(scan_cache.info_weights[a][b] * subject_g);
            }
        }
        M = project_spd_gwas(M, 1e-8);
        Eigen::LDLT<Eigen::MatrixXd> m_solver(M);
        if (m_solver.info() != Eigen::Success) continue;
        Eigen::VectorXd beta_joint = m_solver.solve(score);
        double stat = score.dot(beta_joint);
        Eigen::MatrixXd M_inv = m_solver.solve(Eigen::MatrixXd::Identity(n_basis, n_basis));

        GwasResult res;
        res.chrom = chrom;
        res.snp_id = rsid;
        res.pos = pos;
        res.a1 = alt;
        res.a2 = ref;
        res.freq = freq;
        res.n_obs = (int)keep_subject_indices.size();
        res.beta = beta_joint.size() > 0 ? beta_joint(0) : 0.0;
        res.se = (M_inv.rows() > 0 && M_inv(0, 0) > 0.0) ? std::sqrt(M_inv(0, 0)) : 0.0;
        res.p_value = chisq_prob_local(stat, n_basis);
        res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
        append_result(res, beta_joint, M_inv, stat);
    }
}


void GwasScanner::runRrmScan(
    const std::vector<PhenoReader::PhenoRecord>& records,
    int rrm_order,
    const Eigen::MatrixXf& X_base,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    int n_ind = keep.keep_indices.size();
    int n_obs = records.size();
    int n_basis = rrm_order + 1;
    int rrm_curve_grid = std::max(2, options.rrm_curve_grid);
    int rrm_top = std::max(1, options.rrm_top);
    double rrm_sig_threshold = (options.rrm_sig_threshold > 0.0 && options.rrm_sig_threshold <= 1.0)
        ? options.rrm_sig_threshold : 0.05;

    std::cout << "\n=== Starting Random Regression Model (RRM) GWAS Scan ===\n";
    std::cout << "Observations: " << n_obs << ", Individuals: " << n_ind
              << ", RRM Basis: " << options.rrm_basis << ", Order: " << rrm_order << "\n";
    std::cout << "Dynamic curve grid: " << rrm_curve_grid << ", top SNP curves: " << rrm_top
              << ", window P threshold: " << rrm_sig_threshold << "\n";

    // 1. Build kept-sample lookup and collect longitudinal records by subject
    std::unordered_map<std::string, int> ind_idx_map;
    const auto* fam_ptr = (!is_bgen && !is_pgen) ? &plink_reader.getFamInfo() : nullptr;
    int total_samples = !fids.empty() ? (int)fids.size()
                     : (is_pgen ? (int)pgen_reader.getSampleIds().size()
                                : (is_bgen ? (int)bgen_reader.getSampleIds().size()
                                           : (fam_ptr ? (int)fam_ptr->size() : n_ind)));
    std::vector<std::string> kept_fids(n_ind), kept_iids(n_ind);

    for (int i = 0; i < n_ind; ++i) {
        int orig_idx = keep.keep_indices[i];
        std::string fid, iid;
        if (orig_idx >= 0 && orig_idx < static_cast<int>(fids.size()) && orig_idx < static_cast<int>(iids.size())) {
            fid = fids[orig_idx];
            iid = iids[orig_idx];
        } else if (is_pgen) {
            fid = pgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (is_bgen) {
            fid = bgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (fam_ptr) {
            fid = (*fam_ptr)[orig_idx].fid;
            iid = (*fam_ptr)[orig_idx].iid;
        }
        kept_fids[i] = fid;
        kept_iids[i] = iid;
        ind_idx_map[fid + ":" + iid] = i;
    }

    std::vector<std::vector<double>> subject_times(n_ind);
    std::vector<std::vector<double>> subject_values(n_ind);
    std::vector<double> subject_time_sum(n_ind, 0.0);
    std::vector<int> subject_obs_count(n_ind, 0);
    int valid_obs = 0;
    double tmin = std::numeric_limits<double>::max();
    double tmax = std::numeric_limits<double>::lowest();

    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = ind_idx_map.find(key);
        if (it != ind_idx_map.end() && !std::isnan(rec.value)) {
            int subj = it->second;
            subject_times[subj].push_back(rec.time);
            subject_values[subj].push_back(rec.value);
            subject_time_sum[subj] += rec.time;
            subject_obs_count[subj] += 1;
            if (rec.time < tmin) tmin = rec.time;
            if (rec.time > tmax) tmax = rec.time;
            valid_obs++;
        }
    }

    std::cout << "Matched " << valid_obs << " valid observations to " << n_ind << " genotyped individuals.\n";
    std::cout << "Time range: [" << tmin << ", " << tmax << "]\n";

    if (valid_obs == 0) {
        std::cerr << "Error: No valid observations after matching.\n";
        return;
    }

    // 2. Estimate subject-specific random regression coefficients via within-subject Legendre regression
    Eigen::MatrixXd coef_mat = Eigen::MatrixXd::Constant(n_ind, n_basis, std::numeric_limits<double>::quiet_NaN());
    int n_subjects_with_obs = 0;
    double residual_sse = 0.0;
    for (int subj = 0; subj < n_ind; ++subj) {
        int m = subject_obs_count[subj];
        if (m <= 0) continue;
        n_subjects_with_obs++;

        Eigen::VectorXd time_sub(m), y_sub(m);
        for (int j = 0; j < m; ++j) {
            time_sub(j) = subject_times[subj][j];
            y_sub(j) = subject_values[subj][j];
        }
        Eigen::MatrixXd Phi_sub = buildRrmBasisMatrix(time_sub, tmin, tmax, rrm_order, options.rrm_basis);
        Eigen::MatrixXd gram = Phi_sub.transpose() * Phi_sub;
        double ridge = std::max(1e-8, 1e-8 * gram.trace() / std::max(1, n_basis));
        gram.diagonal().array() += ridge;
        Eigen::VectorXd rhs = Phi_sub.transpose() * y_sub;
        Eigen::VectorXd coef = gram.ldlt().solve(rhs);
        if (!coef.allFinite()) {
            coef = Phi_sub.completeOrthogonalDecomposition().solve(y_sub);
        }
        if (!coef.allFinite()) {
            coef.setZero();
        }
        coef_mat.row(subj) = coef.transpose();
        Eigen::VectorXd resid = y_sub - Phi_sub * coef;
        residual_sse += resid.squaredNorm();
    }

    std::string coef_file = options.out_prefix + ".rrmcoef.txt";
    std::ofstream fcoef(coef_file);
    fcoef << "FID\tIID\tNObs\tMeanTime";
    for (int k = 0; k < n_basis; ++k) fcoef << "\tL" << k;
    fcoef << "\n";
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] <= 0) continue;
        double mean_time = subject_time_sum[subj] / std::max(1, subject_obs_count[subj]);
        fcoef << kept_fids[subj] << "\t" << kept_iids[subj] << "\t" << subject_obs_count[subj] << "\t" << mean_time;
        for (int k = 0; k < n_basis; ++k) fcoef << "\t" << coef_mat(subj, k);
        fcoef << "\n";
    }
    fcoef.close();
    std::cout << "Saved subject-level RRM coefficients to " << coef_file << std::endl;

    std::string meta_file = options.out_prefix + ".rrm.meta.txt";
    std::ofstream fmeta(meta_file);
    fmeta << "Key\tValue\n";
    fmeta << "Mode\tGWAS\n";
    fmeta << "Order\t" << rrm_order << "\n";
    fmeta << "BasisCount\t" << n_basis << "\n";
    fmeta << "TimeMin\t" << tmin << "\n";
    fmeta << "TimeMax\t" << tmax << "\n";
    fmeta << "SubjectCount\t" << n_subjects_with_obs << "\n";
    fmeta << "ObservationCount\t" << valid_obs << "\n";
    fmeta.close();
    std::cout << "Saved RRM metadata to " << meta_file << std::endl;

    DataFilter::KeepInfo basis_keep;
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] > 0) basis_keep.keep_indices.push_back(keep.keep_indices[subj]);
    }
    basis_keep.n_keep = basis_keep.keep_indices.size();
    if (basis_keep.n_keep == 0) {
        std::cerr << "Error: No subjects with valid longitudinal observations.\n";
        return;
    }

    std::vector<int> sample_to_subject_rank(total_samples, -1);
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] > 0) sample_to_subject_rank[keep.keep_indices[subj]] = (int)basis_keep.keep_indices.size(); // temp
    }
    for (int rank = 0; rank < basis_keep.n_keep; ++rank) {
        sample_to_subject_rank[basis_keep.keep_indices[rank]] = rank;
    }

    std::vector<int> obs_to_subject;
    std::vector<int> obs_to_sample;
    std::vector<double> obs_time;
    std::vector<double> obs_value;
    obs_to_subject.reserve(valid_obs);
    obs_to_sample.reserve(valid_obs);
    obs_time.reserve(valid_obs);
    obs_value.reserve(valid_obs);
    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = ind_idx_map.find(key);
        if (it == ind_idx_map.end() || std::isnan(rec.value)) continue;
        int subj = it->second;
        if (subject_obs_count[subj] <= 0) continue;
        int sample_idx = keep.keep_indices[subj];
        int subj_rank = sample_to_subject_rank[sample_idx];
        if (subj_rank < 0) continue;
        obs_to_subject.push_back(subj_rank);
        obs_to_sample.push_back(sample_idx);
        obs_time.push_back(rec.time);
        obs_value.push_back(rec.value);
    }
    int n_obs_rrm = (int)obs_value.size();
    Eigen::VectorXd y_obs(n_obs_rrm);
    Eigen::VectorXd time_obs(n_obs_rrm);
    for (int i = 0; i < n_obs_rrm; ++i) {
        y_obs(i) = obs_value[i];
        time_obs(i) = obs_time[i];
    }
    Eigen::MatrixXd Phi_obs = buildRrmBasisMatrix(time_obs, tmin, tmax, rrm_order, options.rrm_basis);

    Eigen::MatrixXf X_full;
    if (X_base.rows() == total_samples && X_base.cols() > 0) {
        X_full = X_base;
    } else if (X_base.rows() == n_ind && X_base.cols() > 0) {
        X_full = Eigen::MatrixXf::Zero(total_samples, X_base.cols());
        for (int subj = 0; subj < n_ind; ++subj) {
            X_full.row(keep.keep_indices[subj]) = X_base.row(subj);
        }
    } else {
        X_full = Eigen::MatrixXf::Ones(total_samples, 1);
    }

    int extra_basis_cols = std::max(0, n_basis - 1);
    Eigen::MatrixXf X_obs(n_obs_rrm, X_full.cols() + extra_basis_cols);
    for (int i = 0; i < n_obs_rrm; ++i) {
        X_obs.leftCols(X_full.cols()).row(i) = X_full.row(obs_to_sample[i]);
        for (int k = 1; k < n_basis; ++k) {
            X_obs(i, X_full.cols() + (k - 1)) = (float)Phi_obs(i, k);
        }
    }

    Eigen::VectorXd coef_mean = Eigen::VectorXd::Zero(n_basis);
    int coef_n = 0;
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] <= 0) continue;
        coef_mean += coef_mat.row(subj).transpose();
        coef_n++;
    }
    if (coef_n > 0) coef_mean /= coef_n;
    Eigen::MatrixXd K_coeff_init = Eigen::MatrixXd::Zero(n_basis, n_basis);
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] <= 0) continue;
        Eigen::VectorXd c = coef_mat.row(subj).transpose() - coef_mean;
        K_coeff_init.noalias() += c * c.transpose();
    }
    if (coef_n > 1) K_coeff_init /= (coef_n - 1);
    else K_coeff_init.setIdentity();
    K_coeff_init = project_spd_gwas(K_coeff_init, 1e-6);
    Eigen::MatrixXd P_coeff_init = K_coeff_init * 0.8;
    K_coeff_init *= 0.2;

    double sigma2_e_init = residual_sse / std::max(1, n_obs_rrm - coef_n * n_basis);
    if (!(sigma2_e_init > 0.0) || !std::isfinite(sigma2_e_init)) sigma2_e_init = 1.0;

    std::cout << "Running observation-level RRM null model GWAS..." << std::endl;
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.algorithm = options.grm_alg;
    grm_opts.extract_set = extract_ids_;
    grm_opts.exclude_set = exclude_ids_;
    #ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
    #endif

    GrmBuilder global_builder(grm_opts);
    if (is_bgen) {
        bgen_reader.reset();
        global_builder.compute(bgen_reader);
    } else if (is_pgen) {
        pgen_reader.reset();
        global_builder.compute(pgen_reader);
    } else {
        global_builder.compute(plink_reader);
    }
    const Eigen::MatrixXd& G_global_norm = global_builder.getGrm();
    double S_global = (options.grm_alg == 0) ? (double)global_builder.getTotalSnps() : global_builder.getSum2pq();
    auto chroms_set = global_builder.getObservedChroms();
    std::vector<std::string> sorted_chroms = sort_chrom_labels(chroms_set);
    Eigen::MatrixXd G_subject_global = DataFilter::subsetSquare(G_global_norm, basis_keep);

    std::string joint_file = options.out_prefix + ".rrm_joint.assoc.txt";
    std::ofstream fjoint(joint_file);
    fjoint << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA_L0\tSE_L0\tJOINT_CHISQ\tJOINT_P";
    for (int k = 0; k < n_basis; ++k) fjoint << "\tBETA_L" << k;
    for (int k = 0; k < n_basis; ++k) fjoint << "\tSE_L" << k;
    fjoint << "\n";

    std::string null_file = options.out_prefix + ".rrm_joint.null.txt";
    std::ofstream fnull(null_file);
    fnull << "Scope\tSigma2E\tLogL\tAIIter\tNullMode";
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\tK_" << i << "_" << j;
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\tPE_" << i << "_" << j;
    fnull << "\n";

    struct RrmJointRow {
        GwasResult res;
        Eigen::VectorXd beta;
        Eigen::MatrixXd cov;
        double stat = 0.0;
        std::string scope;
    };
    std::vector<RrmJointRow> joint_rows;

    auto append_joint_results = [&](const std::vector<GwasResult>& results,
                                    const std::vector<Eigen::VectorXd>& betas,
                                    const std::vector<Eigen::MatrixXd>& covs,
                                    const std::vector<double>& stats,
                                    const std::string& scope) {
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& res = results[i];
            fjoint << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                   << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                   << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
                   << stats[i] << "\t" << res.p_value;
            for (int k = 0; k < betas[i].size(); ++k) fjoint << "\t" << betas[i](k);
            for (int k = 0; k < covs[i].rows(); ++k) {
                double se = (covs[i](k, k) > 0.0) ? std::sqrt(covs[i](k, k)) : std::numeric_limits<double>::quiet_NaN();
                fjoint << "\t" << se;
            }
            fjoint << "\n";
            RrmJointRow row;
            row.res = res;
            row.beta = betas[i];
            row.cov = covs[i];
            row.stat = stats[i];
            row.scope = scope;
            joint_rows.push_back(std::move(row));
        }
    };

    auto scan_joint_with_grm = [&](const Eigen::MatrixXd& G_subject, const std::string& chrom_filter) {
        int rrm_ai_iter = options.rrm_ai_iter;
        if (rrm_ai_iter < 0) rrm_ai_iter = (n_obs_rrm > 1200) ? 0 : 8;
        if (rrm_ai_iter == 0) {
            std::cout << "RRM null variance: using moment-initialized fast null for "
                      << n_obs_rrm << " observations (AI refinement skipped).\n";
        }
        RrmNullVarianceFit null_fit = estimate_rrm_null_variance(
            G_subject, obs_to_subject, Phi_obs, y_obs, X_obs, K_coeff_init, sigma2_e_init, rrm_ai_iter, P_coeff_init);
        Eigen::MatrixXd V = build_rrm_observation_kernel(G_subject, obs_to_subject, Phi_obs,
                                                         null_fit.K_coeff, null_fit.P_coeff) +
                            null_fit.sigma2_e * Eigen::MatrixXd::Identity(n_obs_rrm, n_obs_rrm);
        RrmNullModelCache cache = build_rrm_null_cache(V, y_obs, X_obs);
        RrmJointScanCache scan_cache = build_rrm_joint_scan_cache(cache, obs_to_subject, Phi_obs, basis_keep.n_keep);
        fnull << (chrom_filter.empty() ? "GLOBAL" : ("LOCO_" + chrom_filter)) << "\t"
              << null_fit.sigma2_e << "\t" << null_fit.logL << "\t" << rrm_ai_iter << "\t"
              << (rrm_ai_iter == 0 ? "moment_fast" : "ai_refined");
        for (int i = 0; i < null_fit.K_coeff.rows(); ++i) {
            for (int j = 0; j < null_fit.K_coeff.cols(); ++j) fnull << "\t" << null_fit.K_coeff(i, j);
        }
        for (int i = 0; i < null_fit.P_coeff.rows(); ++i) {
            for (int j = 0; j < null_fit.P_coeff.cols(); ++j) fnull << "\t" << null_fit.P_coeff(i, j);
        }
        fnull << "\n";

        std::vector<GwasResult> results;
        std::vector<Eigen::VectorXd> betas;
        std::vector<Eigen::MatrixXd> covs;
        std::vector<double> stats;
        auto append_result = [&](const GwasResult& res, const Eigen::VectorXd& beta_joint, const Eigen::MatrixXd& cov_joint, double stat) {
            results.push_back(res);
            betas.push_back(beta_joint);
            covs.push_back(cov_joint);
            stats.push_back(stat);
        };

        if (is_bgen) {
            scan_seq_rrm_joint(
                bgen_reader,
                bgen_reader.getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return bgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage);
                },
                chrom_filter,
                basis_keep.keep_indices,
                scan_cache,
                append_result
            );
        } else if (is_pgen) {
            scan_seq_rrm_joint(
                pgen_reader,
                pgen_reader.getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return pgen_reader.readVariant(chrom, rsid, pos, ref, alt, dosage);
                },
                chrom_filter,
                basis_keep.keep_indices,
                scan_cache,
                append_result
            );
        } else {
            scan_plink_rrm_joint(
                plink_reader,
                chrom_filter,
                basis_keep.keep_indices,
                scan_cache,
                append_result
            );
        }
        append_joint_results(results, betas, covs, stats,
                             chrom_filter.empty() ? "GLOBAL" : ("LOCO_" + chrom_filter));
    };

    if (sorted_chroms.size() <= 1 || options.no_loco) {
        scan_joint_with_grm(G_subject_global, "");
    } else {
        for (const auto& chrom : sorted_chroms) {
            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            if (is_bgen) { bgen_reader.reset(); chr_builder.compute(bgen_reader); }
            else if (is_pgen) { pgen_reader.reset(); chr_builder.compute(pgen_reader); }
            else { chr_builder.compute(plink_reader); }

            double S_chr = (options.grm_alg == 0) ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
            const auto& G_chr_norm = chr_builder.getGrm();
            Eigen::MatrixXd G_chr_sub = DataFilter::subsetSquare(G_chr_norm, basis_keep);
            double S_loco = S_global - S_chr;
            Eigen::MatrixXd G_subject_loco = G_subject_global;
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                double w2 = S_chr / S_loco;
                G_subject_loco = G_subject_global * w1 - G_chr_sub * w2;
            }
            scan_joint_with_grm(G_subject_loco, chrom);
        }
    }
    fjoint.close();
    fnull.close();
    std::cout << "Observation-level RRM joint GWAS saved to " << joint_file << std::endl;
    std::cout << "Observation-level RRM null-model estimates saved to " << null_file << std::endl;

    std::string manifest_file = options.out_prefix + ".rrm.scan_manifest.txt";
    std::ofstream fmanifest(manifest_file);
    fmanifest << "Basis\tAssocFile\tSubjects\tObservations\n";
    std::map<std::string, std::vector<GwasResult>> basis_results;
    std::vector<std::string> basis_names;
    struct RrmSummaryRow {
        GwasResult best_res;
        std::string best_basis;
        double min_p = 1.0;
        double combined_chisq = 0.0;
        double combined_p = 1.0;
        int used_basis = 0;
        std::vector<double> per_basis_p;
        std::vector<double> per_basis_beta;
        std::vector<double> per_basis_se;
    };
    std::vector<RrmSummaryRow> summary_rows;

    for (int k = 0; k < n_basis; ++k) {
        Eigen::VectorXd y_full = Eigen::VectorXd::Zero(total_samples);
        for (int subj = 0; subj < n_ind; ++subj) {
            if (subject_obs_count[subj] > 0 && std::isfinite(coef_mat(subj, k))) {
                y_full(keep.keep_indices[subj]) = coef_mat(subj, k);
            }
        }

        GwasScanner::Options basis_opts = options;
        basis_opts.out_prefix = options.out_prefix + ".L" + std::to_string(k);
        GwasScanner basis_scanner(basis_opts);

        std::vector<GwasResult> basis_scan_results;
        std::cout << "Running RRM coefficient GWAS for L" << k << " on " << basis_keep.n_keep << " subjects...\n";
        basis_scanner.runLocoScan(y_full, X_full, plink_reader, bgen_reader, pgen_reader, is_bgen, is_pgen, basis_keep);
        std::ifstream fassoc(basis_opts.out_prefix + ".assoc.txt");
        std::string line;
        bool first_line = true;
        while (std::getline(fassoc, line)) {
            if (line.empty()) continue;
            if (first_line) { first_line = false; continue; }
            std::stringstream ss(line);
            GwasResult res;
            ss >> res.chrom >> res.snp_id >> res.pos >> res.a1 >> res.a2 >> res.freq >> res.n_obs >> res.beta >> res.se >> res.p_value >> res.log10_p;
            if (!res.snp_id.empty()) basis_scan_results.push_back(res);
        }
        basis_results["L" + std::to_string(k)] = basis_scan_results;
        basis_names.push_back("L" + std::to_string(k));
        fmanifest << "L" << k << "\t" << basis_opts.out_prefix << ".assoc.txt\t"
                  << basis_keep.n_keep << "\t" << valid_obs << "\n";
    }

    fmanifest.close();
    std::cout << "RRM coefficient GWAS manifest saved to " << manifest_file << std::endl;
    std::string summary_file = options.out_prefix + ".rrm.assoc.summary.txt";
    std::ofstream fsum(summary_file);
    fsum << "CHR\tSNP\tBP\tA1\tA2\tBestBasis\tBestBeta\tBestSE\tMinP\tCombinedChiSq\tCombinedP\tNumBasis";
    for (const auto& basis_name : basis_names) fsum << "\tP_" << basis_name;
    fsum << "\n";
    std::set<std::string> all_snps;
    for (const auto& kv : basis_results) {
        for (const auto& res : kv.second) all_snps.insert(res.snp_id);
    }
    for (const auto& snp : all_snps) {
        bool initialized = false;
        GwasResult best_res;
        std::string best_basis;
        double min_p = 1.0;
        double combined_chisq = 0.0;
        int used_basis = 0;
        std::vector<double> per_basis_p(basis_names.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> per_basis_beta(basis_names.size(), std::numeric_limits<double>::quiet_NaN());
        std::vector<double> per_basis_se(basis_names.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t bi = 0; bi < basis_names.size(); ++bi) {
            const auto& res_vec = basis_results[basis_names[bi]];
            for (const auto& res : res_vec) {
                if (res.snp_id == snp) {
                    per_basis_p[bi] = res.p_value;
                    per_basis_beta[bi] = res.beta;
                    per_basis_se[bi] = res.se;
                    break;
                }
            }
        }
        for (const auto& kv : basis_results) {
            for (const auto& res : kv.second) {
                if (res.snp_id != snp || std::isnan(res.p_value)) continue;
                double z = (res.se > 0.0) ? (res.beta / res.se) : 0.0;
                combined_chisq += z * z;
                used_basis++;
                if (!initialized || res.p_value < min_p) {
                    initialized = true;
                    best_res = res;
                    best_basis = kv.first;
                    min_p = res.p_value;
                }
            }
        }
        if (!initialized || used_basis == 0) continue;
        double combined_p = chisq_prob_local(combined_chisq, used_basis);
        fsum << best_res.chrom << "\t" << best_res.snp_id << "\t" << best_res.pos << "\t"
             << best_res.a1 << "\t" << best_res.a2 << "\t" << best_basis << "\t"
             << best_res.beta << "\t" << best_res.se << "\t" << min_p << "\t"
             << combined_chisq << "\t" << combined_p << "\t" << used_basis;
        for (double pval : per_basis_p) fsum << "\t" << pval;
        fsum << "\n";

        RrmSummaryRow row;
        row.best_res = best_res;
        row.best_basis = best_basis;
        row.min_p = min_p;
        row.combined_chisq = combined_chisq;
        row.combined_p = combined_p;
        row.used_basis = used_basis;
        row.per_basis_p = std::move(per_basis_p);
        row.per_basis_beta = std::move(per_basis_beta);
        row.per_basis_se = std::move(per_basis_se);
        summary_rows.push_back(std::move(row));
    }
    fsum.close();

    std::vector<double> unique_times;
    unique_times.reserve(valid_obs);
    for (int subj = 0; subj < n_ind; ++subj) {
        for (double t : subject_times[subj]) unique_times.push_back(t);
    }
    std::sort(unique_times.begin(), unique_times.end());
    unique_times.erase(std::unique(unique_times.begin(), unique_times.end()), unique_times.end());

    std::vector<double> curve_times;
    if ((int)unique_times.size() <= rrm_curve_grid) {
        curve_times = unique_times;
    } else {
        curve_times.resize(rrm_curve_grid);
        double span = tmax - tmin;
        for (int i = 0; i < rrm_curve_grid; ++i) curve_times[i] = tmin + span * i / std::max(1, rrm_curve_grid - 1);
    }

    Eigen::VectorXd time_vec_joint(curve_times.size());
    for (size_t i = 0; i < curve_times.size(); ++i) time_vec_joint((int)i) = curve_times[i];
    Eigen::MatrixXd Phi_grid_joint = buildRrmBasisMatrix(time_vec_joint, tmin, tmax, rrm_order, options.rrm_basis);

    std::sort(joint_rows.begin(), joint_rows.end(), [](const RrmJointRow& a, const RrmJointRow& b) {
        if (std::isfinite(a.res.p_value) && std::isfinite(b.res.p_value) && a.res.p_value != b.res.p_value) return a.res.p_value < b.res.p_value;
        return a.stat > b.stat;
    });

    std::string dynamic_curve_file = options.out_prefix + ".rrm.dynamic_curve.txt";
    std::string dynamic_summary_file = options.out_prefix + ".rrm.dynamic_summary.txt";
    std::ofstream fdyn(dynamic_curve_file);
    std::ofstream fdynsum(dynamic_summary_file);
    fdyn << "Rank\tScope\tCHR\tSNP\tBP\tA1\tA2\tJointChiSq\tJointP\tTime\tBeta\tSE\tP\n";
    fdynsum << "Rank\tScope\tCHR\tSNP\tBP\tA1\tA2\tJointChiSq\tJointP\tBestTime\tBestTimeP"
            << "\tBetaAtBestTime\tSEAtBestTime\tBetaStart\tBetaEnd\tDeltaBeta\tBetaRange"
            << "\tMaxAbsBeta\tEffectRatio\tPeakTime\tPattern"
            << "\tSigThreshold\tSigWindowStart\tSigWindowEnd\tSigGridPoints\tGridPoints\n";

    int top_joint_n = std::min<int>(rrm_top, joint_rows.size());
    for (int rank = 0; rank < top_joint_n; ++rank) {
        const auto& row = joint_rows[rank];
        double best_time = std::numeric_limits<double>::quiet_NaN();
        double best_p = 1.0;
        double best_beta = std::numeric_limits<double>::quiet_NaN();
        double best_se = std::numeric_limits<double>::quiet_NaN();
        double beta_start = std::numeric_limits<double>::quiet_NaN();
        double beta_end = std::numeric_limits<double>::quiet_NaN();
        double max_abs_beta = 0.0;
        double peak_time = std::numeric_limits<double>::quiet_NaN();
        double min_beta = std::numeric_limits<double>::infinity();
        double max_beta = -std::numeric_limits<double>::infinity();
        double sig_window_start = std::numeric_limits<double>::quiet_NaN();
        double sig_window_end = std::numeric_limits<double>::quiet_NaN();
        int sig_grid_points = 0;

        for (int ti = 0; ti < Phi_grid_joint.rows(); ++ti) {
            Eigen::VectorXd phi = Phi_grid_joint.row(ti).transpose();
            double beta_t = phi.dot(row.beta);
            double var_t = (phi.transpose() * row.cov * phi)(0, 0);
            double se_t = (var_t > 0.0 && std::isfinite(var_t)) ? std::sqrt(var_t) : std::numeric_limits<double>::quiet_NaN();
            double p_t = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(se_t) && se_t > 0.0) {
                double z = beta_t / se_t;
                p_t = chisq_prob_local(z * z, 1);
            }
            if (ti == 0) beta_start = beta_t;
            if (ti == Phi_grid_joint.rows() - 1) beta_end = beta_t;
            double abs_beta_t = std::abs(beta_t);
            if (abs_beta_t > max_abs_beta) {
                max_abs_beta = abs_beta_t;
                peak_time = curve_times[ti];
            }
            min_beta = std::min(min_beta, beta_t);
            max_beta = std::max(max_beta, beta_t);
            if (std::isfinite(p_t) && p_t < best_p) {
                best_p = p_t;
                best_time = curve_times[ti];
                best_beta = beta_t;
                best_se = se_t;
            }
            if (std::isfinite(p_t) && p_t < rrm_sig_threshold) {
                if (sig_grid_points == 0) sig_window_start = curve_times[ti];
                sig_window_end = curve_times[ti];
                sig_grid_points++;
            }
            fdyn << (rank + 1) << "\t" << row.scope << "\t" << row.res.chrom << "\t" << row.res.snp_id << "\t"
                 << row.res.pos << "\t" << row.res.a1 << "\t" << row.res.a2 << "\t"
                 << row.stat << "\t" << row.res.p_value << "\t" << curve_times[ti] << "\t"
                 << beta_t << "\t" << se_t << "\t" << p_t << "\n";
        }

        double delta = beta_end - beta_start;
        double range = max_beta - min_beta;
        double scale = std::max(max_abs_beta, 1e-12);
        double effect_ratio = range / scale;
        bool start_active = std::abs(beta_start) >= 0.70 * scale;
        bool end_active = std::abs(beta_end) >= 0.70 * scale;
        bool start_quiet = std::abs(beta_start) <= 0.40 * scale;
        bool end_quiet = std::abs(beta_end) <= 0.40 * scale;
        bool interior_peak = std::isfinite(peak_time) && peak_time > curve_times.front() && peak_time < curve_times.back();
        std::string pattern = "stable";
        if (beta_start * beta_end < 0.0 && std::abs(beta_start) >= 0.20 * scale && std::abs(beta_end) >= 0.20 * scale) {
            pattern = "sign_switch";
        } else if (start_active && end_quiet) {
            pattern = "early_specific";
        } else if (start_quiet && end_active) {
            pattern = "late_specific";
        } else if (interior_peak && start_quiet && end_quiet) {
            pattern = "transient";
        } else if (delta > 0.5 * scale) {
            pattern = "increasing";
        } else if (delta < -0.5 * scale) {
            pattern = "decreasing";
        } else if (range > 0.75 * scale && std::abs(delta) <= 0.35 * scale) {
            pattern = "nonlinear_dynamic";
        } else if (range > 0.5 * scale) {
            pattern = "stage_specific";
        }

        fdynsum << (rank + 1) << "\t" << row.scope << "\t" << row.res.chrom << "\t" << row.res.snp_id << "\t"
                << row.res.pos << "\t" << row.res.a1 << "\t" << row.res.a2 << "\t"
                << row.stat << "\t" << row.res.p_value << "\t" << best_time << "\t" << best_p << "\t"
                << best_beta << "\t" << best_se << "\t" << beta_start << "\t" << beta_end << "\t"
                << delta << "\t" << range << "\t" << max_abs_beta << "\t" << effect_ratio << "\t"
                << peak_time << "\t" << pattern << "\t" << rrm_sig_threshold << "\t"
                << sig_window_start << "\t" << sig_window_end << "\t"
                << sig_grid_points << "\t" << curve_times.size() << "\n";
    }
    fdyn.close();
    fdynsum.close();

    std::sort(summary_rows.begin(), summary_rows.end(), [](const RrmSummaryRow& a, const RrmSummaryRow& b) {
        if (std::isfinite(a.combined_p) && std::isfinite(b.combined_p) && a.combined_p != b.combined_p) return a.combined_p < b.combined_p;
        return a.min_p < b.min_p;
    });

    std::string curve_file = options.out_prefix + ".rrm.top_curve.txt";
    std::ofstream fcurve(curve_file);
    fcurve << "Rank\tCHR\tSNP\tBP\tA1\tA2\tBestBasis\tCombinedP\tTime\tCurveBeta\tApproxSE\tApproxP\n";
    int top_n = std::min<int>(rrm_top, summary_rows.size());
    Eigen::VectorXd time_vec(curve_times.size());
    for (size_t i = 0; i < curve_times.size(); ++i) time_vec((int)i) = curve_times[i];
    Eigen::MatrixXd Phi_grid = buildRrmBasisMatrix(time_vec, tmin, tmax, rrm_order, options.rrm_basis);
    for (int rank = 0; rank < top_n; ++rank) {
        const auto& row = summary_rows[rank];
        for (int ti = 0; ti < Phi_grid.rows(); ++ti) {
            double beta_t = 0.0;
            double var_t = 0.0;
            for (int k = 0; k < n_basis; ++k) {
                double phi = Phi_grid(ti, k);
                if (k < (int)row.per_basis_beta.size() && std::isfinite(row.per_basis_beta[k])) beta_t += phi * row.per_basis_beta[k];
                if (k < (int)row.per_basis_se.size() && std::isfinite(row.per_basis_se[k])) var_t += phi * phi * row.per_basis_se[k] * row.per_basis_se[k];
            }
            double se_t = (var_t > 0.0) ? std::sqrt(var_t) : std::numeric_limits<double>::quiet_NaN();
            double p_t = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(se_t) && se_t > 0.0) {
                double z = beta_t / se_t;
                p_t = chisq_prob_local(z * z, 1);
            }
            fcurve << (rank + 1) << "\t" << row.best_res.chrom << "\t" << row.best_res.snp_id << "\t"
                   << row.best_res.pos << "\t" << row.best_res.a1 << "\t" << row.best_res.a2 << "\t"
                   << row.best_basis << "\t" << row.combined_p << "\t" << curve_times[ti] << "\t"
                   << beta_t << "\t" << se_t << "\t" << p_t << "\n";
        }
    }
    fcurve.close();

    std::cout << "RRM GWAS stage completed: observation-level joint null model + coefficient-based auxiliary scans.\n";
    std::cout << "RRM GWAS summary saved to " << summary_file << std::endl;
    std::cout << "RRM GWAS top-hit curve summary saved to " << curve_file << std::endl;
    std::cout << "RRM dynamic effect curves saved to " << dynamic_curve_file << std::endl;
    std::cout << "RRM dynamic effect summary saved to " << dynamic_summary_file << std::endl;
}

void GwasScanner::runRrmUnbalancedScan(
    const std::vector<PhenoReader::PhenoRecord>& records,
    int rrm_order,
    const Eigen::MatrixXf& X_base,
    const Eigen::MatrixXf& X_tdep,
    const std::vector<std::string>& tdep_names,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    PlinkReader& plink_reader,
    BgenReader& bgen_reader,
    PgenReader& pgen_reader,
    bool is_bgen,
    bool is_pgen,
    const DataFilter::KeepInfo& keep
) {
    int n_ind = keep.keep_indices.size();
    int n_obs = records.size();
    int n_basis = rrm_order + 1;
    int n_tdep = (int)X_tdep.cols();
    int rrm_curve_grid = std::max(2, options.rrm_curve_grid);
    int rrm_top = std::max(1, options.rrm_top);
    double rrm_sig_threshold = (options.rrm_sig_threshold > 0.0 && options.rrm_sig_threshold <= 1.0)
        ? options.rrm_sig_threshold : 0.05;

    std::cout << "\n=== Starting Unbalanced Longitudinal (RRM) GWAS Scan ===\n";
    std::cout << "Observations: " << n_obs << ", Individuals: " << n_ind
              << ", RRM Basis: " << options.rrm_basis << ", Order: " << rrm_order
              << ", Time-dependent covariates: " << n_tdep << "\n";

    std::unordered_map<std::string, int> ind_idx_map;
    const auto* fam_ptr = (!is_bgen && !is_pgen) ? &plink_reader.getFamInfo() : nullptr;
    int total_samples = !fids.empty() ? (int)fids.size()
                     : (is_pgen ? (int)pgen_reader.getSampleIds().size()
                                : (is_bgen ? (int)bgen_reader.getSampleIds().size()
                                           : (fam_ptr ? (int)fam_ptr->size() : n_ind)));
    std::vector<std::string> kept_fids(n_ind), kept_iids(n_ind);
    for (int i = 0; i < n_ind; ++i) {
        int orig_idx = keep.keep_indices[i];
        std::string fid, iid;
        if (orig_idx >= 0 && orig_idx < static_cast<int>(fids.size()) && orig_idx < static_cast<int>(iids.size())) {
            fid = fids[orig_idx];
            iid = iids[orig_idx];
        } else if (is_pgen) {
            fid = pgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (is_bgen) {
            fid = bgen_reader.getSampleIds()[orig_idx];
            iid = fid;
        } else if (fam_ptr) {
            fid = (*fam_ptr)[orig_idx].fid;
            iid = (*fam_ptr)[orig_idx].iid;
        }
        kept_fids[i] = fid;
        kept_iids[i] = iid;
        ind_idx_map[fid + ":" + iid] = i;
    }

    std::vector<std::vector<double>> subject_times(n_ind);
    std::vector<std::vector<double>> subject_values(n_ind);
    std::vector<double> subject_time_sum(n_ind, 0.0);
    std::vector<int> subject_obs_count(n_ind, 0);
    int valid_obs = 0;
    double tmin = std::numeric_limits<double>::max();
    double tmax = std::numeric_limits<double>::lowest();

    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = ind_idx_map.find(key);
        if (it != ind_idx_map.end() && !std::isnan(rec.value)) {
            int subj = it->second;
            subject_times[subj].push_back(rec.time);
            subject_values[subj].push_back(rec.value);
            subject_time_sum[subj] += rec.time;
            subject_obs_count[subj] += 1;
            if (rec.time < tmin) tmin = rec.time;
            if (rec.time > tmax) tmax = rec.time;
            valid_obs++;
        }
    }

    std::cout << "Matched " << valid_obs << " valid observations to " << n_ind << " genotyped individuals.\n";
    std::cout << "Time range: [" << tmin << ", " << tmax << "]\n";
    if (valid_obs == 0) { std::cerr << "Error: No valid observations.\n"; return; }

    std::vector<int> obs_to_subject;
    std::vector<int> obs_to_sample;
    std::vector<double> obs_time;
    std::vector<double> obs_value;
    for (size_t ri = 0; ri < records.size(); ++ri) {
        const auto& rec = records[ri];
        std::string key = rec.fid + ":" + rec.iid;
        auto it = ind_idx_map.find(key);
        if (it == ind_idx_map.end() || std::isnan(rec.value)) continue;
        int subj = it->second;
        obs_to_subject.push_back(subj);
        obs_to_sample.push_back(keep.keep_indices[subj]);
        obs_time.push_back(rec.time);
        obs_value.push_back(rec.value);
    }
    int n_obs_rrm = (int)obs_value.size();
    Eigen::VectorXd y_obs(n_obs_rrm);
    Eigen::VectorXd time_obs(n_obs_rrm);
    for (int i = 0; i < n_obs_rrm; ++i) {
        y_obs(i) = obs_value[i];
        time_obs(i) = obs_time[i];
    }
    Eigen::MatrixXd Phi_obs = buildRrmBasisMatrix(time_obs, tmin, tmax, rrm_order, options.rrm_basis);

    Eigen::MatrixXf X_full;
    if (X_base.rows() == total_samples && X_base.cols() > 0) {
        X_full = X_base;
    } else if (X_base.rows() == n_ind && X_base.cols() > 0) {
        X_full = Eigen::MatrixXf::Zero(total_samples, X_base.cols());
        for (int subj = 0; subj < n_ind; ++subj)
            X_full.row(keep.keep_indices[subj]) = X_base.row(subj);
    } else {
        X_full = Eigen::MatrixXf::Ones(total_samples, 1);
    }

    int n_X_base = X_full.cols();
    int n_X_fix = n_X_base + n_tdep + std::max(0, n_basis - 1);
    Eigen::MatrixXf X_obs(n_obs_rrm, n_X_fix);
    for (int i = 0; i < n_obs_rrm; ++i) {
        X_obs.row(i).leftCols(n_X_base) = X_full.row(obs_to_sample[i]);
        if (n_tdep > 0)
            X_obs.row(i).segment(n_X_base, n_tdep) = X_tdep.row(i);
        for (int k = 1; k < n_basis; ++k)
            X_obs(i, n_X_base + n_tdep + (k - 1)) = (float)Phi_obs(i, k);
    }

    DataFilter::KeepInfo basis_keep;
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] > 0) basis_keep.keep_indices.push_back(keep.keep_indices[subj]);
    }
    basis_keep.n_keep = basis_keep.keep_indices.size();
    if (basis_keep.n_keep == 0) { std::cerr << "Error: No subjects with observations.\n"; return; }

    std::vector<int> sample_to_subject_rank(total_samples, -1);
    for (int rank = 0; rank < basis_keep.n_keep; ++rank)
        sample_to_subject_rank[basis_keep.keep_indices[rank]] = rank;
    std::vector<int> rrm_obs_to_subject(n_obs_rrm);
    for (int i = 0; i < n_obs_rrm; ++i) {
        int sample_idx = obs_to_sample[i];
        int subj_rank = sample_to_subject_rank[sample_idx];
        rrm_obs_to_subject[i] = (subj_rank >= 0) ? subj_rank : 0;
    }

    Eigen::MatrixXd coef_mat = Eigen::MatrixXd::Constant(n_ind, n_basis, std::numeric_limits<double>::quiet_NaN());
    int n_subjects_with_obs = 0;
    double residual_sse = 0.0;
    Eigen::VectorXd coef_mean = Eigen::VectorXd::Zero(n_basis);
    Eigen::MatrixXd bias_sum = Eigen::MatrixXd::Zero(n_basis, n_basis);
    int coef_n = 0;
    for (int subj = 0; subj < n_ind; ++subj) {
        int m = subject_obs_count[subj];
        if (m <= 0) continue;
        n_subjects_with_obs++;
        Eigen::VectorXd time_sub(m), y_sub(m);
        for (int j = 0; j < m; ++j) {
            time_sub(j) = subject_times[subj][j];
            y_sub(j) = subject_values[subj][j];
        }
        Eigen::MatrixXd Phi_sub = buildRrmBasisMatrix(time_sub, tmin, tmax, rrm_order, options.rrm_basis);
        Eigen::MatrixXd gram = Phi_sub.transpose() * Phi_sub;
        double ridge = std::max(1e-8, 1e-8 * gram.trace() / std::max(1, n_basis));
        gram.diagonal().array() += ridge;
        Eigen::LDLT<Eigen::MatrixXd> gram_solver(gram);
        Eigen::VectorXd coef = gram_solver.solve(Phi_sub.transpose() * y_sub);
        if (!coef.allFinite()) coef = Phi_sub.completeOrthogonalDecomposition().solve(y_sub);
        if (!coef.allFinite()) coef.setZero();
        coef_mat.row(subj) = coef.transpose();
        Eigen::VectorXd resid = y_sub - Phi_sub * coef;
        residual_sse += resid.squaredNorm();
        coef_mean += coef;
        Eigen::MatrixXd gram_inv;
        if (gram_solver.info() == Eigen::Success)
            gram_inv = gram_solver.solve(Eigen::MatrixXd::Identity(n_basis, n_basis));
        else
            gram_inv = gram.completeOrthogonalDecomposition().pseudoInverse();
        bias_sum += gram_inv;
        coef_n++;
    }
    if (coef_n > 0) {
        coef_mean /= coef_n;
        bias_sum /= coef_n;
    }

    double sigma2_e_init = residual_sse / std::max(1, n_obs_rrm - coef_n * n_basis);
    if (!(sigma2_e_init > 0.0) || !std::isfinite(sigma2_e_init)) sigma2_e_init = 1.0;

    Eigen::MatrixXd K_coeff_init_raw = Eigen::MatrixXd::Zero(n_basis, n_basis);
    for (int subj = 0; subj < n_ind; ++subj) {
        if (subject_obs_count[subj] <= 0) continue;
        Eigen::VectorXd c = coef_mat.row(subj).transpose() - coef_mean;
        K_coeff_init_raw.noalias() += c * c.transpose();
    }
    if (coef_n > 1) K_coeff_init_raw /= (coef_n - 1);
    else K_coeff_init_raw.setIdentity();
    Eigen::MatrixXd Ktot = project_spd_gwas(K_coeff_init_raw - sigma2_e_init * bias_sum, 1e-6);

    Eigen::MatrixXd P_coeff_init = Ktot * 0.8;
    Eigen::MatrixXd K_coeff_init = Ktot * 0.2;

    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options.min_maf;
    grm_opts.min_hwe = (float)options.min_hwe;
    grm_opts.algorithm = options.grm_alg;
    grm_opts.extract_set = extract_ids_;
    grm_opts.exclude_set = exclude_ids_;
#ifdef _OPENMP
    grm_opts.thread_num = options.n_threads;
#endif

    auto chroms_set = [&]() -> std::set<std::string> {
        std::set<std::string> s;
        if (is_bgen) {
            bgen_reader.loadVariantInfo();
            for (const auto& v : bgen_reader.getVariants()) s.insert(v.chrom);
        } else if (is_pgen) {
            for (const auto& v : pgen_reader.getVariants()) s.insert(v.chrom);
        } else {
            for (const auto& snp : plink_reader.getBimInfo()) s.insert(snp.chrom);
        }
        return s;
    }();
    std::vector<std::string> sorted_chroms = sort_chrom_labels(chroms_set);
    bool use_loco = sorted_chroms.size() > 1 && !options.no_loco;
    bool precise_loco_null = use_loco && !options.rrm_loco_global_null;

    Eigen::MatrixXd G_global_norm;
    GrmBuilder::SplitGrmResult split_grm;
    bool has_split_grm = false;
    bool split_is_subset = false;
    double S_global = 0.0;

    size_t split_n = is_bgen || is_pgen ? (size_t)basis_keep.n_keep : plink_reader.getNumSamples();
    long double split_bytes = (long double)split_n * (long double)split_n *
                              (long double)std::max<size_t>(1, sorted_chroms.size() + 1) *
                              (long double)sizeof(double);
    const long double max_split_bytes = 3.0L * 1024.0L * 1024.0L * 1024.0L;

    if (precise_loco_null && split_bytes <= max_split_bytes) {
        std::cout << "Calculating split GRMs for faster LOCO RRM (estimated memory "
                  << (double)(split_bytes / (1024.0L * 1024.0L)) << " MB)..." << std::endl;
        GrmBuilder split_builder(grm_opts);
        if (is_bgen) {
            bgen_reader.reset();
            split_grm = split_builder.computeSplit(bgen_reader, chroms_set, basis_keep.keep_indices);
            split_is_subset = true;
        } else if (is_pgen) {
            pgen_reader.reset();
            split_grm = split_builder.computeSplit(pgen_reader, chroms_set, basis_keep.keep_indices);
            split_is_subset = true;
        } else {
            split_grm = split_builder.computeSplit(plink_reader, chroms_set);
        }
        G_global_norm = split_grm.global;
        S_global = (options.grm_alg == 0) ? (double)split_grm.global_snp_count : split_grm.global_sum_2pq;
        has_split_grm = true;
    } else {
        if (precise_loco_null) {
            std::cout << "Split GRM memory estimate "
                      << (double)(split_bytes / (1024.0L * 1024.0L))
                      << " MB exceeds limit; using per-chromosome LOCO GRM builds." << std::endl;
        } else if (use_loco && options.rrm_loco_global_null) {
            std::cout << "RRM LOCO global-null mode: reusing the global BLUP null model for chromosome scans." << std::endl;
        }
        GrmBuilder global_builder(grm_opts);
        if (is_bgen) { bgen_reader.reset(); global_builder.compute(bgen_reader); }
        else if (is_pgen) { pgen_reader.reset(); global_builder.compute(pgen_reader); }
        else { global_builder.compute(plink_reader); }
        G_global_norm = global_builder.getGrm();
        S_global = (options.grm_alg == 0) ? (double)global_builder.getTotalSnps() : global_builder.getSum2pq();
        if (!use_loco) {
            chroms_set = global_builder.getObservedChroms();
            sorted_chroms = sort_chrom_labels(chroms_set);
        }
    }
    Eigen::MatrixXd G_subject_global = split_is_subset ? G_global_norm : DataFilter::subsetSquare(G_global_norm, basis_keep);

    std::string joint_file = options.out_prefix + ".rrm_unbal.assoc.txt";
    std::ofstream fjoint(joint_file);
    fjoint << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tJOINT_CHISQ\tJOINT_P";
    for (int k = 0; k < n_basis; ++k) fjoint << "\tBETA_L" << k;
    for (int k = 0; k < n_basis; ++k) fjoint << "\tSE_L" << k;
    if (n_tdep > 0) {
        for (const auto& nm : tdep_names) fjoint << "\tBETA_" << nm;
        for (const auto& nm : tdep_names) fjoint << "\tSE_" << nm;
    }
    fjoint << "\n";

    std::string null_file = options.out_prefix + ".rrm_unbal.null.txt";
    std::ofstream fnull(null_file);
    fnull << "Scope\tSigma2E\tLogL\tAIIter\tNullMode";
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\tK_" << i << "_" << j;
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\tPE_" << i << "_" << j;
    fnull << "\n";

    struct RrmJointRow {
        GwasResult res;
        Eigen::VectorXd beta;
        Eigen::MatrixXd cov;
        double stat = 0.0;
        std::string scope;
    };
    std::vector<RrmJointRow> joint_rows;

    auto append_joint_results = [&](const std::vector<GwasResult>& results,
                                    const std::vector<Eigen::VectorXd>& betas,
                                    const std::vector<Eigen::MatrixXd>& covs,
                                    const std::vector<double>& stats,
                                    const std::string& scope) {
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& res = results[i];
            fjoint << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                   << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                   << res.n_obs << "\t" << stats[i] << "\t" << res.p_value;
            for (int k = 0; k < betas[i].size(); ++k) fjoint << "\t" << betas[i](k);
            for (int k = 0; k < covs[i].rows(); ++k) {
                double se = (covs[i](k, k) > 0.0) ? std::sqrt(covs[i](k, k)) : std::numeric_limits<double>::quiet_NaN();
                fjoint << "\t" << se;
            }
            if (n_tdep > 0 && betas[i].size() > n_basis + n_tdep) {
                int leg_cols = betas[i].size() - n_tdep;
                for (int k = 0; k < n_tdep; ++k) fjoint << "\t" << betas[i](leg_cols + k);
                for (int k = 0; k < n_tdep; ++k) {
                    double se = (covs[i](leg_cols + k, leg_cols + k) > 0.0)
                        ? std::sqrt(covs[i](leg_cols + k, leg_cols + k))
                        : std::numeric_limits<double>::quiet_NaN();
                    fjoint << "\t" << se;
                }
            }
            fjoint << "\n";
            RrmJointRow row;
            row.res = res;
            row.beta = betas[i];
            row.cov = covs[i];
            row.stat = stats[i];
            row.scope = scope;
            joint_rows.push_back(std::move(row));
        }
    };

    int n_keep = basis_keep.n_keep;
    Eigen::MatrixXd Y_blup(n_keep, n_basis);
    for (int i = 0; i < n_keep; ++i) {
        int orig_idx = basis_keep.keep_indices[i];
        Y_blup.row(i) = coef_mat.row(orig_idx) - coef_mean.transpose();
    }

    Eigen::MatrixXd X_blup(n_keep, n_X_base + 1);
    X_blup.col(0).setOnes();
    if (n_X_base > 0) {
        for (int i = 0; i < n_keep; ++i) {
            int orig_idx = basis_keep.keep_indices[i];
            X_blup.row(i).tail(n_X_base) = X_full.row(orig_idx).cast<double>();
        }
    }

    std::cout << "Running BLUP-based RRM GWAS (n_keep=" << n_keep << ", T=" << n_basis << ")..." << std::endl;

    MvEigenVce::Options mv_opts;
    mv_opts.solver = MvEigenVce::VceSolver::AIREML;
    mv_opts.ai_max_iter = std::max(1, options.rrm_ai_iter < 0 ? 20 : options.rrm_ai_iter);
    mv_opts.verbose = false;
    MvEigenVce mv_vce(mv_opts);
    mv_vce.prepare(Y_blup, X_blup, G_subject_global);
    bool vce_ok = mv_vce.runNullModel();

    Eigen::MatrixXd Vg = mv_vce.getVg();
    Eigen::MatrixXd Ve = mv_vce.getVe();

    fnull << "GLOBAL\t" << 0.0 << "\t" << 0.0 << "\t"
          << mv_opts.ai_max_iter << "\t"
          << (mv_opts.solver == MvEigenVce::VceSolver::AIREML ? "ai_refined_blup" : "em_blup");
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\t" << Vg(i, j);
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j) fnull << "\t" << Ve(i, j);
    fnull << "\n";

    MvLmmEngine::Options eng_opts;
    eng_opts.em_mvgwas_test = false;
    eng_opts.use_projection = true;
    eng_opts.center_genotypes = true;
#ifdef _OPENMP
    eng_opts.thread_num = options.n_threads;
#endif
    MvLmmEngine engine(eng_opts);
    engine.prepareFromVce(mv_vce);

    auto write_mv_result = [&](const MvGwasResult& r, const std::string& scope) {
        if (!std::isfinite(r.p_value) || r.p_value < 0.0) return;
        fjoint << r.chrom << "\t" << r.snp_id << "\t" << r.pos << "\t"
               << r.a1 << "\t" << r.a2 << "\t" << r.freq << "\t"
               << n_keep << "\t" << r.wald_stat << "\t" << r.p_value;
        for (int k = 0; k < r.beta.size(); ++k) fjoint << "\t" << r.beta(k);
        for (int k = 0; k < r.se.rows(); ++k)
            fjoint << "\t" << (r.se(k, k) > 0.0 ? std::sqrt(r.se(k, k)) : std::numeric_limits<double>::quiet_NaN());
        fjoint << "\n";
        RrmJointRow row;
        row.res.chrom = r.chrom; row.res.snp_id = r.snp_id; row.res.pos = r.pos;
        row.res.a1 = r.a1; row.res.a2 = r.a2; row.res.freq = r.freq;
        row.res.n_obs = n_keep; row.res.p_value = r.p_value;
        row.beta = r.beta; row.cov = r.se; row.stat = r.wald_stat; row.scope = scope;
        joint_rows.push_back(std::move(row));
    };

    auto scan_with_engine = [&](MvLmmEngine& eng, const std::string& chrom_filter, const std::string& scope) {
        std::vector<MvGwasResult> local;
        if (is_bgen) { bgen_reader.reset(); eng.scan(bgen_reader, local, basis_keep.keep_indices, chrom_filter); }
        else if (is_pgen) { pgen_reader.reset(); eng.scan(pgen_reader, local, basis_keep.keep_indices, chrom_filter); }
        else { eng.scan(plink_reader, local, basis_keep.keep_indices, chrom_filter); }
        for (const auto& r : local) write_mv_result(r, scope);
    };

    auto run_loco_scan = [&](const Eigen::MatrixXd& G_loco, const std::string& chrom) {
        MvEigenVce loco_vce(mv_opts);
        loco_vce.prepare(Y_blup, X_blup, G_loco);
        loco_vce.runNullModel();
        MvLmmEngine loco_engine(eng_opts);
        loco_engine.prepareFromVce(loco_vce);
        scan_with_engine(loco_engine, chrom, "LOCO_" + chrom);
    };

    std::cout << "Scanning with BLUP-based RRM..." << std::endl;
    if (!use_loco) {
        scan_with_engine(engine, "", "GLOBAL");
    } else if (options.rrm_loco_global_null) {
        for (const auto& chrom : sorted_chroms) {
            std::cout << "  LOCO global-null " << chrom << std::endl;
            scan_with_engine(engine, chrom, "LOCO_GLOBAL_NULL_" + chrom);
        }
    } else {
        std::map<std::string, size_t> split_chrom_index;
        if (has_split_grm) {
            for (size_t i = 0; i < split_grm.chrom_labels.size(); ++i) {
                split_chrom_index[split_grm.chrom_labels[i]] = i;
            }
        }
        for (const auto& chrom : sorted_chroms) {
            std::cout << "  LOCO " << chrom << std::endl;
            double S_chr = 0.0;
            Eigen::MatrixXd G_chr_sub;
            if (has_split_grm) {
                auto it = split_chrom_index.find(chrom);
                if (it == split_chrom_index.end()) continue;
                size_t ci = it->second;
                S_chr = (options.grm_alg == 0)
                    ? (double)split_grm.per_chrom_snp_count[ci]
                    : split_grm.per_chrom_sum_2pq[ci];
                G_chr_sub = split_is_subset
                    ? split_grm.per_chrom_raw[ci]
                    : DataFilter::subsetSquare(split_grm.per_chrom_raw[ci], basis_keep);
            } else {
                GrmBuilder::Options chr_opts = grm_opts;
                chr_opts.chrom_filter = chrom;
                GrmBuilder chr_builder(chr_opts);
                if (is_bgen) { bgen_reader.reset(); chr_builder.compute(bgen_reader); }
                else if (is_pgen) { pgen_reader.reset(); chr_builder.compute(pgen_reader); }
                else { chr_builder.compute(plink_reader); }
                S_chr = (options.grm_alg == 0) ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
                G_chr_sub = DataFilter::subsetSquare(chr_builder.getGrm(), basis_keep) * S_chr;
            }
            double S_loco = S_global - S_chr;
            Eigen::MatrixXd G_loco = G_subject_global;
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                G_loco = G_subject_global * w1 - G_chr_sub / S_loco;
            }
            run_loco_scan(G_loco, chrom);
        }
    }

    fjoint.close();
    fnull.close();
    std::cout << "Unbalanced RRM joint GWAS (BLUP-based) saved to " << joint_file << std::endl;
    std::cout << "Unbalanced RRM null-model estimates saved to " << null_file << std::endl;

    std::vector<double> unique_times;
    for (int subj = 0; subj < n_ind; ++subj)
        for (double t : subject_times[subj]) unique_times.push_back(t);
    std::sort(unique_times.begin(), unique_times.end());
    unique_times.erase(std::unique(unique_times.begin(), unique_times.end()), unique_times.end());

    std::vector<double> curve_times;
    if ((int)unique_times.size() <= rrm_curve_grid) {
        curve_times = unique_times;
    } else {
        curve_times.resize(rrm_curve_grid);
        double span = tmax - tmin;
        for (int i = 0; i < rrm_curve_grid; ++i)
            curve_times[i] = tmin + span * i / std::max(1, rrm_curve_grid - 1);
    }

    Eigen::VectorXd time_vec_joint(curve_times.size());
    for (size_t i = 0; i < curve_times.size(); ++i) time_vec_joint((int)i) = curve_times[i];
    Eigen::MatrixXd Phi_grid_joint = buildRrmBasisMatrix(time_vec_joint, tmin, tmax, rrm_order, options.rrm_basis);

    std::sort(joint_rows.begin(), joint_rows.end(), [](const RrmJointRow& a, const RrmJointRow& b) {
        if (std::isfinite(a.res.p_value) && std::isfinite(b.res.p_value) && a.res.p_value != b.res.p_value)
            return a.res.p_value < b.res.p_value;
        return a.stat > b.stat;
    });

    std::string dynamic_curve_file = options.out_prefix + ".rrm_unbal.dynamic_curve.txt";
    std::string dynamic_summary_file = options.out_prefix + ".rrm_unbal.dynamic_summary.txt";
    std::ofstream fdyn(dynamic_curve_file);
    std::ofstream fdynsum(dynamic_summary_file);
    fdyn << "Rank\tScope\tCHR\tSNP\tBP\tA1\tA2\tJointChiSq\tJointP\tTime\tBeta\tSE\tP\n";
    fdynsum << "Rank\tScope\tCHR\tSNP\tBP\tA1\tA2\tJointChiSq\tJointP\tBestTime\tBestTimeP"
            << "\tBetaAtBestTime\tSEAtBestTime\tBetaStart\tBetaEnd\tDeltaBeta\tBetaRange"
            << "\tMaxAbsBeta\tEffectRatio\tPeakTime\tPattern"
            << "\tSigThreshold\tSigWindowStart\tSigWindowEnd\tSigGridPoints\tGridPoints\n";

    int top_joint_n = std::min<int>(rrm_top, joint_rows.size());
    for (int rank = 0; rank < top_joint_n; ++rank) {
        const auto& row = joint_rows[rank];
        double best_time = std::numeric_limits<double>::quiet_NaN();
        double best_p = 1.0;
        double best_beta = std::numeric_limits<double>::quiet_NaN();
        double best_se = std::numeric_limits<double>::quiet_NaN();
        double beta_start = std::numeric_limits<double>::quiet_NaN();
        double beta_end = std::numeric_limits<double>::quiet_NaN();
        double max_abs_beta = 0.0;
        double peak_time = std::numeric_limits<double>::quiet_NaN();
        double min_beta = std::numeric_limits<double>::infinity();
        double max_beta = -std::numeric_limits<double>::infinity();
        double sig_window_start = std::numeric_limits<double>::quiet_NaN();
        double sig_window_end = std::numeric_limits<double>::quiet_NaN();
        int sig_grid_points = 0;

        for (int ti = 0; ti < Phi_grid_joint.rows(); ++ti) {
            Eigen::VectorXd phi = Phi_grid_joint.row(ti).transpose();
            int n_beta_snp = row.beta.size() > n_tdep ? (int)row.beta.size() - n_tdep : (int)row.beta.size();
            Eigen::VectorXd beta_snp = row.beta.head(n_beta_snp);
            Eigen::MatrixXd cov_snp = row.cov.topLeftCorner(n_beta_snp, n_beta_snp);
            double beta_t = phi.dot(beta_snp);
            double var_t = (phi.transpose() * cov_snp * phi)(0, 0);
            double se_t = (var_t > 0.0 && std::isfinite(var_t)) ? std::sqrt(var_t) : std::numeric_limits<double>::quiet_NaN();
            double p_t = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(se_t) && se_t > 0.0) p_t = chisq_prob_local((beta_t / se_t) * (beta_t / se_t), 1);
            if (ti == 0) beta_start = beta_t;
            if (ti == Phi_grid_joint.rows() - 1) beta_end = beta_t;
            double abs_beta_t = std::abs(beta_t);
            if (abs_beta_t > max_abs_beta) { max_abs_beta = abs_beta_t; peak_time = curve_times[ti]; }
            min_beta = std::min(min_beta, beta_t);
            max_beta = std::max(max_beta, beta_t);
            if (std::isfinite(p_t) && p_t < best_p) { best_p = p_t; best_time = curve_times[ti]; best_beta = beta_t; best_se = se_t; }
            if (std::isfinite(p_t) && p_t < rrm_sig_threshold) {
                if (sig_grid_points == 0) sig_window_start = curve_times[ti];
                sig_window_end = curve_times[ti];
                sig_grid_points++;
            }
            fdyn << (rank + 1) << "\t" << row.scope << "\t" << row.res.chrom << "\t" << row.res.snp_id << "\t"
                 << row.res.pos << "\t" << row.res.a1 << "\t" << row.res.a2 << "\t"
                 << row.stat << "\t" << row.res.p_value << "\t" << curve_times[ti] << "\t"
                 << beta_t << "\t" << se_t << "\t" << p_t << "\n";
        }

        double delta = beta_end - beta_start;
        double range = max_beta - min_beta;
        double scale = std::max(max_abs_beta, 1e-12);
        double effect_ratio = range / scale;
        bool start_active = std::abs(beta_start) >= 0.70 * scale;
        bool end_active = std::abs(beta_end) >= 0.70 * scale;
        bool start_quiet = std::abs(beta_start) <= 0.40 * scale;
        bool end_quiet = std::abs(beta_end) <= 0.40 * scale;
        bool interior_peak = std::isfinite(peak_time) && peak_time > curve_times.front() && peak_time < curve_times.back();
        std::string pattern = "stable";
        if (beta_start * beta_end < 0.0 && std::abs(beta_start) >= 0.20 * scale && std::abs(beta_end) >= 0.20 * scale)
            pattern = "sign_switch";
        else if (start_active && end_quiet) pattern = "early_specific";
        else if (start_quiet && end_active) pattern = "late_specific";
        else if (interior_peak && start_quiet && end_quiet) pattern = "transient";
        else if (delta > 0.5 * scale) pattern = "increasing";
        else if (delta < -0.5 * scale) pattern = "decreasing";
        else if (range > 0.75 * scale && std::abs(delta) <= 0.35 * scale) pattern = "nonlinear_dynamic";
        else if (range > 0.5 * scale) pattern = "stage_specific";

        fdynsum << (rank + 1) << "\t" << row.scope << "\t" << row.res.chrom << "\t" << row.res.snp_id << "\t"
                << row.res.pos << "\t" << row.res.a1 << "\t" << row.res.a2 << "\t"
                << row.stat << "\t" << row.res.p_value << "\t" << best_time << "\t" << best_p << "\t"
                << best_beta << "\t" << best_se << "\t" << beta_start << "\t" << beta_end << "\t"
                << delta << "\t" << range << "\t" << max_abs_beta << "\t" << effect_ratio << "\t"
                << peak_time << "\t" << pattern << "\t" << rrm_sig_threshold << "\t"
                << sig_window_start << "\t" << sig_window_end << "\t"
                << sig_grid_points << "\t" << curve_times.size() << "\n";
    }
    fdyn.close();
    fdynsum.close();

    std::cout << "Unbalanced RRM GWAS completed.\n";
    std::cout << "Dynamic effect curves saved to " << dynamic_curve_file << std::endl;
    std::cout << "Dynamic effect summary saved to " << dynamic_summary_file << std::endl;
}


} // namespace cosmic
