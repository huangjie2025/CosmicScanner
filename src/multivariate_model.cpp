#include "multivariate_model.h"
#include "multivariate_helpers.h"
#include "mv_em_update.h"
#include "mv_reml_dense.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <filesystem>
#include <cctype>
#include <numeric>
#include <omp.h>
#include <unsupported/Eigen/KroneckerProduct>
#include "pheno_reader.h"
#include "covar_reader.h"
#include "data_filter.h"
#include "data_loader.h"
#include "grm_builder.h"
#include "mvlmm_engine.h"

#ifndef COSMICSCANNER_VERSION
#define COSMICSCANNER_VERSION "unknown"
#endif
#ifndef COSMICSCANNER_GIT_COMMIT
#define COSMICSCANNER_GIT_COMMIT "unknown"
#endif

namespace cosmic {

static double mv_gamma_q(double a, double x) {
    if (x <= 0.0) return 1.0;
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

double mv_chisq_prob(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return mv_gamma_q(df / 2.0, stat / 2.0);
}

bool mv_fit_gls(const Eigen::MatrixXd& UtY,
                       const Eigen::MatrixXd& UtX,
                       const Eigen::VectorXd& D,
                       const Eigen::MatrixXd& Vg,
                       const Eigen::MatrixXd& Ve,
                       Eigen::VectorXd& beta_all,
                       Eigen::MatrixXd& beta_cov,
                       double* loglik) {
    int n = UtY.rows();
    int n_traits = UtY.cols();
    int n_covars = UtX.cols();
    int pd = n_traits * n_covars;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(pd);

    std::vector<Eigen::MatrixXd> sigma_inv(n);
    std::vector<double> logdet(n, 0.0);
    for (int i = 0; i < n; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        Sigma.diagonal().array() += 1e-8;
        Eigen::LLT<Eigen::MatrixXd> llt(Sigma);
        if (llt.info() != Eigen::Success) return false;
        sigma_inv[i] = llt.solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
        logdet[i] = 2.0 * llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
        Eigen::VectorXd yi = UtY.row(i).transpose();
        for (int a = 0; a < n_covars; ++a) {
            rhs.segment(a * n_traits, n_traits) += UtX(i, a) * (sigma_inv[i] * yi);
            for (int b = 0; b < n_covars; ++b) {
                H.block(a * n_traits, b * n_traits, n_traits, n_traits) += UtX(i, a) * UtX(i, b) * sigma_inv[i];
            }
        }
    }

    H.diagonal().array() += 1e-8;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
    if (ldlt.info() != Eigen::Success) return false;
    beta_cov = ldlt.solve(Eigen::MatrixXd::Identity(pd, pd));
    beta_all = beta_cov * rhs;
    if (!beta_all.allFinite() || !beta_cov.allFinite()) return false;

    if (loglik) {
        double ll = 0.0;
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd fit = Eigen::VectorXd::Zero(n_traits);
            for (int a = 0; a < n_covars; ++a) {
                fit += UtX(i, a) * beta_all.segment(a * n_traits, n_traits);
            }
            Eigen::VectorXd resid = UtY.row(i).transpose() - fit;
            ll += -0.5 * (logdet[i] + (resid.transpose() * sigma_inv[i] * resid)(0, 0));
        }
        *loglik = ll;
    }
    return true;
}

bool mv_em_refit_rotated(const Eigen::MatrixXd& UtY,
                                const Eigen::MatrixXd& UtX,
                                const Eigen::VectorXd& D,
                                int max_iter,
                                double tol,
                                Eigen::MatrixXd& Vg,
                                Eigen::MatrixXd& Ve,
                                const Eigen::MatrixXd* init_vg,
                                const Eigen::MatrixXd* init_ve,
                                int* iterations_out,
                                double* update_norm_out,
                                std::string* stop_reason_out,
                                bool loglik_stop) {
    int n = UtY.rows();
    int d = UtY.cols();
    if (init_vg && init_ve) {
        Vg = *init_vg;
        Ve = *init_ve;
    } else {
        Eigen::MatrixXd centered = UtY.rowwise() - UtY.colwise().mean();
        Eigen::MatrixXd cov = (centered.transpose() * centered) / std::max(1, n - 1);
        Vg = 0.5 * cov;
        Ve = 0.5 * cov;
    }
    Vg = 0.5 * (Vg + Vg.transpose());
    Ve = 0.5 * (Ve + Ve.transpose());
    Vg.diagonal().array() += 1e-8;
    Ve.diagonal().array() += 1e-8;

    if (iterations_out) *iterations_out = 0;
    if (update_norm_out) *update_norm_out = std::numeric_limits<double>::quiet_NaN();
    if (stop_reason_out) *stop_reason_out = "not_started";

    double loglik_old = std::numeric_limits<double>::quiet_NaN();
    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd beta;
        Eigen::MatrixXd beta_cov;
        double loglik_new = std::numeric_limits<double>::quiet_NaN();
        if (!mv_fit_gls(UtY, UtX, D, Vg, Ve, beta, beta_cov,
                        loglik_stop ? &loglik_new : nullptr)) return false;
        if (loglik_stop && iter != 0 && std::isfinite(loglik_new) &&
            std::isfinite(loglik_old) && std::abs(loglik_new - loglik_old) < tol) {
            if (stop_reason_out) *stop_reason_out = "loglik_tolerance";
            return true;
        }
        if (loglik_stop) loglik_old = loglik_new;

        Eigen::MatrixXd Sg = Eigen::MatrixXd::Zero(d, d);
        Eigen::MatrixXd Se = Eigen::MatrixXd::Zero(d, d);
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd fit = Eigen::VectorXd::Zero(d);
            for (int c = 0; c < UtX.cols(); ++c) {
                fit += UtX(i, c) * beta.segment(c * d, d);
            }
            Eigen::VectorXd resid = UtY.row(i).transpose() - fit;
            Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
            Sigma.diagonal().array() += 1e-8;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(Sigma);
            if (ldlt.info() != Eigen::Success) return false;
            Eigen::MatrixXd Pinv = ldlt.solve(Eigen::MatrixXd::Identity(d, d));

            Sg += mv_em_genetic_contribution(D(i), Vg, Pinv, resid);

            Eigen::MatrixXd We = Ve * Pinv;
            Eigen::VectorXd ehat = We * resid;
            Se += ehat * ehat.transpose() + Ve - We * Ve;

            if (UtX.cols() > 0) {
                Eigen::VectorXd xi = UtX.row(i).transpose();
                Eigen::MatrixXd Xi = Eigen::kroneckerProduct(
                    xi.transpose(), Eigen::MatrixXd::Identity(d, d));
                Eigen::MatrixXd beta_unc = Xi * beta_cov * Xi.transpose();
                Eigen::MatrixXd Wg = Vg * D(i) * Pinv;
                if (D(i) > 1e-12) {
                    Sg += Wg * beta_unc * Wg.transpose() / D(i);
                }
                Se += We * beta_unc * We.transpose();
            }
        }

        Eigen::MatrixXd Vg_new = 0.5 * (Sg / n + (Sg / n).transpose());
        Eigen::MatrixXd Ve_new = 0.5 * (Se / n + (Se / n).transpose());
        Vg_new.diagonal().array() += 1e-8;
        Ve_new.diagonal().array() += 1e-8;
        double err = ((Vg_new - Vg).squaredNorm() + (Ve_new - Ve).squaredNorm()) /
                     std::max(1.0, 2.0 * d * d);
        Vg = Vg_new;
        Ve = Ve_new;
        if (!Vg.allFinite() || !Ve.allFinite()) return false;
        if (iterations_out) *iterations_out = iter + 1;
        if (update_norm_out) *update_norm_out = err;
        if (!loglik_stop && err < tol) {
            if (stop_reason_out) *stop_reason_out = "tolerance";
            return true;
        }
    }
    if (iterations_out) *iterations_out = max_iter;
    if (stop_reason_out) *stop_reason_out = "max_iterations";
    return true;
}

double mv_chisq_stat_from_p(double p, int df) {
    if (!(p > 0.0) || p >= 1.0) return 0.0;
    double lo = 0.0;
    double hi = std::max(10.0, (double)df);
    while (mv_chisq_prob(hi, df) > p && hi < 1e6) hi *= 2.0;
    for (int i = 0; i < 80; ++i) {
        double mid = 0.5 * (lo + hi);
        if (mv_chisq_prob(mid, df) > p) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

double mv_gc_correct_p(double p, double lambda, int df) {
    if (!std::isfinite(p) || lambda <= 1.0) return p;
    double stat = mv_chisq_stat_from_p(p, df);
    return mv_chisq_prob(stat / lambda, df);
}

static Eigen::MatrixXd mv_trait_derivative_matrix(int d, int idx) {
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(d, d);
    int k = 0;
    for (int i = 0; i < d; ++i) {
        for (int j = i; j < d; ++j) {
            if (k == idx) {
                M(i, j) = 1.0;
                M(j, i) = 1.0;
                return M;
            }
            ++k;
        }
    }
    return M;
}

static Eigen::MatrixXd mv_load_dense_matrix(const std::string& path, int n) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open matrix file: " + path);
    Eigen::MatrixXd M(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> M(i, j))) throw std::runtime_error("Matrix file has fewer values than expected: " + path);
        }
    }
    return M;
}

static Eigen::MatrixXd mv_parse_square_matrix_arg(const std::string& s, int d, const std::string& name) {
    Eigen::MatrixXd M(d, d);
    std::stringstream ss(s);
    std::string item;
    std::vector<double> vals;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) vals.push_back(std::stod(item));
    }
    if ((int)vals.size() != d * d) throw std::runtime_error(name + " expects d*d comma-separated values");
    for (int i = 0, k = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j, ++k) M(i, j) = vals[k];
    }
    return M;
}

static void mv_write_matrix(std::ostream& os, const Eigen::MatrixXd& M) {
    os << std::setprecision(12);
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            if (i || j) os << ',';
            os << M(i, j);
        }
    }
}

static void mv_add_kron(Eigen::MatrixXd& out,
                        const Eigen::VectorXd& x,
                        const Eigen::MatrixXd& A,
                        double scale = 1.0) {
    const int c = x.size();
    const int d = A.rows();
    for (int i = 0; i < c; ++i) {
        for (int j = 0; j < c; ++j) {
            out.block(i * d, j * d, d, d).noalias() += scale * x(i) * x(j) * A;
        }
    }
}

MvCrtFactors mv_compute_crt_factors(const Eigen::MatrixXd& UtX,
                                           const Eigen::VectorXd& D,
                                           const Eigen::MatrixXd& Vg,
                                           const Eigen::MatrixXd& Ve) {
    MvCrtFactors out;
    const int n = UtX.rows();
    const int c = UtX.cols();
    const int d = Vg.rows();
    const int dc = c * d;
    const int v = d * (d + 1) / 2;
    if (n == 0 || c <= 1 || d == 0 || Ve.rows() != d || D.size() != n) return out;

    std::vector<Eigen::MatrixXd> Sinv(n);
    for (int i = 0; i < n; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        Sigma.diagonal().array() += 1e-8;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Sigma);
        if (ldlt.info() != Eigen::Success) return out;
        Sinv[i] = ldlt.solve(Eigen::MatrixXd::Identity(d, d));
        if (!Sinv[i].allFinite()) return out;
    }

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dc, dc);
    for (int i = 0; i < n; ++i) {
        mv_add_kron(Q, UtX.row(i).transpose(), Sinv[i]);
    }
    Q.diagonal().array() += 1e-10;
    Eigen::LDLT<Eigen::MatrixXd> q_ldlt(Q);
    if (q_ldlt.info() != Eigen::Success) return out;
    Eigen::MatrixXd Qi = q_ldlt.solve(Eigen::MatrixXd::Identity(dc, dc));
    if (!Qi.allFinite()) return out;

    Eigen::MatrixXd snp_block = Qi.block((c - 1) * d, (c - 1) * d, d, d);
    Eigen::LDLT<Eigen::MatrixXd> snp_ldlt(snp_block);
    if (snp_ldlt.info() != Eigen::Success) return out;
    Eigen::MatrixXd Qi_si = snp_ldlt.solve(Eigen::MatrixXd::Identity(d, d));
    if (!Qi_si.allFinite()) return out;

    std::vector<Eigen::MatrixXd> Dpar(v);
    for (int idx = 0; idx < v; ++idx) Dpar[idx] = mv_trait_derivative_matrix(d, idx);

    std::vector<Eigen::MatrixXd> M_g(v, Eigen::MatrixXd::Zero(dc, dc));
    std::vector<Eigen::MatrixXd> M_e(v, Eigen::MatrixXd::Zero(dc, dc));
    std::vector<Eigen::MatrixXd> QiM_g(v), QiM_e(v);
    for (int p = 0; p < v; ++p) {
        for (int i = 0; i < n; ++i) {
            Eigen::MatrixXd A_g = Sinv[i] * (D(i) * Dpar[p]) * Sinv[i];
            Eigen::MatrixXd A_e = Sinv[i] * Dpar[p] * Sinv[i];
            Eigen::VectorXd xi = UtX.row(i).transpose();
            mv_add_kron(M_g[p], xi, A_g);
            mv_add_kron(M_e[p], xi, A_e);
        }
        QiM_g[p] = Qi * M_g[p];
        QiM_e[p] = Qi * M_e[p];
    }

    Eigen::MatrixXd info = Eigen::MatrixXd::Zero(2 * v, 2 * v);
    for (int p1 = 0; p1 < v; ++p1) {
        for (int p2 = p1; p2 < v; ++p2) {
            double gg = 0.0, ge = 0.0, ee = 0.0;
            for (int i = 0; i < n; ++i) {
                Eigen::MatrixXd G1 = D(i) * Dpar[p1];
                Eigen::MatrixXd G2 = D(i) * Dpar[p2];
                const Eigen::MatrixXd& S = Sinv[i];
                gg += 0.5 * (S * G1 * S * G2).trace();
                ge += 0.5 * (S * G1 * S * Dpar[p2]).trace();
                ee += 0.5 * (S * Dpar[p1] * S * Dpar[p2]).trace();
            }
            info(p1, p2) = info(p2, p1) = gg;
            info(p1, p2 + v) = info(p2 + v, p1) = ge;
            info(p1 + v, p2 + v) = info(p2 + v, p1 + v) = ee;
        }
    }
    info.diagonal().array() += 1e-8;
    Eigen::LDLT<Eigen::MatrixXd> info_ldlt(info);
    if (info_ldlt.info() != Eigen::Success) return out;
    Eigen::MatrixXd Hcov = info_ldlt.solve(Eigen::MatrixXd::Identity(2 * v, 2 * v));
    if (!Hcov.allFinite()) return out;

    auto sub_snp = [&](const Eigen::MatrixXd& M) {
        return M.block((c - 1) * d, (c - 1) * d, d, d);
    };

    double B = 0.0, C = 0.0, DD = 0.0;
    for (int v1 = 0; v1 < v; ++v1) {
        Eigen::MatrixXd QiMQi_g1 = QiM_g[v1] * Qi;
        Eigen::MatrixXd QiMQi_e1 = QiM_e[v1] * Qi;
        Eigen::MatrixXd Cg1m = sub_snp(QiMQi_g1) * Qi_si;
        Eigen::MatrixXd Ce1m = sub_snp(QiMQi_e1) * Qi_si;
        double trCg1 = -Cg1m.trace();
        double trCe1 = -Ce1m.trace();

        for (int v2 = v1; v2 < v; ++v2) {
            Eigen::MatrixXd QiMQi_g2 = QiM_g[v2] * Qi;
            Eigen::MatrixXd QiMQi_e2 = QiM_e[v2] * Qi;
            Eigen::MatrixXd Cg2m = sub_snp(QiMQi_g2) * Qi_si;
            Eigen::MatrixXd Ce2m = sub_snp(QiMQi_e2) * Qi_si;
            double trCg2 = -Cg2m.trace();
            double trCe2 = -Ce2m.trace();

            double trCC_gg = (Cg1m * Cg2m).trace();
            double trCC_ge = (Cg1m * Ce2m + Ce1m * Cg2m).trace();
            double trCC_ee = (Ce1m * Ce2m).trace();

            Eigen::MatrixXd qmq_gg = sub_snp(QiM_g[v1] * QiMQi_g2) * Qi_si;
            Eigen::MatrixXd qmq_ge = sub_snp(QiM_g[v1] * QiMQi_e2 + QiM_e[v1] * QiMQi_g2) * Qi_si;
            Eigen::MatrixXd qmq_ee = sub_snp(QiM_e[v1] * QiMQi_e2) * Qi_si;
            double trB_gg = -qmq_gg.trace();
            double trB_ge = -qmq_ge.trace();
            double trB_ee = -qmq_ee.trace();

            Eigen::MatrixXd MM_gg = Eigen::MatrixXd::Zero(dc, dc);
            Eigen::MatrixXd MM_ge = Eigen::MatrixXd::Zero(dc, dc);
            Eigen::MatrixXd MM_ee = Eigen::MatrixXd::Zero(dc, dc);
            for (int i = 0; i < n; ++i) {
                Eigen::VectorXd xi = UtX.row(i).transpose();
                Eigen::MatrixXd G1 = D(i) * Dpar[v1];
                Eigen::MatrixXd G2 = D(i) * Dpar[v2];
                Eigen::MatrixXd E1 = Dpar[v1];
                Eigen::MatrixXd E2 = Dpar[v2];
                mv_add_kron(MM_gg, xi, Sinv[i] * G1 * Sinv[i] * G2 * Sinv[i]);
                mv_add_kron(MM_ge, xi, Sinv[i] * G1 * Sinv[i] * E2 * Sinv[i] +
                                      Sinv[i] * E1 * Sinv[i] * G2 * Sinv[i]);
                mv_add_kron(MM_ee, xi, Sinv[i] * E1 * Sinv[i] * E2 * Sinv[i]);
            }
            trB_gg += (sub_snp(Qi * MM_gg * Qi) * Qi_si).trace();
            trB_ge += 2.0 * (sub_snp(Qi * MM_ge * Qi) * Qi_si).trace();
            trB_ee += (sub_snp(Qi * MM_ee * Qi) * Qi_si).trace();

            double trD_gg = 2.0 * trB_gg;
            double trD_ge = 2.0 * trB_ge;
            double trD_ee = 2.0 * trB_ee;

            double h_gg = Hcov(v1, v2);
            double h_ge = Hcov(v1, v2 + v);
            double h_ee = Hcov(v1 + v, v2 + v);
            double b_add = h_gg * trB_gg + h_ge * trB_ge + h_ee * trB_ee;
            double c_add = h_gg * (trCC_gg + 0.5 * trCg1 * trCg2) +
                           h_ge * (trCC_ge + 0.5 * trCg1 * trCe2 + 0.5 * trCe1 * trCg2) +
                           h_ee * (trCC_ee + 0.5 * trCe1 * trCe2);
            double d_add = h_gg * (trCC_gg + 0.5 * trD_gg) +
                           h_ge * (trCC_ge + 0.5 * trD_ge) +
                           h_ee * (trCC_ee + 0.5 * trD_ee);
            B += b_add; C += c_add; DD += d_add;
            if (v1 != v2) {
                B += b_add; C += c_add; DD += d_add;
            }
        }
    }

    out.a = 2.0 * DD - C;
    out.b = 2.0 * B;
    out.c = C;
    out.ok = std::isfinite(out.a) && std::isfinite(out.b) && std::isfinite(out.c);
    return out;
}

double mv_pcrt(int mode, int df, double p, const MvCrtFactors& crt) {
    if (!crt.ok || !std::isfinite(p) || !(p > 0.0) || p >= 1.0) return p;
    double chisq = mv_chisq_stat_from_p(p, df);
    double chisq_crt = chisq;
    double q = (double)df;
    if (mode == 1) {
        double a = crt.c / (2.0 * q * (q + 2.0));
        double b = 1.0 + (crt.a + crt.b) / (2.0 * q);
        if (std::abs(a) < 1e-12) {
            chisq_crt = chisq / std::max(1e-12, b);
        } else {
            double disc = b * b + 4.0 * a * chisq;
            if (disc >= 0.0) chisq_crt = (-b + std::sqrt(disc)) / (2.0 * a);
        }
    } else if (mode == 2) {
        chisq_crt = chisq / std::max(1e-12, 1.0 + crt.a / (2.0 * q));
    } else {
        chisq_crt = chisq;
    }
    if (!std::isfinite(chisq_crt) || chisq_crt < 0.0) return p;
    return mv_chisq_prob(chisq_crt, df);
}

static bool contains_mv_test(const std::string& tests, const std::string& name) {
    std::stringstream ss(tests);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token == name || token == "all") return true;
    }
    return false;
}

static void mv_apply_gemma_mode(MultivariateModel::Options& opts, int mode) {
    if (mode < 1 || mode > 4) {
        throw std::runtime_error("--mv-gemma-mode must be 1 (Wald), 2 (LRT), 3 (Score), or 4 (all)");
    }
    opts.gemma_mode = mode;
    if (mode == 1) {
        opts.test_wald = true; opts.test_lrt = false; opts.test_score = false;
    } else if (mode == 2) {
        opts.test_wald = false; opts.test_lrt = true; opts.test_score = false;
    } else if (mode == 3) {
        opts.test_wald = false; opts.test_lrt = false; opts.test_score = true;
    } else {
        opts.test_wald = true; opts.test_lrt = true; opts.test_score = true;
    }
}

static std::unordered_set<std::string> mv_load_snp_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open SNP list: " + path);
    std::string id;
    while (in >> id) ids.insert(id);
    return ids;
}

bool mv_scan_filter_passes(const PlinkBimInfo& var,
                                  const MultivariateModel::Options& options) {
    if (!options.scan_chrom.empty() && var.chrom != options.scan_chrom) return false;
    if (options.scan_from_bp >= 0 && var.bp_pos < options.scan_from_bp) return false;
    if (options.scan_to_bp >= 0 && var.bp_pos > options.scan_to_bp) return false;
    if (!options.extract_set.empty() && options.extract_set.find(var.id) == options.extract_set.end()) return false;
    if (!options.exclude_set.empty() && options.exclude_set.find(var.id) != options.exclude_set.end()) return false;
    return true;
}

static std::vector<std::string> mv_split_tsv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, '\t')) fields.push_back(token);
    return fields;
}

static std::vector<std::string> mv_split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) fields.push_back(token);
    return fields;
}

static bool mv_try_parse_double(const std::string& s, double& value) {
    try {
        size_t pos = 0;
        value = std::stod(s, &pos);
        while (pos < s.size() && std::isspace((unsigned char)s[pos])) ++pos;
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

static Eigen::MatrixXd mv_read_numeric_csv_matrix(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open matrix CSV: " + path);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto fields = mv_split_csv_line(line);
        std::vector<double> numeric;
        numeric.reserve(fields.size());
        for (const auto& field : fields) {
            double value = 0.0;
            if (mv_try_parse_double(field, value)) numeric.push_back(value);
        }
        if (numeric.empty()) continue;
        rows.push_back(std::move(numeric));
    }
    if (rows.empty()) throw std::runtime_error("Matrix CSV has no numeric rows: " + path);
    const size_t cols = rows.front().size();
    for (const auto& row : rows) {
        if (row.size() != cols) throw std::runtime_error("Matrix CSV has ragged numeric rows: " + path);
    }
    Eigen::MatrixXd M(rows.size(), cols);
    for (int i = 0; i < (int)rows.size(); ++i) {
        for (int j = 0; j < (int)cols; ++j) M(i, j) = rows[i][j];
    }
    return M;
}

static Eigen::VectorXd mv_load_pheno_column_by_id(const std::string& path,
                                                   int one_based_col,
                                                   const std::vector<std::string>& fids,
                                                   const std::vector<std::string>& iids) {
    if (one_based_col <= 0) throw std::runtime_error("Environment column must be 1-indexed and positive");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open phenotype file for environment column: " + path);
    std::unordered_map<std::string, double> values;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> cols;
        std::string token;
        while (ss >> token) cols.push_back(token);
        if ((int)cols.size() < std::max(2, one_based_col)) continue;
        if (first) {
            first = false;
            std::string c0 = cols[0];
            std::transform(c0.begin(), c0.end(), c0.begin(), [](unsigned char c){ return std::tolower(c); });
            if (c0 == "fid" || c0 == "#fid") continue;
        }
        double value = std::numeric_limits<double>::quiet_NaN();
        if (!mv_try_parse_double(cols[one_based_col - 1], value)) continue;
        values[cols[0] + "\t" + cols[1]] = value;
    }
    Eigen::VectorXd out(fids.size());
    out.setConstant(std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < fids.size(); ++i) {
        auto it = values.find(fids[i] + "\t" + iids[i]);
        if (it != values.end()) out((int)i) = it->second;
    }
    return out;
}

static std::unordered_map<std::string, double> mv_read_score_prepass(const std::string& path) {
    std::unordered_map<std::string, double> scores;
    std::ifstream ifs(path);
    if (!ifs) return scores;
    std::string line;
    if (!std::getline(ifs, line)) return scores;
    auto header = mv_split_tsv_line(line);
    int snp_col = -1, score_col = -1;
    for (int i = 0; i < (int)header.size(); ++i) {
        if (header[i] == "SNP") snp_col = i;
        if (header[i] == "P_SCORE") score_col = i;
    }
    if (snp_col < 0 || score_col < 0) return scores;
    while (std::getline(ifs, line)) {
        auto fields = mv_split_tsv_line(line);
        if ((int)fields.size() <= std::max(snp_col, score_col)) continue;
        try {
            double p = std::stod(fields[score_col]);
            if (std::isfinite(p)) scores[fields[snp_col]] = p;
        } catch (...) {
        }
    }
    return scores;
}

MultivariateModel::Options MultivariateModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.scan_tier = cfg.scan_tier;
    opts.out_prefix = cfg.out_prefix;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.grm_alg = cfg.grm_alg;
    if (!cfg.multivariate_compatible.empty()) {
        if (cfg.multivariate_compatible == "emmvgwas") {
            opts.grm_alg = 6;
            opts.use_covariate_projection = false;
            opts.use_emmvgwas_lambda = true;
            opts.compatibility_mode = "emmvgwas";
        } else {
            throw std::runtime_error("Unknown --mv-compatible: " + cfg.multivariate_compatible + ". Supported value: emmvgwas");
        }
    }
    #ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
    #endif
    opts.no_loco = cfg.no_loco;
    if (!opts.use_emmvgwas_lambda) opts.use_covariate_projection = cfg.use_covariate_projection;
    opts.use_emmv = true;
    opts.use_gemma_refit = false;
    opts.use_gemma_compat = false;
    opts.use_gemma_full = false;
    opts.use_gemma_exact = false;
    opts.use_gemma_gxe = false;
    opts.use_hybrid = false;
    std::string method = cfg.multivariate_method;
    opts.requested_method = method;
    if (!opts.scan_tier.empty() && cfg.multivariate_method_explicit) {
        const std::string expected = opts.scan_tier == "approx" ? "wald"
                                   : opts.scan_tier == "semi" ? "emmv" : "gemma-exact";
        if (method != expected && !(opts.scan_tier == "exact" && method == "gemma-exact-full")) {
            throw std::runtime_error("--scan-tier " + opts.scan_tier + " conflicts with --mv-method " + method);
        }
    }
    if (opts.scan_tier == "approx") method = "wald";
    else if (opts.scan_tier == "semi") method = "emmv";
    else if (opts.scan_tier == "exact") method = "gemma-exact";
    if (opts.scan_tier == "exact" && !cfg.no_loco) {
        throw std::runtime_error("Multivariate --scan-tier exact currently requires --no-loco");
    }
    if (opts.scan_tier == "exact" && opts.use_emmvgwas_lambda) {
        throw std::runtime_error("--scan-tier exact cannot be combined with --mv-compatible emmvgwas; EMmvGWAS is the semi tier");
    }
    if (opts.scan_tier == "approx" && opts.use_emmvgwas_lambda) {
        throw std::runtime_error("--scan-tier approx cannot be combined with --mv-compatible emmvgwas");
    }
    if (method == "wald") {
        opts.use_emmv = false;
    } else if (method == "gemma") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
        opts.deprecated_alias = true;
        std::cerr << "Warning: --mv-method gemma is deprecated experimental GEMMA-style mode; use --mv-method gemma-compat for GEMMA default finite-EM compatibility.\n";
    } else if (method == "gemma-compat" || method == "gemma-default") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
        opts.use_gemma_compat = true;
        if (method == "gemma-default") {
            opts.deprecated_alias = true;
            std::cerr << "Warning: --mv-method gemma-default is an alias; preferred spelling is --mv-method gemma-compat.\n";
        }
    } else if (method == "gemma-full") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
        opts.use_gemma_full = true;
        opts.deprecated_alias = true;
        std::cerr << "Warning: --mv-method gemma-full is deprecated experimental GEMMA-style mode; use --mv-method gemma-compat for GEMMA default finite-EM compatibility.\n";
    } else if (method == "gemma-exact" || method == "gemma-exact-full") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
        opts.use_gemma_full = true;
        opts.use_gemma_exact = true;
        opts.deprecated_alias = true;
        std::cerr << "Warning: --mv-method " << method << " is deprecated experimental GEMMA-style mode and is not a validated fully converged exact REML mode.\n";
    } else if (method == "gemma-gxe" || method == "gxe"
               || method == "gemma-gxe-full" || method == "gxe-full"
               || method == "hybrid") {
        throw std::runtime_error("--mv-method " + method
            + " is unavailable in this focused release. Supported values: wald, emmv, gemma, gemma-compat, gemma-full, gemma-exact-full, eb");
    } else if (method == "eb" || method == "eb-lmm" || method == "lmm-eb") {
        opts.use_eb = true;
        opts.use_emmv = false;
    } else if (method != "emmv") {
        throw std::runtime_error("Unknown --mv-method: " + method + ". Supported values: wald, emmv, gemma, gemma-compat, gemma-full, gemma-exact-full, eb");
    }
    opts.method_name = method;
    if (opts.use_gemma_compat) opts.method_name = "gemma-compat";
    if (opts.use_emmvgwas_lambda) opts.method_name = "emmv";
    if (cfg.multivariate_eb || cfg.eb_lmm) { opts.use_eb = true; opts.use_emmv = false; }
    opts.test_wald = contains_mv_test(cfg.multivariate_tests, "wald");
    opts.test_lrt = contains_mv_test(cfg.multivariate_tests, "lrt");
    opts.test_score = contains_mv_test(cfg.multivariate_tests, "score");
    if (cfg.multivariate_gemma_mode != 4 || cfg.multivariate_tests == "wald,lrt,score") {
        mv_apply_gemma_mode(opts, cfg.multivariate_gemma_mode);
    }
    if (!(cfg.multivariate_pnr >= 0.0 && cfg.multivariate_pnr <= 1.0)) {
        throw std::runtime_error("--mv-pnr must be between 0 and 1");
    }
    opts.pnr = cfg.multivariate_pnr;
    if (cfg.multivariate_emi <= 0 || cfg.multivariate_nri <= 0 || cfg.multivariate_emp <= 0.0 || cfg.multivariate_nrp <= 0.0) {
        throw std::runtime_error("--mv-emi/--mv-nri must be positive and --mv-emp/--mv-nrp must be > 0");
    }
    opts.gemma_emi = cfg.multivariate_emi;
    opts.gemma_nri = cfg.multivariate_nri;
    opts.gemma_emp = cfg.multivariate_emp;
    opts.gemma_nrp = cfg.multivariate_nrp;
    opts.gemma_crt = cfg.multivariate_crt;
    opts.null_only = cfg.multivariate_null_only;
    opts.gemma_force_iterations = cfg.multivariate_gemma_force_iterations;
    opts.k_file = cfg.multivariate_k_file;
    opts.compare_vg = cfg.multivariate_compare_vg;
    opts.compare_ve = cfg.multivariate_compare_ve;
    opts.init_vg = cfg.multivariate_init_vg;
    opts.init_ve = cfg.multivariate_init_ve;
    opts.reml_maxit = std::max(cfg.reml_maxit, opts.use_gemma_exact ? cfg.multivariate_emi : ((opts.use_emmv || opts.use_gemma_refit) ? 1000 : 100));
    opts.reml_tol = (opts.use_emmv || opts.use_gemma_refit) ? std::min(cfg.reml_tol, 1e-8) : cfg.reml_tol;
    opts.refine_top = cfg.multivariate_refine_top;
    opts.refine_p = cfg.multivariate_refine_p;
    opts.refine_discrepancy_top = cfg.multivariate_refine_discrepancy_top;
    opts.refine_discrepancy = cfg.multivariate_refine_discrepancy;
    opts.refine_score_prepass = cfg.multivariate_refine_score_prepass;
    opts.refit_p = cfg.multivariate_refit_p;
    opts.use_aireml = cfg.multivariate_use_aireml;
    opts.em_min_iter = cfg.multivariate_em_min_iter;
    opts.ai_max_iter = cfg.multivariate_ai_max_iter;
    opts.ai_tol = cfg.multivariate_ai_tol;
    opts.output_vbeta = cfg.multivariate_output_vbeta;
    opts.output_emmvgwas = cfg.multivariate_emmvgwas_output;
    opts.emmvgwas_vc_only = cfg.multivariate_emmvgwas_vc_only;
    opts.emmvgwas_G_file = cfg.multivariate_emmvgwas_G_file;
    opts.emmvgwas_R_file = cfg.multivariate_emmvgwas_R_file;
    if (opts.emmvgwas_G_file.empty() != opts.emmvgwas_R_file.empty()) {
        throw std::runtime_error("--mv-emmvgwas-G and --mv-emmvgwas-R must be provided together");
    }
    if (!opts.emmvgwas_G_file.empty() && !opts.use_emmvgwas_lambda) {
        throw std::runtime_error("--mv-emmvgwas-G/--mv-emmvgwas-R require --mv-compatible emmvgwas");
    }
    if (opts.emmvgwas_vc_only && !opts.use_emmvgwas_lambda) {
        throw std::runtime_error("--emmvgwas-vc-only requires --mv-compatible emmvgwas");
    }
    opts.pheno_file = cfg.pheno_file;
    opts.mpheno_list = cfg.mpheno_list;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.env_covar_col = cfg.env_covar_col;
    opts.scan_chrom = cfg.scan_chr;
    opts.scan_from_bp = cfg.scan_from_bp;
    opts.scan_to_bp = cfg.scan_to_bp;
    opts.extract_set = mv_load_snp_id_set(cfg.extract_file);
    opts.exclude_set = mv_load_snp_id_set(cfg.exclude_file);
    return opts;
}

MultivariateModel::MultivariateModel(const Options& opts) : options_(opts) {}

void MultivariateModel::prepare(const std::vector<std::string>& fids,
                                const std::vector<std::string>& iids,
                                int n_samples,
                                PlinkReader& plink_reader,
                                BgenReader& bgen_reader,
                                PgenReader& pgen_reader,
                                bool is_bgen, bool is_pgen) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;

    if (options_.pheno_file.empty() || options_.mpheno_list.empty()) {
        throw std::runtime_error("Multivariate analysis requires --pheno and --mpheno-list.");
    }

    loadPhenotype();
    loadCovariates();
}

void MultivariateModel::loadPhenotype() {
    std::cout << "Loading multivariate phenotypes...\n";
    PhenoReader pheno_reader;
    pheno_reader.loadMatrix(options_.pheno_file, options_.mpheno_list);
    size_t matched = pheno_reader.alignMatrix(fids_, iids_, Y_);
    std::cout << "Matched " << matched << " individuals with phenotype.\n";
    if (matched < 10) throw std::runtime_error("Too few samples.");
}

void MultivariateModel::loadCovariates() {
    Eigen::MatrixXf X_f = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
    X_ = X_f.cast<double>();
    Eigen::VectorXd env_raw;
    if (options_.use_gemma_gxe) {
        if (options_.env_covar_col <= 0) {
            throw std::runtime_error("--mv-method gemma-gxe requires --env-covar <phenotype-column>");
        }
        env_raw = mv_load_pheno_column_by_id(options_.pheno_file, options_.env_covar_col, fids_, iids_);
    }

    keep_.keep_indices.clear();
    keep_.keep_indices.reserve(n_samples_);
    for (int i = 0; i < n_samples_; ++i) {
        bool keep = true;
        for (int t = 0; t < Y_.cols(); ++t) {
            if (!std::isfinite(Y_(i, t))) {
                keep = false;
                break;
            }
        }
        if (keep) {
            for (int j = 0; j < X_.cols(); ++j) {
                if (!std::isfinite(X_(i, j))) {
                    keep = false;
                    break;
                }
            }
        }
        if (keep && options_.use_gemma_gxe && !std::isfinite(env_raw(i))) keep = false;
        if (keep) keep_.keep_indices.push_back(i);
    }
    keep_.n_keep = keep_.keep_indices.size();
    n_keep_ = keep_.n_keep;
    std::cout << "After filtering missing data, " << n_keep_ << " individuals remain.\n";
    if (n_keep_ < 10) throw std::runtime_error("Too few valid samples.");

    Eigen::MatrixXd Y_sub(n_keep_, Y_.cols());
    Eigen::MatrixXd X_sub(n_keep_, X_.cols());
    Eigen::VectorXd env_sub(n_keep_);
    for (int i = 0; i < n_keep_; ++i) {
        Y_sub.row(i) = Y_.row(keep_.keep_indices[i]);
        X_sub.row(i) = X_.row(keep_.keep_indices[i]);
        if (options_.use_gemma_gxe) env_sub(i) = env_raw(keep_.keep_indices[i]);
    }
    Y_ = std::move(Y_sub);
    X_ = std::move(X_sub);
    if (options_.use_gemma_gxe) env_ = std::move(env_sub);
}

void MultivariateModel::run() {
    if (options_.null_only) {
        runNullOnlyDiagnostic();
        return;
    }
    const std::string resolved_tier = options_.scan_tier.empty()
        ? (options_.use_hybrid ? "approx" : (options_.use_gemma_refit ? "exact" : (options_.use_emmv ? "semi" : "approx")))
        : options_.scan_tier;
    const std::string vc_strategy = resolved_tier == "exact" ? "per_snp_full_covariance_refit"
                                  : resolved_tier == "semi" ? "fixed_covariance_ratio_profile_residual"
                                                             : "fixed_null_covariance";
    const std::string test_label = options_.use_gemma_full || options_.use_gemma_exact
        ? "wald,lrt,score"
        : "wald";
    const std::string null_solver = options_.use_gemma_refit
        ? (options_.use_gemma_exact ? "gemma_exact_em" : "gemma_em")
        : (options_.use_emmv ? (options_.use_emmvgwas_lambda ? "emmvgwas_em" : "em_mv") : "fixed_null");
    const int gemma_default_alt_maxit = std::max(1, options_.gemma_emi / 10);
    const double gemma_default_alt_tol = options_.gemma_emp * 10.0;
    const std::string release_status = (options_.use_eb || options_.use_emmvgwas_lambda || options_.use_gemma_refit)
        ? "experimental"
        : "stable";
    {
        std::ofstream meta(options_.out_prefix + ".scan.meta.txt");
        meta << "Item\tValue\n"
             << "release_version\t" << COSMICSCANNER_VERSION << "\n"
             << "git_commit\t" << COSMICSCANNER_GIT_COMMIT << "\n"
             << "release_status\t" << release_status << "\n"
             << "analysis_model\tmultivariate_lmm\n"
             << "model\tmultivariate_lmm\n"
             << "requested_scan_tier\t" << (options_.scan_tier.empty() ? "legacy_method" : options_.scan_tier) << "\n"
             << "resolved_scan_tier\t" << resolved_tier << "\n"
             << "scan_tier\t" << resolved_tier << "\n"
             << "statistical_test\t" << test_label << "\n"
             << "null_solver\t" << null_solver << "\n"
             << "objective_type\tREML\n"
             << "null_optimizer\t" << (options_.use_gemma_refit ? "em_reml_correction+dense_reml_refinement" : null_solver) << "\n"
             << "compatibility_target\t" << (options_.use_gemma_compat ? "gemma_default_finite_em" : "none") << "\n"
             << "per_snp_max_iterations\t" << (options_.use_gemma_compat ? gemma_default_alt_maxit : 0) << "\n"
             << "per_snp_convergence_tolerance\t" << (options_.use_gemma_compat ? gemma_default_alt_tol : 0.0) << "\n"
             << "requested_method\t" << options_.requested_method << "\n"
             << "resolved_method\t" << options_.method_name << "\n"
             << "deprecated_alias\t" << (options_.deprecated_alias ? "true" : "false") << "\n"
             << "forced_wald_iterations\t" << options_.gemma_force_iterations << "\n"
             << "grm_scaling\t" << (options_.k_file.empty() ? "default_vanraden_centered_sum_2pq" : "external_k_file_no_rescaling") << "\n"
             << "external_k_file\t" << (options_.k_file.empty() ? "none" : options_.k_file) << "\n"
             << "loco\t" << (options_.no_loco ? "disabled_no_loco" : "default_loco_when_multiple_chromosomes") << "\n"
             << "fallback\tnone\n"
             << "multivariate_method\t" << options_.method_name << "\n"
             << "compatibility_mode\t" << options_.compatibility_mode << "\n"
             << "snp_vc_strategy\t" << vc_strategy << "\n";
    }
    std::cout << "Multivariate scan tier: " << resolved_tier << " (" << vc_strategy << ")\n";
    mv_validate_output_parent(options_.out_prefix);

    if (options_.use_gemma_refit) {
        runGemmaRefit();
        return;
    }

    std::cout << "Multivariate GWAS method: "
              << (options_.use_eb ? "EB-mvLMM random-SNP Wald" : (options_.use_emmv ? "EMmvGWAS semi-exact Wald" : "fixed-null Wald"))
              << "\n";

    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    grm_opts.thread_num = options_.n_threads;
    #endif

    std::cout << "Calculating Global GRM (Pass 1)...\n";
    GrmBuilder global_builder(grm_opts);
    if (is_pgen_) { pgen_reader_->reset(); global_builder.compute(*pgen_reader_); }
    else if (is_bgen_) { bgen_reader_->reset(); global_builder.compute(*bgen_reader_); }
    else { global_builder.compute(*plink_reader_); }

    const Eigen::MatrixXd& G_global_norm = global_builder.getGrm();
    double S_global = global_builder.getSum2pq();
    auto chroms_set = global_builder.getObservedChroms();
    const bool use_fixed_emmvgwas_vc = !options_.emmvgwas_G_file.empty();
    if (use_fixed_emmvgwas_vc && !options_.no_loco && chroms_set.size() > 1) {
        throw std::runtime_error("--mv-emmvgwas-G/--mv-emmvgwas-R implement EMmvGWAS semi_gwas semantics and require --no-loco");
    }

    std::string out_file = options_.out_prefix + ".assoc.mv.txt";
    {
        auto ofs = mv_open_output_file(out_file);
        eb_effective_tests_ = 0.0;
        eb_tested_markers_ = 0;
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\t" << (options_.use_eb ? "EB_WALD" : "WALD_STAT") << "\tP";
        int n_traits = Y_.cols();
        for (int t = 0; t < n_traits; ++t) {
            ofs << "\t" << (options_.use_eb ? "EB_BETA_" : "BETA_") << (t + 1)
                << "\t" << (options_.use_eb ? "EB_SE_" : "SE_") << (t + 1);
        }
        if (options_.output_vbeta) mv_write_vbeta_header(ofs, n_traits);
        if (options_.use_eb) ofs << "\tMARKER_LAMBDA\tEB_CONFIDENCE";
        ofs << "\n";
    }

    if (options_.no_loco || chroms_set.size() <= 1) {
        if (options_.no_loco) std::cout << "LOCO disabled. Using standard mvLMM...\n";
        else std::cout << "Only " << chroms_set.size() << " chromosome(s) found. Using standard mvLMM (no LOCO)...\n";

        Eigen::MatrixXd G_sub(n_keep_, n_keep_);
        for (int i = 0; i < n_keep_; ++i) {
            for (int j = 0; j < n_keep_; ++j) {
                G_sub(i, j) = G_global_norm(keep_.keep_indices[i], keep_.keep_indices[j]);
            }
        }

        MvEigenVce::Options vce_opts;
        vce_opts.use_projection = options_.use_covariate_projection;
        vce_opts.max_iter = options_.reml_maxit;
        vce_opts.tol = options_.reml_tol;
        vce_opts.solver = options_.use_aireml ? MvEigenVce::VceSolver::AIREML : MvEigenVce::VceSolver::EM;
        vce_opts.em_min_iter = options_.em_min_iter;
        vce_opts.ai_max_iter = options_.ai_max_iter;
        vce_opts.ai_tol = options_.ai_tol;
        vce_opts.use_emmvgwas_lambda = options_.use_emmvgwas_lambda && !use_fixed_emmvgwas_vc;
        MvEigenVce vce(vce_opts);
        vce.prepare(Y_, X_, G_sub);
        if (use_fixed_emmvgwas_vc) {
            Eigen::MatrixXd external_G = mv_read_numeric_csv_matrix(options_.emmvgwas_G_file);
            Eigen::MatrixXd external_R = mv_read_numeric_csv_matrix(options_.emmvgwas_R_file);
            vce.setEmmvgwasFixedComponents(external_G, external_R);
        } else {
            vce.runNullModel();
        }

        if (options_.emmvgwas_vc_only) {
            std::cout << "EMmvGWAS vc_estimation route: writing G/R/beta and skipping SNP scan.\n";
            writeEmmvgwasOutputs(vce, {});
            return;
        }

        MvLmmEngine::Options mv_opts;
        mv_opts.use_projection = options_.use_covariate_projection;
        mv_opts.em_mvgwas_test = options_.use_emmv;
        mv_opts.use_general_sigma_inverse = options_.use_emmvgwas_lambda;
        mv_opts.center_genotypes = !options_.use_emmvgwas_lambda;
        mv_opts.use_dosage_minus_one = options_.use_emmvgwas_lambda;
        mv_opts.eb_lmm = options_.use_eb;
        mv_opts.scan_chrom = options_.scan_chrom;
        mv_opts.scan_from_bp = options_.scan_from_bp;
        mv_opts.scan_to_bp = options_.scan_to_bp;
        mv_opts.extract_set = options_.extract_set;
        mv_opts.exclude_set = options_.exclude_set;
        #ifdef _OPENMP
        mv_opts.thread_num = options_.n_threads;
        #endif
        MvLmmEngine engine(mv_opts);
        engine.prepareFromVce(vce);

        std::vector<MvGwasResult> results;
        if (is_pgen_) { pgen_reader_->reset(); engine.scan(*pgen_reader_, results, keep_.keep_indices); }
        else if (is_bgen_) { bgen_reader_->reset(); engine.scan(*bgen_reader_, results, keep_.keep_indices); }
        else { engine.scan(*plink_reader_, results, keep_.keep_indices); }

        writeResults(out_file, results);
        if (options_.output_emmvgwas) writeEmmvgwasOutputs(vce, results);
    } else {
        std::cout << "Running Multivariate LOCO Scan on " << chroms_set.size() << " chromosomes...\n";
        std::vector<std::string> sorted_chroms(chroms_set.begin(), chroms_set.end());
        std::sort(sorted_chroms.begin(), sorted_chroms.end(), [](const std::string& a, const std::string& b) {
            try { return std::stoi(a) < std::stoi(b); } catch (...) { return a < b; }
        });

        MvLmmEngine::Options mv_opts;
        mv_opts.use_projection = options_.use_covariate_projection;
        mv_opts.em_mvgwas_test = options_.use_emmv;
        mv_opts.use_general_sigma_inverse = options_.use_emmvgwas_lambda;
        mv_opts.center_genotypes = !options_.use_emmvgwas_lambda;
        mv_opts.use_dosage_minus_one = options_.use_emmvgwas_lambda;
        mv_opts.eb_lmm = options_.use_eb;
        #ifdef _OPENMP
        mv_opts.thread_num = options_.n_threads;
        #endif

        for (const std::string& chrom : sorted_chroms) {
            if (!options_.scan_chrom.empty() && chrom != options_.scan_chrom) continue;
            std::cout << "Processing Chromosome " << chrom << "...\n";

            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);

            if (is_bgen_) { bgen_reader_->reset(); chr_builder.compute(*bgen_reader_); }
            else if (is_pgen_) { pgen_reader_->reset(); chr_builder.compute(*pgen_reader_); }
            else { chr_builder.compute(*plink_reader_); }

            double S_chr = chr_builder.getSum2pq();
            const auto& G_chr_norm = chr_builder.getGrm();

            Eigen::MatrixXd G_loco_sub(n_keep_, n_keep_);
            double S_loco = S_global - S_chr;

            if (S_chr > 0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                double w2 = S_chr / S_loco;
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < n_keep_; ++i) {
                    for (int j = 0; j < n_keep_; ++j) {
                        int gi = keep_.keep_indices[i];
                        int gj = keep_.keep_indices[j];
                        G_loco_sub(i, j) = G_global_norm(gi, gj) * w1 - G_chr_norm(gi, gj) * w2;
                    }
                }
            } else {
                for (int i = 0; i < n_keep_; ++i) {
                    for (int j = 0; j < n_keep_; ++j) {
                        G_loco_sub(i, j) = G_global_norm(keep_.keep_indices[i], keep_.keep_indices[j]);
                    }
                }
            }

            MvEigenVce::Options vce_opts;
            vce_opts.use_projection = options_.use_covariate_projection;
            vce_opts.max_iter = options_.reml_maxit;
            vce_opts.tol = options_.reml_tol;
            vce_opts.solver = options_.use_aireml ? MvEigenVce::VceSolver::AIREML : MvEigenVce::VceSolver::EM;
            vce_opts.em_min_iter = options_.em_min_iter;
            vce_opts.ai_max_iter = options_.ai_max_iter;
            vce_opts.ai_tol = options_.ai_tol;
            vce_opts.use_emmvgwas_lambda = options_.use_emmvgwas_lambda;
            MvEigenVce vce(vce_opts);
            vce.prepare(Y_, X_, G_loco_sub);
            vce.runNullModel();

            MvLmmEngine engine(mv_opts);
            engine.prepareFromVce(vce);

            std::vector<MvGwasResult> chrom_results;
            if (is_pgen_) { pgen_reader_->reset(); engine.scan(*pgen_reader_, chrom_results, keep_.keep_indices, chrom); }
            else if (is_bgen_) { bgen_reader_->reset(); engine.scan(*bgen_reader_, chrom_results, keep_.keep_indices, chrom); }
            else { engine.scan(*plink_reader_, chrom_results, keep_.keep_indices, chrom); }

            writeResults(out_file, chrom_results);
        }
    }
    if (options_.use_eb) {
        std::string eb_summary = options_.out_prefix + ".eb.mv.summary.txt";
        writeEbSummary(eb_summary);
        std::cout << "Multivariate EB summary saved to " << eb_summary << std::endl;
    }
    std::cout << "Results saved to " << out_file << std::endl;
}

void MultivariateModel::runNullOnlyDiagnostic() {
    if (options_.k_file.empty()) {
        throw std::runtime_error("--mv-null-only requires --mv-k-file for same-K diagnostics");
    }
    Eigen::MatrixXd K_full = mv_load_dense_matrix(options_.k_file, n_samples_);
    Eigen::MatrixXd K(n_keep_, n_keep_);
    for (int i = 0; i < n_keep_; ++i) {
        for (int j = 0; j < n_keep_; ++j) {
            K(i, j) = K_full(keep_.keep_indices[i], keep_.keep_indices[j]);
        }
    }
    K = 0.5 * (K + K.transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
    Eigen::VectorXd D = es.eigenvalues();
    Eigen::MatrixXd U = es.eigenvectors();
    Eigen::MatrixXd UtY = U.transpose() * Y_;
    Eigen::MatrixXd UtX = U.transpose() * X_;

    MvEigenVce::Options vce_opts;
    vce_opts.use_projection = false;
    vce_opts.max_iter = options_.reml_maxit;
    vce_opts.tol = options_.reml_tol;
    vce_opts.verbose = false;
    vce_opts.use_reml_correction = true;
    vce_opts.use_dense_reml_refinement = true;
    MvEigenVce vce(vce_opts);
    vce.prepare(Y_, X_, K);
    if (!options_.init_vg.empty() || !options_.init_ve.empty()) {
        if (options_.init_vg.empty() || options_.init_ve.empty()) {
            throw std::runtime_error("--mv-init-vg and --mv-init-ve must be provided together");
        }
        vce.setInitialComponents(mv_parse_square_matrix_arg(options_.init_vg, Y_.cols(), "--mv-init-vg"),
                                 mv_parse_square_matrix_arg(options_.init_ve, Y_.cols(), "--mv-init-ve"));
    }
    bool converged = vce.runNullModel();
    Eigen::MatrixXd Vg = vce.getVg();
    Eigen::MatrixXd Ve = vce.getVe();

    MvRemlEval cosmic_eval = mv_eval_reml_dense(UtY, UtX, D, Vg, Ve);

    std::ofstream out(options_.out_prefix + ".mv.null.diag.txt");
    out << "Item\tValue\n";
    out << "n_samples\t" << n_keep_ << "\n";
    out << "n_traits\t" << Y_.cols() << "\n";
    out << "trait_order\t";
    for (size_t i = 0; i < options_.mpheno_list.size(); ++i) {
        if (i) out << ',';
        out << options_.mpheno_list[i];
    }
    out << "\n";
    out << "fixed_effect_rows\t" << X_.rows() << "\n";
    out << "fixed_effect_cols\t" << X_.cols() << "\n";
    out << "k_file\t" << options_.k_file << "\n";
    out << "k_trace\t" << std::setprecision(12) << K.trace() << "\n";
    out << "k_mean_diag\t" << K.diagonal().mean() << "\n";
    out << "k_min_eigen\t" << D.minCoeff() << "\n";
    out << "k_max_eigen\t" << D.maxCoeff() << "\n";
    out << "phenotype_covariance\t";
    Eigen::MatrixXd Yc = Y_.rowwise() - Y_.colwise().mean();
    mv_write_matrix(out, (Yc.transpose() * Yc) / std::max(1, n_keep_ - 1));
    out << "\n";
    out << "cosmic_converged\t" << (converged ? "yes" : "no") << "\n";
    out << "cosmic_iterations\t" << vce.getLastIterations() << "\n";
    out << "cosmic_update_norm\t" << vce.getLastUpdateNorm() << "\n";
    out << "cosmic_internal_loglik\t" << vce.getLastLogLikelihood() << "\n";
    out << "cosmic_Vg\t"; mv_write_matrix(out, Vg); out << "\n";
    out << "cosmic_Ve\t"; mv_write_matrix(out, Ve); out << "\n";
    out << "cosmic_reml_loglik\t" << cosmic_eval.loglik << "\n";
    out << "cosmic_logdet_H\t" << cosmic_eval.logdet_h << "\n";
    out << "cosmic_logdet_XtHinvX\t" << cosmic_eval.logdet_xt_hinv_x << "\n";
    out << "cosmic_quad\t" << cosmic_eval.quad << "\n";
    out << "cosmic_fixed_beta\t";
    for (int i = 0; i < cosmic_eval.beta.size(); ++i) {
        if (i) out << ',';
        out << cosmic_eval.beta(i);
    }
    out << "\n";

    if (!options_.compare_vg.empty() || !options_.compare_ve.empty()) {
        if (options_.compare_vg.empty() || options_.compare_ve.empty()) {
            throw std::runtime_error("--mv-compare-vg and --mv-compare-ve must be provided together");
        }
        Eigen::MatrixXd Vg_cmp = mv_parse_square_matrix_arg(options_.compare_vg, Y_.cols(), "--mv-compare-vg");
        Eigen::MatrixXd Ve_cmp = mv_parse_square_matrix_arg(options_.compare_ve, Y_.cols(), "--mv-compare-ve");
        MvRemlEval cmp = mv_eval_reml_dense(UtY, UtX, D, Vg_cmp, Ve_cmp);
        out << "compare_Vg\t"; mv_write_matrix(out, Vg_cmp); out << "\n";
        out << "compare_Ve\t"; mv_write_matrix(out, Ve_cmp); out << "\n";
        out << "compare_reml_loglik\t" << cmp.loglik << "\n";
        out << "compare_logdet_H\t" << cmp.logdet_h << "\n";
        out << "compare_logdet_XtHinvX\t" << cmp.logdet_xt_hinv_x << "\n";
        out << "compare_quad\t" << cmp.quad << "\n";
        out << "compare_fixed_beta\t";
        for (int i = 0; i < cmp.beta.size(); ++i) {
            if (i) out << ',';
            out << cmp.beta(i);
        }
        out << "\n";
        out << "compare_minus_cosmic_loglik\t" << (cmp.loglik - cosmic_eval.loglik) << "\n";
    }
}


} // namespace cosmic
