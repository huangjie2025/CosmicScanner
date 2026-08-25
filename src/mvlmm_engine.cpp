#include "mvlmm_engine.h"
#include "pgen_reader.h"
#include "bgen_reader.h"
#include <iostream>
#include <cmath>
#include <omp.h>

namespace cosmic {

static double gamma_q(double a, double x) {
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

static double chisq_prob(double stat, int df) {
    if (stat <= 0) return 1.0;
    return gamma_q(df / 2.0, stat / 2.0);
}

static bool mv_variant_passes_filters(const PlinkBimInfo& var,
                                      const std::string& loco_chrom_filter,
                                      const MvLmmEngine::Options& options) {
    if (!loco_chrom_filter.empty() && var.chrom != loco_chrom_filter) return false;
    if (!options.scan_chrom.empty() && var.chrom != options.scan_chrom) return false;
    if (options.scan_from_bp >= 0 && var.bp_pos < options.scan_from_bp) return false;
    if (options.scan_to_bp >= 0 && var.bp_pos > options.scan_to_bp) return false;
    if (!options.extract_set.empty() && options.extract_set.find(var.id) == options.extract_set.end()) return false;
    if (!options.exclude_set.empty() && options.exclude_set.find(var.id) != options.exclude_set.end()) return false;
    return true;
}

static bool mv_has_scan_filters(const std::string& loco_chrom_filter,
                                const MvLmmEngine::Options& options) {
    return !loco_chrom_filter.empty() || !options.scan_chrom.empty() ||
           options.scan_from_bp >= 0 || options.scan_to_bp >= 0 ||
           !options.extract_set.empty() || !options.exclude_set.empty();
}

static Eigen::MatrixXd sanitize_spd_local(const Eigen::MatrixXd& A, double diag_floor = 1e-8) {
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

static bool invert_spd_jitter_local(const Eigen::MatrixXd& A, Eigen::MatrixXd& A_inv) {
    static const double jitters[] = {1e-8, 1e-6, 1e-4, 1e-2};
    Eigen::MatrixXd S = sanitize_spd_local(A, 1e-10);
    for (double jitter : jitters) {
        Eigen::MatrixXd T = S;
        T.diagonal().array() += jitter;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(T);
        if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
            A_inv = ldlt.solve(Eigen::MatrixXd::Identity(T.rows(), T.cols()));
            return A_inv.allFinite();
        }
    }
    return false;
}

struct MvEbFit {
    double lambda = 0.0;
    double confidence = 0.0;
    Eigen::VectorXd beta;
    Eigen::MatrixXd prior_cov;
    Eigen::MatrixXd post_cov;
    double wald = 0.0;
    double p = 1.0;
};

static double mv_logdet_ldlt_local(const Eigen::MatrixXd& M) {
    Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
    if (ldlt.info() != Eigen::Success) return std::numeric_limits<double>::infinity();
    double out = 0.0;
    for (int i = 0; i < ldlt.vectorD().size(); ++i) {
        double v = ldlt.vectorD()(i);
        if (!(v > 0.0) || !std::isfinite(v)) return std::numeric_limits<double>::infinity();
        out += std::log(v);
    }
    return out;
}

static MvEbFit fit_mv_eb_summary(const Eigen::VectorXd& beta_hat,
                                 const Eigen::MatrixXd& beta_cov,
                                 const Eigen::MatrixXd& prior_shape) {
    MvEbFit fit;
    const int q = beta_hat.size();
    fit.beta = Eigen::VectorXd::Zero(q);
    fit.prior_cov = Eigen::MatrixXd::Zero(q, q);
    fit.post_cov = Eigen::MatrixXd::Zero(q, q);
    if (q == 0 || !beta_hat.allFinite()) return fit;
    Eigen::MatrixXd V = sanitize_spd_local(beta_cov, 1e-10);
    Eigen::MatrixXd S = sanitize_spd_local(prior_shape, 1e-8);
    if (S.rows() != q || S.cols() != q) S = Eigen::MatrixXd::Identity(q, q);

    auto eval_ll = [&](double lambda) -> double {
        Eigen::MatrixXd C = V + std::max(0.0, lambda) * S;
        double ld = mv_logdet_ldlt_local(C);
        if (!std::isfinite(ld)) return -std::numeric_limits<double>::infinity();
        Eigen::LDLT<Eigen::MatrixXd> solver(C);
        Eigen::VectorXd x = solver.solve(beta_hat);
        if (!x.allFinite()) return -std::numeric_limits<double>::infinity();
        return -0.5 * (ld + beta_hat.dot(x));
    };

    double best_lambda = 0.0;
    double best_ll = eval_ll(0.0);
    const int grid_n = 41;
    int best_idx = -1;
    std::vector<double> grid(grid_n);
    for (int i = 0; i < grid_n; ++i) {
        double x = -12.0 + 24.0 * i / (double)(grid_n - 1);
        grid[i] = std::exp(x);
        double ll = eval_ll(grid[i]);
        if (ll > best_ll) { best_ll = ll; best_lambda = grid[i]; best_idx = i; }
    }
    if (best_idx >= 0 && best_ll > eval_ll(0.0) + 1e-10) {
        double left = (best_idx == 0) ? std::log(grid[0]) - 2.0 : std::log(grid[best_idx - 1]);
        double right = (best_idx == grid_n - 1) ? std::log(grid.back()) + 2.0 : std::log(grid[best_idx + 1]);
        const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
        double c = right - gr * (right - left);
        double d = left + gr * (right - left);
        double fc = eval_ll(std::exp(c));
        double fd = eval_ll(std::exp(d));
        for (int iter = 0; iter < 60 && std::abs(right - left) > 1e-7; ++iter) {
            if (fc < fd) { left = c; c = d; fc = fd; d = left + gr * (right - left); fd = eval_ll(std::exp(d)); }
            else { right = d; d = c; fd = fc; c = right - gr * (right - left); fc = eval_ll(std::exp(c)); }
        }
        best_lambda = std::exp(0.5 * (left + right));
    }

    if (!(best_lambda > 1e-10) || !std::isfinite(best_lambda)) return fit;
    Eigen::MatrixXd P = best_lambda * S;
    Eigen::MatrixXd C = V + P;
    Eigen::LDLT<Eigen::MatrixXd> c_solver(C);
    if (c_solver.info() != Eigen::Success) return fit;
    Eigen::MatrixXd K = P * c_solver.solve(Eigen::MatrixXd::Identity(q, q));
    fit.lambda = best_lambda;
    fit.prior_cov = P;
    fit.beta = K * beta_hat;
    fit.post_cov = P - K * P;
    fit.post_cov = 0.5 * (fit.post_cov + fit.post_cov.transpose());
    fit.confidence = std::max(0.0, std::min((double)q, K.trace()));
    Eigen::MatrixXd post_inv;
    if (invert_spd_jitter_local(fit.post_cov, post_inv)) {
        fit.wald = fit.beta.transpose() * post_inv * fit.beta;
        if (std::isfinite(fit.wald) && fit.wald >= 0.0) fit.p = chisq_prob(fit.wald, q);
    }
    return fit;
}

MvLmmEngine::MvLmmEngine(const Options& opts) : options(opts) {}

void MvLmmEngine::prepareFromVce(const MvEigenVce& vce) {
    Vg = vce.getVg();
    Ve = vce.getVe();
    D = vce.getEigenvalues();
    W_proj = vce.getProjectedW();
    Uty = vce.getUty();

    n_samples = D.size();
    n_traits = Vg.rows();

    UtX = vce.getUtX();
    n_covars = UtX.cols();
}

// Sparse Null Model Solver using PCG and MC-REML
void MvLmmEngine::processBlock(const Eigen::MatrixXd& S_block,
                               const std::vector<double>& means,
                               size_t global_start_idx,
                               const std::vector<size_t>& global_indices,
                               std::vector<MvGwasResult>& results,
                               const std::vector<struct PlinkBimInfo>& bim,
                               const Eigen::MatrixXd& H_XX_inv,
                               const Eigen::VectorXd& RHS_X,
                               const std::vector<Eigen::MatrixXd>& Sigma_inv,
                               const Eigen::MatrixXd& Cyy_null) {
    size_t current_block_size = S_block.cols();

    // S_block is already rotated by W_proj_t in scan()
    const Eigen::MatrixXd& UtS_block = S_block;

    int pd = n_covars * n_traits;

    // Initialize Accumulators for all SNPs in block
    std::vector<Eigen::MatrixXd> H_XS(current_block_size, Eigen::MatrixXd::Zero(pd, n_traits));
    std::vector<Eigen::MatrixXd> H_SS(current_block_size, Eigen::MatrixXd::Zero(n_traits, n_traits));
    std::vector<Eigen::VectorXd> RHS_S(current_block_size, Eigen::VectorXd::Zero(n_traits));

    // Optimized Loop: Sample-Major
    // Pre-calculate X_Si for sample i if possible, but Si is d x d.
    // Inner loop over k (SNPs) should benefit from cache locality of Si, xi, yi.

    #pragma omp parallel
    {
        // Thread-local accumulators
        std::vector<Eigen::MatrixXd> H_XS_loc(current_block_size, Eigen::MatrixXd::Zero(pd, n_traits));
        std::vector<Eigen::MatrixXd> H_SS_loc(current_block_size, Eigen::MatrixXd::Zero(n_traits, n_traits));
        std::vector<Eigen::VectorXd> RHS_S_loc(current_block_size, Eigen::VectorXd::Zero(n_traits));

        #pragma omp for schedule(static)
        for (int i = 0; i < n_samples; ++i) {
            const Eigen::MatrixXd& Si = Sigma_inv[i];
            Eigen::VectorXd yi = Uty.row(i).transpose(); // d x 1

            // Pre-compute common terms for this sample
            Eigen::VectorXd Si_yi = Si * yi; // d x 1

            for (size_t k = 0; k < current_block_size; ++k) {
                double si = UtS_block(i, k);

                // H_SS += si^2 * Si
                H_SS_loc[k] += (si * si) * Si;

                // RHS_S += si * Si * yi
                RHS_S_loc[k] += si * Si_yi;

                // H_XS
                // Block j: xi(j) * si * Si
                if (n_covars > 0) {
                    Eigen::RowVectorXd xi = UtX.row(i); // 1 x p
                    for (int j = 0; j < n_covars; ++j) {
                        H_XS_loc[k].block(j*n_traits, 0, n_traits, n_traits) += (xi(j) * si) * Si;
                    }
                }
            }
        }

        // Reduction
        #pragma omp critical
        {
            for (size_t k = 0; k < current_block_size; ++k) {
                if (n_covars > 0) {
                    H_XS[k] += H_XS_loc[k];
                }
                H_SS[k] += H_SS_loc[k];
                RHS_S[k] += RHS_S_loc[k];
            }
        }
    }

    // Finalize per SNP (independent)
    #pragma omp parallel for
    for (long long k = 0; k < (long long)current_block_size; ++k) {
        size_t global_idx = global_indices.empty() ? global_start_idx + k : global_indices[k];

        // Schur Complement
        Eigen::MatrixXd Schur;
        Eigen::VectorXd adjusted_RHS;
        if (n_covars > 0) {
            if (options.use_general_sigma_inverse && n_covars == 1) {
                Schur = H_SS[k] - H_XS[k] * H_XX_inv * H_XS[k];
                adjusted_RHS = RHS_S[k] - H_XS[k] * H_XX_inv * RHS_X;
            } else {
                Eigen::MatrixXd H_XS_T = H_XS[k].transpose();
                Schur = H_SS[k] - H_XS_T * H_XX_inv * H_XS[k];
                adjusted_RHS = RHS_S[k] - H_XS_T * H_XX_inv * RHS_X;
            }
        } else {
            Schur = H_SS[k];
            adjusted_RHS = RHS_S[k];
        }

        Eigen::MatrixXd Var_beta;
        Eigen::MatrixXd Schur_inv;
        if (options.use_general_sigma_inverse) {
            Schur_inv = Schur.inverse();
            if (!Schur_inv.allFinite()) continue;
        } else {
            Schur = sanitize_spd_local(Schur, 1e-8);
            if (!invert_spd_jitter_local(Schur, Schur_inv)) continue;
        }
        Eigen::VectorXd beta = Schur_inv * adjusted_RHS;
        if (!beta.allFinite()) continue;

        double wald_stat;
        if (options.em_mvgwas_test) {
            Eigen::MatrixXd Crr = Cyy_null - beta * adjusted_RHS.transpose();
            double df = std::max(1.0, (double)(n_samples - n_covars - 1));
            Eigen::MatrixXd ee = Crr / df;
            Var_beta = Schur_inv * ee;
            Eigen::MatrixXd Var_beta_inv;
            if (options.use_general_sigma_inverse) {
                Var_beta_inv = Var_beta.inverse();
                if (!Var_beta_inv.allFinite()) continue;
                wald_stat = beta.transpose() * Var_beta_inv * beta;
            } else {
                // To compute wald_stat, we need Var_beta^-1 = ee^-1 * Schur.
                Eigen::MatrixXd ee_inv;
                if (!invert_spd_jitter_local(ee, ee_inv)) continue;
                wald_stat = beta.transpose() * ee_inv * Schur * beta;
            }
        } else {
            Var_beta = Schur_inv;
            wald_stat = beta.transpose() * Schur * beta;
        }
        if (!std::isfinite(wald_stat) || wald_stat < 0.0) continue;

        double p_val = chisq_prob(wald_stat, n_traits);
        MvEbFit eb;
        if (options.eb_lmm) {
            eb = fit_mv_eb_summary(beta, Var_beta, Ve);
            beta = eb.beta;
            Var_beta = eb.post_cov;
            wald_stat = eb.wald;
            p_val = eb.p;
        }

        MvGwasResult res;
        if (!bim.empty()) {
            // bim is the block_bim passed from caller, so use local index k
            res.chrom = bim[k].chrom;
            res.snp_id = bim[k].id;
            res.pos = bim[k].bp_pos;
            res.a1 = bim[k].alt;
            res.a2 = bim[k].ref;
        }
        // For PGEN, we might fill these outside or pass them in differently.
        // Currently scan(PgenReader) handles metadata separately.
        // But to unify, we can pass metadata vector.

        res.freq = means[k] / 2.0;
        res.beta = beta;
        res.se = Var_beta;
        res.wald_stat = wald_stat;
        res.p_value = p_val;
        res.p_wald = p_val;
        if (options.eb_lmm) {
            res.eb_marker_lambda = eb.lambda;
            res.eb_confidence = eb.confidence;
            res.eb_prior_cov = eb.prior_cov;
            res.eb_post_cov = eb.post_cov;
        }

        results[global_idx] = res;
    }
}


// --- PGEN/BGEN Scan Support ---

void MvLmmEngine::scan(PgenReader& reader,
                       std::vector<MvGwasResult>& results,
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    size_t m = reader.getNumVariants();
    size_t n = keep_indices.size();

    int n_original = W_proj.rows();
    if (n != (size_t)n_original) {
        std::cerr << "MvLmmEngine::scan: n = " << n << ", n_original = " << n_original << std::endl;
        throw std::runtime_error("MvLmmEngine::scan: Sample size mismatch");
    }

    results.resize(m);

    std::vector<Eigen::MatrixXd> Sigma_inv_local;
    Sigma_inv_local.resize(n_samples);
    #pragma omp parallel for
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        if (options.use_general_sigma_inverse) {
            Sigma_inv_local[i] = Sigma.inverse();
        } else {
            Sigma_inv_local[i] = Sigma.llt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
        }
    }

    int pd = n_covars * n_traits;
    Eigen::MatrixXd H_XX = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd RHS_X = Eigen::VectorXd::Zero(pd);
    Eigen::MatrixXd H_YY = Eigen::MatrixXd::Zero(n_traits, n_traits);

    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Si = Sigma_inv_local[i];
        Eigen::RowVectorXd xi = UtX.row(i);
        Eigen::VectorXd yi = Uty.row(i).transpose();

        H_YY += Si * yi * yi.transpose();
        for (int j = 0; j < n_covars; ++j) {
            for (int k = 0; k < n_covars; ++k) {
                H_XX.block(j*n_traits, k*n_traits, n_traits, n_traits) += xi(j) * xi(k) * Si;
            }
            RHS_X.segment(j*n_traits, n_traits) += xi(j) * (Si * yi);
        }
    }

    Eigen::MatrixXd H_XX_inv;
    Eigen::MatrixXd Cyy_null = H_YY;
    if (pd > 0) {
        if (options.use_general_sigma_inverse) {
            H_XX_inv = H_XX.inverse();
        } else {
            H_XX_inv = H_XX.llt().solve(Eigen::MatrixXd::Identity(pd, pd));
        }
        Eigen::VectorXd B_hat = H_XX_inv * RHS_X;
        if (options.use_general_sigma_inverse && n_covars == 1) {
            Cyy_null = H_YY - B_hat * RHS_X.transpose();
        } else {
            Cyy_null = Eigen::MatrixXd::Zero(n_traits, n_traits);
            for (int i = 0; i < n_samples; ++i) {
                Eigen::VectorXd yi = Uty.row(i).transpose();
                Eigen::VectorXd y_hat = Eigen::VectorXd::Zero(n_traits);
                Eigen::RowVectorXd xi = UtX.row(i);
                for (int j = 0; j < n_covars; ++j) {
                    y_hat += xi(j) * B_hat.segment(j*n_traits, n_traits);
                }
                Eigen::VectorXd resid = yi - y_hat;
                Cyy_null += Sigma_inv_local[i] * resid * resid.transpose();
            }
        }
    } else {
        H_XX_inv = Eigen::MatrixXd::Zero(0, 0);
    }

    bool is_projected = (W_proj.cols() < W_proj.rows());
    std::cout << "Starting Multivariate Scan (PGEN" << (is_projected ? ", Projected" : "") << ") on " << m << " SNPs..." << std::endl;

    size_t processed = 0;
    size_t block_size = 2048;

    Eigen::MatrixXd W_proj_t = W_proj.transpose();

    reader.reset();

    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f;

    while (processed < m) {
        size_t current_block_size = std::min(block_size, m - processed);

        Eigen::MatrixXd S_block_raw(n_original, current_block_size);
        std::vector<double> means(current_block_size);
        std::vector<PlinkBimInfo> block_bim(current_block_size);

        size_t k;
        for (k = 0; k < current_block_size; ++k) {
            bool ok = reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
            if (!ok) break;

            block_bim[k].chrom = chrom;
            block_bim[k].id = rsid;
            block_bim[k].cm_pos = 0;
            block_bim[k].bp_pos = (long long)pos;
            block_bim[k].alt = alt;
            block_bim[k].ref = ref;

            double sum = 0;
            int count = 0;
            for (int i = 0; i < n_original; ++i) {
                double val = (double)dosage_f[keep_indices[i]];
                if (!std::isnan(val) && options.use_dosage_minus_one) val -= 1.0;
                S_block_raw(i, k) = val;
                if (!std::isnan(val)) {
                    sum += val;
                    count++;
                }
            }

            double mean = (count > 0) ? sum / count : 0.0;
            means[k] = mean;

            for (int i = 0; i < n_original; ++i) {
                if (std::isnan(S_block_raw(i, k))) S_block_raw(i, k) = 0.0;
                else if (options.center_genotypes) S_block_raw(i, k) -= mean;
            }
        }

        if (k == 0) break;
        current_block_size = k;

        if (current_block_size < block_size) {
            S_block_raw.conservativeResize(Eigen::NoChange, current_block_size);
            block_bim.resize(current_block_size);
            means.resize(current_block_size);
        }

        Eigen::MatrixXd S_rot = W_proj_t * S_block_raw;

        if (!mv_has_scan_filters(chrom_filter, options)) {
            processBlock(S_rot, means, processed, {}, results, block_bim, H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
        } else {
            std::vector<int> keep_cols;
            std::vector<size_t> global_indices;
            std::vector<double> filtered_means;
            std::vector<PlinkBimInfo> filtered_bim;
            for (size_t ci = 0; ci < current_block_size; ++ci) {
                if (mv_variant_passes_filters(block_bim[ci], chrom_filter, options)) {
                    keep_cols.push_back((int)ci);
                    global_indices.push_back(processed + ci);
                    filtered_means.push_back(means[ci]);
                    filtered_bim.push_back(block_bim[ci]);
                }
            }
            if (!keep_cols.empty()) {
                Eigen::MatrixXd S_filtered(S_rot.rows(), keep_cols.size());
                for (size_t ci = 0; ci < keep_cols.size(); ++ci) {
                    S_filtered.col(ci) = S_rot.col(keep_cols[ci]);
                }
                processBlock(S_filtered, filtered_means, 0, global_indices, results, filtered_bim,
                             H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
            }
        }

        processed += current_block_size;
        std::cout << "\rProcessed " << processed << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << std::endl;
}


void MvLmmEngine::scan(BgenReader& reader,
                       std::vector<MvGwasResult>& results,
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    size_t m = reader.getNumVariants();
    size_t n = keep_indices.size();

    int n_original = W_proj.rows();
    if (n != (size_t)n_original) {
        std::cerr << "MvLmmEngine::scan: n = " << n << ", n_original = " << n_original << std::endl;
        throw std::runtime_error("MvLmmEngine::scan: Sample size mismatch");
    }

    results.resize(m);

    std::vector<Eigen::MatrixXd> Sigma_inv_local;
    Sigma_inv_local.resize(n_samples);
    #pragma omp parallel for
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        if (options.use_general_sigma_inverse) {
            Sigma_inv_local[i] = Sigma.inverse();
        } else {
            Sigma_inv_local[i] = Sigma.llt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
        }
    }

    int pd = n_covars * n_traits;
    Eigen::MatrixXd H_XX = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd RHS_X = Eigen::VectorXd::Zero(pd);
    Eigen::MatrixXd H_YY = Eigen::MatrixXd::Zero(n_traits, n_traits);

    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Si = Sigma_inv_local[i];
        Eigen::RowVectorXd xi = UtX.row(i);
        Eigen::VectorXd yi = Uty.row(i).transpose();

        H_YY += Si * yi * yi.transpose();
        for (int j = 0; j < n_covars; ++j) {
            for (int k = 0; k < n_covars; ++k) {
                H_XX.block(j*n_traits, k*n_traits, n_traits, n_traits) += xi(j) * xi(k) * Si;
            }
            RHS_X.segment(j*n_traits, n_traits) += xi(j) * (Si * yi);
        }
    }

    Eigen::MatrixXd H_XX_inv;
    Eigen::MatrixXd Cyy_null = H_YY;
    if (pd > 0) {
        if (options.use_general_sigma_inverse) {
            H_XX_inv = H_XX.inverse();
        } else {
            H_XX_inv = H_XX.llt().solve(Eigen::MatrixXd::Identity(pd, pd));
        }
        Eigen::VectorXd B_hat = H_XX_inv * RHS_X;
        if (options.use_general_sigma_inverse && n_covars == 1) {
            Cyy_null = H_YY - B_hat * RHS_X.transpose();
        } else {
            Cyy_null = Eigen::MatrixXd::Zero(n_traits, n_traits);
            for (int i = 0; i < n_samples; ++i) {
                Eigen::VectorXd yi = Uty.row(i).transpose();
                Eigen::VectorXd y_hat = Eigen::VectorXd::Zero(n_traits);
                Eigen::RowVectorXd xi = UtX.row(i);
                for (int j = 0; j < n_covars; ++j) {
                    y_hat += xi(j) * B_hat.segment(j*n_traits, n_traits);
                }
                Eigen::VectorXd resid = yi - y_hat;
                Cyy_null += Sigma_inv_local[i] * resid * resid.transpose();
            }
        }
    } else {
        H_XX_inv = Eigen::MatrixXd::Zero(0, 0);
    }

    bool is_projected = (W_proj.cols() < W_proj.rows());
    std::cout << "Starting Multivariate Scan (BGEN" << (is_projected ? ", Projected" : "") << ") on " << m << " SNPs..." << std::endl;

    size_t processed = 0;
    size_t block_size = 2048;

    Eigen::MatrixXd W_proj_t = W_proj.transpose();

    reader.reset();

    std::string chrom, rsid, a1, a2;
    uint32_t pos;
    Eigen::VectorXf dosage_f;

    while (processed < m) {
        size_t current_block_size = std::min(block_size, m - processed);

        Eigen::MatrixXd S_block_raw(n_original, current_block_size);
        std::vector<double> means(current_block_size);
        std::vector<PlinkBimInfo> block_bim(current_block_size);

        size_t k;
        for (k = 0; k < current_block_size; ++k) {
            bool ok = reader.readVariant(chrom, rsid, pos, a1, a2, dosage_f);
            if (!ok) break;

            block_bim[k].chrom = chrom;
            block_bim[k].id = rsid;
            block_bim[k].cm_pos = 0;
            block_bim[k].bp_pos = (long long)pos;
            block_bim[k].alt = a1;
            block_bim[k].ref = a2;

            double sum = 0;
            int count = 0;
            for (int i = 0; i < n_original; ++i) {
                double val = (double)dosage_f[keep_indices[i]];
                if (!std::isnan(val) && options.use_dosage_minus_one) val -= 1.0;
                S_block_raw(i, k) = val;
                if (!std::isnan(val)) {
                    sum += val;
                    count++;
                }
            }

            double mean = (count > 0) ? sum / count : 0.0;
            means[k] = mean;

            for (int i = 0; i < n_original; ++i) {
                if (std::isnan(S_block_raw(i, k))) S_block_raw(i, k) = 0.0;
                else if (options.center_genotypes) S_block_raw(i, k) -= mean;
            }
        }

        if (k == 0) break;
        current_block_size = k;

        if (current_block_size < block_size) {
            S_block_raw.conservativeResize(Eigen::NoChange, current_block_size);
            block_bim.resize(current_block_size);
            means.resize(current_block_size);
        }

        Eigen::MatrixXd S_rot = W_proj_t * S_block_raw;

        if (!mv_has_scan_filters(chrom_filter, options)) {
            processBlock(S_rot, means, processed, {}, results, block_bim, H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
        } else {
            std::vector<int> keep_cols;
            std::vector<size_t> global_indices;
            std::vector<double> filtered_means;
            std::vector<PlinkBimInfo> filtered_bim;
            for (size_t ci = 0; ci < current_block_size; ++ci) {
                if (mv_variant_passes_filters(block_bim[ci], chrom_filter, options)) {
                    keep_cols.push_back((int)ci);
                    global_indices.push_back(processed + ci);
                    filtered_means.push_back(means[ci]);
                    filtered_bim.push_back(block_bim[ci]);
                }
            }
            if (!keep_cols.empty()) {
                Eigen::MatrixXd S_filtered(S_rot.rows(), keep_cols.size());
                for (size_t ci = 0; ci < keep_cols.size(); ++ci) {
                    S_filtered.col(ci) = S_rot.col(keep_cols[ci]);
                }
                processBlock(S_filtered, filtered_means, 0, global_indices, results, filtered_bim,
                             H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
            }
        }

        processed += current_block_size;
        std::cout << "\rProcessed " << processed << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << std::endl;
}



/*
void MvLmmEngine::scan_internal(...) {
    results.resize(m);

    if (use_sparse_grm) {
        // ...
    }
    // ...
}
*/

    // Orphaned code removed


void MvLmmEngine::scan(PlinkReader& reader,
                       std::vector<MvGwasResult>& results,
                       const std::vector<int>& keep_indices,
                       const std::string& chrom_filter) {
    size_t m = reader.getNumSnps();
    size_t n = keep_indices.size();

    int n_original = W_proj.rows();
    if (n != (size_t)n_original) {
        std::cerr << "MvLmmEngine::scan: n = " << n << ", n_original = " << n_original << std::endl;
        throw std::runtime_error("MvLmmEngine::scan: Sample size mismatch");
    }

    results.resize(m);
    const auto& bim = reader.getBimInfo();

    std::vector<Eigen::MatrixXd> Sigma_inv_local;

    // Pre-compute Sigma_inv
    Sigma_inv_local.resize(n_samples);
    #pragma omp parallel for
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        if (options.use_general_sigma_inverse) {
            Sigma_inv_local[i] = Sigma.inverse();
        } else {
            Sigma_inv_local[i] = Sigma.llt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
        }
    }

    // Pre-compute H_XX and H_XX_inv
    int pd = n_covars * n_traits;
    Eigen::MatrixXd H_XX = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd RHS_X = Eigen::VectorXd::Zero(pd);
    Eigen::MatrixXd H_YY = Eigen::MatrixXd::Zero(n_traits, n_traits);

    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd Si = Sigma_inv_local[i];
        Eigen::RowVectorXd xi = UtX.row(i);
        Eigen::VectorXd yi = Uty.row(i).transpose();

        H_YY += Si * yi * yi.transpose();
        for (int j = 0; j < n_covars; ++j) {
            for (int k = 0; k < n_covars; ++k) {
                H_XX.block(j*n_traits, k*n_traits, n_traits, n_traits) += xi(j) * xi(k) * Si;
            }
            RHS_X.segment(j*n_traits, n_traits) += xi(j) * (Si * yi);
        }
    }

    Eigen::MatrixXd H_XX_inv;
    Eigen::MatrixXd Cyy_null = H_YY;
    if (pd > 0) {
        if (options.use_general_sigma_inverse) {
            H_XX_inv = H_XX.inverse();
        } else {
            H_XX_inv = H_XX.llt().solve(Eigen::MatrixXd::Identity(pd, pd));
        }
        Eigen::VectorXd B_hat = H_XX_inv * RHS_X; // pd x 1
        if (options.use_general_sigma_inverse && n_covars == 1) {
            Cyy_null = H_YY - B_hat * RHS_X.transpose();
        } else {
            Cyy_null = Eigen::MatrixXd::Zero(n_traits, n_traits);
            for (int i = 0; i < n_samples; ++i) {
                Eigen::VectorXd yi = Uty.row(i).transpose();
                Eigen::VectorXd y_hat = Eigen::VectorXd::Zero(n_traits);
                Eigen::RowVectorXd xi = UtX.row(i);
                for (int j = 0; j < n_covars; ++j) {
                    y_hat += xi(j) * B_hat.segment(j*n_traits, n_traits);
                }
                Eigen::VectorXd resid = yi - y_hat;
                Cyy_null += Sigma_inv_local[i] * resid * resid.transpose();
            }
        }
    } else {
        H_XX_inv = Eigen::MatrixXd::Zero(0, 0);
    }

    // Scan Loop
    bool is_projected = (W_proj.cols() < W_proj.rows());
    std::cout << "Starting Multivariate Scan (PLINK" << (is_projected ? ", Projected" : "") << ") on " << m << " SNPs..." << std::endl;

    size_t processed = 0;
    size_t block_size = 2048; // Optimized block size

    // U transpose for rotation
    Eigen::MatrixXd W_proj_t = W_proj.transpose(); // W_proj is n_original x n_samples

    while (processed < m) {
        size_t current_block_size = std::min(block_size, m - processed);
        size_t start_idx = processed;
        size_t end_idx = processed + current_block_size;

        if (mv_has_scan_filters(chrom_filter, options)) {
            bool block_has_target = false;
            for (size_t k = 0; k < current_block_size; ++k) {
                if (mv_variant_passes_filters(bim[start_idx + k], chrom_filter, options)) {
                    block_has_target = true;
                    break;
                }
            }
            if (!block_has_target) {
                processed += current_block_size;
                std::cout << "\rProcessed " << processed << "/" << m << " SNPs..." << std::flush;
                continue;
            }
        }

        // 1. Read Block (Raw Doubles)
        Eigen::MatrixXd S_block_raw;
        reader.getSnpBlock(start_idx, end_idx, S_block_raw);

        // 2. Impute and Center
        std::vector<double> means(current_block_size);
        std::vector<struct PlinkBimInfo> block_bim(current_block_size);

        #pragma omp parallel for
        for (long long k = 0; k < (long long)current_block_size; ++k) {
            size_t snp_idx = start_idx + k;
            block_bim[k] = bim[snp_idx];

            double sum = 0;
            int count = 0;
            for(int i=0; i<n_original; ++i) {
                double val = S_block_raw(i, k);
                if (!std::isnan(val) && options.use_dosage_minus_one) val -= 1.0;
                S_block_raw(i, k) = val;
                if (!std::isnan(val)) {
                    sum += val;
                    count++;
                }
            }

            double mean = (count > 0) ? sum / count : 0.0;
            means[k] = mean;

            for(int i=0; i<n_original; ++i) {
                if (std::isnan(S_block_raw(i, k))) S_block_raw(i, k) = 0.0;
                else if (options.center_genotypes) S_block_raw(i, k) -= mean;
            }
        }

        // 3. Rotate Block: S_rot = W_proj' * S_block
        // This is the heavy lifting: Matrix Multiplication
        Eigen::MatrixXd S_rot = W_proj_t * S_block_raw;

        if (!mv_has_scan_filters(chrom_filter, options)) {
            processBlock(S_rot, means, processed, {}, results, block_bim, H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
        } else {
            std::vector<int> keep_cols;
            std::vector<size_t> global_indices;
            std::vector<double> filtered_means;
            std::vector<struct PlinkBimInfo> filtered_bim;
            for (size_t k = 0; k < current_block_size; ++k) {
                if (mv_variant_passes_filters(block_bim[k], chrom_filter, options)) {
                    keep_cols.push_back((int)k);
                    global_indices.push_back(start_idx + k);
                    filtered_means.push_back(means[k]);
                    filtered_bim.push_back(block_bim[k]);
                }
            }
            if (!keep_cols.empty()) {
                Eigen::MatrixXd S_filtered(S_rot.rows(), keep_cols.size());
                for (size_t k = 0; k < keep_cols.size(); ++k) {
                    S_filtered.col(k) = S_rot.col(keep_cols[k]);
                }
                processBlock(S_filtered, filtered_means, 0, global_indices, results, filtered_bim,
                             H_XX_inv, RHS_X, Sigma_inv_local, Cyy_null);
            }
        }

        processed += current_block_size;
        std::cout << "\rProcessed " << processed << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << std::endl;
}

} // namespace cosmic
