#include "lmm_engine.h"
#include "vce.h"
#include "matrix_free_pcg.h"
#include "qc_utils.h"
#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <map>
#include <limits>
#include "rng.h"

namespace cosmic {

// --- Chi-square functions (copied for self-containment) ---
static double gamma_q(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) {
        // Series representation
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
        // Continued fraction representation
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

static void fill_gemma_lmm_tests(GwasResult& res,
                                 double score,
                                 double var_score,
                                 double null_quad,
                                 double n_eff) {
    if (!(var_score > 1e-12) || !(null_quad > 0.0) || !std::isfinite(score) ||
        !std::isfinite(var_score) || !std::isfinite(null_quad)) {
        return;
    }
    const double score_chisq = score * score / var_score;
    const double alt_quad = std::max(null_quad - score_chisq, std::numeric_limits<double>::min());
    const double lrt = std::max(0.0, n_eff * std::log(null_quad / alt_quad));

    res.score_stat = score_chisq;
    res.wald_stat = score_chisq;
    res.lrt_stat = lrt;
    res.p_score = chisq_prob(score_chisq, 1);
    res.p_wald = chisq_prob(score_chisq, 1);
    res.p_lrt = chisq_prob(lrt, 1);
}

static void apply_fixed_ratio_profile(GwasResult& res,
                                      double score,
                                      double var_score,
                                      double null_quad,
                                      double n_eff) {
    if (!(var_score > 1e-12) || !(null_quad > 0.0) || !(n_eff > 1.0)) return;
    const double reduction = score * score / var_score;
    double alt_quad = null_quad - reduction;
    const double tolerance = 1e-10 * null_quad;
    if (alt_quad <= 0.0) {
        if (alt_quad < -tolerance) return;
        alt_quad = std::numeric_limits<double>::min();
    }
    const double profile_scale = alt_quad / std::max(1.0, n_eff - 1.0);
    const double var_beta = profile_scale / var_score;
    if (!(var_beta > 0.0) || !std::isfinite(var_beta)) return;

    res.beta = score / var_score;
    res.se = std::sqrt(var_beta);
    res.wald_stat = res.beta * res.beta / var_beta;
    res.p_wald = chisq_prob(res.wald_stat, 1);
    res.lrt_stat = std::max(0.0, n_eff * std::log(null_quad / alt_quad));
    res.p_lrt = chisq_prob(res.lrt_stat, 1);
    res.p_value = res.p_wald;
    res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
}

struct ProfileRemlFit {
    bool ok = false;
    double logl = -std::numeric_limits<double>::infinity();
    double beta_last = 0.0;
    double se_last = std::numeric_limits<double>::quiet_NaN();
};

static Eigen::MatrixXd append_design_column(const Eigen::MatrixXd& X,
                                            const Eigen::VectorXd& x) {
    Eigen::MatrixXd out(X.rows(), X.cols() + 1);
    if (X.cols() > 0) out.leftCols(X.cols()) = X;
    out.col(X.cols()) = x;
    return out;
}

static ProfileRemlFit profile_lmm_fit(const Eigen::VectorXd& D,
                                      const Eigen::VectorXd& y,
                                      const Eigen::MatrixXd& X,
                                      bool restricted) {
    ProfileRemlFit best;
    const int n = (int)y.size();
    const int p = (int)X.cols();
    const double n_eff = std::max(1.0, restricted ? (double)(n - p) : (double)n);
    if (n <= 0 || D.size() != n || X.rows() != n || n <= p) return best;

    auto eval = [&](double log_eta, Eigen::VectorXd* beta_out,
                    Eigen::MatrixXd* xtx_inv_out, double* sigma2_out) -> double {
        const double eta = std::exp(std::min(20.0, std::max(-20.0, log_eta)));
        Eigen::VectorXd H = eta * D.array() + 1.0;
        for (int i = 0; i < H.size(); ++i) if (!(H[i] > 1e-12)) H[i] = 1e-12;
        Eigen::VectorXd Hinv = H.cwiseInverse();
        Eigen::VectorXd resid = y;
        Eigen::MatrixXd XtHinvX;
        Eigen::MatrixXd XtHinvX_inv;
        Eigen::VectorXd beta;
        double log_det_X = 0.0;

        if (p > 0) {
            Eigen::MatrixXd Xw = X.array().colwise() * Hinv.array().sqrt();
            Eigen::VectorXd yw = y.array() * Hinv.array().sqrt();
            XtHinvX = Xw.transpose() * Xw;
            XtHinvX.diagonal().array() += 1e-10;
            Eigen::LDLT<Eigen::MatrixXd> solver(XtHinvX);
            if (solver.info() == Eigen::Success) {
                XtHinvX_inv = solver.solve(Eigen::MatrixXd::Identity(p, p));
                beta = solver.solve(Xw.transpose() * yw);
                log_det_X = solver.vectorD().array().abs().cwiseMax(1e-300).log().sum();
            } else {
                Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(XtHinvX);
                beta = cod.solve(Xw.transpose() * yw);
                XtHinvX_inv = cod.pseudoInverse();
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(XtHinvX);
                const auto& s = svd.singularValues();
                for (int i = 0; i < s.size(); ++i) log_det_X += std::log(std::max(1e-300, s[i]));
            }
            resid -= X * beta;
        }

        const double quad = (resid.array().square() * Hinv.array()).sum();
        if (!(quad > 0.0) || !std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
        const double log_det_H = H.array().log().sum();
        const double sigma2 = quad / n_eff;
        const double logl = -0.5 * (log_det_H + (restricted ? log_det_X : 0.0) +
                                    n_eff * std::log(sigma2) + n_eff);
        if (beta_out) *beta_out = beta;
        if (xtx_inv_out) *xtx_inv_out = XtHinvX_inv;
        if (sigma2_out) *sigma2_out = sigma2;
        return logl;
    };

    const double lo = -12.0;
    const double hi = 12.0;
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = lo, b = hi;
    double c = b - gr * (b - a);
    double d = a + gr * (b - a);
    double fc = eval(c, nullptr, nullptr, nullptr);
    double fd = eval(d, nullptr, nullptr, nullptr);
    for (int iter = 0; iter < 80; ++iter) {
        if (fc < fd) {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = eval(d, nullptr, nullptr, nullptr);
        } else {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = eval(c, nullptr, nullptr, nullptr);
        }
    }

    const double opt = 0.5 * (a + b);
    Eigen::VectorXd beta;
    Eigen::MatrixXd xtx_inv;
    double sigma2 = 0.0;
    best.logl = eval(opt, &beta, &xtx_inv, &sigma2);
    best.ok = std::isfinite(best.logl);
    if (best.ok && p > 0 && beta.size() == p && xtx_inv.rows() == p) {
        best.beta_last = beta[p - 1];
        const double var = sigma2 * xtx_inv(p - 1, p - 1);
        best.se_last = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
    }
    return best;
}

static void apply_gemma_refit(GwasResult& res,
                              const ProfileRemlFit& null_ml,
                              const ProfileRemlFit& alt_ml,
                              const ProfileRemlFit& alt_reml) {
    if (alt_reml.ok && alt_reml.se_last > 0.0) {
        res.beta = alt_reml.beta_last;
        res.se = alt_reml.se_last;
        res.wald_stat = (res.beta / res.se) * (res.beta / res.se);
        res.p_wald = chisq_prob(res.wald_stat, 1);
        res.p_value = res.p_wald;
        res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
    }
    if (null_ml.ok && alt_ml.ok) {
        res.lrt_stat = std::max(0.0, 2.0 * (alt_ml.logl - null_ml.logl));
        res.p_lrt = chisq_prob(res.lrt_stat, 1);
    }
}
// --------------------------------------------------------

void LmmEngine::prepareBiobank(const Eigen::VectorXd& y,
                               const Eigen::MatrixXf& X,
                               double s2g, double s2e) {

    // Check dimensions
    if (y.size() != X.rows()) throw std::runtime_error("Dimension mismatch between y and X in prepareBiobank");

    // Store VCs
    sigma2_g = s2g;
    sigma2_e = s2e;

    // 1. Project out Covariates from y
    // M = I - X(X'X)^-1 X'
    // y_proj = M * y
    // y_raw = y_proj

    // Convert X to double for precision
    Eigen::MatrixXd X_d = X.cast<double>();
    int n = X_d.rows();
    int p = X_d.cols();

    // QR Decomposition for stable projection
    // X = Q * R
    // P_X = Q * Q' (if Q is thin, n x p)
    // M = I - Q * Q'
    // y_proj = y - Q * (Q' * y)

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(X_d);
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n, p);

    y_raw = y - Q * (Q.transpose() * y);

    // Store Q for later use in PCG (to project G*v)
    // We can store Q in LmmEngine or pass it to scan?
    // Let's store it as X_raw (repurposed) or a new member.
    // Actually, X_raw is float. Let's store Q in X_raw if p is small?
    // Or add a member `Eigen::MatrixXd Q_proj`.
    // Let's add it to class (requires header change) or just use X_raw (which is member).
    // X_raw is MatrixXf. Q is double.
    // Let's resize X_raw to store Q (as float).
    X_raw = Q.cast<float>();

    std::cout << "Biobank Mode Prepared:" << std::endl;
    std::cout << "  Projected covariates (p=" << p << ") out of phenotype." << std::endl;
    std::cout << "  Sigma2_g: " << sigma2_g << ", Sigma2_e: " << sigma2_e << std::endl;
}

double LmmEngine::estimateTraceInverse(MatrixFreePCG& pcg, int n, int k) {
    std::cout << "  Estimating Trace(V^-1) using " << k << " random vectors..." << std::endl;

    RNG rng; // Uses high-res clock seed
    double sum_trace = 0.0;

    // Hutchinson's Trace Estimator
    // Trace(A) = E[z' A z] where z is random vector with E[zz'] = I
    // Usually z is Rademacher (+1/-1) or Normal(0,1)

    for (int i = 0; i < k; ++i) {
        Eigen::VectorXd z(n);
        // Rademacher is slightly better for trace estimation (lower variance)
        for(int j=0; j<n; ++j) {
            z(j) = (rng.normal() > 0) ? 1.0 : -1.0;
        }

        // Solve V * w = z  => w = V^-1 * z
        Eigen::VectorXd w = pcg.solve(z, 1e-3, 50); // Lower precision is fine for trace

        // z' * w = z' * V^-1 * z
        double val = z.dot(w);
        sum_trace += val;

        std::cout << "\r    Iter " << (i+1) << "/" << k << ": Trace=" << (sum_trace/(i+1)) << std::flush;
    }
    std::cout << std::endl;

    return sum_trace / k;
}

void LmmEngine::scanBiobank(PlinkReader& reader,
                            std::vector<GwasResult>& results,
                            const std::vector<int>& keep_indices,
                            const std::string& chrom_filter) {
    // 1. Prepare PCG
    int n = (int)y_raw.size();
    size_t m = reader.getNumSnps();

    if (keep_indices.size() != (size_t)n) throw std::runtime_error("Sample size mismatch in scanBiobank");

    std::cout << "Starting Biobank Scan (Matrix-Free PCG) on " << m << " variants..." << std::endl;
    std::cout << "  Sigma2_g: " << sigma2_g << ", Sigma2_e: " << sigma2_e << std::endl;

    // Define MultOp: v_out = (sigma2_g * G + sigma2_e * I) * v_in
    double inv_M = 1.0 / (double)m;
    auto mult_op = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        Eigen::VectorXd u(m);
        reader.multiply_Zt_v(v_in, u); // u = Z' * v_in

        Eigen::VectorXd w(n);
        reader.multiply_Z_v(u, w);     // w = Z * u = Z * Z' * v_in

        v_out = (sigma2_g * inv_M) * w + sigma2_e * v_in;
    };

    // Preconditioner: Diagonal approximation
    double scalar_precond = 1.0 / (sigma2_g + sigma2_e);
    auto precond_op = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        v_out = v_in * scalar_precond;
    };

    MatrixFreePCG pcg(n);
    pcg.setMultOp(mult_op);
    pcg.setPrecondOp(precond_op);

    // 2. Solve V * x = y
    std::cout << "  Solving V * x = y..." << std::endl;
    Eigen::VectorXd x = pcg.solve(y_raw, 1e-5, 200);

    // 3. Estimate Trace(V^-1)
    double tr_V_inv = estimateTraceInverse(pcg, n, 10);
    double tau = tr_V_inv / n; // Average diagonal element approximation
    std::cout << "  Estimated Trace(V^-1) = " << tr_V_inv << ", tau = " << tau << std::endl;

    // 4. Scan
    std::cout << "  Scanning variants..." << std::endl;
    results.resize(m);

    Eigen::VectorXd scores(m);
    reader.multiply_Zt_v(x, scores);

    const auto& bim = reader.getBimInfo();
    for(size_t i=0; i<m; ++i) {
        if (!options.extract_set.empty() && options.extract_set.find(bim[i].id) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(bim[i].id) != options.exclude_set.end()) continue;
        results[i].chrom = bim[i].chrom;
        results[i].snp_id = bim[i].id;
        results[i].pos = bim[i].bp_pos;
        results[i].a1 = bim[i].alt;
        results[i].a2 = bim[i].ref;
        results[i].n_obs = n;

        // P-value Calculation
        // Score S = g' * V^-1 * y = g' * x
        // Var(S) = g' * V^-1 * g
        // Approx: Var(S) ~ tau * (g'g)
        // Since g is standardized, g'g ~ N (if no missing)
        // Or strictly: g'g = sum( (g_i - mu)^2 / sigma^2 ) = N-1 approx N.
        // So Var(S) ~ tau * N.

        double score = scores[i];
        double var_score = tau * n;

        // Beta approx: beta ~ S / Var(S) ? No.
        // Score Test Statistic T = S^2 / Var(S) ~ Chi2(1)
        // Beta is not directly from Score Test, but usually Beta ~ S / (g' V^-1 g)

        double se = std::sqrt(var_score);
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);

        // Estimate Beta
        // beta = S / Var(S) (This is the One-step approximation)
        double beta = (var_score > 1e-9) ? (score / var_score) : 0.0;
        // SE of beta = 1 / sqrt(Var(S))
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;

        results[i].beta = beta;
        results[i].se = beta_se;
        results[i].p_value = pval;
        results[i].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
    }
}

void LmmEngine::scanBiobank(PgenReader& reader,
                            std::vector<GwasResult>& results,
                            const std::vector<int>& keep_indices,
                            const std::string& chrom_filter) {
    int n = (int)y_raw.size();
    size_t m = reader.getNumVariants();

    if (keep_indices.size() != (size_t)n) throw std::runtime_error("Sample size mismatch in scanBiobank");

    std::cout << "Starting Biobank Scan (Matrix-Free PCG) on PGEN (" << m << " variants)..." << std::endl;
    const int n_total = static_cast<int>(reader.getNumSamples());
    for (int idx : keep_indices) {
        if (idx < 0 || idx >= n_total) {
            throw std::runtime_error("PGEN biobank keep index is out of range.");
        }
    }

    auto expand_to_full_samples = [&](const Eigen::VectorXd& subset) {
        Eigen::VectorXd full = Eigen::VectorXd::Zero(n_total);
        for (int i = 0; i < n; ++i) full(keep_indices[i]) = subset(i);
        return full;
    };

    auto subset_full_samples = [&](const Eigen::VectorXd& full) {
        if (full.size() != n_total) throw std::runtime_error("PGEN biobank full-sample vector size mismatch.");
        Eigen::VectorXd subset(n);
        for (int i = 0; i < n; ++i) subset(i) = full(keep_indices[i]);
        return subset;
    };

    auto multiply_Zt_keep = [&](const Eigen::VectorXd& v_subset, Eigen::VectorXd& out) {
        Eigen::VectorXd v_full = expand_to_full_samples(v_subset);
        reader.multiply_Zt_v(v_full, out, options.thread_num);
    };

    auto multiply_Z_keep = [&](const Eigen::VectorXd& marker_weights, Eigen::VectorXd& out_subset) {
        Eigen::VectorXd out_full;
        reader.multiply_Z_v(marker_weights, out_full, options.thread_num);
        out_subset = subset_full_samples(out_full);
    };

    // 2. Identify Unique Chromosomes and their variant indices
    std::vector<std::string> chroms;
    std::map<std::string, std::vector<size_t>> chrom_map;

    const auto& vars = reader.getVariants();
    for(size_t i=0; i<m; ++i) {
        const std::string& chr = vars[i].chrom;
        if (chrom_map.find(chr) == chrom_map.end()) {
            chroms.push_back(chr);
        }
        chrom_map[chr].push_back(i);
    }

    std::cout << "  Found " << chroms.size() << " chromosomes." << std::endl;

    // Cache Data for LOCO
    if (!reader.isCached()) {
        reader.cacheAllData();
    }

    // 3. Initialize Global Solution (Warm Start)
    // First solve global model to get a good starting point
    // Or solve first chromosome from scratch?
    // Let's solve global first to get V^-1 y estimate.

    std::cout << "  Solving Global Model (for Warm Start)..." << std::endl;

    // Retrieve Q from X_raw
    Eigen::MatrixXd Q_proj = X_raw.cast<double>();

    // Define Global MultOp
    double inv_M = 1.0 / (double)m;
    auto mult_op_global = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        Eigen::VectorXd u(m);
        multiply_Zt_keep(v_in, u);
        Eigen::VectorXd w(n);
        multiply_Z_keep(u, w);

        // Project w onto null space of X
        // w_proj = M * w = w - Q * (Q' * w)
        Eigen::VectorXd w_proj = w - Q_proj * (Q_proj.transpose() * w);

        v_out = (sigma2_g * inv_M) * w_proj + sigma2_e * v_in;
    };

    double scalar_precond = 1.0 / (sigma2_g + sigma2_e);
    auto precond_op = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        v_out = v_in * scalar_precond;
    };

    MatrixFreePCG pcg(n);
    pcg.setMultOp(mult_op_global);
    pcg.setPrecondOp(precond_op);

    Eigen::VectorXd x_global = pcg.solve(y_raw, 1e-4, 100); // Coarse solve

    // 4. LOCO Loop
    results.resize(m);

    for (const auto& chr : chroms) {
        if (!chrom_filter.empty() && chr != chrom_filter) continue;

        std::cout << "  Processing Chromosome " << chr << " (LOCO)..." << std::endl;

        // Define LOCO MultOp: G_loco = G_total - G_chr
        // G_loco * v = (1/M_loco) * Z_loco * Z_loco' * v
        // BUT standard practice uses 1/M_total scaling for all chromosomes?
        // Or 1/M_loco?
        // Usually we assume h2 is fixed. sigma2_g is for total genetic variance.
        // If we remove one chromosome, the genetic variance decreases slightly.
        // But for LOCO, we usually keep sigma2_g fixed or re-estimate.
        // Re-estimating is too slow.
        // We assume sigma2_g is distributed across M SNPs.
        // So G_loco should be scaled by 1/M_loco?
        // Let's stick to 1/M_total for simplicity, or 1/M_loco to keep scale of G=1.
        // Standard: G_loco = (1/M_loco) * sum_{j not in c} Z_j Z_j'

        const std::vector<size_t>& indices = chrom_map[chr];
        size_t m_chr = indices.size();
        size_t m_loco = m - m_chr;
        if (m_loco == 0) continue; // Only 1 chrom?

        double inv_M_loco = 1.0 / (double)m_loco;

        // Masked MultOp
        auto mult_op_loco = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
            Eigen::VectorXd u(m);
            multiply_Zt_keep(v_in, u);

            // Mask out current chromosome
            for (size_t idx : indices) u(idx) = 0.0;

            Eigen::VectorXd w(n);
            multiply_Z_keep(u, w);

            // Re-scale?
            // The sum inside w is sum_{j not in c}.
            // So we multiply by inv_M_loco.
            v_out = (sigma2_g * inv_M_loco) * w + sigma2_e * v_in;
        };

        pcg.setMultOp(mult_op_loco);

        // Solve with Warm Start
        Eigen::VectorXd x_loco = pcg.solve(y_raw, 1e-5, 200, x_global);

        // Estimate Trace for LOCO?
        // Reuse global tau?
        // Global tau = Tr(V^-1)/N.
        // V_loco is similar to V_global.
        // Re-estimating trace for every chromosome is slow (10 more PCG runs).
        // Let's use global trace estimate but scaled?
        // Actually, let's just use the global tau for approximation.
        // Or run a quick trace estimate (k=2?)

        // Quick Trace Estimate (k=5)
        // double tr_V_inv = estimateTraceInverse(pcg, n, 5);
        // double tau = tr_V_inv / n;

        // For now, use global trace estimate from global model (if we computed it).
        // We didn't compute it yet.
        // Let's compute global trace once.
        static double global_tau = 0.0;
        if (global_tau == 0.0) {
            pcg.setMultOp(mult_op_global); // Switch back briefly
            double tr = estimateTraceInverse(pcg, n, 5);
            global_tau = tr / n;
            pcg.setMultOp(mult_op_loco); // Switch back
        }

        // Scan variants on this chromosome
        // Score = Z_j' * x_loco
        // We only need scores for current chromosome.
        // But multiply_Zt_v computes all.
        // We can just extract relevant ones.

        Eigen::VectorXd scores(m);
        multiply_Zt_keep(x_loco, scores);

        for(size_t idx : indices) {
            if (!options.extract_set.empty() && options.extract_set.find(vars[idx].id) == options.extract_set.end()) continue;
            if (!options.exclude_set.empty() && options.exclude_set.find(vars[idx].id) != options.exclude_set.end()) continue;
            results[idx].chrom = vars[idx].chrom;
            results[idx].snp_id = vars[idx].id;
            results[idx].pos = vars[idx].pos;
            results[idx].a1 = vars[idx].alt;
            results[idx].a2 = vars[idx].ref;
            results[idx].n_obs = n;

            double score = scores[idx];
            double var_score = global_tau * n;

            double se = std::sqrt(var_score);
            double stat = (se > 0) ? (score / se) : 0.0;
            double pval = chisq_prob(stat * stat, 1);

            double beta = (var_score > 1e-9) ? (score / var_score) : 0.0;
            double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;

            results[idx].beta = beta;
            results[idx].se = beta_se;
            results[idx].p_value = pval;
            results[idx].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
        }

        // Update Warm Start?
        // x_global = x_loco; // Maybe?
    }
}

void LmmEngine::scanBiobank(BgenReader& reader,
                            std::vector<GwasResult>& results,
                            const std::vector<int>& keep_indices,
                            const std::string& chrom_filter) {
    int n = (int)y_raw.size();
    size_t m = reader.getNumVariants();

    if (keep_indices.size() != (size_t)n) throw std::runtime_error("Sample size mismatch in scanBiobank");

    std::cout << "Starting Biobank Scan (Matrix-Free PCG) on BGEN (" << m << " variants)..." << std::endl;
    std::cout << "  Sigma2_g: " << sigma2_g << ", Sigma2_e: " << sigma2_e << std::endl;

    // Enable BGEN prefetching for global solve
    reader.startPrefetch();

    double inv_M = 1.0 / (double)m;
    auto mult_op = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        Eigen::VectorXd u(m);
        reader.multiply_Zt_v(v_in, u);
        Eigen::VectorXd w(n);
        reader.multiply_Z_v(u, w);
        v_out = (sigma2_g * inv_M) * w + sigma2_e * v_in;
    };

    double scalar_precond = 1.0 / (sigma2_g + sigma2_e);
    auto precond_op = [&](const Eigen::VectorXd& v_in, Eigen::VectorXd& v_out) {
        v_out = v_in * scalar_precond;
    };

    MatrixFreePCG pcg(n);
    pcg.setMultOp(mult_op);
    pcg.setPrecondOp(precond_op);

    std::cout << "  Solving V * x = y..." << std::endl;
    Eigen::VectorXd x = pcg.solve(y_raw, 1e-5, 200);

    // Estimate Trace
    double tr_V_inv = estimateTraceInverse(pcg, n, 10);
    double tau = tr_V_inv / n;
    std::cout << "  Estimated Trace(V^-1) = " << tr_V_inv << ", tau = " << tau << std::endl;

    std::cout << "  Scanning variants..." << std::endl;
    results.resize(m);

    // Enable prefetch for scan
    reader.startPrefetch();

    Eigen::VectorXd scores(m);
    reader.multiply_Zt_v(x, scores);

    // Stop prefetch
    reader.stopPrefetch();

    const auto& vars = reader.getVariants();
    for(size_t i=0; i<m; ++i) {
        if (!options.extract_set.empty() && options.extract_set.find(vars[i].rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(vars[i].rsid) != options.exclude_set.end()) continue;
        if (vars[i].alleles.size() >= 2) {
             results[i].a1 = vars[i].alleles[1]; // Alt
             results[i].a2 = vars[i].alleles[0]; // Ref
        } else {
             results[i].a1 = "NA";
             results[i].a2 = "NA";
        }
        results[i].chrom = vars[i].chrom;
        results[i].snp_id = vars[i].rsid;
        results[i].pos = vars[i].pos;
        results[i].n_obs = n;

        double score = scores[i];
        double var_score = tau * n;

        double se = std::sqrt(var_score);
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);

        double beta = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;

        results[i].beta = beta;
        results[i].se = beta_se;
        results[i].p_value = pval;
        results[i].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
    }
}




LmmEngine::LmmEngine() {}
LmmEngine::LmmEngine(const Options& opts) : options(opts) {}

void LmmEngine::prepareFromVce(const EigenVce& vce) {
    sigma2_g = vce.getSigma2g();
    sigma2_e = vce.getSigma2e();
    null_delta = vce.getDelta();
    D = vce.getEigenvalues();
    W_proj = vce.getProjectedW();
    Uty = vce.getUty();
    UtX = vce.getUtX();
}

namespace {

struct EbMarkerFit {
    double lambda = 0.0;
    double sigma2 = 0.0;
    double beta = 0.0;
    double post_var = 0.0;
    double confidence = 0.0;
    double wald = 0.0;
    double lrt = 0.0;
    double lrt_p = 1.0;
    double p_value = 1.0;
    int conv = 0;
};

static double logdet_ldlt(const Eigen::MatrixXd& M) {
    if (M.rows() == 0) return 0.0;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
    if (ldlt.info() != Eigen::Success) return std::numeric_limits<double>::infinity();
    double out = 0.0;
    const auto& d = ldlt.vectorD();
    for (int i = 0; i < d.size(); ++i) {
        if (!(d(i) > 0.0) || !std::isfinite(d(i))) return std::numeric_limits<double>::infinity();
        out += std::log(d(i));
    }
    return out;
}

static EbMarkerFit fit_eb_marker(const Eigen::VectorXd& D,
                                 const Eigen::VectorXd& Uty,
                                 const Eigen::MatrixXd& UtX,
                                 const Eigen::VectorXd& Utg,
                                 double sigma2_g,
                                 double sigma2_e) {
    EbMarkerFit fit;
    const int n = (int)D.size();
    const int p = (int)UtX.cols();
    const int n_eff = std::max(1, n - p);
    const double ve = std::max(sigma2_e, 1e-12);
    const double eta = std::max(0.0, sigma2_g / ve);

    Eigen::VectorXd R = eta * D.array() + 1.0;
    for (int i = 0; i < R.size(); ++i) R(i) = std::max(R(i), 1e-10);
    Eigen::VectorXd Rinv = R.cwiseInverse();
    Eigen::VectorXd Rinv_z = Rinv.cwiseProduct(Utg);
    const double zRz = Utg.dot(Rinv_z);
    if (!(zRz > 1e-12) || !std::isfinite(zRz)) return fit;

    const double logdetR = R.array().log().sum();
    const double yRy = Uty.dot(Rinv.cwiseProduct(Uty));
    const double zRy = Utg.dot(Rinv.cwiseProduct(Uty));
    Eigen::MatrixXd XtRinvX;
    Eigen::VectorXd XtRinvy;
    Eigen::VectorXd XtRinvz;
    if (p > 0) {
        Eigen::MatrixXd RinvX = UtX.array().colwise() * Rinv.array();
        XtRinvX = UtX.transpose() * RinvX;
        XtRinvy = UtX.transpose() * Rinv.cwiseProduct(Uty);
        XtRinvz = UtX.transpose() * Rinv_z;
    }

    auto eval_loglike = [&](double lambda, double* out_quad = nullptr,
                            Eigen::VectorXd* out_beta = nullptr,
                            double* out_q = nullptr,
                            double* out_zHres = nullptr) -> double {
        lambda = std::max(0.0, lambda);
        const double denom = 1.0 + lambda * zRz;
        if (!(denom > 0.0) || !std::isfinite(denom)) return -std::numeric_limits<double>::infinity();
        const double factor = lambda / denom;
        const double yHy = yRy - factor * zRy * zRy;
        const double q = zRz / denom;
        Eigen::VectorXd XtHiny;
        Eigen::VectorXd XtHinz;
        Eigen::MatrixXd XtHinvX;
        double quad = yHy;
        Eigen::VectorXd beta;
        if (p > 0) {
            XtHiny = XtRinvy - factor * XtRinvz * zRy;
            XtHinz = XtRinvz / denom;
            XtHinvX = XtRinvX - factor * (XtRinvz * XtRinvz.transpose());
            Eigen::LDLT<Eigen::MatrixXd> solver(XtHinvX);
            if (solver.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
            beta = solver.solve(XtHiny);
            if (!beta.allFinite()) return -std::numeric_limits<double>::infinity();
            quad = yHy - XtHiny.dot(beta);
        }
        if (!(quad > 0.0) || !std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
        double logdetX = (p > 0) ? logdet_ldlt(XtHinvX) : 0.0;
        if (!std::isfinite(logdetX)) return -std::numeric_limits<double>::infinity();
        if (out_quad) *out_quad = quad;
        if (out_beta && p > 0) *out_beta = beta;
        if (out_q) *out_q = q;
        if (out_zHres) {
            double val = zRy / denom;
            if (p > 0) val -= XtHinz.dot(beta);
            *out_zHres = val;
        }
        return -0.5 * (logdetR + std::log(denom) + logdetX + n_eff * std::log(quad / (double)n_eff));
    };

    double best_lambda = 0.0;
    const double null_ll = eval_loglike(0.0);
    double best_ll = null_ll;
    const int grid_n = 49;
    std::vector<double> grid;
    grid.reserve(grid_n);
    for (int i = 0; i < grid_n; ++i) {
        double x = -12.0 + (24.0 * i / (double)(grid_n - 1));
        grid.push_back(std::exp(x));
    }
    int best_idx = -1;
    for (int i = 0; i < grid_n; ++i) {
        double ll = eval_loglike(grid[i]);
        if (ll > best_ll) { best_ll = ll; best_lambda = grid[i]; best_idx = i; }
    }
    if (best_idx >= 0 && best_ll > eval_loglike(0.0) + 1e-8) {
        double left = (best_idx == 0) ? std::log(grid[0]) - 2.0 : std::log(grid[best_idx - 1]);
        double right = (best_idx == grid_n - 1) ? std::log(grid.back()) + 2.0 : std::log(grid[best_idx + 1]);
        const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
        double c = right - gr * (right - left);
        double d = left + gr * (right - left);
        double fc = eval_loglike(std::exp(c));
        double fd = eval_loglike(std::exp(d));
        for (int iter = 0; iter < 60 && std::abs(right - left) > 1e-7; ++iter) {
            if (fc < fd) {
                left = c; c = d; fc = fd;
                d = left + gr * (right - left); fd = eval_loglike(std::exp(d));
            } else {
                right = d; d = c; fd = fc;
                c = right - gr * (right - left); fc = eval_loglike(std::exp(c));
            }
        }
        best_lambda = std::exp(0.5 * (left + right));
    }

    if (!(best_lambda > 1e-10) || !std::isfinite(best_lambda)) {
        fit.conv = std::isfinite(null_ll) ? 0 : 1;
        return fit;
    }
    double quad = 0.0, q = 0.0, zHres = 0.0;
    const double final_ll = eval_loglike(best_lambda, &quad, nullptr, &q, &zHres);
    if (!std::isfinite(final_ll)) {
        fit.conv = 1;
        return fit;
    }
    fit.lambda = best_lambda;
    fit.sigma2 = quad / (double)n_eff;
    fit.beta = best_lambda * zHres;
    fit.post_var = fit.sigma2 * std::max(0.0, best_lambda - best_lambda * best_lambda * q);
    fit.confidence = std::min(1.0, std::max(0.0, best_lambda * q));
    if (fit.post_var > 1e-20 && std::isfinite(fit.post_var)) {
        fit.wald = fit.beta * fit.beta / fit.post_var;
        fit.p_value = chisq_prob(fit.wald, 1);
    }
    if (std::isfinite(null_ll)) {
        fit.lrt = std::max(0.0, 2.0 * (final_ll - null_ll));
        fit.lrt_p = chisq_prob(fit.lrt, 1);
    }
    return fit;
}

} // namespace

void LmmEngine::scan(const PlinkReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    int m = reader.getNumSnps();
    int n = keep_indices.size(); // Original sample size
    int n_eig = D.size();        // Projected dimension (n or n-c)
    results.resize(m);

    Eigen::VectorXd V_inv_diag(n_eig);
    for(int i=0; i<n_eig; ++i) V_inv_diag[i] = 1.0 / (sigma2_g * D[i] + sigma2_e);

    Eigen::MatrixXd weighted_UtX;
    if (UtX.cols() > 0) {
        weighted_UtX = UtX.array().colwise() * V_inv_diag.array();
    }

    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    if (UtX.cols() > 0) Xt_Vinv_X = UtX.transpose() * weighted_UtX;

    Eigen::MatrixXd Xt_Vinv_X_inv(UtX.cols(), UtX.cols());
    if (UtX.cols() > 0) {
        Eigen::LLT<Eigen::MatrixXd> lltX(Xt_Vinv_X);
        if(lltX.info() == Eigen::Success) {
            Xt_Vinv_X_inv = lltX.solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
        } else {
            Xt_Vinv_X_inv = Xt_Vinv_X.inverse();
        }
    }

    Eigen::VectorXd Xt_Vinv_y = Eigen::VectorXd::Zero(UtX.cols());
    for(int i=0; i<n_eig; ++i) Xt_Vinv_y += V_inv_diag[i] * (UtX.row(i).transpose() * Uty[i]);

    Eigen::VectorXd beta(UtX.cols());
    if (UtX.cols() > 0) {
        beta = Xt_Vinv_X_inv * Xt_Vinv_y;
    }

    Eigen::VectorXd resid = Uty;
    if (UtX.cols() > 0) {
        resid -= UtX * beta;
    }

    Eigen::VectorXd P_y = resid.cwiseProduct(V_inv_diag);
    const double null_quad = std::max(resid.dot(P_y), std::numeric_limits<double>::min());
    const double n_eff_lmm = std::max(1.0, static_cast<double>(Uty.size() - UtX.cols()));
    ProfileRemlFit null_ml_refit;
    if (options.gemma_refit) null_ml_refit = profile_lmm_fit(D, Uty, UtX, false);

    const auto& bim = reader.getBimInfo();
    const uint8_t* raw_data = reader.getData();
    const size_t bytes_per_snp = reader.getBytesPerSnp();
    const size_t block_size = 512;
    const bool has_weights = (options.obs_weights.size() == n);

    Eigen::MatrixXd G_block(n, block_size);
    std::vector<int> snp_indices;
    std::vector<double> freqs;
    snp_indices.reserve(block_size);
    freqs.reserve(block_size);

    auto flush_block = [&]() {
        const int b = (int)snp_indices.size();
        if (b == 0) return;
        Eigen::MatrixXd G_use = G_block.leftCols(b);
        Eigen::MatrixXd UtG = W_proj.transpose() * G_use;
        Eigen::VectorXd scores = UtG.transpose() * P_y;
        Eigen::MatrixXd VinvG = UtG.array().colwise() * V_inv_diag.array();
        Eigen::VectorXd var_scores = (UtG.array() * VinvG.array()).colwise().sum().transpose();

        if (!options.use_residuals && UtX.cols() > 0) {
            Eigen::MatrixXd Xt_Vinv_G = weighted_UtX.transpose() * UtG;
            Eigen::MatrixXd middle = Xt_Vinv_X_inv * Xt_Vinv_G;
            Eigen::VectorXd term2 = (Xt_Vinv_G.array() * middle.array()).colwise().sum().matrix().transpose();
            var_scores -= term2;
        }

        for (int col = 0; col < b; ++col) {
            const int j = snp_indices[col];
            double var_score = var_scores[col];
            if (!(var_score > 1e-12) || !std::isfinite(var_score)) continue;
            double se = std::sqrt(var_score);
            double score = scores[col];
            double beta_snp = score / var_score;
            double beta_se = 1.0 / se;
            double stat = score / se;
            double pval = chisq_prob(stat * stat, 1);

            results[j].chrom = bim[j].chrom;
            results[j].snp_id = bim[j].id;
            results[j].pos = bim[j].bp_pos;
            results[j].a1 = bim[j].alt;
            results[j].a2 = bim[j].ref;
            results[j].freq = freqs[col];
            results[j].n_obs = n;
            results[j].beta = beta_snp;
            results[j].se = beta_se;
            results[j].p_value = pval;
            results[j].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
            fill_gemma_lmm_tests(results[j], score, var_score, null_quad, n_eff_lmm);
            if (options.profile_scale) {
                apply_fixed_ratio_profile(results[j], score, var_score, null_quad, n_eff_lmm);
            }
            if (options.gemma_refit) {
                Eigen::VectorXd alt_col = UtG.col(col);
                Eigen::MatrixXd alt_X = append_design_column(UtX, alt_col);
                ProfileRemlFit alt_ml_refit = profile_lmm_fit(D, Uty, alt_X, false);
                ProfileRemlFit alt_reml_refit = profile_lmm_fit(D, Uty, alt_X, true);
                apply_gemma_refit(results[j], null_ml_refit, alt_ml_refit, alt_reml_refit);
            }
        }
        snp_indices.clear();
        freqs.clear();
    };

    for (int j = 0; j < m; ++j) {
        if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 && bim[j].bp_pos < options.scan_from_bp) continue;
        if (options.scan_to_bp >= 0 && bim[j].bp_pos > options.scan_to_bp) continue;
        if (!options.extract_set.empty() && options.extract_set.find(bim[j].id) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(bim[j].id) != options.exclude_set.end()) continue;
        if (snp_indices.size() == block_size) flush_block();

        const uint8_t* snp_data = raw_data + 3 + (size_t)j * bytes_per_snp;
        double sum_g = 0.0;
        int valid_n = 0;
        int hom_ref = 0, het = 0, hom_alt = 0;
        int col = (int)snp_indices.size();

        for (int i = 0; i < n; ++i) {
            int sample_idx = keep_indices[i];
            uint8_t byte = snp_data[sample_idx >> 2];
            uint8_t code = (byte >> ((sample_idx & 3) * 2)) & 3;
            if (code != 1) {
                double val = (code == 0) ? 2.0 : (code == 2 ? 1.0 : 0.0);
                G_block(i, col) = val;
                sum_g += val;
                valid_n++;
                if (val == 0.0) hom_ref++;
                else if (val == 1.0) het++;
                else hom_alt++;
            } else {
                G_block(i, col) = std::numeric_limits<double>::quiet_NaN();
            }
        }

        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        if (valid_n < 3 || std::min(freq, 1.0 - freq) < options.min_maf) continue;
        if (options.min_hwe > 0.0 && calculateHWE(hom_ref, het, hom_alt) < options.min_hwe) continue;
        double mean = 2.0 * freq;
        for (int i = 0; i < n; ++i) {
            double val = G_block(i, col);
            G_block(i, col) = std::isnan(val) ? 0.0 : (val - mean);
        }
        if (has_weights) G_block.col(col).array() *= options.obs_weights.array();

        snp_indices.push_back(j);
        freqs.push_back(freq);
    }
    flush_block();
}

void LmmEngine::scanEB(const PlinkReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    const int m = (int)reader.getNumSnps();
    const int n = (int)keep_indices.size();
    results.resize(m);
    const auto& bim = reader.getBimInfo();
    const uint8_t* raw_data = reader.getData();
    const size_t bytes_per_snp = reader.getBytesPerSnp();
    const bool has_weights = (options.obs_weights.size() == n);

    #pragma omp parallel for schedule(dynamic, 16)
    for (int j = 0; j < m; ++j) {
        if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
        if (!options.extract_set.empty() && options.extract_set.find(bim[j].id) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(bim[j].id) != options.exclude_set.end()) continue;

        Eigen::VectorXd g(n);
        const uint8_t* snp_data = raw_data + 3 + (size_t)j * bytes_per_snp;
        double sum_g = 0.0;
        int valid_n = 0;
        for (int i = 0; i < n; ++i) {
            int sample_idx = keep_indices[i];
            uint8_t byte = snp_data[sample_idx >> 2];
            uint8_t code = (byte >> ((sample_idx & 3) * 2)) & 3;
            if (code != 1) {
                double val = (code == 0) ? 2.0 : (code == 2 ? 1.0 : 0.0);
                g[i] = val;
                sum_g += val;
                valid_n++;
            } else {
                g[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        if (valid_n < 3) continue;
        double freq = sum_g / (2.0 * valid_n);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options.min_maf) continue;
        double mean = 2.0 * freq;
        for (int i = 0; i < n; ++i) g[i] = std::isnan(g[i]) ? 0.0 : (g[i] - mean);
        if (has_weights) g.array() *= options.obs_weights.array();

        Eigen::VectorXd Utg = W_proj.transpose() * g;
        EbMarkerFit eb = fit_eb_marker(D, Uty, UtX, Utg, sigma2_g, sigma2_e);

        results[j].chrom = bim[j].chrom;
        results[j].snp_id = bim[j].id;
        results[j].pos = bim[j].bp_pos;
        results[j].a1 = bim[j].alt;
        results[j].a2 = bim[j].ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = eb.beta;
        results[j].se = (eb.post_var > 0.0) ? std::sqrt(eb.post_var) : 0.0;
        results[j].p_value = eb.p_value;
        results[j].log10_p = (eb.p_value > 0.0) ? std::max(0.0, -std::log10(eb.p_value)) : 300.0;
        results[j].eb_marker_lambda = eb.lambda;
        results[j].eb_prior_var = eb.lambda * eb.sigma2;
        results[j].eb_post_var = eb.post_var;
        results[j].eb_confidence = eb.confidence;
        results[j].eb_wald = eb.wald;
        results[j].eb_lrt = eb.lrt;
        results[j].eb_lrt_p = eb.lrt_p;
        results[j].eb_sigma2 = eb.sigma2;
        results[j].eb_conv = eb.conv;
    }
}

// --------------------------------------------------------

void LmmEngine::scan(PgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    int n = keep_indices.size();
    size_t m = reader.getNumVariants();
    results.resize(m);

    int n_eig = D.size();

    // Prepare V inverse diagonal
    Eigen::VectorXd V_inv_diag(n_eig);
    for(int i=0; i<n_eig; ++i) V_inv_diag[i] = 1.0 / (sigma2_g * D[i] + sigma2_e);

    Eigen::MatrixXd weighted_UtX;
    if (UtX.cols() > 0) {
        weighted_UtX = UtX.array().colwise() * V_inv_diag.array();
    }

    // Fixed Effects Part
    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    if (UtX.cols() > 0) Xt_Vinv_X = UtX.transpose() * weighted_UtX;
    Eigen::MatrixXd Xt_Vinv_X_inv;
    if (UtX.cols() > 0) {
        Eigen::LLT<Eigen::MatrixXd> lltX959(Xt_Vinv_X);
        if(lltX959.info() == Eigen::Success) {
            Xt_Vinv_X_inv = lltX959.solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
        } else {
            Xt_Vinv_X_inv = Xt_Vinv_X.inverse();
        }
    }

    Eigen::VectorXd Xt_Vinv_y = Eigen::VectorXd::Zero(UtX.cols());
    for(int i=0; i<n_eig; ++i) Xt_Vinv_y += V_inv_diag[i] * (UtX.row(i).transpose() * Uty[i]);

    Eigen::VectorXd beta_null(UtX.cols());
    if (UtX.cols() > 0) {
        beta_null = Xt_Vinv_X_inv * Xt_Vinv_y;
    }

    Eigen::VectorXd resid = Uty;
    if (UtX.cols() > 0) {
        resid -= UtX * beta_null;
    }

    Eigen::VectorXd P_y = resid.cwiseProduct(V_inv_diag);
    const double null_quad = std::max(resid.dot(P_y), std::numeric_limits<double>::min());
    const double n_eff_lmm = std::max(1.0, static_cast<double>(Uty.size() - UtX.cols()));
    ProfileRemlFit null_ml_refit;
    if (options.gemma_refit) null_ml_refit = profile_lmm_fit(D, Uty, UtX, false);

    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;

    reader.reset();

    // Loop through variants
    for(size_t j=0; j<m; ++j) {
        bool ok = reader.readVariant(chrom, rsid, pos, ref, alt, dosage);
        if (!ok) break;

        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 && pos < static_cast<uint32_t>(options.scan_from_bp)) continue;
        if (options.scan_to_bp >= 0 && pos > static_cast<uint32_t>(options.scan_to_bp)) continue;
        if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;

        // Extract subset and center
        Eigen::VectorXd g(n);
        double sum_g = 0.0;
        int valid_n = 0;

        for(int i=0; i<n; ++i) {
            float val = dosage[keep_indices[i]]; // Pgen dosage is 0..2 or NAN
            g[i] = val;
            if (!std::isnan(val)) {
                sum_g += val;
                valid_n++;
            }
        }

        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        if (valid_n < 3 || std::min(freq, 1.0 - freq) < options.min_maf) continue;
        double mean = 2.0 * freq;

        for(int i=0; i<n; ++i) {
            if (std::isnan(g[i])) g[i] = 0.0; // Mean imputation (centered)
            else g[i] -= mean;
        }

        // Apply observation weights if provided
        if (options.obs_weights.size() == n) {
            g.array() *= options.obs_weights.array();
        }

        // Rotate
        Eigen::VectorXd Utg;
        if (options.use_covariate_projection && W_proj.cols() > 0) {
            Utg = W_proj.transpose() * g;
        } else {
            Utg = W_proj.transpose() * g;
        }

        double score = Utg.dot(P_y);
        double var_score = 0.0;

        Eigen::VectorXd Vinv_g = Utg.cwiseProduct(V_inv_diag);
        double g_Vinv_g = Utg.dot(Vinv_g);

        if (options.use_residuals || UtX.cols() == 0) {
            var_score = g_Vinv_g;
        } else {
            Eigen::VectorXd Xt_Vinv_g = weighted_UtX.transpose() * Utg;
            double term2 = Xt_Vinv_g.transpose() * Xt_Vinv_X_inv * Xt_Vinv_g;
            var_score = g_Vinv_g - term2;
        }

        double se = std::sqrt(var_score);
        double beta_snp = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);

        results[j].chrom = chrom;
        results[j].snp_id = rsid;
        results[j].pos = pos;
        results[j].a1 = alt;
        results[j].a2 = ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = beta_snp;
        results[j].se = beta_se;
        results[j].p_value = pval;
        results[j].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
        fill_gemma_lmm_tests(results[j], score, var_score, null_quad, n_eff_lmm);
        if (options.profile_scale) {
            apply_fixed_ratio_profile(results[j], score, var_score, null_quad, n_eff_lmm);
        }
        if (options.gemma_refit) {
            Eigen::MatrixXd alt_X = append_design_column(UtX, Utg);
            ProfileRemlFit alt_ml_refit = profile_lmm_fit(D, Uty, alt_X, false);
            ProfileRemlFit alt_reml_refit = profile_lmm_fit(D, Uty, alt_X, true);
            apply_gemma_refit(results[j], null_ml_refit, alt_ml_refit, alt_reml_refit);
        }
    }
}

void LmmEngine::scanEB(PgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    const int n = (int)keep_indices.size();
    const size_t m = reader.getNumVariants();
    results.resize(m);
    const bool has_weights = (options.obs_weights.size() == n);
    reader.reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    Eigen::VectorXd g(n);

    for (size_t j = 0; j < m; ++j) {
        if (!reader.readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 && pos < static_cast<uint32_t>(options.scan_from_bp)) continue;
        if (options.scan_to_bp >= 0 && pos > static_cast<uint32_t>(options.scan_to_bp)) continue;
        if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;
        double sum_g = 0.0;
        int valid_n = 0;
        for (int i = 0; i < n; ++i) {
            float val = dosage[keep_indices[i]];
            g[i] = val;
            if (!std::isnan(val)) { sum_g += val; valid_n++; }
        }
        if (valid_n < 3) continue;
        double freq = sum_g / (2.0 * valid_n);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options.min_maf) continue;
        double mean = 2.0 * freq;
        for (int i = 0; i < n; ++i) g[i] = std::isnan(g[i]) ? 0.0 : (g[i] - mean);
        if (has_weights) g.array() *= options.obs_weights.array();

        EbMarkerFit eb = fit_eb_marker(D, Uty, UtX, W_proj.transpose() * g, sigma2_g, sigma2_e);
        results[j].chrom = chrom;
        results[j].snp_id = rsid;
        results[j].pos = pos;
        results[j].a1 = alt;
        results[j].a2 = ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = eb.beta;
        results[j].se = (eb.post_var > 0.0) ? std::sqrt(eb.post_var) : 0.0;
        results[j].p_value = eb.p_value;
        results[j].log10_p = (eb.p_value > 0.0) ? std::max(0.0, -std::log10(eb.p_value)) : 300.0;
        results[j].eb_marker_lambda = eb.lambda;
        results[j].eb_prior_var = eb.lambda * eb.sigma2;
        results[j].eb_post_var = eb.post_var;
        results[j].eb_confidence = eb.confidence;
        results[j].eb_wald = eb.wald;
        results[j].eb_lrt = eb.lrt;
        results[j].eb_lrt_p = eb.lrt_p;
        results[j].eb_sigma2 = eb.sigma2;
        results[j].eb_conv = eb.conv;
    }
}

void LmmEngine::scan(BgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    // Similar to PgenReader but using BgenReader interface
    // Assuming BgenReader has similar readVariant or iterator
    int n = keep_indices.size();
    size_t m = reader.getNumVariants();
    results.resize(m);

    int n_eig = D.size();

    Eigen::VectorXd V_inv_diag(n_eig);
    for(int i=0; i<n_eig; ++i) V_inv_diag[i] = 1.0 / (sigma2_g * D[i] + sigma2_e);

    Eigen::MatrixXd weighted_UtX;
    if (UtX.cols() > 0) {
        weighted_UtX = UtX.array().colwise() * V_inv_diag.array();
    }

    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    if (UtX.cols() > 0) Xt_Vinv_X = UtX.transpose() * weighted_UtX;
    Eigen::MatrixXd Xt_Vinv_X_inv;
    if (UtX.cols() > 0) {
        Eigen::LLT<Eigen::MatrixXd> lltX1058(Xt_Vinv_X);
        if(lltX1058.info() == Eigen::Success) {
            Xt_Vinv_X_inv = lltX1058.solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
        } else {
            Xt_Vinv_X_inv = Xt_Vinv_X.inverse();
        }
    }

    Eigen::VectorXd Xt_Vinv_y = Eigen::VectorXd::Zero(UtX.cols());
    for(int i=0; i<n_eig; ++i) Xt_Vinv_y += V_inv_diag[i] * (UtX.row(i).transpose() * Uty[i]);

    Eigen::VectorXd beta_null(UtX.cols());
    if (UtX.cols() > 0) {
        beta_null = Xt_Vinv_X_inv * Xt_Vinv_y;
    }

    Eigen::VectorXd resid = Uty;
    if (UtX.cols() > 0) {
        resid -= UtX * beta_null;
    }

    Eigen::VectorXd P_y = resid.cwiseProduct(V_inv_diag);
    const double null_quad = std::max(resid.dot(P_y), std::numeric_limits<double>::min());
    const double n_eff_lmm = std::max(1.0, static_cast<double>(Uty.size() - UtX.cols()));
    ProfileRemlFit null_ml_refit;
    if (options.gemma_refit) null_ml_refit = profile_lmm_fit(D, Uty, UtX, false);

    reader.reset();

    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f;

    for(size_t j=0; j<m; ++j) {
        bool ok = reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        if (!ok) break;

        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        if (options.scan_from_bp >= 0 && pos < static_cast<uint32_t>(options.scan_from_bp)) continue;
        if (options.scan_to_bp >= 0 && pos > static_cast<uint32_t>(options.scan_to_bp)) continue;
        if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;

        Eigen::VectorXd g(n);
        double sum_g = 0.0;
        int valid_n = 0;

        for(int i=0; i<n; ++i) {
            float val = dosage_f[keep_indices[i]];
            g[i] = val;
            if (!std::isnan(val)) {
                sum_g += val;
                valid_n++;
            }
        }

        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        if (valid_n < 3 || std::min(freq, 1.0 - freq) < options.min_maf) continue;
        double mean = 2.0 * freq;

        for(int i=0; i<n; ++i) {
            if (std::isnan(g[i])) g[i] = 0.0;
            else g[i] -= mean;
        }

        // Apply observation weights if provided
        if (options.obs_weights.size() == n) {
            g.array() *= options.obs_weights.array();
        }

        Eigen::VectorXd Utg;
        if (options.use_covariate_projection && W_proj.cols() > 0) {
            Utg = W_proj.transpose() * g;
        } else {
            Utg = W_proj.transpose() * g;
        }

        double score = Utg.dot(P_y);
        double var_score = 0.0;

        Eigen::VectorXd Vinv_g = Utg.cwiseProduct(V_inv_diag);
        double g_Vinv_g = Utg.dot(Vinv_g);

        if (options.use_residuals || UtX.cols() == 0) {
            var_score = g_Vinv_g;
        } else {
            Eigen::VectorXd Xt_Vinv_g = weighted_UtX.transpose() * Utg;
            double term2 = Xt_Vinv_g.transpose() * Xt_Vinv_X_inv * Xt_Vinv_g;
            var_score = g_Vinv_g - term2;
        }

        double se = std::sqrt(var_score);
        double beta_snp = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);

        results[j].chrom = chrom;
        results[j].snp_id = rsid;
        results[j].pos = pos;
        results[j].a1 = alt;
        results[j].a2 = ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = beta_snp;
        results[j].se = beta_se;
        results[j].p_value = pval;
        results[j].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
        fill_gemma_lmm_tests(results[j], score, var_score, null_quad, n_eff_lmm);
        if (options.profile_scale) {
            apply_fixed_ratio_profile(results[j], score, var_score, null_quad, n_eff_lmm);
        }
        if (options.gemma_refit) {
            Eigen::MatrixXd alt_X = append_design_column(UtX, Utg);
            ProfileRemlFit alt_ml_refit = profile_lmm_fit(D, Uty, alt_X, false);
            ProfileRemlFit alt_reml_refit = profile_lmm_fit(D, Uty, alt_X, true);
            apply_gemma_refit(results[j], null_ml_refit, alt_ml_refit, alt_reml_refit);
        }
    }
}

void LmmEngine::scanEB(BgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    const int n = (int)keep_indices.size();
    const size_t m = reader.getNumVariants();
    results.resize(m);
    const bool has_weights = (options.obs_weights.size() == n);
    reader.reset();
    std::string chrom, rsid, ref, alt;
    uint32_t pos = 0;
    Eigen::VectorXf dosage;
    Eigen::VectorXd g(n);

    for (size_t j = 0; j < m; ++j) {
        if (!reader.readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        if (!options.extract_set.empty() && options.extract_set.find(rsid) == options.extract_set.end()) continue;
        if (!options.exclude_set.empty() && options.exclude_set.find(rsid) != options.exclude_set.end()) continue;
        double sum_g = 0.0;
        int valid_n = 0;
        for (int i = 0; i < n; ++i) {
            float val = dosage[keep_indices[i]];
            g[i] = val;
            if (!std::isnan(val)) { sum_g += val; valid_n++; }
        }
        if (valid_n < 3) continue;
        double freq = sum_g / (2.0 * valid_n);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < options.min_maf) continue;
        double mean = 2.0 * freq;
        for (int i = 0; i < n; ++i) g[i] = std::isnan(g[i]) ? 0.0 : (g[i] - mean);
        if (has_weights) g.array() *= options.obs_weights.array();

        EbMarkerFit eb = fit_eb_marker(D, Uty, UtX, W_proj.transpose() * g, sigma2_g, sigma2_e);
        results[j].chrom = chrom;
        results[j].snp_id = rsid;
        results[j].pos = pos;
        results[j].a1 = alt;
        results[j].a2 = ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = eb.beta;
        results[j].se = (eb.post_var > 0.0) ? std::sqrt(eb.post_var) : 0.0;
        results[j].p_value = eb.p_value;
        results[j].log10_p = (eb.p_value > 0.0) ? std::max(0.0, -std::log10(eb.p_value)) : 300.0;
        results[j].eb_marker_lambda = eb.lambda;
        results[j].eb_prior_var = eb.lambda * eb.sigma2;
        results[j].eb_post_var = eb.post_var;
        results[j].eb_confidence = eb.confidence;
        results[j].eb_wald = eb.wald;
        results[j].eb_lrt = eb.lrt;
        results[j].eb_lrt_p = eb.lrt_p;
        results[j].eb_sigma2 = eb.sigma2;
        results[j].eb_conv = eb.conv;
    }
}

} // namespace cosmic
