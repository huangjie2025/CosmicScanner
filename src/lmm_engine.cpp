#include "lmm_engine.h"
#include "vce.h"
#include "matrix_free_pcg.h"
#include "qc_utils.h"
#include "sparse_grm_reader.h"
#include <iostream>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <map>
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
        // FIXME: PgenReader doesn't support keep_indices yet, assuming empty or implemented later
        reader.multiply_Zt_v(v_in, u); 
        Eigen::VectorXd w(n);
        reader.multiply_Z_v(u, w);     
        
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
            reader.multiply_Zt_v(v_in, u); 
            
            // Mask out current chromosome
            for (size_t idx : indices) u(idx) = 0.0;
            
            Eigen::VectorXd w(n);
            reader.multiply_Z_v(u, w);     
            
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
        reader.multiply_Zt_v(x_loco, scores);
        
        for(size_t idx : indices) {
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

void LmmEngine::scan(const PlinkReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    int m = reader.getNumSnps();
    int n = keep_indices.size(); // Original sample size
    int n_eig = D.size();        // Projected dimension (n or n-c)
    results.resize(m);
    
    Eigen::VectorXd V_inv_diag(n_eig);
    for(int i=0; i<n_eig; ++i) V_inv_diag[i] = 1.0 / (sigma2_g * D[i] + sigma2_e);
    
    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    Xt_Vinv_X.setZero();
    for(int i=0; i<n_eig; ++i) Xt_Vinv_X += V_inv_diag[i] * (UtX.row(i).transpose() * UtX.row(i));
    
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
    
    SPATest::Options spa_opts;
    spa_opts.pval_threshold = 0.05;
    SPATest spa(spa_opts);
    if (options.use_spa) {
        // Use null model variance for SPA (mu=0 for centered LMM, V_ii = sigma2_g + sigma2_e)
        Eigen::VectorXd null_mu = Eigen::VectorXd::Constant(n, 0.0);
        Eigen::VectorXd null_var = Eigen::VectorXd::Constant(n, sigma2_g + sigma2_e);
        spa.set_null_model(null_mu, null_var);
    }
    
    const auto& bim = reader.getBimInfo();
    #pragma omp parallel for schedule(dynamic, 100) if(options.thread_num > 1)
    for(int j=0; j<m; ++j) {
        if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
        
        Eigen::VectorXf g_raw;
        reader.getSnpGenotypes(j, g_raw);
        
        Eigen::VectorXd g(n);
        double sum_g = 0.0;
        int valid_n = 0;
        for(int i=0; i<n; ++i) {
            g[i] = g_raw[keep_indices[i]];
            if (!std::isnan(g[i])) {
                sum_g += g[i];
                valid_n++;
            }
        }
        
        double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
        double mean = 2.0 * freq;
        
        // Center genotypes (critical for MLMA/Residual mode and Covariate Projection)
        for(int i=0; i<n; ++i) {
            if (std::isnan(g[i])) g[i] = 0.0;
            else g[i] -= mean;
        }
        
        // Apply observation weights if provided (e.g., for Repeatability Model D^{-1/2} transformation)
        if (options.obs_weights.size() == n) {
            g.array() *= options.obs_weights.array();
        }
        
        Eigen::VectorXd Utg;
        if (options.use_covariate_projection && W_proj.cols() > 0) {
            Utg = W_proj.transpose() * g;
        } else {
            // If we don't have projection, we assume W_proj contains standard U
            Utg = W_proj.transpose() * g;
        }
        
        double score = Utg.dot(P_y);
        double var_score = 0.0;
        
        Eigen::VectorXd Vinv_g = Utg.cwiseProduct(V_inv_diag);
        double g_Vinv_g = Utg.dot(Vinv_g);
        
        if (options.use_residuals || UtX.cols() == 0) {
            // GCTA style or Covariate Projection: variance is g' V^-1 g
            // If covariates are already projected out (UtX.cols() == 0), exact variance is also g' V^-1 g
            var_score = g_Vinv_g;
        } else {
            // Exact LMM style: variance is g' P g = g' V^-1 g - g' V^-1 X (X' V^-1 X)^-1 X' V^-1 g
            Eigen::VectorXd Xt_Vinv_g = Eigen::VectorXd::Zero(UtX.cols());
            for(int i=0; i<n_eig; ++i) Xt_Vinv_g += V_inv_diag[i] * (UtX.row(i).transpose() * Utg[i]);
            double term2 = Xt_Vinv_g.transpose() * Xt_Vinv_X_inv * Xt_Vinv_g;
            var_score = g_Vinv_g - term2;
        }
        
        double se = std::sqrt(var_score);
        double beta_snp = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);
        
        if (options.use_spa && pval < spa_opts.pval_threshold) {
            pval = spa.compute_pvalue(g, score, var_score);
        }
        
        results[j].chrom = bim[j].chrom;
        results[j].snp_id = bim[j].id;
        results[j].pos = bim[j].bp_pos;
        results[j].a1 = bim[j].alt;
        results[j].a2 = bim[j].ref;
        results[j].freq = freq;
        results[j].n_obs = n;
        results[j].beta = beta_snp;
        results[j].se = beta_se;
        results[j].p_value = pval;
        results[j].log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
    }
}

// --------------------------------------------------------

void LmmEngine::loadSparseGRM(const std::string& prefix) {
    cosmic::SparseGrmReader::read(prefix, this->sparse_grm);
    this->has_sparse_grm = true;
}

void LmmEngine::scan(PgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter) {
    int n = keep_indices.size();
    size_t m = reader.getNumVariants();
    results.resize(m);
    
    int n_eig = D.size();
    
    // Prepare V inverse diagonal
    Eigen::VectorXd V_inv_diag(n_eig);
    for(int i=0; i<n_eig; ++i) V_inv_diag[i] = 1.0 / (sigma2_g * D[i] + sigma2_e);
    
    // Fixed Effects Part
    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    Xt_Vinv_X.setZero();
    for(int i=0; i<n_eig; ++i) Xt_Vinv_X += V_inv_diag[i] * (UtX.row(i).transpose() * UtX.row(i));
    Eigen::MatrixXd Xt_Vinv_X_inv;
    Eigen::LLT<Eigen::MatrixXd> lltX959(Xt_Vinv_X);
    if(lltX959.info() == Eigen::Success) {
        Xt_Vinv_X_inv = lltX959.solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
    } else {
        Xt_Vinv_X_inv = Xt_Vinv_X.inverse();
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
    
    SPATest::Options spa_opts;
    spa_opts.pval_threshold = 0.05;
    SPATest spa(spa_opts);
    if (options.use_spa) {
        // Use null model variance for SPA (mu=0 for centered LMM, V_ii = sigma2_g + sigma2_e)
        Eigen::VectorXd null_mu = Eigen::VectorXd::Constant(n, 0.0);
        Eigen::VectorXd null_var = Eigen::VectorXd::Constant(n, sigma2_g + sigma2_e);
        spa.set_null_model(null_mu, null_var);
    }
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage;
    
    reader.reset();
    
    // Loop through variants
    for(size_t j=0; j<m; ++j) {
        bool ok = reader.readVariant(chrom, rsid, pos, ref, alt, dosage);
        if (!ok) break;
        
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        
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
        
        if (options.use_residuals) {
            var_score = g_Vinv_g;
        } else {
            Eigen::VectorXd Xt_Vinv_g = Eigen::VectorXd::Zero(UtX.cols());
            int n_eig = D.size();
            for(int i=0; i<n_eig; ++i) Xt_Vinv_g += V_inv_diag[i] * (UtX.row(i).transpose() * Utg[i]);
            double term2 = Xt_Vinv_g.transpose() * Xt_Vinv_X_inv * Xt_Vinv_g;
            var_score = g_Vinv_g - term2;
        }

        double se = std::sqrt(var_score);
        double beta_snp = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);
        
        if (options.use_spa && pval < spa_opts.pval_threshold) {
            pval = spa.compute_pvalue(g, score, var_score);
        }
        
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
    
    Eigen::MatrixXd Xt_Vinv_X(UtX.cols(), UtX.cols());
    Xt_Vinv_X.setZero();
    for(int i=0; i<n_eig; ++i) Xt_Vinv_X += V_inv_diag[i] * (UtX.row(i).transpose() * UtX.row(i));
    Eigen::MatrixXd Xt_Vinv_X_inv;
    Eigen::LLT<Eigen::MatrixXd> lltX1058(Xt_Vinv_X);
    if(lltX1058.info() == Eigen::Success) {
        Xt_Vinv_X_inv = lltX1058.solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
    } else {
        Xt_Vinv_X_inv = Xt_Vinv_X.inverse();
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
    
    SPATest::Options spa_opts;
    spa_opts.pval_threshold = 0.05;
    SPATest spa(spa_opts);
    if (options.use_spa) {
        // Use null model variance for SPA (mu=0 for centered LMM, V_ii = sigma2_g + sigma2_e)
        Eigen::VectorXd null_mu = Eigen::VectorXd::Constant(n, 0.0);
        Eigen::VectorXd null_var = Eigen::VectorXd::Constant(n, sigma2_g + sigma2_e);
        spa.set_null_model(null_mu, null_var);
    }
    
    reader.reset();
    
    std::string chrom, rsid, ref, alt;
    uint32_t pos;
    Eigen::VectorXf dosage_f; 
    
    for(size_t j=0; j<m; ++j) {
        bool ok = reader.readVariant(chrom, rsid, pos, ref, alt, dosage_f);
        if (!ok) break;
        
        if (!chrom_filter.empty() && chrom != chrom_filter) continue;
        
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
        
        if (options.use_residuals) {
            var_score = g_Vinv_g;
        } else {
            Eigen::VectorXd Xt_Vinv_g = Eigen::VectorXd::Zero(UtX.cols());
            for(int i=0; i<n_eig; ++i) Xt_Vinv_g += V_inv_diag[i] * (UtX.row(i).transpose() * Utg[i]);
            double term2 = Xt_Vinv_g.transpose() * Xt_Vinv_X_inv * Xt_Vinv_g;
            var_score = g_Vinv_g - term2;
        }
        
        double se = std::sqrt(var_score);
        double beta_snp = (var_score > 1e-9) ? (score / var_score) : 0.0;
        double beta_se = (var_score > 1e-9) ? (1.0 / se) : NAN;
        double stat = (se > 0) ? (score / se) : 0.0;
        double pval = chisq_prob(stat * stat, 1);
        
        if (options.use_spa && pval < spa_opts.pval_threshold) {
            pval = spa.compute_pvalue(g, score, var_score);
        }
        
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
    }
}

} // namespace cosmic
