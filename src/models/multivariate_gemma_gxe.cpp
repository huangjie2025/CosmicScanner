#include "models/multivariate_model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <omp.h>
#include "grm_builder.h"
#include "mvlmm_engine.h"

namespace cosmic {
namespace {

std::ofstream mv_open_output_file_gxe(const std::string& path, std::ios_base::openmode mode = std::ios::out) {
    std::ofstream ofs(path, mode);
    if (!ofs) throw std::runtime_error("Failed to open output file: " + path);
    return ofs;
}

void mv_write_vbeta_header_gxe(std::ostream& os, int n_traits) {
    for (int i = 0; i < n_traits; ++i) {
        for (int j = i; j < n_traits; ++j) os << "\tVBETA_" << (i + 1) << "_" << (j + 1);
    }
}

void mv_write_vbeta_values_gxe(std::ostream& os, const Eigen::MatrixXd& vbeta, int n_traits) {
    for (int i = 0; i < n_traits; ++i) {
        for (int j = i; j < n_traits; ++j) {
            double value = (i < vbeta.rows() && j < vbeta.cols()) ? vbeta(i, j) : std::numeric_limits<double>::quiet_NaN();
            os << "\t" << value;
        }
    }
}

double mv_gamma_q_gxe(double a, double x) {
    if (x <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0; del *= x / ap; sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7) return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
        }
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1e-30, d = 1.0 / b, h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a); b += 2.0;
            d = an * d + b; if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c; if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d; double del = d * c; h *= del;
            if (std::abs(del - 1.0) < 1e-7) return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
        }
    }
    return 0.0;
}

double mv_chisq_prob_gxe(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return mv_gamma_q_gxe(df / 2.0, stat / 2.0);
}

bool mv_fit_gls_gxe(const Eigen::MatrixXd& UtY,
                    const Eigen::MatrixXd& UtX,
                    const Eigen::VectorXd& D,
                    const Eigen::MatrixXd& Vg,
                    const Eigen::MatrixXd& Ve,
                    Eigen::VectorXd& beta_all,
                    Eigen::MatrixXd& beta_cov,
                    double* loglik = nullptr) {
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
            for (int b = 0; b < n_covars; ++b) H.block(a * n_traits, b * n_traits, n_traits, n_traits) += UtX(i, a) * UtX(i, b) * sigma_inv[i];
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
            for (int a = 0; a < n_covars; ++a) fit += UtX(i, a) * beta_all.segment(a * n_traits, n_traits);
            Eigen::VectorXd resid = UtY.row(i).transpose() - fit;
            ll += -0.5 * (logdet[i] + (resid.transpose() * sigma_inv[i] * resid)(0, 0));
        }
        *loglik = ll;
    }
    return true;
}

bool mv_em_refit_rotated_gxe(const Eigen::MatrixXd& UtY,
                             const Eigen::MatrixXd& UtX,
                             const Eigen::VectorXd& D,
                             int max_iter,
                             double tol,
                             Eigen::MatrixXd& Vg,
                             Eigen::MatrixXd& Ve,
                             const Eigen::MatrixXd* init_vg = nullptr,
                             const Eigen::MatrixXd* init_ve = nullptr) {
    int n = UtY.rows();
    int d = UtY.cols();
    if (init_vg && init_ve) { Vg = *init_vg; Ve = *init_ve; }
    else {
        Eigen::MatrixXd centered = UtY.rowwise() - UtY.colwise().mean();
        Eigen::MatrixXd cov = (centered.transpose() * centered) / std::max(1, n - 1);
        Vg = 0.5 * cov; Ve = 0.5 * cov;
    }
    Vg = 0.5 * (Vg + Vg.transpose()); Ve = 0.5 * (Ve + Ve.transpose());
    Vg.diagonal().array() += 1e-8; Ve.diagonal().array() += 1e-8;
    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd beta; Eigen::MatrixXd beta_cov;
        if (!mv_fit_gls_gxe(UtY, UtX, D, Vg, Ve, beta, beta_cov, nullptr)) return false;
        Eigen::MatrixXd Sg = Eigen::MatrixXd::Zero(d, d), Se = Eigen::MatrixXd::Zero(d, d);
        for (int i = 0; i < n; ++i) {
            Eigen::VectorXd fit = Eigen::VectorXd::Zero(d);
            for (int c = 0; c < UtX.cols(); ++c) fit += UtX(i, c) * beta.segment(c * d, d);
            Eigen::VectorXd resid = UtY.row(i).transpose() - fit;
            Eigen::MatrixXd Sigma = D(i) * Vg + Ve; Sigma.diagonal().array() += 1e-8;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(Sigma); if (ldlt.info() != Eigen::Success) return false;
            Eigen::MatrixXd Pinv = ldlt.solve(Eigen::MatrixXd::Identity(d, d));
            Eigen::MatrixXd Wg = Vg * D(i) * Pinv; Eigen::VectorXd uhat = Wg * resid;
            Sg += uhat * uhat.transpose() + Vg * D(i) - Wg * Vg * D(i);
            Eigen::MatrixXd We = Ve * Pinv; Eigen::VectorXd ehat = We * resid;
            Se += ehat * ehat.transpose() + Ve - We * Ve;
        }
        Eigen::MatrixXd Vg_new = 0.5 * (Sg / n + (Sg / n).transpose());
        Eigen::MatrixXd Ve_new = 0.5 * (Se / n + (Se / n).transpose());
        Vg_new.diagonal().array() += 1e-8; Ve_new.diagonal().array() += 1e-8;
        double err = ((Vg_new - Vg).squaredNorm() + (Ve_new - Ve).squaredNorm()) / std::max(1.0, 2.0 * d * d);
        Vg = Vg_new; Ve = Ve_new;
        if (!Vg.allFinite() || !Ve.allFinite()) return false;
        if (err < tol) return true;
    }
    return true;
}

bool mv_scan_filter_passes_gxe(const PlinkBimInfo& var, const MultivariateModel::Options& options) {
    if (!options.scan_chrom.empty() && var.chrom != options.scan_chrom) return false;
    if (options.scan_from_bp >= 0 && var.bp_pos < options.scan_from_bp) return false;
    if (options.scan_to_bp >= 0 && var.bp_pos > options.scan_to_bp) return false;
    if (!options.extract_set.empty() && options.extract_set.find(var.id) == options.extract_set.end()) return false;
    if (!options.exclude_set.empty() && options.exclude_set.find(var.id) != options.exclude_set.end()) return false;
    return true;
}

} // namespace

void MultivariateModel::runGemmaGxe() {
    if (!options_.no_loco) throw std::runtime_error("--mv-method gemma-gxe currently requires --no-loco");
    std::cout << "Multivariate GWAS method: "
              << (options_.use_gemma_full ? "GEMMA-style multivariate GxE Wald/LRT/Score" : "GEMMA-style multivariate GxE interaction Wald (fast null)")
              << "\n";

    const int n = n_keep_;
    const int q = Y_.cols();
    Eigen::VectorXd env = env_;
    double emean = env.mean();
    double esd = std::sqrt((env.array() - emean).square().mean());
    if (esd > 1e-12) env = (env.array() - emean) / esd;
    Eigen::MatrixXd X_env(n, X_.cols() + 1);
    X_env.leftCols(X_.cols()) = X_;
    X_env.col(X_.cols()) = env;

    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    grm_opts.thread_num = options_.n_threads;
    #endif
    GrmBuilder grm_builder(grm_opts);
    if (is_pgen_) { pgen_reader_->reset(); grm_builder.compute(*pgen_reader_); }
    else if (is_bgen_) { bgen_reader_->reset(); grm_builder.compute(*bgen_reader_); }
    else { grm_builder.compute(*plink_reader_); }
    const Eigen::MatrixXd& G_global = grm_builder.getGrm();
    Eigen::MatrixXd G_sub(n, n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) G_sub(i, j) = G_global(keep_.keep_indices[i], keep_.keep_indices[j]);

    MvEigenVce::Options vce_opts;
    vce_opts.use_projection = false;
    vce_opts.max_iter = options_.reml_maxit;
    vce_opts.tol = options_.reml_tol;
    vce_opts.solver = options_.use_aireml ? MvEigenVce::VceSolver::AIREML : MvEigenVce::VceSolver::EM;
    vce_opts.em_min_iter = options_.em_min_iter;
    vce_opts.ai_max_iter = options_.ai_max_iter;
    vce_opts.ai_tol = options_.ai_tol;
    MvEigenVce vce(vce_opts);
    vce.prepare(Y_, X_env, G_sub);
    vce.runNullModel();

    Eigen::VectorXd D = vce.getEigenvalues();
    Eigen::MatrixXd U = vce.getProjectedW();
    Eigen::MatrixXd UtY = vce.getUty();
    Eigen::MatrixXd UtXenv = vce.getUtX();
    Eigen::MatrixXd Vg = vce.getVg();
    Eigen::MatrixXd Ve = vce.getVe();
    std::vector<Eigen::MatrixXd> Sigma_inv(n);
    for (int i = 0; i < n; ++i) {
        Eigen::MatrixXd Sigma = D(i) * Vg + Ve;
        Sigma.diagonal().array() += 1e-8;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Sigma);
        if (ldlt.info() != Eigen::Success) throw std::runtime_error("GEMMA-GxE Sigma factorization failed");
        Sigma_inv[i] = ldlt.solve(Eigen::MatrixXd::Identity(q, q));
    }

    const auto& bim = plink_reader_->getBimInfo();
    const size_t m = is_pgen_ ? pgen_reader_->getNumVariants() : is_bgen_ ? bgen_reader_->getNumVariants() : plink_reader_->getNumSnps();
    std::string out_file = options_.out_prefix + ".assoc.mv.gxe.txt";
    auto ofs = mv_open_output_file_gxe(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tWALD_GXE\tP_GXE";
    if (options_.use_gemma_full) ofs << "\tP_WALD_GXE\tP_LRT_GXE\tP_SCORE_GXE";
    for (int t = 0; t < q; ++t) ofs << "\tBETA_GXE_" << (t + 1) << "\tSE_GXE_" << (t + 1);
    if (options_.output_vbeta) mv_write_vbeta_header_gxe(ofs, q);
    ofs << "\n";

    Eigen::MatrixXd U_t = U.transpose();
    const size_t block_size = 512;
    if (is_pgen_) pgen_reader_->reset();
    if (is_bgen_) bgen_reader_->reset();
    for (size_t start = 0; start < m; start += block_size) {
        size_t end = std::min(m, start + block_size);
        Eigen::MatrixXd raw;
        size_t current_block = end - start;
        std::vector<PlinkBimInfo> block_bim(current_block);
        if (is_pgen_) {
            raw = Eigen::MatrixXd(n_samples_, current_block);
            Eigen::VectorXf dosage_f; std::string chrom, rsid, ref, alt; uint32_t pos = 0; size_t k = 0;
            for (; k < current_block; ++k) {
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = (double)dosage_f[i];
                block_bim[k].chrom = chrom; block_bim[k].id = rsid; block_bim[k].cm_pos = 0.0; block_bim[k].bp_pos = (long long)pos; block_bim[k].alt = alt; block_bim[k].ref = ref;
            }
            if (k == 0) break; current_block = k; raw.conservativeResize(Eigen::NoChange, current_block); block_bim.resize(current_block);
        } else if (is_bgen_) {
            raw = Eigen::MatrixXd(n_samples_, current_block);
            Eigen::VectorXf dosage_f; std::string chrom, rsid, a1, a2; uint32_t pos = 0; size_t k = 0;
            for (; k < current_block; ++k) {
                if (!bgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = (double)dosage_f[i];
                block_bim[k].chrom = chrom; block_bim[k].id = rsid; block_bim[k].cm_pos = 0.0; block_bim[k].bp_pos = (long long)pos; block_bim[k].alt = a1; block_bim[k].ref = a2;
            }
            if (k == 0) break; current_block = k; raw.conservativeResize(Eigen::NoChange, current_block); block_bim.resize(current_block);
        } else {
            plink_reader_->getSnpBlock(start, end, raw);
            for (size_t k = 0; k < current_block; ++k) block_bim[k] = bim[start + k];
        }
        for (size_t k = 0; k < current_block; ++k) {
            const auto& snp = block_bim[k];
            if (!mv_scan_filter_passes_gxe(snp, options_)) continue;
            Eigen::VectorXd g(n); double sum = 0.0; int count = 0;
            for (int i = 0; i < n; ++i) { double val = raw(keep_.keep_indices[i], k); g(i) = val; if (std::isfinite(val)) { sum += val; ++count; } }
            if (count == 0) continue;
            double mean = sum / (double)count;
            for (int i = 0; i < n; ++i) { if (!std::isfinite(g(i))) g(i) = mean; if (mean > 1.0) g(i) = 2.0 - g(i); }
            Eigen::VectorXd z = g.array() * env.array();
            Eigen::VectorXd Utg = U_t * g;
            Eigen::VectorXd Utz = U_t * z;
            const int c = UtXenv.cols() + 1;
            const int pd = c * q;
            Eigen::MatrixXd H_CC = Eigen::MatrixXd::Zero(pd, pd);
            Eigen::VectorXd RHS_C = Eigen::VectorXd::Zero(pd);
            Eigen::MatrixXd H_CZ = Eigen::MatrixXd::Zero(pd, q);
            Eigen::MatrixXd H_ZZ = Eigen::MatrixXd::Zero(q, q);
            Eigen::VectorXd RHS_Z = Eigen::VectorXd::Zero(q);
            for (int i = 0; i < n; ++i) {
                const Eigen::MatrixXd& Si = Sigma_inv[i];
                Eigen::VectorXd yi = UtY.row(i).transpose();
                Eigen::VectorXd ci(c);
                ci.head(UtXenv.cols()) = UtXenv.row(i).transpose(); ci(c - 1) = Utg(i);
                double zi = Utz(i); Eigen::VectorXd Siy = Si * yi;
                for (int a = 0; a < c; ++a) {
                    RHS_C.segment(a * q, q) += ci(a) * Siy;
                    H_CZ.block(a * q, 0, q, q) += ci(a) * zi * Si;
                    for (int b = 0; b < c; ++b) H_CC.block(a * q, b * q, q, q) += ci(a) * ci(b) * Si;
                }
                RHS_Z += zi * Siy; H_ZZ += zi * zi * Si;
            }
            H_CC.diagonal().array() += 1e-10;
            Eigen::LDLT<Eigen::MatrixXd> cc_ldlt(H_CC); if (cc_ldlt.info() != Eigen::Success) continue;
            Eigen::MatrixXd HCC_inv_HCZ = cc_ldlt.solve(H_CZ);
            Eigen::VectorXd HCC_inv_RHSC = cc_ldlt.solve(RHS_C);
            Eigen::MatrixXd Schur = H_ZZ - H_CZ.transpose() * HCC_inv_HCZ;
            Eigen::VectorXd adj = RHS_Z - H_CZ.transpose() * HCC_inv_RHSC;
            Schur = 0.5 * (Schur + Schur.transpose()).eval(); Schur.diagonal().array() += 1e-10;
            Eigen::LDLT<Eigen::MatrixXd> s_ldlt(Schur); if (s_ldlt.info() != Eigen::Success) continue;
            Eigen::VectorXd beta = s_ldlt.solve(adj); if (!beta.allFinite()) continue;
            Eigen::MatrixXd vbeta = s_ldlt.solve(Eigen::MatrixXd::Identity(q, q));
            double stat = beta.dot(Schur * beta); if (!std::isfinite(stat) || stat < 0.0) continue;
            double p_score = mv_chisq_prob_gxe(stat, q);
            double p_wald = p_score;
            double p_lrt = std::numeric_limits<double>::quiet_NaN();
            if (options_.use_gemma_full) {
                Eigen::MatrixXd UtX_main(UtXenv.rows(), UtXenv.cols() + 1); UtX_main.leftCols(UtXenv.cols()) = UtXenv; UtX_main.col(UtXenv.cols()) = Utg;
                Eigen::MatrixXd UtX_full(UtXenv.rows(), UtXenv.cols() + 2); UtX_full.leftCols(UtXenv.cols()) = UtXenv; UtX_full.col(UtXenv.cols()) = Utg; UtX_full.col(UtXenv.cols() + 1) = Utz;
                Eigen::MatrixXd Vg_main, Ve_main, Vg_full, Ve_full;
                int em_iter = std::max(1, options_.reml_maxit / 10); double em_tol = std::max(options_.reml_tol * 10.0, 1e-7);
                if (mv_em_refit_rotated_gxe(UtY, UtX_main, D, em_iter, em_tol, Vg_main, Ve_main, &Vg, &Ve)) {
                    Eigen::VectorXd beta_main_all; Eigen::MatrixXd cov_main; double ll_main = std::numeric_limits<double>::quiet_NaN();
                    mv_fit_gls_gxe(UtY, UtX_main, D, Vg_main, Ve_main, beta_main_all, cov_main, &ll_main);
                    if (mv_em_refit_rotated_gxe(UtY, UtX_full, D, em_iter, em_tol, Vg_full, Ve_full, &Vg_main, &Ve_main)) {
                        Eigen::VectorXd beta_full_all; Eigen::MatrixXd cov_full; double ll_full = std::numeric_limits<double>::quiet_NaN();
                        if (mv_fit_gls_gxe(UtY, UtX_full, D, Vg_full, Ve_full, beta_full_all, cov_full, &ll_full)) {
                            int offset = (UtX_full.cols() - 1) * q;
                            beta = beta_full_all.segment(offset, q); vbeta = cov_full.block(offset, offset, q, q);
                            Eigen::MatrixXd vbeta_inv = vbeta.ldlt().solve(Eigen::MatrixXd::Identity(q, q));
                            stat = (beta.transpose() * vbeta_inv * beta)(0, 0); if (!std::isfinite(stat) || stat < 0.0) continue;
                            p_wald = mv_chisq_prob_gxe(stat, q);
                            if (std::isfinite(ll_main) && std::isfinite(ll_full)) p_lrt = mv_chisq_prob_gxe(std::max(0.0, 2.0 * (ll_full - ll_main)), q);
                        }
                    }
                }
            }
            double p = options_.use_gemma_full ? (options_.gemma_mode == 2 ? p_lrt : (options_.gemma_mode == 3 ? p_score : p_wald)) : p_score;
            ofs << snp.chrom << "\t" << snp.id << "\t" << snp.bp_pos << "\t" << snp.alt << "\t" << snp.ref << "\t" << mean / 2.0 << "\t" << stat << "\t" << p;
            if (options_.use_gemma_full) ofs << "\t" << p_wald << "\t" << p_lrt << "\t" << p_score;
            for (int t = 0; t < q; ++t) { double se = vbeta(t, t) > 0.0 ? std::sqrt(vbeta(t, t)) : std::numeric_limits<double>::quiet_NaN(); ofs << "\t" << beta(t) << "\t" << se; }
            if (options_.output_vbeta) mv_write_vbeta_values_gxe(ofs, vbeta, q);
            ofs << "\n";
        }
        std::cout << "\rProcessed " << std::min(start + current_block, m) << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << "\nResults saved to " << out_file << std::endl;
}

} // namespace cosmic
