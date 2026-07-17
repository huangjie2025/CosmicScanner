#include "multivariate_model.h"
#include "multivariate_helpers.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <omp.h>
#include "grm_builder.h"
#include "mvlmm_engine.h"

namespace cosmic {

static Eigen::MatrixXd mv_refit_load_dense_matrix(const std::string& path, int n) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open matrix file: " + path);
    Eigen::MatrixXd M(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(in >> M(i, j))) {
                throw std::runtime_error("Matrix file has fewer values than expected: " + path);
            }
        }
    }
    return M;
}

void MultivariateModel::runGemmaRefit() {
    mv_validate_output_parent(options_.out_prefix);

    if (!options_.no_loco) {
        std::cout << "Multivariate GEMMA-style refit currently runs no-LOCO; use --no-loco for explicit validation.\n";
    }

    std::cout << "Multivariate GWAS method: "
              << (options_.use_gemma_exact ? "GEMMA exact-style mvLMM full alternative refit"
                                            : (options_.use_gemma_full ? "GEMMA-compatible mvLMM Wald/LRT/Score" : "GEMMA-compatible mvLMM Wald"))
              << " (pnr=" << options_.pnr << ", mode=" << options_.gemma_mode << ")\n";
    if (options_.use_gemma_exact) {
        std::cout << "  exact controls: emi=" << options_.gemma_emi
                  << ", nri=" << options_.gemma_nri
                  << ", emp=" << options_.gemma_emp
                  << ", nrp=" << options_.gemma_nrp << "\n";
    }
    if (options_.gemma_crt) {
        std::cout << "  CRT: GEMMA-style Edgeworth/PCRT correction enabled for p < "
                  << options_.pnr << "\n";
    }

    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    grm_opts.thread_num = options_.n_threads;
    #endif

    Eigen::MatrixXd G_sub(n_keep_, n_keep_);
    if (!options_.k_file.empty()) {
        std::cout << "Using external K matrix for GEMMA-style refit: " << options_.k_file << "\n";
        Eigen::MatrixXd K_full = mv_refit_load_dense_matrix(options_.k_file, n_samples_);
        for (int i = 0; i < n_keep_; ++i) {
            for (int j = 0; j < n_keep_; ++j) {
                G_sub(i, j) = K_full(keep_.keep_indices[i], keep_.keep_indices[j]);
            }
        }
        G_sub = 0.5 * (G_sub + G_sub.transpose());
    } else if (refit_grm_override_) {
        std::cout << "Using LOCO GRM override for GEMMA-style refit"
                  << (refit_chrom_filter_.empty() ? "" : (" (chr " + refit_chrom_filter_ + ")"))
                  << "...\n";
        G_sub = *refit_grm_override_;
    } else {
        std::cout << "Calculating Global GRM for GEMMA-style refit...\n";
        GrmBuilder global_builder(grm_opts);
        if (is_pgen_) { pgen_reader_->reset(); global_builder.compute(*pgen_reader_); }
        else if (is_bgen_) { bgen_reader_->reset(); global_builder.compute(*bgen_reader_); }
        else { global_builder.compute(*plink_reader_); }
        const Eigen::MatrixXd& G_global = global_builder.getGrm();

        for (int i = 0; i < n_keep_; ++i) {
            for (int j = 0; j < n_keep_; ++j) {
                G_sub(i, j) = G_global(keep_.keep_indices[i], keep_.keep_indices[j]);
            }
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G_sub);
    Eigen::VectorXd D = es.eigenvalues();
    for (int i = 0; i < D.size(); ++i) if (D(i) < 0.0) D(i) = 0.0;
    Eigen::MatrixXd U = es.eigenvectors();
    Eigen::MatrixXd UtY = U.transpose() * Y_;
    Eigen::MatrixXd UtX_base = U.transpose() * X_;

    MvEigenVce::Options null_vce_opts;
    null_vce_opts.use_projection = false;
    null_vce_opts.max_iter = options_.reml_maxit;
    null_vce_opts.tol = options_.reml_tol;
    null_vce_opts.verbose = false;
    null_vce_opts.use_reml_correction = true;
    null_vce_opts.use_dense_reml_refinement = true;
    if (options_.use_aireml) {
        null_vce_opts.solver = MvEigenVce::VceSolver::AIREML;
        null_vce_opts.em_min_iter = options_.em_min_iter;
        null_vce_opts.ai_max_iter = options_.ai_max_iter;
        null_vce_opts.ai_tol = options_.ai_tol;
    }
    MvEigenVce null_vce(null_vce_opts);
    null_vce.prepare(Y_, X_, G_sub);
    null_vce.runNullModel();
    Eigen::MatrixXd Vg_null = null_vce.getVg();
    Eigen::MatrixXd Ve_null = null_vce.getVe();
    Eigen::VectorXd beta_null;
    Eigen::MatrixXd cov_null;
    double loglik_null = std::numeric_limits<double>::quiet_NaN();
    mv_fit_gls(UtY, UtX_base, D, Vg_null, Ve_null, beta_null, cov_null, &loglik_null);

    std::string out_file = options_.out_prefix + (options_.use_hybrid ? ".refined.assoc.mv.txt" : ".assoc.mv.txt");
    if (!refit_append_mode_) {
        auto ofs = mv_open_output_file(out_file);
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tWALD_STAT\tP";
        if (options_.use_gemma_full) {
            ofs << "\tP_WALD\tP_LRT\tP_SCORE\tP_WALD_GC\tP_LRT_GC\tP_SCORE_GC";
            if (options_.gemma_crt) {
                ofs << "\tP_WALD_CRT\tP_LRT_CRT\tP_SCORE_CRT";
            }
        }
        for (int t = 0; t < Y_.cols(); ++t) ofs << "\tBETA_" << (t + 1) << "\tSE_" << (t + 1);
        if (options_.output_vbeta) mv_write_vbeta_header(ofs, Y_.cols());
        ofs << "\n";
    }

    const auto& bim = plink_reader_->getBimInfo();
    const size_t m = is_pgen_ ? pgen_reader_->getNumVariants()
                   : is_bgen_ ? bgen_reader_->getNumVariants()
                              : plink_reader_->getNumSnps();
    const size_t block_size = 256;
    std::vector<MvGwasResult> results;
    results.reserve(m);

    if (is_pgen_) pgen_reader_->reset();
    if (is_bgen_) bgen_reader_->reset();

    for (size_t start = 0; start < m; start += block_size) {
        size_t end = std::min(start + block_size, m);
        Eigen::MatrixXd raw;
        size_t current_block = end - start;
        std::vector<PlinkBimInfo> block_bim(current_block);
        if (is_pgen_) {
            raw = Eigen::MatrixXd(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, rsid, ref, alt;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = (double)dosage_f[i];
                block_bim[k].chrom = chrom;
                block_bim[k].id = rsid;
                block_bim[k].cm_pos = 0.0;
                block_bim[k].bp_pos = (long long)pos;
                block_bim[k].alt = alt;
                block_bim[k].ref = ref;
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else if (is_bgen_) {
            raw = Eigen::MatrixXd(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, rsid, a1, a2;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!bgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = (double)dosage_f[i];
                block_bim[k].chrom = chrom;
                block_bim[k].id = rsid;
                block_bim[k].cm_pos = 0.0;
                block_bim[k].bp_pos = (long long)pos;
                block_bim[k].alt = a1;
                block_bim[k].ref = a2;
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else {
            plink_reader_->getSnpBlock(start, end, raw);
            for (size_t k = 0; k < current_block; ++k) block_bim[k] = bim[start + k];
        }

        std::vector<MvGwasResult> block_results(current_block);
        std::vector<char> block_valid(current_block, 0);

        #pragma omp parallel for schedule(dynamic)
        for (long long kk_ll = 0; kk_ll < (long long)current_block; ++kk_ll) {
            size_t kk = (size_t)kk_ll;
            size_t snp_idx = start + kk;
            if (!refine_mask_.empty() && !refine_mask_[snp_idx]) continue;
            if (!mv_scan_filter_passes(block_bim[kk], options_)) continue;
            if (!refit_chrom_filter_.empty() && block_bim[kk].chrom != refit_chrom_filter_) continue;
            Eigen::VectorXd x(n_keep_);
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < n_keep_; ++i) {
                double v = raw(keep_.keep_indices[i], kk);
                if (!std::isnan(v)) {
                    sum += v;
                    ++count;
                }
            }
            double mean = count > 0 ? sum / count : 0.0;
            for (int i = 0; i < n_keep_; ++i) {
                double v = raw(keep_.keep_indices[i], kk);
                x(i) = std::isnan(v) ? mean : v;
            }

            Eigen::MatrixXd X_aug(n_keep_, X_.cols() + 1);
            X_aug.leftCols(X_.cols()) = X_;
            X_aug.col(X_.cols()) = x;

            Eigen::MatrixXd UtX_aug(UtX_base.rows(), UtX_base.cols() + 1);
            UtX_aug.leftCols(UtX_base.cols()) = UtX_base;
            UtX_aug.col(UtX_base.cols()) = U.transpose() * x;

            int n_traits = Y_.cols();
            int n_covars = UtX_aug.cols();
            int offset = (n_covars - 1) * n_traits;

            Eigen::VectorXd beta;
            Eigen::MatrixXd Vbeta;
            double wald_stat = std::numeric_limits<double>::quiet_NaN();
            double p_wald = std::numeric_limits<double>::quiet_NaN();
            double p_lrt = std::numeric_limits<double>::quiet_NaN();
            double p_score_val = std::numeric_limits<double>::quiet_NaN();
            double p_wald_crt = std::numeric_limits<double>::quiet_NaN();
            double p_lrt_crt = std::numeric_limits<double>::quiet_NaN();
            double p_score_crt = std::numeric_limits<double>::quiet_NaN();

            Eigen::VectorXd beta_score_all;
            Eigen::MatrixXd cov_score_all;
            bool have_score = mv_fit_gls(UtY, UtX_aug, D, Vg_null, Ve_null,
                                         beta_score_all, cov_score_all, nullptr);
            if (!have_score) continue;
            Eigen::VectorXd beta_score = beta_score_all.segment(offset, n_traits);
            Eigen::MatrixXd Vbeta_score = cov_score_all.block(offset, offset, n_traits, n_traits);
            Eigen::MatrixXd Vs_inv = Vbeta_score.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
            double score_stat = (beta_score.transpose() * Vs_inv * beta_score)(0, 0);
            double p_score_calc = std::numeric_limits<double>::quiet_NaN();
            if (std::isfinite(score_stat) && score_stat >= 0.0) p_score_calc = mv_chisq_prob(score_stat, n_traits);
            if (options_.test_score) p_score_val = p_score_calc;
            p_score_crt = p_score_val;

            beta = beta_score;
            Vbeta = Vbeta_score;
            wald_stat = score_stat;

            if (!options_.use_gemma_full) {
                Eigen::MatrixXd Vg_wald, Ve_wald;
                int snp_refit_maxit = std::max(30, options_.reml_maxit / 20);
                double snp_refit_tol = std::max(options_.reml_tol * 10.0, 1e-7);
                if (!mv_em_refit_rotated(UtY, UtX_aug, D, snp_refit_maxit, snp_refit_tol,
                                         Vg_wald, Ve_wald, &Vg_null, &Ve_null)) continue;
                Eigen::VectorXd beta_all;
                Eigen::MatrixXd Hinv;
                if (!mv_fit_gls(UtY, UtX_aug, D, Vg_wald, Ve_wald, beta_all, Hinv, nullptr)) continue;
                beta = beta_all.segment(offset, n_traits);
                Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                Eigen::MatrixXd Vbeta_inv = Vbeta.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                wald_stat = (beta.transpose() * Vbeta_inv * beta)(0, 0);
                if (!std::isfinite(wald_stat) || wald_stat < 0.0) continue;
                p_wald = mv_chisq_prob(wald_stat, n_traits);
                p_wald_crt = p_wald;
            } else {
                if (options_.gemma_crt && options_.test_score &&
                    std::isfinite(p_score_val) && p_score_val < options_.pnr) {
                    MvCrtFactors crt = mv_compute_crt_factors(UtX_aug, D, Vg_null, Ve_null);
                    p_score_crt = mv_pcrt(3, n_traits, p_score_val, crt);
                }
                if (options_.test_lrt) {
                    Eigen::MatrixXd Vg_lrt, Ve_lrt;
                    int em_iter = options_.use_gemma_exact ? std::max(1, options_.gemma_emi)
                                                           : std::max(1, options_.reml_maxit / 10);
                    double em_tol = options_.use_gemma_exact ? options_.gemma_emp
                                                             : std::max(options_.reml_tol * 10.0, 1e-7);
                    if (mv_em_refit_rotated(UtY, UtX_aug, D, em_iter, em_tol,
                                            Vg_lrt, Ve_lrt, &Vg_null, &Ve_null)) {
                        Eigen::VectorXd beta_all;
                        Eigen::MatrixXd Hinv;
                        double loglik_alt = std::numeric_limits<double>::quiet_NaN();
                        if (mv_fit_gls(UtY, UtX_aug, D, Vg_lrt, Ve_lrt, beta_all, Hinv, &loglik_alt) &&
                            std::isfinite(loglik_alt) && std::isfinite(loglik_null)) {
                            double lrt = std::max(0.0, 2.0 * (loglik_alt - loglik_null));
                            p_lrt = mv_chisq_prob(lrt, n_traits);
                            p_lrt_crt = p_lrt;
                            if (std::isfinite(p_lrt) && p_lrt < options_.pnr) {
                                Eigen::MatrixXd Vg_lrt2, Ve_lrt2;
                                int nr_like_iter = options_.use_gemma_exact ? std::max(1, options_.gemma_nri)
                                                                            : std::max(1, options_.ai_max_iter);
                                double nr_like_tol = options_.use_gemma_exact ? options_.gemma_nrp
                                                                              : options_.ai_tol;
                                if (mv_em_refit_rotated(UtY, UtX_aug, D, nr_like_iter, nr_like_tol,
                                                        Vg_lrt2, Ve_lrt2, &Vg_lrt, &Ve_lrt) &&
                                    mv_fit_gls(UtY, UtX_aug, D, Vg_lrt2, Ve_lrt2, beta_all, Hinv, &loglik_alt) &&
                                    std::isfinite(loglik_alt)) {
                                    lrt = std::max(0.0, 2.0 * (loglik_alt - loglik_null));
                                    p_lrt = mv_chisq_prob(lrt, n_traits);
                                    Vg_lrt = Vg_lrt2;
                                    Ve_lrt = Ve_lrt2;
                                    p_lrt_crt = p_lrt;
                                }
                            }
                            if (options_.gemma_crt && std::isfinite(p_lrt) && p_lrt < options_.pnr) {
                                MvCrtFactors crt = mv_compute_crt_factors(UtX_aug, D, Vg_lrt, Ve_lrt);
                                p_lrt_crt = mv_pcrt(2, n_traits, p_lrt, crt);
                            }
                            if (!options_.test_wald) {
                                beta = beta_all.segment(offset, n_traits);
                                Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                            }
                        }
                    }
                }

                if (options_.test_wald) {
                    Eigen::MatrixXd Vg_wald, Ve_wald;
                    int em_iter = options_.use_gemma_exact ? std::max(1, options_.gemma_emi)
                                                           : std::max(1, options_.reml_maxit / 10);
                    double em_tol = options_.use_gemma_exact ? options_.gemma_emp
                                                             : std::max(options_.reml_tol * 10.0, 1e-7);
                    if (!mv_em_refit_rotated(UtY, UtX_aug, D, em_iter, em_tol,
                                             Vg_wald, Ve_wald, &Vg_null, &Ve_null)) continue;
                    Eigen::VectorXd beta_all;
                    Eigen::MatrixXd Hinv;
                    if (!mv_fit_gls(UtY, UtX_aug, D, Vg_wald, Ve_wald, beta_all, Hinv, nullptr)) continue;
                    beta = beta_all.segment(offset, n_traits);
                    Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                    Eigen::MatrixXd Vbeta_inv = Vbeta.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                    wald_stat = (beta.transpose() * Vbeta_inv * beta)(0, 0);
                    if (!std::isfinite(wald_stat) || wald_stat < 0.0) continue;
                    p_wald = mv_chisq_prob(wald_stat, n_traits);
                    p_wald_crt = p_wald;
                    if (std::isfinite(p_wald) && p_wald < options_.pnr) {
                        Eigen::MatrixXd Vg_wald2, Ve_wald2;
                        int nr_like_iter = options_.use_gemma_exact ? std::max(1, options_.gemma_nri)
                                                                    : std::max(1, options_.ai_max_iter);
                        double nr_like_tol = options_.use_gemma_exact ? options_.gemma_nrp
                                                                      : options_.ai_tol;
                        if (mv_em_refit_rotated(UtY, UtX_aug, D, nr_like_iter, nr_like_tol,
                                                Vg_wald2, Ve_wald2, &Vg_wald, &Ve_wald) &&
                            mv_fit_gls(UtY, UtX_aug, D, Vg_wald2, Ve_wald2, beta_all, Hinv, nullptr)) {
                            beta = beta_all.segment(offset, n_traits);
                            Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                            Vbeta_inv = Vbeta.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                            wald_stat = (beta.transpose() * Vbeta_inv * beta)(0, 0);
                            if (std::isfinite(wald_stat) && wald_stat >= 0.0) p_wald = mv_chisq_prob(wald_stat, n_traits);
                            p_wald_crt = p_wald;
                            Vg_wald = Vg_wald2;
                            Ve_wald = Ve_wald2;
                        }
                    }
                    if (options_.gemma_crt && std::isfinite(p_wald) && p_wald < options_.pnr) {
                        MvCrtFactors crt = mv_compute_crt_factors(UtX_aug, D, Vg_wald, Ve_wald);
                        p_wald_crt = mv_pcrt(1, n_traits, p_wald, crt);
                    }
                }
            }

            MvGwasResult res;
            res.chrom = block_bim[kk].chrom;
            res.snp_id = block_bim[kk].id;
            res.pos = block_bim[kk].bp_pos;
            res.a1 = block_bim[kk].alt;
            res.a2 = block_bim[kk].ref;
            res.freq = mean / 2.0;
            res.beta = beta;
            res.se = Vbeta;
            res.wald_stat = wald_stat;
            if (options_.use_gemma_full && options_.gemma_mode == 2) res.p_value = options_.gemma_crt ? p_lrt_crt : p_lrt;
            else if (options_.use_gemma_full && options_.gemma_mode == 3) res.p_value = options_.gemma_crt ? p_score_crt : p_score_val;
            else res.p_value = options_.gemma_crt && std::isfinite(p_wald_crt) ? p_wald_crt
                               : (std::isfinite(p_wald) ? p_wald : (std::isfinite(p_lrt) ? p_lrt : p_score_val));
            res.p_wald = p_wald;
            res.p_lrt = p_lrt;
            res.p_score = p_score_val;
            res.p_wald_crt = p_wald_crt;
            res.p_lrt_crt = p_lrt_crt;
            res.p_score_crt = p_score_crt;
            block_results[kk] = std::move(res);
            block_valid[kk] = 1;
        }

        for (size_t kk = 0; kk < current_block; ++kk) {
            if (block_valid[kk]) results.push_back(std::move(block_results[kk]));
        }
        std::cout << "\rProcessed " << end << "/" << m << " SNPs..." << std::flush;
    }
    std::cout << "\n";
    if (options_.use_gemma_full) {
        auto median_lambda = [&](double MvGwasResult::*field) {
            std::vector<double> stats;
            stats.reserve(results.size());
            for (const auto& res : results) {
                double p = res.*field;
                if (std::isfinite(p) && p > 0.0 && p < 1.0) {
                    stats.push_back(mv_chisq_stat_from_p(p, Y_.cols()));
                }
            }
            if (stats.empty()) return 1.0;
            std::sort(stats.begin(), stats.end());
            double median = stats[stats.size() / 2];
            double expected = mv_chisq_stat_from_p(0.5, Y_.cols());
            if (!(expected > 0.0)) return 1.0;
            return std::max(1.0, median / expected);
        };
        double lambda_wald = median_lambda(&MvGwasResult::p_wald);
        double lambda_lrt = median_lambda(&MvGwasResult::p_lrt);
        double lambda_score = median_lambda(&MvGwasResult::p_score);

        auto ofs = mv_open_output_file(out_file, std::ios::out | std::ios::app);
        for (const auto& res : results) {
            if (std::isnan(res.p_value)) continue;
            double p_wald_gc = mv_gc_correct_p(res.p_wald, lambda_wald, Y_.cols());
            double p_lrt_gc = mv_gc_correct_p(res.p_lrt, lambda_lrt, Y_.cols());
            double p_score_gc = mv_gc_correct_p(res.p_score, lambda_score, Y_.cols());
            ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
                << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
                << res.wald_stat << "\t" << res.p_value << "\t"
                << res.p_wald << "\t" << res.p_lrt << "\t" << res.p_score << "\t"
                << p_wald_gc << "\t" << p_lrt_gc << "\t" << p_score_gc;
            if (options_.gemma_crt) {
                ofs << "\t" << res.p_wald_crt << "\t" << res.p_lrt_crt << "\t" << res.p_score_crt;
            }
            for (int t = 0; t < Y_.cols(); ++t) {
                ofs << "\t" << res.beta(t) << "\t" << std::sqrt(std::max(0.0, res.se(t, t)));
            }
            if (options_.output_vbeta) mv_write_vbeta_values(ofs, res.se, Y_.cols());
            ofs << "\n";
        }
    } else {
        writeResults(out_file, results);
    }
    std::cout << "Results saved to " << out_file << std::endl;
}


} // namespace cosmic
