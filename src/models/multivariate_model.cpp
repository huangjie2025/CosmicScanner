#include "models/multivariate_model.h"
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
#include <omp.h>
#include "pheno_reader.h"
#include "covar_reader.h"
#include "data_filter.h"
#include "data_loader.h"
#include "grm_builder.h"
#include "mvlmm_engine.h"

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

static double mv_chisq_prob(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return mv_gamma_q(df / 2.0, stat / 2.0);
}

static bool mv_fit_gls(const Eigen::MatrixXd& UtY,
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

static bool mv_em_refit_rotated(const Eigen::MatrixXd& UtY,
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

    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd beta;
        Eigen::MatrixXd beta_cov;
        if (!mv_fit_gls(UtY, UtX, D, Vg, Ve, beta, beta_cov, nullptr)) return false;

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

            Eigen::MatrixXd Wg = Vg * D(i) * Pinv;
            Eigen::VectorXd uhat = Wg * resid;
            Sg += uhat * uhat.transpose() + Vg * D(i) - Wg * Vg * D(i);

            Eigen::MatrixXd We = Ve * Pinv;
            Eigen::VectorXd ehat = We * resid;
            Se += ehat * ehat.transpose() + Ve - We * Ve;
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
        if (err < tol) return true;
    }
    return true;
}

static double mv_chisq_stat_from_p(double p, int df) {
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

static double mv_gc_correct_p(double p, double lambda, int df) {
    if (!std::isfinite(p) || lambda <= 1.0) return p;
    double stat = mv_chisq_stat_from_p(p, df);
    return mv_chisq_prob(stat / lambda, df);
}

static bool contains_mv_test(const std::string& tests, const std::string& name) {
    std::stringstream ss(tests);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token == name || token == "all") return true;
    }
    return false;
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

static bool mv_scan_filter_passes(const PlinkBimInfo& var,
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
    opts.out_prefix = cfg.out_prefix;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.grm_alg = cfg.grm_alg;
    if (!cfg.multivariate_compatible.empty()) {
        if (cfg.multivariate_compatible == "emmvgwas") {
            opts.grm_alg = 6;
            opts.use_covariate_projection = false;
            opts.use_emmvgwas_lambda = true;
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
    opts.use_gemma_full = false;
    opts.use_hybrid = false;
    if (cfg.multivariate_method == "wald") {
        opts.use_emmv = false;
    } else if (cfg.multivariate_method == "gemma") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
    } else if (cfg.multivariate_method == "gemma-full") {
        opts.use_emmv = false;
        opts.use_gemma_refit = true;
        opts.use_gemma_full = true;
    } else if (cfg.multivariate_method == "hybrid") {
        opts.use_emmv = true;
        opts.use_hybrid = true;
    } else if (cfg.multivariate_method == "eb" || cfg.multivariate_method == "eb-lmm" || cfg.multivariate_method == "lmm-eb") {
        opts.use_eb = true;
        opts.use_emmv = false;
    } else if (cfg.multivariate_method != "emmv") {
        throw std::runtime_error("Unknown --mv-method: " + cfg.multivariate_method + ". Supported values: emmv, wald, eb, gemma, gemma-full, hybrid");
    }
    if (cfg.multivariate_eb || cfg.eb_lmm) { opts.use_eb = true; opts.use_emmv = false; }
    opts.test_wald = contains_mv_test(cfg.multivariate_tests, "wald");
    opts.test_lrt = contains_mv_test(cfg.multivariate_tests, "lrt");
    opts.test_score = contains_mv_test(cfg.multivariate_tests, "score");
    opts.reml_maxit = std::max(cfg.reml_maxit, (opts.use_emmv || opts.use_gemma_refit) ? 1000 : 100);
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
    opts.pheno_file = cfg.pheno_file;
    opts.mpheno_list = cfg.mpheno_list;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
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
        if (keep) keep_.keep_indices.push_back(i);
    }
    keep_.n_keep = keep_.keep_indices.size();
    n_keep_ = keep_.n_keep;
    std::cout << "After filtering missing data, " << n_keep_ << " individuals remain.\n";
    if (n_keep_ < 10) throw std::runtime_error("Too few valid samples.");

    Eigen::MatrixXd Y_sub(n_keep_, Y_.cols());
    Eigen::MatrixXd X_sub(n_keep_, X_.cols());
    for (int i = 0; i < n_keep_; ++i) {
        Y_sub.row(i) = Y_.row(keep_.keep_indices[i]);
        X_sub.row(i) = X_.row(keep_.keep_indices[i]);
    }
    Y_ = std::move(Y_sub);
    X_ = std::move(X_sub);
}

void MultivariateModel::run() {
    if (options_.use_hybrid) {
        runHybrid();
        return;
    }

    if (options_.use_gemma_refit) {
        runGemmaRefit();
        return;
    }

    std::cout << "Multivariate GWAS method: "
              << (options_.use_eb ? "experimental EB-mvLMM random-SNP Wald" : (options_.use_emmv ? "EMmvGWAS semi-exact Wald" : "fixed-null Wald"))
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

    std::string out_file = options_.out_prefix + ".assoc.mv.txt";
    {
        std::ofstream ofs(out_file);
        eb_effective_tests_ = 0.0;
        eb_tested_markers_ = 0;
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\t" << (options_.use_eb ? "EB_WALD" : "WALD_STAT") << "\tP";
        int n_traits = Y_.cols();
        for (int t = 0; t < n_traits; ++t) {
            ofs << "\t" << (options_.use_eb ? "EB_BETA_" : "BETA_") << (t + 1)
                << "\t" << (options_.use_eb ? "EB_SE_" : "SE_") << (t + 1);
        }
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
        vce_opts.use_emmvgwas_lambda = options_.use_emmvgwas_lambda;
        MvEigenVce vce(vce_opts);
        vce.prepare(Y_, X_, G_sub);
        vce.runNullModel();

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

void MultivariateModel::runHybrid() {
    if (is_bgen_ || is_pgen_) {
        throw std::runtime_error("--mv-method hybrid currently supports PLINK --bfile input only.");
    }
    std::cout << "Multivariate GWAS method: "
              << (options_.no_loco ? "hybrid EMmvGWAS scan + GEMMA-full top-SNP refinement"
                                   : "LOCO-aware hybrid EMmvGWAS scan + GEMMA-full top-SNP refinement")
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

    std::cout << "Hybrid stage 1: " << (options_.no_loco ? "EMmvGWAS full scan" : "LOCO-aware EMmvGWAS full scan") << "...\n";
    GrmBuilder global_builder(grm_opts);
    global_builder.compute(*plink_reader_);
    const Eigen::MatrixXd& G_global = global_builder.getGrm();
    double S_global = global_builder.getSum2pq();
    auto chroms_set = global_builder.getObservedChroms();

    Eigen::MatrixXd G_global_sub(n_keep_, n_keep_);
    for (int i = 0; i < n_keep_; ++i) {
        for (int j = 0; j < n_keep_; ++j) {
            G_global_sub(i, j) = G_global(keep_.keep_indices[i], keep_.keep_indices[j]);
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
    MvLmmEngine::Options mv_opts;
    mv_opts.use_projection = options_.use_covariate_projection;
    mv_opts.em_mvgwas_test = true;
    mv_opts.use_general_sigma_inverse = options_.use_emmvgwas_lambda;
    mv_opts.center_genotypes = !options_.use_emmvgwas_lambda;
    mv_opts.use_dosage_minus_one = options_.use_emmvgwas_lambda;
    mv_opts.scan_chrom = options_.scan_chrom;
    mv_opts.scan_from_bp = options_.scan_from_bp;
    mv_opts.scan_to_bp = options_.scan_to_bp;
    mv_opts.extract_set = options_.extract_set;
    mv_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    mv_opts.thread_num = options_.n_threads;
    #endif
    std::string out_file = options_.out_prefix + ".assoc.mv.txt";
    {
        std::ofstream ofs(out_file);
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tWALD_STAT\tP";
        for (int t = 0; t < Y_.cols(); ++t) ofs << "\tBETA_" << (t + 1) << "\tSE_" << (t + 1);
        ofs << "\n";
    }

    std::vector<MvGwasResult> results(plink_reader_->getNumSnps());
    std::map<std::string, Eigen::MatrixXd> loco_refit_grms;
    if (options_.no_loco || chroms_set.size() <= 1) {
        MvEigenVce vce(vce_opts);
        vce.prepare(Y_, X_, G_global_sub);
        vce.runNullModel();

        MvLmmEngine engine(mv_opts);
        engine.prepareFromVce(vce);
        engine.scan(*plink_reader_, results, keep_.keep_indices);
        writeResults(out_file, results);
    } else {
        std::vector<std::string> sorted_chroms(chroms_set.begin(), chroms_set.end());
        std::sort(sorted_chroms.begin(), sorted_chroms.end(), [](const std::string& a, const std::string& b) {
            try { return std::stoi(a) < std::stoi(b); } catch (...) { return a < b; }
        });

        for (const std::string& chrom : sorted_chroms) {
            if (!options_.scan_chrom.empty() && chrom != options_.scan_chrom) continue;
            std::cout << "Hybrid LOCO screening chromosome " << chrom << "...\n";

            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            chr_builder.compute(*plink_reader_);

            double S_chr = chr_builder.getSum2pq();
            double S_loco = S_global - S_chr;
            const auto& G_chr = chr_builder.getGrm();
            Eigen::MatrixXd G_loco_sub = G_global_sub;
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                double w2 = S_chr / S_loco;
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < n_keep_; ++i) {
                    for (int j = 0; j < n_keep_; ++j) {
                        int gi = keep_.keep_indices[i];
                        int gj = keep_.keep_indices[j];
                        G_loco_sub(i, j) = G_global(gi, gj) * w1 - G_chr(gi, gj) * w2;
                    }
                }
            }
            loco_refit_grms[chrom] = G_loco_sub;

            MvEigenVce vce(vce_opts);
            vce.prepare(Y_, X_, G_loco_sub);
            vce.runNullModel();

            MvLmmEngine engine(mv_opts);
            engine.prepareFromVce(vce);
            std::vector<MvGwasResult> chrom_results;
            engine.scan(*plink_reader_, chrom_results, keep_.keep_indices, chrom);
            writeResults(out_file, chrom_results);
            for (size_t i = 0; i < chrom_results.size() && i < results.size(); ++i) {
                if (!chrom_results[i].snp_id.empty() && std::isfinite(chrom_results[i].p_value)) {
                    results[i] = std::move(chrom_results[i]);
                }
            }
        }
    }

    std::vector<std::pair<double, size_t>> valid_by_p;
    valid_by_p.reserve(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        if (std::isfinite(results[i].p_value) && !results[i].snp_id.empty()) {
            valid_by_p.emplace_back(results[i].p_value, i);
        }
    }
    std::sort(valid_by_p.begin(), valid_by_p.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    std::unordered_map<std::string, double> prepass_score_p;
    size_t score_prepass_candidates = 0;
    if (options_.refine_score_prepass > 0 && options_.refine_discrepancy_top > 0 && !valid_by_p.empty()) {
        score_prepass_candidates = std::min(valid_by_p.size(), (size_t)options_.refine_score_prepass);
        refine_mask_.assign(plink_reader_->getNumSnps(), 0);
        for (size_t j = 0; j < score_prepass_candidates; ++j) {
            refine_mask_[valid_by_p[j].second] = 1;
        }

        bool old_full = options_.use_gemma_full;
        bool old_refit = options_.use_gemma_refit;
        bool old_no_loco = options_.no_loco;
        double old_refit_p = options_.refit_p;
        std::string old_out_prefix = options_.out_prefix;
        options_.use_gemma_full = true;
        options_.use_gemma_refit = true;
        options_.no_loco = true;
        options_.refit_p = 0.0;
        options_.out_prefix = old_out_prefix + ".score_prepass";

        std::cout << "Hybrid Stage 1.5 Score-only prepass on top " << score_prepass_candidates << " SNPs...\n";
        if (!loco_refit_grms.empty()) {
            bool first_refit = true;
            for (auto& entry : loco_refit_grms) {
                refit_grm_override_ = &entry.second;
                refit_chrom_filter_ = entry.first;
                refit_append_mode_ = !first_refit;
                runGemmaRefit();
                first_refit = false;
            }
            refit_grm_override_ = nullptr;
            refit_chrom_filter_.clear();
            refit_append_mode_ = false;
        } else {
            runGemmaRefit();
        }

        prepass_score_p = mv_read_score_prepass(options_.out_prefix + ".refined.assoc.mv.txt");
        options_.out_prefix = old_out_prefix;
        options_.use_gemma_full = old_full;
        options_.use_gemma_refit = old_refit;
        options_.no_loco = old_no_loco;
        options_.refit_p = old_refit_p;
        refine_mask_.clear();
        std::cout << "Hybrid Score-only prepass produced " << prepass_score_p.size() << " score p-values.\n";
    }

    std::vector<char> selected(results.size(), 0), selected_by_p(results.size(), 0), selected_by_top(results.size(), 0), selected_by_discrepancy(results.size(), 0);
    size_t candidates_by_p = 0;
    size_t candidates_by_top = 0;
    bool p_threshold_enabled = options_.refine_p < 1.0;
    for (const auto& item : valid_by_p) {
        if (p_threshold_enabled && item.first <= options_.refine_p) {
            selected[item.second] = 1;
            selected_by_p[item.second] = 1;
            ++candidates_by_p;
        }
    }
    if (options_.refine_top > 0) {
        size_t top_n = std::min(valid_by_p.size(), (size_t)options_.refine_top);
        for (size_t j = 0; j < top_n; ++j) {
            selected[valid_by_p[j].second] = 1;
            selected_by_top[valid_by_p[j].second] = 1;
        }
        candidates_by_top = top_n;
    }

    auto neg_log10_p = [](double p) {
        if (!std::isfinite(p) || p <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        return -std::log10(std::max(p, 1e-300));
    };
    std::vector<std::pair<double, size_t>> valid_by_discrepancy;
    valid_by_discrepancy.reserve(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        if (results[i].snp_id.empty()) continue;
        double p_wald = std::isfinite(results[i].p_wald) ? results[i].p_wald : results[i].p_value;
        double p_score = results[i].p_score;
        auto prepass_it = prepass_score_p.find(results[i].snp_id);
        if (prepass_it != prepass_score_p.end()) p_score = prepass_it->second;
        double lw = neg_log10_p(p_wald);
        double ls = neg_log10_p(p_score);
        if (std::isfinite(lw) && std::isfinite(ls)) {
            double discrepancy = std::abs(ls - lw);
            if (discrepancy >= options_.refine_discrepancy) {
                valid_by_discrepancy.emplace_back(discrepancy, i);
            }
        }
    }
    std::sort(valid_by_discrepancy.begin(), valid_by_discrepancy.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    size_t candidates_by_discrepancy = 0;
    if (options_.refine_discrepancy_top > 0) {
        size_t disc_n = std::min(valid_by_discrepancy.size(), (size_t)options_.refine_discrepancy_top);
        for (size_t j = 0; j < disc_n; ++j) {
            selected[valid_by_discrepancy[j].second] = 1;
            selected_by_discrepancy[valid_by_discrepancy[j].second] = 1;
        }
        candidates_by_discrepancy = disc_n;
    }

    std::vector<std::pair<double, size_t>> candidates;
    candidates.reserve(valid_by_p.size());
    for (const auto& item : valid_by_p) {
        if (selected[item.second]) candidates.push_back(item);
    }
    size_t candidates_before_top = candidates_by_p;

    refine_mask_.assign(plink_reader_->getNumSnps(), 0);
    std::string snplist = options_.out_prefix + ".refine.snplist";
    std::string summary_file = options_.out_prefix + ".hybrid.summary.tsv";
    {
        std::ofstream ofs(snplist);
        for (const auto& item : candidates) {
            refine_mask_[item.second] = 1;
            ofs << results[item.second].snp_id << "\t" << item.first << "\n";
        }
    }
    {
        std::ofstream ofs(summary_file);
        ofs << "stage1_mode\t" << ((options_.no_loco || chroms_set.size() <= 1) ? "global" : "loco") << "\n";
        ofs << "screening_method\tEMmvGWAS semi-exact Wald\n";
        ofs << "refinement_method\t" << (loco_refit_grms.empty() ? "GEMMA-full Wald/LRT/Score"
                                                                 : "per-chrom LOCO GEMMA-full Wald/LRT/Score") << "\n";
        ofs << "refine_p\t" << options_.refine_p << "\n";
        ofs << "refine_p_threshold_enabled\t" << (p_threshold_enabled ? "yes" : "no") << "\n";
        ofs << "refine_top\t" << options_.refine_top << "\n";
        ofs << "refine_discrepancy\t" << options_.refine_discrepancy << "\n";
        ofs << "refine_discrepancy_top\t" << options_.refine_discrepancy_top << "\n";
        ofs << "refine_score_prepass\t" << options_.refine_score_prepass << "\n";
        ofs << "score_prepass_candidates\t" << score_prepass_candidates << "\n";
        ofs << "score_prepass_scores\t" << prepass_score_p.size() << "\n";
        ofs << "discrepancy_available\t" << (!prepass_score_p.empty() || !valid_by_discrepancy.empty() ? "yes" : "no") << "\n";
        ofs << "candidates_by_p\t" << candidates_by_p << "\n";
        ofs << "candidates_by_top\t" << candidates_by_top << "\n";
        ofs << "candidates_by_discrepancy\t" << candidates_by_discrepancy << "\n";
        ofs << "candidates_before_top_filter\t" << candidates_before_top << "\n";
        ofs << "candidates_selected\t" << candidates.size() << "\n";
        ofs << "snplist\t" << snplist << "\n";
    }
    std::cout << "Hybrid selected " << candidates.size() << " SNPs for GEMMA-full refinement: " << snplist << "\n";
    if (candidates.empty()) {
        std::cout << "Hybrid refinement skipped because no SNP passed --mv-refine-p.\n";
        return;
    }

    bool old_full = options_.use_gemma_full;
    bool old_refit = options_.use_gemma_refit;
    bool old_no_loco = options_.no_loco;
    options_.use_gemma_full = true;
    options_.use_gemma_refit = true;
    options_.no_loco = true;

    if (!loco_refit_grms.empty()) {
        bool first_refit = true;
        for (auto& entry : loco_refit_grms) {
            const std::string& chrom = entry.first;
            refit_grm_override_ = &entry.second;
            refit_chrom_filter_ = chrom;
            refit_append_mode_ = !first_refit;
            std::cout << "Hybrid Stage 2 LOCO GEMMA-full refinement chromosome " << chrom << "...\n";
            runGemmaRefit();
            first_refit = false;
        }
        refit_grm_override_ = nullptr;
        refit_chrom_filter_.clear();
        refit_append_mode_ = false;
    } else {
        runGemmaRefit();
    }

    options_.use_gemma_full = old_full;
    options_.use_gemma_refit = old_refit;
    options_.no_loco = old_no_loco;
    refine_mask_.clear();
}

void MultivariateModel::runGemmaRefit() {
    if (!options_.no_loco) {
        std::cout << "Multivariate GEMMA-style refit currently runs no-LOCO; use --no-loco for explicit validation.\n";
    }

    std::cout << "Multivariate GWAS method: "
              << (options_.use_gemma_full ? "GEMMA-full MVP per-SNP refit (Wald/LRT/Score)" : "GEMMA-style per-SNP variance-component refit (Wald)")
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

    Eigen::MatrixXd G_sub(n_keep_, n_keep_);
    if (refit_grm_override_) {
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
        std::ofstream ofs(out_file);
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tWALD_STAT\tP";
        if (options_.use_gemma_full) {
            ofs << "\tP_WALD\tP_LRT\tP_SCORE\tP_WALD_GC\tP_LRT_GC\tP_SCORE_GC";
        }
        for (int t = 0; t < Y_.cols(); ++t) ofs << "\tBETA_" << (t + 1) << "\tSE_" << (t + 1);
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

            if (options_.use_gemma_full && options_.refit_p >= 1.0) {
                // Fast path: null Vg, Ve for Wald/LRT/Score (matches GEMMA mvLMM)
                Eigen::VectorXd beta_all;
                Eigen::MatrixXd Hinv;
                double loglik_alt = std::numeric_limits<double>::quiet_NaN();
                if (!mv_fit_gls(UtY, UtX_aug, D, Vg_null, Ve_null,
                                beta_all, Hinv, &loglik_alt)) continue;
                beta = beta_all.segment(offset, n_traits);
                Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                Eigen::MatrixXd Vbeta_inv = Vbeta.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                wald_stat = (beta.transpose() * Vbeta_inv * beta)(0, 0);
                if (!std::isfinite(wald_stat) || wald_stat < 0.0) continue;
                p_wald = mv_chisq_prob(wald_stat, n_traits);
                p_score_val = options_.test_score ? p_wald : p_score_val;
                if (options_.test_lrt && std::isfinite(loglik_alt) && std::isfinite(loglik_null)) {
                    double lrt = std::max(0.0, 2.0 * (loglik_alt - loglik_null));
                    p_lrt = mv_chisq_prob(lrt, n_traits);
                }
            } else {
                // Refit path: Score pre-screen + optional EM refit
                Eigen::VectorXd beta_score_all;
                Eigen::MatrixXd cov_score_all;
                bool need_score = options_.use_gemma_full &&
                                  (options_.test_score || options_.refit_p < 1.0);
                bool have_score = need_score && mv_fit_gls(UtY, UtX_aug, D, Vg_null, Ve_null,
                                                           beta_score_all, cov_score_all, nullptr);
                if (have_score) {
                    Eigen::VectorXd beta_s = beta_score_all.segment(offset, n_traits);
                    Eigen::MatrixXd Vs = cov_score_all.block(offset, offset, n_traits, n_traits);
                    Eigen::MatrixXd Vs_inv = Vs.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                    double stat = (beta_s.transpose() * Vs_inv * beta_s)(0, 0);
                    if (std::isfinite(stat) && stat >= 0.0) p_score_val = mv_chisq_prob(stat, n_traits);
                }

                bool do_full_refit = !options_.use_gemma_full ||
                                     (std::isfinite(p_score_val) && p_score_val <= options_.refit_p);

                if (!do_full_refit) {
                    if (!have_score) continue;
                    beta = beta_score_all.segment(offset, n_traits);
                    Vbeta = cov_score_all.block(offset, offset, n_traits, n_traits);
                } else {
                    Eigen::MatrixXd Vg, Ve;
                    int snp_refit_maxit = std::max(30, options_.reml_maxit / 20);
                    double snp_refit_tol = std::max(options_.reml_tol * 10.0, 1e-7);
                    if (!mv_em_refit_rotated(UtY, UtX_aug, D, snp_refit_maxit, snp_refit_tol,
                                             Vg, Ve, &Vg_null, &Ve_null)) continue;
                    Eigen::VectorXd beta_all;
                    Eigen::MatrixXd Hinv;
                    double loglik_alt = std::numeric_limits<double>::quiet_NaN();
                    if (!mv_fit_gls(UtY, UtX_aug, D, Vg, Ve, beta_all, Hinv, &loglik_alt)) continue;
                    beta = beta_all.segment(offset, n_traits);
                    Vbeta = Hinv.block(offset, offset, n_traits, n_traits);
                    Eigen::MatrixXd Vbeta_inv = Vbeta.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
                    wald_stat = (beta.transpose() * Vbeta_inv * beta)(0, 0);
                    if (!std::isfinite(wald_stat) || wald_stat < 0.0) continue;
                    p_wald = mv_chisq_prob(wald_stat, n_traits);
                    p_score_val = options_.test_score ? p_score_val : p_score_val;
                    if (options_.use_gemma_full && options_.test_lrt && std::isfinite(loglik_alt) && std::isfinite(loglik_null)) {
                        double lrt = std::max(0.0, 2.0 * (loglik_alt - loglik_null));
                        p_lrt = mv_chisq_prob(lrt, n_traits);
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
            res.p_value = std::isfinite(p_wald) ? p_wald : p_score_val;
            res.p_wald = p_wald;
            res.p_lrt = p_lrt;
            res.p_score = p_score_val;
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

        std::ofstream ofs(out_file, std::ios::app);
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
            for (int t = 0; t < Y_.cols(); ++t) {
                ofs << "\t" << res.beta(t) << "\t" << std::sqrt(std::max(0.0, res.se(t, t)));
            }
            ofs << "\n";
        }
    } else {
        writeResults(out_file, results);
    }
    std::cout << "Results saved to " << out_file << std::endl;
}

void MultivariateModel::writeResults(const std::string& out_file, const std::vector<MvGwasResult>& results) {
    int n_traits = Y_.cols();
    std::ofstream ofs(out_file, std::ios::app);
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.wald_stat << "\t" << res.p_value;
        for (int t = 0; t < n_traits; ++t) {
            double beta = res.beta(t);
            double se = std::sqrt(res.se(t, t));
            ofs << "\t" << beta << "\t" << se;
        }
        if (options_.use_eb) {
            double eb_confidence = std::isfinite(res.eb_confidence) ? res.eb_confidence : 0.0;
            eb_confidence = std::max(0.0, std::min((double)n_traits, eb_confidence));
            eb_effective_tests_ += eb_confidence;
            eb_tested_markers_++;
            ofs << "\t" << res.eb_marker_lambda << "\t" << eb_confidence;
        }
        ofs << "\n";
    }
}

void MultivariateModel::writeEbSummary(const std::string& out_file) const {
    std::ofstream ofs(out_file);
    if (!ofs) throw std::runtime_error("Cannot write multivariate EB summary: " + out_file);
    ofs << std::setprecision(12);
    double m_eff = std::max(eb_effective_tests_, 1.0);
    ofs << "Metric\tValue\n";
    ofs << "method\tEB-mvLMM random-SNP scan (Phi=lambda*Ve)\n";
    ofs << "tested_markers\t" << eb_tested_markers_ << "\n";
    ofs << "traits\t" << Y_.cols() << "\n";
    ofs << "effective_tests\t" << eb_effective_tests_ << "\n";
    ofs << "effective_tests_per_marker\t" << (eb_tested_markers_ > 0 ? eb_effective_tests_ / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "bonferroni_0.05\t" << (eb_tested_markers_ > 0 ? 0.05 / (double)eb_tested_markers_ : 0.0) << "\n";
    ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    ofs << "note\tExperimental summary-level EB shrinkage; validate Type I error before production use.\n";
}

} // namespace cosmic
