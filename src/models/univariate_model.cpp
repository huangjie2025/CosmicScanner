#include "models/univariate_model.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "gwas_scanner.h"
#include "qc_utils.h"

namespace cosmic {

namespace {

std::unordered_set<std::string> load_lm_snp_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (ifs >> id) ids.insert(id);
    if (ids.empty()) throw std::runtime_error("SNP list is empty: " + path);
    return ids;
}

bool lm_region_passes(const PlinkBimInfo& snp,
                      const std::string& chr,
                      long long from_bp,
                      long long to_bp) {
    if (!chr.empty() && snp.chrom != chr) return false;
    if (from_bp >= 0 && snp.bp_pos < from_bp) return false;
    if (to_bp >= 0 && snp.bp_pos > to_bp) return false;
    return true;
}

double lm_gamma_q(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x == 0.0) return 1.0;
    const double eps = 1e-12;
    const double fpmin = 1e-300;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < 200; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * eps) {
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / fpmin;
        double d = 1.0 / std::max(b, fpmin);
        double h = d;
        for (int i = 1; i <= 200; ++i) {
            double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < fpmin) d = fpmin;
            c = b + an / c;
            if (std::abs(c) < fpmin) c = fpmin;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < eps) {
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
            }
        }
    }
    return 0.0;
}

double lm_betacf(double a, double b, double x) {
    const int maxit = 200;
    const double eps = 3e-14;
    const double fpmin = 1e-300;
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < fpmin) d = fpmin;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= maxit; ++m) {
        int m2 = 2 * m;
        double aa = static_cast<double>(m) * (b - m) * x /
                    ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x /
             ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < fpmin) d = fpmin;
        c = 1.0 + aa / c;
        if (std::abs(c) < fpmin) c = fpmin;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) <= eps) break;
    }
    return h;
}

double lm_reg_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                         a * std::log(x) + b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * lm_betacf(a, b, x) / a;
    }
    return 1.0 - bt * lm_betacf(b, a, 1.0 - x) / b;
}

double lm_fdist_q(double f, double d1, double d2) {
    if (!std::isfinite(f) || f <= 0.0) return 1.0;
    double x = d2 / (d2 + d1 * f);
    return std::clamp(lm_reg_incomplete_beta(d2 / 2.0, d1 / 2.0, x), 0.0, 1.0);
}

double lm_chisq_q(double stat, int df) {
    if (!std::isfinite(stat) || stat <= 0.0) return 1.0;
    return std::clamp(lm_gamma_q(df / 2.0, stat / 2.0), 0.0, 1.0);
}

} // namespace

static EigenVce::NullSolver parse_lmm_null_solver(const std::string& value) {
    if (value == "brent") return EigenVce::NullSolver::Brent;
    if (value == "idul") return EigenVce::NullSolver::IDUL;
    if (value == "idul_plus" || value == "idul+") return EigenVce::NullSolver::IDULPlus;
    if (value == "aidul") return EigenVce::NullSolver::AIDUL;
    throw std::runtime_error("Unknown --lmm-null-solver: " + value + " (expected brent, idul, idul_plus, or aidul)");
}

UnivariateModel::Options UnivariateModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.out_ssf = cfg.out_ssf;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.fast_lmm = cfg.fast_lmm;
    opts.scan_tier = cfg.scan_tier;
    if (!opts.scan_tier.empty() && cfg.fast_lmm && opts.scan_tier != "approx") {
        throw std::runtime_error("--fast-lmm is a legacy fixed-null alias and conflicts with --scan-tier " + opts.scan_tier);
    }
    if (!opts.scan_tier.empty() && cfg.lmm_gemma_refit && opts.scan_tier != "exact") {
        throw std::runtime_error("--lmm-gemma-refit conflicts with --scan-tier " + opts.scan_tier);
    }
    if (opts.scan_tier == "semi") {
        opts.profile_scale = true;
        opts.lmm_gemma_refit = false;
    } else if (opts.scan_tier == "exact") {
        opts.lmm_gemma_refit = true;
    } else if (opts.scan_tier == "approx") {
        opts.lmm_gemma_refit = false;
    }
    opts.lmm_gemma_mode = cfg.lmm_gemma_mode;
    if (opts.scan_tier.empty()) opts.lmm_gemma_refit = cfg.lmm_gemma_refit;
    opts.linear_model = cfg.linear_model;
    opts.linear_model_mode = cfg.linear_model_mode;
    opts.eb_lmm = cfg.eb_lmm;
    opts.n_threads = 1;
    #ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
    #endif
    opts.no_loco = cfg.no_loco;
    opts.grm_alg = cfg.grm_alg;
    opts.lmm_null_solver = parse_lmm_null_solver(cfg.lmm_null_solver);
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.biobank_mode = cfg.biobank_mode;
    opts.use_covariate_projection = cfg.use_covariate_projection;
    opts.use_spa = cfg.use_spa;
    opts.regenie_block_size = cfg.block_size;
    opts.regenie_model_snps = cfg.regenie_model_snps;
    opts.regenie_write_pred = cfg.regenie_write_pred;
    opts.regenie_use_pred = cfg.regenie_use_pred;
    opts.regenie_step = cfg.regenie_step;
    opts.scan_chr = cfg.scan_chr;
    opts.scan_from_bp = cfg.scan_from_bp;
    opts.scan_to_bp = cfg.scan_to_bp;
    opts.run_skat = cfg.run_skat;
    opts.run_burden = cfg.run_burden;
    opts.set_file = cfg.set_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.mpheno = cfg.mpheno;
    return opts;
}

UnivariateModel::UnivariateModel(const Options& opts) : options_(opts) {}

void UnivariateModel::prepare(const std::vector<std::string>& fids,
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

    y_ = DataLoader::loadPhenotype(options_.pheno_file, options_.mpheno, fids_, iids_, n_samples_, plink_reader, is_bgen_, is_pgen_);
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering missing data, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples.");
}

void UnivariateModel::run() {
    if (options_.linear_model) {
        runLinearModel();
        return;
    }
    if (options_.lmm_gemma_mode < 0 || options_.lmm_gemma_mode > 4) {
        throw std::runtime_error("--lmm-gemma-mode must be 1, 2, 3, or 4");
    }
    if (options_.lmm_gemma_mode > 0 && options_.biobank_mode) {
        throw std::runtime_error("--lmm-gemma-mode is for exact/fast LMM scans; it is not supported with --biobank-mode");
    }
    if (options_.lmm_gemma_refit && options_.eb_lmm) {
        throw std::runtime_error("--lmm-gemma-refit cannot be combined with EB-LMM");
    }
    if (!options_.scan_tier.empty() && (options_.eb_lmm || options_.biobank_mode)) {
        throw std::runtime_error("--scan-tier is currently supported for the standard Gaussian single-trait LMM only");
    }

    GwasScanner::Options opts;
    opts.out_prefix = options_.out_prefix;
    opts.out_ssf = options_.out_ssf;
    opts.min_maf = options_.min_maf;
    opts.min_hwe = options_.min_hwe;
    opts.fast_lmm = options_.fast_lmm;
    opts.scan_tier = options_.scan_tier;
    opts.profile_scale = options_.profile_scale;
    opts.lmm_gemma_mode = options_.lmm_gemma_mode;
    opts.lmm_gemma_refit = options_.lmm_gemma_refit;
    opts.eb_lmm = options_.eb_lmm;
    #ifdef _OPENMP
    opts.n_threads = options_.n_threads;
    #endif
    opts.no_loco = options_.no_loco;
    opts.grm_alg = options_.grm_alg;
    opts.lmm_null_solver = options_.lmm_null_solver;
    opts.reml_maxit = options_.reml_maxit;
    opts.reml_tol = options_.reml_tol;
    opts.biobank_mode = options_.biobank_mode;
    opts.use_covariate_projection = options_.use_covariate_projection;
    opts.use_spa = options_.use_spa;
    opts.run_skat = options_.run_skat;
    opts.run_burden = options_.run_burden;
    opts.set_file = options_.set_file;
    opts.regenie_block_size = options_.regenie_block_size;
    opts.regenie_model_snps = options_.regenie_model_snps;
    opts.regenie_write_pred = options_.regenie_write_pred;
    opts.regenie_use_pred = options_.regenie_use_pred;
    opts.regenie_step = options_.regenie_step;
    opts.scan_chr = options_.scan_chr;
    opts.scan_from_bp = options_.scan_from_bp;
    opts.scan_to_bp = options_.scan_to_bp;
    opts.extract_file = options_.extract_file;
    opts.exclude_file = options_.exclude_file;

    GwasScanner scanner(opts);

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
    if (options_.run_skat || options_.run_burden) {
        scanner.runRareVariantScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else if (options_.biobank_mode) {
        scanner.runBiobankScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else
#endif
    {
        scanner.runLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    }
}

void UnivariateModel::runLinearModel() {
    if (options_.linear_model_mode < 1 || options_.linear_model_mode > 4) {
        throw std::runtime_error("--lm-mode must be 1, 2, 3, or 4");
    }

    const auto extract_set = load_lm_snp_id_set(options_.extract_file);
    const auto exclude_set = load_lm_snp_id_set(options_.exclude_file);
    const int n = keep_.n_keep;
    const int p = X_.cols();
    const int df = n - p - 1;
    if (df <= 0) {
        throw std::runtime_error("Too few samples for GEMMA-style LM after covariates and SNP effect.");
    }

    Eigen::VectorXd y(n);
    Eigen::MatrixXd W(n, p);
    for (int i = 0; i < n; ++i) {
        const int idx = keep_.keep_indices[i];
        y[i] = y_[idx];
        for (int j = 0; j < p; ++j) W(i, j) = static_cast<double>(X_(idx, j));
    }

    Eigen::MatrixXd WtW_inv;
    Eigen::VectorXd Wty;
    double yPwy = y.squaredNorm();
    if (p > 0) {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(W.transpose() * W);
        if (ldlt.info() != Eigen::Success) throw std::runtime_error("LM covariate matrix factorization failed.");
        WtW_inv = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
        Wty = W.transpose() * y;
        yPwy -= Wty.dot(WtW_inv * Wty);
    }
    if (!std::isfinite(yPwy) || yPwy <= 0.0) {
        throw std::runtime_error("LM residual phenotype variance is non-positive after covariate projection.");
    }

    std::ofstream ofs(options_.out_prefix + ".lm.assoc.txt");
    if (!ofs) throw std::runtime_error("Cannot write LM association output: " + options_.out_prefix + ".lm.assoc.txt");
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tBETA\tSE\tP";
    if (options_.linear_model_mode == 4) ofs << "\tP_WALD\tP_LRT\tP_SCORE";
    ofs << "\n";
    ofs << std::setprecision(12);

    const size_t m = is_pgen_ ? pgen_reader_->getNumVariants()
                   : is_bgen_ ? bgen_reader_->getNumVariants()
                              : plink_reader_->getNumSnps();
    const auto& bim = plink_reader_->getBimInfo();
    if (is_pgen_) pgen_reader_->reset();
    if (is_bgen_) bgen_reader_->reset();

    size_t tested = 0;
    bool hwe_dosage_warning = false;
    const size_t block_size = 512;
    for (size_t start = 0; start < m; start += block_size) {
        const size_t end = std::min(m, start + block_size);
        size_t current_block = end - start;
        Eigen::MatrixXd raw;
        std::vector<PlinkBimInfo> block_bim(current_block);

        if (is_pgen_) {
            raw.resize(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, rsid, ref, alt;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = dosage_f[i];
                block_bim[k] = {chrom, rsid, 0.0, static_cast<long long>(pos), alt, ref};
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else if (is_bgen_) {
            raw.resize(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, rsid, a1, a2;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!bgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = dosage_f[i];
                block_bim[k] = {chrom, rsid, 0.0, static_cast<long long>(pos), a1, a2};
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else {
            plink_reader_->getSnpBlock(start, end, raw);
            for (size_t k = 0; k < current_block; ++k) block_bim[k] = bim[start + k];
        }

        for (size_t k = 0; k < current_block; ++k) {
            const auto& snp = block_bim[k];
            if (!lm_region_passes(snp, options_.scan_chr, options_.scan_from_bp, options_.scan_to_bp)) continue;
            if (!extract_set.empty() && !extract_set.count(snp.id)) continue;
            if (!exclude_set.empty() && exclude_set.count(snp.id)) continue;

            Eigen::VectorXd g(n);
            double sum = 0.0;
            int count = 0;
            long long hom0 = 0, het = 0, hom2 = 0;
            bool discrete_genotypes = true;
            for (int i = 0; i < n; ++i) {
                double val = raw(keep_.keep_indices[i], k);
                if (!std::isfinite(val)) continue;
                sum += val;
                ++count;
                double rounded = std::round(val);
                if (std::abs(val - rounded) <= 1e-6 && rounded >= 0.0 && rounded <= 2.0) {
                    if (rounded < 0.5) ++hom0;
                    else if (rounded < 1.5) ++het;
                    else ++hom2;
                } else {
                    discrete_genotypes = false;
                }
            }
            if (count == 0) continue;
            const double mean = sum / count;
            const double freq = mean / 2.0;
            const double maf = std::min(freq, 1.0 - freq);
            if (maf < options_.min_maf) continue;
            if (options_.min_hwe > 0.0) {
                if (discrete_genotypes && count == n) {
                    if (calculateHWE(hom0, het, hom2) < options_.min_hwe) continue;
                } else if (!hwe_dosage_warning) {
                    std::cerr << "Warning: --lm applies --hwe only to fully observed hard-call genotypes; "
                              << "dosage/missing variants are not HWE-filtered.\n";
                    hwe_dosage_warning = true;
                }
            }
            for (int i = 0; i < n; ++i) {
                double val = raw(keep_.keep_indices[i], k);
                g[i] = std::isfinite(val) ? val : mean;
            }

            double xPwx = g.squaredNorm();
            double xPwy = g.dot(y);
            if (p > 0) {
                Eigen::VectorXd Wtg = W.transpose() * g;
                xPwx -= Wtg.dot(WtW_inv * Wtg);
                xPwy -= Wtg.dot(WtW_inv * Wty);
            }
            if (!std::isfinite(xPwx) || xPwx <= 1e-12) continue;
            double yPxy = yPwy - xPwy * xPwy / xPwx;
            if (!std::isfinite(yPxy) || yPxy <= 0.0) yPxy = std::numeric_limits<double>::min();

            const double beta = xPwy / xPwx;
            const double se_wald = std::sqrt(yPxy / (static_cast<double>(df) * xPwx));
            const double se_score = std::sqrt(yPwy / (static_cast<double>(n) * xPwx));
            const double f_wald = (se_wald > 0.0) ? (beta * beta / (se_wald * se_wald)) : 0.0;
            const double f_score = (se_score > 0.0) ? (beta * beta / (se_score * se_score)) : 0.0;
            const double p_wald = lm_fdist_q(f_wald, 1.0, static_cast<double>(df));
            const double lrt = static_cast<double>(n) * (std::log(yPwy) - std::log(yPxy));
            const double p_lrt = lm_chisq_q(lrt, 1);
            const double p_score = lm_fdist_q(f_score, 1.0, static_cast<double>(df));

            double se = se_wald;
            double pval = p_wald;
            if (options_.linear_model_mode == 2) pval = p_lrt;
            else if (options_.linear_model_mode == 3) {
                pval = p_score;
                se = se_score;
            }

            ofs << snp.chrom << '\t' << snp.id << '\t' << snp.bp_pos << '\t'
                << snp.alt << '\t' << snp.ref << '\t' << freq << '\t'
                << beta << '\t' << se << '\t' << pval;
            if (options_.linear_model_mode == 4) {
                ofs << '\t' << p_wald << '\t' << p_lrt << '\t' << p_score;
            }
            ofs << '\n';
            ++tested;
        }
    }

    std::cout << "GEMMA-style LM scan complete. Tested " << tested
              << " variants. Results written to " << options_.out_prefix
              << ".lm.assoc.txt\n";
}

} // namespace cosmic
