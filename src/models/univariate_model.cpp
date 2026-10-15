#include "models/univariate_model.h"
#include "genetics/genotype_encoding.h"
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
#include "analysis_dataset.h"
#include "covar_reader.h"
#include "pheno_reader.h"
#include "validation.h"

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
    opts.out_gz = cfg.out_gz;
    opts.manhattan_tsv = cfg.manhattan_tsv;
    opts.qq_tsv = cfg.qq_tsv;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.bgen_dom_max_missing = cfg.bgen_dom_max_missing;
    opts.bgen_dom_min_info = cfg.bgen_dom_min_info;
    opts.bgen_dom_entropy_warning = cfg.bgen_dom_entropy_warning;
    opts.fast_lmm = cfg.fast_lmm;
    opts.fastlmm = cfg.fastlmm;
    opts.fastlmm_trace_threshold = cfg.fastlmm_trace_threshold;
    opts.fastlmm_min_rank = cfg.fastlmm_min_rank;
    opts.fastlmm_max_rank = cfg.fastlmm_max_rank;
    opts.fastlmm_matrix_free = cfg.fastlmm_matrix_free;
    opts.fastlmm_no_fallback = cfg.fastlmm_no_fallback;
    opts.lmm_low_rank = cfg.lmm_low_rank;
    opts.lmm_low_rank_oversample = cfg.lmm_low_rank_oversample;
    opts.lmm_low_rank_seed = cfg.lmm_low_rank_seed;
    opts.lmm_low_rank_error_tol = cfg.lmm_low_rank_error_tol;
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
    opts.snp_model = cfg.snp_model;
    opts.gwas_model = cfg.gwas_model;
    opts.rare_method = cfg.rare_method;
    opts.farmcpu = cfg.farmcpu;
    opts.blink = cfg.blink;
    opts.firth_logistic = cfg.firth_logistic;
    opts.saige_style = cfg.saige_style;
    opts.logistic_test = cfg.logistic_test;
    opts.logistic_glmm = cfg.logistic_glmm;
    opts.conditional_joint = cfg.conditional_joint;
    opts.farmcpu_max_iterations = cfg.farmcpu_max_iterations;
    opts.farmcpu_pseudo_qtn = cfg.farmcpu_pseudo_qtn;
    opts.condition_snp_file = cfg.condition_snp_file;
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
    opts.rare_max_maf = cfg.rare_max_maf;
    opts.rare_min_mac = cfg.rare_min_mac;
    opts.set_file = cfg.set_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.mpheno = cfg.mpheno;
    opts.sparse_cutoff = cfg.sparse_cutoff;
    return opts;
}

UnivariateModel::UnivariateModel(const Options& opts) : options_(opts) {}

bool UnivariateModel::usesAnalysisDatasetContract() const {
    return !options_.linear_model &&
           !options_.eb_lmm &&
           !options_.biobank_mode &&
           !options_.run_skat &&
           !options_.run_burden;
}

AnalysisDataset UnivariateModel::buildAnalysisDatasetForStandardLmm(PlinkReader& plink_reader,
                                                                    BgenReader& bgen_reader,
                                                                    PgenReader& pgen_reader,
                                                                    VcfReader& vcf_reader) const {
    AnalysisDataset dataset;
    const auto sample_keys = sample_keys_from_fid_iid(fids_, iids_);

    GenotypeSourceMetadata genotype_metadata;
    genotype_metadata.sample_order = sample_keys;
    genotype_metadata.n_samples = sample_keys.size();
    genotype_metadata.n_markers = generic_reader_ ? generic_reader_->getNumVariants()
                                : is_vcf_ ? vcf_reader.getNumVariants()
                                : is_pgen_ ? pgen_reader_->getNumVariants()
                                : is_bgen_ ? bgen_reader_->getNumVariants()
                                           : plink_reader.getNumSnps();
    genotype_metadata.format = generic_reader_ ? GenotypeStorageFormat::ExternalMatrix
                             : is_vcf_ ? GenotypeStorageFormat::Vcf
                             : is_pgen_ ? GenotypeStorageFormat::PlinkPgen
                             : is_bgen_ ? GenotypeStorageFormat::Bgen
                                        : GenotypeStorageFormat::PlinkBed;
    genotype_metadata.capabilities.supports_streaming = true;
    genotype_metadata.capabilities.supports_random_access = !is_bgen_;
    genotype_metadata.capabilities.supports_prefetch = is_bgen_;
    genotype_metadata.capabilities.materializes_full_matrix = false;
    dataset.setGenotypeSource(std::make_shared<MetadataOnlyGenotypeSource>(std::move(genotype_metadata)));

    Eigen::VectorXd y_aligned(n_samples_);
    if (options_.pheno_file.empty()) {
        if (is_bgen_ || is_pgen_ || generic_reader_) {
            throw std::runtime_error("BGEN/PGEN/HapMap/text input requires --pheno file.");
        }
        const auto& fam = plink_reader.getFamInfo();
        for (int i = 0; i < n_samples_; ++i) y_aligned(i) = fam[static_cast<std::size_t>(i)].phenotype;
    } else {
        PhenoReader pheno_reader;
        pheno_reader.load(options_.pheno_file, options_.mpheno);
        pheno_reader.align(fids_, iids_, y_aligned);
    }

    std::vector<double> y_values;
    std::vector<bool> y_missing;
    y_values.reserve(static_cast<std::size_t>(y_aligned.size()));
    y_missing.reserve(static_cast<std::size_t>(y_aligned.size()));
    for (int i = 0; i < y_aligned.size(); ++i) {
        y_values.push_back(y_aligned(i));
        y_missing.push_back(std::isnan(y_aligned(i)));
    }
    PhenotypeTable phenotype = phenotype_table_from_records(sample_keys, y_values, {"trait"}, y_missing);

    Eigen::MatrixXf covariates;
    std::vector<SampleKey> covariate_keys;
    std::vector<std::string> covariate_names;
    bool has_covariates = false;
    int total_covariate_cols = 0;
    if (!options_.qcovar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.load(options_.qcovar_file);
        Eigen::MatrixXf qcovariates;
        cov_reader.align(fids_, iids_, qcovariates);
        total_covariate_cols += qcovariates.cols();
        covariates.conservativeResize(n_samples_, total_covariate_cols);
        covariates.rightCols(qcovariates.cols()) = qcovariates;
        for (int c = 0; c < qcovariates.cols(); ++c) covariate_names.push_back("qcovar" + std::to_string(c + 1));
        has_covariates = true;
    }
    if (!options_.covar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.loadCategorical(options_.covar_file);
        Eigen::MatrixXf ccovariates;
        cov_reader.align(fids_, iids_, ccovariates);
        const int old_cols = total_covariate_cols;
        total_covariate_cols += ccovariates.cols();
        covariates.conservativeResize(n_samples_, total_covariate_cols);
        covariates.middleCols(old_cols, ccovariates.cols()) = ccovariates;
        for (int c = 0; c < ccovariates.cols(); ++c) covariate_names.push_back("covar" + std::to_string(c + 1));
        has_covariates = true;
    }

    if (has_covariates) {
        std::vector<double> cov_values;
        std::vector<bool> cov_missing;
        cov_values.reserve(static_cast<std::size_t>(covariates.rows() * covariates.cols()));
        cov_missing.reserve(static_cast<std::size_t>(covariates.rows() * covariates.cols()));
        for (int r = 0; r < covariates.rows(); ++r) {
            for (int c = 0; c < covariates.cols(); ++c) {
                cov_values.push_back(covariates(r, c));
                cov_missing.push_back(std::isnan(covariates(r, c)));
            }
        }
        covariate_keys = sample_keys;
        dataset.setCovariateTable(covariate_table_from_numeric(covariate_keys, covariate_names, cov_values, cov_missing));
    }

    dataset.setPhenotypeTable(std::move(phenotype));
    dataset.setSampleRegistry(SampleRegistry::fromGenotypeOrder(
        sample_keys, sample_keys, has_covariates ? covariate_keys : std::vector<SampleKey>{}));
    ValidationReport construction_report;
    validate_analysis_dataset_prerun(
        construction_report, dataset, "CosmicScanner standard LMM",
        true, true, has_covariates, true, true);
    construction_report.throw_if_fail();
    return dataset;
}

void UnivariateModel::applyAnalysisDatasetForStandardLmm(const AnalysisDataset& dataset) {
    const auto y_values = phenotype_values_in_analysis_order(dataset);
    y_.resize(static_cast<int>(y_values.size()));
    for (int i = 0; i < y_.size(); ++i) y_(i) = y_values[static_cast<std::size_t>(i)];

    std::size_t cov_cols = 0;
    const auto cov_values = covariate_values_in_analysis_order(dataset, &cov_cols);
    X_.resize(static_cast<int>(dataset.samples().size()), static_cast<int>(cov_cols + 1));
    X_.col(0).setOnes();
    for (std::size_t r = 0; r < dataset.samples().size(); ++r) {
        for (std::size_t c = 0; c < cov_cols; ++c) {
            X_(static_cast<int>(r), static_cast<int>(c + 1)) =
                static_cast<float>(cov_values[r * cov_cols + c]);
        }
    }
    if (cov_cols > 0) {
        std::cout << "Using " << (cov_cols + 1) << " fixed-effect columns (including intercept).\n";
    }
}

void UnivariateModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              VcfReader& vcf_reader,
                               bool is_bgen, bool is_pgen, bool is_vcf,
                               SequentialGenotypeMatrix* generic_reader) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    vcf_reader_ = &vcf_reader;
    generic_reader_ = generic_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;
    is_vcf_ = is_vcf;

    if (usesAnalysisDatasetContract()) {
        auto dataset = buildAnalysisDatasetForStandardLmm(plink_reader, *bgen_reader_, *pgen_reader_, *vcf_reader_);
        ValidationReport report;
        validate_analysis_dataset_prerun(
            report, dataset, "CosmicScanner standard LMM",
            true, true,
            !options_.qcovar_file.empty() || !options_.covar_file.empty(),
            true, true);
        std::cout << report.summary(false);
        report.throw_if_fail();
        applyAnalysisDatasetForStandardLmm(dataset);
        std::cout << dataset.diagnosticSummary("production-contract mode");
        analysis_dataset_ = dataset;
        has_analysis_dataset_ = true;
        validation_report_.merge(report);
    } else {
        y_ = DataLoader::loadPhenotype(options_.pheno_file, options_.mpheno, fids_, iids_, n_samples_, plink_reader,
                                       is_bgen_ || is_vcf_ || generic_reader_, is_pgen_);
        X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
    }

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering missing data, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples.");
}

void UnivariateModel::run() {
    if (!options_.logistic_test.empty() || options_.logistic_glmm) {
        runLogistic();
        return;
    }
    if (options_.saige_style) {
        runSaigeStyle();
        return;
    }
    if (options_.farmcpu || options_.blink || options_.firth_logistic ||
        !options_.condition_snp_file.empty()) {
        if (generic_reader_) throw std::runtime_error("Multi-locus models are not supported for HapMap/text input");
        runMultiLocusModel();
        return;
    }
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
    if ((options_.lmm_low_rank > 0 || options_.fastlmm) &&
        (options_.eb_lmm || options_.biobank_mode || options_.profile_scale ||
          options_.lmm_gemma_refit || options_.use_spa || is_bgen_ || is_pgen_ || is_vcf_ || generic_reader_)) {
        throw std::runtime_error("Low-rank LMM currently requires BED input and the standard fixed-null Gaussian LMM");
    }

    GwasScanner::Options opts;
    opts.out_prefix = options_.out_prefix;
    opts.out_ssf = options_.out_ssf;
    opts.min_maf = options_.min_maf;
    opts.min_hwe = options_.min_hwe;
    opts.bgen_dom_max_missing = options_.bgen_dom_max_missing;
    opts.bgen_dom_min_info = options_.bgen_dom_min_info;
    opts.bgen_dom_entropy_warning = options_.bgen_dom_entropy_warning;
    opts.snp_model = options_.snp_model;
    opts.gwas_model_md = options_.gwas_model;
    opts.fast_lmm = options_.fast_lmm;
    opts.fastlmm = options_.fastlmm;
    opts.fastlmm_trace_threshold = options_.fastlmm_trace_threshold;
    opts.fastlmm_min_rank = options_.fastlmm_min_rank;
    opts.fastlmm_max_rank = options_.fastlmm_max_rank;
    opts.lmm_low_rank = options_.lmm_low_rank;
    opts.lmm_low_rank_oversample = options_.lmm_low_rank_oversample;
    opts.lmm_low_rank_seed = options_.lmm_low_rank_seed;
    opts.lmm_low_rank_error_tol = options_.lmm_low_rank_error_tol;
    opts.fastlmm_matrix_free = options_.fastlmm_matrix_free;
    opts.fastlmm_no_fallback = options_.fastlmm_no_fallback;
    opts.out_gz = options_.out_gz;
    opts.manhattan_tsv = options_.manhattan_tsv;
    opts.qq_tsv = options_.qq_tsv;
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
    opts.rare_method = options_.rare_method;
    opts.run_burden = options_.run_burden;
    opts.rare_max_maf = options_.rare_max_maf;
    opts.rare_min_mac = options_.rare_min_mac;
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

    if (options_.run_skat || options_.run_burden) {
        if (generic_reader_) throw std::runtime_error("Rare-variant tests are not supported for HapMap/text input");
        scanner.runRareVariantScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, *vcf_reader_,
                                   is_bgen_, is_pgen_, is_vcf_, keep_);
    } else if (options_.biobank_mode) {
        if (generic_reader_) throw std::runtime_error("REGENIE-style biobank mode is not supported for HapMap/text input");
        scanner.runBiobankScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else if (options_.gwas_model == "additive_dominance") {
        if (is_bgen_) {
            throw std::runtime_error(
                "--gwas-model additive_dominance is not supported for BGEN input: the BGEN "
                "dosage-only backend does not expose P(AB) or hardcall genotype classes, so it "
                "cannot represent the production dominance model. Convert to PLINK BED or PGEN "
                "hardcalls, or use the default additive model.");
        }
        if (is_vcf_ || generic_reader_) {
            throw std::runtime_error("--gwas-model additive_dominance currently requires PLINK BED or PGEN hardcall input");
        }
        scanner.runAdLocoScan(y_, X_, *plink_reader_, *pgen_reader_, is_pgen_, keep_);
    } else if (options_.gwas_model == "dominance") {
        if (is_vcf_ || generic_reader_) {
            throw std::runtime_error("--gwas-model dominance requires PLINK BED, PGEN hardcall, or BGEN probability input");
        }
        scanner.runDominanceLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_,
                                     is_bgen_, is_pgen_, keep_);
    } else {
        scanner.runLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, *vcf_reader_,
                            is_bgen_, is_pgen_, is_vcf_, keep_, generic_reader_);
    }
    validation_report_.merge(scanner.validationReport());
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
    const bool exact_hwe_available = !is_bgen_;
    {
        std::ofstream meta(options_.out_prefix + ".scan.meta.txt");
        meta << "Item\tValue\n"
             << "analysis_model\tunivariate_lm\n"
             << "SNP_MODEL\t" << options_.snp_model << "\n"
             << "HWE_APPLIED\t" << (options_.min_hwe > 0.0 && exact_hwe_available ? "true" : "false") << "\n"
             << "HWE_METHOD\t" << (options_.min_hwe > 0.0 && exact_hwe_available ? "EXACT_HARDCALL" : "NONE") << "\n"
             << "HWE_REASON_IF_SKIPPED\t" << (options_.min_hwe <= 0.0 ? "FILTER_DISABLED" : exact_hwe_available ? "NONE" : "CONTINUOUS_DOSAGE_WITHOUT_PROBABILITIES") << "\n";
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

    const size_t m = generic_reader_ ? generic_reader_->getNumVariants()
                   : is_vcf_ ? vcf_reader_->getNumVariants()
                   : is_pgen_ ? pgen_reader_->getNumVariants()
                   : is_bgen_ ? bgen_reader_->getNumVariants()
                              : plink_reader_->getNumSnps();
    const auto& bim = plink_reader_->getBimInfo();
    if (is_pgen_) pgen_reader_->reset();
    if (is_bgen_) bgen_reader_->reset();
    if (is_vcf_) vcf_reader_->reset();
    if (generic_reader_) generic_reader_->reset();

    size_t tested = 0;
    const size_t block_size = 512;
    for (size_t start = 0; start < m; start += block_size) {
        const size_t end = std::min(m, start + block_size);
        size_t current_block = end - start;
        Eigen::MatrixXd raw;
        std::vector<PlinkBimInfo> block_bim(current_block);

        if (generic_reader_) {
            raw.resize(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, id, ref, alt;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!generic_reader_->readVariant(chrom, id, pos, ref, alt, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = dosage_f[i];
                block_bim[k] = {chrom, id, 0.0, static_cast<long long>(pos), alt, ref};
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else if (is_vcf_) {
            raw.resize(n_samples_, current_block);
            Eigen::VectorXf dosage_f;
            std::string chrom, id, ref, alt;
            uint32_t pos = 0;
            size_t k = 0;
            for (; k < current_block; ++k) {
                if (!vcf_reader_->readVariant(chrom, id, pos, ref, alt, dosage_f)) break;
                for (int i = 0; i < n_samples_; ++i) raw(i, k) = dosage_f[i];
                block_bim[k] = {chrom, id, 0.0, static_cast<long long>(pos), alt, ref};
            }
            if (k == 0) break;
            current_block = k;
            raw.conservativeResize(Eigen::NoChange, current_block);
            block_bim.resize(current_block);
        } else if (is_pgen_) {
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

            Eigen::VectorXd marker_raw(n);
            for (int i = 0; i < n; ++i) {
                marker_raw[i] = raw(keep_.keep_indices[i], k);
            }
            const auto prepared = genetics::prepareAssociationMarker(
                marker_raw, options_.snp_model, options_.min_maf, options_.min_hwe,
                exact_hwe_available, false);
            if (!prepared.passes_maf || !prepared.passes_hwe) continue;
            const Eigen::VectorXd& g = prepared.values;
            const double freq = prepared.allele_frequency;

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
