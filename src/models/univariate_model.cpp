#include "models/univariate_model.h"
#include <iostream>
#include <stdexcept>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "gwas_scanner.h"

namespace cosmic {

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
    GwasScanner::Options opts;
    opts.out_prefix = options_.out_prefix;
    opts.out_ssf = options_.out_ssf;
    opts.min_maf = options_.min_maf;
    opts.min_hwe = options_.min_hwe;
    opts.fast_lmm = options_.fast_lmm;
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

    if (options_.run_skat || options_.run_burden) {
        scanner.runRareVariantScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else if (options_.biobank_mode) {
        scanner.runBiobankScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    } else {
        scanner.runLocoScan(y_, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_, keep_);
    }
}

} // namespace cosmic
