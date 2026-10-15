#include "models/survival_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"

namespace cosmic {

SurvivalModel::Options SurvivalModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.grm_alg = cfg.grm_alg;
    opts.n_threads = 1;
#ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
#endif
    opts.no_loco = cfg.no_loco;
    opts.use_spa = cfg.use_spa;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    return opts;
}

SurvivalModel::SurvivalModel(const Options& opts) : options_(opts) {}

void SurvivalModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              VcfReader& vcf_reader,
                              bool is_bgen, bool is_pgen, bool is_vcf) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    vcf_reader_ = &vcf_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;
    is_vcf_ = is_vcf;

    loadPhenotype();
    loadCovariates();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for survival model.");
}

void SurvivalModel::loadPhenotype() {
    std::string file = options_.pheno_file;
    if (file.empty()) throw std::runtime_error("Survival model requires --pheno with FID IID TIME EVENT columns.");

    std::ifstream ifs(file);
    if (!ifs.is_open()) throw std::runtime_error("Cannot open: " + file);

    std::string line;
    std::vector<double> y_vals(n_samples_, 0.0), t_vals, e_vals;

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<std::string> cols;
        std::string token;
        while (iss >> token) cols.push_back(token);
        if (cols.size() < (size_t)options_.event_col) continue;

        std::string fid = cols[0], iid = cols[1];
        int idx = -1;
        for (int j = 0; j < n_samples_; ++j) {
            if (fids_[j] == fid && iids_[j] == iid) { idx = j; break; }
        }
        if (idx < 0) continue;

        t_vals.push_back(std::stod(cols[options_.time_col - 1]));
        e_vals.push_back(std::stod(cols[options_.event_col - 1]));
        y_vals[idx] = 1.0;
    }
    ifs.close();

    y_.resize(n_samples_);
    for (int i = 0; i < n_samples_; ++i) y_(i) = y_vals[i];
    survival_time_ = Eigen::Map<Eigen::VectorXd>(t_vals.data(), t_vals.size());
    event_flag_ = Eigen::Map<Eigen::VectorXd>(e_vals.data(), e_vals.size());

    std::cout << "Loaded " << e_vals.size() << " survival records ("
              << (int)event_flag_.sum() << " events).\n";
}

void SurvivalModel::loadCovariates() {
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
}

void SurvivalModel::runCoxNull() {
    std::cout << "Running Cox null model...\n";
}

void SurvivalModel::martingaleResiduals(Eigen::VectorXd& M, Eigen::VectorXd& eta) {
    int n = keep_.n_keep;
    M.resize(n);
    eta = Eigen::VectorXd::Ones(n);

    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;

    Eigen::VectorXd cum_haz = Eigen::VectorXd::Zero(n);
    double cum_sum = 0.0;
    for (int i = n - 1; i >= 0; --i) {
        int idx = order[i];
        cum_sum += event_flag_(idx) / (n - i);
        cum_haz(idx) = cum_sum;
    }

    for (int i = 0; i < n; ++i) {
        M(i) = event_flag_(i) - cum_haz(i);
    }
}

void SurvivalModel::scanGwas() {
    int n = keep_.n_keep;

    Eigen::VectorXd M, eta;
    martingaleResiduals(M, eta);

    std::cout << "Running GWAS on martingale residuals...\n";

    GwasScanner::Options gopts;
    gopts.out_prefix = options_.out_prefix;
    gopts.min_maf = options_.min_maf;
    gopts.min_hwe = options_.min_hwe;
    gopts.fast_lmm = true;
    gopts.n_threads = options_.n_threads;
    gopts.no_loco = options_.no_loco;
    gopts.use_spa = options_.use_spa;

    GwasScanner scanner(gopts);
    scanner.runLocoScan(M, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, *vcf_reader_,
                        is_bgen_, is_pgen_, is_vcf_, keep_);
}

void SurvivalModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_survival_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) return;
    ofs << "sigma2_g\tsigma2_e\n";
    ofs << sigma2_g_ << "\t" << sigma2_e_ << "\n";
    ofs.close();
    std::cout << "Survival null model written to " << fname << "\n";
}

void SurvivalModel::run() {
    runCoxNull();
    scanGwas();
    writeNullSummary();
}

} // namespace cosmic
