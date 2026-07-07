#include "models/gemma_vc_model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <omp.h>

#include "data_loader.h"
#include "eigen_vce.h"
#include "grm_builder.h"

namespace cosmic {
namespace {

std::unordered_set<std::string> load_snp_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (in >> id) ids.insert(id);
    return ids;
}

Eigen::MatrixXd residualize(const Eigen::VectorXd& y, const Eigen::MatrixXf& Xf) {
    Eigen::MatrixXd X = Xf.cast<double>();
    if (X.cols() == 0) return y;
    Eigen::MatrixXd XtX = X.transpose() * X;
    XtX.diagonal().array() += 1e-10;
    Eigen::VectorXd beta = XtX.ldlt().solve(X.transpose() * y);
    return y - X * beta;
}

EigenVce::NullSolver parse_solver(const std::string& value, int vc_mode) {
    if (vc_mode == 3) return EigenVce::NullSolver::AIDUL;
    if (value == "brent") return EigenVce::NullSolver::Brent;
    if (value == "idul") return EigenVce::NullSolver::IDUL;
    if (value == "idul_plus" || value == "idul+") return EigenVce::NullSolver::IDULPlus;
    if (value == "aidul") return EigenVce::NullSolver::AIDUL;
    throw std::runtime_error("Unknown --lmm-null-solver for --vc: " + value);
}

} // namespace

GemmaVcModel::Options GemmaVcModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.mpheno = cfg.mpheno;
    opts.mode = cfg.vc_mode;
    opts.grm_alg = cfg.grm_alg;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.use_covariate_projection = cfg.use_covariate_projection;
    opts.lmm_null_solver = cfg.lmm_null_solver;
    #ifdef _OPENMP
    opts.n_threads = cfg.num_threads > 0 ? cfg.num_threads : omp_get_max_threads();
    #else
    opts.n_threads = 1;
    #endif
    return opts;
}

GemmaVcModel::GemmaVcModel(const Options& opts) : options_(opts) {}

void GemmaVcModel::prepare(const std::vector<std::string>& fids,
                           const std::vector<std::string>& iids,
                           int n_samples,
                           PlinkReader& plink_reader,
                           BgenReader& bgen_reader,
                           PgenReader& pgen_reader,
                           bool is_bgen,
                           bool is_pgen) {
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
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for --vc");
}

Eigen::MatrixXd GemmaVcModel::buildGrm() const {
    GrmBuilder::Options opts;
    opts.min_maf = static_cast<float>(options_.min_maf);
    opts.min_hwe = static_cast<float>(options_.min_hwe);
    opts.algorithm = options_.grm_alg;
    opts.extract_set = load_snp_set(options_.extract_file);
    opts.exclude_set = load_snp_set(options_.exclude_file);
    opts.thread_num = options_.n_threads;

    GrmBuilder builder(opts);
    if (is_bgen_) {
        bgen_reader_->reset();
        builder.compute(*bgen_reader_);
    } else if (is_pgen_) {
        pgen_reader_->reset();
        builder.compute(*pgen_reader_);
    } else {
        builder.compute(*plink_reader_);
    }
    return DataFilter::subsetSquare(builder.getGrm(), keep_);
}

GemmaVcModel::VcResult GemmaVcModel::runHe(const Eigen::VectorXd& y,
                                           const Eigen::MatrixXf& X,
                                           const Eigen::MatrixXd& G) const {
    Eigen::VectorXd r = residualize(y, X);
    const int n = static_cast<int>(r.size());
    const double phen_var = std::max(0.0, r.squaredNorm() / std::max(1, n - static_cast<int>(X.cols())));
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double kij = G(i, j);
            num += kij * r(i) * r(j);
            den += kij * kij;
        }
    }
    double vg = (den > 1e-20) ? num / den : 0.0;
    if (!std::isfinite(vg)) vg = 0.0;
    vg = std::max(0.0, vg);
    const double mean_diag = G.diagonal().mean();
    double ve = phen_var - vg * mean_diag;
    if (!std::isfinite(ve) || ve < 0.0) ve = std::max(0.0, phen_var - vg);

    VcResult res;
    res.method = "HE";
    res.sigma2_g = vg;
    res.sigma2_e = ve;
    res.pve = (vg + ve > 0.0) ? vg / (vg + ve) : 0.0;
    res.delta = (vg > 0.0) ? ve / vg : std::numeric_limits<double>::infinity();
    res.n = n;
    res.n_covars = X.cols();
    return res;
}

GemmaVcModel::VcResult GemmaVcModel::runReml(const Eigen::VectorXd& y,
                                             const Eigen::MatrixXf& X,
                                             const Eigen::MatrixXd& G) const {
    EigenVce::Options vopts;
    vopts.use_covariate_projection = options_.use_covariate_projection;
    vopts.null_solver = parse_solver(options_.lmm_null_solver, options_.mode);
    vopts.max_iter = options_.reml_maxit;
    vopts.tol = options_.reml_tol;
    EigenVce vce(vopts);
    vce.prepare(y, X, G);
    vce.runNullModel();

    VcResult res;
    res.method = options_.mode == 3 ? "AI_REML" : "REML";
    res.sigma2_g = vce.getSigma2g();
    res.sigma2_e = vce.getSigma2e();
    res.pve = (res.sigma2_g + res.sigma2_e > 0.0)
        ? res.sigma2_g / (res.sigma2_g + res.sigma2_e) : 0.0;
    res.delta = vce.getDelta();
    res.n = static_cast<int>(y.size());
    res.n_covars = X.cols();
    return res;
}

void GemmaVcModel::writeResult(const VcResult& res) const {
    std::ofstream vc(options_.out_prefix + ".vc.txt");
    if (!vc) throw std::runtime_error("Cannot write --vc output");
    vc << "component\tvariance\tse\tpve\tse_pve\n";
    vc << std::scientific << std::setprecision(10);
    vc << "Vg\t" << res.sigma2_g << "\tNA\t" << res.pve << "\tNA\n";
    vc << "Ve\t" << res.sigma2_e << "\tNA\t" << (1.0 - res.pve) << "\tNA\n";

    std::ofstream summary(options_.out_prefix + ".vc.summary.txt");
    if (!summary) throw std::runtime_error("Cannot write --vc summary output");
    summary << "method\t" << res.method << "\n";
    summary << "n\t" << res.n << "\n";
    summary << "n_covars\t" << res.n_covars << "\n";
    summary << "n_snps\t" << res.n_snps << "\n";
    summary << "sigma2_g\t" << std::scientific << std::setprecision(10) << res.sigma2_g << "\n";
    summary << "sigma2_e\t" << res.sigma2_e << "\n";
    summary << "pve\t" << res.pve << "\n";
    summary << "delta\t" << res.delta << "\n";
}

void GemmaVcModel::run() {
    if (options_.mode < 1 || options_.mode > 3) {
        throw std::runtime_error("--vc must be 1, 2, or 3");
    }
    std::cout << "Running GEMMA-style variance component estimation (--vc "
              << options_.mode << ")...\n";
    Eigen::MatrixXd G = buildGrm();
    Eigen::VectorXd y = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X = DataFilter::subsetRows(X_, keep_);

    VcResult res = (options_.mode == 1) ? runHe(y, X, G) : runReml(y, X, G);
    res.n_snps = static_cast<int>(is_bgen_ ? bgen_reader_->getNumVariants()
                         : is_pgen_ ? pgen_reader_->getNumVariants()
                                    : plink_reader_->getNumSnps());
    writeResult(res);
    std::cout << "GEMMA-style VC output written to " << options_.out_prefix
              << ".vc.txt and " << options_.out_prefix << ".vc.summary.txt\n";
}

} // namespace cosmic
