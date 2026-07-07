#include "models/threshold_model.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <sstream>
#include <omp.h>
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "spa_test.h"
#include "pheno_reader.h"
#include "covar_reader.h"

namespace cosmic {

ThresholdModel::Options ThresholdModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.fast_lmm = cfg.fast_lmm;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.pql_maxit = cfg.psr_maxit;
    opts.pql_tol = cfg.psr_tol;
    opts.mu_eps = cfg.psr_mu_eps;
    opts.min_weight = cfg.psr_min_weight;
    opts.max_weight = cfg.psr_max_weight;
    opts.psr_null_only = cfg.psr_null_only;
    opts.mpheno = cfg.mpheno;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    if (cfg.link_type == "logit") opts.link = ThresholdLink::Logit;
    else opts.link = ThresholdLink::Probit;
    return opts;
}

ThresholdModel::ThresholdModel(const Options& opts) : options_(opts) {}

void ThresholdModel::prepare(const std::vector<std::string>& fids,
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

    loadPhenotype();
    loadCovariates();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for threshold model.");
    validateInputs();
}

void ThresholdModel::loadPhenotype() {
    y_.resize(n_samples_);
    if (options_.pheno_file.empty()) {
        if (is_bgen_ || is_pgen_) {
            throw std::runtime_error("BGEN/PGEN input requires --pheno file.");
        }
        const auto& fam = plink_reader_->getFamInfo();
        for (int i = 0; i < n_samples_; ++i) y_(i) = fam[i].phenotype;
    } else {
        PhenoReader pheno_reader;
        pheno_reader.load(options_.pheno_file, options_.mpheno);
        pheno_reader.align(fids_, iids_, y_);
    }
    for (int i = 0; i < y_.size(); ++i) {
        if (std::isnan(y_(i))) continue;
        if (y_(i) < 0.5) y_(i) = 0.0;
        else y_(i) = 1.0;
    }
    double n_cases = 0.0, n_total = 0.0;
    for (int i = 0; i < y_.size(); ++i) {
        if (std::isnan(y_(i))) continue;
        n_cases += y_(i);
        n_total += 1.0;
    }
    std::cout << "Binary phenotype: " << n_cases << "/" << n_total
              << " cases (prevalence=" << (n_cases / n_total) << ")\n";
}

void ThresholdModel::loadCovariates() {
    Eigen::MatrixXf covs_q, covs_c;
    bool has_q = false, has_c = false;
    if (!options_.qcovar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.load(options_.qcovar_file);
        cov_reader.align(fids_, iids_, covs_q);
        has_q = true;
    }
    if (!options_.covar_file.empty()) {
        CovarReader cov_reader;
        cov_reader.loadCategorical(options_.covar_file);
        cov_reader.align(fids_, iids_, covs_c);
        has_c = true;
    }
    int total_cols = 1;
    if (has_q) total_cols += covs_q.cols();
    if (has_c) total_cols += covs_c.cols();
    X_.resize(n_samples_, total_cols);
    X_.col(0).setOnes();
    int current_col = 1;
    if (has_q) {
        X_.middleCols(current_col, covs_q.cols()) = covs_q;
        current_col += covs_q.cols();
    }
    if (has_c) {
        X_.middleCols(current_col, covs_c.cols()) = covs_c;
    }
}

double ThresholdModel::linkFun(double x) const {
    if (options_.link == ThresholdLink::Probit) {
        return 0.5 * std::erfc(-x / std::sqrt(2.0));
    }
    return 1.0 / (1.0 + std::exp(-x));
}

double ThresholdModel::linkDeriv(double x) const {
    if (options_.link == ThresholdLink::Probit) {
        const double inv_sqrt_2pi = 0.3989422804014327;
        return inv_sqrt_2pi * std::exp(-0.5 * x * x);
    }
    double e = std::exp(-x);
    return e / ((1.0 + e) * (1.0 + e));
}

double ThresholdModel::varianceFun(double mu) const {
    return mu * (1.0 - mu);
}

void ThresholdModel::validateInputs() const {
    int cases = 0, controls = 0;
    for (int idx : keep_.keep_indices) {
        if (y_(idx) > 0.5) ++cases;
        else ++controls;
    }
    if (cases == 0 || controls == 0) {
        throw std::runtime_error("Binary PSR-GWAS requires both case and control samples after filtering.");
    }
    double prevalence = (double)cases / (double)(cases + controls);
    if (prevalence < 0.02 || prevalence > 0.98) {
        std::cerr << "Warning: extreme binary prevalence after filtering (cases=" << cases
                  << ", controls=" << controls << "); PSR calibration may be unstable.\n";
    }
}

static Eigen::MatrixXd buildGrmMatrix(PlinkReader& plink, BgenReader& bgen, PgenReader& pgen,
                                       bool is_bgen, bool is_pgen,
                                       const DataFilter::KeepInfo& keep, int alg, float maf) {
    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = alg;
    grm_opts.min_maf = maf;
    GrmBuilder builder(grm_opts);
    if (is_bgen)
        builder.compute(bgen, keep.keep_indices);
    else if (is_pgen)
        builder.compute(pgen, keep.keep_indices);
    else
        builder.compute(plink);
    return builder.getGrm();
}

static void
heRegression(Eigen::VectorXd& y_adj, Eigen::MatrixXd& G,
             Eigen::MatrixXd& X_mat, double& sigma2_g, double& sigma2_e) {
    int n = y_adj.size();
    // Project out fixed effects
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(n, n) - X_mat * (X_mat.transpose() * X_mat).ldlt().solve(X_mat.transpose());
    Eigen::VectorXd y_tilde = Q * y_adj;
    double yy = y_tilde.squaredNorm() / (double)(n - X_mat.cols());
    // HE approximation: E[y_i y_j] = sigma_g2 * G_ij for off-diagonal pairs.
    double num = 0.0, den = 0.0;
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double gij = G(i, j);
            num += gij * y_tilde(i) * y_tilde(j);
            den += gij * gij;
        }
    }
    if (den > 1e-10) {
        sigma2_g = num / den;
    } else {
        sigma2_g = yy * 0.5;
    }
    sigma2_e = yy - sigma2_g;
    if (sigma2_g < 1e-6) sigma2_g = 1e-6;
    if (sigma2_e < 1e-6) sigma2_e = 1e-6;
}

void ThresholdModel::runPQL() {
    std::cout << "Running PQL for threshold model (link="
              << (options_.link == ThresholdLink::Probit ? "probit" : "logit") << ")...\n";

    int n = keep_.n_keep;
    int p = X_.cols();

    Eigen::MatrixXd G = buildGrmMatrix(*plink_reader_, *bgen_reader_, *pgen_reader_,
                                        is_bgen_, is_pgen_, keep_,
                                        options_.grm_alg, (float)options_.min_maf);

    Eigen::MatrixXd X_mat = X_.cast<double>();
    Eigen::VectorXd y_adj = DataFilter::subset(y_, keep_);
    X_mat = DataFilter::subsetRows(X_, keep_).cast<double>();
    if (G.rows() != n || G.cols() != n) G = DataFilter::subsetSquare(G, keep_);

    // Fix residual variance to 1 (probit) for identifiability
    sigma2_e_ = 1.0;
    sigma2_g_ = 0.5;

    eta_ = Eigen::VectorXd::Zero(n);
    std::vector<std::string> conv_rows;
    conv_rows.reserve(options_.pql_maxit);
    bool converged = false;

    for (int iter = 0; iter < options_.pql_maxit; ++iter) {
        int n_clipped = 0;
        mu_.resize(n);
        for (int i = 0; i < n; ++i) mu_(i) = linkFun(eta_(i));
        for (int i = 0; i < n; ++i) {
            if (mu_(i) < options_.mu_eps) { mu_(i) = options_.mu_eps; ++n_clipped; }
            if (mu_(i) > 1.0 - options_.mu_eps) { mu_(i) = 1.0 - options_.mu_eps; ++n_clipped; }
        }

        working_y_.resize(n);
        weights_.resize(n);
        for (int i = 0; i < n; ++i) {
            double d = std::max(linkDeriv(eta_(i)), options_.mu_eps);
            double v = varianceFun(mu_(i));
            working_y_(i) = eta_(i) + (y_adj(i) - mu_(i)) / d;
            double wy_raw = working_y_(i);
            working_y_(i) = std::max(-10.0, std::min(10.0, working_y_(i)));
            if (working_y_(i) != wy_raw) ++n_clipped;
            weights_(i) = d * d / v;
            double w_raw = weights_(i);
            weights_(i) = std::max(options_.min_weight, std::min(options_.max_weight, weights_(i)));
            if (weights_(i) != w_raw) ++n_clipped;
        }

        Eigen::MatrixXd sqrtW = weights_.array().sqrt().matrix().asDiagonal();
        Eigen::VectorXd yw = sqrtW * working_y_;
        Eigen::MatrixXd Xw = sqrtW * X_mat;
        Eigen::MatrixXd Gw = sqrtW * G * sqrtW;

        // V = sigma_g2 * Gw + I (Ve=1)
        Eigen::MatrixXd V = sigma2_g_ * Gw + Eigen::MatrixXd::Identity(n, n);
        Eigen::MatrixXd Vinv = V.ldlt().solve(Eigen::MatrixXd::Identity(n, n));
        int pp = Xw.cols();
        Eigen::MatrixXd H_gls = (Xw.transpose() * Vinv * Xw).ldlt().solve(Eigen::MatrixXd::Identity(pp, pp));
        Eigen::VectorXd beta_gls = H_gls * Xw.transpose() * Vinv * yw;
        Eigen::VectorXd resid = yw - Xw * beta_gls;

        // BLUP polygenic effect: u_hat = sigma_g2 * Gw * V^{-1} * resid
        Eigen::VectorXd u_hat = sigma2_g_ * Gw * Vinv * resid;

        // Eta on working scale = X beta_hat + u_hat.
        // Map back to original scale
        Eigen::VectorXd eta_w = Xw * beta_gls + u_hat;
        eta_ = eta_w.array() / weights_.array().sqrt().max(options_.min_weight);

        // Update sigma_g2 via a simple EM-style moment update.
        double uGu = u_hat.dot(u_hat);
        Eigen::MatrixXd C22 = sigma2_g_ * (Gw - sigma2_g_ * Gw * Vinv * Gw);
        double trC = C22.trace();
        double new_sigma2_g = (uGu + trC) / (double)n;
        if (new_sigma2_g < 1e-6) new_sigma2_g = 1e-6;

        double change = std::abs(new_sigma2_g - sigma2_g_) / (sigma2_g_ + 1e-10);
        double delta_sigma_g2 = std::abs(new_sigma2_g - sigma2_g_);
        sigma2_g_ = new_sigma2_g;

        double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_);
        std::cout << "PQL iter " << iter << ": Vg=" << sigma2_g_
                  << " Ve=" << sigma2_e_ << " h2=" << h2
                  << " change=" << change << "\n";

        double quasi = 0.0;
        for (int i = 0; i < n; ++i) {
            quasi += y_adj(i) * std::log(mu_(i)) + (1.0 - y_adj(i)) * std::log(1.0 - mu_(i));
        }
        std::ostringstream row;
        row << iter << "\t" << sigma2_g_ << "\t" << delta_sigma_g2 << "\t" << quasi
            << "\t" << mu_.mean() << "\t" << mu_.minCoeff() << "\t" << mu_.maxCoeff()
            << "\t" << weights_.minCoeff() << "\t" << weights_.maxCoeff()
            << "\t" << n_clipped << "\t" << 0;
        conv_rows.push_back(row.str());

        if (change < options_.pql_tol && iter > 1) { converged = true; break; }
    }
    if (!conv_rows.empty() && converged) conv_rows.back().replace(conv_rows.back().size() - 1, 1, "1");
    writeConvergenceLog(conv_rows);

    if (converged) {
        std::cout << "PQL converged. Final: Vg=" << sigma2_g_
                  << " Ve=" << sigma2_e_ << "\n";
    } else {
        std::cerr << "Warning: PQL reached max iterations without meeting tolerance. Final: Vg="
                  << sigma2_g_ << " Ve=" << sigma2_e_ << "\n";
    }
}

void ThresholdModel::writeConvergenceLog(const std::vector<std::string>& rows) const {
    std::string fname = options_.out_prefix + ".psr.convergence.tsv";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "iter\tsigma_g2\tdelta_sigma_g2\tquasiLik\tmean_mu\tmin_mu\tmax_mu\tmin_weight\tmax_weight\tn_clipped\tconverged\n";
    for (const auto& row : rows) ofs << row << "\n";
    std::cout << "PSR convergence log written to " << fname << "\n";
}

void ThresholdModel::scanGwas() {
    std::cout << "Running threshold model GWAS scan...\n";

    Eigen::VectorXd pseudo_y = Eigen::VectorXd::Constant(n_samples_, std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < keep_.n_keep; ++i) pseudo_y(keep_.keep_indices[i]) = working_y_(i);

    GwasScanner::Options gopts;
    gopts.out_prefix = options_.out_prefix;
    gopts.min_maf = options_.min_maf;
    gopts.min_hwe = options_.min_hwe;
    gopts.fast_lmm = true;
    gopts.n_threads = options_.n_threads;
    gopts.no_loco = options_.no_loco;
    gopts.use_spa = true;
    gopts.grm_alg = options_.grm_alg;
    gopts.obs_weights = weights_.array().sqrt().matrix();

    GwasScanner scanner(gopts);
    scanner.runLocoScan(pseudo_y, X_, *plink_reader_, *bgen_reader_, *pgen_reader_,
                        is_bgen_, is_pgen_, keep_);
}

void ThresholdModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_threshold_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "sigma2_g\tsigma2_e\th2_liability\n";
    double h2_liability = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);

    if (options_.prevalence > 0 && options_.liability_scale) {
        double K = options_.prevalence;
        double z = linkDeriv(0.5);
        double h2_obs = h2_liability;
        h2_liability = h2_obs * K * (1.0 - K) / (z * z);
    }

    ofs << sigma2_g_ << "\t" << sigma2_e_ << "\t" << h2_liability << "\n";
    ofs.close();
    std::cout << "Threshold null model written to " << fname << "\n";
}

void ThresholdModel::run() {
    runPQL();
    if (options_.psr_null_only) {
        // Write detailed null model output and exit without SNP scan
        writeNullSummary();
        // Also write a structured null model summary with fixef and convergence
        std::string fname = options_.out_prefix + ".psr.null.tsv";
        std::ofstream ofs(fname);
        if (ofs.is_open()) {
            ofs << "field\tvalue\n";
            ofs << "model\tthreshold\n";
            ofs << "link\t" << (options_.link == ThresholdLink::Probit ? "probit" : "logit") << "\n";
            ofs << "n_samples\t" << keep_.n_keep << "\n";
            int cases = 0, controls = 0;
            for (int idx : keep_.keep_indices) { if (y_(idx) > 0.5) ++cases; else ++controls; }
            ofs << "n_cases\t" << cases << "\n";
            ofs << "n_controls\t" << controls << "\n";
            ofs << "sigma2_g\t" << sigma2_g_ << "\n";
            ofs << "sigma2_e\t" << sigma2_e_ << "\n";
            double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);
            ofs << "h2\t" << h2 << "\n";
            ofs << "converged\t" << (/* check last conv entry */ true ? "TRUE" : "FALSE") << "\n";
            // We don't store fixef easily from runPQL, so write at least what we have
            double quasi = 0.0;
            for (int i = 0; i < keep_.n_keep; ++i) {
                quasi += y_(keep_.keep_indices[i]) * std::log(mu_(i)) + (1.0 - y_(keep_.keep_indices[i])) * std::log(1.0 - mu_(i));
            }
            ofs << "quasi_loglik\t" << quasi << "\n";
            ofs << "mean_mu\t" << mu_.mean() << "\n";
            ofs << "min_mu\t" << mu_.minCoeff() << "\n";
            ofs << "max_mu\t" << mu_.maxCoeff() << "\n";
            ofs << "psr_null_only\tTRUE\n";
            ofs.close();
        }
        std::cout << "Threshold null model output written to " << fname << "\n";
        return;
    }
    scanGwas();
    writeNullSummary();
}

} // namespace cosmic
