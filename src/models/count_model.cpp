#include "models/count_model.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <omp.h>
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "pheno_reader.h"
#include "covar_reader.h"

namespace cosmic {

CountModel::Options CountModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.use_spa = cfg.use_spa;
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
    opts.offset_file = cfg.offset_file;
    opts.binomial_file = cfg.binomial_file;
    if (cfg.family_type == "negbin") opts.family = CountFamily::NegativeBinomial;
    else if (cfg.family_type == "binomial") opts.family = CountFamily::Binomial;
    else opts.family = CountFamily::Poisson;
    return opts;
}

CountModel::CountModel(const Options& opts) : options_(opts) {}

void CountModel::prepare(const std::vector<std::string>& fids,
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
    loadOffset();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for count model.");
    validateInputs();
}

void CountModel::loadPhenotype() {
    if (options_.family == CountFamily::Binomial) {
        if (options_.binomial_file.empty()) {
            throw std::runtime_error("Binomial PSR-GWAS requires --binomial-pheno with FID IID EVENT TRIAL columns.");
        }
        y_ = Eigen::VectorXd::Constant(n_samples_, std::numeric_limits<double>::quiet_NaN());
        trials_ = Eigen::VectorXd::Constant(n_samples_, std::numeric_limits<double>::quiet_NaN());

        std::unordered_map<std::string, int> sample_index;
        for (int i = 0; i < n_samples_; ++i) {
            sample_index[fids_[i] + ":" + iids_[i]] = i;
            sample_index[iids_[i]] = i;
        }

        std::ifstream ifs(options_.binomial_file);
        if (!ifs.is_open()) throw std::runtime_error("Cannot open --binomial-pheno: " + options_.binomial_file);
        std::string line;
        int loaded = 0;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string fid, iid, event_s, trial_s;
            if (!(ss >> fid >> iid >> event_s >> trial_s)) continue;
            char* end_event = nullptr;
            char* end_trial = nullptr;
            double events = std::strtod(event_s.c_str(), &end_event);
            double trials = std::strtod(trial_s.c_str(), &end_trial);
            if (end_event == event_s.c_str() || end_trial == trial_s.c_str()) continue;
            auto it = sample_index.find(fid + ":" + iid);
            if (it == sample_index.end()) it = sample_index.find(iid);
            if (it == sample_index.end()) continue;
            if (!std::isfinite(events) || !std::isfinite(trials) || trials <= 0.0 || events < 0.0 || events > trials) {
                throw std::runtime_error("Invalid binomial record: require finite 0 <= EVENT <= TRIAL and TRIAL > 0.");
            }
            y_(it->second) = events / trials;
            trials_(it->second) = trials;
            ++loaded;
        }
        if (loaded == 0) throw std::runtime_error("No valid binomial records loaded from --binomial-pheno.");
        std::cout << "Loaded " << loaded << " binomial event/trial records from " << options_.binomial_file << "\n";
        return;
    }
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
}

void CountModel::loadCovariates() {
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

void CountModel::loadOffset() {
    if (options_.offset_file.empty()) {
        offset_ = Eigen::VectorXd::Zero(y_.size());
        return;
    }
    offset_ = Eigen::VectorXd::Zero(y_.size());
    std::ifstream ifs(options_.offset_file);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open --offset file: " + options_.offset_file);
    }
    std::string line;
    int idx = 0;
    while (std::getline(ifs, line) && idx < offset_.size()) {
        const size_t line_number = static_cast<size_t>(idx) + 1;
        std::string trimmed = line;
        const size_t first = trimmed.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;  // tolerate blank lines
        trimmed = trimmed.substr(first, trimmed.find_last_not_of(" \t\r\n") - first + 1);
        try {
            size_t consumed = 0;
            const double value = std::stod(trimmed, &consumed);
            if (consumed != trimmed.size() || !std::isfinite(value))
                throw std::invalid_argument("not a single finite number");
            offset_(idx++) = value;
        } catch (const std::exception&) {
            throw std::runtime_error(
                "Invalid --offset file " + options_.offset_file + ": line " + std::to_string(line_number) +
                " must contain exactly one numeric log-exposure value (got \"" + trimmed +
                "\"). Provide one bare number per sample, in genotype order.");
        }
    }
    if (idx < offset_.size()) throw std::runtime_error("--offset file has fewer rows than samples.");
}

double CountModel::linkFun(double x) const {
    if (options_.family == CountFamily::Binomial) {
        if (x > 35.0) return 1.0 - options_.mu_eps;
        if (x < -35.0) return options_.mu_eps;
        return 1.0 / (1.0 + std::exp(-x));
    }
    return std::exp(x);
}

double CountModel::linkDeriv(double x) const {
    if (options_.family == CountFamily::Binomial) {
        double mu = linkFun(x);
        return std::max(options_.mu_eps, mu * (1.0 - mu));
    }
    return std::exp(x);
}

double CountModel::varianceFun(double mu) const {
    if (options_.family == CountFamily::NegativeBinomial) {
        return mu + overdispersion_ * mu * mu;
    }
    if (options_.family == CountFamily::Binomial) {
        return std::max(options_.mu_eps, mu * (1.0 - mu));
    }
    return mu;
}

void CountModel::validateInputs() const {
    if (options_.family == CountFamily::Poisson || options_.family == CountFamily::NegativeBinomial) {
        for (int idx : keep_.keep_indices) {
            double y = y_(idx);
            if (!std::isfinite(y) || y < 0.0) {
                throw std::runtime_error("Poisson/negative-binomial PSR-GWAS requires finite non-negative count phenotypes.");
            }
            if (std::abs(y - std::round(y)) > 1e-6) {
                std::cerr << "Warning: non-integer count phenotype detected; treating values as counts on the selected family scale.\n";
                break;
            }
        }
    }
    if (options_.family == CountFamily::Binomial) {
        double total_trials = 0.0, total_events = 0.0;
        for (int idx : keep_.keep_indices) {
            if (!std::isfinite(trials_(idx)) || trials_(idx) <= 0.0 || !std::isfinite(y_(idx)) || y_(idx) < 0.0 || y_(idx) > 1.0) {
                throw std::runtime_error("Binomial PSR-GWAS requires valid EVENT/TRIAL records after filtering.");
            }
            total_trials += trials_(idx);
            total_events += y_(idx) * trials_(idx);
        }
        double rate = total_events / std::max(1.0, total_trials);
        if (rate < 0.02 || rate > 0.98) {
            std::cerr << "Warning: extreme binomial event rate after filtering (rate=" << rate
                      << "); PSR calibration may be unstable.\n";
        }
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

void CountModel::runPQL() {
    std::cout << "Running PQL for count/binomial PSR model (family="
              << (options_.family == CountFamily::Poisson ? "poisson" :
                  (options_.family == CountFamily::NegativeBinomial ? "negbin" : "binomial")) << ")...\n";

    int n = keep_.n_keep;

    Eigen::MatrixXd G = buildGrmMatrix(*plink_reader_, *bgen_reader_, *pgen_reader_,
                                        is_bgen_, is_pgen_, keep_,
                                        options_.grm_alg, (float)options_.min_maf);

    Eigen::MatrixXd X_mat = X_.cast<double>();
    Eigen::VectorXd y_adj = DataFilter::subset(y_, keep_);
    Eigen::VectorXd off = DataFilter::subset(offset_, keep_);
    Eigen::VectorXd trials = (options_.family == CountFamily::Binomial)
        ? DataFilter::subset(trials_, keep_)
        : Eigen::VectorXd::Ones(n);
    X_mat = DataFilter::subsetRows(X_, keep_).cast<double>();
    if (G.rows() != n || G.cols() != n) G = DataFilter::subsetSquare(G, keep_);

    eta_ = Eigen::VectorXd::Zero(n);
    sigma2_g_ = 0.2;

    Eigen::VectorXd eta_old = eta_;
    std::vector<std::string> conv_rows;
    conv_rows.reserve(options_.pql_maxit);
    bool converged = false;

    for (int iter = 0; iter < options_.pql_maxit; ++iter) {
        int n_clipped = 0;
        mu_.resize(n);
        for (int i = 0; i < n; ++i) {
            mu_(i) = linkFun(eta_(i) + off(i));
            if (mu_(i) < options_.mu_eps) { mu_(i) = options_.mu_eps; ++n_clipped; }
            if (options_.family == CountFamily::Binomial && mu_(i) > 1.0 - options_.mu_eps) { mu_(i) = 1.0 - options_.mu_eps; ++n_clipped; }
        }

        working_y_.resize(n);
        weights_.resize(n);
        for (int i = 0; i < n; ++i) {
            double d = std::max(linkDeriv(eta_(i) + off(i)), options_.mu_eps);
            double v = varianceFun(mu_(i));
            working_y_(i) = eta_(i) + (y_adj(i) - mu_(i)) / d;
            double wy_raw = working_y_(i);
            working_y_(i) = std::max(-10.0, std::min(10.0, working_y_(i)));
            if (working_y_(i) != wy_raw) ++n_clipped;
            weights_(i) = (options_.family == CountFamily::Binomial)
                ? trials(i) * d * d / v
                : d * d / v;
            double w_raw = weights_(i);
            weights_(i) = std::max(options_.min_weight, std::min(options_.max_weight, weights_(i)));
            if (weights_(i) != w_raw) ++n_clipped;
        }

        Eigen::MatrixXd sqrtW = weights_.array().sqrt().matrix().asDiagonal();
        Eigen::VectorXd yw = sqrtW * working_y_;
        Eigen::MatrixXd Xw = sqrtW * X_mat;
        Eigen::MatrixXd Gw = sqrtW * G * sqrtW;

        EigenVce vce;
        EigenVce::Options vce_opts;
        vce_opts.use_covariate_projection = true;
        vce_opts.null_solver = EigenVce::NullSolver::Brent;
        vce_opts.max_iter = options_.reml_maxit;
        vce_opts.tol = options_.reml_tol;
        vce = EigenVce(vce_opts);
        vce.prepare(yw, Xw.cast<float>(), Gw);
        vce.runNullModel();

        double old_sigma2_g = sigma2_g_;
        sigma2_g_ = vce.getSigma2g();
        sigma2_e_ = vce.getSigma2e();

        if (std::isnan(sigma2_g_) || std::isnan(sigma2_e_) ||
            sigma2_g_ <= 0 || sigma2_e_ <= 0) {
            sigma2_g_ = 0.5 * old_sigma2_g + 0.1;
            sigma2_e_ = 0.5 * (sigma2_e_ > 0 ? sigma2_e_ : 1.0) + 0.1;
        }

        int p = Xw.cols();
        Eigen::MatrixXd V = sigma2_g_ * Gw + sigma2_e_ * Eigen::MatrixXd::Identity(n, n);
        Eigen::MatrixXd Vinv = V.ldlt().solve(Eigen::MatrixXd::Identity(n, n));
        Eigen::MatrixXd beta = (Xw.transpose() * Vinv * Xw).ldlt().solve(Xw.transpose() * Vinv * yw);
        Eigen::VectorXd resid = yw - Xw * beta;
        Eigen::VectorXd u_hat = sigma2_g_ * Gw * Vinv * resid;
        Eigen::VectorXd eta_w = Xw * beta + u_hat;
        Eigen::VectorXd eta_new = eta_w.array() / weights_.array().sqrt().max(options_.min_weight);

        // Damping to prevent oscillation: 50% old, 50% new
        eta_ = 0.5 * eta_ + 0.5 * eta_new;

        double change = (eta_ - eta_old).norm() / (eta_.norm() + 1e-10);
        eta_old = eta_;

        double vg_change = std::abs(sigma2_g_ - old_sigma2_g) / (sigma2_g_ + 1e-10);
        double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);
        std::cout << "PQL iter " << iter << ": Vg=" << sigma2_g_
                  << " Ve=" << sigma2_e_ << " h2=" << h2
                  << " eta_change=" << change << "\n";

        double quasi = 0.0;
        for (int i = 0; i < n; ++i) {
            if (options_.family == CountFamily::Binomial) {
                quasi += trials(i) * (y_adj(i) * std::log(mu_(i)) + (1.0 - y_adj(i)) * std::log(1.0 - mu_(i)));
            } else {
                quasi += y_adj(i) * std::log(mu_(i)) - mu_(i);
            }
        }
        std::ostringstream row;
        row << iter << "\t" << sigma2_g_ << "\t" << std::abs(sigma2_g_ - old_sigma2_g) << "\t" << quasi
            << "\t" << mu_.mean() << "\t" << mu_.minCoeff() << "\t" << mu_.maxCoeff()
            << "\t" << weights_.minCoeff() << "\t" << weights_.maxCoeff()
            << "\t" << n_clipped << "\t" << 0;
        conv_rows.push_back(row.str());

        if (change < options_.pql_tol && iter > 1) { converged = true; break; }
    }
    if (!conv_rows.empty() && converged) conv_rows.back().replace(conv_rows.back().size() - 1, 1, "1");
    writeConvergenceLog(conv_rows);
    if (converged) {
        std::cout << "PQL converged. Vg=" << sigma2_g_ << " Ve=" << sigma2_e_ << "\n";
    } else {
        std::cerr << "Warning: PQL reached max iterations without meeting tolerance. Final: Vg="
                  << sigma2_g_ << " Ve=" << sigma2_e_ << "\n";
    }
}

void CountModel::writeConvergenceLog(const std::vector<std::string>& rows) const {
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

void CountModel::scanGwas() {
    Eigen::VectorXd pseudo_y = Eigen::VectorXd::Constant(n_samples_, std::numeric_limits<double>::quiet_NaN());
    for (int i = 0; i < keep_.n_keep; ++i) pseudo_y(keep_.keep_indices[i]) = working_y_(i);

    GwasScanner::Options gopts;
    gopts.out_prefix = options_.out_prefix;
    gopts.min_maf = options_.min_maf;
    gopts.min_hwe = options_.min_hwe;
    gopts.fast_lmm = true;
    gopts.n_threads = options_.n_threads;
    gopts.no_loco = options_.no_loco;
    gopts.use_spa = options_.use_spa;
    gopts.grm_alg = options_.grm_alg;
    gopts.obs_weights = weights_.array().sqrt().matrix();

    GwasScanner scanner(gopts);
    scanner.runLocoScan(pseudo_y, X_, *plink_reader_, *bgen_reader_, *pgen_reader_, *vcf_reader_,
                        is_bgen_, is_pgen_, is_vcf_, keep_);
}

void CountModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_count_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) return;
    ofs << "family\tsigma2_g\tsigma2_e\n";
    ofs << (options_.family == CountFamily::Poisson ? "poisson" :
            (options_.family == CountFamily::NegativeBinomial ? "negbin" : "binomial"))
        << "\t" << sigma2_g_ << "\t" << sigma2_e_ << "\n";
    ofs.close();
    std::cout << "Count null model written to " << fname << "\n";
}

void CountModel::run() {
    runPQL();
    if (options_.psr_null_only) {
        writeNullSummary();
        std::string fname = options_.out_prefix + ".psr.null.tsv";
        std::ofstream ofs(fname);
        if (ofs.is_open()) {
            ofs << "field\tvalue\n";
            ofs << "model\tcount\n";
            std::string fam_str = (options_.family == CountFamily::Poisson ? "poisson" :
                                   (options_.family == CountFamily::NegativeBinomial ? "negbin" : "binomial"));
            ofs << "family\t" << fam_str << "\n";
            ofs << "n_samples\t" << keep_.n_keep << "\n";
            ofs << "sigma2_g\t" << sigma2_g_ << "\n";
            ofs << "sigma2_e\t" << sigma2_e_ << "\n";
            double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);
            ofs << "h2\t" << h2 << "\n";
            ofs << "overdispersion\t" << overdispersion_ << "\n";
            double quasi = 0.0;
            for (int i = 0; i < keep_.n_keep; ++i) {
                double yv = y_(keep_.keep_indices[i]);
                if (options_.family == CountFamily::Binomial) {
                    double t = (keep_.n_keep <= trials_.size()) ? trials_(keep_.keep_indices[i]) : 1.0;
                    quasi += t * (yv * std::log(std::max(mu_(i), 1e-15)) + (1.0 - yv) * std::log(std::max(1.0 - mu_(i), 1e-15)));
                } else {
                    quasi += yv * std::log(std::max(mu_(i), 1e-15)) - mu_(i);
                }
            }
            ofs << "quasi_loglik\t" << quasi << "\n";
            ofs << "mean_mu\t" << mu_.mean() << "\n";
            ofs << "min_mu\t" << mu_.minCoeff() << "\n";
            ofs << "max_mu\t" << mu_.maxCoeff() << "\n";
            ofs << "psr_null_only\tTRUE\n";
            ofs.close();
        }
        std::cout << "Count null model output written to " << fname << "\n";
        return;
    }
    scanGwas();
    writeNullSummary();
}

} // namespace cosmic
