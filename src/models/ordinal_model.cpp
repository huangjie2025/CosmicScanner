#include "models/ordinal_model.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <omp.h>

#include "models/threshold_model.h"
#include "grm_builder.h"
#include "pheno_reader.h"
#include "covar_reader.h"

namespace cosmic {

OrdinalModel::Options OrdinalModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.link_type = cfg.link_type;
    opts.pql_maxit = cfg.psr_maxit;
    opts.pql_tol = cfg.psr_tol;
    opts.mu_eps = cfg.psr_mu_eps;
    opts.min_weight = cfg.psr_min_weight;
    opts.max_weight = cfg.psr_max_weight;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.ordinal_mode = cfg.ordinal_mode;
    opts.ordinal_stacked_engine = cfg.ordinal_stacked_engine;
    opts.ordinal_vc_method = cfg.ordinal_vc_method;
    opts.ordinal_sigma2_g = cfg.ordinal_sigma2_g;
    opts.ordinal_vc_grid_min = cfg.ordinal_vc_grid_min;
    opts.ordinal_vc_grid_max = cfg.ordinal_vc_grid_max;
    opts.ordinal_vc_grid_points = cfg.ordinal_vc_grid_points;
    opts.ordinal_allow_fallback = cfg.ordinal_allow_fallback;
    opts.ordinal_test = cfg.ordinal_test;
    opts.pheno_file = cfg.pheno_file;
    opts.ordinal_file = cfg.ordinal_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    return opts;
}

OrdinalModel::OrdinalModel(const Options& opts) : options_(opts) {}

void OrdinalModel::prepare(const std::vector<std::string>& fids,
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
    loadOrdinalPhenotype();
    loadCovariates();
    validateInputs();
}

void OrdinalModel::loadOrdinalPhenotype() {
    std::string path = !options_.ordinal_file.empty() ? options_.ordinal_file : options_.pheno_file;
    if (path.empty()) {
        throw std::runtime_error("Ordinal PSR-GWAS requires --ordinal-pheno or --pheno with FID IID ORDINAL columns.");
    }

    y_ = Eigen::VectorXd::Constant(n_samples_, std::numeric_limits<double>::quiet_NaN());
    std::unordered_map<std::string, int> sample_index;
    for (int i = 0; i < n_samples_; ++i) {
        sample_index[fids_[i] + ":" + iids_[i]] = i;
        sample_index[iids_[i]] = i;
    }

    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open ordinal phenotype file: " + path);
    std::string line;
    int loaded = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string fid, iid, value_s;
        if (!(ss >> fid >> iid >> value_s)) continue;
        char* end = nullptr;
        double value = std::strtod(value_s.c_str(), &end);
        if (end == value_s.c_str() || !std::isfinite(value)) continue;
        auto it = sample_index.find(fid + ":" + iid);
        if (it == sample_index.end()) it = sample_index.find(iid);
        if (it == sample_index.end()) continue;
        y_(it->second) = value;
        ++loaded;
    }
    if (loaded == 0) throw std::runtime_error("No valid ordinal records loaded.");

    std::set<double> unique_levels;
    for (int i = 0; i < y_.size(); ++i) {
        if (!std::isfinite(y_(i))) continue;
        unique_levels.insert(y_(i));
        level_counts_[y_(i)] += 1;
    }
    levels_.assign(unique_levels.begin(), unique_levels.end());
    std::cout << "Loaded " << loaded << " ordinal records with " << levels_.size() << " levels.\n";
}

void OrdinalModel::loadCovariates() {
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
    int total_cols = 0;
    if (has_q) total_cols += covs_q.cols();
    if (has_c) total_cols += covs_c.cols();
    if (total_cols == 0) { has_covariates_ = false; return; }
    has_covariates_ = true;
    X_covariates_.resize(n_samples_, total_cols);
    int col = 0;
    if (has_q) {
        X_covariates_.middleCols(col, covs_q.cols()) = covs_q.cast<double>();
        col += covs_q.cols();
    }
    if (has_c) {
        X_covariates_.middleCols(col, covs_c.cols()) = covs_c.cast<double>();
    }
}

void OrdinalModel::validateInputs() const {
    if (levels_.size() < 3) {
        throw std::runtime_error("Ordinal PSR-GWAS requires at least 3 ordered categories after filtering.");
    }
    for (double level : levels_) {
        auto it = level_counts_.find(level);
        int count = (it == level_counts_.end()) ? 0 : it->second;
        if (count < 5) {
            std::cerr << "Warning: ordinal level " << level << " has only " << count
                      << " samples; PSR calibration may be unstable.\n";
        }
    }
}

std::string OrdinalModel::writeThresholdPhenotype(double threshold, int threshold_index) const {
    std::string path = options_.out_prefix + ".ordinal.t" + std::to_string(threshold_index) + ".pheno";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write temporary ordinal threshold phenotype: " + path);
    for (int i = 0; i < n_samples_; ++i) {
        if (!std::isfinite(y_(i))) continue;
        double binary = (y_(i) <= threshold) ? 1.0 : 0.0;
        out << fids_[i] << "\t" << iids_[i] << "\t" << binary << "\n";
    }
    return path;
}

void OrdinalModel::writeThresholdSummary() const {
    std::string path = options_.out_prefix + ".ordinal.thresholds.tsv";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write ordinal thresholds: " + path);
    out << "threshold_index\tthreshold_value\tn_leq\tn_gt\tordinal_mode\n";
    for (size_t k = 0; k + 1 < levels_.size(); ++k) {
        double threshold = levels_[k];
        int n_leq = 0, n_gt = 0;
        for (int i = 0; i < y_.size(); ++i) {
            if (!std::isfinite(y_(i))) continue;
            if (y_(i) <= threshold) ++n_leq;
            else ++n_gt;
        }
        out << (k + 1) << "\t" << threshold << "\t" << n_leq << "\t" << n_gt
            << "\texperimental_cumulative_binary\n";
    }
}

void OrdinalModel::writeCombinedConvergence(const std::vector<std::string>& threshold_prefixes) const {
    std::ofstream out(options_.out_prefix + ".psr.convergence.tsv");
    out << "threshold_index\titer\tsigma_g2\tdelta_sigma_g2\tquasiLik\tmean_mu\tmin_mu\tmax_mu\tmin_weight\tmax_weight\tn_clipped\tconverged\n";
    for (size_t k = 0; k < threshold_prefixes.size(); ++k) {
        std::ifstream in(threshold_prefixes[k] + ".psr.convergence.tsv");
        std::string line;
        bool header = true;
        while (std::getline(in, line)) {
            if (header) { header = false; continue; }
            if (!line.empty()) out << (k + 1) << "\t" << line << "\n";
        }
    }
}

void OrdinalModel::writeCombinedGwas(const std::vector<std::string>& threshold_prefixes, double runtime_sec) const {
    struct BestRow {
        double p = 2.0;
        double beta = std::numeric_limits<double>::quiet_NaN();
        double se = std::numeric_limits<double>::quiet_NaN();
        double log10p = std::numeric_limits<double>::quiet_NaN();
        int threshold_index = 0;
    };
    std::map<std::string, BestRow> best;
    bool all_converged = true;

    for (size_t k = 0; k < threshold_prefixes.size(); ++k) {
        std::ifstream conv(threshold_prefixes[k] + ".psr.convergence.tsv");
        std::string line, last;
        std::getline(conv, line);
        while (std::getline(conv, line)) if (!line.empty()) last = line;
        if (last.empty() || last.substr(last.find_last_of('\t') + 1) != "1") all_converged = false;

        std::ifstream in(threshold_prefixes[k] + ".assoc.txt");
        if (!in.is_open()) throw std::runtime_error("Cannot open ordinal threshold association file: " + threshold_prefixes[k] + ".assoc.txt");
        std::string header;
        std::getline(in, header);
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string chr, snp, bp, a1, a2, freq, n, beta_s, se_s, p_s, logp_s;
            ss >> chr >> snp >> bp >> a1 >> a2 >> freq >> n >> beta_s >> se_s >> p_s >> logp_s;
            double p = std::stod(p_s);
            auto& row = best[snp];
            if (p < row.p) {
                row.p = p;
                row.beta = std::stod(beta_s);
                row.se = std::stod(se_s);
                row.log10p = std::stod(logp_s);
                row.threshold_index = static_cast<int>(k + 1);
            }
        }
    }

    std::ofstream out(options_.out_prefix + ".gwas.tsv");
    out << "trait_type\tordinal_mode\tsnp\teffect\tse\twald_or_stat\tp_value\tlog10p\tn_levels\tn_thresholds\tbest_threshold\tconverged\tstatus_class\truntime\n";
    std::string status = all_converged ? "OK" : "CONVERGED_WITH_WARNING";
    for (const auto& kv : best) {
        const auto& row = kv.second;
        double wald = (row.se > 0.0) ? (row.beta / row.se) * (row.beta / row.se) : std::numeric_limits<double>::quiet_NaN();
        out << "ordinal\texperimental_cumulative_binary\t" << kv.first << "\t"
            << row.beta << "\t" << row.se << "\t" << wald << "\t" << row.p << "\t"
            << row.log10p << "\t" << levels_.size() << "\t" << threshold_prefixes.size() << "\t"
            << row.threshold_index << "\t" << (all_converged ? "TRUE" : "FALSE") << "\t"
            << status << "\t" << runtime_sec << "\n";
    }
}

// --- Stacked Ordinal PSR Implementation ---

Eigen::MatrixXd OrdinalModel::buildGrmMatrix() const {
    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = (float)options_.min_maf;
    GrmBuilder builder(grm_opts);
    if (is_bgen_)
        builder.compute(*bgen_reader_);
    else if (is_pgen_)
        builder.compute(*pgen_reader_);
    else
        builder.compute(*plink_reader_);
    return builder.getGrm();
}

void OrdinalModel::buildStackedData(int n_thresh, int n_stacked,
                                    Eigen::VectorXd& z_stack,
                                    Eigen::MatrixXd& X_stack,
                                    Eigen::VectorXi& row2indiv) const {
    int n = n_samples_;
    z_stack.resize(n_stacked);
    X_stack.resize(n_stacked, n_thresh);
    row2indiv.resize(n_stacked);
    X_stack.setZero();

    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n_thresh; ++k) {
            int r = i * n_thresh + k;
            z_stack(r) = (y_(i) <= levels_[k]) ? 1.0 : 0.0;
            row2indiv(r) = i;
            X_stack(r, k) = 1.0;
        }
    }

    if (has_covariates_) {
        int n_cov = X_covariates_.cols();
        Eigen::MatrixXd X_ext(n_stacked, n_thresh + n_cov);
        X_ext.leftCols(n_thresh) = X_stack;
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n_thresh; ++k) {
                int r = i * n_thresh + k;
                X_ext.row(r).tail(n_cov) = X_covariates_.row(i);
            }
        }
        X_stack = X_ext;
    }
}

void OrdinalModel::writeStackedThresholds(int n_thresh) const {
    std::string path = options_.out_prefix + ".ordinal.thresholds.tsv";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write stacked ordinal thresholds: " + path);
    out << "threshold_index\tthreshold_value\tn_leq\tn_gt\tordinal_mode\n";
    for (int k = 0; k < n_thresh; ++k) {
        double threshold = levels_[k];
        int n_leq = 0, n_gt = 0;
        for (int i = 0; i < y_.size(); ++i) {
            if (!std::isfinite(y_(i))) continue;
            if (y_(i) <= threshold) ++n_leq;
            else ++n_gt;
        }
        out << (k + 1) << "\t" << threshold << "\t" << n_leq << "\t" << n_gt
            << "\tstacked_psr\n";
    }
}

void OrdinalModel::writeStackedConvergence(const std::vector<std::string>& rows) const {
    std::string fname = options_.out_prefix + ".psr.convergence.tsv";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "iter\tsigma_g2\tdelta_sigma_g2\tquasiLik\tmean_mu\tmin_mu\tmax_mu\tmin_weight\tmax_weight\tn_clipped\tordinal_stacked_engine\tvc_method\tvc_update\tconverged\n";
    for (const auto& row : rows) ofs << row << "\n";
    std::cout << "Stacked PSR convergence log written to " << fname << "\n";
}

void OrdinalModel::writeStackedBlup(const Eigen::VectorXd& u_hat,
                                    const Eigen::VectorXd& sigma2_g) const {
    std::string path = options_.out_prefix + ".ordinal.blup.tsv";
    std::ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("Cannot write " + path);
    out << "FID\tIID\tblup\tsigma2_g\n";
    for (int i = 0; i < n_samples_; ++i) {
        out << fids_[i] << "\t" << iids_[i] << "\t"
            << u_hat(i) << "\t" << sigma2_g(i) << "\n";
    }
    std::cout << "Stacked PSR BLUP written to " << path << "\n";
}

void OrdinalModel::writeStackedVc(const std::vector<std::string>& vc_rows) const {
    std::string path = options_.out_prefix + ".ordinal.vc.tsv";
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "Warning: cannot write " << path << "\n";
        return;
    }
    out << "iter\tsigma2_g\treml_loglik\tboundary\tgrid_min\tgrid_max\tgrid_points\tn_grid_valid\tjitter_used\tjitter_value\tstatus\n";
    for (const auto& row : vc_rows) out << row << "\n";
    std::cout << "Stacked PSR VC log written to " << path << "\n";
}

OrdinalModel::GridRemlResult OrdinalModel::estimateSigma2GridReml(
    const Eigen::VectorXd& y_tilde,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& K_stack,
    const Eigen::VectorXd& w,
    double current_sigma2) const
{
    int N = y_tilde.size();
    int p = X.cols();
    GridRemlResult result;
    result.sigma2_g = current_sigma2;
    result.reml_loglik = -std::numeric_limits<double>::infinity();
    result.boundary = false;
    result.n_valid = 0;

    Eigen::MatrixXd R_inv = w.asDiagonal();
    int n_points = options_.ordinal_vc_grid_points;
    double g_min = std::max(options_.ordinal_vc_grid_min, 0.0);
    double g_max = std::max(options_.ordinal_vc_grid_max, g_min + 1e-6);
    double step = (g_max - g_min) / (std::max(n_points - 1, 1));

    // Precompute X'R^{-1} and y'R^{-1} terms that don't depend on sigma2_g
    Eigen::MatrixXd RtX = R_inv * X;
    Eigen::MatrixXd XtR_1_X = X.transpose() * RtX;
    Eigen::VectorXd R_inv_y = R_inv * y_tilde;
    Eigen::VectorXd XtR_1_y = X.transpose() * R_inv_y;

    // Precompute K_stack quantities that don't depend on sigma2_g
    Eigen::MatrixXd KR = K_stack * R_inv;
    Eigen::MatrixXd KX = K_stack * X;
    Eigen::VectorXd Ky = K_stack * R_inv_y;

    for (int gi = 0; gi < n_points; ++gi) {
        double g = g_min + gi * step;
        if (n_points == 1) g = (g_min + g_max) / 2.0;

        double jitter = 0.0;
        bool jittered = false;
        Eigen::MatrixXd V = g * K_stack;
        for (int i = 0; i < N; ++i) V(i, i) += 1.0 / w(i);

        Eigen::LDLT<Eigen::MatrixXd> ldlt;
        ldlt.compute(V);
        bool ok = (ldlt.info() == Eigen::Success);
        if (!ok) {
            // Try jitter
            for (double j : {1e-8, 1e-6, 1e-4}) {
                for (int i = 0; i < N; ++i) V(i, i) += j;
                ldlt.compute(V);
                if (ldlt.info() == Eigen::Success) { jitter = j; jittered = true; ok = true; break; }
                for (int i = 0; i < N; ++i) V(i, i) -= j;
            }
        }
        if (!ok) continue;

        ++result.n_valid;

        // log|V|
        Eigen::VectorXd d = ldlt.vectorD();
        double log_det_V = 0.0;
        for (int i = 0; i < d.size(); ++i) {
            if (d(i) > 0) log_det_V += std::log(d(i));
            else { ok = false; break; }
        }
        if (!ok) continue;

        // X'V^{-1}X = X'R^{-1}X - g * (KX)' V^{-1} (KX) ... actually use ldlt
        // V = g*K + R^{-1}, so V^{-1} = R - g*R*K*V^{-1} (using Woodbury)
        // Simpler: compute X'V^{-1}X directly
        Eigen::MatrixXd Vinv_X = ldlt.solve(X);
        Eigen::MatrixXd XtVinvX = X.transpose() * Vinv_X;

        // log|X'V^{-1}X|
        Eigen::LDLT<Eigen::MatrixXd> ldlt_xtvx;
        ldlt_xtvx.compute(XtVinvX);
        Eigen::VectorXd d2 = ldlt_xtvx.vectorD();
        double log_det_XtVX = 0.0;
        for (int i = 0; i < d2.size(); ++i) {
            if (d2(i) > 0) log_det_XtVX += std::log(d2(i));
            else { ok = false; break; }
        }
        if (!ok) continue;

        // y'Py = y'V^{-1}y - y'V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}y
        Eigen::VectorXd Vinv_y = ldlt.solve(y_tilde);
        double ytPy = y_tilde.dot(Vinv_y) - Vinv_y.transpose() * X * ldlt_xtvx.solve(X.transpose() * Vinv_y);

        double reml_ll = -0.5 * (log_det_V + log_det_XtVX + ytPy);
        if (std::isfinite(reml_ll) && reml_ll > result.reml_loglik) {
            result.sigma2_g = g;
            result.reml_loglik = reml_ll;
            result.boundary = (g <= options_.ordinal_vc_grid_min + 1e-10);
            result.jitter_used = jittered;
            result.jitter_value = jitter;
            result.status = "OK";
        }
    }

    if (result.n_valid == 0) {
        throw std::runtime_error("Grid REML failed: all sigma2_g candidates produced non-finite V");
    }
    return result;
}

void OrdinalModel::runStackedPsr() {
    auto start = std::chrono::steady_clock::now();
    int C = static_cast<int>(levels_.size());
    int n_thresh = C - 1;
    int n = n_samples_;
    int N = n * n_thresh;

    if (N > 5000) {
        std::cout << "Stacked PSR: " << N << " rows (n=" << n << ", thresholds=" << n_thresh
                  << "). Dense first pass.\n";
    }
    std::cout << "Running stacked PSR with " << C << " levels, "
              << n_thresh << " thresholds, " << N << " stacked rows.\n";

    // 1. Build stacked data
    Eigen::VectorXd z;
    Eigen::MatrixXd X;
    Eigen::VectorXi row2indiv;
    buildStackedData(n_thresh, N, z, X, row2indiv);
    int p = X.cols();

    // 2. Write thresholds and build GRM
    writeStackedThresholds(n_thresh);
    Eigen::MatrixXd G = buildGrmMatrix();

    Eigen::MatrixXd K_stack(N, N);
    for (int r1 = 0; r1 < N; ++r1) {
        int i1 = row2indiv(r1);
        for (int r2 = 0; r2 < N; ++r2) {
            int i2 = row2indiv(r2);
            K_stack(r1, r2) = G(i1, i2);
        }
    }

    // 3. Initialize alpha from cumulative proportions (logit scale)
    Eigen::VectorXd theta(p);
    theta.setZero();
    for (int k = 0; k < n_thresh; ++k) {
        double prop = 0.0;
        for (int i = 0; i < n; ++i) {
            if (y_(i) <= levels_[k]) prop += 1.0;
        }
        prop /= n;
        double p_clamped = std::max(std::min(prop, 1.0 - 1e-10), 1e-10);
        theta(k) = std::log(p_clamped / (1.0 - p_clamped));
    }

    // 4. Determine engine
    bool use_dense = (options_.ordinal_stacked_engine == "dense-zkz");
    if (use_dense && N > 2000) {
        throw std::runtime_error("dense-zkz engine requires N_stack <= 2000 (current N="
                                 + std::to_string(N) + "). Use fixed-gls engine instead.");
    }
    std::string engine_label = use_dense ? "dense_zkz_blup_first_pass" : "fixed_vc_weighted_gls_first_pass";
    std::cout << "Stacked PSR engine: " << options_.ordinal_stacked_engine
              << " (label: " << engine_label << ")\n";

    if (options_.ordinal_test != "wald" && options_.ordinal_test != "score" &&
        options_.ordinal_test != "robust-score" && options_.ordinal_test != "all") {
        throw std::runtime_error("Invalid --ordinal-test: " + options_.ordinal_test +
                                 " (expected wald, score, robust-score, or all)");
    }

    double sigma2_g = options_.ordinal_sigma2_g;
    std::string vc_method = options_.ordinal_vc_method;
    bool vc_update = false;
    if (vc_method == "grid-reml") {
        if (!use_dense) {
            throw std::runtime_error("grid-reml VC update requires --ordinal-stacked-engine dense-zkz");
        }
        vc_update = true;
        std::cout << "Stacked PSR VC method: grid-reml (sigma2_g range [" << options_.ordinal_vc_grid_min
                  << ", " << options_.ordinal_vc_grid_max << "], " << options_.ordinal_vc_grid_points << " points)\n";
    } else {
        std::cout << "Stacked PSR VC method: fixed (sigma2_g=" << sigma2_g << ")\n";
    }
    Eigen::VectorXd eta = X * theta;

    // 5. PSR/PQL iteration
    bool converged = false;
    std::vector<std::string> conv_rows;
    std::vector<std::string> vc_rows;
    conv_rows.reserve(options_.pql_maxit);
    Eigen::VectorXd w(N), y_tilde(N);
    Eigen::LDLT<Eigen::MatrixXd> ldlt_V;

    double prev_sigma2_g = sigma2_g;
    double prev_quasi = -std::numeric_limits<double>::infinity();

    Eigen::VectorXd mu(N);
    for (int iter = 0; iter < options_.pql_maxit; ++iter) {
        int n_clipped = 0;
        for (int i = 0; i < N; ++i) {
            mu(i) = 1.0 / (1.0 + std::exp(-eta(i)));
            if (mu(i) < options_.mu_eps) { mu(i) = options_.mu_eps; ++n_clipped; }
            if (mu(i) > 1.0 - options_.mu_eps) { mu(i) = 1.0 - options_.mu_eps; ++n_clipped; }
        }
        for (int i = 0; i < N; ++i) {
            double d = std::max(mu(i) * (1.0 - mu(i)), options_.mu_eps);
            y_tilde(i) = eta(i) + (z(i) - mu(i)) / d;
            double raw = y_tilde(i);
            y_tilde(i) = std::max(-10.0, std::min(10.0, y_tilde(i)));
            if (y_tilde(i) != raw) ++n_clipped;
            w(i) = std::max(options_.min_weight, std::min(options_.max_weight, d));
        }

        // VC update (grid-reml) before GLS
        if (vc_update) {
            GridRemlResult gr = estimateSigma2GridReml(y_tilde, X, K_stack, w, sigma2_g);
            double dsg = std::abs(gr.sigma2_g - prev_sigma2_g) / (gr.sigma2_g + 1e-10);
            sigma2_g = gr.sigma2_g;
            prev_sigma2_g = sigma2_g;
            std::ostringstream vc_row;
            vc_row << iter << "\t" << gr.sigma2_g << "\t" << gr.reml_loglik
                   << "\t" << (gr.boundary ? "TRUE" : "FALSE")
                   << "\t" << options_.ordinal_vc_grid_min
                   << "\t" << options_.ordinal_vc_grid_max
                   << "\t" << options_.ordinal_vc_grid_points
                   << "\t" << gr.n_valid
                   << "\t" << (gr.jitter_used ? "TRUE" : "FALSE")
                   << "\t" << gr.jitter_value
                   << "\t" << gr.status;
            vc_rows.push_back(vc_row.str());
        }

        Eigen::VectorXd beta_gls(p);
        if (use_dense) {
            Eigen::MatrixXd V = sigma2_g * K_stack;
            for (int i = 0; i < N; ++i) V(i, i) += 1.0 / w(i);
            ldlt_V.compute(V);
            if (ldlt_V.info() != Eigen::Success) {
                throw std::runtime_error("LDLT of V failed in PSR iteration " + std::to_string(iter));
            }
            Eigen::MatrixXd Vinv_X = ldlt_V.solve(X);
            Eigen::MatrixXd XtVinvX = X.transpose() * Vinv_X;
            XtVinvX.diagonal().array() += 1e-6;
            Eigen::VectorXd Vinv_y = ldlt_V.solve(y_tilde);
            Eigen::VectorXd XtViny = X.transpose() * Vinv_y;
            Eigen::MatrixXd H = XtVinvX.ldlt().solve(Eigen::MatrixXd::Identity(p, p));
            beta_gls = H * XtViny;
        } else {
            Eigen::MatrixXd sqrtW = w.array().sqrt().matrix().asDiagonal();
            Eigen::VectorXd yw = sqrtW * y_tilde;
            Eigen::MatrixXd Xw = sqrtW * X;
            Eigen::MatrixXd XtX = Xw.transpose() * Xw;
            XtX.diagonal().array() += 1e-6;
            Eigen::MatrixXd H = XtX.ldlt().solve(Eigen::MatrixXd::Identity(p, p));
            beta_gls = H * Xw.transpose() * yw;
        }

        if (iter > 0) {
            double damp = 0.5;
            beta_gls = damp * beta_gls + (1.0 - damp) * theta;
        }
        for (int j = 0; j < p; ++j) {
            beta_gls(j) = std::max(-10.0, std::min(10.0, beta_gls(j)));
        }
        eta = X * beta_gls;
        theta = beta_gls;

        if (eta.hasNaN()) {
            throw std::runtime_error("eta has NaN after GLS update in stacked PSR");
        }
        double change = 1.0;
        if (iter > 0) {
            Eigen::VectorXd eta_prev = eta;
            change = (eta - eta_prev).norm() / (eta.norm() + 1e-10);
        }

        double quasi = 0.0;
        for (int i = 0; i < N; ++i) {
            double m = std::max(std::min(mu(i), 1.0 - 1e-15), 1e-15);
            quasi += z(i) * std::log(m) + (1.0 - z(i)) * std::log(1.0 - m);
        }

        double delta_sg = (iter == 0) ? 0.0
            : std::abs(sigma2_g - prev_sigma2_g) / (sigma2_g + 1e-10);

        std::ostringstream row;
        row << iter << "\t" << sigma2_g << "\t" << delta_sg
            << "\t" << quasi
            << "\t" << mu.mean() << "\t" << mu.minCoeff() << "\t" << mu.maxCoeff()
            << "\t" << w.minCoeff() << "\t" << w.maxCoeff()
            << "\t" << n_clipped
            << "\t" << engine_label
            << "\t" << vc_method
            << "\t" << (vc_update ? "TRUE" : "FALSE")
            << "\t" << 0;
        conv_rows.push_back(row.str());

        double dq = std::abs(quasi - prev_quasi) / (std::abs(quasi) + 1e-10);
        prev_quasi = quasi;
        prev_sigma2_g = sigma2_g;

        if (change < options_.pql_tol && iter > 0) { converged = true; break; }
    }
    if (!conv_rows.empty() && converged) {
        conv_rows.back().replace(conv_rows.back().size() - 1, 1, "1");
    }
    writeStackedConvergence(conv_rows);
    if (!vc_rows.empty()) {
        vc_rows.back().replace(vc_rows.back().rfind('\t') + 1, std::string::npos, "DONE");
        writeStackedVc(vc_rows);
    }

    std::cout << "Stacked PSR (" << engine_label << ") "
              << (converged ? "converged" : "finished")
              << ". Final Vg=" << sigma2_g
              << ", vc_method=" << vc_method
              << ", vc_update=" << (vc_update ? "TRUE" : "FALSE") << "\n";

    // 6. Precompute for SNP scan
    Eigen::VectorXd beta_hat(p);
    Eigen::VectorXd Py(N);
    // fixed-gls precomputed quantities (used in scan, declared here for scope)
    Eigen::VectorXd yw_scaled(N);
    Eigen::MatrixXd XtVinvX_inv(p, p);
    Eigen::MatrixXd Xw_scaled(N, p);
    Eigen::VectorXd XtViny_fixed(p);
    Eigen::MatrixXd Vinv_X;
    if (use_dense) {
        Vinv_X = ldlt_V.solve(X);
        Eigen::VectorXd Vinv_y = ldlt_V.solve(y_tilde);
        Eigen::MatrixXd XtVinvX = X.transpose() * Vinv_X;
        XtVinvX_inv = XtVinvX.ldlt().solve(Eigen::MatrixXd::Identity(p, p));
        Eigen::VectorXd XtViny = X.transpose() * Vinv_y;
        beta_hat = XtVinvX_inv * XtViny;
        Py = Vinv_y - Vinv_X * (XtVinvX_inv * XtViny);
    } else {
        Eigen::MatrixXd sqrtW = w.array().sqrt().matrix().asDiagonal();
        Xw_scaled = sqrtW * X;
        yw_scaled = sqrtW * y_tilde;
        Eigen::MatrixXd XtVinvX = Xw_scaled.transpose() * Xw_scaled;
        XtVinvX_inv = XtVinvX.ldlt().solve(Eigen::MatrixXd::Identity(p, p));
        XtViny_fixed = Xw_scaled.transpose() * yw_scaled;
        beta_hat = XtVinvX_inv * XtViny_fixed;
        Py = yw_scaled - Xw_scaled * (XtVinvX_inv * XtViny_fixed);
    }

    // 6b. Score test precomputation
    z_stack_ = z;
    mu_ = mu;
    w_ = w;
    row2indiv_ = row2indiv;

    bool do_score_diag = (options_.ordinal_test != "wald");
    bool do_raw = do_score_diag;
    bool do_cluster = (options_.ordinal_test == "robust-score" || options_.ordinal_test == "all");

    Eigen::VectorXd r_resid_fixed;
    Eigen::VectorXd Pr_dense;
    Eigen::VectorXd Vinv_r;
    Eigen::VectorXd XtVinv_r;
    std::ofstream score_out;

    if (do_score_diag) {
        Eigen::VectorXd r_raw = z - mu;
        if (use_dense) {
            Vinv_r = ldlt_V.solve(r_raw);
            XtVinv_r = X.transpose() * Vinv_r;
            Pr_dense = Vinv_r - Vinv_X * (XtVinvX_inv * XtVinv_r);
        } else {
            Eigen::MatrixXd sqrtW = w.array().sqrt().matrix().asDiagonal();
            Eigen::VectorXd rw = sqrtW * r_raw;
            Eigen::VectorXd Xw_t_rw = Xw_scaled.transpose() * rw;
            Eigen::VectorXd br = XtVinvX_inv * Xw_t_rw;
            r_resid_fixed = rw - Xw_scaled * br;
        }

        std::string score_path = options_.out_prefix + ".ordinal.score.tsv";
        score_out.open(score_path);
        if (!score_out.is_open()) throw std::runtime_error("Cannot write " + score_path);
        score_out << "trait_type\tordinal_mode\tordinal_implementation\t"
                  << "ordinal_stacked_engine\t"
                  << "ordinal_vc_method\tvc_update\tsigma2_g\t"
                  << "n_levels\tn_thresholds\t"
                  << "chr\tsnp\tbp\ta1\ta2\tfreq\t"
                  << "n_stack\tn_clusters\t"
                  << "wald_stat\twald_p\twald_log10p\t"
                  << "psr_score_stat\tpsr_score_p\tpsr_score_log10p\t"
                  << "raw_score_stat\traw_score_p\traw_score_log10p\t"
                  << "cluster_score_stat\tcluster_score_p\tcluster_score_log10p\t"
                  << "abs_wald_minus_psr_score\twald_psr_score_equivalent\t"
                  << "raw_score_projection\t"
                  << "converged\tstatus_class\tstatus\terror_message\n";
    }

    int n_snps = 0;
    if (is_pgen_) n_snps = static_cast<int>(pgen_reader_->getNumVariants());
    else if (is_bgen_) n_snps = static_cast<int>(bgen_reader_->getNumVariants());
    else n_snps = static_cast<int>(plink_reader_->getNumSnps());

    std::string out_path = options_.out_prefix + ".gwas.tsv";
    std::ofstream out(out_path);
    if (!out.is_open()) throw std::runtime_error("Cannot write " + out_path);
    out << "trait_type\tordinal_mode\tordinal_implementation\tordinal_fallback_used\t"
        << "vc_method\tvc_update\tsigma2_g\t"
        << "n_levels\tn_thresholds\t"
        << "chr\tsnp\tbp\ta1\ta2\tfreq\teffect\tse\twald_or_stat\tp_value\tlog10p\tassociation_test\t"
        << "converged\tstatus_class\truntime\terror_message\n";
    std::string status_str = converged ? "OK" : "CONVERGED_WITH_WARNING";

    int n_skipped_lowmaf = 0;
    double runtime_sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    for (int si = 0; si < n_snps; ++si) {
        std::string chr, snp_id, a1, a2;
        int bp = 0;
        double freq = 0.0;
        Eigen::VectorXd dosage_vec;
        bool read_ok = false;

        if (is_pgen_) {
            std::string rsid;
            uint32_t pos = 0;
            Eigen::VectorXf dosage_f;
            read_ok = pgen_reader_->readVariant(chr, rsid, pos, a1, a2, dosage_f);
            if (read_ok) {
                snp_id = rsid;
                bp = static_cast<int>(pos);
                dosage_vec = dosage_f.cast<double>();
            }
        } else if (is_bgen_) {
            std::string rsid;
            uint32_t pos = 0;
            Eigen::VectorXf dosage_f;
            read_ok = bgen_reader_->readVariant(chr, rsid, pos, a1, a2, dosage_f);
            if (read_ok) {
                snp_id = rsid;
                bp = static_cast<int>(pos);
                dosage_vec = dosage_f.cast<double>();
            }
        } else {
            const auto& bim = plink_reader_->getBimInfo();
            chr = bim[si].chrom;
            snp_id = bim[si].id;
            bp = static_cast<int>(bim[si].bp_pos);
            a1 = bim[si].alt;
            a2 = bim[si].ref;
            Eigen::VectorXf geno_f;
            size_t n_read = plink_reader_->getSnpGenotypes(static_cast<size_t>(si), geno_f);
            read_ok = (n_read > 0);
            if (read_ok) dosage_vec = geno_f.cast<double>();
        }

        if (!read_ok) { ++n_skipped_lowmaf; continue; }

        double sum_g = 0.0, n_g = 0.0;
        for (int i = 0; i < dosage_vec.size(); ++i) {
            if (std::isfinite(dosage_vec(i))) { sum_g += dosage_vec(i); n_g += 1.0; }
        }
        if (n_g < 1.0) { ++n_skipped_lowmaf; continue; }
        freq = sum_g / (2.0 * n_g);
        if (freq < options_.min_maf || freq > 1.0 - options_.min_maf) { ++n_skipped_lowmaf; continue; }

        // Build stacked genotype
        Eigen::VectorXd g_stacked(N);
        for (int i = 0; i < n; ++i) {
            double g = (i < dosage_vec.size() && std::isfinite(dosage_vec(i))) ? dosage_vec(i) : 0.0;
            for (int k = 0; k < n_thresh; ++k) {
                g_stacked(i * n_thresh + k) = g;
            }
        }

        // Wald test and score diagnostics
        double gtPg = 0.0, gtPy = 0.0;
        Eigen::VectorXd Vg, XtVg;
        Eigen::VectorXd gw, Xw_t_gw;
        if (use_dense) {
            Vg = ldlt_V.solve(g_stacked);
            XtVg = X.transpose() * Vg;
            gtPg = g_stacked.dot(Vg) - XtVg.dot(XtVinvX_inv * XtVg);
            gtPy = g_stacked.dot(Py);
        } else {
            Eigen::MatrixXd sqrtW = w.array().sqrt().matrix().asDiagonal();
            gw = sqrtW * g_stacked;
            Xw_t_gw = Xw_scaled.transpose() * gw;
            gtPg = gw.squaredNorm() - Xw_t_gw.dot(XtVinvX_inv * Xw_t_gw);
            gtPy = gw.dot(yw_scaled) - Xw_t_gw.dot(XtVinvX_inv * XtViny_fixed);
        }

        if (gtPg <= 1e-12) { ++n_skipped_lowmaf; continue; }

        double gamma = gtPy / gtPg;
        double se = std::sqrt(1.0 / gtPg);
        double wald_stat = gtPy * gtPy / gtPg;
        double wald = wald_stat;
        double z = std::abs(gtPy) / std::sqrt(gtPg);
        double p_val = std::erfc(z / std::sqrt(2.0));
        double log10p = (p_val > 0.0) ? -std::log10(p_val) : 15.0;

        // Score diagnostics
        double psr_score_stat = wald_stat;
        double psr_score_p = p_val;
        double psr_score_log10p = log10p;

        double raw_score_stat = -1, raw_score_p = -1, raw_score_log10p = -1;
        double cluster_score_stat = -1, cluster_score_p = -1, cluster_score_log10p = -1;
        double abs_wald_minus_psr = std::abs(wald_stat - psr_score_stat);
        bool wald_psr_equivalent = (abs_wald_minus_psr < 1e-10);

        if (do_raw) {
            double U_raw;
            if (use_dense) {
                double gtPr = g_stacked.dot(Vinv_r) - XtVg.dot(XtVinvX_inv * XtVinv_r);
                U_raw = gtPr;
            } else {
                U_raw = gw.dot(r_resid_fixed);
            }
            raw_score_stat = (U_raw * U_raw) / gtPg;
            if (raw_score_stat > 0) {
                double z_raw = std::sqrt(raw_score_stat);
                raw_score_p = std::erfc(z_raw / std::sqrt(2.0));
                raw_score_log10p = (raw_score_p > 0) ? -std::log10(raw_score_p) : 15.0;
            }
        }

        if (do_cluster) {
            std::vector<double> S(n, 0.0);
            if (use_dense) {
                Eigen::VectorXd bg_dense = XtVinvX_inv * XtVg;
                Eigen::VectorXd Pg = Vg - Vinv_X * bg_dense;
                for (int r = 0; r < N; ++r) {
                    S[row2indiv_(r)] += Pg(r) * Pr_dense(r);
                }
            } else {
                Eigen::VectorXd bg = XtVinvX_inv * Xw_t_gw;
                for (int r = 0; r < N; ++r) {
                    double P_gw_r = gw(r) - Xw_scaled.row(r).dot(bg);
                    S[row2indiv_(r)] += P_gw_r * r_resid_fixed(r);
                }
            }
            double U_cl = 0, V_cl = 0;
            for (int i = 0; i < n; ++i) {
                U_cl += S[i];
                V_cl += S[i] * S[i];
            }
            cluster_score_stat = (V_cl > 0) ? (U_cl * U_cl / V_cl) : -1;
            if (cluster_score_stat > 0) {
                double z_cl = std::sqrt(cluster_score_stat);
                cluster_score_p = std::erfc(z_cl / std::sqrt(2.0));
                cluster_score_log10p = (cluster_score_p > 0) ? -std::log10(cluster_score_p) : 15.0;
            }
        }

        std::string association_test = "wald";
        double output_stat = wald_stat;
        double output_p = p_val;
        double output_log10p = log10p;
        if (options_.ordinal_test == "score") {
            association_test = "psr_score";
            output_stat = psr_score_stat;
            output_p = psr_score_p;
            output_log10p = psr_score_log10p;
        } else if (options_.ordinal_test == "robust-score" && cluster_score_stat > 0) {
            association_test = "cluster_score";
            output_stat = cluster_score_stat;
            output_p = cluster_score_p;
            output_log10p = cluster_score_log10p;
        }

        out << "ordinal\tstacked_psr\t" << engine_label << "\tFALSE\t"
            << vc_method << "\t" << (vc_update ? "TRUE" : "FALSE") << "\t" << sigma2_g << "\t"
            << C << "\t" << n_thresh << "\t"
            << chr << "\t" << snp_id << "\t" << bp << "\t"
            << a1 << "\t" << a2 << "\t" << freq << "\t"
            << gamma << "\t" << se << "\t" << output_stat << "\t" << output_p << "\t" << output_log10p << "\t"
            << association_test << "\t"
            << (converged ? "TRUE" : "FALSE") << "\t" << status_str << "\t"
            << runtime_sec << "\t\n";

        if (do_score_diag) {
            std::string raw_proj = use_dense ? "dense_v" : "wls";
            score_out << "ordinal\tstacked_psr\t" << engine_label << "\t"
                      << options_.ordinal_stacked_engine << "\t"
                      << vc_method << "\t" << (vc_update ? "TRUE" : "FALSE") << "\t" << sigma2_g << "\t"
                      << C << "\t" << n_thresh << "\t"
                      << chr << "\t" << snp_id << "\t" << bp << "\t"
                      << a1 << "\t" << a2 << "\t" << freq << "\t"
                       << N << "\t" << n << "\t"
                      << wald_stat << "\t" << p_val << "\t" << log10p << "\t"
                      << psr_score_stat << "\t" << psr_score_p << "\t" << psr_score_log10p << "\t"
                      << raw_score_stat << "\t" << raw_score_p << "\t" << raw_score_log10p << "\t"
                      << cluster_score_stat << "\t" << cluster_score_p << "\t" << cluster_score_log10p << "\t"
                      << abs_wald_minus_psr << "\t"
                      << (wald_psr_equivalent ? "TRUE" : "FALSE") << "\t"
                      << raw_proj << "\t"
                      << (converged ? "TRUE" : "FALSE") << "\t" << status_str << "\t"
                      << "OK\t\n";
        }
    }

    if (n_skipped_lowmaf > 0) {
        std::cout << "Stacked PSR: skipped " << n_skipped_lowmaf << " variants (low MAF or missing).\n";
    }
    std::cout << "Stacked PSR results written to " << out_path << "\n";

    if (do_score_diag) {
        score_out.close();
        std::string score_path = options_.out_prefix + ".ordinal.score.tsv";
        std::cout << "Stacked PSR score diagnostics written to " << score_path << "\n";
    }

    // 7. BLUP output (dense-zkz only)
    if (use_dense && converged) {
        Eigen::VectorXd resid = y_tilde - X * beta_hat;
        Eigen::VectorXd Vinv_resid = ldlt_V.solve(resid);
        Eigen::VectorXd z_vinv_resid = Eigen::VectorXd::Zero(n);
        for (int r = 0; r < N; ++r) {
            z_vinv_resid(row2indiv(r)) += Vinv_resid(r);
        }
        Eigen::VectorXd u_hat = sigma2_g * G.selfadjointView<Eigen::Lower>() * z_vinv_resid;
        Eigen::VectorXd sigma2_g_vec = Eigen::VectorXd::Constant(n, sigma2_g);
        writeStackedBlup(u_hat, sigma2_g_vec);
    }
}

void OrdinalModel::run() {
    auto start = std::chrono::steady_clock::now();
    bool fallback_used = false;

    if (options_.ordinal_mode == "stacked-psr") {
        try {
            runStackedPsr();
            auto end = std::chrono::steady_clock::now();
            double runtime = std::chrono::duration<double>(end - start).count();
            std::cout << "Ordinal PSR-GWAS stacked-psr results written to "
                      << options_.out_prefix << ".gwas.tsv (runtime=" << runtime << "s)\n";
            return;
        } catch (const std::exception& e) {
            if (options_.ordinal_allow_fallback) {
                std::cerr << "WARNING: stacked-psr failed (" << e.what()
                          << "). Falling back to experimental cumulative-binary.\n";
                fallback_used = true;
            } else {
                throw std::runtime_error(
                    std::string("stacked-psr failed: ") + e.what()
                    + ". Use --ordinal-allow-fallback to fall back to experimental cumulative-binary.");
            }
        }
    }

    // Experimental cumulative-binary path (unchanged)
    writeThresholdSummary();
    std::vector<std::string> threshold_prefixes;
    for (size_t k = 0; k + 1 < levels_.size(); ++k) {
        double threshold = levels_[k];
        std::string pheno = writeThresholdPhenotype(threshold, static_cast<int>(k + 1));
        std::string prefix = options_.out_prefix + ".ordinal.t" + std::to_string(k + 1);
        threshold_prefixes.push_back(prefix);

        ThresholdModel::Options opts;
        opts.out_prefix = prefix;
        opts.min_maf = options_.min_maf;
        opts.min_hwe = options_.min_hwe;
        opts.grm_alg = options_.grm_alg;
        opts.n_threads = options_.n_threads;
        opts.no_loco = options_.no_loco;
        opts.fast_lmm = options_.fast_lmm;
        opts.link = (options_.link_type == "probit") ? ThresholdLink::Probit : ThresholdLink::Logit;
        opts.pql_maxit = options_.pql_maxit;
        opts.pql_tol = options_.pql_tol;
        opts.mu_eps = options_.mu_eps;
        opts.min_weight = options_.min_weight;
        opts.max_weight = options_.max_weight;
        opts.reml_maxit = options_.reml_maxit;
        opts.reml_tol = options_.reml_tol;
        opts.mpheno = 1;
        opts.pheno_file = pheno;
        opts.qcovar_file = options_.qcovar_file;
        opts.covar_file = options_.covar_file;

        std::string mode_label = fallback_used ? "stacked-psr (fallback to cumulative-binary)" : "experimental cumulative-binary";
        std::cout << "Ordinal cumulative threshold " << (k + 1) << "/" << (levels_.size() - 1)
                  << " (y <= " << threshold << ") using " << mode_label << " PSR.\n";
        ThresholdModel model(opts);
        model.prepare(fids_, iids_, n_samples_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_);
        model.run();
    }
    auto end = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    writeCombinedConvergence(threshold_prefixes);
    writeCombinedGwas(threshold_prefixes, runtime);
    std::cout << "Ordinal PSR-GWAS "
              << (fallback_used ? "stacked-psr (fallback) " : "experimental cumulative-binary ")
              << "results written to " << options_.out_prefix << ".gwas.tsv\n";
}

} // namespace cosmic
