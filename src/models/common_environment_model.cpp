#include "models/common_environment_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <set>
#include <chrono>
#include <iomanip>
#include <stdexcept>
#include <omp.h>
#include "pheno_reader.h"
#include "covar_reader.h"
#include "grm_builder.h"
#include "repeatability_helpers.h"
#include "qc_utils.h"
#include "matrix_io.h"

namespace cosmic {

namespace {
// Extract the submatrix of G corresponding to kept subjects (rank >= 0).
// Ranks are assigned in ascending sample order, so index order is preserved.
Eigen::MatrixXd subsetByRank(const Eigen::MatrixXd& G, const std::vector<int>& subject_rank) {
    std::vector<int> keep;
    for (int s = 0; s < (int)subject_rank.size(); ++s) {
        if (subject_rank[s] >= 0) keep.push_back(s);
    }
    Eigen::MatrixXd out(keep.size(), keep.size());
    for (size_t i = 0; i < keep.size(); ++i) {
        for (size_t j = 0; j < keep.size(); ++j) {
            out(i, j) = G(keep[i], keep[j]);
        }
    }
    return out;
}

void reportGroupDiagnostics(int n_skipped, const std::vector<int>& group_record_count) {
    if (n_skipped > 0) {
        std::cout << "Skipped " << n_skipped
                  << " phenotype records with no common-environment group "
                     "(excluded from model).\n";
    }
    int n_groups = (int)group_record_count.size();
    if (n_groups == 0) return;
    int gmin = *std::min_element(group_record_count.begin(), group_record_count.end());
    int gmax = *std::max_element(group_record_count.begin(), group_record_count.end());
    double gmean = std::accumulate(group_record_count.begin(), group_record_count.end(), 0.0) / n_groups;
    int n_small = 0;
    for (int c : group_record_count) if (c < 2) ++n_small;
    std::cout << "Group sizes: min=" << gmin << " max=" << gmax
              << " mean=" << gmean << " (" << n_small << " groups with <2 records)\n";
}
} // namespace

CommonEnvironmentModel::Options CommonEnvironmentModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.tpheno = cfg.tpheno_set ? cfg.tpheno : 0;
    opts.reml_vg = cfg.reml_vg;
    opts.reml_vc = cfg.reml_vc;
    opts.reml_ve = cfg.reml_ve;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.common_env_file = cfg.common_env_file;
    opts.group_file_record_level = cfg.common_env_record_level;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    if (cfg.num_threads > 0) {
        opts.n_threads = cfg.num_threads;
        omp_set_num_threads(cfg.num_threads);
    }
    return opts;
}

CommonEnvironmentModel::CommonEnvironmentModel(const Options& opts) : options_(opts) {}

// ---------------------------------------------------------------------------
// prepare
// ---------------------------------------------------------------------------

void CommonEnvironmentModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              VcfReader& vcf_reader,
                              bool is_bgen, bool is_pgen, bool is_vcf) {
    fids_ = fids; iids_ = iids; n_samples_ = n_samples;
    is_bgen_ = is_bgen; is_pgen_ = is_pgen;
    is_vcf_ = is_vcf;
    plink_reader_ = &plink_reader; bgen_reader_ = &bgen_reader; pgen_reader_ = &pgen_reader;
    vcf_reader_ = &vcf_reader;

    loadPhenotype();
    loadCommonEnvMapping();
    loadCovariates();
    buildGRM();
}

void CommonEnvironmentModel::loadPhenotype() {
    if (options_.pheno_file.empty()) {
        throw std::runtime_error(
            "Common-environment model requires --pheno file (long format: FID IID [TIME] VALUE)");
    }
    std::cout << "Loading phenotype records for common-environment model...\n";
    PhenoReader pheno_reader;
    auto records = (options_.tpheno > 0)
        ? pheno_reader.loadLongitudinalRecords(options_.pheno_file, options_.tpheno, options_.mpheno)
        : pheno_reader.loadAllRecords(options_.pheno_file, options_.mpheno);

    std::unordered_map<std::string, int> fid_iid_to_gen_idx;
    for (int i = 0; i < n_samples_; ++i) {
        fid_iid_to_gen_idx[fids_[i] + ":" + iids_[i]] = i;
    }

    gen_records_.clear();
    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = fid_iid_to_gen_idx.find(key);
        if (it == fid_iid_to_gen_idx.end()) continue;
        GenRecord gr;
        gr.aid = it->second + 1;
        gr.y = rec.value;
        gr.time = rec.time;
        gr.idstr = key;
        gen_records_.push_back(std::move(gr));
    }
    n_recs_kept_ = (int)gen_records_.size();
    std::cout << "Kept " << n_recs_kept_ << " phenotype records matching genotype samples.\n";
    if (n_recs_kept_ < 10) throw std::runtime_error("Too few records for common-environment model.");

    std::set<int> subjects;
    for (const auto& rec : gen_records_) subjects.insert(rec.aid - 1);
    n_subjects_kept_ = (int)subjects.size();
}

void CommonEnvironmentModel::loadCommonEnvMapping() {
    if (options_.common_env_file.empty()) {
        throw std::runtime_error(
            "Common-environment model requires --common-env-file "
            "(FID IID GROUP_ID per subject, or per record with --common-env-record-level)");
    }
    std::ifstream ifs(options_.common_env_file);
    if (!ifs.is_open()) {
        throw std::runtime_error("Cannot open common-environment file: " + options_.common_env_file);
    }

    std::string line;
    int line_no = 0;
    int n_skipped_unknown_subject = 0;
    int n_skipped_missing = 0;
    int n_header_skipped = 0;

    std::set<std::string> known_subjects;
    for (int s = 0; s < n_samples_; ++s) {
        known_subjects.insert(fids_[s] + ":" + iids_[s]);
    }

    std::vector<std::vector<std::string>> group_lines;
    if (!options_.group_file_record_level) {
        // Subject-level mapping: FID IID GROUP, one line per subject.
        std::unordered_map<std::string, int> subject_to_group;
        std::unordered_map<std::string, int> label_to_idx;
        while (std::getline(ifs, line)) {
            ++line_no;
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::vector<std::string> tok;
            std::string t;
            while (iss >> t) tok.push_back(t);
            if (tok.size() < 3) {
                if (line_no == 1) { ++n_header_skipped; continue; }
                ++n_skipped_missing;
                std::cerr << "Warning: common-env line " << line_no
                          << " has <3 columns; skipped.\n";
                continue;
            }
            std::string key = tok[0] + ":" + tok[1];
            if (known_subjects.find(key) == known_subjects.end()) {
                ++n_skipped_unknown_subject;
                continue;
            }
            // Assign the group level by GROUP LABEL (not by subject).
            std::string gid = tok[2];
            auto gidx = label_to_idx.find(gid);
            if (gidx == label_to_idx.end()) {
                gidx = label_to_idx.emplace(gid, (int)group_labels_.size()).first;
                group_labels_.push_back(gid);
            }
            if (!subject_to_group.emplace(key, gidx->second).second) {
                std::cerr << "Warning: duplicate common-env entry for " << key
                          << " (line " << line_no << "); keeping first.\n";
            }
        }
        if (n_header_skipped > 0) {
            std::cout << "Skipped header line in common-env file.\n";
        }
        if (n_skipped_unknown_subject > 0) {
            std::cout << "Skipped " << n_skipped_unknown_subject
                      << " common-env entries for subjects absent from genotype data.\n";
        }
        if (group_labels_.empty()) {
            throw std::runtime_error("Common-environment file contains no valid group lines.");
        }
        n_groups_ = (int)group_labels_.size();
        obs_to_group_.clear();
        obs_to_group_.reserve(n_recs_kept_);
        for (int r = 0; r < n_recs_kept_; ++r) {
            auto it = subject_to_group.find(gen_records_[r].idstr);
            if (it == subject_to_group.end()) {
                ++n_skipped_missing;
                obs_to_group_.push_back(-1);
            } else {
                obs_to_group_.push_back(it->second);
            }
        }
    } else {
        // Record-level mapping: FID IID GROUP, one line per phenotype record
        // (same order as the phenotype file's genotype-matched records).
        // Optional TIME column: FID IID TIME GROUP.
        while (std::getline(ifs, line)) {
            ++line_no;
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::vector<std::string> tok;
            std::string t;
            while (iss >> t) tok.push_back(t);
            if (tok.size() < 3) {
                if (line_no == 1) { ++n_header_skipped; continue; }
                ++n_skipped_missing;
                std::cerr << "Warning: common-env record line " << line_no
                          << " has <3 columns; skipped.\n";
                continue;
            }
            group_lines.push_back(std::move(tok));
        }
        int n_mismatch = 0;
        obs_to_group_.clear();
        obs_to_group_.reserve(n_recs_kept_);
        std::unordered_map<std::string, int> label_to_idx;
        for (int r = 0; r < n_recs_kept_; ++r) {
            const auto& rec = gen_records_[r];
            if (r < (int)group_lines.size()) {
                const auto& tok = group_lines[r];
                std::string file_id = tok[0] + ":" + tok[1];
                if (file_id != rec.idstr) {
                    ++n_mismatch;
                    std::cerr << "Warning: common-env record line " << (r + 1)
                              << " ID '" << file_id << "' does not match phenotype record '"
                              << rec.idstr << "'\n";
                    obs_to_group_.push_back(-1);
                    continue;
                }
                std::string gid = tok[2];
                auto gidx = label_to_idx.find(gid);
                if (gidx == label_to_idx.end()) {
                    gidx = label_to_idx.emplace(gid, (int)group_labels_.size()).first;
                    group_labels_.push_back(gid);
                }
                obs_to_group_.push_back(gidx->second);
            } else {
                ++n_skipped_missing;
                obs_to_group_.push_back(-1);
            }
        }
        if (n_mismatch > 0) {
            std::cerr << "Warning: " << n_mismatch << "/" << n_recs_kept_
                      << " record-level common-env lines mismatched phenotype order.\n";
        }
        n_groups_ = (int)group_labels_.size();
        if (n_groups_ == 0) {
            throw std::runtime_error("Common-environment file contains no valid group assignments.");
        }
    }

    if (n_skipped_missing > 0) {
        std::cout << n_skipped_missing << " phenotype records have no common-environment group; "
                     "they are excluded from the model.\n";
    }

    // Drop records without a group (shared by both mapping modes)
    std::vector<GenRecord> kept;
    std::vector<int> kept_obs_to_group;
    kept.reserve(n_recs_kept_);
    kept_obs_to_group.reserve(n_recs_kept_);
    for (int r = 0; r < n_recs_kept_; ++r) {
        if (obs_to_group_[r] >= 0) {
            kept.push_back(gen_records_[r]);
            kept_obs_to_group.push_back(obs_to_group_[r]);
        }
    }
    int dropped = n_recs_kept_ - (int)kept.size();
    if (dropped > 0) std::cout << "Dropped " << dropped << " records without group assignment.\n";
    gen_records_ = std::move(kept);
    obs_to_group_ = std::move(kept_obs_to_group);
    n_recs_kept_ = (int)gen_records_.size();
    if (n_recs_kept_ < 10) throw std::runtime_error("Too few records with group assignment.");

    std::vector<int> group_record_count(n_groups_, 0);
    for (int g : obs_to_group_) ++group_record_count[g];
    reportGroupDiagnostics(0, group_record_count);
}

void CommonEnvironmentModel::loadCovariates() {
    Eigen::MatrixXf covs_q, covs_c;
    bool has_q = !options_.qcovar_file.empty();
    bool has_c = !options_.covar_file.empty();
    if (has_q) {
        CovarReader cov_reader;
        cov_reader.load(options_.qcovar_file);
        cov_reader.align(fids_, iids_, covs_q);
    }
    if (has_c) {
        CovarReader cov_reader;
        cov_reader.loadCategorical(options_.covar_file);
        cov_reader.align(fids_, iids_, covs_c);
    }

    int p = 1 + (has_q ? covs_q.cols() : 0) + (has_c ? covs_c.cols() : 0);
    fd_.p = p;
    fd_.rows.resize(n_recs_kept_);
    for (int i = 0; i < n_recs_kept_; ++i) {
        int gi = gen_records_[i].aid - 1;
        fd_.rows[i].push_back({0, 1.0});
        int col = 1;
        if (has_q) {
            for (int k = 0; k < covs_q.cols(); ++k) {
                float v = covs_q(gi, k);
                if (!std::isnan(v)) fd_.rows[i].push_back({col + k, (double)v});
            }
            col += covs_q.cols();
        }
        if (has_c) {
            for (int k = 0; k < covs_c.cols(); ++k) {
                float v = covs_c(gi, k);
                if (!std::isnan(v)) fd_.rows[i].push_back({col + k, (double)v});
            }
        }
    }

    X_design_.resize(n_recs_kept_, p);
    for (int i = 0; i < n_recs_kept_; ++i) {
        X_design_.row(i).setZero();
        for (const auto& kv : fd_.rows[i]) {
            if (kv.first >= 0 && kv.first < p) X_design_(i, kv.first) = kv.second;
        }
    }

    // Drop aliased fixed-effect columns (rank-deficient design).
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X_design_);
    qr.setThreshold(1e-10);
    const int rank = qr.rank();
    if (rank < p) {
        Eigen::VectorXi piv = qr.colsPermutation().indices();
        std::vector<int> keep_cols;
        for (int i = 0; i < rank; ++i) keep_cols.push_back(piv[i]);
        std::sort(keep_cols.begin(), keep_cols.end());
        std::vector<int> old_to_new(p, -1);
        for (int i = 0; i < (int)keep_cols.size(); ++i) old_to_new[keep_cols[i]] = i;
        Eigen::MatrixXd X_new(n_recs_kept_, rank);
        for (int j = 0; j < rank; ++j) X_new.col(j) = X_design_.col(keep_cols[j]);
        X_design_ = std::move(X_new);
        for (auto& row : fd_.rows) {
            std::vector<std::pair<int, double>> new_row;
            for (const auto& kv : row) {
                if (kv.first >= 0 && kv.first < (int)old_to_new.size() && old_to_new[kv.first] >= 0) {
                    new_row.push_back({old_to_new[kv.first], kv.second});
                }
            }
            row.swap(new_row);
        }
        fd_.p = rank;
        std::cout << "Dropped " << (p - rank) << " aliased fixed-effect column(s).\n";
    }

    y_obs_.resize(n_recs_kept_);
    for (int i = 0; i < n_recs_kept_; ++i) y_obs_(i) = gen_records_[i].y;
}

void CommonEnvironmentModel::buildGRM() {
    std::cout << "Calculating genomic GRM for common-environment model...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
#ifdef _OPENMP
    grm_opts.thread_num = omp_get_max_threads();
#endif
    GrmBuilder builder(grm_opts);
    if (is_pgen_) { pgen_reader_->reset(); builder.compute(*pgen_reader_); }
    else if (is_bgen_) { bgen_reader_->reset(); builder.compute(*bgen_reader_); }
    else { builder.compute(*plink_reader_); }
    G_subject_ = builder.getGrm().cast<double>();
    std::cout << "GRM built: " << G_subject_.rows() << "x" << G_subject_.cols() << "\n";
}

RandomEffectSpec CommonEnvironmentModel::additiveSpec() const {
    RandomEffectSpec spec;
    spec.name = "additive";
    spec.kernel = RandomEffectKernel::GRM;
    spec.obs_to_level = obs_to_sample_;
    spec.n_levels = n_subjects_kept_;
    spec.vc_index = 0;
    return spec;
}

RandomEffectSpec CommonEnvironmentModel::commonEnvSpec() const {
    RandomEffectSpec spec;
    spec.name = "common";
    spec.kernel = RandomEffectKernel::Identity;
    spec.obs_to_level = obs_to_group_;
    spec.n_levels = n_groups_;
    spec.vc_index = 1;
    return spec;
}

// ---------------------------------------------------------------------------
// VCE: record-level REML, 3 components (sigma2_a, sigma2_c, sigma2_e)
// ---------------------------------------------------------------------------

double CommonEnvironmentModel::recordLogL(const Eigen::MatrixXd& G_used,
                                          double vg, double vc, double ve) {
    if (!(vg >= 0.0) || !(vc >= 0.0) || !(ve > 0.0)) {
        return -std::numeric_limits<double>::infinity();
    }
    if (!std::isfinite(vg) || !std::isfinite(vc) || !std::isfinite(ve)) {
        return -std::numeric_limits<double>::infinity();
    }
    if ((int)obs_to_sample_.size() != n_recs_kept_ || (int)obs_to_group_.size() != n_recs_kept_) {
        return -std::numeric_limits<double>::infinity();
    }
    std::vector<RandomEffectSpec> specs = {additiveSpec(), commonEnvSpec()};
    try {
        RecordLevelNullModel nm;
        if (!nm.prepare(y_obs_, X_design_, specs, {vg, vc}, ve, &G_used)) {
            return -std::numeric_limits<double>::infinity();
        }
        double log_det_v = nm.logDetV();
        double log_det_x = 0.0;
        // logL = -0.5 (log|V| + log|X'V^-1X| + y'Py + (n-p)log(2pi))
        // log|X'V^-1X| via factorization: recompute cheaply from VinvX.
        Eigen::MatrixXd XtVinvX = X_design_.transpose() * nm.VinvX();
        Eigen::LDLT<Eigen::MatrixXd> xldlt(XtVinvX);
        if (xldlt.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
        log_det_x = (xldlt.vectorD().array().abs().cwiseMax(1e-300).log()).sum();
        const int n = n_recs_kept_;
        const int p = (int)X_design_.cols();
        double logL = -0.5 * (log_det_v + log_det_x + nm.nullQuad() + (n - p) * std::log(2.0 * M_PI));
        return std::isfinite(logL) ? logL : -std::numeric_limits<double>::infinity();
    } catch (const std::exception&) {
        return -std::numeric_limits<double>::infinity();
    }
}

bool CommonEnvironmentModel::runVCE(const Eigen::MatrixXd& G_used) {
    std::cout << "Running record-level 3-component REML VCE "
              << "(sigma2_a, sigma2_c, sigma2_e)...\n";
    auto t0 = std::chrono::high_resolution_clock::now();

    const bool fixed_mode = (options_.reml_vg >= 0.0 && options_.reml_vc >= 0.0 && options_.reml_ve >= 0.0);
    if (fixed_mode) {
        v_g_ = options_.reml_vg; v_c_ = options_.reml_vc; v_e_ = options_.reml_ve;
        std::cout << "Using user-provided variance components (skipping VCE):\n"
                  << "V(a)=" << v_g_ << " V(c)=" << v_c_ << " V(e)=" << v_e_ << "\n";
        return true;
    }

    // OLS-based initialization
    Eigen::VectorXd beta_ols = X_design_.colPivHouseholderQr().solve(y_obs_);
    Eigen::VectorXd resid = y_obs_ - X_design_ * beta_ols;
    double total_var = resid.squaredNorm() / std::max(1, n_recs_kept_ - (int)X_design_.cols());
    if (!(total_var > 0.0) || !std::isfinite(total_var)) {
        double my = y_obs_.mean();
        total_var = (y_obs_.array() - my).square().sum() / std::max(1, n_recs_kept_ - 1);
    }
    std::cout << "  OLS residual variance = " << total_var << "\n";

    // Moment-based between-group estimate for one of the starts
    double between_group_var = 0.0;
    {
        std::vector<double> g_sum(n_groups_, 0.0), g_cnt(n_groups_, 0.0);
        for (int r = 0; r < n_recs_kept_; ++r) {
            if (obs_to_group_[r] < 0 || obs_to_group_[r] >= n_groups_) continue;
            g_sum[obs_to_group_[r]] += resid[r];
            g_cnt[obs_to_group_[r]] += 1.0;
        }
        double ss_between = 0.0, sum_n2 = 0.0;
        int g_used = 0;
        for (int g = 0; g < n_groups_; ++g) {
            if (g_cnt[g] <= 0.0) continue;
            ++g_used;
            double mg = g_sum[g] / g_cnt[g];
            ss_between += g_cnt[g] * mg * mg;
            sum_n2 += g_cnt[g] * g_cnt[g];
        }
        double kbar = (n_recs_kept_ - sum_n2 / std::max(1.0, (double)n_recs_kept_)) /
                      std::max(1.0, (double)g_used - 1.0);
        double msb = ss_between / std::max(1.0, (double)g_used - 1.0);
        double msw = std::max((resid.squaredNorm() - ss_between) / std::max(1.0, (double)n_recs_kept_ - (double)g_used), 1e-10);
        between_group_var = std::max(0.0, (msb - msw) / std::max(1e-12, kbar));
    }

    std::vector<std::pair<double, double>> inits = {
        {0.30, 0.20}, {0.20, 0.30}, {0.40, 0.10},
        {0.10, 0.20}, {0.15, 0.35}, {0.05, 0.10},
    };
    if (between_group_var > 1e-8) {
        double vg0 = std::max(total_var * 0.3, 1e-4);
        double vc0 = std::min(2.0 * between_group_var, total_var * 0.5);
        inits.push_back({vg0, vc0});
    }

    double best_logL = -std::numeric_limits<double>::infinity();
    double best_vg = total_var * 0.3, best_vc = total_var * 0.15, best_ve = total_var * 0.55;
    int best_status = 999;

    for (size_t ti = 0; ti < inits.size(); ++ti) {
        double lv[3] = {
            std::log(std::max(total_var * inits[ti].first, 1e-5)),
            std::log(std::max(total_var * inits[ti].second, 1e-5)),
            std::log(std::max(total_var * (1.0 - inits[ti].first - inits[ti].second), 1e-5)),
        };
        double vg = std::exp(lv[0]), vc = std::exp(lv[1]), ve = std::exp(lv[2]);
        double cur = recordLogL(G_used, vg, vc, ve);
        if (!std::isfinite(cur)) continue;
        double step = 0.4;
        int iter = 0;
        int max_iter = std::max(20, options_.reml_maxit);
        for (; iter < max_iter && step > 1e-5; ++iter) {
            bool improved = false;
            for (int d = 0; d < 3; ++d) {
                for (double dir : {-1.0, 1.0}) {
                    double cand[3] = {lv[0], lv[1], lv[2]};
                    cand[d] += dir * step;
                    double logL = recordLogL(G_used, std::exp(cand[0]), std::exp(cand[1]), std::exp(cand[2]));
                    if (std::isfinite(logL) && logL > cur + 1e-8) {
                        lv[0] = cand[0]; lv[1] = cand[1]; lv[2] = cand[2];
                        cur = logL;
                        improved = true;
                    }
                }
            }
            if (!improved) step *= 0.5;
        }
        vg = std::exp(lv[0]); vc = std::exp(lv[1]); ve = std::exp(lv[2]);
        double total = vg + vc + ve;
        int status = 0;
        if (!std::isfinite(cur)) status = 3;
        else if (vg / total < 1e-4) status = 1;
        else if (ve / total < 0.05) status = 2;
        std::cout << "  init #" << ti << " (vg_f=" << inits[ti].first
                  << ", vc_f=" << inits[ti].second << "): V(a)=" << vg
                  << " V(c)=" << vc << " V(e)=" << ve
                  << " logL=" << cur << " status=" << status
                  << " iters=" << (iter + 1) << "\n";
        if (status >= 0 && std::isfinite(cur) &&
            (status < best_status || (status == best_status && cur > best_logL))) {
            best_logL = cur; best_status = status;
            best_vg = vg; best_vc = vc; best_ve = ve;
        }
    }

    if (!std::isfinite(best_logL)) {
        std::cerr << "ERROR: All common-environment VCE starts failed.\n";
        best_vg = std::max(total_var * 0.3, 1e-4);
        best_vc = std::max(between_group_var, 1e-4);
        best_ve = std::max(total_var * 0.55, 1e-4);
        best_status = 3;
    }

    v_g_ = best_vg; v_c_ = best_vc; v_e_ = best_ve;
    vce_converged_ = (best_status <= 1);
    vce_iters_ = 0; // reported per-init; global flag only
    vce_best_logL_ = best_logL;

    double total = v_g_ + v_c_ + v_e_;
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "VCE done in " << std::chrono::duration<double>(t1 - t0).count() << "s\n";
    std::cout << "V(a)  = " << v_g_ << "\n"
              << "V(c)  = " << v_c_ << "\n"
              << "V(e)  = " << v_e_ << "\n"
              << "V(p)  = " << total << "\n"
              << "h2    = " << v_g_ / total << "\n"
              << "ce_fraction = " << v_c_ / total
              << " (common-environment fraction)\n";
    return std::isfinite(best_logL);
}

// ---------------------------------------------------------------------------
// Fixed-null SNP scan
// ---------------------------------------------------------------------------

void CommonEnvironmentModel::scanWithGrm(const Eigen::MatrixXd& G_used,
                                         const std::string& assoc_file,
                                         const std::string& chrom_filter) {
    std::vector<RandomEffectSpec> specs = {additiveSpec(), commonEnvSpec()};
    RecordLevelNullModel nm;
    if (!nm.prepare(y_obs_, X_design_, specs, {v_g_, v_c_}, v_e_, &G_used)) {
        throw std::runtime_error("Common-environment null model factorization failed at scan stage.");
    }
    std::cout << "  Scan: fixed-null record-level GLS (V fixed, n=" << n_recs_kept_
              << ", p=" << X_design_.cols() << ")\n";

    struct VariantMeta {
        std::string chrom, snp_id, a1, a2;
        long long pos;
    };

    const int block_size = 256;
    Eigen::MatrixXd raw_block(n_recs_kept_, block_size);
    std::vector<VariantMeta> metas;
    metas.reserve(block_size);
    int filled = 0;

    auto flush_block = [&]() {
        if (filled == 0) return;
        const int bsz = filled;
        Eigen::MatrixXd G_block(n_recs_kept_, bsz);
        std::vector<double> freqs(bsz, 0.0);
        std::vector<int> valid_ns(bsz, 0);
        for (int c = 0; c < bsz; ++c) {
            double sum_g = 0.0;
            int valid_n = 0;
            long long h1 = 0, ht = 0, h2 = 0;
            for (int r = 0; r < n_recs_kept_; ++r) {
                float val = (float)raw_block(r, c);
                G_block(r, c) = val;
                if (!std::isnan(val)) {
                    sum_g += val; ++valid_n;
                    int gt = (int)std::llround(val);
                    if (gt == 0) ++h2; else if (gt == 1) ++ht; else ++h1;
                }
            }
            freqs[c] = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
            valid_ns[c] = valid_n;
            double maf = std::min(freqs[c], 1.0 - freqs[c]);
            bool pass = valid_n > 0 && maf >= options_.min_maf;
            if (pass && options_.min_hwe > 0.0) {
                pass = calculateHWE(h1, ht, h2) >= options_.min_hwe;
            }
            double mean = 2.0 * freqs[c];
            for (int r = 0; r < n_recs_kept_; ++r) {
                if (!pass) G_block(r, c) = 0.0;
                else if (std::isnan(G_block(r, c))) G_block(r, c) = 0.0;
                else G_block(r, c) -= mean;
            }
            if (!pass) { freqs[c] = 0.0; valid_ns[c] = 0; }
        }

        Eigen::MatrixXd VinvG = nm.applyVinv(G_block);
        Eigen::VectorXd scores = G_block.transpose() * nm.Py();
        Eigen::MatrixXd xvX = G_block.transpose() * nm.VinvX();
        Eigen::MatrixXd xvX_inv = xvX * nm.XtVinvXInv();
        Eigen::VectorXd term2_all = xvX.cwiseProduct(xvX_inv).rowwise().sum();
        Eigen::VectorXd gVg = G_block.cwiseProduct(VinvG).colwise().sum();

        for (int c = 0; c < bsz; ++c) {
            if (valid_ns[c] <= 0) continue;
            const double var_score = gVg(c) - term2_all(c);
            GwasResult res;
            res.chrom = metas[c].chrom; res.snp_id = metas[c].snp_id; res.pos = metas[c].pos;
            res.a1 = metas[c].a1; res.a2 = metas[c].a2;
            res.freq = freqs[c];
            res.n_obs = valid_ns[c];
            if (var_score > 1e-12 && std::isfinite(scores(c))) {
                res.beta = scores(c) / var_score;
                res.se = std::sqrt(1.0 / var_score);
                const double stat = res.beta / res.se;
                res.p_value = std::erfc(std::abs(stat) / std::sqrt(2.0));
                res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
            }
            appendAssocResult(assoc_file, res, false);
        }
        metas.clear();
        filled = 0;
    };

    auto add_variant = [&](const Eigen::VectorXf& dosage, const std::string& chrom,
                           const std::string& snp_id, long long pos,
                           const std::string& a1, const std::string& a2) {
        if (!chrom_filter.empty() && chrom != chrom_filter) return;
        if (filled == block_size) flush_block();
        for (int r = 0; r < n_recs_kept_; ++r) {
            float v = dosage[gen_records_[r].aid - 1];
            raw_block(r, filled) = std::isnan(v)
                ? std::numeric_limits<double>::quiet_NaN() : (double)v;
        }
        metas.push_back({chrom, snp_id, a1, a2, pos});
        ++filled;
    };

    if (is_bgen_) {
        bgen_reader_->reset();
        std::string chrom, rsid, a1, a2;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        while (bgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage)) {
            add_variant(dosage, chrom, rsid, (long long)pos, a1, a2);
        }
        flush_block();
        return;
    }
    if (is_pgen_) {
        pgen_reader_->reset();
        std::string chrom, rsid, a1, a2;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        while (pgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage)) {
            add_variant(dosage, chrom, rsid, (long long)pos, a1, a2);
        }
        flush_block();
        return;
    }

    // PLINK path: random-access, parallel genotype fetch.
    const auto& bim = plink_reader_->getBimInfo();
    const int m = (int)bim.size();
    for (int start = 0; start < m; start += block_size) {
        int end = std::min(m, start + block_size);
        std::vector<int> snp_indices;
        snp_indices.reserve(end - start);
        for (int j = start; j < end; ++j) {
            if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
            snp_indices.push_back(j);
        }
        if (snp_indices.empty()) continue;
        int bsz = (int)snp_indices.size();
        if (filled > 0) flush_block();
        std::vector<std::vector<float>> tmp(bsz);
#pragma omp parallel for schedule(dynamic, 8)
        for (int c = 0; c < bsz; ++c) {
            Eigen::VectorXf g_raw;
            plink_reader_->getSnpGenotypes(snp_indices[c], g_raw);
            tmp[c].assign(g_raw.data(), g_raw.data() + g_raw.size());
        }
        for (int c = 0; c < bsz; ++c) {
            const auto& snp = bim[snp_indices[c]];
            add_variant(Eigen::Map<const Eigen::VectorXf>(tmp[c].data(), (int)tmp[c].size()),
                        snp.chrom, snp.id, snp.bp_pos, snp.alt, snp.ref);
        }
    }
    flush_block();
}

// ---------------------------------------------------------------------------
// run / helpers
// ---------------------------------------------------------------------------

void CommonEnvironmentModel::run() {
    std::cout << "\n--- Common-Environment GWAS Model ---\n";
    std::cout << "Model: y = Xb + x_j beta_j + Z_a a + Z_c c + e\n"
              << "       a ~ N(0, K sigma2_a), c ~ N(0, I_g sigma2_c)\n"
              << "  (" << n_groups_ << " environment group levels)\n";

    // Build obs_to_sample_ (record -> subject rank over kept subjects)
    std::set<int> subject_set;
    for (int r = 0; r < n_recs_kept_; ++r) subject_set.insert(gen_records_[r].aid - 1);
    n_subjects_kept_ = (int)subject_set.size();
    std::vector<int> subject_rank(n_samples_, -1);
    {
        int rank = 0;
        for (int s : subject_set) subject_rank[s] = rank++;
    }
    obs_to_sample_.resize(n_recs_kept_);
    for (int r = 0; r < n_recs_kept_; ++r) obs_to_sample_[r] = subject_rank[gen_records_[r].aid - 1];

    std::cout << "Data: " << n_recs_kept_ << " records, " << n_subjects_kept_
              << " subjects, " << n_groups_ << " groups, "
              << X_design_.cols() << " fixed-effect columns (incl. intercept).\n";

    if (is_bgen_ || is_pgen_) {
        std::cout << "Loading variant info for chromosome labels...\n";
    }
    std::set<std::string> chroms_set;
    if (is_bgen_) {
        bgen_reader_->reset();
        bgen_reader_->loadVariantInfo();
        for (const auto& v : bgen_reader_->getVariants()) chroms_set.insert(v.chrom);
    } else if (is_pgen_) {
        for (const auto& v : pgen_reader_->getVariants()) chroms_set.insert(v.chrom);
    } else {
        for (const auto& snp : plink_reader_->getBimInfo()) chroms_set.insert(snp.chrom);
    }
    bool use_loco = chroms_set.size() > 1 && !options_.no_loco;

    std::string assoc_file = options_.out_prefix + ".assoc.txt";
    {
        std::ofstream ofs(assoc_file);
        writeAssocHeader(assoc_file, false);
    }

    const std::string null_file = options_.out_prefix + ".common_env.null.txt";
    Eigen::MatrixXd G_used = subsetByRank(G_subject_, subject_rank);
    if (!use_loco) {
        runVCE(G_used);
        scanWithGrm(G_used, assoc_file, "");
        writeNullSummary(null_file, "record_level_gls");
    } else {
        GrmBuilder::Options grm_opts;
        grm_opts.min_maf = (float)options_.min_maf;
        grm_opts.min_hwe = (float)options_.min_hwe;
        grm_opts.algorithm = options_.grm_alg;
#ifdef _OPENMP
        grm_opts.thread_num = omp_get_max_threads();
#endif
        GrmBuilder split_builder(grm_opts);
        GrmBuilder::SplitGrmResult split;
        if (is_bgen_) { bgen_reader_->reset(); split = split_builder.computeSplit(*bgen_reader_, chroms_set); }
        else if (is_pgen_) { pgen_reader_->reset(); split = split_builder.computeSplit(*pgen_reader_, chroms_set); }
        else { split = split_builder.computeSplit(*plink_reader_, chroms_set); }

        double S_global = (options_.grm_alg == 0) ? (double)split.global_snp_count
                                                  : split.global_sum_2pq;
        std::string last_used_mode = "record_level_gls_loco";
        for (size_t ci = 0; ci < split.chrom_labels.size(); ++ci) {
            const std::string& chrom = split.chrom_labels[ci];
            std::cout << "  LOCO chr " << chrom << "...\n";
            double S_chr = (options_.grm_alg == 0)
                ? (double)split.per_chrom_snp_count[ci]
                : split.per_chrom_sum_2pq[ci];
            double S_loco = S_global - S_chr;
            Eigen::MatrixXd G_chr = subsetByRank(split.per_chrom_raw[ci], subject_rank);
            if (S_chr > 0.0 && S_loco > 1e-6) {
                Eigen::MatrixXd G_loco = G_used * (S_global / S_loco) - G_chr / S_loco;
                G_loco = 0.5 * (G_loco + G_loco.transpose()).eval();
                runVCE(G_loco);
                scanWithGrm(G_loco, assoc_file, chrom);
            } else {
                runVCE(G_used);
                scanWithGrm(G_used, assoc_file, chrom);
            }
            writeNullSummary(null_file, last_used_mode);
        }
    }
    std::cout << "Common-environment GWAS results saved to " << assoc_file << "\n";
}

void CommonEnvironmentModel::writeNullSummary(const std::string& filename,
                                              const std::string& scan_mode) const {
    std::ofstream ofs(filename);
    double total = v_g_ + v_c_ + v_e_;
    ofs << "Item\tValue\n"
        << "n_records\t" << n_recs_kept_ << "\n"
        << "n_subjects\t" << n_subjects_kept_ << "\n"
        << "n_groups\t" << n_groups_ << "\n"
        << "sigma2_a\t" << v_g_ << "\n"
        << "sigma2_common\t" << v_c_ << "\n"
        << "sigma2_e\t" << v_e_ << "\n"
        << "h2\t" << ((total > 0.0) ? (v_g_ / total) : 0.0) << "\n"
        << "common_environment_fraction\t" << ((total > 0.0) ? (v_c_ / total) : 0.0) << "\n"
        << "converged\t" << (vce_converged_ ? "yes" : "no") << "\n"
        << "vce_best_logL\t" << vce_best_logL_ << "\n"
        << "solver\trecord_level_dense_ldlt" << "\n"
        << "scan_mode\t" << scan_mode << "\n";
    ofs.close();
    std::cout << "Common-environment null summary written to " << filename << "\n";
}

// ---------------------------------------------------------------------------
// Small private helpers (anonymous namespace, top of file)
// ---------------------------------------------------------------------------

} // namespace cosmic