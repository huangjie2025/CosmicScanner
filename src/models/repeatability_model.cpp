#include "repeatability_model.h"
#include "repeatability_helpers.h"
#include "vce.h"
#include "stcg_vce.h"
#include "logger.h"
#include "matrix_io.h"
#include "matrix_free_pcg.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <sstream>
#include <sys/resource.h>
#include <omp.h>
#include "grm_builder.h"
#include "pheno_reader.h"
#include "covar_reader.h"
#include "data_filter.h"
#include "qc_utils.h"
#include "vce.h"
#include "eigen_vce.h"
#include "lmm_engine.h"

extern "C" {
    void openblas_set_num_threads(int num_threads) __attribute__((weak));
}

namespace cosmic {

static EigenVce::NullSolver parse_repeatability_lmm_null_solver(const std::string& value) {
    if (value == "brent") return EigenVce::NullSolver::Brent;
    if (value == "idul") return EigenVce::NullSolver::IDUL;
    if (value == "idul_plus" || value == "idul+") return EigenVce::NullSolver::IDULPlus;
    if (value == "aidul") return EigenVce::NullSolver::AIDUL;
    throw std::runtime_error("Unknown --lmm-null-solver: " + value + " (expected brent, idul, idul_plus, or aidul)");
}

static std::unordered_set<std::string> load_repeatability_snp_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (ifs >> id) ids.insert(id);
    if (ids.empty()) throw std::runtime_error("SNP list is empty: " + path);
    return ids;
}

RepeatabilityModel::Options RepeatabilityModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.out_ssf = cfg.out_ssf;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.grm_alg = cfg.grm_alg;
    opts.n_threads = 1;
    #ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
    #endif
    opts.no_loco = cfg.no_loco;
    opts.fast_lmm = cfg.fast_lmm;
    opts.scan_tier = cfg.scan_tier;
    opts.eb_lmm = cfg.eb_lmm;
    opts.mixed_precision = cfg.repeatability_mixed_precision;
    opts.refit_validation = cfg.repeatability_refit_validation;
    if (!opts.scan_tier.empty() && cfg.repeatability_refit_validation && opts.scan_tier != "exact") {
        throw std::runtime_error("--repeatability-refit-validation conflicts with --scan-tier " + opts.scan_tier);
    }
    if (opts.scan_tier == "exact") opts.refit_validation = true;
    opts.refit_limit = cfg.repeatability_refit_limit;
    opts.grm_cache_prefix = cfg.repeatability_grm_cache_prefix;
    opts.weight_file = cfg.repeatability_weight_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.extract_set = load_repeatability_snp_id_set(cfg.extract_file);
    opts.exclude_set = load_repeatability_snp_id_set(cfg.exclude_file);
    opts.use_covariate_projection = cfg.use_covariate_projection;
    opts.lmm_null_solver = parse_repeatability_lmm_null_solver(cfg.lmm_null_solver);
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.reml_vg = cfg.reml_vg;
    opts.reml_vpe = cfg.reml_vpe;
    opts.reml_ve = cfg.reml_ve;
    opts.pheno_file = cfg.pheno_file;
    opts.mpheno = cfg.mpheno;
    opts.tpheno = cfg.tpheno_set ? cfg.tpheno : 0;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.use_spa = cfg.use_spa;
    // Thread count override
    if (cfg.num_threads > 0) {
        opts.n_threads = cfg.num_threads;
        omp_set_num_threads(cfg.num_threads);
    }
    // VCE refinement mode
    std::string vrm = cfg.vce_refine_mode;
    std::transform(vrm.begin(), vrm.end(), vrm.begin(), ::tolower);
    if (vrm == "full-record-exact") opts.vce_refine_mode = VceRefineMode::FullRecordExact;
    else if (vrm == "subject-level-exact") opts.vce_refine_mode = VceRefineMode::SubjectLevelExact;
    else if (vrm == "fast") opts.vce_refine_mode = VceRefineMode::Fast;
    else opts.vce_refine_mode = VceRefineMode::Auto;
    // IASBLUP-style repeatability mode
    std::string rm = cfg.repeatability_mode;
    std::transform(rm.begin(), rm.end(), rm.begin(), ::tolower);
    if (rm == "exact") opts.repeatability_mode = RepeatabilityMode::Exact;
    else if (rm == "fast") opts.repeatability_mode = RepeatabilityMode::Fast;
    else if (rm == "stcg") opts.repeatability_mode = RepeatabilityMode::STCG;
    else opts.repeatability_mode = RepeatabilityMode::Auto;
    if (!opts.scan_tier.empty() && opts.repeatability_mode != RepeatabilityMode::Exact) {
        throw std::runtime_error("--scan-tier for repeatability GWAS requires --repeatability-mode exact");
    }
    return opts;
}

namespace {

static std::string dbg_num(double x) {
    if (std::isnan(x)) return "null";
    if (!std::isfinite(x)) return (x > 0.0) ? "1e999" : "-1e999";
    std::ostringstream os;
    os << std::setprecision(15) << x;
    return os.str();
}

template <typename Derived>
static int dbg_nonfinite_count(const Eigen::MatrixBase<Derived>& m) {
    int count = 0;
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (!std::isfinite((double)m(i, j))) ++count;
        }
    }
    return count;
}

template <typename Derived>
static std::pair<double, double> dbg_diag_minmax(const Eigen::MatrixBase<Derived>& m) {
    if (m.rows() == 0 || m.cols() == 0) return {std::numeric_limits<double>::quiet_NaN(),
                                                std::numeric_limits<double>::quiet_NaN()};
    double mn = std::numeric_limits<double>::infinity();
    double mx = -std::numeric_limits<double>::infinity();
    int d = std::min<int>(m.rows(), m.cols());
    for (int i = 0; i < d; ++i) {
        double v = (double)m(i, i);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
    }
    return {mn, mx};
}

static int drop_aliased_fixed_effects(Eigen::MatrixXf& X, FixedDesignG& fd) {
    if (X.cols() <= 1) return 0;
    Eigen::MatrixXd Xd = X.cast<double>();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Xd);
    qr.setThreshold(1e-10);
    const int rank = qr.rank();
    if (rank >= X.cols()) return 0;

    Eigen::VectorXi piv = qr.colsPermutation().indices();
    std::vector<int> keep_cols;
    keep_cols.reserve(rank);
    for (int i = 0; i < rank; ++i) keep_cols.push_back(piv[i]);
    std::sort(keep_cols.begin(), keep_cols.end());

    std::vector<int> old_to_new(X.cols(), -1);
    for (int i = 0; i < (int)keep_cols.size(); ++i) old_to_new[keep_cols[i]] = i;

    Eigen::MatrixXf X_new(X.rows(), rank);
    for (int j = 0; j < rank; ++j) X_new.col(j) = X.col(keep_cols[j]);
    X = std::move(X_new);

    for (auto& row : fd.rows) {
        std::vector<std::pair<int, double>> new_row;
        new_row.reserve(row.size());
        for (const auto& kv : row) {
            if (kv.first < 0 || kv.first >= (int)old_to_new.size()) continue;
            int new_idx = old_to_new[kv.first];
            if (new_idx >= 0) new_row.push_back({new_idx, kv.second});
        }
        row.swap(new_row);
    }

    if (!fd.names.empty() && (int)fd.names.size() == (int)old_to_new.size()) {
        std::vector<std::string> new_names;
        new_names.reserve(rank);
        for (int idx : keep_cols) new_names.push_back(fd.names[idx]);
        fd.names.swap(new_names);
    }
    fd.p = rank;
    return (int)old_to_new.size() - rank;
}

static void set_blas_threads(int n) {
    if (openblas_set_num_threads) openblas_set_num_threads(n);
}

static double gamma_q(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0; del *= x / ap; sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7)
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
        }
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1e-30, d = 1.0 / b, h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a); b += 2.0;
            d = an * d + b; if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c; if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d; double del = d * c; h *= del;
            if (std::abs(del - 1.0) < 1e-7)
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
        }
    }
    return 0.0;
}

static double chisq_prob(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return gamma_q(df / 2.0, stat / 2.0);
}

// Algorithm 0 (Yang/GCTA default), 3 (Yang et al), 5 (TrueVar) apply per-SNP
// standardization (divide by sqrt(2pq) or sqrt(var)), so each SNP contributes
// equally to the GRM trace. LOCO scaling uses SNP count for these.
// Algorithms 1 (VanRaden), 2 (Zeng), 4 (Vitezica) only center, so LOCO
// scaling uses sum_2pq to account for heterogeneous per-SNP variance.

} // namespace

RepeatabilityModel::RepeatabilityModel(const Options& opts) : options_(opts) {}

void RepeatabilityModel::prepare(const std::vector<std::string>& fids,
                                 const std::vector<std::string>& iids,
                                 int n_samples,
                                 PlinkReader& plink_reader,
                                 BgenReader& bgen_reader,
                                 PgenReader& pgen_reader,
                                 bool is_bgen, bool is_pgen) {
    fids_ = fids; iids_ = iids; n_samples_ = n_samples;
    is_bgen_ = is_bgen; is_pgen_ = is_pgen;
    plink_reader_ = &plink_reader; bgen_reader_ = &bgen_reader; pgen_reader_ = &pgen_reader;

    loadPhenotype();
    loadCovariates();
    auto grm_start = std::chrono::high_resolution_clock::now();
    buildGRM();
    auto grm_end = std::chrono::high_resolution_clock::now();
    grm_time_sec_ = std::chrono::duration<double>(grm_end - grm_start).count();
    buildMMEComponents();
    auto vce_start = std::chrono::high_resolution_clock::now();
    runVCE();
    auto vce_end = std::chrono::high_resolution_clock::now();
    vce_time_sec_ = std::chrono::duration<double>(vce_end - vce_start).count();
    std::cout << "Repeatability prepare timings: GRM=" << std::fixed << std::setprecision(3)
              << grm_time_sec_ << "s, VCE=" << vce_time_sec_ << "s\n";
}

void RepeatabilityModel::loadPhenotype() {
    if (options_.pheno_file.empty()) throw std::runtime_error("Repeatability model requires --pheno file (long format).");
    std::cout << "Loading long-format phenotype data...\n";
    PhenoReader pheno_reader;
    auto records = (options_.tpheno > 0)
        ? pheno_reader.loadLongitudinalRecords(options_.pheno_file, options_.tpheno, options_.mpheno)
        : pheno_reader.loadAllRecords(options_.pheno_file, options_.mpheno);

    std::unordered_map<std::string, int> fid_iid_to_gen_idx;
    for (int i = 0; i < n_samples_; ++i) fid_iid_to_gen_idx[fids_[i] + ":" + iids_[i]] = i;

    gen_records_.clear();
    record_weights_.resize(0);
    std::unordered_map<std::string, std::deque<double>> weight_by_id;
    std::unordered_map<std::string, double> weight_by_id_time;
    bool has_weight_file = !options_.weight_file.empty();
    if (has_weight_file) {
        std::ifstream wfs(options_.weight_file);
        if (!wfs) throw std::runtime_error("Cannot open repeatability weight file: " + options_.weight_file);
        std::string line;
        int line_no = 0;
        while (std::getline(wfs, line)) {
            ++line_no;
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            std::vector<std::string> tok;
            std::string t;
            while (iss >> t) tok.push_back(t);
            if (tok.size() < 3) continue;
            try {
                std::string id_key = tok[0] + ":" + tok[1];
                if (tok.size() >= 4) {
                    double time = std::stod(tok[2]);
                    double weight = std::stod(tok[3]);
                    if (weight > 0.0 && std::isfinite(weight)) {
                        weight_by_id_time[id_key + ":" + std::to_string(time)] = weight;
                    }
                } else {
                    double weight = std::stod(tok[2]);
                    if (weight > 0.0 && std::isfinite(weight)) weight_by_id[id_key].push_back(weight);
                }
            } catch (...) {
                if (line_no == 1) continue;
                throw std::runtime_error("Invalid repeatability weight line " + std::to_string(line_no) + ": " + line);
            }
        }
    }
    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = fid_iid_to_gen_idx.find(key);
        if (it != fid_iid_to_gen_idx.end()) {
            GenRecord gr;
            gr.aid = it->second + 1;
            gr.y = rec.value;
            gr.time = rec.time;
            gr.idstr = key;
            gen_records_.push_back(gr);
            double weight = 1.0;
            if (has_weight_file) {
                bool found = false;
                if (!std::isnan(rec.time)) {
                    auto wt = weight_by_id_time.find(key + ":" + std::to_string(rec.time));
                    if (wt != weight_by_id_time.end()) { weight = wt->second; found = true; }
                }
                if (!found) {
                    auto wi = weight_by_id.find(key);
                    if (wi != weight_by_id.end() && !wi->second.empty()) {
                        weight = wi->second.front();
                        wi->second.pop_front();
                        found = true;
                    }
                }
                if (!found) throw std::runtime_error("Missing repeatability weight for record: " + key);
            }
            record_weights_.conservativeResize(record_weights_.size() + 1);
            record_weights_[record_weights_.size() - 1] = weight;
        }
    }
    n_recs_kept_ = (int)gen_records_.size();
    std::cout << "Kept " << n_recs_kept_ << " phenotype records matching genotype samples.\n";
    if (record_weights_.size() != n_recs_kept_) record_weights_ = Eigen::VectorXd::Ones(n_recs_kept_);
    if (has_weight_file) std::cout << "Loaded record-specific residual weights for repeatability model.\n";
    if (n_recs_kept_ < 10) throw std::runtime_error("Too few records.");
}

void RepeatabilityModel::loadCovariates() {
    Eigen::MatrixXf covs_q, covs_c;
    bool has_q = false, has_c = false;
    bool has_record_time = false;
    for (const auto& rec : gen_records_) {
        if (!std::isnan(rec.time)) { has_record_time = true; break; }
    }

    std::vector<std::string> rec_fids, rec_iids;
    std::vector<double> rec_times;
    if (has_record_time) {
        rec_fids.reserve(gen_records_.size());
        rec_iids.reserve(gen_records_.size());
        rec_times.reserve(gen_records_.size());
        for (const auto& rec : gen_records_) {
            size_t pos = rec.idstr.find(':');
            rec_fids.push_back(pos == std::string::npos ? "" : rec.idstr.substr(0, pos));
            rec_iids.push_back(pos == std::string::npos ? rec.idstr : rec.idstr.substr(pos + 1));
            rec_times.push_back(rec.time);
        }
    }

    if (!options_.qcovar_file.empty()) {
        CovarReader cov_reader;
        if (has_record_time) {
            cov_reader.loadWithTime(options_.qcovar_file, 1);
            cov_reader.alignRecords(rec_fids, rec_iids, rec_times, covs_q);
        } else {
            cov_reader.load(options_.qcovar_file);
            cov_reader.align(fids_, iids_, covs_q);
        }
        has_q = true;
    }
    if (!options_.covar_file.empty()) {
        CovarReader cov_reader;
        if (has_record_time) {
            cov_reader.loadCategoricalWithTime(options_.covar_file, 1);
            cov_reader.alignRecords(rec_fids, rec_iids, rec_times, covs_c);
        } else {
            cov_reader.loadCategorical(options_.covar_file);
            cov_reader.align(fids_, iids_, covs_c);
        }
        has_c = true;
    }

    int n_cov_cols = 1;
    if (has_q) n_cov_cols += covs_q.cols();
    if (has_c) n_cov_cols += covs_c.cols();

    fd_.p = n_cov_cols;
    fd_.rows.resize(n_recs_kept_);
    for (int i = 0; i < n_recs_kept_; ++i) {
        int gen_idx = gen_records_[i].aid - 1;
        fd_.rows[i].push_back({0, 1.0});
        int current_col = 1;
        if (has_q) {
            for (int k = 0; k < covs_q.cols(); ++k) {
                float val = has_record_time ? covs_q(i, k) : covs_q(gen_idx, k);
                if (!std::isnan(val)) fd_.rows[i].push_back({current_col + k, (double)val});
            }
            current_col += covs_q.cols();
        }
        if (has_c) {
            for (int k = 0; k < covs_c.cols(); ++k) {
                float val = has_record_time ? covs_c(i, k) : covs_c(gen_idx, k);
                if (!std::isnan(val)) fd_.rows[i].push_back({current_col + k, (double)val});
            }
        }
    }

    X_obs_.resize(n_recs_kept_, fd_.p);
    y_obs_.resize(n_recs_kept_);
    obs_to_sample_.resize(n_recs_kept_);
    for (int i = 0; i < n_recs_kept_; ++i) {
        y_obs_(i) = gen_records_[i].y;
        for (int j = 0; j < fd_.p; ++j) X_obs_(i, j) = 0.0f;
        for (const auto& kv : fd_.rows[i]) {
            if (kv.first >= 0 && kv.first < fd_.p) X_obs_(i, kv.first) = (float)kv.second;
        }
        obs_to_sample_[i] = gen_records_[i].aid - 1;
    }

    int dropped_aliased_cols = drop_aliased_fixed_effects(X_obs_, fd_);
    if (dropped_aliased_cols > 0) {
        LOG_WARN("Repeatability fixed-effect design matrix is rank-deficient; dropped "
                 << dropped_aliased_cols << " aliased columns before VCE/scan");
    }

    // #region debug-point A:design-matrix
    {
        Eigen::MatrixXd X_dbg = X_obs_.cast<double>();
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_dbg(X_dbg);
        int zero_cols = 0;
        int intercept_only_rows = 0;
        for (int j = 0; j < X_dbg.cols(); ++j) {
            bool any_nonzero = false;
            for (int i = 0; i < X_dbg.rows(); ++i) {
                if (X_dbg(i, j) != 0.0f) {
                    any_nonzero = true;
                    break;
                }
            }
            if (!any_nonzero) ++zero_cols;
        }
        for (int i = 0; i < X_dbg.rows(); ++i) {
            bool only_intercept = (X_dbg.cols() > 0 && X_dbg(i, 0) == 1.0);
            for (int j = 1; j < X_dbg.cols(); ++j) {
                if (X_dbg(i, j) != 0.0f) {
                    only_intercept = false;
                    break;
                }
            }
            if (only_intercept) ++intercept_only_rows;
        }
        int nonfinite_y = 0;
        for (int i = 0; i < y_obs_.size(); ++i) if (!std::isfinite(y_obs_[i])) ++nonfinite_y;
        std::ostringstream data;
        data << "{"
             << "\"n_records\":" << n_recs_kept_
             << ",\"n_subjects\":" << n_samples_
             << ",\"x_cols\":" << fd_.p
             << ",\"qcov_cols\":" << (has_q ? covs_q.cols() : 0)
             << ",\"ccov_cols\":" << (has_c ? covs_c.cols() : 0)
             << ",\"qcov_rows\":" << (has_q ? covs_q.rows() : 0)
             << ",\"ccov_rows\":" << (has_c ? covs_c.rows() : 0)
             << ",\"x_rank\":" << qr_dbg.rank()
             << ",\"aliased_cols_dropped\":" << dropped_aliased_cols
             << ",\"zero_cols\":" << zero_cols
             << ",\"intercept_only_rows\":" << intercept_only_rows
             << ",\"x_nonfinite\":" << dbg_nonfinite_count(X_dbg)
             << ",\"y_nonfinite\":" << nonfinite_y
             << "}";
        std::cerr << "Repeatability design matrix diagnostics: "
                  << data.str() << std::endl;
    }
}

void RepeatabilityModel::buildGRM() {
    std::cout << "Calculating GRM...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    #ifdef _OPENMP
    grm_opts.thread_num = omp_get_max_threads();
    #endif

    {
        GrmBuilder builder(grm_opts);
        if (is_pgen_) { pgen_reader_->reset(); builder.compute(*pgen_reader_); }
        else if (is_bgen_) { bgen_reader_->reset(); builder.compute(*bgen_reader_); }
        else { builder.compute(*plink_reader_); }

        Eigen::MatrixXd G = builder.getGrm().cast<double>();
        G_subject_ = G;
        const bool needs_grm_inverse =
            options_.repeatability_mode == RepeatabilityMode::STCG ||
            (options_.repeatability_mode == RepeatabilityMode::Auto && n_recs_kept_ > 80000) ||
            options_.vce_refine_mode == VceRefineMode::FullRecordExact ||
            options_.refit_validation;
        if (needs_grm_inverse) {
            std::cout << "Inverting GRM for full-record/STCG/refit path...\n";
            G.diagonal().array() += 0.001;
            Eigen::LLT<Eigen::MatrixXd> lltG(G);
            if (lltG.info() == Eigen::Success) {
                Ginv_ = lltG.solve(Eigen::MatrixXd::Identity(G.rows(), G.cols()));
            } else {
                Ginv_ = G.inverse();
            }
        } else {
            Ginv_.resize(0, 0);
        }
    }
    set_blas_threads(options_.n_threads);
}

void RepeatabilityModel::buildMMEComponents() {
    if (Ginv_.size() == 0) {
        components_.clear();
        return;
    }
    Ginv_adapter_ = std::make_unique<DenseMatrixAdapter>(Ginv_);
    RandomComponent comp_u;
    comp_u.Qinv = Ginv_adapter_.get();
    comp_u.id_map.resize(n_recs_kept_);
    for (int i = 0; i < n_recs_kept_; ++i) comp_u.id_map[i] = gen_records_[i].aid - 1;

    I_mat_.resize(n_samples_, n_samples_); I_mat_.setIdentity();
    I_adapter_ = std::make_unique<SparseMatrixAdapter>(I_mat_);
    RandomComponent comp_pe;
    comp_pe.Qinv = I_adapter_.get();
    comp_pe.id_map = comp_u.id_map;

    components_ = {comp_u, comp_pe};

    // For STCG mode, attach genotype matrix operator for matrix-free VCE
    if (options_.repeatability_mode == RepeatabilityMode::STCG ||
        (options_.repeatability_mode == RepeatabilityMode::Auto && n_recs_kept_ > 80000)) {
        if (!is_bgen_ && !is_pgen_ && plink_reader_) {
            // PlinkReader inherits GenotypeMatrix — direct attachment
            components_[0].geno_mat = plink_reader_;
            LOG_INFO("STCG mode: attached PlinkReader as GenotypeMatrix operator");
        }
        // BGEN/PGEN don't implement GenotypeMatrix yet; STCG_VCE will fall back
        // to GRM-based AI_REML in runVCE()
    }
}

void RepeatabilityModel::runVCE() {
    bool skip_vce = (options_.reml_vg >= 0.0 && options_.reml_vpe >= 0.0 && options_.reml_ve >= 0.0);
    if (skip_vce) {
        v_g_ = options_.reml_vg; v_pe_ = options_.reml_vpe; v_e_ = options_.reml_ve;
        std::cout << "Using user-provided variance components (skip VCE):\n";
        std::cout << "V(G)  = " << v_g_ << "\nV(pe) = " << v_pe_ << "\nV(e)  = " << v_e_ << "\n";
        return;
    }

    std::cout << "Setting up VCE (Repeatability Model)...\n";

    // Use matrix-free STCG VCE for large PLINK repeatability jobs. This avoids
    // fragile dense subject-level eigendecomposition while keeping exact scan
    // available after variance components are estimated.
    bool use_stcg_vce = (!is_bgen_ && !is_pgen_ && plink_reader_ &&
                         (options_.repeatability_mode == RepeatabilityMode::STCG ||
                          (options_.repeatability_mode == RepeatabilityMode::Auto && n_recs_kept_ > 80000)));

    if (use_stcg_vce) {
        // STCG_VCE: matrix-free VCE using stochastic trace + PCG
        std::cout << "  Using STCG_VCE (matrix-free) for null model...\n";
        Eigen::MatrixXd Xd = X_obs_.cast<double>();
        Eigen::VectorXd beta_ols = Xd.colPivHouseholderQr().solve(y_obs_);
        Eigen::VectorXd resid_ols = y_obs_ - Xd * beta_ols;
        double total_var_obs = resid_ols.squaredNorm() / std::max(1, n_recs_kept_ - (int)X_obs_.cols());

        double vg_init = std::max(total_var_obs * 0.5, 1e-4);
        double vpe_init = std::max(total_var_obs * 0.3, 1e-4);
        double ve_init = std::max(total_var_obs * 0.2, 1e-4);

        VCEConfig stcg_cfg;
        stcg_cfg.verbose = true;
        stcg_cfg.max_iter = options_.reml_maxit;
        stcg_cfg.tol = options_.reml_tol;
        stcg_cfg.vce_mode = "stcg";
        stcg_cfg.set_algorithm_from_string(stcg_cfg.vce_mode);
        stcg_cfg.print_report = true;

        STCG_VCE stcg_vce(gen_records_, fd_, components_, y_obs_, stcg_cfg);
        stcg_vce.setInitialVariances({vg_init, vpe_init}, ve_init);
        stcg_vce.solve();

        v_e_ = stcg_vce.getVarE();
        auto vars = stcg_vce.getVarsU();
        if (vars.size() >= 2) {
            v_g_ = vars[0]; v_pe_ = vars[1];
        }

        double total_var = v_e_ + v_g_ + v_pe_;
        std::cout << "V(e)  = " << v_e_ << "\nV(G)  = " << v_g_ << "\nV(pe) = " << v_pe_
                  << "\nV(p)  = " << total_var << "\nh2    = " << v_g_ / total_var
                  << "\nr2    = " << (v_g_ + v_pe_) / total_var << " (Repeatability)\n";
        return;
    }

    // Standard AI_REML path (Exact/Fast modes, or STCG without geno_mat)
    // The repeatability model VCE has two stable local optima:
    //   A) Correct: Ve dominates (~70% of total), Vg and Vpe share the remainder
    //   B) Spurious: Vg+Vpe absorb most variance, Ve is squeezed to <10%
    // The solver can converge to either depending on initialization.
    Eigen::MatrixXd Xd_init = X_obs_.cast<double>();
    Eigen::VectorXd beta_init = Xd_init.colPivHouseholderQr().solve(y_obs_);
    Eigen::VectorXd resid_init = y_obs_ - Xd_init * beta_init;
    double weight_sum_init = (record_weights_.size() == n_recs_kept_) ? record_weights_.sum() : (double)n_recs_kept_;
    double total_var_obs = (record_weights_.size() == n_recs_kept_)
        ? (record_weights_.array() * resid_init.array().square()).sum() / std::max(1.0, weight_sum_init - (double)X_obs_.cols())
        : resid_init.squaredNorm() / std::max(1, n_recs_kept_ - (int)X_obs_.cols());
    if (!std::isfinite(total_var_obs) || total_var_obs <= 0.0) {
        double mean_y = y_obs_.mean();
        total_var_obs = (y_obs_.array() - mean_y).square().sum() / std::max(1, n_recs_kept_ - 1);
    }
    std::cout << "  OLS residual variance = " << total_var_obs << " (X cols=" << X_obs_.cols() << ", n_recs=" << n_recs_kept_ << ")\n";
    const Eigen::MatrixXd Xd_vce = X_obs_.cast<double>();
    std::unique_ptr<Eigen::MatrixXd> V_exact_buf;
    double subject_vce_best_logL = -std::numeric_limits<double>::infinity();
    int subject_vce_best_status = 999;

    const int full_record_refine_max_records = 1000;
    const int full_record_refine_max_subjects = 300;

    {
        std::vector<int> keep_subject_indices = obs_to_sample_;
        std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
        keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
        const int q = (int)keep_subject_indices.size();

        std::vector<int> sample_to_rank(n_samples_, -1);
        for (int i = 0; i < q; ++i) sample_to_rank[keep_subject_indices[i]] = i;

        Eigen::VectorXd D = Eigen::VectorXd::Zero(q);
        Eigen::VectorXd y_sum = Eigen::VectorXd::Zero(q);
        Eigen::MatrixXd X_sum = Eigen::MatrixXd::Zero(q, Xd_vce.cols());
        for (int i = 0; i < n_recs_kept_; ++i) {
            int r = sample_to_rank[obs_to_sample_[i]];
            double wi = (record_weights_.size() == n_recs_kept_) ? record_weights_[i] : 1.0;
            D[r] += wi;
            y_sum[r] += wi * y_obs_[i];
            X_sum.row(r) += wi * Xd_vce.row(i);
        }
        double log_det_weight = 0.0;
        for (int i = 0; i < n_recs_kept_; ++i) {
            double wi = (record_weights_.size() == n_recs_kept_) ? record_weights_[i] : 1.0;
            if (!(wi > 0.0) || !std::isfinite(wi)) return;
            log_det_weight += std::log(wi);
        }

        Eigen::MatrixXd G_used(q, q);
        for (int i = 0; i < q; ++i)
            for (int j = 0; j < q; ++j)
                G_used(i, j) = G_subject_(keep_subject_indices[i], keep_subject_indices[j]);

        std::cout << "  Subject-level VCE (eigendecomp of " << q << "x" << q << " GRM)...\n";
        auto t_eig_start = std::chrono::high_resolution_clock::now();
        G_used = 0.5 * (G_used + G_used.transpose()).eval();
        Eigen::MatrixXd G_used_clean = G_used;

        auto dense_ldlt_logL = [&](double vg, double vpe, double ve) -> double {
            auto dense_fail = [&](const char* reason) -> double {
                (void)reason;
                return -std::numeric_limits<double>::infinity();
            };
            if (!(vg >= 0.0) || !(vpe >= 0.0) || !(ve > 0.0)) return -std::numeric_limits<double>::infinity();
            const int n = n_recs_kept_;
            const int p = Xd_vce.cols();

            auto factorize_spd_bent = [](Eigen::MatrixXd& A, Eigen::LLT<Eigen::MatrixXd>& solver) -> bool {
                A = 0.5 * (A + A.transpose()).eval();
                const double scale = std::max(1.0, A.diagonal().cwiseAbs().mean());
                for (double rel : {0.0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0}) {
                    Eigen::MatrixXd B = A;
                    if (rel > 0.0) B.diagonal().array() += rel * scale;
                    solver.compute(B);
                    if (solver.info() == Eigen::Success) {
                        A.swap(B);
                        return true;
                    }
                }
                return false;
            };
            auto llt_logdet = [](const Eigen::LLT<Eigen::MatrixXd>& solver) -> double {
                Eigen::MatrixXd L = solver.matrixL();
                return 2.0 * L.diagonal().array().log().sum();
            };

            Eigen::MatrixXd K = vg * G_used;
            K.diagonal().array() += vpe;
            Eigen::LLT<Eigen::MatrixXd> k_solver;
            if (!factorize_spd_bent(K, k_solver)) return dense_fail("K_factorize");

            Eigen::MatrixXd Iq = Eigen::MatrixXd::Identity(q, q);
            Eigen::MatrixXd Kinv = k_solver.solve(Iq);
            Eigen::MatrixXd M = Kinv;
            M.diagonal().array() += D.array() / ve;
            Eigen::LLT<Eigen::MatrixXd> m_solver;
            if (!factorize_spd_bent(M, m_solver)) return dense_fail("M_factorize");

            Eigen::VectorXd Minv_y_sum = m_solver.solve(y_sum / ve);
            Eigen::MatrixXd Minv_X_sum(q, p);
            for (int j = 0; j < p; ++j) Minv_X_sum.col(j) = m_solver.solve(X_sum.col(j) / ve);
            if (!Minv_y_sum.allFinite() || !Minv_X_sum.allFinite()) return dense_fail("M_solve_nonfinite");

            Eigen::VectorXd Vinv_y(n);
            Eigen::MatrixXd VinvX(n, p);
            for (int i = 0; i < n; ++i) {
                int r = sample_to_rank[obs_to_sample_[i]];
                double wi = (record_weights_.size() == n) ? record_weights_[i] : 1.0;
                Vinv_y[i] = wi * (y_obs_[i] - Minv_y_sum[r]) / ve;
                for (int j = 0; j < p; ++j) VinvX(i, j) = wi * (Xd_vce(i, j) - Minv_X_sum(r, j)) / ve;
            }

            Eigen::MatrixXd XtVinvX = Xd_vce.transpose() * VinvX;
            Eigen::LLT<Eigen::MatrixXd> x_solver;
            if (!factorize_spd_bent(XtVinvX, x_solver)) return dense_fail("XtVinvX_factorize");

            double log_det_v = n * std::log(ve) - log_det_weight +
                               llt_logdet(k_solver) + llt_logdet(m_solver);
            double log_det_xt = llt_logdet(x_solver);
            Eigen::VectorXd beta = x_solver.solve(Xd_vce.transpose() * Vinv_y);
            Eigen::VectorXd Py = Vinv_y - VinvX * beta;
            double yPy = y_obs_.dot(Py);
            double logL = -0.5 * (log_det_v + log_det_xt + yPy + (n - p) * std::log(2.0 * M_PI));
            return std::isfinite(logL) ? logL : dense_fail("logL_nonfinite");
        };

        auto run_dense_ldlt_fallback = [&]() -> bool {
            std::cerr << "WARN: Falling back to dense LDLT repeatability VCE.\n";
            std::vector<std::pair<double,double>> inits = {
                {0.5, 0.3}, {0.3, 0.2}, {0.2, 0.1}, {0.1, 0.3}, {0.05, 0.25}
            };
            double best_logL = -std::numeric_limits<double>::infinity();
            double best_vg = 0.0, best_vpe = 0.0, best_ve = total_var_obs;
            for (const auto& init : inits) {
                double lv[3] = {
                    std::log(std::max(total_var_obs * init.first, 1e-8)),
                    std::log(std::max(total_var_obs * init.second, 1e-8)),
                    std::log(std::max(total_var_obs * (1.0 - init.first - init.second), 1e-8))
                };
                double current = dense_ldlt_logL(std::exp(lv[0]), std::exp(lv[1]), std::exp(lv[2]));
                for (double step : {0.4, 0.2}) {
                    bool improved = true;
                    int local_iter = 0;
                    while (improved && local_iter++ < 4) {
                        improved = false;
                        for (int d = 0; d < 3; ++d) {
                            for (double dir : {-1.0, 1.0}) {
                                double cand[3] = {lv[0], lv[1], lv[2]};
                                cand[d] += dir * step;
                                double logL = dense_ldlt_logL(std::exp(cand[0]), std::exp(cand[1]), std::exp(cand[2]));
                                if (std::isfinite(logL) && logL > current + 1e-8) {
                                    lv[0] = cand[0]; lv[1] = cand[1]; lv[2] = cand[2];
                                    current = logL;
                                    improved = true;
                                }
                            }
                        }
                    }
                }
                double vg = std::exp(lv[0]);
                double vpe = std::exp(lv[1]);
                double ve = std::exp(lv[2]);
                std::cout << "    dense fallback init (vg_f=" << init.first << ",vpe_f=" << init.second
                          << "): Vg=" << vg << " Vpe=" << vpe << " Ve=" << ve
                          << " logL=" << current << "\n";
                if (std::isfinite(current) && current > best_logL) {
                    best_logL = current; best_vg = vg; best_vpe = vpe; best_ve = ve;
                }
            }
            if (!std::isfinite(best_logL)) return false;
            v_g_ = best_vg; v_pe_ = best_vpe; v_e_ = best_ve;
            subject_vce_best_logL = best_logL;
            subject_vce_best_status = 0;
            return true;
        };

        auto run_moment_fallback = [&]() {
            std::cerr << "WARN: Falling back to residual moment repeatability VCE.\n";
            Eigen::VectorXd r_sum = Eigen::VectorXd::Zero(q);
            Eigen::VectorXd r_sumsq = Eigen::VectorXd::Zero(q);
            Eigen::VectorXd r_n = Eigen::VectorXd::Zero(q);
            for (int i = 0; i < n_recs_kept_; ++i) {
                int rr = sample_to_rank[obs_to_sample_[i]];
                double e = resid_init[i];
                r_sum[rr] += e;
                r_sumsq[rr] += e * e;
                r_n[rr] += 1.0;
            }
            double grand = resid_init.mean();
            double within_ss = 0.0;
            double between_ss = 0.0;
            double sum_n2 = 0.0;
            for (int i = 0; i < q; ++i) {
                if (r_n[i] <= 0.0) continue;
                double mean_i = r_sum[i] / r_n[i];
                within_ss += r_sumsq[i] - r_sum[i] * mean_i;
                between_ss += r_n[i] * (mean_i - grand) * (mean_i - grand);
                sum_n2 += r_n[i] * r_n[i];
            }
            double msw = within_ss / std::max(1.0, (double)n_recs_kept_ - (double)q);
            double msb = between_ss / std::max(1.0, (double)q - 1.0);
            double kbar = ((double)n_recs_kept_ - sum_n2 / std::max(1.0, (double)n_recs_kept_)) / std::max(1.0, (double)q - 1.0);
            double repeat = (msb - msw) / std::max(1e-12, msb + (kbar - 1.0) * msw);
            repeat = std::max(0.0, std::min(0.95, repeat));
            double total = std::max(total_var_obs, 1e-8);
            v_g_ = std::max(total * repeat * 0.25, 1e-8);
            v_pe_ = std::max(total * repeat * 0.75, 1e-8);
            v_e_ = std::max(total * (1.0 - repeat), 1e-8);
            subject_vce_best_logL = -std::numeric_limits<double>::infinity();
            subject_vce_best_status = 1;
            std::cout << "    moment fallback: Vg=" << v_g_ << " Vpe=" << v_pe_
                      << " Ve=" << v_e_ << " repeatability=" << repeat << "\n";
        };

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;
        set_blas_threads(1);
        eig_solver.compute(G_used);
        if (eig_solver.info() != Eigen::Success) {
            for (double rel : {0.0, 1e-10, 1e-8, 1e-6, 1e-4}) {
                Eigen::MatrixXd retry = G_used_clean;
                if (rel > 0.0) {
                    double scale = std::max(1.0, retry.diagonal().cwiseAbs().mean());
                    retry.diagonal().array() += rel * scale;
                }
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> retry_solver;
                retry_solver.compute(retry);
                if (retry_solver.info() == Eigen::Success) {
                    eig_solver = std::move(retry_solver);
                    if (rel > 0.0) std::cerr << "WARN: GRM eigendecomposition succeeded after diagonal jitter " << rel << ".\n";
                    break;
                }
            }
        }
        set_blas_threads(options_.n_threads);
        if (eig_solver.info() != Eigen::Success) {
            std::cerr << "ERROR: GRM eigendecomposition failed.\n";
            G_used = G_used_clean;
            if (!run_dense_ldlt_fallback()) run_moment_fallback();
            auto t_vce_end = std::chrono::high_resolution_clock::now();
            double vce_seconds = std::chrono::duration<double>(t_vce_end - t_eig_start).count();
            std::cout << "  Fallback VCE done in " << vce_seconds << "s\n";
            return;
        }
        Eigen::VectorXd eig_vals = eig_solver.eigenvalues();
        Eigen::MatrixXd eig_vecs = eig_solver.eigenvectors();
        auto t_eig_end = std::chrono::high_resolution_clock::now();
        std::cout << "  GRM eigendecomposition done in " << std::chrono::duration<double>(t_eig_end - t_eig_start).count() << "s\n";

        double d_mean = D.mean();

        auto fast_subject_logL = [&](double vg, double vpe, double ve) -> double {
            const int n = n_recs_kept_;
            const int p = Xd_vce.cols();
            auto debug_fast_fail = [&](const char* reason, const Eigen::MatrixXd* mat, const Eigen::VectorXd* diag) {
                static int debug_fast_fail_count = 0;
                if (debug_fast_fail_count >= 12) return;
                ++debug_fast_fail_count;
                std::ostringstream data;
                data << "{"
                     << "\"reason\":\"" << reason << "\""
                     << ",\"vg\":" << dbg_num(vg)
                     << ",\"vpe\":" << dbg_num(vpe)
                     << ",\"ve\":" << dbg_num(ve)
                     << ",\"n_records\":" << n_recs_kept_
                     << ",\"n_subjects\":" << q
                     << ",\"x_cols\":" << p
                     << ",\"D_min\":" << dbg_num(D.size() ? D.minCoeff() : std::numeric_limits<double>::quiet_NaN())
                     << ",\"D_max\":" << dbg_num(D.size() ? D.maxCoeff() : std::numeric_limits<double>::quiet_NaN());
                if (mat) {
                    auto [diag_min, diag_max] = dbg_diag_minmax(*mat);
                    data << ",\"mat_rows\":" << mat->rows()
                         << ",\"mat_cols\":" << mat->cols()
                         << ",\"mat_nonfinite\":" << dbg_nonfinite_count(*mat)
                         << ",\"mat_diag_min\":" << dbg_num(diag_min)
                         << ",\"mat_diag_max\":" << dbg_num(diag_max);
                }
                if (diag && diag->size() > 0) {
                    data << ",\"diag_min\":" << dbg_num(diag->minCoeff())
                         << ",\"diag_max\":" << dbg_num(diag->maxCoeff())
                         << ",\"diag_nonpos\":" << (int)((diag->array() <= 0.0).count());
                }
                data << "}";
                std::cerr << "Repeatability subject VCE diagnostics: "
                          << data.str() << std::endl;
            };
            if (!(vg >= 0.0) || !(vpe >= 0.0) || !(ve > 0.0)) return -std::numeric_limits<double>::infinity();

            Eigen::VectorXd kappa = vg * eig_vals.array() + vpe;
            if ((kappa.array() <= 0.0).any()) {
                debug_fast_fail("kappa_nonpositive", nullptr, &kappa);
                return -std::numeric_limits<double>::infinity();
            }
            Eigen::VectorXd inv_kappa = kappa.cwiseInverse();

            Eigen::VectorXd Kinv_diag(q);
            for (int i = 0; i < q; ++i) {
                double s = 0.0;
                for (int j = 0; j < q; ++j) s += eig_vecs(i, j) * eig_vecs(i, j) * inv_kappa[j];
                Kinv_diag[i] = s;
            }

            Eigen::VectorXd M_diag = Kinv_diag + D / ve;
            Eigen::VectorXd M_diag_inv = M_diag.cwiseInverse();

            auto M_mult = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
                Eigen::VectorXd Utx = eig_vecs.transpose() * x;
                Eigen::VectorXd Kinv_x = eig_vecs * (inv_kappa.asDiagonal() * Utx);
                return Kinv_x + D.cwiseProduct(x) / ve;
            };

            auto pcg_solve = [&](const Eigen::VectorXd& b, int max_iter = 50) -> Eigen::VectorXd {
                const double b_norm = b.norm();
                if (b_norm == 0.0) return Eigen::VectorXd::Zero(b.size());
                Eigen::VectorXd x = M_diag_inv.cwiseProduct(b);
                Eigen::VectorXd Ax = M_mult(x);
                Eigen::VectorXd r = b - Ax;
                if (r.norm() < 1e-10 * b_norm) return x;
                Eigen::VectorXd z = M_diag_inv.cwiseProduct(r);
                Eigen::VectorXd p = z;
                double rsold = r.dot(z);
                if (!std::isfinite(rsold) || rsold < 0.0) {
                    return Eigen::VectorXd::Constant(b.size(), std::numeric_limits<double>::quiet_NaN());
                }

                for (int it = 0; it < max_iter; ++it) {
                    Eigen::VectorXd Ap = M_mult(p);
                    double denom = p.dot(Ap);
                    if (!std::isfinite(denom) || denom <= 0.0) {
                        return Eigen::VectorXd::Constant(b.size(), std::numeric_limits<double>::quiet_NaN());
                    }
                    double alpha = rsold / denom;
                    x += alpha * p;
                    r -= alpha * Ap;
                    if (r.norm() < 1e-10 * b_norm) break;
                    z = M_diag_inv.cwiseProduct(r);
                    double rsnew = r.dot(z);
                    if (!std::isfinite(rsnew) || rsnew < 0.0) {
                        return Eigen::VectorXd::Constant(b.size(), std::numeric_limits<double>::quiet_NaN());
                    }
                    p = z + (rsnew / rsold) * p;
                    rsold = rsnew;
                }
                return x;
            };

            Eigen::VectorXd Minv_y_sum = pcg_solve(y_sum / ve);
            if (!Minv_y_sum.allFinite()) return -std::numeric_limits<double>::infinity();
            Eigen::MatrixXd Minv_X_sum(q, p);
            for (int j = 0; j < p; ++j) {
                Minv_X_sum.col(j) = pcg_solve(X_sum.col(j) / ve);
                if (!Minv_X_sum.col(j).allFinite()) return -std::numeric_limits<double>::infinity();
            }

            Eigen::VectorXd Vinv_y(n);
            Eigen::MatrixXd VinvX(n, p);
            for (int i = 0; i < n; ++i) {
                int r = sample_to_rank[obs_to_sample_[i]];
                double wi = (record_weights_.size() == n) ? record_weights_[i] : 1.0;
                Vinv_y[i] = wi * (y_obs_[i] - Minv_y_sum[r]) / ve;
                for (int j = 0; j < p; ++j) VinvX(i, j) = wi * (Xd_vce(i, j) - Minv_X_sum(r, j)) / ve;
            }

            Eigen::MatrixXd XtVinvX = Xd_vce.transpose() * VinvX;
            Eigen::LDLT<Eigen::MatrixXd> x_solver;
            if (!factorize_with_jitter(XtVinvX, x_solver)) {
                debug_fast_fail("XtVinvX_factorize_fail", &XtVinvX, nullptr);
                return -std::numeric_limits<double>::infinity();
            }
            Eigen::VectorXd Dx = x_solver.vectorD();
            if ((Dx.array() <= 0.0).any()) {
                debug_fast_fail("XtVinvX_nonpositive_diag", &XtVinvX, &Dx);
                return -std::numeric_limits<double>::infinity();
            }

            double log_det_K = kappa.array().log().sum();
            double log_det_M;
            if (q <= 6000) {
                Eigen::MatrixXd M_exact = eig_vecs * inv_kappa.asDiagonal() * eig_vecs.transpose();
                for (int k = 0; k < q; ++k) M_exact(k, k) += D[k] / ve;
                Eigen::LDLT<Eigen::MatrixXd> M_ldlt(M_exact);
                if (M_ldlt.info() != Eigen::Success) {
                    debug_fast_fail("M_exact_factorize_fail", &M_exact, nullptr);
                    return -std::numeric_limits<double>::infinity();
                }
                Eigen::VectorXd Dm = M_ldlt.vectorD();
                if ((Dm.array() <= 0.0).any()) {
                    debug_fast_fail("M_exact_nonpositive_diag", &M_exact, &Dm);
                    return -std::numeric_limits<double>::infinity();
                }
                log_det_M = Dm.array().log().sum();
            } else {
                Eigen::VectorXd approx_M_eigvals = inv_kappa.array() + d_mean / ve;
                log_det_M = approx_M_eigvals.array().log().sum();
            }
            double log_det_v = n * std::log(ve) - log_det_weight + log_det_K + log_det_M;
            double log_det_xt = Dx.array().log().sum();
            Eigen::VectorXd beta = x_solver.solve(Xd_vce.transpose() * Vinv_y);
            Eigen::VectorXd Py = Vinv_y - VinvX * beta;
            double yPy = y_obs_.dot(Py);
            double logL = -0.5 * (log_det_v + log_det_xt + yPy + (n - p) * std::log(2.0 * M_PI));
            if (!std::isfinite(logL)) {
                debug_fast_fail("logL_nonfinite", &XtVinvX, &Dx);
                return -std::numeric_limits<double>::infinity();
            }
            return logL;
        };

        auto optimize_vce = [&](double vg0, double vpe0, double ve0) -> std::pair<int, double> {
            double lv[3] = {std::log(std::max(vg0, 1e-8)), std::log(std::max(vpe0, 1e-8)), std::log(std::max(ve0, 1e-8))};
            double current = fast_subject_logL(std::exp(lv[0]), std::exp(lv[1]), std::exp(lv[2]));
            if (!std::isfinite(current)) return {-1, -std::numeric_limits<double>::infinity()};

            double step = 0.3;
            int max_iter = std::max(1, options_.reml_maxit);
            for (int iter = 0; iter < max_iter && step > 1e-5; ++iter) {
                bool improved = false;
                for (int d = 0; d < 3; ++d) {
                    for (double dir : {-1.0, 1.0}) {
                        double cand[3] = {lv[0], lv[1], lv[2]};
                        cand[d] += dir * step;
                        double logL = fast_subject_logL(std::exp(cand[0]), std::exp(cand[1]), std::exp(cand[2]));
                        if (std::isfinite(logL) && logL > current + 1e-8) {
                            lv[0] = cand[0]; lv[1] = cand[1]; lv[2] = cand[2];
                            current = logL;
                            improved = true;
                        }
                    }
                }
                if (!improved) step *= 0.5;
            }

            v_g_ = std::exp(lv[0]); v_pe_ = std::exp(lv[1]); v_e_ = std::exp(lv[2]);
            double total = v_g_ + v_pe_ + v_e_;
            int status = 0;
            if (!std::isfinite(current) || total <= 0.0) status = 3;
            else if (v_g_ / total < 0.001) status = 1;
            else if (v_e_ / total < 0.1) status = 2;
            return {status, current};
        };

        std::vector<std::pair<double,double>> sub_inits = {
            {0.5, 0.3}, {0.3, 0.2}, {0.2, 0.1}
        };
        if (options_.vce_refine_mode == VceRefineMode::FullRecordExact ||
            options_.vce_refine_mode == VceRefineMode::SubjectLevelExact) {
            sub_inits.insert(sub_inits.end(), {{0.1, 0.3}, {0.4, 0.4}, {0.15, 0.15}, {0.05, 0.4}, {0.6, 0.1}});
        }
        double best_logL = -std::numeric_limits<double>::infinity();
        int best_status = 999;
        double best_vg = 0, best_vpe = 0, best_ve = total_var_obs;
        auto t_vce_start = std::chrono::high_resolution_clock::now();
        for (size_t ti = 0; ti < sub_inits.size(); ++ti) {
            auto& [vg_f, vpe_f] = sub_inits[ti];
            double vg0 = std::max(total_var_obs * vg_f, 1e-4);
            double vpe0 = std::max(total_var_obs * vpe_f, 1e-4);
            double ve0 = std::max(total_var_obs * (1.0 - vg_f - vpe_f), 1e-4);
            auto [status, logL] = optimize_vce(vg0, vpe0, ve0);
            std::cout << "    init #" << ti << " (vg_f=" << vg_f << ",vpe_f=" << vpe_f
                      << "): Vg=" << v_g_ << " Vpe=" << v_pe_ << " Ve=" << v_e_
                      << " logL=" << logL << " status=" << status << "\n";
            if (status >= 0 && std::isfinite(logL) &&
                (status < best_status || (status == best_status && logL > best_logL))) {
                best_logL = logL; best_status = status;
                best_vg = v_g_; best_vpe = v_pe_; best_ve = v_e_;
            }
        }
        auto t_vce_end = std::chrono::high_resolution_clock::now();
        double vce_seconds = std::chrono::duration<double>(t_vce_end - t_vce_start).count();
        std::cout << "  Subject-level VCE done in " << vce_seconds << "s\n";

        if (best_status == 999 || !std::isfinite(best_logL)) {
            std::cerr << "WARN: All subject-level repeatability VCE starts failed.\n";
            if (!run_dense_ldlt_fallback()) run_moment_fallback();
            return;
        }

        v_g_ = best_vg; v_pe_ = best_vpe; v_e_ = best_ve;
        subject_vce_best_logL = best_logL;
        subject_vce_best_status = best_status;
    }

    auto exact_repeatability_logL = [&](double vg, double vpe, double ve) -> double {
        if (!(vg >= 0.0) || !(vpe >= 0.0) || !(ve > 0.0)) return -std::numeric_limits<double>::infinity();
        const int n = n_recs_kept_;
        const int p = Xd_vce.cols();
        if (!V_exact_buf) V_exact_buf = std::make_unique<Eigen::MatrixXd>(n, n);
        Eigen::MatrixXd& V_exact = *V_exact_buf;
        for (int i = 0; i < n; ++i) {
            int si = obs_to_sample_[i];
            for (int j = 0; j <= i; ++j) {
                int sj = obs_to_sample_[j];
                double value = vg * G_subject_(si, sj);
                if (si == sj) value += vpe;
                if (i == j) {
                    double wi = (record_weights_.size() == n) ? record_weights_[i] : 1.0;
                    value += ve / wi;
                }
                V_exact(i, j) = value;
                V_exact(j, i) = value;
            }
        }

        Eigen::LDLT<Eigen::MatrixXd> v_solver;
        if (!factorize_with_jitter(V_exact, v_solver)) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd Dv = v_solver.vectorD();
        if ((Dv.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();
        double log_det_v = Dv.array().log().sum();

        Eigen::MatrixXd VinvX = v_solver.solve(Xd_vce);
        Eigen::VectorXd Vinvy = v_solver.solve(y_obs_);
        Eigen::MatrixXd XtVinvX = Xd_vce.transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> x_solver;
        if (!factorize_with_jitter(XtVinvX, x_solver)) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd Dx = x_solver.vectorD();
        if ((Dx.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();
        double log_det_xt = Dx.array().log().sum();
        Eigen::VectorXd beta = x_solver.solve(Xd_vce.transpose() * Vinvy);
        Eigen::VectorXd Py = Vinvy - VinvX * beta;
        double yPy = y_obs_.dot(Py);
        double logL = -0.5 * (log_det_v + log_det_xt + yPy + (n - p) * std::log(2.0 * M_PI));
        return std::isfinite(logL) ? logL : -std::numeric_limits<double>::infinity();
    };

    auto subject_repeatability_logL = [&](double vg, double vpe, double ve) -> double {
        if (!(vg >= 0.0) || !(vpe >= 0.0) || !(ve > 0.0)) return -std::numeric_limits<double>::infinity();
        const int n = n_recs_kept_;
        const int p = X_obs_.cols();

        std::vector<int> keep_subject_indices = obs_to_sample_;
        std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
        keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
        const int q = (int)keep_subject_indices.size();
        if (q <= 0 || n <= p) return -std::numeric_limits<double>::infinity();

        std::vector<int> sample_to_rank(n_samples_, -1);
        for (int i = 0; i < q; ++i) sample_to_rank[keep_subject_indices[i]] = i;

        Eigen::VectorXd D = Eigen::VectorXd::Zero(q);
        Eigen::VectorXd y_sum = Eigen::VectorXd::Zero(q);
        Eigen::MatrixXd X_sum = Eigen::MatrixXd::Zero(q, p);
        double log_det_weight = 0.0;
        for (int i = 0; i < n; ++i) {
            int r = sample_to_rank[obs_to_sample_[i]];
            if (r < 0) return -std::numeric_limits<double>::infinity();
            double wi = (record_weights_.size() == n) ? record_weights_[i] : 1.0;
            if (!(wi > 0.0) || !std::isfinite(wi)) return -std::numeric_limits<double>::infinity();
            D[r] += wi;
            y_sum[r] += wi * y_obs_[i];
            X_sum.row(r) += wi * Xd_vce.row(i);
            log_det_weight += std::log(wi);
        }
        if ((D.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();

        Eigen::MatrixXd G_used(q, q);
        for (int i = 0; i < q; ++i) {
            for (int j = 0; j < q; ++j) G_used(i, j) = G_subject_(keep_subject_indices[i], keep_subject_indices[j]);
        }
        Eigen::MatrixXd K = vg * G_used;
        K.diagonal().array() += vpe;

        Eigen::LDLT<Eigen::MatrixXd> k_solver;
        if (!factorize_with_jitter(K, k_solver)) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd Dk = k_solver.vectorD();
        if ((Dk.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();
        Eigen::MatrixXd M = k_solver.solve(Eigen::MatrixXd::Identity(q, q));
        M.diagonal().array() += D.array() / ve;
        Eigen::LDLT<Eigen::MatrixXd> m_solver;
        if (!factorize_with_jitter(M, m_solver)) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd Dm = m_solver.vectorD();
        if ((Dm.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();

        Eigen::VectorXd corr_y = m_solver.solve(y_sum / ve);
        Eigen::MatrixXd corr_X(q, p);
        for (int j = 0; j < p; ++j) corr_X.col(j) = m_solver.solve(X_sum.col(j) / ve);

        Eigen::VectorXd Vinvy(n);
        Eigen::MatrixXd VinvX(n, p);
        for (int i = 0; i < n; ++i) {
            int r = sample_to_rank[obs_to_sample_[i]];
            double wi = (record_weights_.size() == n) ? record_weights_[i] : 1.0;
            Vinvy[i] = wi * (y_obs_[i] - corr_y[r]) / ve;
            for (int j = 0; j < p; ++j) VinvX(i, j) = wi * (Xd_vce(i, j) - corr_X(r, j)) / ve;
        }

        Eigen::MatrixXd XtVinvX = Xd_vce.transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> x_solver;
        if (!factorize_with_jitter(XtVinvX, x_solver)) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd Dx = x_solver.vectorD();
        if ((Dx.array() <= 0.0).any()) return -std::numeric_limits<double>::infinity();

        double log_det_v = n * std::log(ve) - log_det_weight + Dk.array().log().sum() + Dm.array().log().sum();
        double log_det_xt = Dx.array().log().sum();
        Eigen::VectorXd beta = x_solver.solve(Xd_vce.transpose() * Vinvy);
        Eigen::VectorXd Py = Vinvy - VinvX * beta;
        double yPy = y_obs_.dot(Py);
        double logL = -0.5 * (log_det_v + log_det_xt + yPy + (n - p) * std::log(2.0 * M_PI));
        return std::isfinite(logL) ? logL : -std::numeric_limits<double>::infinity();
    };

    auto refine_subject_reml = [&](double& vg, double& vpe, double& ve, double& best_logL) -> bool {
        std::vector<int> keep_subject_indices = obs_to_sample_;
        std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
        keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
        if (keep_subject_indices.empty()) return false;

        double lv[3] = {std::log(std::max(vg, 1e-8)), std::log(std::max(vpe, 1e-8)), std::log(std::max(ve, 1e-8))};
        double current = subject_repeatability_logL(std::exp(lv[0]), std::exp(lv[1]), std::exp(lv[2]));
        if (!std::isfinite(current)) return false;
        double step = 0.5;
        int max_iter = std::max(1, options_.reml_maxit);
        for (int iter = 0; iter < max_iter && step > 1e-4; ++iter) {
            bool improved = false;
            for (int d = 0; d < 3; ++d) {
                for (double dir : {-1.0, 1.0}) {
                    double cand[3] = {lv[0], lv[1], lv[2]};
                    cand[d] += dir * step;
                    double logL = subject_repeatability_logL(std::exp(cand[0]), std::exp(cand[1]), std::exp(cand[2]));
                    if (std::isfinite(logL) && logL > current + 1e-7) {
                        lv[0] = cand[0]; lv[1] = cand[1]; lv[2] = cand[2];
                        current = logL;
                        improved = true;
                    }
                }
            }
            if (!improved) step *= 0.5;
        }
        vg = std::exp(lv[0]);
        vpe = std::exp(lv[1]);
        ve = std::exp(lv[2]);
        best_logL = current;
        return true;
    };

    auto refine_full_record_reml = [&](double& vg, double& vpe, double& ve, double& best_logL) -> bool {
        std::vector<int> keep_subject_indices = obs_to_sample_;
        std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
        keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
        if (options_.vce_refine_mode != VceRefineMode::FullRecordExact &&
            (n_recs_kept_ > full_record_refine_max_records ||
             keep_subject_indices.size() > (size_t)full_record_refine_max_subjects)) return false;

        const double min_var = std::max(total_var_obs * 1e-10, 1e-10);
        const double base_logL = exact_repeatability_logL(vg, vpe, ve);
        if (!std::isfinite(base_logL)) return false;
        struct StartPoint { double vg, vpe, ve; };
        std::vector<StartPoint> starts = {
            {std::max(vg, min_var), std::max(vpe, min_var), std::max(ve, min_var)},
            {std::max(1.0, min_var), std::max(1.0, min_var), std::max(8.0, min_var)},
            {std::max(0.5, min_var), std::max(0.5, min_var), std::max(10.0, min_var)},
            {std::max(2.0, min_var), std::max(2.0, min_var), std::max(6.0, min_var)},
            {std::max(0.1, min_var), std::max(3.0, min_var), std::max(8.0, min_var)},
            {std::max(total_var_obs * 0.10, min_var), std::max(total_var_obs * 0.10, min_var), std::max(total_var_obs * 0.80, min_var)},
            {std::max(total_var_obs * 0.05, min_var), std::max(total_var_obs * 0.20, min_var), std::max(total_var_obs * 0.75, min_var)},
            {std::max(total_var_obs * 0.20, min_var), std::max(total_var_obs * 0.05, min_var), std::max(total_var_obs * 0.75, min_var)},
            {std::max(total_var_obs * 0.15, min_var), std::max(total_var_obs * 0.25, min_var), std::max(total_var_obs * 0.60, min_var)},
        };

        double global_best = base_logL;
        double best_lv[3] = {std::log(std::max(vg, min_var)), std::log(std::max(vpe, min_var)), std::log(std::max(ve, min_var))};
        for (const auto& start : starts) {
            double lv[3] = {std::log(start.vg), std::log(start.vpe), std::log(start.ve)};
            double current = exact_repeatability_logL(std::exp(lv[0]), std::exp(lv[1]), std::exp(lv[2]));
            if (!std::isfinite(current)) continue;

            double step = 1.0;
            int max_iter = std::max(1, options_.reml_maxit);
            for (int iter = 0; iter < max_iter && step > 5e-5; ++iter) {
                bool improved = false;
                for (int d = 0; d < 3; ++d) {
                    for (double dir : {-1.0, 1.0}) {
                        double cand[3] = {lv[0], lv[1], lv[2]};
                        cand[d] += dir * step;
                        double cvg = std::max(std::exp(cand[0]), min_var);
                        double cvpe = std::max(std::exp(cand[1]), min_var);
                        double cve = std::max(std::exp(cand[2]), min_var);
                        double logL = exact_repeatability_logL(cvg, cvpe, cve);
                        if (std::isfinite(logL) && logL > current + 1e-8) {
                            lv[0] = std::log(cvg); lv[1] = std::log(cvpe); lv[2] = std::log(cve);
                            current = logL;
                            improved = true;
                        }
                    }
                }
                if (!improved) step *= 0.5;
            }

            if (current > global_best) {
                global_best = current;
                best_lv[0] = lv[0]; best_lv[1] = lv[1]; best_lv[2] = lv[2];
            }
        }

        if (!std::isfinite(global_best) || global_best <= base_logL + 1e-7) return false;
        vg = std::exp(best_lv[0]);
        vpe = std::exp(best_lv[1]);
        ve = std::exp(best_lv[2]);
        best_logL = global_best;
        return true;
    };

    auto use_full_record_vce = [&]() -> bool {
        if (options_.vce_refine_mode == VceRefineMode::FullRecordExact) return true;
        return false;
    };

    auto use_full_record_logL = [&]() -> bool {
        if (options_.vce_refine_mode == VceRefineMode::FullRecordExact) return true;
        return false;
    };

    auto use_full_record_refine = [&]() -> bool {
        if (options_.vce_refine_mode == VceRefineMode::FullRecordExact) return true;
        return false;
    };

    auto run_vce_from_values = [&](double vg_init, double vpe_init, double ve_init,
                                   int max_iter, double tol) -> std::pair<int, double> {
        if (!use_full_record_vce()) {
            return {0, 0.0};
        }
        VCEConfig vce_cfg; vce_cfg.verbose = true; vce_cfg.max_iter = max_iter;
        vce_cfg.tol = tol; vce_cfg.vce_mode = "em";
        vce_cfg.set_algorithm_from_string(vce_cfg.vce_mode);
        vce_cfg.print_report = false;
        AI_REML solver(gen_records_, fd_, components_, y_obs_, vce_cfg);
        solver.setInitialVariances({vg_init, vpe_init}, ve_init);
        solver.solve();
        v_e_ = solver.getVarE();
        auto vars = solver.getVarsU();
        if (vars.size() < 2) return {-1, -std::numeric_limits<double>::infinity()};
        v_g_ = vars[0]; v_pe_ = vars[1];
        double total = v_g_ + v_pe_ + v_e_;
        double logL = use_full_record_logL()
            ? exact_repeatability_logL(v_g_, v_pe_, v_e_)
            : subject_repeatability_logL(v_g_, v_pe_, v_e_);
        int status = 0;
        if (!std::isfinite(logL) || total <= 0.0) status = 3;
        else if (v_g_ / total < 0.001) status = 1;
        else if (v_e_ / total < 0.1) status = 2;
        return {status, logL};
    };

    auto run_single_vce = [&](double vg_frac, double vpe_frac) -> std::pair<int, double> {
        double vg_init = std::max(total_var_obs * vg_frac, 1e-4);
        double vpe_init = std::max(total_var_obs * vpe_frac, 1e-4);
        double ve_init  = std::max(total_var_obs * (1.0 - vg_frac - vpe_frac), 1e-4);
        return run_vce_from_values(vg_init, vpe_init, ve_init, options_.reml_maxit, options_.reml_tol);
    };

    // Multi-start VCE: try multiple initializations, pick the one with highest logL
    struct TrialResult { double vg, vpe, ve, logL; int status; };
    std::vector<TrialResult> trials;

    if (use_full_record_vce()) {
        std::vector<std::pair<double,double>> inits = {
            {0.5, 0.3},
            {0.7, 0.2},
            {0.2, 0.2},
            {0.3, 0.5},
            {0.4, 0.1},
            {0.1, 0.1},
            {0.05, 0.05},
        };
        for (auto& [vg_f, vpe_f] : inits) {
            auto [status, logL] = run_single_vce(vg_f, vpe_f);
            trials.push_back({v_g_, v_pe_, v_e_, logL, status});
        }
    } else {
        trials.push_back({v_g_, v_pe_, v_e_, subject_vce_best_logL, subject_vce_best_status});
    }

    // Select best: prefer valid biological/numerical regimes, then exact REML logL.
    int best_idx = 0;
    double best_logL = -std::numeric_limits<double>::infinity();
    for (int preferred_status : {0, 1, 2}) {
        for (int i = 0; i < (int)trials.size(); ++i) {
            auto& t = trials[i];
            if (t.status == preferred_status && std::isfinite(t.logL) && t.logL > best_logL) {
                best_logL = t.logL;
                best_idx = i;
            }
        }
        if (best_logL > -std::numeric_limits<double>::infinity()) break;
    }
    // If no finite model-specific logL exists, pick the least pathological trial.
    if (best_logL == -std::numeric_limits<double>::infinity()) {
        int best_status = 999;
        for (int i = 0; i < (int)trials.size(); ++i) {
            if (trials[i].status < best_status) {
                best_status = trials[i].status;
                best_logL = trials[i].logL;
                best_idx = i;
            }
        }
    }

    v_g_ = trials[best_idx].vg;
    v_pe_ = trials[best_idx].vpe;
    v_e_ = trials[best_idx].ve;

    {
        double start_vg = v_g_, start_vpe = v_pe_, start_ve = v_e_;
        double refined_logL = best_logL;
        if (options_.vce_refine_mode == VceRefineMode::Fast || options_.vce_refine_mode == VceRefineMode::Auto) {
            LOG_INFO("Repeatability VCE refine mode=" << (options_.vce_refine_mode == VceRefineMode::Auto ? "auto" : "fast")
                     << ", using subject-level Woodbury VCE without extra exact refinement");
        } else if (use_full_record_refine() && refine_full_record_reml(v_g_, v_pe_, v_e_, refined_logL)) {
            best_logL = refined_logL;
            LOG_INFO("Repeatability full-record exact REML refinement accepted (logL=" << refined_logL << ")");
        } else if (options_.vce_refine_mode != VceRefineMode::FullRecordExact &&
                   refine_subject_reml(v_g_, v_pe_, v_e_, refined_logL)) {
            best_logL = refined_logL;
            LOG_INFO("Repeatability subject-level exact REML refinement accepted (logL=" << refined_logL << ")");
        } else if (options_.vce_refine_mode == VceRefineMode::FullRecordExact) {
            LOG_INFO("Repeatability full-record exact REML refinement did not improve, keeping best EM estimate");
        } else {
            v_g_ = start_vg;
            v_pe_ = start_vpe;
            v_e_ = start_ve;
            LOG_INFO("Repeatability subject-level exact REML refinement skipped");
        }
    }

    LOG_INFO("Multi-start VCE: best init #" << best_idx << " (logL=" << best_logL
             << ", status=" << trials[best_idx].status << ")");

    double total_var = v_e_ + v_g_ + v_pe_;
    std::cout << "V(e)  = " << v_e_ << "\nV(G)  = " << v_g_ << "\nV(pe) = " << v_pe_
              << "\nV(p)  = " << total_var << "\nh2    = " << v_g_ / total_var
              << "\nr2    = " << (v_g_ + v_pe_) / total_var << " (Repeatability)\n";
}

void RepeatabilityModel::run() {
    // IASBLUP-style strategy selection
    RepeatabilityMode mode = options_.repeatability_mode;
    if (mode == RepeatabilityMode::Auto) {
        // Auto-select based on data size
        if (n_recs_kept_ > 80000) {
            mode = RepeatabilityMode::STCG;
            LOG_INFO("Repeatability auto-mode: n_records=" << n_recs_kept_ << " > 80K -> STCG (matrix-free)");
        } else if (options_.fast_lmm || n_recs_kept_ > 20000) {
            mode = RepeatabilityMode::Fast;
            LOG_INFO("Repeatability auto-mode: n_records=" << n_recs_kept_ << " -> Fast (aggregated mean)");
        } else {
            mode = RepeatabilityMode::Exact;
            LOG_INFO("Repeatability auto-mode: n_records=" << n_recs_kept_ << " -> Exact (record-level GLS)");
        }
    }
    if (options_.eb_lmm && mode != RepeatabilityMode::Fast) {
        std::cout << "Repeatability EB-LMM uses the validated fast aggregated subject-level scan; switching mode to fast.\n";
        mode = RepeatabilityMode::Fast;
    }

    // Print solver report for repeatability
    SolverReport report;
    report.task = TaskType::Repeatability;
    report.mme_dim = n_recs_kept_;
    report.n_components = 3;  // G + PE + E
    switch (mode) {
        case RepeatabilityMode::Exact:
            report.strategy = VceStrategy::MMEBased;
            report.solver = SolverType::SparseLDLT;
            report.trace_mode = TraceMode::Exact;
            report.strategy_reason = "R1: Record-level exact GLS";
            break;
        case RepeatabilityMode::Fast:
            report.strategy = VceStrategy::VBased;
            report.solver = SolverType::DenseVEigen;
            report.trace_mode = TraceMode::Exact;
            report.strategy_reason = "R2: Aggregated mean + Eigen null";
            break;
        case RepeatabilityMode::STCG:
            report.strategy = VceStrategy::MatrixFree;
            report.solver = SolverType::STCG;
            report.trace_mode = TraceMode::Hutch;
            report.strategy_reason = "R3: Matrix-free STCG";
            break;
        default: break;
    }
    print_solver_report(report);

    switch (mode) {
        case RepeatabilityMode::Exact:
            scanExact();
            break;
        case RepeatabilityMode::Fast:
            scanFast();
            break;
        case RepeatabilityMode::STCG:
            scanSTCG();
            break;
        default:
            scanExact();
            break;
    }
}

void RepeatabilityModel::writeNullSummary(const std::string& filename, const std::string& scan_mode,
                                          const std::vector<std::pair<std::string, std::string>>& extra_fields) const {
    std::ofstream ofs(filename);
    bool extra_has_grm_cache_status = false;
    for (const auto& [key, value] : extra_fields) {
        if (key == "grm_cache_status") {
            extra_has_grm_cache_status = true;
            break;
        }
    }
    double total_var = v_g_ + v_pe_ + v_e_;
    std::vector<char> observed_subjects(n_samples_, 0);
    int n_subjects_kept = 0;
    for (int sample_index : obs_to_sample_) {
        if (sample_index >= 0 && sample_index < n_samples_ && !observed_subjects[sample_index]) {
            observed_subjects[sample_index] = 1;
            ++n_subjects_kept;
        }
    }
    ofs << "Item\tValue\n"
        << "n_records\t" << n_recs_kept_ << "\n"
        << "n_subjects\t" << n_subjects_kept << "\n"
        << "Vg\t" << v_g_ << "\nVpe\t" << v_pe_ << "\nVe\t" << v_e_ << "\nVp\t" << total_var << "\n"
        << "h2\t" << ((total_var > 0.0) ? (v_g_ / total_var) : 0.0) << "\n"
        << "r2\t" << ((total_var > 0.0) ? ((v_g_ + v_pe_) / total_var) : 0.0) << "\n"
        << "kernel_precision\t" << (options_.mixed_precision ? "float32_mixed" : "float64_exact") << "\n"
        << "grm_cache_prefix\t" << (options_.grm_cache_prefix.empty() ? "NA" : options_.grm_cache_prefix) << "\n";
    if (!extra_has_grm_cache_status) {
        ofs << "grm_cache_status\t" << (options_.grm_cache_prefix.empty() ? "disabled" : "not_used") << "\n";
    }
    ofs
        << "vce_refine_mode\t"
        << (options_.vce_refine_mode == VceRefineMode::FullRecordExact ? "full-record-exact" :
            options_.vce_refine_mode == VceRefineMode::SubjectLevelExact ? "subject-level-exact" :
            options_.vce_refine_mode == VceRefineMode::Fast ? "fast" : "auto") << "\n"
        << "time_grm_sec\t" << grm_time_sec_ << "\n"
        << "time_vce_sec\t" << vce_time_sec_ << "\n"
        << "time_scan_vbuild_sec\t" << scan_vbuild_time_sec_ << "\n"
        << "time_snp_scan_sec\t" << scan_snp_time_sec_ << "\n"
        << "scan_mode\t" << scan_mode << "\n";
    for (const auto& [key, value] : extra_fields) {
        ofs << key << "\t" << value << "\n";
    }
}

} // namespace cosmic
