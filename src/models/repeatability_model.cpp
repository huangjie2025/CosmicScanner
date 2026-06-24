#include "repeatability_model.h"
#include "vce.h"
#include "stcg_vce.h"
#include "logger.h"
#include "matrix_io.h"
#include "matrix_free_pcg.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <deque>
#include <sstream>
#include <sys/resource.h>
#include <omp.h>
#include <filesystem>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "grm_builder.h"
#include "pheno_reader.h"
#include "covar_reader.h"
#include "data_filter.h"
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
    opts.eb_lmm = cfg.eb_lmm;
    opts.mixed_precision = cfg.repeatability_mixed_precision;
    opts.refit_validation = cfg.repeatability_refit_validation;
    opts.refit_limit = cfg.repeatability_refit_limit;
    opts.grm_cache_prefix = cfg.repeatability_grm_cache_prefix;
    opts.weight_file = cfg.repeatability_weight_file;
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
    return opts;
}

namespace {

// #region debug-point helper:repeatability-woodbury
static std::string dbg_read_env_value(const std::string& file_path,
                                      const std::string& key,
                                      const std::string& fallback) {
    std::ifstream ifs(file_path);
    if (!ifs) return fallback;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.rfind(key + "=", 0) == 0) return line.substr(key.size() + 1);
    }
    return fallback;
}

static std::string dbg_json_escape(std::string s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (char c : s) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            case '\'': out += " "; break;
            default: out += c; break;
        }
    }
    return out;
}

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

static void dbg_report_event(const std::string& hypothesis_id,
                             const std::string& location,
                             const std::string& msg,
                             const std::string& data_json) {
    const std::string env_file = ".dbg/repeatability-woodbury.env";
    const std::string server_url = dbg_read_env_value(env_file, "DEBUG_SERVER_URL", "http://127.0.0.1:7777/event");
    const std::string session_id = dbg_read_env_value(env_file, "DEBUG_SESSION_ID", "repeatability-woodbury");
    std::ostringstream payload;
    payload << "{\"sessionId\":\"" << dbg_json_escape(session_id)
            << "\",\"runId\":\"pre-fix\""
            << ",\"hypothesisId\":\"" << dbg_json_escape(hypothesis_id)
            << "\",\"location\":\"" << dbg_json_escape(location)
            << "\",\"msg\":\"" << dbg_json_escape(msg)
            << "\",\"data\":" << data_json << "}";
    std::string cmd = "curl -s -X POST '" + server_url
        + "' -H 'Content-Type: application/json' --data '"
        + payload.str() + "' >/dev/null 2>&1";
    (void)std::system(cmd.c_str());
}
// #endregion

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

static double wald_pvalue_1df(double z_score) {
    return std::erfc(std::abs(z_score) / std::sqrt(2.0));
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

static std::vector<std::string> sortChroms(const std::set<std::string>& chroms_set) {
    std::vector<std::string> chroms(chroms_set.begin(), chroms_set.end());
    std::sort(chroms.begin(), chroms.end(), [](const std::string& a, const std::string& b) {
        try { return std::stoi(a) < std::stoi(b); } catch (...) { return a < b; }
    });
    return chroms;
}

// Algorithm 0 (Yang/GCTA default), 3 (Yang et al), 5 (TrueVar) apply per-SNP
// standardization (divide by sqrt(2pq) or sqrt(var)), so each SNP contributes
// equally to the GRM trace. LOCO scaling uses SNP count for these.
// Algorithms 1 (VanRaden), 2 (Zeng), 4 (Vitezica) only center, so LOCO
// scaling uses sum_2pq to account for heterogeneous per-SNP variance.
static bool uses_per_snp_scaling(int alg) {
    return alg == 0 || alg == 3 || alg == 5;
}

static void writeAssocHeader(const std::string& filename, bool out_ssf, bool eb_lmm = false) {
    std::ofstream ofs(filename);
    if (eb_lmm && !out_ssf) {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tEB_BETA\tEB_SE\tEB_WALD\tP\tLOG10P"
            << "\tMARKER_LAMBDA\tPRIOR_VAR\tPOST_VAR\tEB_CONFIDENCE\n";
    } else if (out_ssf) {
        ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }
}

static void appendAssocResult(const std::string& filename, const GwasResult& res, bool out_ssf, bool eb_lmm = false) {
    std::ofstream ofs(filename, std::ios::app);
    if (std::isnan(res.p_value)) return;
    if (eb_lmm && !out_ssf) {
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.eb_wald << "\t" << res.p_value << "\t" << res.log10_p << "\t"
            << res.eb_marker_lambda << "\t" << res.eb_prior_var << "\t"
            << res.eb_post_var << "\t" << res.eb_confidence << "\n";
    } else if (out_ssf) {
        ofs << res.chrom << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.p_value << "\t" << res.snp_id << "\n";
    } else {
        ofs << res.chrom << "\t" << res.snp_id << "\t" << res.pos << "\t"
            << res.a1 << "\t" << res.a2 << "\t" << res.freq << "\t"
            << res.n_obs << "\t" << res.beta << "\t" << res.se << "\t"
            << res.p_value << "\t" << res.log10_p << "\n";
    }
}

static void appendAssocResults(const std::string& filename, const std::vector<GwasResult>& results, bool out_ssf, bool eb_lmm = false) {
    for (const auto& res : results) appendAssocResult(filename, res, out_ssf, eb_lmm);
}

static void accumulateEbResults(const std::vector<GwasResult>& results, double& effective_tests, long long& tested_markers) {
    for (const auto& res : results) {
        if (std::isnan(res.p_value)) continue;
        effective_tests += std::max(0.0, std::min(1.0, res.eb_confidence));
        tested_markers++;
    }
}

static void writeRepeatabilityEbSummary(const std::string& filename, double effective_tests, long long tested_markers) {
    std::ofstream ofs(filename);
    if (!ofs) throw std::runtime_error("Cannot write repeatability EB summary: " + filename);
    ofs << std::setprecision(12);
    double m_eff = std::max(effective_tests, 1.0);
    ofs << "Metric\tValue\n";
    ofs << "method\tRepeatability EB-LMM random-SNP scan (fast aggregated subject-level)\n";
    ofs << "tested_markers\t" << tested_markers << "\n";
    ofs << "effective_tests\t" << effective_tests << "\n";
    ofs << "effective_tests_per_marker\t" << (tested_markers > 0 ? effective_tests / (double)tested_markers : 0.0) << "\n";
    ofs << "bonferroni_0.05\t" << (tested_markers > 0 ? 0.05 / (double)tested_markers : 0.0) << "\n";
    ofs << "eb_mbonf_0.05\t" << (0.05 / m_eff) << "\n";
    ofs << "note\tExperimental repeatability EB uses aggregated weighted subject means; exact record-level EB is not used in this mode.\n";
}

struct NullModelCache {
    Eigen::MatrixXd VinvX;
    Eigen::MatrixXd XtVinvX_inv;
    Eigen::VectorXd P_y;
    Eigen::MatrixXd Minv;
    Eigen::VectorXd Wz;
    Eigen::MatrixXd H_sub;
    Eigen::VectorXd y_sub_proj;
    Eigen::MatrixXd X_sub_proj;
    Eigen::MatrixXd K_subject;
    Eigen::VectorXd Dz;
    Eigen::MatrixXd X_sub_sum;
    Eigen::VectorXd Py_sub_sum;
    double ve = 1.0;
};

template <typename IndexToMeta, typename ResultAppender>
static void processExactSubjectBlock(const Eigen::MatrixXd& G_block,
                                     const std::vector<double>& freqs,
                                     const std::vector<int>& valid_ns,
                                     IndexToMeta index_to_meta,
                                     const NullModelCache& cache,
                                     int n_subjects_keep,
                                     ResultAppender append_result) {
    if (G_block.cols() == 0) return;
    Eigen::MatrixXd rhs = G_block;
    for (int s = 0; s < n_subjects_keep; ++s) rhs.row(s) *= cache.Wz[s] / cache.ve;
    Eigen::MatrixXd correction = cache.Minv * rhs;
    Eigen::MatrixXd VinvG = (G_block - correction) / cache.ve;
    Eigen::RowVectorXd scores = G_block.transpose() * cache.Py_sub_sum;
    Eigen::MatrixXd XtVinvG = cache.X_sub_sum.transpose() * VinvG;

    for (int c = 0; c < G_block.cols(); ++c) {
        double g_Vinv_g = (cache.Wz.array() * G_block.col(c).array() * VinvG.col(c).array()).sum();
        Eigen::VectorXd xvg = XtVinvG.col(c);
        double term2 = xvg.transpose() * cache.XtVinvX_inv * xvg;
        double var_score = g_Vinv_g - term2;

        GwasResult res;
        index_to_meta(c, res);
        res.freq = freqs[c];
        res.n_obs = n_subjects_keep;
        if (valid_ns[c] > 0 && var_score > 1e-9) {
            double se = std::sqrt(var_score);
            double stat = scores[c] / se;
            res.beta = scores[c] / var_score;
            res.se = 1.0 / se;
            res.p_value = chisq_prob(stat * stat, 1);
            res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
        }
        append_result(res);
    }
}

static NullModelCache buildNullModelCache(const Eigen::MatrixXd& V,
                                          const Eigen::VectorXd& y,
                                          const Eigen::MatrixXf& X) {
    NullModelCache cache;
    Eigen::LDLT<Eigen::MatrixXd> v_solver;
    if (!factorize_with_jitter(V, v_solver)) {
        throw std::runtime_error("Repeatability scan failed to factorize covariance matrix V.");
    }
    Eigen::MatrixXd Xd = X.cast<double>();
    cache.VinvX = v_solver.solve(Xd);
    Eigen::VectorXd Vinv_y = v_solver.solve(y);
    Eigen::MatrixXd XtVinvX = Xd.transpose() * cache.VinvX;
    Eigen::LDLT<Eigen::MatrixXd> x_solver(XtVinvX);
    if (x_solver.info() == Eigen::Success) {
        cache.XtVinvX_inv = x_solver.solve(Eigen::MatrixXd::Identity(XtVinvX.rows(), XtVinvX.cols()));
    } else {
        cache.XtVinvX_inv = XtVinvX.inverse();
    }
    Eigen::VectorXd beta_null = cache.XtVinvX_inv * (Xd.transpose() * Vinv_y);
    cache.P_y = Vinv_y - cache.VinvX * beta_null;
    return cache;
}

static std::string sanitizeChromLabel(std::string chrom) {
    for (char& c : chrom) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
    }
    return chrom;
}

static std::string getCacheFilePath(const std::string& prefix, const std::string& label) {
    return prefix + "." + label + ".f32grm.bin";
}

static std::string getCacheMetaPath(const std::string& prefix, const std::string& label) {
    return prefix + "." + label + ".meta.txt";
}

struct MappedFloat32Matrix {
    int fd = -1;
    size_t mapped_bytes = 0;
    void* mapped_ptr = MAP_FAILED;
    const float* data = nullptr;
    uint32_t rows = 0, cols = 0;

    ~MappedFloat32Matrix() { close(); }
    void close() {
        if (mapped_ptr != MAP_FAILED) { munmap(mapped_ptr, mapped_bytes); mapped_ptr = MAP_FAILED; }
        if (fd >= 0) { ::close(fd); fd = -1; }
        data = nullptr; mapped_bytes = 0; rows = 0; cols = 0;
    }
    bool open(const std::string& file) {
        close();
        fd = ::open(file.c_str(), O_RDONLY);
        if (fd < 0) return false;
        struct stat st;
        if (fstat(fd, &st) != 0 || st.st_size < (off_t)(sizeof(uint32_t) * 2)) { close(); return false; }
        mapped_bytes = (size_t)st.st_size;
        mapped_ptr = mmap(nullptr, mapped_bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapped_ptr == MAP_FAILED) { close(); return false; }
        const uint32_t* header = reinterpret_cast<const uint32_t*>(mapped_ptr);
        rows = header[0]; cols = header[1];
        if (rows == 0 || cols == 0) { close(); return false; }
        size_t expected = sizeof(uint32_t) * 2 + sizeof(float) * (size_t)rows * (size_t)cols;
        if (mapped_bytes < expected) { close(); return false; }
        data = reinterpret_cast<const float*>(header + 2);
        return true;
    }
};

static bool readMatrixFloat32Mmap(const std::string& file, Eigen::MatrixXf& out) {
    MappedFloat32Matrix mapped;
    if (!mapped.open(file)) return false;
    Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> view(mapped.data, (int)mapped.rows, (int)mapped.cols);
    out = view;
    return true;
}

static void writeMatrixFloat32(const std::string& file, const Eigen::MatrixXf& mat) {
    std::filesystem::path p(file);
    if (!p.parent_path().empty()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(file, std::ios::binary);
    if (!ofs) throw std::runtime_error("Cannot write GRM cache file: " + file);
    uint32_t rows = (uint32_t)mat.rows(), cols = (uint32_t)mat.cols();
    ofs.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    ofs.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    ofs.write(reinterpret_cast<const char*>(mat.data()), sizeof(float) * rows * cols);
}

static bool readCacheMeta(const std::string& file, double* scaling_factor, std::set<std::string>* observed_chroms) {
    std::ifstream ifs(file);
    if (!ifs) return false;
    std::string key;
    while (ifs >> key) {
        if (key == "scale" && scaling_factor) { ifs >> *scaling_factor; }
        else if (key == "chroms" && observed_chroms) {
            observed_chroms->clear();
            std::string chrom_line; std::getline(ifs, chrom_line);
            std::istringstream iss(chrom_line); std::string chrom;
            while (iss >> chrom) observed_chroms->insert(chrom);
        } else { std::string dummy; std::getline(ifs, dummy); }
    }
    return true;
}

static void writeCacheMeta(const std::string& file, double scaling_factor, const std::set<std::string>& observed_chroms) {
    std::filesystem::path p(file);
    if (!p.parent_path().empty()) std::filesystem::create_directories(p.parent_path());
    std::ofstream ofs(file);
    if (!ofs) throw std::runtime_error("Cannot write GRM cache meta file: " + file);
    ofs << "scale " << scaling_factor << "\nchroms";
    for (const auto& chrom : observed_chroms) ofs << " " << chrom;
    ofs << "\n";
}

template <typename ResultAppender>
static void scanPlinkExact(const PlinkReader& reader,
                           const std::string& chrom_filter,
                           const std::vector<int>& keep_subject_indices,
                           const std::vector<int>& obs_to_subject,
                           const NullModelCache& cache,
                            int n_subjects_keep,
                            ResultAppender append_result) {
    const auto& bim = reader.getBimInfo();
    int m = (int)bim.size();
    const int block_size = 1024;
    for (int start = 0; start < m; start += block_size) {
        int end = std::min(m, start + block_size);
        std::vector<int> snp_indices;
        snp_indices.reserve(end - start);
        for (int j = start; j < end; ++j) {
            if (chrom_filter.empty() || bim[j].chrom == chrom_filter) snp_indices.push_back(j);
        }
        if (snp_indices.empty()) continue;

        int bsz = (int)snp_indices.size();
        Eigen::MatrixXd G_block(n_subjects_keep, bsz);
        std::vector<double> freqs(bsz, 0.0);
        std::vector<int> valid_ns(bsz, 0);

        #pragma omp parallel for schedule(dynamic, 16)
        for (int c = 0; c < bsz; ++c) {
            Eigen::VectorXf g_raw;
            int j = snp_indices[c];
            reader.getSnpGenotypes(j, g_raw);
            double sum_g = 0.0;
            int valid_n = 0;
            for (int s = 0; s < n_subjects_keep; ++s) {
                float val = g_raw[keep_subject_indices[s]];
                G_block(s, c) = val;
                if (!std::isnan(val)) { sum_g += val; ++valid_n; }
            }
            freqs[c] = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
            valid_ns[c] = valid_n;
            double mean = 2.0 * freqs[c];
            for (int s = 0; s < n_subjects_keep; ++s) {
                if (std::isnan(G_block(s, c))) G_block(s, c) = 0.0;
                else G_block(s, c) -= mean;
            }
        }

        processExactSubjectBlock(G_block, freqs, valid_ns,
            [&](int c, GwasResult& res) {
                int j = snp_indices[c];
                res.chrom = bim[j].chrom; res.snp_id = bim[j].id; res.pos = bim[j].bp_pos;
                res.a1 = bim[j].alt; res.a2 = bim[j].ref;
            }, cache, n_subjects_keep, append_result);
    }
}

template <typename Reader, typename ReadFn, typename ResultAppender>
static void scanSequentialExact(Reader& reader, size_t m, ReadFn read_variant,
                                const std::string& chrom_filter,
                                const std::vector<int>& keep_subject_indices,
                                const std::vector<int>& obs_to_subject,
                                const NullModelCache& cache,
                                int n_subjects_keep,
                                ResultAppender append_result) {
    reader.reset();
    const size_t CHUNK_SIZE = 1024;
    struct SnpData {
        std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
    };
    std::vector<SnpData> chunk_data; chunk_data.reserve(CHUNK_SIZE);
    size_t variants_read = 0;
    while (variants_read < m) {
        chunk_data.clear();
        bool hit_eof = false;
        for (size_t i = 0; i < CHUNK_SIZE && variants_read < m; ++i) {
            SnpData d; bool ok = read_variant(d.chrom, d.rsid, d.pos, d.ref, d.alt, d.dosage);
            if (!ok) { hit_eof = true; break; }
            variants_read++;
            if (!chrom_filter.empty() && d.chrom != chrom_filter) continue;
            chunk_data.push_back(std::move(d));
        }
        if (chunk_data.empty()) {
            if (hit_eof) break;
            continue;
        }
        int bsz = (int)chunk_data.size();
        Eigen::MatrixXd G_block(n_subjects_keep, bsz);
        std::vector<double> freqs(bsz, 0.0);
        std::vector<int> valid_ns(bsz, 0);
        #pragma omp parallel for schedule(dynamic, 16)
        for (int i = 0; i < bsz; ++i) {
            const auto& d = chunk_data[i];
            double sum_g = 0.0;
            int valid_n = 0;
            for (int s = 0; s < n_subjects_keep; ++s) {
                float val = d.dosage[keep_subject_indices[s]];
                G_block(s, i) = val;
                if (!std::isnan(val)) { sum_g += val; ++valid_n; }
            }
            freqs[i] = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
            valid_ns[i] = valid_n;
            double mean = 2.0 * freqs[i];
            for (int s = 0; s < n_subjects_keep; ++s) {
                if (std::isnan(G_block(s, i))) G_block(s, i) = 0.0;
                else G_block(s, i) -= mean;
            }
        }
        processExactSubjectBlock(G_block, freqs, valid_ns,
            [&](int c, GwasResult& res) {
                const auto& d = chunk_data[c];
                res.chrom = d.chrom; res.snp_id = d.rsid; res.pos = d.pos;
                res.a1 = d.alt; res.a2 = d.ref;
            }, cache, n_subjects_keep, append_result);
    }
}

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
        dbg_report_event("A", "repeatability_model.cpp:loadCovariates",
                         "[DEBUG] repeatability design matrix summary", data.str());
    }
    // #endregion
}

void RepeatabilityModel::buildGRM() {
    std::cout << "Calculating GRM...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.algorithm = options_.grm_alg;
    #ifdef _OPENMP
    grm_opts.thread_num = omp_get_max_threads();
    #endif

    {
        GrmBuilder builder(grm_opts);
        if (is_pgen_) { pgen_reader_->reset(); builder.compute(*pgen_reader_); }
        else if (is_bgen_) { bgen_reader_->reset(); builder.compute(*bgen_reader_); }
        else { builder.compute(*plink_reader_); }

        std::cout << "Inverting GRM...\n";
        Eigen::MatrixXd G = builder.getGrm().cast<double>();
        G_subject_ = G;
        G.diagonal().array() += 0.001;
        Eigen::LLT<Eigen::MatrixXd> lltG(G);
        if (lltG.info() == Eigen::Success) {
            Ginv_ = lltG.solve(Eigen::MatrixXd::Identity(G.rows(), G.cols()));
        } else {
            Ginv_ = G.inverse();
        }
    }
    set_blas_threads(options_.n_threads);
}

void RepeatabilityModel::buildMMEComponents() {
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

    // Check if STCG_VCE path is available (geno_mat attached)
    bool use_stcg_vce = false;

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
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(G_used);
        if (eig_solver.info() != Eigen::Success) {
            std::cerr << "ERROR: GRM eigendecomposition failed.\n";
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
                dbg_report_event("B", "repeatability_model.cpp:fast_subject_logL",
                                 "[DEBUG] repeatability subject VCE failure", data.str());
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
                Eigen::VectorXd x = M_diag_inv.cwiseProduct(b);
                Eigen::VectorXd Ax = M_mult(x);
                Eigen::VectorXd r = b - Ax;
                Eigen::VectorXd z = M_diag_inv.cwiseProduct(r);
                Eigen::VectorXd p = z;
                double rsold = r.dot(z);

                for (int it = 0; it < max_iter; ++it) {
                    Eigen::VectorXd Ap = M_mult(p);
                    double alpha = rsold / p.dot(Ap);
                    x += alpha * p;
                    r -= alpha * Ap;
                    if (r.norm() < 1e-10 * b.norm()) break;
                    z = M_diag_inv.cwiseProduct(r);
                    double rsnew = r.dot(z);
                    p = z + (rsnew / rsold) * p;
                    rsold = rsnew;
                }
                return x;
            };

            Eigen::VectorXd Minv_y_sum = pcg_solve(y_sum / ve);
            Eigen::MatrixXd Minv_X_sum(q, p);
            for (int j = 0; j < p; ++j) {
                Minv_X_sum.col(j) = pcg_solve(X_sum.col(j) / ve);
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
            if (status < best_status || (status == best_status && logL > best_logL)) {
                best_logL = logL; best_status = status;
                best_vg = v_g_; best_vpe = v_pe_; best_ve = v_e_;
            }
        }
        auto t_vce_end = std::chrono::high_resolution_clock::now();
        double vce_seconds = std::chrono::duration<double>(t_vce_end - t_vce_start).count();
        std::cout << "  Subject-level VCE done in " << vce_seconds << "s\n";

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
        std::cout << "Repeatability EB-LMM is experimental and currently uses the fast aggregated subject-level scan; switching mode to fast.\n";
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

void RepeatabilityModel::scanFast() {
    auto scan_total_start = std::chrono::high_resolution_clock::now();
    std::cout << "\nRunning Repeatability GWAS Scan (fast aggregated mean model)...\n";

    struct SubjectData { double sum_y = 0.0, count = 0.0; std::vector<double> sum_X; };
    std::vector<SubjectData> subj_data(n_samples_);
    std::vector<bool> subj_present(n_samples_, false);
    for (int i = 0; i < n_recs_kept_; ++i) {
        int gen_idx = gen_records_[i].aid - 1;
        if (gen_idx < 0 || gen_idx >= n_samples_) continue;
        subj_present[gen_idx] = true; subj_data[gen_idx].sum_y += gen_records_[i].y;
        subj_data[gen_idx].count += 1.0;
        if (subj_data[gen_idx].sum_X.empty()) subj_data[gen_idx].sum_X.resize(fd_.p - 1, 0.0);
        for (const auto& kv : fd_.rows[i])
            if (kv.first > 0 && kv.first - 1 < (int)subj_data[gen_idx].sum_X.size())
                subj_data[gen_idx].sum_X[kv.first - 1] += kv.second;
    }

    std::vector<int> agg_keep_indices; std::vector<double> y_agg_vec, n_counts;
    std::vector<std::vector<double>> X_agg_vec;
    for (int i = 0; i < n_samples_; ++i) {
        if (!subj_present[i]) continue;
        agg_keep_indices.push_back(i); double n_i = subj_data[i].count;
        n_counts.push_back(n_i); y_agg_vec.push_back(subj_data[i].sum_y / n_i);
        std::vector<double> x_row; x_row.push_back(1.0);
        for (double val : subj_data[i].sum_X) x_row.push_back(val / n_i);
        X_agg_vec.push_back(std::move(x_row));
    }

    int n_agg = (int)agg_keep_indices.size();
    std::cout << "Aggregated to " << n_agg << " unique subjects for fast repeatability scan.\n";

    Eigen::VectorXd y_agg(n_agg); Eigen::MatrixXf X_agg(n_agg, fd_.p); Eigen::VectorXd obs_weights(n_agg);
    double n_mean = 0.0, n_var = 0.0, n_min = 1e30, n_max = 0.0;
    for (int i = 0; i < n_agg; ++i) {
        y_agg(i) = y_agg_vec[i];
        for (int j = 0; j < fd_.p; ++j) X_agg(i, j) = (float)X_agg_vec[i][j];
        double cnt = n_counts[i]; n_mean += cnt; n_min = std::min(n_min, cnt); n_max = std::max(n_max, cnt);
        double var_resid = v_pe_ + v_e_ / cnt;
        obs_weights(i) = (var_resid > 1e-12) ? (1.0 / std::sqrt(var_resid)) : 1.0;
    }
    if (n_agg > 0) { n_mean /= n_agg; for (double c : n_counts) n_var += (c - n_mean) * (c - n_mean); n_var /= n_agg; }
    double mean_w = obs_weights.mean(); if (mean_w > 0.0) obs_weights /= mean_w;

    std::cout << "Repeated Measures Balance Check:\n  Mean records/subject: " << n_mean
              << "\n  Min/Max records: " << n_min << " / " << n_max
              << "\n  Variance of records: " << n_var << "\n";
    if (n_var > 0.1) std::cout << "WARNING: Fast repeatability mode is approximate for unbalanced data.\n";

    DataFilter::KeepInfo agg_keep; agg_keep.keep_indices = agg_keep_indices; agg_keep.n_keep = n_agg;

    std::string agg_pheno_file = options_.out_prefix + ".agg.pheno";
    { std::ofstream ofs(agg_pheno_file); ofs << "FID\tIID\tPHENO_MEAN\tN_RECS\n";
      for (size_t i = 0; i < agg_keep_indices.size(); ++i) {
          int idx = agg_keep_indices[i]; std::string fid = "0", iid = "0";
          if (is_pgen_) { iid = pgen_reader_->getSampleIds()[idx]; fid = iid; }
          else if (is_bgen_) { iid = bgen_reader_->getSampleIds()[idx]; fid = iid; }
          else { fid = plink_reader_->getFamInfo()[idx].fid; iid = plink_reader_->getFamInfo()[idx].iid; }
          ofs << fid << "\t" << iid << "\t" << y_agg(i) << "\t" << n_counts[i] << "\n";
      }
    }

    std::vector<std::pair<std::string, std::string>> fast_summary_fields = {
        {"mean_records", std::to_string(n_mean)},
        {"min_records", std::to_string(n_min)},
        {"max_records", std::to_string(n_max)}
    };
    writeNullSummary(options_.out_prefix + ".repeatability.null.txt", "fast_repeatability_mean_model", fast_summary_fields);

    Eigen::VectorXd y_w = y_agg; Eigen::MatrixXf X_w = X_agg;
    for (int i = 0; i < n_agg; ++i) { double w = obs_weights(i); y_w(i) *= w; X_w.row(i) *= (float)w; }

    LmmEngine::Options lmm_fast_opts;
    lmm_fast_opts.fast_mode = true; lmm_fast_opts.min_maf = options_.min_maf;
    lmm_fast_opts.min_hwe = options_.min_hwe; lmm_fast_opts.use_covariate_projection = options_.use_covariate_projection;
    lmm_fast_opts.use_spa = options_.use_spa; lmm_fast_opts.thread_num = options_.n_threads;

    std::string assoc_file = options_.out_prefix + ".assoc.txt";
    writeAssocHeader(assoc_file, options_.out_ssf, options_.eb_lmm);
    double eb_effective_tests = 0.0;
    long long eb_tested_markers = 0;

    GrmBuilder::Options grm_opts; grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe; grm_opts.algorithm = options_.grm_alg;
    #ifdef _OPENMP
    grm_opts.thread_num = omp_get_max_threads();
    #endif

    GrmBuilder global_builder(grm_opts);
    if (is_bgen_) { bgen_reader_->reset(); global_builder.compute(*bgen_reader_); }
    else if (is_pgen_) { pgen_reader_->reset(); global_builder.compute(*pgen_reader_); }
    else { global_builder.compute(*plink_reader_); }
    auto chroms_set_fast = global_builder.getObservedChroms();

    if (!options_.no_loco && chroms_set_fast.size() > 1) {
        std::cout << "Fast repeatability GWAS: LOCO scan...\n";
        double S_global_fast = (uses_per_snp_scaling(options_.grm_alg)) ? (double)global_builder.getTotalSnps() : global_builder.getSum2pq();
        const auto& G_global = global_builder.getGrm();
        auto sorted_chroms = sortChroms(chroms_set_fast);

        for (int chr_idx = 0; chr_idx < (int)sorted_chroms.size(); ++chr_idx) {
            const auto& chrom = sorted_chroms[chr_idx];
            std::cout << "  Fast LOCO chromosome " << chrom << "...\n";

            GrmBuilder::Options chr_opts = grm_opts;
            chr_opts.chrom_filter = chrom;
            GrmBuilder chr_builder(chr_opts);
            if (is_bgen_) { bgen_reader_->reset(); chr_builder.compute(*bgen_reader_); }
            else if (is_pgen_) { pgen_reader_->reset(); chr_builder.compute(*pgen_reader_); }
            else { chr_builder.compute(*plink_reader_); }

            double S_chr = (uses_per_snp_scaling(options_.grm_alg)) ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
            const auto& G_chr_norm = chr_builder.getGrm();
            double S_loco = S_global_fast - S_chr;
            int n_sub = agg_keep.n_keep;
            Eigen::MatrixXd G_loco_sub(n_sub, n_sub);
            if (S_chr > 0.0 && S_loco > 1e-6) {
                double w1 = S_global_fast / S_loco;
                double w2 = S_chr / S_loco;
                #pragma omp parallel for collapse(2)
                for (int i = 0; i < n_sub; ++i) for (int j = 0; j < n_sub; ++j) {
                    int gi = agg_keep.keep_indices[i], gj = agg_keep.keep_indices[j];
                    G_loco_sub(i, j) = G_global(gi, gj) * w1 - G_chr_norm(gi, gj) * w2;
                }
            } else {
                G_loco_sub = DataFilter::subsetSquare(G_global, agg_keep);
            }
            G_loco_sub.array().colwise() *= obs_weights.array();
            G_loco_sub.array().rowwise() *= obs_weights.transpose().array();

            EigenVce::Options vce_fast_opts; vce_fast_opts.use_covariate_projection = options_.use_covariate_projection;
            vce_fast_opts.null_solver = options_.lmm_null_solver;
            vce_fast_opts.max_iter = options_.reml_maxit;
            vce_fast_opts.tol = options_.reml_tol;
            EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_loco_sub); vce_fast.runNullModel();
            LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

            std::vector<GwasResult> chrom_results;
            if (is_bgen_) { bgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*bgen_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*bgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else if (is_pgen_) { pgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*pgen_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*pgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else { if (options_.eb_lmm) engine_fast.scanEB(*plink_reader_, chrom_results, agg_keep.keep_indices, chrom); else engine_fast.scan(*plink_reader_, chrom_results, agg_keep.keep_indices, chrom); }

            appendAssocResults(assoc_file, chrom_results, options_.out_ssf, options_.eb_lmm);
            if (options_.eb_lmm) accumulateEbResults(chrom_results, eb_effective_tests, eb_tested_markers);
        }
    } else {
        std::cout << "Fast repeatability GWAS: global scan on " << n_agg << " aggregated subjects.\n";
        Eigen::MatrixXd G_fast_sub = DataFilter::subsetSquare(G_subject_, agg_keep);
        G_fast_sub.array().colwise() *= obs_weights.array();
        G_fast_sub.array().rowwise() *= obs_weights.transpose().array();

        EigenVce::Options vce_fast_opts; vce_fast_opts.use_covariate_projection = options_.use_covariate_projection;
        vce_fast_opts.null_solver = options_.lmm_null_solver;
        vce_fast_opts.max_iter = options_.reml_maxit;
        vce_fast_opts.tol = options_.reml_tol;
        EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_fast_sub); vce_fast.runNullModel();
        LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

        std::vector<GwasResult> results_fast;
        if (is_bgen_) { bgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*bgen_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*bgen_reader_, results_fast, agg_keep.keep_indices); }
        else if (is_pgen_) { pgen_reader_->reset(); if (options_.eb_lmm) engine_fast.scanEB(*pgen_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*pgen_reader_, results_fast, agg_keep.keep_indices); }
        else { if (options_.eb_lmm) engine_fast.scanEB(*plink_reader_, results_fast, agg_keep.keep_indices); else engine_fast.scan(*plink_reader_, results_fast, agg_keep.keep_indices); }
        appendAssocResults(assoc_file, results_fast, options_.out_ssf, options_.eb_lmm);
        if (options_.eb_lmm) accumulateEbResults(results_fast, eb_effective_tests, eb_tested_markers);
    }
    auto scan_total_end = std::chrono::high_resolution_clock::now();
    scan_snp_time_sec_ = std::chrono::duration<double>(scan_total_end - scan_total_start).count();
    fast_summary_fields.push_back({"time_fast_total_sec", std::to_string(scan_snp_time_sec_)});
    writeNullSummary(options_.out_prefix + ".repeatability.null.txt", "fast_repeatability_mean_model", fast_summary_fields);
    if (options_.eb_lmm) {
        std::string eb_summary = options_.out_prefix + ".repeatability.eb.summary.txt";
        writeRepeatabilityEbSummary(eb_summary, eb_effective_tests, eb_tested_markers);
        std::cout << "Repeatability EB summary saved to " << eb_summary << std::endl;
    }
    std::cout << "Fast repeatability GWAS results saved to " << assoc_file << std::endl;
}

void RepeatabilityModel::scanExact() {
    scan_vbuild_time_sec_ = 0.0;
    scan_snp_time_sec_ = 0.0;
    std::cout << "\nRunning Repeatability GWAS Scan (exact record-level GLS)...\n";

    std::vector<int> keep_subject_indices = obs_to_sample_;
    std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
    keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
    int n_subjects_keep = (int)keep_subject_indices.size();
    int n_obs = n_recs_kept_;

    std::vector<int> sample_to_subject_rank(n_samples_, -1);
    for (int i = 0; i < n_subjects_keep; ++i) sample_to_subject_rank[keep_subject_indices[i]] = i;
    std::vector<int> obs_to_subject(n_obs);
    for (int i = 0; i < n_obs; ++i) obs_to_subject[i] = sample_to_subject_rank[obs_to_sample_[i]];

    DataFilter::KeepInfo keep_subjects; keep_subjects.keep_indices = keep_subject_indices; keep_subjects.n_keep = n_subjects_keep;

    std::cout << "Repeatability GWAS: building subject-level GRM...\n";
    if (options_.mixed_precision) std::cout << "Repeatability GWAS: mixed precision enabled.\n";

    GrmBuilder::Options grm_opts; grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe; grm_opts.algorithm = options_.grm_alg;
    grm_opts.thread_num = options_.n_threads;
    std::string global_cache_status = options_.grm_cache_prefix.empty() ? "disabled" : "not_used";

    // Pre-declare variables shared between lambda and LOCO path
    bool use_single_pass_loco = false;
    GrmBuilder::SplitGrmResult split_res;
    std::vector<Eigen::MatrixXf> per_chrom_raw_subject;
    std::vector<double> chr_scaling;
    std::vector<std::string> chr_cache_st;

    auto compute_or_load_subject_grm = [&](const std::string& chrom_filter,
                                           double* scaling_factor_out,
                                           std::set<std::string>* observed_chroms_out,
                                           std::string* cache_status_out) -> Eigen::MatrixXd {
        use_single_pass_loco = options_.grm_cache_prefix.empty()
                               && !is_bgen_ && !options_.no_loco && options_.grm_alg != 2;

        if (use_single_pass_loco && chrom_filter.empty()) {
            GrmBuilder split_builder(grm_opts);
            if (is_bgen_) {
                const auto& variants = bgen_reader_->getVariants();
                std::set<std::string> all_chroms;
                for (const auto& v : variants) all_chroms.insert(v.chrom);
                bgen_reader_->reset();
                split_res = split_builder.computeSplit(*bgen_reader_, all_chroms, keep_subject_indices);
            } else if (is_pgen_) {
                const auto& variants = pgen_reader_->getVariants();
                std::set<std::string> all_chroms;
                for (const auto& v : variants) all_chroms.insert(v.chrom);
                pgen_reader_->reset();
                split_res = split_builder.computeSplit(*pgen_reader_, all_chroms, keep_subject_indices);
            } else {
                const auto& bim = plink_reader_->getBimInfo();
                std::set<std::string> all_chroms;
                for (const auto& b : bim) all_chroms.insert(b.chrom);
                split_res = split_builder.computeSplit(*plink_reader_, all_chroms);
            }

            double Sg = (uses_per_snp_scaling(options_.grm_alg))
                ? (double)split_res.global_snp_count : split_res.global_sum_2pq;
            Eigen::MatrixXd G_full = (is_bgen_ || is_pgen_)
                ? split_res.global : DataFilter::subsetSquare(split_res.global, keep_subjects);

            if (observed_chroms_out) *observed_chroms_out = split_builder.getObservedChroms();
            if (scaling_factor_out) *scaling_factor_out = Sg;
            if (cache_status_out) *cache_status_out = "single_pass";

            std::set<std::string> obs = split_builder.getObservedChroms();
            std::vector<std::string> local_chroms = sortChroms(obs);
            size_t nc = local_chroms.size();
            per_chrom_raw_subject.assign(nc, Eigen::MatrixXf());
            chr_scaling.assign(nc, 0.0);
            chr_cache_st.assign(nc, "single_pass");

            for (size_t i = 0; i < nc; ++i) {
                auto it = std::find(split_res.chrom_labels.begin(), split_res.chrom_labels.end(), local_chroms[i]);
                if (it != split_res.chrom_labels.end()) {
                    size_t k = std::distance(split_res.chrom_labels.begin(), it);
                    per_chrom_raw_subject[i] = (is_bgen_ || is_pgen_)
                        ? split_res.per_chrom_raw[k].cast<float>()
                        : DataFilter::subsetSquare(split_res.per_chrom_raw[k], keep_subjects).cast<float>();
                    chr_scaling[i] = (uses_per_snp_scaling(options_.grm_alg))
                        ? (double)split_res.per_chrom_snp_count[k] : split_res.per_chrom_sum_2pq[k];
                }
            }

            return G_full;
        }
        set_blas_threads(1);
        std::string cache_label = chrom_filter.empty() ? "global" : ("chr_" + sanitizeChromLabel(chrom_filter));
        std::string cache_file = options_.grm_cache_prefix.empty() ? "" : getCacheFilePath(options_.grm_cache_prefix, cache_label);
        std::string cache_meta = options_.grm_cache_prefix.empty() ? "" : getCacheMetaPath(options_.grm_cache_prefix, cache_label);

        if (!cache_file.empty()) {
            Eigen::MatrixXf cached_f32;
            if (readMatrixFloat32Mmap(cache_file, cached_f32)) {
                bool meta_ok = true;
                if (!cache_meta.empty()) meta_ok = readCacheMeta(cache_meta, scaling_factor_out, observed_chroms_out);
                if (meta_ok) {
                    if (cache_status_out) *cache_status_out = "mmap_load";
                    std::cout << "Repeatability GWAS: loaded GRM cache via mmap from " << cache_file << std::endl;
                    return cached_f32.cast<double>();
                }
                if (cache_status_out) *cache_status_out = "meta_miss_recompute";
            }
        }

        GrmBuilder local_builder(grm_opts);
        if (!chrom_filter.empty()) { GrmBuilder::Options lo = grm_opts; lo.chrom_filter = chrom_filter; local_builder = GrmBuilder(lo); }
        if (is_bgen_) { bgen_reader_->reset(); local_builder.compute(*bgen_reader_); }
        else if (is_pgen_) { pgen_reader_->reset(); local_builder.compute(*pgen_reader_); }
        else { local_builder.compute(*plink_reader_); }

        if (scaling_factor_out) *scaling_factor_out = (uses_per_snp_scaling(options_.grm_alg)) ? (double)local_builder.getTotalSnps() : local_builder.getSum2pq();
        if (observed_chroms_out) *observed_chroms_out = local_builder.getObservedChroms();

        Eigen::MatrixXd G_subject = DataFilter::subsetSquare(local_builder.getGrm(), keep_subjects);
        if (!cache_file.empty()) {
            writeMatrixFloat32(cache_file, G_subject.cast<float>());
            writeCacheMeta(cache_meta, scaling_factor_out ? *scaling_factor_out : 0.0, local_builder.getObservedChroms());
            if (cache_status_out) *cache_status_out = "computed_save";
            std::cout << "Repeatability GWAS: saved GRM cache to " << cache_file << std::endl;
        } else if (cache_status_out) { *cache_status_out = "disabled"; }
        return G_subject;
    };

    std::set<std::string> observed_chroms; double S_global = 0.0;
    Eigen::MatrixXd G_subject_global = compute_or_load_subject_grm("", &S_global, &observed_chroms, &global_cache_status);
    Eigen::MatrixXf G_subject_global_f32;
    if (options_.mixed_precision) G_subject_global_f32 = G_subject_global.cast<float>();
    std::vector<std::string> chroms = sortChroms(observed_chroms);
    bool loco_enabled = !options_.no_loco && chroms.size() > 1;

    std::string out_file = options_.out_prefix + ".assoc.txt";
    std::string null_file = options_.out_prefix + ".repeatability.null.txt";
    std::string cache_summary_file = options_.out_prefix + ".repeatability.cache.txt";
    std::vector<std::pair<std::string, std::string>> cache_entries;
    cache_entries.push_back({"global", global_cache_status});
    writeAssocHeader(out_file, options_.out_ssf);
    writeNullSummary(null_file, loco_enabled ? "exact_repeatability_loco" : "exact_repeatability_global",
                     {{"grm_cache_status", global_cache_status}});

    auto scan_with_subject_grm = [&](const Eigen::MatrixXd& G_subject, const std::string& chrom_filter) {
        auto v_start = std::chrono::high_resolution_clock::now();
        Eigen::MatrixXd K = v_g_ * G_subject + v_pe_ * Eigen::MatrixXd::Identity(n_subjects_keep, n_subjects_keep);
        Eigen::LDLT<Eigen::MatrixXd> k_solver;
        if (!factorize_with_jitter(K, k_solver)) {
            throw std::runtime_error("Repeatability GWAS failed to factorize subject-level K matrix.");
        }

        Eigen::VectorXd Wz = Eigen::VectorXd::Zero(n_subjects_keep);
        Eigen::VectorXd y_sub = Eigen::VectorXd::Zero(n_subjects_keep);
        Eigen::MatrixXd X_sub = Eigen::MatrixXd::Zero(n_subjects_keep, X_obs_.cols());
        for (int i = 0; i < n_obs; ++i) {
            int subj = obs_to_subject[i];
            double wi = (record_weights_.size() == n_obs) ? record_weights_[i] : 1.0;
            Wz[subj] += wi;
            y_sub[subj] += wi * y_obs_[i];
            X_sub.row(subj) += wi * X_obs_.row(i).cast<double>();
        }
        Eigen::MatrixXd M = k_solver.solve(Eigen::MatrixXd::Identity(n_subjects_keep, n_subjects_keep));
        M.diagonal().array() += Wz.array() / v_e_;
        Eigen::LDLT<Eigen::MatrixXd> m_solver;
        if (!factorize_with_jitter(M, m_solver)) {
            // #region debug-point C:woodbury-factorization
            {
                auto [k_diag_min, k_diag_max] = dbg_diag_minmax(K);
                auto [m_diag_min, m_diag_max] = dbg_diag_minmax(M);
                std::ostringstream data;
                data << "{"
                     << "\"chrom\":\"" << (chrom_filter.empty() ? "global" : chrom_filter) << "\""
                     << ",\"vg\":" << dbg_num(v_g_)
                     << ",\"vpe\":" << dbg_num(v_pe_)
                     << ",\"ve\":" << dbg_num(v_e_)
                     << ",\"n_subjects\":" << n_subjects_keep
                     << ",\"n_records\":" << n_obs
                     << ",\"Wz_min\":" << dbg_num(Wz.size() ? Wz.minCoeff() : std::numeric_limits<double>::quiet_NaN())
                     << ",\"Wz_max\":" << dbg_num(Wz.size() ? Wz.maxCoeff() : std::numeric_limits<double>::quiet_NaN())
                     << ",\"K_nonfinite\":" << dbg_nonfinite_count(K)
                     << ",\"M_nonfinite\":" << dbg_nonfinite_count(M)
                     << ",\"K_diag_min\":" << dbg_num(k_diag_min)
                     << ",\"K_diag_max\":" << dbg_num(k_diag_max)
                     << ",\"M_diag_min\":" << dbg_num(m_diag_min)
                     << ",\"M_diag_max\":" << dbg_num(m_diag_max)
                     << "}";
                dbg_report_event("C", "repeatability_model.cpp:scanExact",
                                 "[DEBUG] repeatability woodbury factorization failure", data.str());
            }
            // #endregion
            throw std::runtime_error("Repeatability GWAS failed to factorize Woodbury subject matrix.");
        }
        Eigen::MatrixXd Minv = m_solver.solve(Eigen::MatrixXd::Identity(n_subjects_keep, n_subjects_keep));

        auto woodbury_vinv_vec = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            Eigen::VectorXd ztrinv = Eigen::VectorXd::Zero(n_subjects_keep);
            for (int i = 0; i < n_obs; ++i) {
                double wi = (record_weights_.size() == n_obs) ? record_weights_[i] : 1.0;
                ztrinv[obs_to_subject[i]] += wi * v[i] / v_e_;
            }
            Eigen::VectorXd correction = Minv * ztrinv;
            Eigen::VectorXd out(n_obs);
            for (int i = 0; i < n_obs; ++i) {
                double wi = (record_weights_.size() == n_obs) ? record_weights_[i] : 1.0;
                out[i] = wi * (v[i] - correction[obs_to_subject[i]]) / v_e_;
            }
            return out;
        };

        Eigen::VectorXd Vinv_y = woodbury_vinv_vec(y_obs_);
        Eigen::MatrixXd VinvX(n_obs, X_obs_.cols());
        for (int j = 0; j < X_obs_.cols(); ++j) VinvX.col(j) = woodbury_vinv_vec(X_obs_.col(j).cast<double>());
        Eigen::MatrixXd XtVinvX = X_obs_.cast<double>().transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> x_solver;
        Eigen::MatrixXd XtVinvX_inv;
        if (factorize_with_jitter(XtVinvX, x_solver)) {
            XtVinvX_inv = x_solver.solve(Eigen::MatrixXd::Identity(XtVinvX.rows(), XtVinvX.cols()));
        } else {
            // Last resort: use invert_psd_matrix from libCosmicCore (eigenvalue correction)
            XtVinvX_inv = invert_psd_matrix(XtVinvX, 1e-4);
            LOG_WARN("XtVinvX required eigenvalue-corrected inversion (near-singular design matrix)");
        }
        Eigen::VectorXd beta_null = XtVinvX_inv * (X_obs_.cast<double>().transpose() * Vinv_y);
        Eigen::VectorXd P_y = Vinv_y - VinvX * beta_null;

        NullModelCache cache; cache.VinvX = VinvX; cache.XtVinvX_inv = XtVinvX_inv; cache.P_y = P_y;
        cache.Minv = Minv;
        cache.Wz = Wz;
        cache.X_sub_sum = X_sub;
        cache.Py_sub_sum = Eigen::VectorXd::Zero(n_subjects_keep);
        cache.ve = v_e_;
        for (int i = 0; i < n_obs; ++i) cache.Py_sub_sum[obs_to_subject[i]] += P_y[i];

        auto v_end = std::chrono::high_resolution_clock::now();
        auto scan_start = std::chrono::high_resolution_clock::now();
        set_blas_threads(1);

        std::vector<GwasResult> chunk_results;
        auto append_result = [&](const GwasResult& res) { chunk_results.push_back(res); };

        if (is_bgen_) {
            scanSequentialExact(*bgen_reader_, bgen_reader_->getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage);
                }, chrom_filter, keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
        } else if (is_pgen_) {
            scanSequentialExact(*pgen_reader_, pgen_reader_->getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage);
                }, chrom_filter, keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
        } else {
            scanPlinkExact(*plink_reader_, chrom_filter, keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
        }

        auto scan_end = std::chrono::high_resolution_clock::now();
        double v_time = std::chrono::duration<double>(v_end - v_start).count();
        double scan_time = std::chrono::duration<double>(scan_end - scan_start).count();
        #pragma omp critical
        {
            scan_vbuild_time_sec_ += v_time;
            scan_snp_time_sec_ += scan_time;
            std::cout << "  " << (chrom_filter.empty() ? "Global" : "Chr " + chrom_filter) << " scan timings:"
                      << " V_build=" << std::fixed << std::setprecision(3) << v_time << "s, SNP_scan=" << scan_time << "s\n";
            appendAssocResults(out_file, chunk_results, options_.out_ssf);
        }
    };

    if (!loco_enabled) {
        std::cout << "Repeatability GWAS: running exact global scan on " << n_subjects_keep << " subjects and " << n_obs << " records.\n";
        scan_with_subject_grm(options_.mixed_precision ? G_subject_global_f32.cast<double>() : G_subject_global, "");
    } else if (use_single_pass_loco) {
        // Data already computed in compute_or_load_subject_grm during global GRM call
        std::cout << "Repeatability GWAS: running exact LOCO scan on " << chroms.size() << " chromosomes (single-pass GRM).\n";
        #pragma omp parallel for schedule(dynamic, 1) num_threads(options_.n_threads)
        for (size_t i = 0; i < chroms.size(); ++i) {
            const auto& chrom = chroms[i];
            #pragma omp critical
            { std::cout << "  Starting LOCO chromosome " << chrom << std::endl; cache_entries.push_back({"chr_" + sanitizeChromLabel(chrom), "single_pass"}); }
            double S_loco = S_global - chr_scaling[i];
            Eigen::MatrixXd G_global_d = options_.mixed_precision ? G_subject_global_f32.cast<double>() : G_subject_global;
            Eigen::MatrixXd G_loco = G_global_d;
            if (chr_scaling[i] > 0.0 && S_loco > 1e-6) {
                G_loco = (G_global_d * S_global - per_chrom_raw_subject[i].cast<double>()) / S_loco;
            }
            scan_with_subject_grm(G_loco, chrom);
        }
    } else {
        std::cout << "Repeatability GWAS: running exact LOCO scan on " << chroms.size() << " chromosomes.\n";
        std::vector<Eigen::MatrixXf> per_chrom_norm_subject;
        std::vector<double> chr_scaling2(chroms.size(), 0.0);
        std::vector<std::string> chr_cache_st2(chroms.size(), "disabled");

        #ifdef _OPENMP
        if (options_.n_threads > 1) {
            #pragma omp parallel for num_threads(options_.n_threads) schedule(dynamic, 1)
            for (size_t i = 0; i < chroms.size(); ++i) {
                double S_chr = 0.0; std::string status = "disabled";
                Eigen::MatrixXd G_chr = compute_or_load_subject_grm(chroms[i], &S_chr, nullptr, &status);
                #pragma omp critical
                { per_chrom_norm_subject.push_back(G_chr.cast<float>()); chr_scaling2[i] = S_chr; chr_cache_st2[i] = status; }
            }
        } else
        #endif
        {
            for (size_t i = 0; i < chroms.size(); ++i) {
                std::cout << "  Precomputing GRM for chromosome " << chroms[i] << std::endl;
                double S_chr = 0.0; std::string status = "disabled";
                Eigen::MatrixXd G_chr = compute_or_load_subject_grm(chroms[i], &S_chr, nullptr, &status);
                per_chrom_norm_subject.push_back(G_chr.cast<float>()); chr_scaling2[i] = S_chr; chr_cache_st2[i] = status;
            }
        }

        bool can_parallelize = (!is_bgen_ && !is_pgen_ && options_.n_threads > 1);
        #pragma omp parallel for schedule(dynamic, 1) if(can_parallelize) num_threads(options_.n_threads)
        for (size_t i = 0; i < chroms.size(); ++i) {
            const auto& chrom = chroms[i];
            #pragma omp critical
            { std::cout << "  Starting LOCO chromosome " << chrom << std::endl; cache_entries.push_back({"chr_" + sanitizeChromLabel(chrom), chr_cache_st2[i]}); }
            double S_loco = S_global - chr_scaling2[i];
            Eigen::MatrixXd G_subject_loco = G_subject_global;
            if (chr_scaling2[i] > 0.0 && S_loco > 1e-6) {
                double w1 = S_global / S_loco;
                double w2 = chr_scaling2[i] / S_loco;
                G_subject_loco = G_subject_global * w1 - per_chrom_norm_subject[i].cast<double>() * w2;
            }
            scan_with_subject_grm(G_subject_loco, chrom);
        }
    }

    { std::ofstream ofs(cache_summary_file); ofs << "Label\tStatus\n"; for (const auto& e : cache_entries) ofs << e.first << "\t" << e.second << "\n"; }
    writeNullSummary(null_file, loco_enabled ? "exact_repeatability_loco" : "exact_repeatability_global",
                     {{"grm_cache_status", global_cache_status}});
    std::cout << "Repeatability GWAS results saved to " << out_file << std::endl;
    std::cout << "Repeatability null summary saved to " << null_file << std::endl;
    std::cout << "Repeatability cache summary saved to " << cache_summary_file << std::endl;

    if (options_.refit_validation) scanRefitValidation();
}

void RepeatabilityModel::scanRefitValidation() {
    if (is_bgen_ || is_pgen_) throw std::runtime_error("Repeatability refit validation currently supports PLINK input only.");
    std::cout << "\nRunning Repeatability GWAS validation (per-SNP refit; small-data gold standard)...\n";

    const auto& bim = plink_reader_->getBimInfo();
    int n_variants = (int)bim.size();
    int n_obs = n_recs_kept_;
    int n_subjects = n_samples_;
    int variants_to_run = (options_.refit_limit > 0) ? std::min(n_variants, options_.refit_limit) : n_variants;

    std::string out_file = options_.out_prefix + ".refit.assoc.txt";
    std::string meta_file = options_.out_prefix + ".refit.meta.txt";
    writeAssocHeader(out_file, options_.out_ssf);

    { std::ofstream meta(meta_file); meta << "Item\tValue\nmethod\trepeatability_per_snp_refit_validation\nvariant_reader\tplink\n"
        << "variants_total\t" << n_variants << "\nvariants_attempted\t" << variants_to_run
        << "\nn_records\t" << n_obs << "\nn_subjects\t" << n_subjects
        << "\ninit_Vg\t" << v_g_ << "\ninit_Vpe\t" << v_pe_ << "\ninit_Ve\t" << v_e_ << "\n"; }

    Eigen::VectorXf g_raw; int written = 0, skipped_no_info = 0, skipped_solver = 0, skipped_nonpositive = 0;
    for (int j = 0; j < variants_to_run; ++j) {
        plink_reader_->getSnpGenotypes(j, g_raw);
        if (g_raw.size() == 0) { skipped_no_info++; continue; }
        double sum_g = 0.0; int valid_n = 0;
        for (int s = 0; s < n_subjects; ++s) { float val = g_raw[s]; if (!std::isnan(val)) { sum_g += val; valid_n++; } }
        if (valid_n == 0) { skipped_no_info++; continue; }
        double mean_dosage = sum_g / valid_n; double freq = sum_g / (2.0 * valid_n);
        if (freq <= 1e-8 || freq >= 1.0 - 1e-8) { skipped_no_info++; continue; }

        FixedDesignG fd_snp = fd_; const int snp_col = fd_.p; fd_snp.p = fd_.p + 1; fd_snp.rows.resize(fd_.rows.size());
        for (size_t i = 0; i < fd_.rows.size(); ++i) {
            fd_snp.rows[i] = fd_.rows[i]; int sample_idx = obs_to_sample_[i];
            double dosage = mean_dosage;
            if (sample_idx >= 0 && sample_idx < g_raw.size() && !std::isnan(g_raw[sample_idx])) dosage = g_raw[sample_idx];
            fd_snp.rows[i].push_back({snp_col, dosage});
        }

        VCEConfig vce_cfg; vce_cfg.verbose = false; vce_cfg.max_iter = options_.reml_maxit; vce_cfg.tol = options_.reml_tol; vce_cfg.vce_mode = "em";
        vce_cfg.set_algorithm_from_string(vce_cfg.vce_mode);
        AI_REML refit_solver(gen_records_, fd_snp, components_, y_obs_, vce_cfg);
        refit_solver.setInitialVars(v_e_, {v_g_, v_pe_}); refit_solver.solve();
        const auto& vars = refit_solver.getVarsU();
        if (vars.size() < 2) { skipped_solver++; continue; }
        double vg = vars[0], vpe = vars[1], ve = refit_solver.getVarE();
        if (vg <= 0.0 || vpe < 0.0 || ve <= 0.0) { skipped_nonpositive++; continue; }

        Eigen::MatrixXd X_full(n_obs, fd_snp.p); X_full.leftCols(fd_.p) = X_obs_.cast<double>();
        for (int i = 0; i < n_obs; ++i) {
            int sample_idx = obs_to_sample_[i]; double dosage = mean_dosage;
            if (sample_idx >= 0 && sample_idx < g_raw.size() && !std::isnan(g_raw[sample_idx])) dosage = g_raw[sample_idx];
            X_full(i, snp_col) = dosage;
        }

        Eigen::MatrixXd V(n_obs, n_obs);
        for (int i = 0; i < n_obs; ++i) { int si = obs_to_sample_[i]; for (int k = 0; k <= i; ++k) { int sk = obs_to_sample_[k];
            double value = vg * G_subject_(si, sk); if (si == sk) value += vpe; if (i == k) value += ve; V(i, k) = value; V(k, i) = value; } }

        Eigen::LDLT<Eigen::MatrixXd> v_solver;
        if (!factorize_with_jitter(V, v_solver)) { skipped_solver++; continue; }
        Eigen::MatrixXd VinvX = v_solver.solve(X_full); Eigen::VectorXd Vinvy = v_solver.solve(y_obs_);
        Eigen::MatrixXd XtVinvX = X_full.transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> x_solver;
        if (!factorize_with_jitter(XtVinvX, x_solver)) { skipped_solver++; continue; }
        Eigen::MatrixXd XtVinvX_inv = x_solver.solve(Eigen::MatrixXd::Identity(XtVinvX.rows(), XtVinvX.cols()));
        Eigen::VectorXd beta = XtVinvX_inv * (X_full.transpose() * Vinvy);
        double beta_j = beta(snp_col), var_beta = XtVinvX_inv(snp_col, snp_col);
        if (!(var_beta > 0.0)) { skipped_solver++; continue; }
        double se_j = std::sqrt(var_beta), z = beta_j / se_j, p = wald_pvalue_1df(z);

        GwasResult res; res.chrom = bim[j].chrom; res.snp_id = bim[j].id; res.pos = bim[j].bp_pos;
        res.a1 = bim[j].alt; res.a2 = bim[j].ref; res.freq = freq; res.n_obs = n_subjects;
        res.beta = beta_j; res.se = se_j; res.p_value = p; res.log10_p = (p > 0.0) ? -std::log10(p) : 300.0;
        appendAssocResult(out_file, res, options_.out_ssf); written++;
    }

    { std::ofstream meta(meta_file, std::ios::app); meta << "variants_written\t" << written << "\nskipped_no_info\t" << skipped_no_info
        << "\nskipped_solver\t" << skipped_solver << "\nskipped_nonpositive\t" << skipped_nonpositive << "\n"; }
    std::cout << "Repeatability refit validation results saved to " << out_file << "\n";
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
    ofs << "Item\tValue\n"
        << "n_records\t" << n_recs_kept_ << "\n"
        << "n_subjects\t" << n_samples_ << "\n"
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

void RepeatabilityModel::scanSTCG() {
    std::cout << "\nRunning Repeatability GWAS Scan (STCG PCG-based)...\n";

    int n_obs = n_recs_kept_;
    std::vector<int> keep_subject_indices = obs_to_sample_;
    std::sort(keep_subject_indices.begin(), keep_subject_indices.end());
    keep_subject_indices.erase(std::unique(keep_subject_indices.begin(), keep_subject_indices.end()), keep_subject_indices.end());
    int n_subjects_keep = (int)keep_subject_indices.size();

    std::vector<int> sample_to_subject_rank(n_samples_, -1);
    for (int i = 0; i < n_subjects_keep; ++i) sample_to_subject_rank[keep_subject_indices[i]] = i;
    std::vector<int> obs_to_subject(n_obs);
    for (int i = 0; i < n_obs; ++i) obs_to_subject[i] = sample_to_subject_rank[obs_to_sample_[i]];

    std::cout << "  Records: " << n_obs << ", Subjects: " << n_subjects_keep << "\n";
    std::cout << "  V(G)=" << v_g_ << " V(PE)=" << v_pe_ << " V(E)=" << v_e_ << "\n";

    std::vector<std::vector<int>> subj_to_rec(n_subjects_keep);
    for (int i = 0; i < n_obs; ++i) subj_to_rec[obs_to_subject[i]].push_back(i);

    std::vector<int> n_obs_per_subj(n_subjects_keep, 0);
    for (int i = 0; i < n_obs; ++i) n_obs_per_subj[obs_to_subject[i]]++;

    auto scan_stcg_with_grm = [&](const Eigen::MatrixXd& G_loco,
                                  const std::string& chrom_filter,
                                  std::ofstream& ofs) {
        int q = n_subjects_keep;
        int p = X_obs_.cols();
        double vg = v_g_, vpe = v_pe_, ve = v_e_;

        Eigen::VectorXd Zty(q);
        Eigen::MatrixXd ZtX(q, p);
        Zty.setZero(); ZtX.setZero();
        for (int i = 0; i < n_obs; ++i) {
            int s = obs_to_subject[i];
            Zty[s] += y_obs_[i];
            for (int j = 0; j < p; ++j) ZtX(s, j) += X_obs_(i, j);
        }

        Eigen::VectorXd D_diag(q);
        for (int s = 0; s < q; ++s) D_diag[s] = n_obs_per_subj[s];

        Eigen::MatrixXd K = vg * G_loco + vpe * Eigen::MatrixXd::Identity(q, q);
        Eigen::LDLT<Eigen::MatrixXd> K_ldlt(K);
        Eigen::MatrixXd Kinv_Zty = K_ldlt.solve(Zty);
        Eigen::MatrixXd Kinv_ZtX = K_ldlt.solve(ZtX);

        Eigen::MatrixXd M = K_ldlt.solve(Eigen::MatrixXd::Identity(q, q));
        for (int s = 0; s < q; ++s) M(s, s) += D_diag[s] / ve;
        Eigen::LDLT<Eigen::MatrixXd> M_ldlt(M);

        auto woodbury_Vinv_vec = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
            Eigen::VectorXd Ztv(q);
            Ztv.setZero();
            for (int i = 0; i < n_obs; ++i) Ztv[obs_to_subject[i]] += v[i];
            Eigen::VectorXd Minv_Ztv = M_ldlt.solve(Ztv);
            Eigen::VectorXd out(n_obs);
            for (int i = 0; i < n_obs; ++i) {
                out[i] = v[i] / ve - Minv_Ztv[obs_to_subject[i]] / (ve * ve);
            }
            return out;
        };

        Eigen::VectorXd Vinv_y = woodbury_Vinv_vec(y_obs_);
        Eigen::MatrixXd Vinv_X(n_obs, p);
        for (int j = 0; j < p; ++j) {
            Eigen::VectorXd X_col = X_obs_.col(j).cast<double>();
            Vinv_X.col(j) = woodbury_Vinv_vec(X_col);
        }
        Eigen::MatrixXd XtVinvX = X_obs_.cast<double>().transpose() * Vinv_X;
        Eigen::LDLT<Eigen::MatrixXd> x_solver(XtVinvX);
        Eigen::VectorXd beta;
        if (x_solver.info() == Eigen::Success) {
            Eigen::VectorXd Xty = X_obs_.cast<double>().transpose() * Vinv_y;
            beta = x_solver.solve(Xty);
        } else {
            beta = Eigen::VectorXd::Zero(p);
        }
        Eigen::VectorXd resid = y_obs_ - X_obs_.cast<double>() * beta;
        Eigen::VectorXd P_y = woodbury_Vinv_vec(resid);

        if (chrom_filter.empty()) {
            std::cout << "  Null model (Woodbury) completed (global). Beta[0]=" << (beta.size() > 0 ? beta[0] : 0.0) << "\n";
        } else {
            std::cout << "  Null model (Woodbury) completed (LOCO chr " << chrom_filter << ").\n";
        }

        auto process_snp = [&](const std::string& chrom, const std::string& rsid, uint32_t pos,
                               const std::string& ref, const std::string& alt,
                               const Eigen::VectorXf& dosage) {
            Eigen::VectorXd g(n_obs);
            double sum_g = 0.0;
            int valid_n = 0;
            for (int i = 0; i < n_obs; ++i) {
                float val = dosage[obs_to_sample_[i]];
                g[i] = val;
                if (!std::isnan(val)) { sum_g += val; valid_n++; }
                else g[i] = 0.0;
            }
            if (valid_n < 10) return;
            double maf = sum_g / (2.0 * valid_n);
            if (maf > 0.5) maf = 1.0 - maf;
            if (maf < options_.min_maf) return;
            double g_mean = sum_g / valid_n;
            g.array() -= g_mean;

            double score = g.dot(P_y);
            Eigen::VectorXd Vinv_g = woodbury_Vinv_vec(g);
            Eigen::VectorXd XtVinvX_XtVinv_g = x_solver.solve(X_obs_.cast<double>().transpose() * Vinv_g);
            Eigen::VectorXd P_g = Vinv_g - Vinv_X * XtVinvX_XtVinv_g;
            double score_var = g.dot(P_g);
            if (score_var < 1e-15) return;
            double z = score / std::sqrt(score_var);
            double pval = std::erfc(std::abs(z) / std::sqrt(2.0));
            double snp_beta = score / score_var;
            double se = 1.0 / std::sqrt(score_var);
            double log10p = (pval > 0.0 && pval <= 1.0) ? -std::log10(pval) : (pval > 1.0 ? 0.0 : 300.0);

            ofs << chrom << "\t" << rsid << "\t" << pos << "\t"
                << alt << "\t" << ref << "\t" << maf << "\t" << n_subjects_keep << "\t"
                << snp_beta << "\t" << se << "\t" << pval << "\t" << log10p << "\n";
        };

        if (is_bgen_) {
            bgen_reader_->reset();
            size_t m = bgen_reader_->getNumVariants();
            for (size_t j = 0; j < m; ++j) {
                std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
                if (!bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                process_snp(chrom, rsid, pos, ref, alt, dosage);
            }
        } else if (is_pgen_) {
            pgen_reader_->reset();
            size_t m = pgen_reader_->getNumVariants();
            for (size_t j = 0; j < m; ++j) {
                std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                process_snp(chrom, rsid, pos, ref, alt, dosage);
            }
        } else {
            const auto& bim = plink_reader_->getBimInfo();
            int m = plink_reader_->getNumSnps();
            for (int j = 0; j < m; ++j) {
                if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
                Eigen::VectorXf dosage;
                plink_reader_->getSnpGenotypes(j, dosage);
                process_snp(bim[j].chrom, bim[j].id, bim[j].bp_pos, bim[j].ref, bim[j].alt, dosage);
            }
        }
    };

    // Determine LOCO mode
    GrmBuilder::Options grm_opts; grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe; grm_opts.algorithm = options_.grm_alg;
    grm_opts.thread_num = options_.n_threads;

    bool loco_enabled = !options_.no_loco;
    std::set<std::string> observed_chroms;

    if (is_bgen_) {
        const auto& variants = bgen_reader_->getVariants();
        for (const auto& v : variants) observed_chroms.insert(v.chrom);
    } else if (is_pgen_) {
        const auto& variants = pgen_reader_->getVariants();
        for (const auto& v : variants) observed_chroms.insert(v.chrom);
    } else {
        const auto& bim = plink_reader_->getBimInfo();
        for (const auto& b : bim) observed_chroms.insert(b.chrom);
    }
    std::vector<std::string> chroms = sortChroms(observed_chroms);
    if (chroms.size() <= 1) loco_enabled = false;

    std::string out_file = options_.out_prefix + ".stcg.assoc.txt";
    std::string null_file = options_.out_prefix + ".stcg.null.summary.txt";
    std::ofstream ofs(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";

    if (!loco_enabled) {
        std::cout << "  Running global GRM scan (LOCO disabled or single chromosome).\n";
        scan_stcg_with_grm(G_subject_, "", ofs);
        writeNullSummary(null_file, "STCG",
            {{"V_G", std::to_string(v_g_)}, {"V_PE", std::to_string(v_pe_)},
             {"V_E", std::to_string(v_e_)}, {"n_records", std::to_string(n_obs)},
             {"n_subjects", std::to_string(n_subjects_keep)}, {"loco", "false"}});
    } else {
        std::cout << "  Running LOCO scan on " << chroms.size() << " chromosomes.\n";

        // Compute global GRM and per-chromosome GRMs for LOCO
        double S_global = 0.0;
        bool use_single_pass = options_.grm_cache_prefix.empty() && options_.grm_alg != 2;
        GrmBuilder::SplitGrmResult split_res;
        std::vector<Eigen::MatrixXf> per_chrom_raw;
        std::vector<double> chr_scaling;
        DataFilter::KeepInfo ki; ki.keep_indices = keep_subject_indices; ki.n_keep = keep_subject_indices.size();

        if (use_single_pass) {
            GrmBuilder split_builder(grm_opts);
            if (is_bgen_) {
                bgen_reader_->reset();
                split_res = split_builder.computeSplit(*bgen_reader_, observed_chroms, keep_subject_indices);
            } else if (is_pgen_) {
                pgen_reader_->reset();
                split_res = split_builder.computeSplit(*pgen_reader_, observed_chroms, keep_subject_indices);
            } else {
                split_res = split_builder.computeSplit(*plink_reader_, observed_chroms);
            }
            S_global = (uses_per_snp_scaling(options_.grm_alg))
                ? (double)split_res.global_snp_count : split_res.global_sum_2pq;

            per_chrom_raw.assign(chroms.size(), Eigen::MatrixXf());
            chr_scaling.assign(chroms.size(), 0.0);
            for (size_t i = 0; i < chroms.size(); ++i) {
                auto it = std::find(split_res.chrom_labels.begin(), split_res.chrom_labels.end(), chroms[i]);
                if (it != split_res.chrom_labels.end()) {
                    size_t k = std::distance(split_res.chrom_labels.begin(), it);
                    per_chrom_raw[i] = (is_bgen_ || is_pgen_)
                        ? split_res.per_chrom_raw[k].cast<float>()
                        : DataFilter::subsetSquare(split_res.per_chrom_raw[k], ki).cast<float>();
                    chr_scaling[i] = (uses_per_snp_scaling(options_.grm_alg))
                        ? (double)split_res.per_chrom_snp_count[k] : split_res.per_chrom_sum_2pq[k];
                }
            }
        } else {
            // Compute global GRM from the already-built G_subject_
            S_global = (uses_per_snp_scaling(options_.grm_alg))
                ? (double)(is_bgen_ ? bgen_reader_->getNumVariants() : is_pgen_ ? pgen_reader_->getNumVariants() : plink_reader_->getNumSnps())
                : 1.0;

            per_chrom_raw.assign(chroms.size(), Eigen::MatrixXf());
            chr_scaling.assign(chroms.size(), 0.0);
            for (size_t i = 0; i < chroms.size(); ++i) {
                GrmBuilder chr_builder(grm_opts);
                GrmBuilder::Options lo = grm_opts; lo.chrom_filter = chroms[i];
                chr_builder = GrmBuilder(lo);
                if (is_bgen_) { bgen_reader_->reset(); chr_builder.compute(*bgen_reader_); }
                else if (is_pgen_) { pgen_reader_->reset(); chr_builder.compute(*pgen_reader_); }
                else { chr_builder.compute(*plink_reader_); }
                per_chrom_raw[i] = DataFilter::subsetSquare(chr_builder.getGrm(), ki).cast<float>();
                chr_scaling[i] = (uses_per_snp_scaling(options_.grm_alg))
                    ? (double)chr_builder.getTotalSnps() : chr_builder.getSum2pq();
            }
        }

        Eigen::MatrixXd G_global = (use_single_pass)
            ? ((is_bgen_ || is_pgen_)
                ? split_res.global
                : DataFilter::subsetSquare(split_res.global, ki))
            : G_subject_;

        // Run LOCO: for each chromosome, compute LOCO GRM and scan
        for (size_t i = 0; i < chroms.size(); ++i) {
            const auto& chrom = chroms[i];
            std::cout << "  STCG LOCO chromosome " << chrom << "...\n";

            double S_loco = S_global - chr_scaling[i];
            Eigen::MatrixXd G_loco = G_global;
            if (chr_scaling[i] > 0.0 && S_loco > 1e-6) {
                G_loco = (G_global * S_global - per_chrom_raw[i].cast<double>()) / S_loco;
            }
            scan_stcg_with_grm(G_loco, chrom, ofs);
        }

        writeNullSummary(null_file, "STCG_LOCO",
            {{"V_G", std::to_string(v_g_)}, {"V_PE", std::to_string(v_pe_)},
             {"V_E", std::to_string(v_e_)}, {"n_records", std::to_string(n_obs)},
             {"n_subjects", std::to_string(n_subjects_keep)}, {"loco", "true"},
             {"n_chromosomes", std::to_string(chroms.size())}});
    }

    ofs.close();
    std::cout << "STCG repeatability GWAS results saved to " << out_file << "\n";
}

} // namespace cosmic
