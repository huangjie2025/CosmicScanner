#include "repeatability_model.h"
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
#include <unordered_map>
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
    opts.mixed_precision = cfg.repeatability_mixed_precision;
    opts.refit_validation = cfg.repeatability_refit_validation;
    opts.refit_limit = cfg.repeatability_refit_limit;
    opts.grm_cache_prefix = cfg.repeatability_grm_cache_prefix;
    opts.use_covariate_projection = cfg.use_covariate_projection;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.reml_vg = cfg.reml_vg;
    opts.reml_vpe = cfg.reml_vpe;
    opts.reml_ve = cfg.reml_ve;
    opts.pheno_file = cfg.pheno_file;
    opts.mpheno = cfg.mpheno;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.use_spa = cfg.use_spa;
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

static void set_blas_threads(int n) {
    if (openblas_set_num_threads) openblas_set_num_threads(n);
}

static double wald_pvalue_1df(double z_score) {
    return std::erfc(std::abs(z_score) / std::sqrt(2.0));
}

static bool factorize_with_jitter(const Eigen::MatrixXd& matrix,
                                  Eigen::LDLT<Eigen::MatrixXd>& solver) {
    const double jitter_values[] = {0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
    for (double jitter : jitter_values) {
        Eigen::MatrixXd work = matrix;
        if (jitter > 0.0) work.diagonal().array() += jitter;
        solver.compute(work);
        if (solver.info() == Eigen::Success) {
            // Verify the solution is not degenerate
            Eigen::VectorXd diag = solver.vectorD();
            if (diag.cwiseAbs().minCoeff() > 1e-15) return true;
        }
    }
    return false;
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

static void writeAssocHeader(const std::string& filename, bool out_ssf) {
    std::ofstream ofs(filename);
    if (out_ssf) {
        ofs << "chromosome\tbase_pair_location\teffect_allele\tother_allele\teffect_allele_frequency\tsample_size\tbeta\tstandard_error\tp_value\tvariant_id\n";
    } else {
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }
}

static void appendAssocResult(const std::string& filename, const GwasResult& res, bool out_ssf) {
    std::ofstream ofs(filename, std::ios::app);
    if (std::isnan(res.p_value)) return;
    if (out_ssf) {
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

static void appendAssocResults(const std::string& filename, const std::vector<GwasResult>& results, bool out_ssf) {
    for (const auto& res : results) appendAssocResult(filename, res, out_ssf);
}

struct NullModelCache {
    Eigen::MatrixXd VinvX;
    Eigen::MatrixXd XtVinvX_inv;
    Eigen::VectorXd P_y;
    Eigen::MatrixXd H_sub;
    Eigen::VectorXd y_sub_proj;
    Eigen::MatrixXd X_sub_proj;
};

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
    int n_obs = (int)obs_to_subject.size();
    const auto& bim = reader.getBimInfo();
    int m = (int)bim.size();
    std::vector<GwasResult> chunk_results(m);
    std::vector<bool> valid(m, false);

    #pragma omp parallel
    {
        Eigen::VectorXf g_raw;
        Eigen::VectorXd subject_g(n_subjects_keep);
        #pragma omp for schedule(dynamic, 64)
        for (int j = 0; j < m; ++j) {
            if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
            reader.getSnpGenotypes(j, g_raw);
            double sum_g = 0.0; int valid_n = 0;
            for (int s = 0; s < n_subjects_keep; ++s) {
                float val = g_raw[keep_subject_indices[s]];
                subject_g[s] = val;
                if (!std::isnan(val)) { sum_g += val; valid_n++; }
            }
            double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
            double mean = 2.0 * freq;
            for (int s = 0; s < n_subjects_keep; ++s) {
                if (std::isnan(subject_g[s])) subject_g[s] = 0.0;
                else subject_g[s] -= mean;
            }
            Eigen::VectorXd v = cache.H_sub * subject_g;
            double g_Vinv_g = subject_g.dot(v);
            double score = subject_g.dot(cache.y_sub_proj);
            Eigen::VectorXd XtVinv_g = cache.X_sub_proj.transpose() * subject_g;
            double term2 = XtVinv_g.transpose() * cache.XtVinvX_inv * XtVinv_g;
            double var_score = g_Vinv_g - term2;

            GwasResult res;
            res.chrom = bim[j].chrom; res.snp_id = bim[j].id; res.pos = bim[j].bp_pos;
            res.a1 = bim[j].alt; res.a2 = bim[j].ref; res.freq = freq; res.n_obs = n_subjects_keep;
            if (var_score > 1e-9) {
                double se = std::sqrt(var_score); double stat = score / se;
                res.beta = score / var_score; res.se = 1.0 / se;
                res.p_value = chisq_prob(stat * stat, 1);
                res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
            }
            chunk_results[j] = res; valid[j] = true;
        }
    }
    for (int j = 0; j < m; ++j) if (valid[j]) append_result(chunk_results[j]);
}

template <typename Reader, typename ReadFn, typename ResultAppender>
static void scanSequentialExact(Reader& reader, size_t m, ReadFn read_variant,
                                const std::string& chrom_filter,
                                const std::vector<int>& keep_subject_indices,
                                const std::vector<int>& obs_to_subject,
                                const NullModelCache& cache,
                                int n_subjects_keep,
                                ResultAppender append_result) {
    int n_obs = (int)obs_to_subject.size();
    reader.reset();
    const size_t CHUNK_SIZE = 1024;
    struct SnpData {
        std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
    };
    std::vector<SnpData> chunk_data; chunk_data.reserve(CHUNK_SIZE);
    size_t variants_read = 0;
    while (variants_read < m) {
        chunk_data.clear();
        for (size_t i = 0; i < CHUNK_SIZE && variants_read < m; ++i) {
            SnpData d; bool ok = read_variant(d.chrom, d.rsid, d.pos, d.ref, d.alt, d.dosage);
            if (!ok) break; variants_read++;
            if (!chrom_filter.empty() && d.chrom != chrom_filter) continue;
            chunk_data.push_back(std::move(d));
        }
        if (chunk_data.empty()) break;
        int csz = (int)chunk_data.size();
        std::vector<GwasResult> chunk_results(csz);
        std::vector<bool> valid(csz, false);
        #pragma omp parallel
        {
            Eigen::VectorXd subject_g(n_subjects_keep);
            #pragma omp for schedule(dynamic, 64)
            for (int i = 0; i < csz; ++i) {
                const auto& d = chunk_data[i];
                double sum_g = 0.0; int valid_n = 0;
                for (int s = 0; s < n_subjects_keep; ++s) {
                    float val = d.dosage[keep_subject_indices[s]];
                    subject_g[s] = val;
                    if (!std::isnan(val)) { sum_g += val; valid_n++; }
                }
                double freq = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
                double mean = 2.0 * freq;
                for (int s = 0; s < n_subjects_keep; ++s) {
                    if (std::isnan(subject_g[s])) subject_g[s] = 0.0;
                    else subject_g[s] -= mean;
                }
                Eigen::VectorXd v = cache.H_sub * subject_g;
                double g_Vinv_g = subject_g.dot(v);
                double score = subject_g.dot(cache.y_sub_proj);
                Eigen::VectorXd XtVinv_g = cache.X_sub_proj.transpose() * subject_g;
                double term2 = XtVinv_g.transpose() * cache.XtVinvX_inv * XtVinv_g;
                double var_score = g_Vinv_g - term2;
                GwasResult res;
                res.chrom = d.chrom; res.snp_id = d.rsid; res.pos = d.pos;
                res.a1 = d.alt; res.a2 = d.ref; res.freq = freq; res.n_obs = n_subjects_keep;
                if (var_score > 1e-9) {
                    double se = std::sqrt(var_score); double stat = score / se;
                    res.beta = score / var_score; res.se = 1.0 / se;
                    res.p_value = chisq_prob(stat * stat, 1);
                    res.log10_p = (res.p_value > 0.0) ? -std::log10(res.p_value) : 300.0;
                }
                chunk_results[i] = res; valid[i] = true;
            }
        }
        for (int i = 0; i < csz; ++i) if (valid[i]) append_result(chunk_results[i]);
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
    buildGRM();
    buildMMEComponents();
    runVCE();
}

void RepeatabilityModel::loadPhenotype() {
    if (options_.pheno_file.empty()) throw std::runtime_error("Repeatability model requires --pheno file (long format).");
    std::cout << "Loading long-format phenotype data...\n";
    PhenoReader pheno_reader;
    auto records = pheno_reader.loadAllRecords(options_.pheno_file, options_.mpheno);

    std::unordered_map<std::string, int> fid_iid_to_gen_idx;
    for (int i = 0; i < n_samples_; ++i) fid_iid_to_gen_idx[fids_[i] + ":" + iids_[i]] = i;

    gen_records_.clear();
    for (const auto& rec : records) {
        std::string key = rec.fid + ":" + rec.iid;
        auto it = fid_iid_to_gen_idx.find(key);
        if (it != fid_iid_to_gen_idx.end()) {
            GenRecord gr;
            gr.aid = it->second + 1;
            gr.y = rec.value;
            gr.idstr = key;
            gen_records_.push_back(gr);
        }
    }
    n_recs_kept_ = (int)gen_records_.size();
    std::cout << "Kept " << n_recs_kept_ << " phenotype records matching genotype samples.\n";
    if (n_recs_kept_ < 10) throw std::runtime_error("Too few records.");
}

void RepeatabilityModel::loadCovariates() {
    Eigen::MatrixXf covs_q, covs_c;
    bool has_q = false, has_c = false;
    if (!options_.qcovar_file.empty()) {
        CovarReader cov_reader; cov_reader.load(options_.qcovar_file);
        cov_reader.align(fids_, iids_, covs_q); has_q = true;
    }
    if (!options_.covar_file.empty()) {
        CovarReader cov_reader; cov_reader.loadCategorical(options_.covar_file);
        cov_reader.align(fids_, iids_, covs_c); has_c = true;
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
                float val = covs_q(gen_idx, k);
                if (!std::isnan(val)) fd_.rows[i].push_back({current_col + k, (double)val});
            }
            current_col += covs_q.cols();
        }
        if (has_c) {
            for (int k = 0; k < covs_c.cols(); ++k) {
                float val = covs_c(gen_idx, k);
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
    bool use_stcg_vce = (options_.repeatability_mode == RepeatabilityMode::STCG ||
                         (options_.repeatability_mode == RepeatabilityMode::Auto && n_recs_kept_ > 80000))
                        && components_[0].geno_mat != nullptr;

    if (use_stcg_vce) {
        // STCG_VCE: matrix-free VCE using stochastic trace + PCG
        std::cout << "  Using STCG_VCE (matrix-free) for null model...\n";
        double total_var_obs = y_obs_.squaredNorm() / n_recs_kept_;

        double vg_init = std::max(total_var_obs * 0.5, 1e-4);
        double vpe_init = std::max(total_var_obs * 0.3, 1e-4);
        double ve_init = std::max(total_var_obs * 0.2, 1e-4);

        VCEConfig stcg_cfg;
        stcg_cfg.verbose = true;
        stcg_cfg.max_iter = options_.reml_maxit;
        stcg_cfg.tol = options_.reml_tol;
        stcg_cfg.vce_mode = "stcg";
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
    double total_var_obs = y_obs_.squaredNorm() / n_recs_kept_;

    auto run_single_vce = [&](double vg_frac, double vpe_frac) -> std::pair<int, double> {
        double vg_init = std::max(total_var_obs * vg_frac, 1e-4);
        double vpe_init = std::max(total_var_obs * vpe_frac, 1e-4);
        double ve_init  = std::max(total_var_obs * (1.0 - vg_frac - vpe_frac), 1e-4);
        VCEConfig vce_cfg; vce_cfg.verbose = true; vce_cfg.max_iter = options_.reml_maxit;
        vce_cfg.tol = options_.reml_tol; vce_cfg.vce_mode = "em";
        vce_cfg.print_report = false;  // Suppress report for multi-start trials
        AI_REML solver(gen_records_, fd_, components_, y_obs_, vce_cfg);
        solver.setInitialVariances({vg_init, vpe_init}, ve_init);
        solver.solve();
        v_e_ = solver.getVarE();
        auto vars = solver.getVarsU();
        if (vars.size() < 2) return {-1, -std::numeric_limits<double>::infinity()};
        v_g_ = vars[0]; v_pe_ = vars[1];
        double logL = solver.getLogLikelihood();
        double total = v_g_ + v_pe_ + v_e_;
        // 0 = acceptable, 1 = Vg collapsed, 2 = Ve squeezed (Vg+Vpe regime)
        int status = 0;
        if (v_g_ / total < 0.001) status = 1;
        else if (v_e_ / total < 0.1) status = 2;
        return {status, logL};
    };

    // Multi-start VCE: try multiple initializations, pick the one with highest logL
    struct TrialResult { double vg, vpe, ve, logL; int status; };
    std::vector<TrialResult> trials;
    std::vector<std::pair<double,double>> inits = {
        {0.5, 0.3},   // Standard split
        {0.7, 0.2},   // High genetic init
        {0.2, 0.2},   // High residual init
        {0.3, 0.5},   // High PE init
        {0.4, 0.1},   // Low PE init
    };

    for (auto& [vg_f, vpe_f] : inits) {
        auto [status, logL] = run_single_vce(vg_f, vpe_f);
        trials.push_back({v_g_, v_pe_, v_e_, logL, status});
    }

    // Select best: prefer status==0, then highest logL
    int best_idx = 0;
    double best_logL = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < (int)trials.size(); ++i) {
        auto& t = trials[i];
        if (t.status == 0 && t.logL > best_logL) {
            best_logL = t.logL;
            best_idx = i;
        }
    }
    // If no status==0 found, pick highest logL regardless
    if (best_logL == -std::numeric_limits<double>::infinity()) {
        for (int i = 0; i < (int)trials.size(); ++i) {
            if (trials[i].logL > best_logL) {
                best_logL = trials[i].logL;
                best_idx = i;
            }
        }
    }

    v_g_ = trials[best_idx].vg;
    v_pe_ = trials[best_idx].vpe;
    v_e_ = trials[best_idx].ve;

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

    writeNullSummary(options_.out_prefix + ".repeatability.null.txt", "fast_repeatability_mean_model",
                     {{"mean_records", std::to_string(n_mean)},
                      {"min_records", std::to_string(n_min)},
                      {"max_records", std::to_string(n_max)}});

    Eigen::VectorXd y_w = y_agg; Eigen::MatrixXf X_w = X_agg;
    for (int i = 0; i < n_agg; ++i) { double w = obs_weights(i); y_w(i) *= w; X_w.row(i) *= (float)w; }

    LmmEngine::Options lmm_fast_opts;
    lmm_fast_opts.fast_mode = true; lmm_fast_opts.min_maf = options_.min_maf;
    lmm_fast_opts.min_hwe = options_.min_hwe; lmm_fast_opts.use_covariate_projection = options_.use_covariate_projection;
    lmm_fast_opts.use_spa = options_.use_spa; lmm_fast_opts.thread_num = options_.n_threads;

    std::string assoc_file = options_.out_prefix + ".assoc.txt";
    writeAssocHeader(assoc_file, options_.out_ssf);

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
            EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_loco_sub); vce_fast.runNullModel();
            LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

            std::vector<GwasResult> chrom_results;
            if (is_bgen_) { bgen_reader_->reset(); engine_fast.scan(*bgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else if (is_pgen_) { pgen_reader_->reset(); engine_fast.scan(*pgen_reader_, chrom_results, agg_keep.keep_indices, chrom); }
            else { engine_fast.scan(*plink_reader_, chrom_results, agg_keep.keep_indices, chrom); }

            appendAssocResults(assoc_file, chrom_results, options_.out_ssf);
        }
    } else {
        std::cout << "Fast repeatability GWAS: global scan on " << n_agg << " aggregated subjects.\n";
        Eigen::MatrixXd G_fast_sub = DataFilter::subsetSquare(G_subject_, agg_keep);
        G_fast_sub.array().colwise() *= obs_weights.array();
        G_fast_sub.array().rowwise() *= obs_weights.transpose().array();

        EigenVce::Options vce_fast_opts; vce_fast_opts.use_covariate_projection = options_.use_covariate_projection;
        EigenVce vce_fast(vce_fast_opts); vce_fast.prepare(y_w, X_w, G_fast_sub); vce_fast.runNullModel();
        LmmEngine engine_fast(lmm_fast_opts); engine_fast.prepareFromVce(vce_fast);

        std::vector<GwasResult> results_fast;
        if (is_bgen_) { bgen_reader_->reset(); engine_fast.scan(*bgen_reader_, results_fast, agg_keep.keep_indices); }
        else if (is_pgen_) { pgen_reader_->reset(); engine_fast.scan(*pgen_reader_, results_fast, agg_keep.keep_indices); }
        else { engine_fast.scan(*plink_reader_, results_fast, agg_keep.keep_indices); }
        appendAssocResults(assoc_file, results_fast, options_.out_ssf);
    }
    std::cout << "Fast repeatability GWAS results saved to " << assoc_file << std::endl;
}

void RepeatabilityModel::scanExact() {
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
                               && !options_.no_loco && options_.grm_alg != 2;

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
        Eigen::VectorXd Dz = Eigen::VectorXd::Zero(n_subjects_keep);
        for (int i = 0; i < n_obs; ++i) Dz[obs_to_subject[i]] += 1.0;

        Eigen::MatrixXd H_inv = K;
        for (int i = 0; i < n_subjects_keep; ++i) if (Dz[i] > 0) H_inv(i, i) += v_e_ / Dz[i];
        Eigen::LLT<Eigen::MatrixXd> llt(H_inv);
        Eigen::MatrixXd H_sub = llt.solve(Eigen::MatrixXd::Identity(n_subjects_keep, n_subjects_keep));

        Eigen::VectorXd y_sub = Eigen::VectorXd::Zero(n_subjects_keep);
        Eigen::MatrixXd X_sub = Eigen::MatrixXd::Zero(n_subjects_keep, X_obs_.cols());
        for (int i = 0; i < n_obs; ++i) { int subj = obs_to_subject[i]; y_sub[subj] += y_obs_[i]; X_sub.row(subj) += X_obs_.row(i).cast<double>(); }

        Eigen::VectorXd u_y = H_sub * (K * y_sub);
        for (int i = 0; i < n_subjects_keep; ++i) if (Dz[i] > 0) u_y[i] /= Dz[i];
        Eigen::MatrixXd U_X = H_sub * (K * X_sub);
        for (int i = 0; i < n_subjects_keep; ++i) if (Dz[i] > 0) U_X.row(i) /= Dz[i];

        Eigen::VectorXd Vinv_y(n_obs); Eigen::MatrixXd VinvX(n_obs, X_obs_.cols());
        for (int i = 0; i < n_obs; ++i) {
            int subj = obs_to_subject[i];
            Vinv_y[i] = (y_obs_[i] - u_y[subj]) / v_e_;
            VinvX.row(i) = (X_obs_.row(i).cast<double>() - U_X.row(subj)) / v_e_;
        }
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

        Eigen::VectorXd y_sub_bar = Eigen::VectorXd::Zero(n_subjects_keep);
        Eigen::MatrixXd X_sub_bar = Eigen::MatrixXd::Zero(n_subjects_keep, X_obs_.cols());
        for (int i = 0; i < n_subjects_keep; ++i) if (Dz[i] > 0) { y_sub_bar[i] = y_sub[i] / Dz[i]; X_sub_bar.row(i) = X_sub.row(i) / Dz[i]; }

        NullModelCache cache; cache.VinvX = VinvX; cache.XtVinvX_inv = XtVinvX_inv; cache.P_y = P_y;
        cache.H_sub = H_sub; cache.y_sub_proj = H_sub * y_sub_bar - H_sub * X_sub_bar * beta_null;
        cache.X_sub_proj = H_sub * X_sub_bar;

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
    double total_var = v_g_ + v_pe_ + v_e_;
    ofs << "Item\tValue\n"
        << "n_records\t" << n_recs_kept_ << "\n"
        << "n_subjects\t" << n_samples_ << "\n"
        << "Vg\t" << v_g_ << "\nVpe\t" << v_pe_ << "\nVe\t" << v_e_ << "\nVp\t" << total_var << "\n"
        << "h2\t" << ((total_var > 0.0) ? (v_g_ / total_var) : 0.0) << "\n"
        << "r2\t" << ((total_var > 0.0) ? ((v_g_ + v_pe_) / total_var) : 0.0) << "\n"
        << "kernel_precision\t" << (options_.mixed_precision ? "float32_mixed" : "float64_exact") << "\n"
        << "grm_cache_prefix\t" << (options_.grm_cache_prefix.empty() ? "NA" : options_.grm_cache_prefix) << "\n"
        << "grm_cache_status\t"
        << (options_.grm_cache_prefix.empty() ? "disabled" : "not_used") << "\n"
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

    // Build PE grouping: for each subject, list of record indices
    std::vector<std::vector<int>> subj_to_rec(n_subjects_keep);
    for (int i = 0; i < n_obs; ++i) subj_to_rec[obs_to_subject[i]].push_back(i);

    // Build V operator: V = vg*ZGZ' + vpe*ZZ' + ve*I
    // For standardized GRM: ZGZ' has diagonal ≈ vg, off-diagonal = vg*G[i,j]
    // ZZ' has 1 for same-subject records, 0 otherwise
    // We use the GRM (G_subject_) for the genetic component
    auto apply_V = [&](const Eigen::VectorXd& v, Eigen::VectorXd& out) {
        out.setZero(n_obs);
        // ve * v
        out.noalias() += v_e_ * v;
        // vpe * ZZ' * v: for each subject, scatter sum to all records
        for (int s = 0; s < n_subjects_keep; ++s) {
            const auto& recs = subj_to_rec[s];
            double sum_v = 0.0;
            for (int r : recs) sum_v += v[r];
            for (int r : recs) out[r] += v_pe_ * sum_v;
        }
        // vg * ZGZ' * v: genetic component using GRM
        // First gather v from records to subjects (mean per subject)
        // Then apply G, then scatter back
        if (G_subject_.size() > 0) {
            // Gather: v_subj[s] = sum of v over records for subject s
            Eigen::VectorXd v_subj = Eigen::VectorXd::Zero(n_subjects_keep);
            for (int s = 0; s < n_subjects_keep; ++s) {
                for (int r : subj_to_rec[s]) v_subj[s] += v[r];
            }
            // Apply G: Gv = G * v_subj
            Eigen::VectorXd Gv = G_subject_ * v_subj;
            // Scatter back: for each record of subject s, add vg * Gv[s]
            for (int s = 0; s < n_subjects_keep; ++s) {
                for (int r : subj_to_rec[s]) out[r] += v_g_ * Gv[s];
            }
        } else {
            // No GRM available: use diagonal approximation (G_ii ≈ 1)
            for (int s = 0; s < n_subjects_keep; ++s) {
                for (int r : subj_to_rec[s]) out[r] += v_g_ * v[r];
            }
        }
    };

    // Build Jacobi preconditioner: V_diag^{-1}
    // V_ii = vg + vpe * n_obs_i + ve
    std::vector<int> n_obs_per_subj(n_subjects_keep, 0);
    for (int i = 0; i < n_obs; ++i) n_obs_per_subj[obs_to_subject[i]]++;

    Eigen::VectorXd V_diag_inv(n_obs);
    for (int i = 0; i < n_obs; ++i) {
        double V_ii = v_g_ + v_pe_ * n_obs_per_subj[obs_to_subject[i]] + v_e_;
        V_diag_inv[i] = 1.0 / std::max(V_ii, 1e-10);
    }

    auto apply_precond = [&](const Eigen::VectorXd& r, Eigen::VectorXd& z) {
        z = r.cwiseProduct(V_diag_inv);
    };

    // Compute null model: V^{-1} * y and V^{-1} * X using PCG
    MatrixFreePCG pcg(n_obs);
    pcg.setMultOp(apply_V);
    pcg.setPrecondOp(apply_precond);
    pcg.setQuiet(true);

    // Compute V^{-1} * y
    Eigen::VectorXd Vinv_y = pcg.solve(y_obs_, 1e-8, 500);

    // Compute V^{-1} * X columns for projection
    Eigen::MatrixXd Vinv_X(n_obs, X_obs_.cols());
    for (int j = 0; j < X_obs_.cols(); ++j) {
        Eigen::VectorXd X_col = X_obs_.col(j).cast<double>();
        Vinv_X.col(j) = pcg.solve(X_col, 1e-8, 500);
    }

    // Compute beta = (X'V^{-1}X)^{-1} X'V^{-1}y
    Eigen::MatrixXd XtVinvX = X_obs_.cast<double>().transpose() * Vinv_X;
    Eigen::LDLT<Eigen::MatrixXd> x_solver(XtVinvX);
    Eigen::VectorXd beta;
    if (x_solver.info() == Eigen::Success) {
        Eigen::VectorXd Xty = X_obs_.cast<double>().transpose() * Vinv_y;
        beta = x_solver.solve(Xty);
    } else {
        beta = Eigen::VectorXd::Zero(X_obs_.cols());
    }

    // Compute P * y = V^{-1} * (y - X*beta)
    Eigen::VectorXd resid = y_obs_ - X_obs_.cast<double>() * beta;
    Eigen::VectorXd P_y = pcg.solve(resid, 1e-8, 500);

    std::cout << "  Null model PCG completed. Beta[0]=" << (beta.size() > 0 ? beta[0] : 0.0) << "\n";

    // Open output file
    std::string out_file = options_.out_prefix + ".stcg.assoc.txt";
    std::ofstream ofs(out_file);
    ofs << "SNP\tChrom\tBP\tA1\tA2\tMAF\tScore\tScoreVar\tZ\tP\n";

    // Expand genotype from subject-level to record-level and compute score
    auto process_snp = [&](const std::string& chrom, const std::string& rsid, uint32_t pos,
                           const std::string& ref, const std::string& alt,
                           const Eigen::VectorXf& dosage) {
        // Build record-level genotype
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

        // Center genotype
        double g_mean = sum_g / valid_n;
        g.array() -= g_mean;

        // Score = g' * P * y
        double score = g.dot(P_y);

        // Score variance = g' * P * g via PCG
        // P * g = V^{-1} * (g - X * (X'V^{-1}X)^{-1} X'V^{-1}g)
        Eigen::VectorXd Vinv_g = pcg.solve(g, 1e-8, 200);
        Eigen::VectorXd XtVinvX_XtVinv_g = x_solver.solve(X_obs_.cast<double>().transpose() * Vinv_g);
        Eigen::VectorXd P_g = Vinv_g - Vinv_X * XtVinvX_XtVinv_g;
        double score_var = g.dot(P_g);

        if (score_var < 1e-15) return;

        double z = score / std::sqrt(score_var);
        double pval = 2.0 * std::erfc(std::abs(z) / std::sqrt(2.0));

        ofs << rsid << "\t" << chrom << "\t" << pos << "\t"
            << alt << "\t" << ref << "\t" << maf << "\t"
            << score << "\t" << score_var << "\t" << z << "\t" << pval << "\n";
    };

    // Sequential scan using readVariant for all formats
    if (is_bgen_) {
        bgen_reader_->reset();
        size_t m = bgen_reader_->getNumVariants();
        for (size_t j = 0; j < m; ++j) {
            std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
            if (!bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
            process_snp(chrom, rsid, pos, ref, alt, dosage);
        }
    } else if (is_pgen_) {
        pgen_reader_->reset();
        size_t m = pgen_reader_->getNumVariants();
        for (size_t j = 0; j < m; ++j) {
            std::string chrom, rsid, ref, alt; uint32_t pos; Eigen::VectorXf dosage;
            if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
            process_snp(chrom, rsid, pos, ref, alt, dosage);
        }
    } else {
        // PLINK format
        const auto& bim = plink_reader_->getBimInfo();
        int m = plink_reader_->getNumSnps();
        for (int j = 0; j < m; ++j) {
            Eigen::VectorXf dosage;
            plink_reader_->getSnpGenotypes(j, dosage);
            process_snp(bim[j].chrom, bim[j].id, bim[j].bp_pos, bim[j].ref, bim[j].alt, dosage);
        }
    }

    ofs.close();
    std::cout << "STCG repeatability GWAS results saved to " << out_file << "\n";

    writeNullSummary(options_.out_prefix + ".stcg.null.summary.txt", "STCG",
        {{"V_G", std::to_string(v_g_)}, {"V_PE", std::to_string(v_pe_)},
         {"V_E", std::to_string(v_e_)}, {"n_records", std::to_string(n_obs)},
         {"n_subjects", std::to_string(n_subjects_keep)}});
}

} // namespace cosmic
