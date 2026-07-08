#include "repeatability_model.h"
#include "repeatability_helpers.h"

#include "data_filter.h"
#include "grm_builder.h"
#include "logger.h"
#include "matrix_io.h"
#include "qc_utils.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <omp.h>

namespace cosmic {
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

struct NullModelCache {
    Eigen::MatrixXd VinvX;
    Eigen::MatrixXd XtVinvX_inv;
    Eigen::VectorXd P_y;
    Eigen::MatrixXd Minv;
    Eigen::VectorXd Wz;
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
        if (valid_ns[c] <= 0) continue;
        double g_Vinv_g = (cache.Wz.array() * G_block.col(c).array() * VinvG.col(c).array()).sum();
        Eigen::VectorXd xvg = XtVinvG.col(c);
        double term2 = xvg.transpose() * cache.XtVinvX_inv * xvg;
        double var_score = g_Vinv_g - term2;

        GwasResult res;
        index_to_meta(c, res);
        res.freq = freqs[c];
        res.n_obs = valid_ns[c];
        if (var_score > 1e-9) {
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

template <typename ResultAppender>
static void scanPlinkExact(const PlinkReader& reader,
                           const std::string& chrom_filter,
                           double min_maf,
                           double min_hwe,
                           const std::unordered_set<std::string>& extract_set,
                           const std::unordered_set<std::string>& exclude_set,
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
            if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
            if (!extract_set.empty() && extract_set.find(bim[j].id) == extract_set.end()) continue;
            if (!exclude_set.empty() && exclude_set.find(bim[j].id) != exclude_set.end()) continue;
            snp_indices.push_back(j);
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
            long long n_hom1 = 0, n_het = 0, n_hom2 = 0;
            for (int s = 0; s < n_subjects_keep; ++s) {
                float val = g_raw[keep_subject_indices[s]];
                G_block(s, c) = val;
                if (!std::isnan(val)) {
                    sum_g += val; ++valid_n;
                    int gt = (int)std::llround(val);
                    if (gt == 0) ++n_hom2;
                    else if (gt == 1) ++n_het;
                    else if (gt == 2) ++n_hom1;
                }
            }
            freqs[c] = (valid_n > 0) ? (sum_g / (2.0 * valid_n)) : 0.0;
            valid_ns[c] = valid_n;
            double maf = std::min(freqs[c], 1.0 - freqs[c]);
            bool pass = valid_n > 0 && maf >= min_maf;
            if (pass && min_hwe > 0.0) pass = calculateHWE(n_hom1, n_het, n_hom2) >= min_hwe;
            double mean = 2.0 * freqs[c];
            for (int s = 0; s < n_subjects_keep; ++s) {
                if (!pass) G_block(s, c) = 0.0;
                else if (std::isnan(G_block(s, c))) G_block(s, c) = 0.0;
                else G_block(s, c) -= mean;
            }
            if (!pass) { freqs[c] = 0.0; valid_ns[c] = 0; }
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
                                 double min_maf,
                                 const std::unordered_set<std::string>& extract_set,
                                 const std::unordered_set<std::string>& exclude_set,
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
            if (!extract_set.empty() && extract_set.find(d.rsid) == extract_set.end()) continue;
            if (!exclude_set.empty() && exclude_set.find(d.rsid) != exclude_set.end()) continue;
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
            double maf = std::min(freqs[i], 1.0 - freqs[i]);
            bool pass = valid_n > 0 && maf >= min_maf;
            double mean = 2.0 * freqs[i];
            for (int s = 0; s < n_subjects_keep; ++s) {
                if (!pass) G_block(s, i) = 0.0;
                else if (std::isnan(G_block(s, i))) G_block(s, i) = 0.0;
                else G_block(s, i) -= mean;
            }
            if (!pass) { freqs[i] = 0.0; valid_ns[i] = 0; }
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
    grm_opts.extract_set = options_.extract_set;
    grm_opts.exclude_set = options_.exclude_set;
    grm_opts.thread_num = options_.n_threads;
    std::string global_cache_status = options_.grm_cache_prefix.empty() ? "disabled" : "not_used";

    bool use_single_pass_loco = false;
    GrmBuilder::SplitGrmResult split_res;
    std::vector<Eigen::MatrixXf> per_chrom_raw_subject;
    std::vector<double> chr_scaling;

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
            std::cerr << "Repeatability exact scan Woodbury factorization diagnostics: "
                      << data.str() << std::endl;
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

        std::vector<GwasResult> chunk_results;
        auto append_result = [&](const GwasResult& res) { chunk_results.push_back(res); };

        if (is_bgen_) {
            scanSequentialExact(*bgen_reader_, bgen_reader_->getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage);
                }, chrom_filter, options_.min_maf, options_.extract_set, options_.exclude_set,
                keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
        } else if (is_pgen_) {
            scanSequentialExact(*pgen_reader_, pgen_reader_->getNumVariants(),
                [&](std::string& chrom, std::string& rsid, uint32_t& pos, std::string& ref, std::string& alt, Eigen::VectorXf& dosage) {
                    return pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage);
                }, chrom_filter, options_.min_maf, options_.extract_set, options_.exclude_set,
                keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
        } else {
            scanPlinkExact(*plink_reader_, chrom_filter, options_.min_maf, options_.min_hwe,
                           options_.extract_set, options_.exclude_set,
                           keep_subject_indices, obs_to_subject, cache, n_subjects_keep, append_result);
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
        std::vector<Eigen::MatrixXf> per_chrom_norm_subject(chroms.size());
        std::vector<double> chr_scaling2(chroms.size(), 0.0);
        std::vector<std::string> chr_cache_st2(chroms.size(), "disabled");

        for (size_t i = 0; i < chroms.size(); ++i) {
            std::cout << "  Precomputing GRM for chromosome " << chroms[i] << std::endl;
            double S_chr = 0.0; std::string status = "disabled";
            Eigen::MatrixXd G_chr = compute_or_load_subject_grm(chroms[i], &S_chr, nullptr, &status);
            per_chrom_norm_subject[i] = G_chr.cast<float>();
            chr_scaling2[i] = S_chr;
            chr_cache_st2[i] = status;
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

} // namespace cosmic
