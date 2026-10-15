#include "gwas_scanner.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <omp.h>

namespace cosmic {

namespace {

double biobank_gamma_q(double a, double x) {
    if (x < 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;
        for (int i = 0; i < 100; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::abs(del) < std::abs(sum) * 1e-7) {
                return 1.0 - sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
            }
        }
    } else {
        double b = x + 1.0 - a;
        double c = 1.0 / 1e-30;
        double d = 1.0 / b;
        double h = d;
        for (int i = 1; i <= 100; ++i) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::abs(d) < 1e-30) d = 1e-30;
            c = b + an / c;
            if (std::abs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (std::abs(del - 1.0) < 1e-7) {
                return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
            }
        }
    }
    return 0.0;
}

double biobank_chisq_prob(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return biobank_gamma_q(df / 2.0, stat / 2.0);
}

Eigen::MatrixXd build_biobank_covariates(const Eigen::MatrixXf& X,
                                         const Eigen::VectorXd& loco_pred,
                                         bool include_loco) {
    const int n = X.rows();
    const int p = X.cols() + (include_loco ? 1 : 0);
    Eigen::MatrixXd C(n, p);
    if (X.cols() > 0) C.leftCols(X.cols()) = X.cast<double>();
    if (include_loco) C.col(p - 1) = loco_pred;
    return C;
}

bool variant_in_scan_window(const std::string& chrom,
                            long long pos,
                            const std::string& chrom_filter,
                            long long from_bp,
                            long long to_bp) {
    if (!chrom_filter.empty() && chrom != chrom_filter) return false;
    if (from_bp >= 0 && pos < from_bp) return false;
    if (to_bp >= 0 && pos > to_bp) return false;
    return true;
}

void scan_dosage_biobank_variant(const Eigen::VectorXf& dosage,
                                 const std::vector<int>& keep_indices,
                                 const Eigen::MatrixXd& Q,
                                 const Eigen::VectorXd& y_resid,
                                 double y_resid_sq,
                                 int df,
                                 const std::string& chrom,
                                 const std::string& rsid,
                                 uint32_t pos,
                                 const std::string& ref,
                                 const std::string& alt,
                                 double min_maf,
                                 std::vector<GwasResult>& results) {
    const int n = (int)keep_indices.size();
    Eigen::VectorXd g(n);
    double sum_g = 0.0;
    int valid_n = 0;
    for (int i = 0; i < n; ++i) {
        float val = dosage[keep_indices[i]];
        g[i] = val;
        if (!std::isnan(val)) {
            sum_g += val;
            valid_n++;
        }
    }
    if (valid_n < 3) return;
    double freq = sum_g / (2.0 * valid_n);
    double maf = std::min(freq, 1.0 - freq);
    if (maf < min_maf) return;
    double mean = 2.0 * freq;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(g[i])) g[i] = mean;
        g[i] -= mean;
    }
    if (Q.cols() > 0) g -= Q * (Q.transpose() * g);
    double gg = g.squaredNorm();
    if (!(gg > 1e-12) || !std::isfinite(gg)) return;
    double gy = g.dot(y_resid);
    double beta = gy / gg;
    double rss = y_resid_sq - (gy * gy / gg);
    if (rss < 0.0 && rss > -1e-8 * y_resid_sq) rss = 0.0;
    double sigma2 = rss / (double)df;
    if (!(sigma2 > 0.0) || !std::isfinite(sigma2)) return;
    double se = std::sqrt(sigma2 / gg);
    double stat = beta / se;
    double pval = biobank_chisq_prob(stat * stat, 1);

    GwasResult res;
    res.chrom = chrom;
    res.snp_id = rsid;
    res.pos = pos;
    res.a1 = ref;
    res.a2 = alt;
    res.freq = freq;
    res.n_obs = n;
    res.beta = beta;
    res.se = se;
    res.p_value = pval;
    res.log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
    results.push_back(std::move(res));
}

} // namespace

Eigen::MatrixXd orthonormal_covariates(const Eigen::MatrixXd& C) {
    if (C.cols() == 0) return Eigen::MatrixXd(C.rows(), 0);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(C);
    int rank = qr.rank();
    if (rank <= 0) return Eigen::MatrixXd(C.rows(), 0);
    Eigen::MatrixXd Q_full = qr.householderQ() * Eigen::MatrixXd::Identity(C.rows(), rank);
    return Q_full;
}

void write_biobank_predictions(const std::string& path,
                               const PlinkReader& reader,
                               const std::vector<int>& keep_indices,
                               const std::vector<Eigen::VectorXd>& preds) {
    if (path.empty()) return;
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Cannot write prediction cache: " + path);
    ofs << "FID\tIID";
    for (int c = 0; c <= 22; ++c) ofs << "\tPRED" << c;
    ofs << "\n";
    const auto& fam = reader.getFamInfo();
    for (int i = 0; i < (int)keep_indices.size(); ++i) {
        int idx = keep_indices[i];
        ofs << fam[idx].fid << "\t" << fam[idx].iid;
        for (int c = 0; c <= 22; ++c) {
            double val = (c < (int)preds.size() && preds[c].size() == (int)keep_indices.size()) ? preds[c][i] : 0.0;
            ofs << "\t" << std::setprecision(12) << val;
        }
        ofs << "\n";
    }
}

void write_biobank_predictions(const std::string& path,
                               const std::vector<std::string>& sample_ids,
                               const std::vector<int>& keep_indices,
                               const std::vector<Eigen::VectorXd>& preds) {
    if (path.empty()) return;
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Cannot write prediction cache: " + path);
    ofs << "FID\tIID";
    for (int c = 0; c <= 22; ++c) ofs << "\tPRED" << c;
    ofs << "\n";
    for (int i = 0; i < (int)keep_indices.size(); ++i) {
        int idx = keep_indices[i];
        const std::string& sid = sample_ids[idx];
        size_t sep = sid.find(':');
        if (sep == std::string::npos) sep = sid.find('_');
        if (sep != std::string::npos) {
            ofs << sid.substr(0, sep) << "\t" << sid.substr(sep + 1);
        } else {
            ofs << sid << "\t" << sid;
        }
        for (int c = 0; c <= 22; ++c) {
            double val = (c < (int)preds.size() && preds[c].size() == (int)keep_indices.size()) ? preds[c][i] : 0.0;
            ofs << "\t" << std::setprecision(12) << val;
        }
        ofs << "\n";
    }
}

std::vector<Eigen::VectorXd> read_biobank_predictions(const std::string& path,
                                                      const PlinkReader& reader,
                                                      const std::vector<int>& keep_indices) {
    std::vector<Eigen::VectorXd> preds(23, Eigen::VectorXd::Zero((int)keep_indices.size()));
    if (path.empty()) return preds;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open prediction cache: " + path);

    std::string header;
    std::getline(ifs, header);
    std::unordered_map<std::string, std::vector<double>> cache;
    std::string fid, iid;
    while (ifs >> fid >> iid) {
        std::vector<double> vals(23, 0.0);
        for (int c = 0; c <= 22; ++c) ifs >> vals[c];
        cache[fid + ":" + iid] = std::move(vals);
    }

    const auto& fam = reader.getFamInfo();
    for (int i = 0; i < (int)keep_indices.size(); ++i) {
        int idx = keep_indices[i];
        auto it = cache.find(fam[idx].fid + ":" + fam[idx].iid);
        if (it == cache.end()) throw std::runtime_error("Prediction cache missing sample: " + fam[idx].fid + ":" + fam[idx].iid);
        for (int c = 0; c <= 22; ++c) preds[c][i] = it->second[c];
    }
    return preds;
}

std::vector<Eigen::VectorXd> read_biobank_predictions_by_sample_id(const std::string& path,
                                                                   const std::vector<std::string>& sample_ids,
                                                                   const std::vector<int>& keep_indices) {
    std::vector<Eigen::VectorXd> preds(23, Eigen::VectorXd::Zero((int)keep_indices.size()));
    if (path.empty()) return preds;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open prediction cache: " + path);

    std::string header;
    std::getline(ifs, header);
    std::unordered_map<std::string, std::vector<double>> cache_by_pair;
    std::unordered_map<std::string, std::vector<double>> cache_by_iid;
    std::string fid, iid;
    while (ifs >> fid >> iid) {
        std::vector<double> vals(23, 0.0);
        for (int c = 0; c <= 22; ++c) ifs >> vals[c];
        cache_by_iid[iid] = vals;
        cache_by_pair[fid + ":" + iid] = std::move(vals);
    }

    for (int i = 0; i < (int)keep_indices.size(); ++i) {
        int idx = keep_indices[i];
        const std::string& sid = sample_ids[idx];
        const std::vector<double>* vals = nullptr;

        auto it_iid = cache_by_iid.find(sid);
        if (it_iid != cache_by_iid.end()) vals = &it_iid->second;

        if (!vals) {
            auto it_pair = cache_by_pair.find(sid);
            if (it_pair != cache_by_pair.end()) vals = &it_pair->second;
        }

        if (!vals) {
            auto it_pair = cache_by_pair.find(sid + ":" + sid);
            if (it_pair != cache_by_pair.end()) vals = &it_pair->second;
        }

        if (!vals) {
            size_t sep = sid.find('_');
            if (sep == std::string::npos) sep = sid.find(':');
            if (sep != std::string::npos) {
                std::string sfid = sid.substr(0, sep);
                std::string siid = sid.substr(sep + 1);
                auto it_iid2 = cache_by_iid.find(siid);
                if (it_iid2 != cache_by_iid.end()) vals = &it_iid2->second;
                if (!vals) {
                    auto it_pair2 = cache_by_pair.find(sfid + ":" + siid);
                    if (it_pair2 != cache_by_pair.end()) vals = &it_pair2->second;
                }
            }
        }

        if (!vals) throw std::runtime_error("Prediction cache missing sample ID: " + sid);
        for (int c = 0; c <= 22; ++c) preds[c][i] = (*vals)[c];
    }
    return preds;
}

void scan_pgen_biobank_linear_parallel(const PgenReader& reader,
                                       uint32_t start_idx,
                                       uint32_t count,
                                       const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& X,
                                       const Eigen::VectorXd& loco_pred,
                                       bool include_loco,
                                       const std::vector<int>& keep_indices,
                                       const std::string& chrom_filter,
                                       long long from_bp,
                                       long long to_bp,
                                       double min_maf,
                                       const std::unordered_set<std::string>& extract_ids,
                                       const std::unordered_set<std::string>& exclude_ids,
                                       int n_threads,
                                       std::vector<GwasResult>& results) {
    const int n = (int)keep_indices.size();
    Eigen::MatrixXd C = build_biobank_covariates(X, loco_pred, include_loco);
    Eigen::MatrixXd Q = orthonormal_covariates(C);
    Eigen::VectorXd y_resid = y;
    if (Q.cols() > 0) y_resid -= Q * (Q.transpose() * y_resid);
    double y_resid_sq = y_resid.squaredNorm();
    const int df = std::max(1, n - (int)Q.cols() - 1);

    int threads = std::max(1, std::min(n_threads, (int)count));
    threads = std::min(threads, 32);
    results.clear();
    results.reserve(count);
    const uint32_t batch_size = 4096;
    uint32_t end_idx = start_idx + count;
    for (uint32_t batch_start = start_idx; batch_start < end_idx; batch_start += batch_size) {
        uint32_t batch_count = std::min<uint32_t>(batch_size, end_idx - batch_start);
        std::vector<std::vector<GwasResult>> thread_results(threads);

        #pragma omp parallel num_threads(threads)
        {
        int tid = omp_get_thread_num();
        int nt = omp_get_num_threads();
        uint32_t local_start = batch_start + (uint32_t)((uint64_t)batch_count * tid / nt);
        uint32_t local_end = batch_start + (uint32_t)((uint64_t)batch_count * (tid + 1) / nt);
        std::vector<GwasResult> local_results;
        local_results.reserve(local_end - local_start);

        PgenReader local_reader;
        local_reader.openSharedMetadata(reader);
        local_reader.seekVariant(local_start);
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        for (uint32_t vi = local_start; vi < local_end; vi++) {
            if (!local_reader.readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
            if (!variant_in_scan_window(chrom, pos, chrom_filter, from_bp, to_bp)) continue;
            if (!extract_ids.empty() && extract_ids.find(rsid) == extract_ids.end()) continue;
            if (!exclude_ids.empty() && exclude_ids.find(rsid) != exclude_ids.end()) continue;
            scan_dosage_biobank_variant(dosage, keep_indices, Q, y_resid, y_resid_sq, df,
                                        chrom, rsid, pos, ref, alt, min_maf, local_results);
        }
        local_reader.releaseVariantRange(local_start, local_end - local_start);
        thread_results[tid] = std::move(local_results);
        }

        size_t total = results.size();
        for (const auto& local : thread_results) total += local.size();
        results.reserve(total);
        for (auto& local : thread_results) {
            results.insert(results.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
        }
    }
}

void scan_bgen_biobank_linear_parallel(const BgenReader& reader,
                                       uint32_t start_idx,
                                       uint32_t count,
                                       const Eigen::VectorXd& y,
                                       const Eigen::MatrixXf& X,
                                       const Eigen::VectorXd& loco_pred,
                                       bool include_loco,
                                       const std::vector<int>& keep_indices,
                                       const std::string& chrom_filter,
                                       long long from_bp,
                                       long long to_bp,
                                       double min_maf,
                                       const std::unordered_set<std::string>& extract_ids,
                                       const std::unordered_set<std::string>& exclude_ids,
                                       int n_threads,
                                       std::vector<GwasResult>& results) {
    const int n = (int)keep_indices.size();
    Eigen::MatrixXd C = build_biobank_covariates(X, loco_pred, include_loco);
    Eigen::MatrixXd Q = orthonormal_covariates(C);
    Eigen::VectorXd y_resid = y;
    if (Q.cols() > 0) y_resid -= Q * (Q.transpose() * y_resid);
    double y_resid_sq = y_resid.squaredNorm();
    const int df = std::max(1, n - (int)Q.cols() - 1);

    int threads = std::max(1, std::min(n_threads, (int)count));
    results.clear();
    results.reserve(count);
    const uint32_t batch_size = 4096;
    uint32_t end_idx = start_idx + count;
    for (uint32_t batch_start = start_idx; batch_start < end_idx; batch_start += batch_size) {
        uint32_t batch_count = std::min<uint32_t>(batch_size, end_idx - batch_start);
        std::vector<std::vector<GwasResult>> thread_results(threads);

        #pragma omp parallel num_threads(threads)
        {
        int tid = omp_get_thread_num();
        int nt = omp_get_num_threads();
        uint32_t local_start = batch_start + (uint32_t)((uint64_t)batch_count * tid / nt);
        uint32_t local_end = batch_start + (uint32_t)((uint64_t)batch_count * (tid + 1) / nt);
        std::vector<GwasResult> local_results;
        local_results.reserve(local_end - local_start);

        BgenReader local_reader;
        local_reader.openSharedMetadata(reader);
        local_reader.seekVariant(local_start);
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        Eigen::VectorXf dosage;
        for (uint32_t vi = local_start; vi < local_end; vi++) {
            if (!local_reader.readVariant(chrom, rsid, pos, ref, alt, dosage)) break;
            if (!variant_in_scan_window(chrom, pos, chrom_filter, from_bp, to_bp)) continue;
            if (!extract_ids.empty() && extract_ids.find(rsid) == extract_ids.end()) continue;
            if (!exclude_ids.empty() && exclude_ids.find(rsid) != exclude_ids.end()) continue;
            scan_dosage_biobank_variant(dosage, keep_indices, Q, y_resid, y_resid_sq, df,
                                        chrom, rsid, pos, ref, alt, min_maf, local_results);
        }
        thread_results[tid] = std::move(local_results);
        }

        size_t total = results.size();
        for (const auto& local : thread_results) total += local.size();
        results.reserve(total);
        for (auto& local : thread_results) {
            results.insert(results.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
        }
    }
}

void scan_plink_biobank_linear(const PlinkReader& reader,
                               const Eigen::VectorXd& y,
                               const Eigen::MatrixXf& X,
                               const Eigen::VectorXd& loco_pred,
                               bool include_loco,
                               const std::vector<int>& keep_indices,
                               const std::string& chrom_filter,
                               long long from_bp,
                               long long to_bp,
                               double min_maf,
                               const std::unordered_set<std::string>& extract_ids,
                               const std::unordered_set<std::string>& exclude_ids,
                               int n_threads,
                               std::vector<GwasResult>& results) {
    const int n = (int)keep_indices.size();
    const int m = (int)reader.getNumSnps();
    const auto& bim = reader.getBimInfo();
    const uint8_t* raw_data = reader.getData();
    const size_t bytes_per_snp = reader.getBytesPerSnp();

    Eigen::MatrixXd C = build_biobank_covariates(X, loco_pred, include_loco);
    Eigen::MatrixXd Q = orthonormal_covariates(C);
    Eigen::VectorXd y_resid = y;
    if (Q.cols() > 0) y_resid -= Q * (Q.transpose() * y_resid);
    double y_resid_sq = y_resid.squaredNorm();

    const int df = std::max(1, n - (int)Q.cols() - 1);
    results.clear();
    results.reserve(chrom_filter.empty() ? m : m / 22 + 1024);

    int threads = std::max(1, n_threads);
    const int snp_batch_size = 4096;
    for (int batch_start = 0; batch_start < m; batch_start += snp_batch_size) {
        int batch_end = std::min(m, batch_start + snp_batch_size);
        std::vector<std::vector<GwasResult>> thread_results(threads);
        #pragma omp parallel num_threads(threads)
        {
        int tid = omp_get_thread_num();
        int nt = omp_get_num_threads();
        int j_start = batch_start + (int)((long long)(batch_end - batch_start) * tid / nt);
        int j_end = batch_start + (int)((long long)(batch_end - batch_start) * (tid + 1) / nt);
        std::vector<GwasResult> local_results;
        local_results.reserve((size_t)(j_end - j_start));
        Eigen::VectorXd g(n);
        for (int j = j_start; j < j_end; ++j) {
        if (!variant_in_scan_window(bim[j].chrom, bim[j].bp_pos, chrom_filter, from_bp, to_bp)) continue;
        if (!extract_ids.empty() && extract_ids.find(bim[j].id) == extract_ids.end()) continue;
        if (!exclude_ids.empty() && exclude_ids.find(bim[j].id) != exclude_ids.end()) continue;
        const uint8_t* snp_data = raw_data + 3 + (size_t)j * bytes_per_snp;
        double sum_g = 0.0;
        int valid_n = 0;
        for (int i = 0; i < n; ++i) {
            int sample_idx = keep_indices[i];
            uint8_t byte = snp_data[sample_idx >> 2];
            uint8_t code = (byte >> ((sample_idx & 3) * 2)) & 3;
            if (code != 1) {
                double val = (code == 0) ? 2.0 : (code == 2 ? 1.0 : 0.0);
                g[i] = val;
                sum_g += val;
                valid_n++;
            } else {
                g[i] = std::numeric_limits<double>::quiet_NaN();
            }
        }
        if (valid_n < 3) continue;
        double freq = sum_g / (2.0 * valid_n);
        double maf = std::min(freq, 1.0 - freq);
        if (maf < min_maf) continue;
        double mean = 2.0 * freq;
        for (int i = 0; i < n; ++i) {
            if (std::isnan(g[i])) g[i] = mean;
            g[i] -= mean;
        }
        if (Q.cols() > 0) g -= Q * (Q.transpose() * g);
        double gg = g.squaredNorm();
        if (!(gg > 1e-12) || !std::isfinite(gg)) continue;
        double gy = g.dot(y_resid);
        double beta = gy / gg;
        double rss = y_resid_sq - (gy * gy / gg);
        if (rss < 0.0 && rss > -1e-8 * y_resid_sq) rss = 0.0;
        double sigma2 = rss / (double)df;
        if (!(sigma2 > 0.0) || !std::isfinite(sigma2)) continue;
        double se = std::sqrt(sigma2 / gg);
        double stat = beta / se;
        double pval = biobank_chisq_prob(stat * stat, 1);

        GwasResult res;
        res.chrom = bim[j].chrom;
        res.snp_id = bim[j].id;
        res.pos = bim[j].bp_pos;
        res.a1 = bim[j].alt;
        res.a2 = bim[j].ref;
        res.freq = freq;
        res.n_obs = n;
        res.beta = beta;
        res.se = se;
        res.p_value = pval;
        res.log10_p = (pval > 0) ? -std::log10(pval) : 300.0;
        local_results.push_back(std::move(res));
        }
        thread_results[tid] = std::move(local_results);
        }
        size_t total = results.size();
        for (const auto& local : thread_results) total += local.size();
        results.reserve(total);
        for (auto& local : thread_results) {
            results.insert(results.end(), std::make_move_iterator(local.begin()), std::make_move_iterator(local.end()));
        }
        reader.releaseSnpRange((size_t)batch_start, (size_t)(batch_end - batch_start));
    }
}

} // namespace cosmic
