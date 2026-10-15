#include "gemma_summary_tools.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "grm_builder.h"

namespace cosmic {
namespace {

std::unordered_set<std::string> load_ids(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (in >> id) ids.insert(id);
    return ids;
}

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::string normalized = line;
    for (char& c : normalized) {
        if (c == ',') c = ' ';
    }
    std::stringstream ss(normalized);
    std::vector<std::string> out;
    std::string x;
    while (ss >> x) out.push_back(x);
    return out;
}

int find_col(const std::unordered_map<std::string, int>& cols,
             const std::vector<std::string>& names) {
    for (const auto& n : names) {
        auto it = cols.find(upper(n));
        if (it != cols.end()) return it->second;
    }
    return -1;
}

std::vector<double> read_numbers(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open numeric file: " + path);
    std::vector<double> values;
    double x;
    while (in >> x) values.push_back(x);
    if (values.empty()) throw std::runtime_error("Empty numeric file: " + path);
    return values;
}

std::vector<std::string> read_prefix_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open prefix list: " + path);
    std::vector<std::string> prefixes;
    std::string p;
    while (in >> p) prefixes.push_back(p);
    if (prefixes.empty()) throw std::runtime_error("Empty prefix list: " + path);
    return prefixes;
}

double chisq1_from_p(double p) {
    if (!(p > 0.0) || p >= 1.0) return 0.0;
    double lo = 0.0;
    double hi = 1.0;
    auto surv = [](double x) {
        return std::erfc(std::sqrt(std::max(0.0, x) / 2.0));
    };
    while (surv(hi) > p && hi < 1e6) hi *= 2.0;
    for (int i = 0; i < 100; ++i) {
        double mid = 0.5 * (lo + hi);
        if (surv(mid) > p) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

struct SummaryQ {
    double q = 0.0;
    double vq = 0.0;
    double snp_weight_sum = 0.0;
    int ni_total = 0;
    int usable = 0;
    int total = 0;
};

struct CategorySpec {
    int ncat = 1;
    std::unordered_map<std::string, int> cat;
    std::unordered_map<std::string, std::vector<double>> wcat;
    std::unordered_map<std::string, double> wsnp;
};

struct SummaryQMulti {
    Eigen::VectorXd q;
    Eigen::MatrixXd vq;
    Eigen::VectorXd s;
    int ni_total = 0;
    int usable = 0;
    int total = 0;
};

struct BetaRecord {
    std::string rs;
    int cat = 0;
    int n = 0;
    double weight = 1.0;
    double z2 = 0.0;
};

struct CorRecord {
    std::string rs;
    int n = 0;
    double bp = 0.0;
    double cm = 0.0;
    double var = 0.0;
    std::vector<double> cor;
    bool in_beta = false;
    int beta_n = 0;
    int cat = 0;
    double weight = 1.0;
};

struct CorSummary {
    std::vector<CorRecord> records;
    std::unordered_map<std::string, size_t> index;
};

SummaryQ compute_summary_q(const GwasConfig& cfg) {
    std::ifstream in(cfg.beta_file);
    if (!in) throw std::runtime_error("Cannot open beta file: " + cfg.beta_file);
    std::string header_line;
    if (!std::getline(in, header_line)) throw std::runtime_error("Empty beta file: " + cfg.beta_file);
    auto header = split_ws(header_line);
    std::unordered_map<std::string, int> cols;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;

    const int z_col = find_col(cols, {"Z", "ZSCORE", "Z_SCORE"});
    const int beta_col = find_col(cols, {"BETA", "B"});
    const int se_col = find_col(cols, {"SE", "SE_BETA", "SEBETA"});
    const int chisq_col = find_col(cols, {"CHISQ", "CHISQ_STAT", "CHISQUARE", "WALD", "LRT", "SCORE"});
    const int p_col = find_col(cols, {"P", "PVALUE", "P_VALUE", "P-VALUE"});
    const int n_col = find_col(cols, {"N", "N_TOTAL", "NTOTAL"});
    const int nmis_col = find_col(cols, {"N_MIS", "NMIS", "N_MISS"});
    const int nobs_col = find_col(cols, {"N_OBS", "NOBS"});
    if (z_col < 0 && (beta_col < 0 || se_col < 0) && chisq_col < 0 && p_col < 0) {
        throw std::runtime_error("beta file needs Z, or BETA+SE, or CHISQ, or P columns");
    }
    if (n_col < 0 && (nmis_col < 0 || nobs_col < 0)) {
        throw std::runtime_error("beta file needs N, or N_MIS+N_OBS columns");
    }

    std::vector<double> z2;
    std::vector<int> ns;
    std::string line;
    SummaryQ res;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < header.size()) continue;
        ++res.total;
        double val = 0.0;
        if (z_col >= 0) {
            double z = std::stod(f[z_col]);
            val = z * z;
        } else if (beta_col >= 0 && se_col >= 0) {
            double se = std::stod(f[se_col]);
            if (se == 0.0) continue;
            double z = std::stod(f[beta_col]) / se;
            val = z * z;
        } else if (chisq_col >= 0) {
            val = std::stod(f[chisq_col]);
        } else if (p_col >= 0) {
            val = chisq1_from_p(std::stod(f[p_col]));
        }
        int n = 0;
        if (n_col >= 0) n = std::stoi(f[n_col]);
        else n = std::stoi(f[nmis_col]) + std::stoi(f[nobs_col]);
        if (n <= 0 || !std::isfinite(val)) continue;
        z2.push_back(val);
        ns.push_back(n);
        res.ni_total = std::max(res.ni_total, n);
    }
    if (z2.empty()) throw std::runtime_error("found no usable summary-stat rows");

    double vec_q = 0.0;
    double vec_s = 0.0;
    for (size_t i = 0; i < z2.size(); ++i) {
        vec_q += (z2[i] - 1.0) / static_cast<double>(ns[i]);
        vec_s += 1.0;
    }
    res.q = vec_q / vec_s;
    res.snp_weight_sum = vec_s;
    res.usable = static_cast<int>(z2.size());

    const int n_block = std::max(1, cfg.summary_block_count);
    std::vector<double> block_q(n_block, 0.0), block_s(n_block, 0.0);
    const int per_block = std::max(1, static_cast<int>(std::floor(z2.size() / static_cast<double>(n_block))));
    int b = 0, t = 0;
    for (size_t i = 0; i < z2.size(); ++i) {
        block_q[b] += (z2[i] - 1.0);
        block_s[b] += 1.0;
        if (b < n_block - 1) {
            if (t < per_block - 1) ++t;
            else { ++b; t = 0; }
        }
    }
    std::vector<double> jack;
    for (int k = 0; k < n_block; ++k) {
        if (block_s[k] != 0.0 && vec_s != block_s[k]) {
            jack.push_back((vec_q - block_q[k]) / (vec_s - block_s[k]));
        }
    }
    double mean = 0.0;
    for (double v : jack) mean += v;
    if (!jack.empty()) mean /= static_cast<double>(jack.size());
    for (double& v : jack) {
        v -= mean;
        res.vq += v * v;
    }
    if (jack.size() > 1) res.vq = res.vq / static_cast<double>(jack.size()) * static_cast<double>(jack.size() - 1);
    return res;
}

bool is_marker_col(const std::string& name) {
    const std::string u = upper(name);
    return u == "RS" || u == "SNP" || u == "ID" || u == "CHR" || u == "CHROM" ||
           u == "PS" || u == "POS" || u == "BP" || u == "POSITION" || u == "CM" ||
           u == "A1" || u == "ALLELE1" || u == "A0" || u == "A2" || u == "ALLELE0";
}

CategorySpec read_category_spec(const GwasConfig& cfg) {
    CategorySpec spec;
    auto read_one_cat = [&](const std::string& category_file) {
        std::ifstream in(category_file);
        if (!in) throw std::runtime_error("Cannot open category file: " + category_file);
        std::string header_line;
        if (!std::getline(in, header_line)) throw std::runtime_error("Empty category file: " + category_file);
        auto header = split_ws(header_line);
        std::unordered_map<std::string, int> cols;
        for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;
        const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
        const int chr_col = find_col(cols, {"CHR", "CHROM"});
        const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
        if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) {
            throw std::runtime_error("category file needs RS/SNP/ID, or CHR+PS/POS columns");
        }
        std::vector<int> cat_cols;
        for (int i = 0; i < static_cast<int>(header.size()); ++i) {
            if (!is_marker_col(header[i])) cat_cols.push_back(i);
        }
        if (cat_cols.empty()) throw std::runtime_error("category file has no category columns");
        spec.ncat = std::max(spec.ncat, static_cast<int>(cat_cols.size()));
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_ws(line);
            if (f.size() < header.size()) continue;
            std::string rs = rs_col >= 0 ? f[rs_col] : (f[chr_col] + ":" + f[pos_col]);
            for (int c = 0; c < spec.ncat; ++c) {
                const double v = std::stod(f[cat_cols[c]]);
                if (std::abs(v - 1.0) < 1e-12) {
                    spec.cat[rs] = c;
                    break;
                }
            }
        }
    };
    if (!cfg.category_file.empty()) {
        read_one_cat(cfg.category_file);
    }
    if (!cfg.multi_category_file.empty()) {
        for (const auto& path : read_prefix_list(cfg.multi_category_file)) {
            read_one_cat(path);
        }
    }
    if (!cfg.category_weight_file.empty()) {
        std::ifstream in(cfg.category_weight_file);
        if (!in) throw std::runtime_error("Cannot open wcat file: " + cfg.category_weight_file);
        std::string header_line;
        if (!std::getline(in, header_line)) throw std::runtime_error("Empty wcat file: " + cfg.category_weight_file);
        auto header = split_ws(header_line);
        std::unordered_map<std::string, int> cols;
        for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;
        const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
        const int chr_col = find_col(cols, {"CHR", "CHROM"});
        const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
        if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) {
            throw std::runtime_error("wcat file needs RS/SNP/ID, or CHR+PS/POS columns");
        }
        std::vector<int> weight_cols;
        for (int i = 0; i < static_cast<int>(header.size()); ++i) {
            if (!is_marker_col(header[i])) weight_cols.push_back(i);
        }
        if (spec.ncat == 1 && cfg.category_file.empty() && cfg.multi_category_file.empty()) {
            spec.ncat = static_cast<int>(weight_cols.size());
        }
        if (static_cast<int>(weight_cols.size()) != spec.ncat) {
            throw std::runtime_error("wcat category column count does not match --cat");
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto f = split_ws(line);
            if (f.size() < header.size()) continue;
            std::string rs = rs_col >= 0 ? f[rs_col] : (f[chr_col] + ":" + f[pos_col]);
            std::vector<double> w(spec.ncat, 0.0);
            for (int c = 0; c < spec.ncat; ++c) w[c] = std::stod(f[weight_cols[c]]);
            spec.wcat[rs] = std::move(w);
        }
    }
    if (!cfg.snp_weight_summary_file.empty()) {
        std::ifstream in(cfg.snp_weight_summary_file);
        if (!in) throw std::runtime_error("Cannot open wsnp file: " + cfg.snp_weight_summary_file);
        std::string rs;
        double w = 0.0;
        while (in >> rs >> w) spec.wsnp[rs] = w;
    }
    return spec;
}

int category_for_snp(const CategorySpec& spec, const std::string& rs) {
    if (spec.cat.empty()) return 0;
    auto it = spec.cat.find(rs);
    if (it == spec.cat.end()) return -1;
    return it->second;
}

double weight_for_snp(const CategorySpec& spec, const std::string& rs, int cat, double fallback = 1.0) {
    double w = fallback;
    auto wit = spec.wcat.find(rs);
    if (wit != spec.wcat.end() && cat >= 0 && cat < static_cast<int>(wit->second.size())) {
        w *= wit->second[cat];
    }
    auto sit = spec.wsnp.find(rs);
    if (sit != spec.wsnp.end()) w *= sit->second;
    return w;
}

std::vector<BetaRecord> read_beta_records(const GwasConfig& cfg, const CategorySpec& spec) {
    std::ifstream in(cfg.beta_file);
    if (!in) throw std::runtime_error("Cannot open beta file: " + cfg.beta_file);
    std::string header_line;
    if (!std::getline(in, header_line)) throw std::runtime_error("Empty beta file: " + cfg.beta_file);
    auto header = split_ws(header_line);
    std::unordered_map<std::string, int> cols;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;
    const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
    const int chr_col = find_col(cols, {"CHR", "CHROM"});
    const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
    const int z_col = find_col(cols, {"Z", "ZSCORE", "Z_SCORE"});
    const int beta_col = find_col(cols, {"BETA", "B"});
    const int se_col = find_col(cols, {"SE", "SE_BETA", "SEBETA"});
    const int chisq_col = find_col(cols, {"CHISQ", "CHISQ_STAT", "CHISQUARE", "WALD", "LRT", "SCORE"});
    const int p_col = find_col(cols, {"P", "PVALUE", "P_VALUE", "P-VALUE"});
    const int n_col = find_col(cols, {"N", "N_TOTAL", "NTOTAL"});
    const int nmis_col = find_col(cols, {"N_MIS", "NMIS", "N_MISS"});
    const int nobs_col = find_col(cols, {"N_OBS", "NOBS"});
    const int ncase_col = find_col(cols, {"NCASE", "N_CASE"});
    const int ncontrol_col = find_col(cols, {"NCONTROL", "N_CONTROL"});
    if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) throw std::runtime_error("beta file needs SNP ID columns");
    if (z_col < 0 && (beta_col < 0 || se_col < 0) && chisq_col < 0 && p_col < 0) {
        throw std::runtime_error("beta file needs Z, or BETA+SE, or CHISQ, or P columns");
    }
    if (n_col < 0 && (nmis_col < 0 || nobs_col < 0) && (ncase_col < 0 || ncontrol_col < 0)) {
        throw std::runtime_error("beta file needs N, N_MIS+N_OBS, or NCASE+NCONTROL columns");
    }

    std::vector<BetaRecord> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < header.size()) continue;
        const std::string rs = rs_col >= 0 ? f[rs_col] : (f[chr_col] + ":" + f[pos_col]);
        const int cat = category_for_snp(spec, rs);
        if (cat < 0) continue;
        double z2 = 0.0;
        if (z_col >= 0) {
            double z = std::stod(f[z_col]);
            z2 = z * z;
        } else if (beta_col >= 0 && se_col >= 0) {
            double se = std::stod(f[se_col]);
            if (se == 0.0) continue;
            double z = std::stod(f[beta_col]) / se;
            z2 = z * z;
        } else if (chisq_col >= 0) {
            z2 = std::stod(f[chisq_col]);
        } else {
            z2 = chisq1_from_p(std::stod(f[p_col]));
        }
        int n = 0;
        if (n_col >= 0) n = std::stoi(f[n_col]);
        else if (nmis_col >= 0 && nobs_col >= 0) n = std::stoi(f[nmis_col]) + std::stoi(f[nobs_col]);
        else n = std::stoi(f[ncase_col]) + std::stoi(f[ncontrol_col]);
        const double w = weight_for_snp(spec, rs, cat, 1.0);
        if (n <= 0 || w == 0.0 || !std::isfinite(z2) || !std::isfinite(w)) continue;
        rows.push_back({rs, cat, n, w, z2});
    }
    if (rows.empty()) throw std::runtime_error("found no usable summary-stat rows");
    return rows;
}

SummaryQMulti compute_summary_q_multi(const GwasConfig& cfg, const CategorySpec& spec) {
    auto rows = read_beta_records(cfg, spec);
    SummaryQMulti out;
    out.q = Eigen::VectorXd::Zero(spec.ncat);
    out.vq = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
    out.s = Eigen::VectorXd::Zero(spec.ncat);
    Eigen::VectorXd raw_q = Eigen::VectorXd::Zero(spec.ncat);
    Eigen::VectorXi n_snps = Eigen::VectorXi::Zero(spec.ncat);
    out.total = static_cast<int>(rows.size());
    for (const auto& r : rows) {
        raw_q(r.cat) += (r.z2 - 1.0) * r.weight / static_cast<double>(r.n);
        out.s(r.cat) += r.weight;
        n_snps(r.cat) += 1;
        out.ni_total = std::max(out.ni_total, r.n);
        ++out.usable;
    }
    for (int c = 0; c < spec.ncat; ++c) {
        if (out.s(c) != 0.0) out.q(c) = raw_q(c) / out.s(c);
    }

    const int n_block = std::max(1, cfg.summary_block_count);
    for (int l = 0; l < spec.ncat; ++l) {
        const int per_block = static_cast<int>(std::floor(n_snps(l) / static_cast<double>(n_block)));
        if (per_block <= 0) continue;
        std::vector<Eigen::VectorXd> block_q(n_block, Eigen::VectorXd::Zero(spec.ncat));
        std::vector<Eigen::VectorXd> block_s(n_block, Eigen::VectorXd::Zero(spec.ncat));
        int b = 0, t = 0;
        for (const auto& r : rows) {
            block_q[b](r.cat) += (r.z2 - 1.0) * r.weight;
            block_s[b](r.cat) += r.weight;
            if (r.cat == l && b < n_block - 1) {
                if (t < per_block - 1) ++t;
                else { ++b; t = 0; }
            }
        }
        std::vector<Eigen::VectorXd> jack(n_block, Eigen::VectorXd::Zero(spec.ncat));
        Eigen::VectorXd mean = Eigen::VectorXd::Zero(spec.ncat);
        Eigen::VectorXi count = Eigen::VectorXi::Zero(spec.ncat);
        for (int k = 0; k < n_block; ++k) {
            for (int i = 0; i < spec.ncat; ++i) {
                if (block_s[k](i) != 0.0 && out.s(i) != block_s[k](i)) {
                    jack[k](i) = (raw_q(i) - block_q[k](i)) / (out.s(i) - block_s[k](i));
                    mean(i) += jack[k](i);
                    count(i) += 1;
                }
            }
        }
        for (int i = 0; i < spec.ncat; ++i) if (count(i) > 0) mean(i) /= static_cast<double>(count(i));
        for (int k = 0; k < n_block; ++k) {
            for (int i = 0; i < spec.ncat; ++i) if (jack[k](i) != 0.0) jack[k](i) -= mean(i);
        }
        for (int i = 0; i < spec.ncat; ++i) {
            double d = 0.0;
            int n = 0;
            for (int k = 0; k < n_block; ++k) {
                if (jack[k](l) != 0.0 && jack[k](i) != 0.0) {
                    d += jack[k](l) * jack[k](i);
                    ++n;
                }
            }
            if (n != 0) d = d / static_cast<double>(n) * static_cast<double>(n - 1);
            out.vq(i, l) += d;
            out.vq(l, i) = out.vq(i, l);
        }
    }
    for (int i = 0; i < spec.ncat; ++i) {
        for (int j = i + 1; j < spec.ncat; ++j) {
            out.vq(i, j) *= 0.5;
            out.vq(j, i) = out.vq(i, j);
        }
    }
    return out;
}

void write_vector_file(const std::string& path, const Eigen::VectorXd& v) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write vector file: " + path);
    out << std::setprecision(10);
    for (int i = 0; i < v.size(); ++i) out << v(i) << "\n";
}

void write_matrix_file(const std::string& path, const Eigen::MatrixXd& m) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write matrix file: " + path);
    out << std::setprecision(10);
    for (int i = 0; i < m.rows(); ++i) {
        for (int j = 0; j < m.cols(); ++j) {
            if (j) out << '\t';
            out << m(i, j);
        }
        out << "\n";
    }
}

CorSummary read_cor_summary(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open cor file: " + path);
    std::string header_line;
    if (!std::getline(in, header_line)) throw std::runtime_error("Empty cor file: " + path);
    auto header = split_ws(header_line);
    std::unordered_map<std::string, int> cols;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;

    const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
    const int chr_col = find_col(cols, {"CHR", "CHROM"});
    const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
    const int cm_col = find_col(cols, {"CM"});
    const int n_col = find_col(cols, {"N", "N_TOTAL", "NTOTAL"});
    const int nmis_col = find_col(cols, {"N_MIS", "NMIS", "N_MISS"});
    const int nobs_col = find_col(cols, {"N_OBS", "NOBS"});
    const int af_col = find_col(cols, {"AF", "MAF", "F", "ALLELE_FREQ", "ALLELE_FREQUENCY"});
    const int var_col = find_col(cols, {"VAR"});
    const int cor_col = find_col(cols, {"COR", "R"});
    if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) {
        throw std::runtime_error("cor file needs RS, or CHR+PS/POS columns");
    }
    if (n_col < 0 && (nmis_col < 0 || nobs_col < 0)) {
        throw std::runtime_error("cor file needs N, or N_MIS+N_OBS columns");
    }
    if (var_col < 0 && af_col < 0) {
        throw std::runtime_error("cor file needs VAR or AF column");
    }
    if (cor_col < 0) throw std::runtime_error("cor file needs COR/R column");
    if (cor_col != static_cast<int>(header.size()) - 1) {
        throw std::runtime_error("cor file COR/R column must be the last fixed column");
    }

    CorSummary out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < header.size() - 1) continue;
        CorRecord rec;
        if (rs_col >= 0) rec.rs = f[rs_col];
        else rec.rs = f[chr_col] + ":" + f[pos_col];
        if (out.index.count(rec.rs)) {
            throw std::runtime_error("duplicate SNP in cor file: " + rec.rs);
        }
        if (n_col >= 0) rec.n = std::stoi(f[n_col]);
        else rec.n = std::stoi(f[nmis_col]) + std::stoi(f[nobs_col]);
        if (pos_col >= 0) rec.bp = std::stod(f[pos_col]);
        if (cm_col >= 0) rec.cm = std::stod(f[cm_col]);
        if (var_col >= 0) {
            rec.var = std::stod(f[var_col]);
        } else {
            const double af = std::stod(f[af_col]);
            rec.var = 2.0 * af * (1.0 - af);
        }
        for (size_t k = static_cast<size_t>(cor_col); k < f.size(); ++k) {
            const std::string val = upper(f[k]);
            if (val == "NA" || val == "NAN" || val == "-NAN") {
                rec.cor.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                rec.cor.push_back(std::stod(f[k]));
            }
        }
        if (rec.n <= 0 || rec.var <= 0.0 || !std::isfinite(rec.var)) continue;
        out.index[rec.rs] = out.records.size();
        out.records.push_back(std::move(rec));
    }
    if (out.records.empty()) throw std::runtime_error("cor file has no usable SNP rows");
    return out;
}

SummaryQ compute_summary_q_with_cor(const GwasConfig& cfg, CorSummary& cor) {
    std::ifstream in(cfg.beta_file);
    if (!in) throw std::runtime_error("Cannot open beta file: " + cfg.beta_file);
    std::string header_line;
    if (!std::getline(in, header_line)) throw std::runtime_error("Empty beta file: " + cfg.beta_file);
    auto header = split_ws(header_line);
    std::unordered_map<std::string, int> cols;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;

    const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
    const int chr_col = find_col(cols, {"CHR", "CHROM"});
    const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
    const int z_col = find_col(cols, {"Z", "ZSCORE", "Z_SCORE"});
    const int beta_col = find_col(cols, {"BETA", "B"});
    const int se_col = find_col(cols, {"SE", "SE_BETA", "SEBETA"});
    const int chisq_col = find_col(cols, {"CHISQ", "CHISQ_STAT", "CHISQUARE", "WALD", "LRT", "SCORE"});
    const int p_col = find_col(cols, {"P", "PVALUE", "P_VALUE", "P-VALUE"});
    const int n_col = find_col(cols, {"N", "N_TOTAL", "NTOTAL"});
    const int nmis_col = find_col(cols, {"N_MIS", "NMIS", "N_MISS"});
    const int nobs_col = find_col(cols, {"N_OBS", "NOBS"});
    const int af_col = find_col(cols, {"AF", "MAF", "F", "ALLELE_FREQ", "ALLELE_FREQUENCY"});
    const int var_col = find_col(cols, {"VAR"});
    if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) {
        throw std::runtime_error("beta file needs RS/SNP/ID, or CHR+PS/POS columns for --beta + --cor");
    }
    if (z_col < 0 && (beta_col < 0 || se_col < 0) && chisq_col < 0 && p_col < 0) {
        throw std::runtime_error("beta file needs Z, or BETA+SE, or CHISQ, or P columns");
    }
    if (n_col < 0 && (nmis_col < 0 || nobs_col < 0)) {
        throw std::runtime_error("beta file needs N, or N_MIS+N_OBS columns");
    }

    SummaryQ res;
    double qvar_diag = 0.0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < header.size()) continue;
        ++res.total;

        std::string rs = rs_col >= 0 ? f[rs_col] : (f[chr_col] + ":" + f[pos_col]);
        auto it = cor.index.find(rs);
        if (it == cor.index.end()) continue;
        CorRecord& rec = cor.records[it->second];
        if (rec.in_beta) throw std::runtime_error("duplicate SNP in beta file: " + rs);

        double z2 = 0.0;
        if (z_col >= 0) {
            double z = std::stod(f[z_col]);
            z2 = z * z;
        } else if (beta_col >= 0 && se_col >= 0) {
            double se = std::stod(f[se_col]);
            if (se == 0.0) continue;
            double z = std::stod(f[beta_col]) / se;
            z2 = z * z;
        } else if (chisq_col >= 0) {
            z2 = std::stod(f[chisq_col]);
        } else if (p_col >= 0) {
            z2 = chisq1_from_p(std::stod(f[p_col]));
        }
        int n = n_col >= 0 ? std::stoi(f[n_col]) : (std::stoi(f[nmis_col]) + std::stoi(f[nobs_col]));
        if (n <= 0 || !std::isfinite(z2)) continue;

        double var_x = rec.var;
        if (var_col >= 0) var_x = std::stod(f[var_col]);
        else if (af_col >= 0) {
            const double af = std::stod(f[af_col]);
            var_x = 2.0 * af * (1.0 - af);
        }
        if (var_x <= 0.0 || !std::isfinite(var_x)) continue;

        rec.in_beta = true;
        rec.beta_n = n;
        rec.var = var_x;
        res.q += (z2 - 1.0) * var_x / static_cast<double>(n);
        res.snp_weight_sum += var_x;
        qvar_diag += var_x * var_x / (static_cast<double>(n) * static_cast<double>(n));
        res.ni_total = std::max(res.ni_total, n);
        ++res.usable;
    }
    if (res.usable == 0) throw std::runtime_error("found no beta/cor overlapping SNPs");
    res.vq = 2.0 * qvar_diag;
    return res;
}

struct SummaryS {
    double s = 0.0;
    double svar = 0.0;
    double qvar_ld = 0.0;
    int ns_pair = 0;
};

struct DirectSummaryMulti {
    Eigen::VectorXd q;
    Eigen::MatrixXd vq;
    Eigen::VectorXd s;
    Eigen::MatrixXd S;
    Eigen::MatrixXd Svar;
    int ni_total = 0;
    int overlap = 0;
    int ns_pair = 0;
};

struct PrefixSummaryData {
    Eigen::MatrixXd Vq;
    Eigen::VectorXd q;
    Eigen::VectorXd study_s;
    Eigen::MatrixXd S;
    Eigen::MatrixXd Svar;
    Eigen::VectorXd ref_s;
    double n_study = 0.0;
    double n_ref = 0.0;
};

SummaryS compute_s_from_cor(const CorSummary& cor) {
    SummaryS out;
    int ns_test = 0;
    for (size_t i = 0; i < cor.records.size(); ++i) {
        const auto& r1 = cor.records[i];
        if (!r1.in_beta) continue;
        ++ns_test;
        const double var1 = r1.var;
        const int n1 = r1.n;
        const double diag = var1 * var1;
        out.s += (1.0 - 1.0 / static_cast<double>(n1)) * diag;
        out.svar += diag * diag / (static_cast<double>(n1) * static_cast<double>(n1));

        for (size_t k = 0; k < r1.cor.size(); ++k) {
            const size_t j = i + k + 1;
            if (j >= cor.records.size()) break;
            const auto& r2 = cor.records[j];
            if (!r2.in_beta) continue;
            const double corr = r1.cor[k];
            if (!std::isfinite(corr)) continue;
            const double var2 = r2.var;
            const int n2 = r2.n;
            const int n12 = std::min(n1, n2);
            const double d1 = corr * corr - 1.0 / static_cast<double>(n12);
            const double d2 = var1 * var2;
            const double d3 = corr * corr / (static_cast<double>(r1.beta_n) * static_cast<double>(r2.beta_n));
            out.qvar_ld += 2.0 * d3 * d2;
            out.s += 2.0 * d1 * d2;
            out.svar += 2.0 * d2 * d2 / (static_cast<double>(n12) * static_cast<double>(n12));
            ++out.ns_pair;
        }
    }
    if (ns_test == 0) throw std::runtime_error("cor file has no SNPs overlapping beta file");
    if (out.ns_pair > 0) {
        out.svar = 2.0 * out.svar * static_cast<double>(ns_test) * static_cast<double>(ns_test) /
                   (2.0 * static_cast<double>(out.ns_pair));
    } else {
        out.svar *= 2.0;
    }
    if (std::abs(out.s) < 1e-20) throw std::runtime_error("beta/cor summary S is numerically zero");
    return out;
}

DirectSummaryMulti compute_direct_multi(const GwasConfig& cfg, const CategorySpec& spec) {
    CorSummary cor = read_cor_summary(cfg.cor_file);
    DirectSummaryMulti out;
    out.q = Eigen::VectorXd::Zero(spec.ncat);
    out.vq = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
    out.s = Eigen::VectorXd::Zero(spec.ncat);
    out.S = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
    out.Svar = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);

    std::ifstream in(cfg.beta_file);
    if (!in) throw std::runtime_error("Cannot open beta file: " + cfg.beta_file);
    std::string header_line;
    if (!std::getline(in, header_line)) throw std::runtime_error("Empty beta file: " + cfg.beta_file);
    auto header = split_ws(header_line);
    std::unordered_map<std::string, int> cols;
    for (int i = 0; i < static_cast<int>(header.size()); ++i) cols[upper(header[i])] = i;
    const int rs_col = find_col(cols, {"RS", "SNP", "ID"});
    const int chr_col = find_col(cols, {"CHR", "CHROM"});
    const int pos_col = find_col(cols, {"PS", "POS", "BP", "POSITION"});
    const int z_col = find_col(cols, {"Z", "ZSCORE", "Z_SCORE"});
    const int beta_col = find_col(cols, {"BETA", "B"});
    const int se_col = find_col(cols, {"SE", "SE_BETA", "SEBETA"});
    const int chisq_col = find_col(cols, {"CHISQ", "CHISQ_STAT", "CHISQUARE", "WALD", "LRT", "SCORE"});
    const int p_col = find_col(cols, {"P", "PVALUE", "P_VALUE", "P-VALUE"});
    const int n_col = find_col(cols, {"N", "N_TOTAL", "NTOTAL"});
    const int nmis_col = find_col(cols, {"N_MIS", "NMIS", "N_MISS"});
    const int nobs_col = find_col(cols, {"N_OBS", "NOBS"});
    const int af_col = find_col(cols, {"AF", "MAF", "F", "ALLELE_FREQ", "ALLELE_FREQUENCY"});
    const int var_col = find_col(cols, {"VAR"});
    if (rs_col < 0 && (chr_col < 0 || pos_col < 0)) throw std::runtime_error("beta file needs SNP ID columns");
    if (z_col < 0 && (beta_col < 0 || se_col < 0) && chisq_col < 0 && p_col < 0) {
        throw std::runtime_error("beta file needs Z, or BETA+SE, or CHISQ, or P columns");
    }
    if (n_col < 0 && (nmis_col < 0 || nobs_col < 0)) throw std::runtime_error("beta file needs N or N_MIS+N_OBS columns");

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < header.size()) continue;
        const std::string rs = rs_col >= 0 ? f[rs_col] : (f[chr_col] + ":" + f[pos_col]);
        auto it = cor.index.find(rs);
        if (it == cor.index.end()) continue;
        const int cat = category_for_snp(spec, rs);
        if (cat < 0) continue;
        CorRecord& rec = cor.records[it->second];
        if (rec.in_beta) throw std::runtime_error("duplicate SNP in beta file: " + rs);

        double z2 = 0.0;
        if (z_col >= 0) {
            double z = std::stod(f[z_col]);
            z2 = z * z;
        } else if (beta_col >= 0 && se_col >= 0) {
            double se = std::stod(f[se_col]);
            if (se == 0.0) continue;
            double z = std::stod(f[beta_col]) / se;
            z2 = z * z;
        } else if (chisq_col >= 0) {
            z2 = std::stod(f[chisq_col]);
        } else {
            z2 = chisq1_from_p(std::stod(f[p_col]));
        }
        int n = n_col >= 0 ? std::stoi(f[n_col]) : (std::stoi(f[nmis_col]) + std::stoi(f[nobs_col]));
        if (n <= 0 || !std::isfinite(z2)) continue;
        double var_x = rec.var;
        if (var_col >= 0) var_x = std::stod(f[var_col]);
        else if (af_col >= 0) {
            const double af = std::stod(f[af_col]);
            var_x = 2.0 * af * (1.0 - af);
        }
        const double weight = weight_for_snp(spec, rs, cat, var_x);
        if (weight == 0.0 || !std::isfinite(weight)) continue;
        rec.in_beta = true;
        rec.beta_n = n;
        rec.cat = cat;
        rec.weight = weight;
        out.q(cat) += (z2 - 1.0) * weight / static_cast<double>(n);
        out.s(cat) += weight;
        out.vq(cat, cat) += 2.0 * weight * weight / (static_cast<double>(n) * static_cast<double>(n));
        out.ni_total = std::max(out.ni_total, n);
        ++out.overlap;
    }
    if (out.overlap == 0) throw std::runtime_error("found no beta/cor/category overlapping SNPs");

    int ns_test = 0;
    for (size_t i = 0; i < cor.records.size(); ++i) {
        const auto& r1 = cor.records[i];
        if (!r1.in_beta) continue;
        ++ns_test;
        const int c1 = r1.cat;
        const double w1 = r1.weight;
        const int n1 = r1.n;
        out.S(c1, c1) += (1.0 - 1.0 / static_cast<double>(n1)) * w1 * w1;
        out.Svar(c1, c1) += w1 * w1 * w1 * w1 / (static_cast<double>(n1) * static_cast<double>(n1));
        for (size_t k = 0; k < r1.cor.size(); ++k) {
            const size_t j = i + k + 1;
            if (j >= cor.records.size()) break;
            const auto& r2 = cor.records[j];
            if (!r2.in_beta) continue;
            const double corr = r1.cor[k];
            if (!std::isfinite(corr)) continue;
            const int c2 = r2.cat;
            const int n12 = std::min(r1.n, r2.n);
            const double d1 = corr * corr - 1.0 / static_cast<double>(n12);
            const double d2 = w1 * r2.weight;
            const double d3 = corr * corr / (static_cast<double>(r1.beta_n) * static_cast<double>(r2.beta_n));
            if (c1 == c2) {
                out.vq(c1, c1) += 2.0 * d3 * d2;
                out.S(c1, c1) += 2.0 * d1 * d2;
                out.Svar(c1, c1) += 2.0 * d2 * d2 / (static_cast<double>(n12) * static_cast<double>(n12));
            } else {
                out.S(c1, c2) += d1 * d2;
                out.S(c2, c1) += d1 * d2;
                out.Svar(c1, c2) += d2 * d2 / (static_cast<double>(n12) * static_cast<double>(n12));
                out.Svar(c2, c1) += d2 * d2 / (static_cast<double>(n12) * static_cast<double>(n12));
            }
            ++out.ns_pair;
        }
    }
    if (out.ns_pair > 0) {
        out.Svar *= 2.0 * static_cast<double>(ns_test) * static_cast<double>(ns_test) /
                    (2.0 * static_cast<double>(out.ns_pair));
    } else {
        out.Svar *= 2.0;
    }
    return out;
}

void write_summary_vc_outputs(const GwasConfig& cfg,
                              const Eigen::MatrixXd& Vq,
                              const Eigen::MatrixXd& S,
                              const Eigen::MatrixXd& Svar,
                              const Eigen::VectorXd& q,
                              const Eigen::VectorXd& s,
                              double n_study,
                              const std::string& method,
                              const std::string& meta_extra) {
    if (S.rows() != S.cols() || S.rows() != q.size() || s.size() != q.size()) {
        throw std::runtime_error("summary VC matrix dimensions do not agree");
    }
    std::string solver = "full-piv-lu";
    Eigen::FullPivLU<Eigen::MatrixXd> lu(S);
    Eigen::MatrixXd Si;
    if (lu.isInvertible()) {
        Si = lu.inverse();
    } else {
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::VectorXd singular = svd.singularValues();
        const double tol = std::max(S.rows(), S.cols()) * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, singular.size() > 0 ? singular.maxCoeff() : 0.0);
        Eigen::VectorXd inv = singular;
        for (int i = 0; i < inv.size(); ++i) inv(i) = inv(i) > tol ? 1.0 / inv(i) : 0.0;
        Si = svd.matrixV() * inv.asDiagonal() * svd.matrixU().transpose();
        solver = "svd-pseudoinverse";
    }
    Eigen::VectorXd pve = Si * q;
    Eigen::VectorXd sigma = Eigen::VectorXd::Zero(pve.size());
    for (int i = 0; i < pve.size(); ++i) {
        sigma(i) = std::abs(s(i)) > 0.0 ? pve(i) / s(i) : 0.0;
    }
    Eigen::MatrixXd var_mat = Eigen::MatrixXd::Zero(S.rows(), S.cols());
    for (int i = 0; i < S.rows(); ++i) {
        for (int j = 0; j < S.cols(); ++j) {
            var_mat(i, j) = Svar(i, j) * pve(i) * pve(j) + Vq(i, j) / (n_study * n_study);
        }
    }
    Eigen::MatrixXd cov_pve = Si * var_mat * Si.transpose();
    const double total_pve = pve.sum();
    const double total_var = cov_pve.sum();
    const double total_se = std::sqrt(std::max(0.0, total_var));
    const double total_s = s.sum();

    std::ofstream out(cfg.out_prefix + ".summary.vc.txt");
    if (!out) throw std::runtime_error("Cannot write summary-stat --vc output");
    out << "component\tpve\tse_pve\tsigma2_per_snp\tse_sigma2_per_snp\tenrichment\tse_enrichment\n";
    out << std::scientific << std::setprecision(10);
    for (int i = 0; i < pve.size(); ++i) {
        const double se = std::sqrt(std::max(0.0, cov_pve(i, i)));
        const double se_sigma = s(i) != 0.0 ? se / s(i) : 0.0;
        const double expected = (total_s != 0.0 && total_pve != 0.0) ? (s(i) / total_s * total_pve) : 0.0;
        const double enrich = expected != 0.0 ? pve(i) / expected : 1.0;
        out << "Vg" << (i + 1) << '\t' << pve(i) << '\t' << se << '\t'
            << sigma(i) << '\t' << se_sigma << '\t' << enrich << "\tNA\n";
    }
    out << "total_pve\t" << total_pve << '\t' << total_se << "\tNA\tNA\tNA\tNA\n";

    std::ofstream meta(cfg.out_prefix + ".summary.vc.meta.txt");
    if (!meta) throw std::runtime_error("Cannot write summary-stat --vc metadata");
    meta << "method\t" << method << "\n";
    meta << "S_solver\t" << solver << "\n";
    meta << "n_components\t" << q.size() << "\n";
    meta << "n_study\t" << std::setprecision(10) << n_study << "\n";
    meta << meta_extra;
    for (int i = 0; i < q.size(); ++i) meta << "q" << (i + 1) << '\t' << q(i) << "\n";
}

PrefixSummaryData read_study_prefix(const std::string& prefix) {
    auto vq = read_numbers(prefix + ".Vq.txt");
    auto qv = read_numbers(prefix + ".q.txt");
    auto size = read_numbers(prefix + ".size.txt");
    if (size.size() < 2 || qv.empty()) throw std::runtime_error("invalid study prefix: " + prefix);
    const int nvc = static_cast<int>(size.size()) - 1;
    if (static_cast<int>(qv.size()) != nvc || static_cast<int>(vq.size()) != nvc * nvc) {
        throw std::runtime_error("study prefix matrix dimensions do not match: " + prefix);
    }
    PrefixSummaryData d;
    d.Vq = Eigen::MatrixXd::Zero(nvc, nvc);
    d.q = Eigen::VectorXd::Zero(nvc);
    d.study_s = Eigen::VectorXd::Zero(nvc);
    for (int i = 0; i < nvc; ++i) {
        d.q(i) = qv[i];
        d.study_s(i) = size[i];
        for (int j = 0; j < nvc; ++j) d.Vq(i, j) = vq[i * nvc + j];
    }
    d.n_study = size.back();
    return d;
}

PrefixSummaryData read_ref_prefix(const std::string& prefix) {
    auto s = read_numbers(prefix + ".S.txt");
    auto size = read_numbers(prefix + ".size.txt");
    if (size.size() < 2) throw std::runtime_error("invalid ref prefix: " + prefix);
    const int nvc = static_cast<int>(size.size()) - 1;
    if (static_cast<int>(s.size()) != 2 * nvc * nvc) {
        throw std::runtime_error("ref prefix matrix dimensions do not match: " + prefix);
    }
    PrefixSummaryData d;
    d.S = Eigen::MatrixXd::Zero(nvc, nvc);
    d.Svar = Eigen::MatrixXd::Zero(nvc, nvc);
    d.ref_s = Eigen::VectorXd::Zero(nvc);
    for (int i = 0; i < nvc; ++i) {
        d.ref_s(i) = size[i];
        for (int j = 0; j < nvc; ++j) {
            d.S(i, j) = s[i * nvc + j];
            d.Svar(i, j) = s[nvc * nvc + i * nvc + j];
        }
    }
    d.n_ref = size.back();
    return d;
}

PrefixSummaryData read_or_merge_studies(const GwasConfig& cfg) {
    if (!cfg.multi_study_file.empty()) {
        auto prefixes = read_prefix_list(cfg.multi_study_file);
        PrefixSummaryData acc;
        bool first = true;
        for (const auto& p : prefixes) {
            PrefixSummaryData cur = read_study_prefix(p);
            if (first) {
                acc.Vq = Eigen::MatrixXd::Zero(cur.Vq.rows(), cur.Vq.cols());
                acc.q = Eigen::VectorXd::Zero(cur.q.size());
                acc.study_s = Eigen::VectorXd::Zero(cur.study_s.size());
                first = false;
            }
            if (cur.q.size() != acc.q.size()) throw std::runtime_error("mstudy prefixes have inconsistent dimensions");
            for (int i = 0; i < cur.q.size(); ++i) {
                acc.q(i) += cur.q(i) * cur.study_s(i);
                acc.study_s(i) += cur.study_s(i);
                for (int j = 0; j < cur.q.size(); ++j) {
                    acc.Vq(i, j) += cur.Vq(i, j) * cur.study_s(i) * cur.study_s(j);
                }
            }
            acc.n_study = std::max(acc.n_study, cur.n_study);
        }
        for (int i = 0; i < acc.q.size(); ++i) {
            if (acc.study_s(i) == 0.0) continue;
            acc.q(i) /= acc.study_s(i);
            for (int j = 0; j < acc.q.size(); ++j) {
                if (acc.study_s(j) != 0.0) acc.Vq(i, j) /= (acc.study_s(i) * acc.study_s(j));
            }
        }
        return acc;
    }
    if (cfg.study_prefix.empty()) throw std::runtime_error("summary-stat route requires --study/-study or --mstudy/-mstudy");
    return read_study_prefix(cfg.study_prefix);
}

PrefixSummaryData read_or_merge_refs(const GwasConfig& cfg) {
    if (!cfg.multi_ref_file.empty()) {
        auto prefixes = read_prefix_list(cfg.multi_ref_file);
        PrefixSummaryData acc;
        bool first = true;
        for (const auto& p : prefixes) {
            PrefixSummaryData cur = read_ref_prefix(p);
            if (first) {
                acc.S = Eigen::MatrixXd::Zero(cur.S.rows(), cur.S.cols());
                acc.Svar = Eigen::MatrixXd::Zero(cur.Svar.rows(), cur.Svar.cols());
                acc.ref_s = Eigen::VectorXd::Zero(cur.ref_s.size());
                first = false;
            }
            if (cur.ref_s.size() != acc.ref_s.size()) throw std::runtime_error("mref prefixes have inconsistent dimensions");
            acc.ref_s += cur.ref_s;
            acc.n_ref = std::max(acc.n_ref, cur.n_ref);
            for (int i = 0; i < cur.ref_s.size(); ++i) {
                for (int j = 0; j < cur.ref_s.size(); ++j) {
                    acc.S(i, j) += cur.S(i, j) * cur.ref_s(i) * cur.ref_s(j);
                    acc.Svar(i, j) += cur.Svar(i, j) * cur.ref_s(i) * cur.ref_s(j) * cur.ref_s(i) * cur.ref_s(j);
                }
            }
        }
        for (int i = 0; i < acc.ref_s.size(); ++i) {
            if (acc.ref_s(i) == 0.0) continue;
            for (int j = 0; j < acc.ref_s.size(); ++j) {
                if (acc.ref_s(j) == 0.0) continue;
                acc.S(i, j) /= (acc.ref_s(i) * acc.ref_s(j));
                acc.Svar(i, j) /= (acc.ref_s(i) * acc.ref_s(j) * acc.ref_s(i) * acc.ref_s(j));
            }
        }
        return acc;
    }
    if (cfg.ref_prefix.empty()) throw std::runtime_error("summary-stat route requires --ref/-ref or --mref/-mref");
    return read_ref_prefix(cfg.ref_prefix);
}

void center_matrix(Eigen::MatrixXd& G) {
    const int n = G.rows();
    Eigen::VectorXd row_mean = G.rowwise().mean();
    double grand = row_mean.mean();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            G(i, j) = G(i, j) - row_mean(i) - row_mean(j) + grand;
        }
    }
}

double scale_mean_diag(Eigen::MatrixXd& G) {
    const double mean_diag = G.diagonal().mean();
    if (std::abs(mean_diag) > 1e-20) G.array() /= mean_diag;
    return mean_diag;
}

double comp_s_one(const Eigen::MatrixXd& A, const Eigen::MatrixXd& K, int n_cvt) {
    const int n = A.rows();
    double tr_AK = 0.0, sum_A = 0.0, sum_K = 0.0, sum_AK = 0.0, tr_A = 0.0, tr_K = 0.0;
    for (int l = 0; l < n; ++l) {
        double s_A = 0.0, s_K = 0.0;
        for (int k = 0; k < n; ++k) {
            const double a = A(l, k);
            const double b = K(l, k);
            s_A += a;
            s_K += b;
            tr_AK += a * b;
            sum_A += a;
            sum_K += b;
            if (l == k) {
                tr_A += a;
                tr_K += b;
            }
        }
        sum_AK += s_A * s_K;
    }
    sum_A /= static_cast<double>(n);
    sum_K /= static_cast<double>(n);
    sum_AK /= static_cast<double>(n);
    tr_A -= sum_A;
    tr_K -= sum_K;
    if (tr_A == 0.0 || tr_K == 0.0 || n == n_cvt) return 0.0;
    const double d = tr_AK - 2.0 * sum_AK + sum_A * sum_K;
    return d / (tr_A * tr_K) - 1.0 / static_cast<double>(n - n_cvt);
}

} // namespace

GemmaSummaryTools::GemmaSummaryTools(const GwasConfig& cfg) : cfg_(cfg) {}

void GemmaSummaryTools::runGs(PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              bool is_bgen,
                              bool is_pgen) {
    if (cfg_.gemma_s_mode != 1) {
        throw std::runtime_error("--gs currently supports GEMMA mode 1");
    }

    std::cout << "Calculating GEMMA-style S matrix...\n";
    if (!cfg_.category_file.empty() || !cfg_.multi_category_file.empty() ||
        !cfg_.category_weight_file.empty() || !cfg_.snp_weight_summary_file.empty()) {
        CategorySpec spec = read_category_spec(cfg_);
        std::vector<std::unordered_set<std::string>> cat_sets(spec.ncat);
        for (const auto& kv : spec.cat) {
            if (kv.second >= 0 && kv.second < spec.ncat) cat_sets[kv.second].insert(kv.first);
        }
        if (spec.cat.empty() && !spec.wcat.empty()) {
            for (const auto& kv : spec.wcat) {
                for (int c = 0; c < spec.ncat && c < static_cast<int>(kv.second.size()); ++c) {
                    if (kv.second[c] != 0.0) cat_sets[c].insert(kv.first);
                }
            }
        }
        if (spec.cat.empty() && spec.wcat.empty() && !spec.wsnp.empty()) {
            for (const auto& kv : spec.wsnp) if (kv.second != 0.0) cat_sets[0].insert(kv.first);
        }
        for (int c = 0; c < spec.ncat; ++c) {
            if (cat_sets[c].empty()) {
                throw std::runtime_error("--gs category " + std::to_string(c + 1) + " has no SNPs");
            }
        }

        const auto user_extract = load_ids(cfg_.extract_file);
        const auto user_exclude = load_ids(cfg_.exclude_file);
        std::vector<Eigen::MatrixXd> Kc(spec.ncat);
        Eigen::VectorXd size = Eigen::VectorXd::Zero(spec.ncat + 1);
        int n_ref = 0;
        for (int c = 0; c < spec.ncat; ++c) {
            GrmBuilder::Options opts;
            opts.min_maf = static_cast<float>(cfg_.min_maf);
            opts.min_hwe = static_cast<float>(cfg_.min_hwe);
            opts.algorithm = cfg_.grm_alg;
            if (user_extract.empty()) {
                opts.extract_set = cat_sets[c];
            } else {
                for (const auto& id : cat_sets[c]) if (user_extract.count(id)) opts.extract_set.insert(id);
            }
            opts.exclude_set = user_exclude;
            #ifdef _OPENMP
            opts.thread_num = cfg_.num_threads > 0 ? cfg_.num_threads : omp_get_max_threads();
            #endif

            GrmBuilder builder(opts);
            if (is_bgen) {
                bgen_reader.reset();
                builder.compute(bgen_reader);
            } else if (is_pgen) {
                pgen_reader.reset();
                builder.compute(pgen_reader);
            } else {
                builder.compute(plink_reader);
            }
            Kc[c] = builder.getGrm();
            center_matrix(Kc[c]);
            scale_mean_diag(Kc[c]);
            size(c) = builder.getTotalSnps();
            n_ref = Kc[c].rows();
        }
        size(spec.ncat) = n_ref;

        Eigen::MatrixXd S = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
        Eigen::MatrixXd Svar = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
        for (int i = 0; i < spec.ncat; ++i) {
            for (int j = i; j < spec.ncat; ++j) {
                S(i, j) = comp_s_one(Kc[i], Kc[j], 1);
                if (!std::isfinite(S(i, j))) S(i, j) = 0.0;
                S(j, i) = S(i, j);
            }
        }

        // Leave-one-sample jackknife approximation of GEMMA Svar for each component pair.
        std::vector<Eigen::MatrixXd> jack;
        jack.reserve(std::max(0, n_ref));
        for (int drop = 0; drop < n_ref; ++drop) {
            std::vector<Eigen::MatrixXd> subc(spec.ncat);
            for (int c = 0; c < spec.ncat; ++c) {
                subc[c].resize(n_ref - 1, n_ref - 1);
                int rr = 0;
                for (int r = 0; r < n_ref; ++r) {
                    if (r == drop) continue;
                    int cc = 0;
                    for (int col = 0; col < n_ref; ++col) {
                        if (col == drop) continue;
                        subc[c](rr, cc++) = Kc[c](r, col);
                    }
                    ++rr;
                }
                center_matrix(subc[c]);
                scale_mean_diag(subc[c]);
            }
            Eigen::MatrixXd s_drop = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
            for (int i = 0; i < spec.ncat; ++i) {
                for (int j = i; j < spec.ncat; ++j) {
                    s_drop(i, j) = comp_s_one(subc[i], subc[j], 1);
                    if (!std::isfinite(s_drop(i, j))) s_drop(i, j) = 0.0;
                    s_drop(j, i) = s_drop(i, j);
                }
            }
            jack.push_back(std::move(s_drop));
        }
        if (jack.size() > 1) {
            Eigen::MatrixXd mean = Eigen::MatrixXd::Zero(spec.ncat, spec.ncat);
            for (const auto& m : jack) mean += m;
            mean /= static_cast<double>(jack.size());
            for (const auto& m : jack) {
                Eigen::MatrixXd d = m - mean;
                Svar.array() += d.array().square();
            }
            Svar *= (static_cast<double>(jack.size() - 1) / static_cast<double>(jack.size()));
        }
        Eigen::MatrixXd stacked(2 * spec.ncat, spec.ncat);
        stacked.block(0, 0, spec.ncat, spec.ncat) = S;
        stacked.block(spec.ncat, 0, spec.ncat, spec.ncat) = Svar;
        write_matrix_file(cfg_.out_prefix + ".S.txt", stacked);
        write_vector_file(cfg_.out_prefix + ".size.txt", size);
        std::ofstream snpf(cfg_.out_prefix + ".snps.txt");
        if (!snpf) throw std::runtime_error("Cannot write GEMMA-style snps output");
        for (int c = 0; c < spec.ncat; ++c) {
            for (const auto& id : cat_sets[c]) snpf << id << "\tVg" << (c + 1) << "\n";
        }
        std::cout << "GEMMA-style multi-category S outputs written to " << cfg_.out_prefix
                  << ".S.txt, " << cfg_.out_prefix << ".size.txt, and "
                  << cfg_.out_prefix << ".snps.txt (" << spec.ncat << " categories)\n";
        return;
    }

    GrmBuilder::Options opts;
    opts.min_maf = static_cast<float>(cfg_.min_maf);
    opts.min_hwe = static_cast<float>(cfg_.min_hwe);
    opts.algorithm = cfg_.grm_alg;
    opts.extract_set = load_ids(cfg_.extract_file);
    opts.exclude_set = load_ids(cfg_.exclude_file);
    #ifdef _OPENMP
    opts.thread_num = cfg_.num_threads > 0 ? cfg_.num_threads : omp_get_max_threads();
    #endif

    GrmBuilder builder(opts);
    if (is_bgen) {
        bgen_reader.reset();
        builder.compute(bgen_reader);
    } else if (is_pgen) {
        pgen_reader.reset();
        builder.compute(pgen_reader);
    } else {
        builder.compute(plink_reader);
    }

    Eigen::MatrixXd K = builder.getGrm();
    center_matrix(K);
    scale_mean_diag(K);
    Eigen::MatrixXd A = K;
    const double s = comp_s_one(A, K, 1);

    // Leave-one-sample jackknife approximation of GEMMA Svar for the one-component case.
    std::vector<double> loo;
    loo.reserve(K.rows());
    for (int drop = 0; drop < K.rows(); ++drop) {
        Eigen::MatrixXd sub(K.rows() - 1, K.cols() - 1);
        int rr = 0;
        for (int i = 0; i < K.rows(); ++i) {
            if (i == drop) continue;
            int cc = 0;
            for (int j = 0; j < K.cols(); ++j) {
                if (j == drop) continue;
                sub(rr, cc++) = K(i, j);
            }
            ++rr;
        }
        center_matrix(sub);
        scale_mean_diag(sub);
        loo.push_back(comp_s_one(sub, sub, 1));
    }
    double mean = 0.0;
    for (double v : loo) mean += v;
    mean /= std::max<size_t>(1, loo.size());
    double svar = 0.0;
    for (double v : loo) svar += (v - mean) * (v - mean);
    if (loo.size() > 1) svar *= (static_cast<double>(loo.size() - 1) / static_cast<double>(loo.size()));

    std::ofstream sf(cfg_.out_prefix + ".S.txt");
    std::ofstream sizef(cfg_.out_prefix + ".size.txt");
    std::ofstream snpf(cfg_.out_prefix + ".snps.txt");
    if (!sf || !sizef || !snpf) throw std::runtime_error("Cannot write GEMMA-style S matrix outputs");
    sf << std::setprecision(10) << s << "\n" << svar << "\n";
    sizef << std::setprecision(10) << builder.getTotalSnps() << "\n" << K.rows() << "\n";
    if (is_bgen) {
        bgen_reader.loadVariantInfo();
        for (const auto& v : bgen_reader.getVariants()) snpf << v.rsid << "\n";
    } else if (is_pgen) {
        for (const auto& v : pgen_reader.getVariants()) snpf << v.id << "\n";
    } else {
        for (const auto& v : plink_reader.getBimInfo()) snpf << v.id << "\n";
    }

    std::cout << "GEMMA-style S outputs written to " << cfg_.out_prefix
              << ".S.txt, " << cfg_.out_prefix << ".size.txt, and "
              << cfg_.out_prefix << ".snps.txt\n";
}

void GemmaSummaryTools::runGq() {
    if (cfg_.gemma_q_mode != 1) {
        throw std::runtime_error("--gq currently supports GEMMA mode 1");
    }
    if (cfg_.beta_file.empty()) {
        throw std::runtime_error("--gq requires --beta/-beta <file>");
    }
    if (!cfg_.category_file.empty() || !cfg_.multi_category_file.empty() ||
        !cfg_.category_weight_file.empty() || !cfg_.snp_weight_summary_file.empty()) {
        CategorySpec spec = read_category_spec(cfg_);
        SummaryQMulti sq = compute_summary_q_multi(cfg_, spec);
        write_matrix_file(cfg_.out_prefix + ".Vq.txt", sq.vq);
        write_vector_file(cfg_.out_prefix + ".q.txt", sq.q);
        Eigen::VectorXd size(sq.s.size() + 1);
        for (int i = 0; i < sq.s.size(); ++i) size(i) = sq.s(i);
        size(sq.s.size()) = sq.ni_total;
        write_vector_file(cfg_.out_prefix + ".size.txt", size);
        std::cout << "GEMMA-style multi-category q outputs written to " << cfg_.out_prefix
                  << ".Vq.txt, " << cfg_.out_prefix << ".q.txt, and "
                  << cfg_.out_prefix << ".size.txt (" << sq.usable
                  << " usable rows; " << spec.ncat << " categories)\n";
        return;
    }
    SummaryQ sq = compute_summary_q(cfg_);

    std::ofstream vqf(cfg_.out_prefix + ".Vq.txt");
    std::ofstream qf(cfg_.out_prefix + ".q.txt");
    std::ofstream sizef(cfg_.out_prefix + ".size.txt");
    if (!vqf || !qf || !sizef) throw std::runtime_error("Cannot write GEMMA-style q outputs");
    vqf << std::setprecision(10) << sq.vq << "\n";
    qf << std::setprecision(10) << sq.q << "\n";
    sizef << std::setprecision(10) << sq.snp_weight_sum << "\n" << sq.ni_total << "\n";
    std::cout << "GEMMA-style q outputs written to " << cfg_.out_prefix
              << ".Vq.txt, " << cfg_.out_prefix << ".q.txt, and "
              << cfg_.out_prefix << ".size.txt (" << sq.usable
              << " usable rows; " << sq.total << " total rows)\n";
}

void GemmaSummaryTools::runSummaryVc() {
    PrefixSummaryData study = read_or_merge_studies(cfg_);
    PrefixSummaryData ref = read_or_merge_refs(cfg_);
    if (study.q.size() != ref.ref_s.size()) {
        throw std::runtime_error("summary-stat --vc study/ref dimensions do not match");
    }
    std::ostringstream extra;
    extra << "study_prefix\t" << cfg_.study_prefix << "\n";
    extra << "ref_prefix\t" << cfg_.ref_prefix << "\n";
    extra << "mstudy_file\t" << cfg_.multi_study_file << "\n";
    extra << "mref_file\t" << cfg_.multi_ref_file << "\n";
    extra << "n_ref\t" << ref.n_ref << "\n";
    extra << "ref_weight_sum\t" << ref.ref_s.sum() << "\n";
    write_summary_vc_outputs(cfg_, study.Vq, ref.S, ref.Svar, study.q, study.study_s, study.n_study,
                             cfg_.multi_study_file.empty() && cfg_.multi_ref_file.empty()
                                 ? "study_ref_multi_component"
                                 : "mstudy_mref_multi_component",
                             extra.str());

    std::cout << "GEMMA-style summary-stat VC output written to " << cfg_.out_prefix
              << ".summary.vc.txt and " << cfg_.out_prefix << ".summary.vc.meta.txt\n";
}

void GemmaSummaryTools::runBetaCorVc() {
    if (cfg_.beta_file.empty() || cfg_.cor_file.empty()) {
        throw std::runtime_error("summary-stat --vc with direct inputs requires --beta/-beta and --cor/-cor");
    }
    if (!cfg_.category_file.empty() || !cfg_.multi_category_file.empty() ||
        !cfg_.category_weight_file.empty() || !cfg_.snp_weight_summary_file.empty()) {
        CategorySpec spec = read_category_spec(cfg_);
        DirectSummaryMulti dm = compute_direct_multi(cfg_, spec);
        std::ostringstream extra;
        extra << "beta_file\t" << cfg_.beta_file << "\n";
        extra << "cor_file\t" << cfg_.cor_file << "\n";
        extra << "overlap_snp_count\t" << dm.overlap << "\n";
        extra << "ld_pair_count\t" << dm.ns_pair << "\n";
        extra << "category_file\t" << cfg_.category_file << "\n";
        extra << "mcat_file\t" << cfg_.multi_category_file << "\n";
        extra << "wcat_file\t" << cfg_.category_weight_file << "\n";
        extra << "wsnp_file\t" << cfg_.snp_weight_summary_file << "\n";
        write_summary_vc_outputs(cfg_, dm.vq, dm.S, dm.Svar, dm.q, dm.s,
                                 static_cast<double>(dm.ni_total),
                                 "beta_cor_multi_component", extra.str());
        Eigen::MatrixXd combined_s(dm.S.rows() * 2, dm.S.cols());
        combined_s.topRows(dm.S.rows()) = dm.S;
        combined_s.bottomRows(dm.S.rows()) = dm.Svar;
        write_matrix_file(cfg_.out_prefix + ".S.txt", combined_s);
        write_matrix_file(cfg_.out_prefix + ".Svar.txt", dm.Svar);
        write_matrix_file(cfg_.out_prefix + ".Vq.txt", dm.vq);
        write_vector_file(cfg_.out_prefix + ".q.txt", dm.q);
        Eigen::VectorXd size(dm.s.size() + 1);
        for (int i = 0; i < dm.s.size(); ++i) size(i) = dm.s(i);
        size(dm.s.size()) = dm.ni_total;
        write_vector_file(cfg_.out_prefix + ".size.txt", size);
        std::cout << "GEMMA-style multi-category beta/cor summary-stat VC output written to "
                  << cfg_.out_prefix << ".summary.vc.txt and "
                  << cfg_.out_prefix << ".summary.vc.meta.txt (" << dm.overlap
                  << " overlapping SNPs; " << dm.ns_pair << " LD pairs; "
                  << spec.ncat << " categories)\n";
        return;
    }
    CorSummary cor = read_cor_summary(cfg_.cor_file);
    SummaryQ sq = compute_summary_q_with_cor(cfg_, cor);
    SummaryS ss = compute_s_from_cor(cor);
    sq.vq += 2.0 * ss.qvar_ld;

    const double pve = sq.q / ss.s;
    const double sigma2_per_snp = pve / sq.snp_weight_sum;
    const double var_pve_component = ss.svar * pve * pve +
                                     sq.vq / (static_cast<double>(sq.ni_total) * static_cast<double>(sq.ni_total));
    const double se_pve = std::sqrt(std::max(0.0, var_pve_component)) / std::abs(ss.s);
    const double se_sigma2 = se_pve / sq.snp_weight_sum;

    std::ofstream out(cfg_.out_prefix + ".summary.vc.txt");
    if (!out) throw std::runtime_error("Cannot write beta/cor summary-stat --vc output");
    out << "component\tpve\tse_pve\tsigma2_per_snp\tse_sigma2_per_snp\tenrichment\tse_enrichment\n";
    out << std::scientific << std::setprecision(10)
        << "Vg\t" << pve << '\t' << se_pve << '\t'
        << sigma2_per_snp << '\t' << se_sigma2 << "\t1.0000000000e+00\tNA\n";
    out << "total_pve\t" << pve << '\t' << se_pve << "\tNA\tNA\tNA\tNA\n";

    std::ofstream meta(cfg_.out_prefix + ".summary.vc.meta.txt");
    if (!meta) throw std::runtime_error("Cannot write beta/cor summary-stat --vc metadata");
    meta << "method\tbeta_cor_single_component\n";
    meta << "beta_file\t" << cfg_.beta_file << "\n";
    meta << "cor_file\t" << cfg_.cor_file << "\n";
    meta << "n_study\t" << std::setprecision(10) << sq.ni_total << "\n";
    meta << "overlap_snp_count\t" << sq.usable << "\n";
    meta << "study_snp_weight_sum\t" << sq.snp_weight_sum << "\n";
    meta << "ld_pair_count\t" << ss.ns_pair << "\n";
    meta << "S\t" << ss.s << "\n";
    meta << "Svar\t" << ss.svar << "\n";
    meta << "q\t" << sq.q << "\n";
    meta << "Vq\t" << sq.vq << "\n";

    std::cout << "GEMMA-style beta/cor summary-stat VC output written to "
              << cfg_.out_prefix << ".summary.vc.txt and "
              << cfg_.out_prefix << ".summary.vc.meta.txt (" << sq.usable
              << " overlapping SNPs; " << ss.ns_pair << " LD pairs)\n";
}

void GemmaSummaryTools::runSummaryCi() {
    PrefixSummaryData study = read_or_merge_studies(cfg_);
    PrefixSummaryData ref = read_or_merge_refs(cfg_);
    if (study.q.size() != ref.ref_s.size()) {
        throw std::runtime_error("summary-stat --ci study/ref dimensions do not match");
    }
    const int nvc = static_cast<int>(study.q.size());
    Eigen::FullPivLU<Eigen::MatrixXd> lu(ref.S);
    if (!lu.isInvertible()) throw std::runtime_error("summary-stat --ci S matrix is singular");
    Eigen::MatrixXd Si = lu.inverse();
    Eigen::VectorXd pve = Si * study.q;
    Eigen::MatrixXd var_mat = Eigen::MatrixXd::Zero(nvc, nvc);
    for (int i = 0; i < nvc; ++i) {
        for (int j = 0; j < nvc; ++j) {
            var_mat(i, j) = ref.Svar(i, j) * pve(i) * pve(j) +
                            study.Vq(i, j) / (study.n_study * study.n_study);
        }
    }
    Eigen::MatrixXd cov_pve = Si * var_mat * Si.transpose();

    std::ofstream out(cfg_.out_prefix + ".summary.ci.txt");
    if (!out) throw std::runtime_error("Cannot write summary-stat --ci output");
    out << "component\tpve\tse_pve\tci95_lower\tci95_upper\n";
    out << std::scientific << std::setprecision(10);
    for (int i = 0; i < nvc; ++i) {
        const double se_pve = std::sqrt(std::max(0.0, cov_pve(i, i)));
        const double lo = pve(i) - 1.96 * se_pve;
        const double hi = pve(i) + 1.96 * se_pve;
        out << "Vg" << (i + 1) << '\t' << pve(i) << '\t' << se_pve << '\t' << lo << '\t' << hi << "\n";
    }
    const double total = pve.sum();
    const double total_se = std::sqrt(std::max(0.0, cov_pve.sum()));
    out << "total_pve\t" << total << '\t' << total_se << '\t'
        << total - 1.96 * total_se << '\t' << total + 1.96 * total_se << "\n";
    std::cout << "GEMMA-style summary-stat CI output written to "
              << cfg_.out_prefix << ".summary.ci.txt\n";
}

} // namespace cosmic
