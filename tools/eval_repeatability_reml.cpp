#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

std::vector<std::string> split_ws(const std::string& s) {
    std::stringstream ss(s);
    std::vector<std::string> out;
    std::string x;
    while (ss >> x) out.push_back(x);
    return out;
}

std::string id_key(const std::string& fid, const std::string& iid) {
    return fid + ":" + iid;
}

std::string record_key(const std::string& fid, const std::string& iid, const std::string& time) {
    return fid + ":" + iid + ":" + time;
}

int require_col(const std::vector<std::string>& header, const std::string& name, const std::string& file) {
    for (int i = 0; i < (int)header.size(); ++i) {
        if (header[i] == name) return i;
    }
    throw std::runtime_error("Missing column '" + name + "' in " + file);
}

struct BedData {
    MatrixXd X;
    std::vector<std::string> ids;
};

BedData read_bed(const std::string& prefix) {
    std::ifstream fam(prefix + ".fam");
    if (!fam) throw std::runtime_error("Cannot open FAM: " + prefix + ".fam");
    std::string line;
    std::vector<std::string> ids;
    while (std::getline(fam, line)) {
        auto v = split_ws(line);
        if (v.size() >= 2) ids.push_back(id_key(v[0], v[1]));
    }
    const int n = (int)ids.size();

    std::ifstream bim(prefix + ".bim");
    if (!bim) throw std::runtime_error("Cannot open BIM: " + prefix + ".bim");
    int m = 0;
    while (std::getline(bim, line)) ++m;

    std::ifstream bed(prefix + ".bed", std::ios::binary);
    if (!bed) throw std::runtime_error("Cannot open BED: " + prefix + ".bed");
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(bed)), {});
    if (bytes.size() < 3 || bytes[0] != 0x6c || bytes[1] != 0x1b || bytes[2] != 0x01) {
        throw std::runtime_error("Unsupported BED file (expected SNP-major PLINK BED)");
    }

    MatrixXd X(n, m);
    const int bytes_per_snp = (n + 3) / 4;
    size_t pos = 3;
    for (int j = 0; j < m; ++j) {
        std::vector<double> col(n, std::numeric_limits<double>::quiet_NaN());
        for (int b = 0; b < bytes_per_snp; ++b) {
            unsigned char c = bytes[pos++];
            for (int k = 0; k < 4; ++k) {
                int i = b * 4 + k;
                if (i >= n) break;
                int lo = (c >> (2 * k)) & 1;
                int hi = (c >> (2 * k + 1)) & 1;
                if (lo == 0 && hi == 0) col[i] = 0.0;
                else if (lo == 1 && hi == 1) col[i] = 2.0;
                else if (lo == 0 && hi == 1) col[i] = 1.0;
            }
        }
        double sum = 0.0;
        int cnt = 0;
        for (double x : col) {
            if (std::isfinite(x)) { sum += x; ++cnt; }
        }
        double mean = cnt > 0 ? sum / cnt : 0.0;
        for (int i = 0; i < n; ++i) X(i, j) = std::isfinite(col[i]) ? col[i] : mean;
    }
    return {X, ids};
}

MatrixXd make_vanraden_grm(const MatrixXd& X) {
    const int n = X.rows();
    const int m = X.cols();
    std::vector<int> keep;
    std::vector<double> p(m);
    double denom = 0.0;
    for (int j = 0; j < m; ++j) {
        p[j] = X.col(j).mean() / 2.0;
        if (p[j] > 0.0 && p[j] < 1.0) {
            keep.push_back(j);
            denom += 2.0 * p[j] * (1.0 - p[j]);
        }
    }
    if (keep.empty() || denom <= 0.0) throw std::runtime_error("No valid SNPs for GRM");
    MatrixXd M(n, (int)keep.size());
    for (int jj = 0; jj < (int)keep.size(); ++jj) {
        int j = keep[jj];
        M.col(jj) = X.col(j).array() - 2.0 * p[j];
    }
    return (M * M.transpose()) / denom;
}

struct ModelData {
    VectorXd y;
    MatrixXd X;
    std::vector<int> obs_to_subject;
};

ModelData read_model_data(const std::string& pheno_file,
                          const std::string& covar_file,
                          const std::string& trait,
                          const std::string& covar_name,
                          const std::vector<std::string>& ids) {
    std::map<std::string, std::string> covar_by_record;
    std::ifstream cf(covar_file);
    if (!cf) throw std::runtime_error("Cannot open covar: " + covar_file);
    std::string line;
    if (!std::getline(cf, line)) throw std::runtime_error("Empty covar file: " + covar_file);
    auto ch = split_ws(line);
    int c_fid = require_col(ch, "FID", covar_file);
    int c_iid = require_col(ch, "IID", covar_file);
    int c_time = require_col(ch, "TIME", covar_file);
    int c_cov = require_col(ch, covar_name, covar_file);
    std::set<std::string> levels_set;
    while (std::getline(cf, line)) {
        auto v = split_ws(line);
        if ((int)v.size() <= std::max({c_fid, c_iid, c_time, c_cov})) continue;
        std::string key = record_key(v[c_fid], v[c_iid], v[c_time]);
        covar_by_record[key] = v[c_cov];
        levels_set.insert(v[c_cov]);
    }
    std::vector<std::string> levels(levels_set.begin(), levels_set.end());

    std::map<std::string, int> id_index;
    for (int i = 0; i < (int)ids.size(); ++i) id_index[ids[i]] = i;

    std::ifstream pf(pheno_file);
    if (!pf) throw std::runtime_error("Cannot open pheno: " + pheno_file);
    if (!std::getline(pf, line)) throw std::runtime_error("Empty pheno file: " + pheno_file);
    auto ph = split_ws(line);
    int p_fid = require_col(ph, "FID", pheno_file);
    int p_iid = require_col(ph, "IID", pheno_file);
    int p_time = require_col(ph, "TIME", pheno_file);
    int p_trait = require_col(ph, trait, pheno_file);

    std::vector<double> yv;
    std::vector<int> subj;
    std::vector<std::string> cov_levels;
    while (std::getline(pf, line)) {
        auto v = split_ws(line);
        if ((int)v.size() <= std::max({p_fid, p_iid, p_time, p_trait})) continue;
        auto id_it = id_index.find(id_key(v[p_fid], v[p_iid]));
        if (id_it == id_index.end()) continue;
        auto cov_it = covar_by_record.find(record_key(v[p_fid], v[p_iid], v[p_time]));
        if (cov_it == covar_by_record.end()) continue;
        yv.push_back(std::stod(v[p_trait]));
        subj.push_back(id_it->second);
        cov_levels.push_back(cov_it->second);
    }
    if (yv.empty()) throw std::runtime_error("No matched phenotype/covariate records");

    const int n = (int)yv.size();
    const int p = 1 + std::max(0, (int)levels.size() - 1);
    VectorXd y(n);
    MatrixXd X = MatrixXd::Zero(n, p);
    for (int i = 0; i < n; ++i) {
        y[i] = yv[i];
        X(i, 0) = 1.0;
        for (int l = 1; l < (int)levels.size(); ++l) {
            if (cov_levels[i] == levels[l]) X(i, l) = 1.0;
        }
    }
    return {y, X, subj};
}

double exact_reml_logl(const VectorXd& y, const MatrixXd& X, const std::vector<int>& obs_to_subject,
                       const MatrixXd& G, double vg, double vpe, double ve) {
    if (!(vg >= 0.0) || !(vpe >= 0.0) || !(ve > 0.0)) return -std::numeric_limits<double>::infinity();
    const int n = y.size();
    const int p = X.cols();
    MatrixXd V(n, n);
    for (int i = 0; i < n; ++i) {
        int si = obs_to_subject[i];
        for (int j = 0; j <= i; ++j) {
            int sj = obs_to_subject[j];
            double value = vg * G(si, sj);
            if (si == sj) value += vpe;
            if (i == j) value += ve;
            V(i, j) = V(j, i) = value;
        }
    }
    Eigen::LDLT<MatrixXd> v_ldlt(V);
    if (v_ldlt.info() != Eigen::Success || (v_ldlt.vectorD().array() <= 0.0).any()) {
        return -std::numeric_limits<double>::infinity();
    }
    double logdet_v = v_ldlt.vectorD().array().log().sum();
    MatrixXd VinvX = v_ldlt.solve(X);
    VectorXd Vinvy = v_ldlt.solve(y);
    MatrixXd XtVinvX = X.transpose() * VinvX;
    Eigen::LDLT<MatrixXd> x_ldlt(XtVinvX);
    if (x_ldlt.info() != Eigen::Success || (x_ldlt.vectorD().array() <= 0.0).any()) {
        return -std::numeric_limits<double>::infinity();
    }
    double logdet_x = x_ldlt.vectorD().array().log().sum();
    VectorXd beta = x_ldlt.solve(X.transpose() * Vinvy);
    VectorXd Py = Vinvy - VinvX * beta;
    double yPy = y.dot(Py);
    return -0.5 * (logdet_v + logdet_x + yPy + (n - p) * std::log(2.0 * M_PI));
}

std::map<std::string, double> read_null_file(const std::string& file) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("Cannot open null file: " + file);
    std::string line;
    std::getline(in, line);
    std::map<std::string, double> values;
    while (std::getline(in, line)) {
        auto v = split_ws(line);
        if (v.size() >= 2) {
            try { values[v[0]] = std::stod(v[1]); } catch (...) {}
        }
    }
    return values;
}

void usage() {
    std::cerr << "Usage: eval_repeatability_reml --bfile PREFIX --pheno FILE --covar FILE "
              << "[--trait TNB] [--covar-name PARITY] (--vg V --vpe V --ve V | --null FILE) [--label NAME]\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::string bfile, pheno, covar, null_file, label = "input";
        std::string trait = "TNB", covar_name = "PARITY";
        double vg = std::numeric_limits<double>::quiet_NaN();
        double vpe = std::numeric_limits<double>::quiet_NaN();
        double ve = std::numeric_limits<double>::quiet_NaN();
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            auto next = [&]() -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("Missing value for " + a);
                return argv[++i];
            };
            if (a == "--bfile") bfile = next();
            else if (a == "--pheno") pheno = next();
            else if (a == "--covar") covar = next();
            else if (a == "--trait") trait = next();
            else if (a == "--covar-name") covar_name = next();
            else if (a == "--null") null_file = next();
            else if (a == "--label") label = next();
            else if (a == "--vg") vg = std::stod(next());
            else if (a == "--vpe") vpe = std::stod(next());
            else if (a == "--ve") ve = std::stod(next());
            else if (a == "--help" || a == "-h") { usage(); return 0; }
            else throw std::runtime_error("Unknown argument: " + a);
        }
        if (bfile.empty() || pheno.empty() || covar.empty()) { usage(); return 2; }
        if (!null_file.empty()) {
            auto vals = read_null_file(null_file);
            if (vals.count("Vg")) vg = vals["Vg"];
            if (vals.count("Vpe")) vpe = vals["Vpe"];
            if (vals.count("Ve")) ve = vals["Ve"];
            if (label == "input") label = null_file;
        }
        if (!std::isfinite(vg) || !std::isfinite(vpe) || !std::isfinite(ve)) {
            throw std::runtime_error("Provide --vg/--vpe/--ve or --null FILE containing Vg/Vpe/Ve");
        }

        BedData bed = read_bed(bfile);
        MatrixXd G = make_vanraden_grm(bed.X);
        ModelData data = read_model_data(pheno, covar, trait, covar_name, bed.ids);
        double logl = exact_reml_logl(data.y, data.X, data.obs_to_subject, G, vg, vpe, ve);
        double vp = vg + vpe + ve;

        std::cout << "label\tn_records\tn_subjects\tn_fixed\tlogL\tVg\tVpe\tVe\th2\tr2\n";
        std::cout << label << '\t' << data.y.size() << '\t' << bed.ids.size() << '\t' << data.X.cols() << '\t'
                  << logl << '\t' << vg << '\t' << vpe << '\t' << ve << '\t'
                  << (vg / vp) << '\t' << ((vg + vpe) / vp) << '\n';
        return std::isfinite(logl) ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "eval_repeatability_reml error: " << e.what() << "\n";
        return 1;
    }
}
