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

void OrdinalModel::run() {
    auto start = std::chrono::steady_clock::now();
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

        std::cout << "Ordinal cumulative threshold " << (k + 1) << "/" << (levels_.size() - 1)
                  << " (y <= " << threshold << ") using experimental cumulative-binary PSR.\n";
        ThresholdModel model(opts);
        model.prepare(fids_, iids_, n_samples_, *plink_reader_, *bgen_reader_, *pgen_reader_, is_bgen_, is_pgen_);
        model.run();
    }
    auto end = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration<double>(end - start).count();
    writeCombinedConvergence(threshold_prefixes);
    writeCombinedGwas(threshold_prefixes, runtime);
    std::cout << "Ordinal PSR-GWAS experimental cumulative-binary results written to "
              << options_.out_prefix << ".gwas.tsv\n";
}

} // namespace cosmic
