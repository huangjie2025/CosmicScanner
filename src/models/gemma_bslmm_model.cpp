#include "models/gemma_bslmm_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <unordered_set>

#include "data_loader.h"

namespace cosmic {

namespace {

std::unordered_set<std::string> load_id_set(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (in >> id) ids.insert(id);
    return ids;
}

double sigmoid(double x) {
    if (x >= 0.0) {
        double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    double z = std::exp(x);
    return z / (1.0 + z);
}

std::vector<std::string> split_ws(const std::string& line) {
    std::stringstream ss(line);
    std::vector<std::string> out;
    std::string x;
    while (ss >> x) out.push_back(x);
    return out;
}

void write_bslmm_meta(const std::string& out_prefix,
                      int mode,
                      const std::string& implementation,
                      int n_samples,
                      int n_snps,
                      double h2,
                      double h_min,
                      double h_max,
                      double rho_min,
                      double rho_max,
                      double p_log_min,
                      double p_log_max,
                      double gamma_mean,
                      int s_min,
                      int s_max,
                      int burnin,
                      int sampling,
                      int rpace,
                      int wpace,
                      long long seed) {
    std::ofstream meta(out_prefix + ".bslmm.meta.txt");
    if (!meta) throw std::runtime_error("Cannot write BSLMM metadata");
    meta << std::setprecision(12);
    meta << "mode\t" << mode << "\n";
    meta << "implementation\t" << implementation << "\n";
    meta << "n_samples\t" << n_samples << "\n";
    meta << "n_snps\t" << n_snps << "\n";
    meta << "h2\t" << h2 << "\n";
    meta << "h_min\t" << h_min << "\n";
    meta << "h_max\t" << h_max << "\n";
    meta << "rho_min\t" << rho_min << "\n";
    meta << "rho_max\t" << rho_max << "\n";
    meta << "pmin_log10_pi\t" << p_log_min << "\n";
    meta << "pmax_log10_pi\t" << p_log_max << "\n";
    meta << "pi_min\t" << std::pow(10.0, p_log_min) << "\n";
    meta << "pi_max\t" << std::pow(10.0, p_log_max) << "\n";
    meta << "gamma_mean\t" << gamma_mean << "\n";
    meta << "s_min\t" << s_min << "\n";
    meta << "s_max\t" << s_max << "\n";
    meta << "burnin\t" << burnin << "\n";
    meta << "sampling\t" << sampling << "\n";
    meta << "rpace\t" << rpace << "\n";
    meta << "wpace\t" << wpace << "\n";
    meta << "seed\t" << seed << "\n";
}

} // namespace

GemmaBslmmModel::Options GemmaBslmmModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.mpheno = cfg.mpheno;
    opts.mode = cfg.bslmm_mode;
    opts.top_sparse = cfg.bslmm_top;
    opts.h2 = cfg.bslmm_h2;
    opts.min_maf = cfg.min_maf;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.bf_file = cfg.bslmm_bf_file;
    opts.hyp_file = cfg.bslmm_hyp_file;
    opts.h_grid = cfg.bslmm_hgrid;
    opts.rho_grid = cfg.bslmm_rho_grid;
    opts.h_min = cfg.bslmm_h_min;
    opts.h_max = cfg.bslmm_h_max;
    opts.rho_min = cfg.bslmm_rho_min;
    opts.rho_max = cfg.bslmm_rho_max;
    opts.p_min = cfg.bslmm_p_min;
    opts.p_max = cfg.bslmm_p_max;
    opts.p_min_set = cfg.bslmm_p_min_set;
    opts.p_max_set = cfg.bslmm_p_max_set;
    opts.gamma_mean = cfg.bslmm_gamma_mean;
    opts.h_scale = cfg.bslmm_h_scale;
    opts.rho_scale = cfg.bslmm_rho_scale;
    opts.p_scale = cfg.bslmm_p_scale;
    opts.mh_steps = cfg.bslmm_mh_steps;
    opts.s_min = cfg.bslmm_s_min;
    opts.rpace = cfg.bslmm_rpace;
    opts.wpace = cfg.bslmm_wpace;
    opts.burnin = cfg.bslmm_burnin;
    opts.sampling = cfg.bslmm_sampling;
    opts.seed = cfg.bslmm_seed;
    return opts;
}

GemmaBslmmModel::GemmaBslmmModel(const Options& options) : options_(options) {}

void GemmaBslmmModel::prepare(const std::vector<std::string>& fids,
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

    if (options_.mode < 1 || options_.mode > 5) {
        throw std::runtime_error("--bslmm supports modes 1, 2, 3, 4, and 5");
    }
    if (options_.mode == 5) {
        if (options_.bf_file.empty()) throw std::runtime_error("--bslmm 5 requires --bf/-bf");
        if (options_.hyp_file.empty()) throw std::runtime_error("--bslmm 5 requires --hyp/-hyp");
        return;
    }
    if (options_.pheno_file.empty()) throw std::runtime_error("--bslmm requires --pheno");
    if (!(options_.h2 > 0.0 && options_.h2 < 1.0)) {
        throw std::runtime_error("--bslmm-h2 must be in (0,1)");
    }
    if (options_.mode == 4 && (options_.h_grid <= 0 || options_.rho_grid <= 0)) {
        throw std::runtime_error("--bslmm 4 requires positive --bslmm-hgrid and --bslmm-rho-grid");
    }
    if (!(options_.h_min > 0.0 && options_.h_max < 1.0 && options_.h_min <= options_.h_max)) {
        throw std::runtime_error("BSLMM/DAP requires 0 < -hmin <= -hmax < 1");
    }
    if (!(options_.rho_min > 0.0 && options_.rho_max < 1.0 && options_.rho_min <= options_.rho_max)) {
        throw std::runtime_error("BSLMM/DAP requires 0 < -rmin <= -rmax < 1");
    }
    if (!(options_.p_min <= options_.p_max && options_.p_max <= 0.0)) {
        throw std::runtime_error("BSLMM/DAP requires GEMMA-style -pmin <= -pmax <= 0, where p is log10(pi)");
    }
    if (!(options_.h_scale > 0.0 && options_.rho_scale > 0.0 && options_.p_scale > 0.0)) {
        throw std::runtime_error("BSLMM proposal scales -hscale/-rscale/-pscale must be positive");
    }
    if (options_.s_min < 0 || options_.top_sparse < options_.s_min) {
        throw std::runtime_error("BSLMM requires 0 <= -smin <= -smax/--bslmm-top");
    }

    y_ = DataLoader::loadPhenotype(options_.pheno_file, options_.mpheno, fids_, iids_, n_samples_,
                                   *plink_reader_, is_bgen_, is_pgen_);
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for --bslmm");
    if (options_.mode == 3) {
        for (int i = 0; i < keep_.n_keep; ++i) {
            double v = y_[keep_.keep_indices[i]];
            if (!(std::abs(v) < 1e-8 || std::abs(v - 1.0) < 1e-8)) {
                throw std::runtime_error("--bslmm 3 requires 0/1 phenotypes");
            }
        }
    }
}

void GemmaBslmmModel::loadDesign(Eigen::MatrixXd& Z, std::vector<VariantMeta>& variants) const {
    const auto extract = load_id_set(options_.extract_file);
    const auto exclude = load_id_set(options_.exclude_file);
    std::vector<Eigen::VectorXd> cols;
    const int n = keep_.n_keep;

    auto add_variant = [&](const VariantMeta& meta, const Eigen::VectorXf& dosage_f) {
        if (!extract.empty() && !extract.count(meta.id)) return;
        if (!exclude.empty() && exclude.count(meta.id)) return;
        Eigen::VectorXd g(n);
        double sum = 0.0;
        int valid = 0;
        for (int i = 0; i < n; ++i) {
            double v = dosage_f[keep_.keep_indices[i]];
            g[i] = v;
            if (std::isfinite(v)) {
                sum += v;
                ++valid;
            }
        }
        if (valid == 0) return;
        VariantMeta m = meta;
        m.mean = sum / valid;
        m.freq = m.mean / 2.0;
        double maf = std::min(m.freq, 1.0 - m.freq);
        if (maf < options_.min_maf) return;
        for (int i = 0; i < n; ++i) if (!std::isfinite(g[i])) g[i] = m.mean;
        g.array() -= m.mean;
        double ss = g.squaredNorm();
        if (!(ss > 1e-12)) return;
        m.sd = std::sqrt(ss / std::max(1, n - 1));
        g.array() /= m.sd;
        variants.push_back(m);
        cols.push_back(std::move(g));
    };

    if (is_pgen_) {
        pgen_reader_->reset();
        Eigen::VectorXf dosage;
        std::string chrom, rsid, ref, alt;
        uint32_t pos = 0;
        while (pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage)) {
            add_variant({chrom, rsid, (long long)pos, alt, ref}, dosage);
        }
    } else if (is_bgen_) {
        bgen_reader_->reset();
        Eigen::VectorXf dosage;
        std::string chrom, rsid, a1, a2;
        uint32_t pos = 0;
        while (bgen_reader_->readVariant(chrom, rsid, pos, a1, a2, dosage)) {
            add_variant({chrom, rsid, (long long)pos, a1, a2}, dosage);
        }
    } else {
        const auto& bim = plink_reader_->getBimInfo();
        const size_t block_size = 512;
        for (size_t start = 0; start < bim.size(); start += block_size) {
            size_t end = std::min(bim.size(), start + block_size);
            Eigen::MatrixXf raw;
            plink_reader_->getSnpBlock(start, end, raw);
            for (size_t k = 0; k < end - start; ++k) {
                const auto& snp = bim[start + k];
                Eigen::VectorXf dosage = raw.col((int)k);
                add_variant({snp.chrom, snp.id, snp.bp_pos, snp.alt, snp.ref}, dosage);
            }
        }
    }

    Z.resize(n, (int)cols.size());
    for (int j = 0; j < (int)cols.size(); ++j) Z.col(j) = cols[j];
}

void GemmaBslmmModel::run() {
    if (options_.mode == 5) {
        runDapFineMap();
        return;
    }
    Eigen::MatrixXd Z;
    std::vector<VariantMeta> variants;
    loadDesign(Z, variants);
    const int n = Z.rows();
    const int m = Z.cols();
    if (m == 0) throw std::runtime_error("--bslmm found no usable SNPs after filters");
    const double p_log_min = options_.p_min_set ? options_.p_min : -std::log10(std::max(1, m));
    const double p_log_max = options_.p_max_set ? options_.p_max : 0.0;
    if (!(p_log_min <= p_log_max && p_log_max <= 0.0)) {
        throw std::runtime_error("BSLMM/DAP requires GEMMA-style -pmin <= -pmax <= 0 after resolving defaults");
    }

    Eigen::VectorXd y(n);
    Eigen::MatrixXd C(n, X_.cols());
    for (int i = 0; i < n; ++i) {
        int idx = keep_.keep_indices[i];
        y[i] = y_[idx];
        for (int j = 0; j < X_.cols(); ++j) C(i, j) = X_(idx, j);
    }
    const double y_mean = y.mean();
    Eigen::VectorXd y_work = y;
    if (C.cols() > 0) {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(C.transpose() * C);
        Eigen::VectorXd b = ldlt.solve(C.transpose() * y);
        y_work = y - C * b;
    }
    if (options_.mode == 3) {
        double p = std::min(0.999, std::max(0.001, y.mean()));
        y_work = y.array() - p;
    } else {
        y_work.array() -= y_work.mean();
    }

    const double lambda = (double)m * (1.0 - options_.h2) / options_.h2;
    Eigen::MatrixXd gram = Z.transpose() * Z;
    gram.diagonal().array() += lambda;
    Eigen::VectorXd beta_std = gram.ldlt().solve(Z.transpose() * y_work);
    Eigen::VectorXd alpha_raw(m);
    for (int j = 0; j < m; ++j) alpha_raw[j] = beta_std[j] / variants[j].sd;

    if (options_.mode == 4) {
        const double y_var = std::max(1e-12, y_work.squaredNorm() / std::max(1, n - 1));
        std::vector<std::array<double, 5>> hyper;
        hyper.reserve(options_.h_grid * options_.rho_grid);
        const int n_grid = options_.h_grid * options_.rho_grid;
        for (int hi = 0; hi < options_.h_grid; ++hi) {
            const double h = options_.h_grid == 1 ? options_.h2
                : options_.h_min + (options_.h_max - options_.h_min) * static_cast<double>(hi) / static_cast<double>(options_.h_grid - 1);
            for (int ri = 0; ri < options_.rho_grid; ++ri) {
                const double rho = options_.rho_grid == 1 ? 0.5
                    : options_.rho_min + (options_.rho_max - options_.rho_min) * static_cast<double>(ri) / static_cast<double>(options_.rho_grid - 1);
                const double pi = std::pow(10.0, 0.5 * (p_log_min + p_log_max));
                const double expected_sparse = std::max(1.0, options_.gamma_mean > 0.0 ? options_.gamma_mean : pi * static_cast<double>(m));
                const double sa2 = h * rho / ((1.0 - h) * expected_sparse);
                const double sb2 = h * (1.0 - rho) / ((1.0 - h) * static_cast<double>(m));
                hyper.push_back({h, rho, sa2, sb2, 1.0 / static_cast<double>(n_grid)});
            }
        }

        Eigen::VectorXd zy = Z.transpose() * y_work;
        std::ofstream hyp(options_.out_prefix + ".hyp.txt");
        std::ofstream bf(options_.out_prefix + ".bf.txt");
        if (!hyp || !bf) throw std::runtime_error("Cannot write DAP outputs");
        hyp << "h\trho\tsa2\tsb2\tweight\n";
        hyp << std::scientific << std::setprecision(6);
        for (const auto& h : hyper) hyp << h[0] << '\t' << h[1] << '\t' << h[2] << '\t' << h[3] << '\t' << h[4] << '\n';
        bf << "chr\trs\tps\tn_miss";
        for (int k = 0; k < static_cast<int>(hyper.size()); ++k) bf << "\tBF" << (k + 1);
        bf << "\n" << std::scientific << std::setprecision(6);
        for (int j = 0; j < m; ++j) {
            const double zscore = zy[j] / std::sqrt(std::max(1e-12, y_var * Z.col(j).squaredNorm()));
            const double z2 = zscore * zscore;
            bf << variants[j].chrom << '\t' << variants[j].id << '\t' << variants[j].pos << "\t0";
            for (const auto& h : hyper) {
                const double prior_var = std::max(1e-10, h[2] + h[3]);
                const double r = prior_var / (prior_var + 1.0 / std::max(1.0, Z.col(j).squaredNorm()));
                const double abf = std::sqrt(std::max(1e-12, 1.0 - r)) * std::exp(0.5 * r * z2);
                bf << '\t' << abf;
            }
            bf << '\n';
        }
        std::cout << "GEMMA-style BSLMM-DAP mode 4 complete. Outputs: "
                  << options_.out_prefix << ".bf.txt and .hyp.txt\n"
                  << "  Note: this is a deterministic Wakefield-style ABF approximation, not GEMMA's exact GSL DAP engine.\n";
        write_bslmm_meta(options_.out_prefix, options_.mode,
                         "deterministic-wakefield-abf-dap-approximation",
                         n, m, options_.h2, options_.h_min, options_.h_max,
                         options_.rho_min, options_.rho_max, p_log_min, p_log_max,
                         options_.gamma_mean, options_.s_min, options_.top_sparse,
                         options_.burnin, options_.sampling, options_.rpace, options_.wpace,
                         options_.seed);
        return;
    }

    std::vector<int> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return std::abs(beta_std[a]) > std::abs(beta_std[b]);
    });
    int top = (options_.mode == 2) ? 0 : std::min(std::max(options_.s_min, options_.top_sparse), m);
    std::vector<char> gamma(m, 0);
    for (int i = 0; i < top; ++i) gamma[order[i]] = 1;

    Eigen::VectorXd beta_sparse = Eigen::VectorXd::Zero(m);
    if (top > 0) {
        Eigen::MatrixXd Zs(n, top);
        for (int i = 0; i < top; ++i) Zs.col(i) = Z.col(order[i]);
        Eigen::MatrixXd gs = Zs.transpose() * Zs;
        gs.diagonal().array() += std::max(1e-8, 0.1 * lambda);
        Eigen::VectorXd bs = gs.ldlt().solve(Zs.transpose() * y_work);
        for (int i = 0; i < top; ++i) beta_sparse[order[i]] = bs[i] / variants[order[i]].sd;
    }

    Eigen::VectorXd pred = Eigen::VectorXd::Constant(n, options_.mode == 3 ? 0.0 : y_mean);
    for (int j = 0; j < m; ++j) pred += Z.col(j) * beta_std[j];
    if (top > 0) {
        for (int j = 0; j < m; ++j) pred += Z.col(j) * (beta_sparse[j] * variants[j].sd * 0.25);
    }
    if (options_.mode == 3) {
        for (int i = 0; i < n; ++i) pred[i] = sigmoid(pred[i]);
    }

    {
        std::ofstream out(options_.out_prefix + ".param.txt");
        if (!out) throw std::runtime_error("Cannot write BSLMM param output");
        out << "chr\trs\tps\tn_miss\talpha\tbeta\tgamma\n";
        out << std::scientific << std::setprecision(6);
        for (int j = 0; j < m; ++j) {
            out << variants[j].chrom << '\t' << variants[j].id << '\t' << variants[j].pos
                << "\t0\t" << alpha_raw[j] << '\t' << beta_sparse[j] << '\t' << (int)gamma[j] << '\n';
        }
    }
    {
        std::ofstream out(options_.out_prefix + ".bv.txt");
        if (!out) throw std::runtime_error("Cannot write BSLMM bv output");
        out << std::scientific << std::setprecision(6);
        int k = 0;
        for (int i = 0; i < n_samples_; ++i) {
            if (k < keep_.n_keep && keep_.keep_indices[k] == i) out << pred[k++] << '\n';
            else out << "NA\n";
        }
    }
    {
        std::ofstream out(options_.out_prefix + ".hyp.txt");
        if (!out) throw std::runtime_error("Cannot write BSLMM hyp output");
        double var_y = (y.array() - y.mean()).square().sum() / std::max(1, n - 1);
        double var_p = (pred.array() - pred.mean()).square().sum() / std::max(1, n - 1);
        double pve = var_y > 0.0 ? std::min(1.0, std::max(0.0, var_p / var_y)) : 0.0;
        out << "h\tpve\trho\tpge\tpi\tn_gamma\n";
        const double pi_min = std::pow(10.0, p_log_min);
        const double pi_max = std::pow(10.0, p_log_max);
        const double pi = top > 0 ? std::min(pi_max, std::max(pi_min, (double)top / (double)m)) : 0.0;
        out << std::scientific << std::setprecision(6)
            << options_.h2 << '\t' << pve << '\t'
            << (options_.mode == 2 ? 1.0 : 0.5) << '\t'
            << (top > 0 ? 0.5 : 0.0) << '\t'
            << pi << '\t'
            << top << '\n';
    }

    std::cout << "GEMMA-style BSLMM mode " << options_.mode << " complete. Used "
              << m << " SNPs and " << n << " individuals. Outputs: "
              << options_.out_prefix << ".param.txt, .bv.txt, .hyp.txt\n";
    if (options_.mode != 2) {
        std::cout << "  Note: mode " << options_.mode
                  << " is a deterministic sparse+ridge approximation, not GEMMA's full MCMC sampler.\n";
    }
    write_bslmm_meta(options_.out_prefix, options_.mode,
                     options_.mode == 2 ? "deterministic-ridge-gblup" : "deterministic-sparse-ridge-approximation",
                     n, m, options_.h2, options_.h_min, options_.h_max,
                     options_.rho_min, options_.rho_max, p_log_min, p_log_max,
                     options_.gamma_mean, options_.s_min, options_.top_sparse,
                     options_.burnin, options_.sampling, options_.rpace, options_.wpace,
                     options_.seed);
}

void GemmaBslmmModel::runDapFineMap() const {
    if (!(options_.p_min <= options_.p_max && options_.p_max <= 0.0)) {
        throw std::runtime_error("BSLMM/DAP requires GEMMA-style -pmin <= -pmax <= 0, where p is log10(pi)");
    }
    std::ifstream hyp_in(options_.hyp_file);
    if (!hyp_in) throw std::runtime_error("Cannot open --hyp file: " + options_.hyp_file);
    std::string line;
    if (!std::getline(hyp_in, line)) throw std::runtime_error("Empty --hyp file: " + options_.hyp_file);
    std::vector<double> hyp_weight;
    std::vector<std::array<double, 5>> hyper;
    while (std::getline(hyp_in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.size() < 5) continue;
        std::array<double, 5> h{};
        for (int i = 0; i < 5; ++i) h[i] = std::stod(f[i]);
        hyper.push_back(h);
        hyp_weight.push_back(std::max(0.0, h[4]));
    }
    if (hyper.empty()) throw std::runtime_error("--hyp file has no usable hyperparameter rows");
    double wsum = std::accumulate(hyp_weight.begin(), hyp_weight.end(), 0.0);
    if (wsum <= 0.0) {
        std::fill(hyp_weight.begin(), hyp_weight.end(), 1.0 / static_cast<double>(hyp_weight.size()));
    } else {
        for (double& w : hyp_weight) w /= wsum;
    }

    std::ifstream bf_in(options_.bf_file);
    if (!bf_in) throw std::runtime_error("Cannot open --bf file: " + options_.bf_file);
    if (!std::getline(bf_in, line)) throw std::runtime_error("Empty --bf file: " + options_.bf_file);
    auto header = split_ws(line);
    bool gemma_bf = header.size() >= 5 && (header[0] == "chr" || header[0] == "CHR");
    std::vector<std::string> rsids;
    std::vector<double> scores;
    while (std::getline(bf_in, line)) {
        if (line.empty()) continue;
        auto f = split_ws(line);
        if (f.empty()) continue;
        const size_t rs_col = gemma_bf ? 1 : 0;
        const size_t bf_start = gemma_bf ? 4 : 2;
        if (f.size() <= bf_start) continue;
        rsids.push_back(f[rs_col]);
        double score = 0.0;
        const size_t nbf = f.size() - bf_start;
        for (size_t k = 0; k < nbf; ++k) {
            const double bf = std::max(0.0, std::stod(f[bf_start + k]));
            const double w = k < hyp_weight.size() ? hyp_weight[k] : (1.0 / static_cast<double>(nbf));
            score += w * bf;
        }
        scores.push_back(score);
    }
    if (rsids.empty()) throw std::runtime_error("--bf file has no usable SNP rows");
    double total = std::accumulate(scores.begin(), scores.end(), 0.0);
    std::vector<double> pip(scores.size(), 0.0);
    if (total > 0.0) {
        const double p_log_min = options_.p_min_set ? options_.p_min : -std::log10(std::max<size_t>(1, scores.size()));
        const double p_log_max = options_.p_max_set ? options_.p_max : 0.0;
        const double pi_mid = std::pow(10.0, 0.5 * (p_log_min + p_log_max));
        const double prior_gamma = options_.gamma_mean > 0.0 ? options_.gamma_mean : pi_mid * static_cast<double>(scores.size());
        const double expected_gamma = std::min<double>(
            std::max<double>(options_.s_min, std::min<double>(options_.top_sparse, prior_gamma)),
            scores.size());
        for (size_t i = 0; i < scores.size(); ++i) {
            pip[i] = std::min(1.0, scores[i] / total * expected_gamma);
        }
    }

    std::ofstream gamma(options_.out_prefix + ".gamma.txt");
    std::ofstream hyp(options_.out_prefix + ".hyp.txt");
    std::ofstream coef(options_.out_prefix + ".coef.txt");
    if (!gamma || !hyp || !coef) throw std::runtime_error("Cannot write DAP fine-mapping outputs");
    hyp << "h\trho\tsa2\tsb2\tweight\n";
    hyp << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < hyper.size(); ++i) {
        hyp << hyper[i][0] << '\t' << hyper[i][1] << '\t' << hyper[i][2] << '\t'
            << hyper[i][3] << '\t' << hyp_weight[i] << '\n';
    }
    gamma << "rs\tgamma\n" << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < rsids.size(); ++i) gamma << rsids[i] << '\t' << pip[i] << '\n';
    coef << "coef\n" << std::scientific << std::setprecision(6);
    for (double p : pip) coef << p << '\n';

    std::cout << "GEMMA-style BSLMM-DAP mode 5 complete. Outputs: "
              << options_.out_prefix << ".gamma.txt, .hyp.txt, .coef.txt\n"
              << "  Note: this is a deterministic posterior-weight summary over BF/HYP grids, not GEMMA's full DAP EM annotation model.\n";
    const double p_log_min = options_.p_min_set ? options_.p_min : -std::log10(std::max<size_t>(1, scores.size()));
    const double p_log_max = options_.p_max_set ? options_.p_max : 0.0;
    write_bslmm_meta(options_.out_prefix, options_.mode,
                     "deterministic-posterior-weight-dap-summary",
                     0, static_cast<int>(scores.size()), options_.h2, options_.h_min, options_.h_max,
                     options_.rho_min, options_.rho_max, p_log_min, p_log_max,
                     options_.gamma_mean, options_.s_min, options_.top_sparse,
                     options_.burnin, options_.sampling, options_.rpace, options_.wpace,
                     options_.seed);
}

} // namespace cosmic
