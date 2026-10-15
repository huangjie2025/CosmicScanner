#include "models/gemma_bslmm_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
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
    opts.thin = cfg.bslmm_thin;
    opts.seed = cfg.bslmm_seed;
    opts.mcmc = cfg.bslmm_mcmc;
    opts.bvsr = cfg.bvsr;
    return opts;
}

GemmaBslmmModel::GemmaBslmmModel(const Options& options) : options_(options) {}

void GemmaBslmmModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              VcfReader& vcf_reader,
                              bool is_bgen, bool is_pgen, bool is_vcf) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    vcf_reader_ = &vcf_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;
    is_vcf_ = is_vcf;

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
    if ((options_.mcmc || options_.bvsr) &&
        (options_.burnin < 0 || options_.sampling < 1 || options_.thin < 1)) {
        throw std::runtime_error("BSLMM/BVSR MCMC requires non-negative burn-in and positive samples/thinning");
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

    if (options_.mcmc || options_.bvsr) {
        runPosteriorMcmc(Z, variants, y_work);
        return;
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

void GemmaBslmmModel::runPosteriorMcmc(const Eigen::MatrixXd& Z,
                                        const std::vector<VariantMeta>& variants,
                                        const Eigen::VectorXd& y_work) const {
    const int n = Z.rows();
    const int m = Z.cols();
    const bool use_polygenic = !options_.bvsr;
    const std::uint64_t seed = options_.seed == 0 ? 5489u : static_cast<std::uint64_t>(options_.seed);
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    const double y_variance = std::max(1e-8, y_work.squaredNorm() / std::max(1, n - 1));
    const double pi_min = options_.p_min_set ? std::pow(10.0, options_.p_min)
                                             : 1.0 / static_cast<double>(m);
    const double pi_max = options_.p_max_set ? std::pow(10.0, options_.p_max) : 1.0;
    double pi = std::min(pi_max, std::max(pi_min,
        options_.gamma_mean > 0.0 ? options_.gamma_mean / static_cast<double>(m)
                                  : std::min(0.1, 10.0 / static_cast<double>(m))));
    double residual_variance = 0.5 * y_variance;
    double polygenic_variance = use_polygenic ? 0.25 * y_variance : 0.0;
    double slab_variance = y_variance / std::max(1.0, pi * static_cast<double>(m));

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd residual = y_work;
    std::vector<unsigned char> included(m, 0);
    std::vector<int> marker_order(m);
    std::iota(marker_order.begin(), marker_order.end(), 0);

    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    std::vector<int> positive_eigenvalues;
    if (use_polygenic) {
        Eigen::MatrixXd kinship = (Z * Z.transpose()) / static_cast<double>(m);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(kinship);
        if (eigensolver.info() != Eigen::Success)
            throw std::runtime_error("BSLMM MCMC polygenic GRM eigendecomposition failed");
        eigenvalues = eigensolver.eigenvalues().cwiseMax(0.0);
        eigenvectors = eigensolver.eigenvectors();
        for (int i = 0; i < n; ++i)
            if (eigenvalues[i] > 1e-10) positive_eigenvalues.push_back(i);
    }

    Eigen::VectorXd beta_sum = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd pip_sum = Eigen::VectorXd::Zero(m);
    double pve_sum = 0.0;
    double pge_sum = 0.0;
    double pi_sum = 0.0;
    double residual_variance_sum = 0.0;
    double polygenic_variance_sum = 0.0;
    double pve_first_half = 0.0;
    double pve_second_half = 0.0;
    int first_half_count = 0;
    int second_half_count = 0;
    long long inclusion_changes = 0;
    long long inclusion_updates = 0;
    int retained = 0;

    auto inverse_gamma = [&](double shape, double scale) {
        std::gamma_distribution<double> gamma(shape, 1.0);
        return scale / std::max(1e-12, gamma(rng));
    };
    auto sample_beta = [&](double a, double b) {
        std::gamma_distribution<double> ga(a, 1.0), gb(b, 1.0);
        const double x = ga(rng);
        const double y = gb(rng);
        return x / std::max(1e-12, x + y);
    };

    const int total_iterations = options_.burnin + options_.sampling;
    for (int iteration = 0; iteration < total_iterations; ++iteration) {
        if (use_polygenic) {
            residual += u;
            Eigen::VectorXd transformed = eigenvectors.transpose() * residual;
            Eigen::VectorXd sampled = Eigen::VectorXd::Zero(n);
            for (int index : positive_eigenvalues) {
                const double prior_variance = std::max(1e-12, polygenic_variance * eigenvalues[index]);
                const double posterior_variance = 1.0 /
                    (1.0 / prior_variance + 1.0 / residual_variance);
                const double posterior_mean = posterior_variance * transformed[index] / residual_variance;
                sampled[index] = posterior_mean + std::sqrt(posterior_variance) * normal(rng);
            }
            u = eigenvectors * sampled;
            residual -= u;
        }

        std::shuffle(marker_order.begin(), marker_order.end(), rng);
        int included_count = 0;
        double beta_square_sum = 0.0;
        for (int marker : marker_order) {
            residual.noalias() += Z.col(marker) * beta[marker];
            const double precision = Z.col(marker).squaredNorm() / residual_variance + 1.0 / slab_variance;
            const double posterior_variance = 1.0 / std::max(1e-12, precision);
            const double posterior_mean = posterior_variance * Z.col(marker).dot(residual) / residual_variance;
            const double log_bayes_factor = 0.5 *
                (std::log(posterior_variance / slab_variance) +
                 posterior_mean * posterior_mean / posterior_variance);
            const double log_prior_odds = std::log(std::max(1e-12, pi)) -
                                          std::log(std::max(1e-12, 1.0 - pi));
            const double log_odds = std::max(-40.0, std::min(40.0, log_prior_odds + log_bayes_factor));
            const double inclusion_probability = 1.0 / (1.0 + std::exp(-log_odds));
            const bool old_state = included[marker] != 0;
            included[marker] = uniform(rng) < inclusion_probability;
            beta[marker] = included[marker]
                ? posterior_mean + std::sqrt(posterior_variance) * normal(rng)
                : 0.0;
            residual.noalias() -= Z.col(marker) * beta[marker];
            inclusion_changes += old_state != (included[marker] != 0);
            ++inclusion_updates;
            if (included[marker]) {
                ++included_count;
                beta_square_sum += beta[marker] * beta[marker];
            }
        }

        slab_variance = inverse_gamma(2.0 + 0.5 * included_count,
                                      0.01 * y_variance + 0.5 * beta_square_sum);
        pi = std::min(pi_max, std::max(pi_min,
            sample_beta(1.0 + included_count, 1.0 + m - included_count)));
        residual_variance = inverse_gamma(2.0 + 0.5 * n,
                                          0.1 * y_variance + 0.5 * residual.squaredNorm());
        if (use_polygenic && !positive_eigenvalues.empty()) {
            const Eigen::VectorXd transformed_u = eigenvectors.transpose() * u;
            double quadratic = 0.0;
            for (int index : positive_eigenvalues)
                quadratic += transformed_u[index] * transformed_u[index] / eigenvalues[index];
            polygenic_variance = inverse_gamma(
                2.0 + 0.5 * positive_eigenvalues.size(),
                0.1 * y_variance + 0.5 * quadratic);
        }

        if (iteration < options_.burnin ||
            ((iteration - options_.burnin) % options_.thin) != 0) continue;
        ++retained;
        beta_sum += beta;
        for (int marker = 0; marker < m; ++marker) pip_sum[marker] += included[marker] != 0;
        const Eigen::VectorXd sparse = Z * beta;
        const Eigen::VectorXd genetic = u + sparse;
        const double genetic_variance = genetic.squaredNorm() / std::max(1, n - 1);
        const double sparse_variance = sparse.squaredNorm() / std::max(1, n - 1);
        const double pve = genetic_variance /
                           std::max(1e-12, genetic_variance + residual_variance);
        const double pge = genetic_variance > 1e-12 ? sparse_variance / genetic_variance : 0.0;
        pve_sum += pve;
        pge_sum += std::min(1.0, std::max(0.0, pge));
        pi_sum += pi;
        residual_variance_sum += residual_variance;
        polygenic_variance_sum += polygenic_variance;
        if (retained <= (options_.sampling / options_.thin + 1) / 2) {
            pve_first_half += pve;
            ++first_half_count;
        } else {
            pve_second_half += pve;
            ++second_half_count;
        }
    }
    if (retained == 0) throw std::runtime_error("BSLMM/BVSR MCMC retained no posterior samples");

    const Eigen::VectorXd posterior_beta = beta_sum / static_cast<double>(retained);
    const Eigen::VectorXd pip = pip_sum / static_cast<double>(retained);
    const std::string effects_path = options_.out_prefix +
        (options_.bvsr ? ".bvsr.effects.tsv" : ".bslmm.effects.tsv");
    std::ofstream effects(effects_path);
    if (!effects) throw std::runtime_error("Cannot write BSLMM/BVSR posterior effects");
    effects << "SNP\tA1\tA2\tPOSTERIOR_BETA\tPIP\n" << std::setprecision(12);
    for (int marker = 0; marker < m; ++marker) {
        effects << variants[marker].id << '\t' << variants[marker].a1 << '\t'
                << variants[marker].a2 << '\t'
                << posterior_beta[marker] / variants[marker].sd << '\t'
                << pip[marker] << '\n';
    }

    const double first_half_mean = first_half_count > 0 ? pve_first_half / first_half_count : 0.0;
    const double second_half_mean = second_half_count > 0 ? pve_second_half / second_half_count : first_half_mean;
    std::ofstream posterior(options_.out_prefix + ".bslmm.posterior.tsv");
    if (!posterior) throw std::runtime_error("Cannot write BSLMM/BVSR posterior summary");
    posterior << "METRIC\tVALUE\n" << std::setprecision(12);
    posterior << "MODEL\t" << (options_.bvsr ? "BVSR" : "BSLMM") << '\n';
    posterior << "PVE\t" << pve_sum / retained << '\n';
    posterior << "PGE\t" << pge_sum / retained << '\n';
    posterior << "PI\t" << pi_sum / retained << '\n';
    posterior << "RESIDUAL_VARIANCE\t" << residual_variance_sum / retained << '\n';
    posterior << "POLYGENIC_VARIANCE\t" << polygenic_variance_sum / retained << '\n';
    posterior << "ITERATIONS\t" << total_iterations << '\n';
    posterior << "BURNIN\t" << options_.burnin << '\n';
    posterior << "SAMPLING_ITERATIONS\t" << options_.sampling << '\n';
    posterior << "THIN\t" << options_.thin << '\n';
    posterior << "SAMPLES_RETAINED\t" << retained << '\n';
    posterior << "INCLUSION_CHANGE_RATE\t"
              << static_cast<double>(inclusion_changes) / std::max<long long>(1, inclusion_updates) << '\n';
    posterior << "PVE_SPLIT_HALF_DIFFERENCE\t" << std::abs(first_half_mean - second_half_mean) << '\n';
    posterior << "SEED\t" << seed << '\n';

    std::cout << (options_.bvsr ? "BVSR" : "BSLMM")
              << " MCMC complete: retained " << retained << " posterior samples. Outputs: "
              << effects_path << " and " << options_.out_prefix << ".bslmm.posterior.tsv\n";
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
