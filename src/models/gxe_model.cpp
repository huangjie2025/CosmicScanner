#include "models/gxe_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <set>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "eigen_vce.h"

namespace {

double gamma_q_gxe(double a, double x) {
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

double chisq_prob_gxe(double stat, int df) {
    if (stat <= 0.0) return 1.0;
    return gamma_q_gxe(df / 2.0, stat / 2.0);
}

} // anonymous namespace

namespace cosmic {

GxEModel::Options GxEModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.use_spa = cfg.use_spa;
    opts.rrm_order = cfg.rrm_order;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.tpheno = cfg.tpheno;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.env_covar_col = cfg.env_covar_col;
    return opts;
}

GxEModel::GxEModel(const Options& opts) : options_(opts) {}

void GxEModel::prepare(const std::vector<std::string>& fids,
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

    loadPhenotypeAndEnvironment();
    loadCovariates();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for GxE model.");
}

void GxEModel::loadPhenotypeAndEnvironment() {
    std::string file = options_.pheno_file;
    if (file.empty()) throw std::runtime_error("GxE model requires --pheno.");

    y_ = DataLoader::loadPhenotype(file, options_.mpheno, fids_, iids_, n_samples_, *plink_reader_, is_bgen_, is_pgen_);

    if (options_.env_covar_col > 0) {
        env_covar_ = Eigen::VectorXd::Zero(n_samples_);
        std::ifstream ifs(file);
        std::string line;
        int idx = 0;
        while (std::getline(ifs, line) && idx < n_samples_) {
            std::istringstream iss(line);
            std::vector<std::string> cols;
            std::string token;
            while (iss >> token) cols.push_back(token);
            if ((int)cols.size() > options_.env_covar_col - 1) {
                env_covar_(idx) = std::stod(cols[options_.env_covar_col - 1]);
            }
            idx++;
        }
    } else {
        env_covar_ = Eigen::VectorXd::Ones(n_samples_);
    }

    // Standardize environment covariate (mean 0, unit variance)
    // This is essential for the reaction norm model: K1 = G ⊙ (ee') would
    // otherwise be on a wildly different scale than K0 = G, causing ill-conditioned
    // V and REML divergence. Standardization is standard practice and the 2-df
    // joint test is invariant to linear rescaling of the environment.
    if (options_.env_covar_col > 0 && n_samples_ > 1) {
        double emean = env_covar_.head(n_samples_).mean();
        double evar = (env_covar_.head(n_samples_).array() - emean).square().mean();
        double esd = std::sqrt(evar);
        if (esd > 1e-10) {
            env_covar_.head(n_samples_) = (env_covar_.head(n_samples_).array() - emean) / esd;
        }
        std::cout << "  Environment covariate standardized (mean=0, sd=1) from raw mean="
                  << emean << " sd=" << esd << "\n";
    }
}

void GxEModel::loadCovariates() {
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
}

void GxEModel::runVCE() {
    std::cout << "Running GxE reaction norm VCE (AI-REML)...\n";

    int n = keep_.n_keep;
    Eigen::MatrixXd X_mat = X_.cast<double>();
    Eigen::VectorXd y_adj = y_.segment(0, n);
    Eigen::VectorXd env = env_covar_.segment(0, n);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = (float)options_.min_maf;
    GrmBuilder grm_builder(grm_opts);
    if (is_bgen_)
        grm_builder.compute(*bgen_reader_, keep_.keep_indices);
    else if (is_pgen_)
        grm_builder.compute(*pgen_reader_, keep_.keep_indices);
    else
        grm_builder.compute(*plink_reader_);
    Eigen::MatrixXd G = grm_builder.getGrm();

    // K0 = G, K1 = G ⊙ (ee')
    Eigen::MatrixXd eeT = env * env.transpose();
    Eigen::MatrixXd K0 = G;
    Eigen::MatrixXd K1(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            K1(i, j) = G(i, j) * eeT(i, j);

    sigma2_g0_ = 0.3;
    sigma2_g1_ = 0.1;
    sigma2_e_ = 0.6;

    // Initialize σ²_e from phenotypic variance for better starting point
    double vy = (y_adj.array() - y_adj.mean()).square().mean();
    sigma2_e_ = std::max(0.3, vy * 0.6);
    sigma2_g0_ = vy * 0.2;
    sigma2_g1_ = vy * 0.05;

    double prev_logL = -std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < options_.reml_maxit; ++iter) {
        // Save current values for convergence check
        double g0_prev = sigma2_g0_, g1_prev = sigma2_g1_, e_prev = sigma2_e_;

        Eigen::MatrixXd V = sigma2_g0_ * K0 + sigma2_g1_ * K1 + sigma2_e_ * I;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
        if (ldlt.info() != Eigen::Success) {
            // V not PD — add jitter and retry once
            V.diagonal().array() += 1e-6 * vy;
            ldlt.compute(V);
            if (ldlt.info() != Eigen::Success) break;
        }
        Eigen::MatrixXd Vinv = ldlt.solve(I);
        Eigen::MatrixXd VX = Vinv * X_mat;
        Eigen::MatrixXd XtVX = X_mat.transpose() * VX;
        Eigen::LDLT<Eigen::MatrixXd> xldlt(XtVX);
        if (xldlt.info() != Eigen::Success) break;
        Eigen::MatrixXd VX_VXt_inv = xldlt.solve(VX.transpose());
        Eigen::MatrixXd P = Vinv - VX * VX_VXt_inv;
        Eigen::VectorXd Py = P * y_adj;

        // V_i = ∂V/∂σ²_i: V_0 = K0, V_1 = K1, V_2 = I
        // AI_ij = 0.5 · y' · P · V_i · P · V_j · P · y
        // AI_ij = 0.5 · (PV_iPy)' · V_j · Py  ==  0.5 · (PK_iPy)' · (V_j · Py)
        Eigen::VectorXd PK0Py = P * (K0 * Py);
        Eigen::VectorXd PK1Py = P * (K1 * Py);
        Eigen::VectorXd PPy = P * Py;

        Eigen::Matrix3d AI;
        AI(0, 0) = 0.5 * (K0 * Py).dot(PK0Py);
        AI(1, 1) = 0.5 * (K1 * Py).dot(PK1Py);
        AI(2, 2) = 0.5 * Py.dot(PPy);
        AI(0, 1) = AI(1, 0) = 0.5 * (K0 * Py).dot(PK1Py);
        AI(0, 2) = AI(2, 0) = 0.5 * (K0 * Py).dot(PPy);
        AI(1, 2) = AI(2, 1) = 0.5 * (K1 * Py).dot(PPy);

        // Gradient: g_i = -0.5 · tr(P·V_i) + 0.5 · y'·P·V_i·P·y
        Eigen::Vector3d grad;
        grad(0) = -0.5 * (P * K0).trace() + 0.5 * (K0 * Py).dot(Py);
        grad(1) = -0.5 * (P * K1).trace() + 0.5 * (K1 * Py).dot(Py);
        grad(2) = -0.5 * P.trace() + 0.5 * Py.squaredNorm();

        // Solve: AI · Δ = grad
        Eigen::Vector3d delta = AI.ldlt().solve(grad);

        // Step-halving line search
        double old_logL = -0.5 * (ldlt.vectorD().array().log().sum()
                          + xldlt.vectorD().array().log().sum()
                          + y_adj.dot(Py));
        bool accepted = false;
        for (int step = 0; step < 12; ++step) {
            double g0_new = sigma2_g0_ + delta(0);
            double g1_new = sigma2_g1_ + delta(1);
            double e_new = sigma2_e_ + delta(2);
            if (g0_new > 1e-6 && g1_new > 1e-6 && e_new > 1e-6) {
                Eigen::MatrixXd V_new = g0_new * K0 + g1_new * K1 + e_new * I;
                Eigen::LDLT<Eigen::MatrixXd> new_ldlt(V_new);
                if (new_ldlt.info() == Eigen::Success) {
                    Eigen::MatrixXd new_Vinv = new_ldlt.solve(I);
                    Eigen::MatrixXd new_VX = new_Vinv * X_mat;
                    Eigen::MatrixXd new_XtVX = X_mat.transpose() * new_VX;
                    Eigen::LDLT<Eigen::MatrixXd> new_xldlt(new_XtVX);
                    if (new_xldlt.info() == Eigen::Success) {
                        Eigen::MatrixXd new_C = new_xldlt.solve(new_VX.transpose());
                        Eigen::VectorXd new_Py = (new_Vinv - new_VX * new_C) * y_adj;
                        double new_logL = -0.5 * (new_ldlt.vectorD().array().log().sum()
                                          + new_xldlt.vectorD().array().log().sum()
                                          + y_adj.dot(new_Py));
                        if (new_logL >= old_logL - 1e-8) {
                            sigma2_g0_ = g0_new; sigma2_g1_ = g1_new; sigma2_e_ = e_new;
                            accepted = true;
                            break;
                        }
                    }
                }
            }
            delta *= 0.5;
        }
        if (!accepted) {
            // Fall back to EM update with ratio-based step for stability
            double yPK0y = Py.dot(K0 * Py), yPK1y = Py.dot(K1 * Py), yPPy = Py.squaredNorm();
            double tr_PK0 = (P * K0).trace(), tr_PK1 = (P * K1).trace(), tr_P = P.trace();
            // EM ratio update: σ²_new = σ²_old · (y'PK_iPy / tr_PK_i)
            // This is more stable than direct division when traces are near zero
            double ratio_g0 = (tr_PK0 > 1e-8) ? (yPK0y / (tr_PK0 * sigma2_g0_ + 1e-10)) : 1.0;
            double ratio_g1 = (tr_PK1 > 1e-8) ? (yPK1y / (tr_PK1 * sigma2_g1_ + 1e-10)) : 1.0;
            double ratio_e  = (tr_P   > 1e-8) ? (yPPy  / (tr_P   * sigma2_e_ + 1e-10)) : 1.0;
            // Damp the EM step to prevent oscillation
            double damp = 0.5;
            sigma2_g0_ *= (1.0 - damp + damp * std::max(0.0, ratio_g0));
            sigma2_g1_ *= (1.0 - damp + damp * std::max(0.0, ratio_g1));
            sigma2_e_  *= (1.0 - damp + damp * std::max(0.0, ratio_e));
        }

        // Boundary clamping: use small floor proportional to vy
        double floor_vc = std::max(1e-6, 1e-4 * vy);
        if (sigma2_g0_ < floor_vc) sigma2_g0_ = floor_vc;
        if (sigma2_g1_ < floor_vc) sigma2_g1_ = floor_vc;
        if (sigma2_e_  < floor_vc) sigma2_e_  = floor_vc;

        if (iter % 5 == 0 || iter == options_.reml_maxit - 1) {
            std::cout << "GxE iter " << iter << ": Vg0=" << sigma2_g0_
                      << " Vg1=" << sigma2_g1_ << " Ve=" << sigma2_e_ << "\n";
        }

        // Convergence: check relative change in all 3 parameters
        double cur_logL = -0.5 * (ldlt.vectorD().array().log().sum()
                        + xldlt.vectorD().array().log().sum()
                        + y_adj.dot(Py));
        if (iter > 3) {
            double chg = std::abs(sigma2_g0_ - g0_prev) / (g0_prev + 1e-10)
                       + std::abs(sigma2_g1_ - g1_prev) / (g1_prev + 1e-10)
                       + std::abs(sigma2_e_  - e_prev)  / (e_prev  + 1e-10);
            if (chg < 1e-7) break;
            // Also break if logL hasn't improved significantly
            if (std::abs(cur_logL - prev_logL) < 1e-6 && iter > 5) break;
        }
        prev_logL = cur_logL;
    }

    std::cout << "Final GxE VCE: Vg0=" << sigma2_g0_
              << " Vg1=" << sigma2_g1_ << " Ve=" << sigma2_e_ << "\n";
}

void GxEModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_gxe_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) return;
    ofs << "sigma2_g0\tsigma2_g1\tsigma2_e\n";
    ofs << sigma2_g0_ << "\t" << sigma2_g1_ << "\t" << sigma2_e_ << "\n";
    ofs.close();
    std::cout << "GxE null model written to " << fname << "\n";
}

void GxEModel::run() {
    runVCE();
    scanGwas();
    writeNullSummary();
}

// Helper: 2-df Wald test for GxE scan
struct GxeSnpTestResult {
    double beta_g = 0.0;
    double se_g = 0.0;
    double beta_gxe = 0.0;
    double se_gxe = 0.0;
    double joint_chisq = 0.0;
    double joint_p = 1.0;
};

static GxeSnpTestResult gxe_wald_2df(
    double score_g, double score_gxe,
    double m_gg, double m_g_gxe, double m_gxe_gxe)
{
    GxeSnpTestResult res;
    Eigen::Matrix2d M;
    M << m_gg, m_g_gxe, m_g_gxe, m_gxe_gxe;
    Eigen::LDLT<Eigen::Matrix2d> ldlt(M);
    if (ldlt.info() != Eigen::Success) {
        M.diagonal().array() += 1e-8;
        ldlt.compute(M);
        if (ldlt.info() != Eigen::Success) return res;
    }
    Eigen::Vector2d score(score_g, score_gxe);
    Eigen::Vector2d beta = ldlt.solve(score);
    if (!beta.allFinite()) return res;
    Eigen::Matrix2d M_inv = ldlt.solve(Eigen::Matrix2d::Identity());
    res.beta_g = beta(0);
    res.se_g = (M_inv(0,0) > 0) ? std::sqrt(M_inv(0,0)) : std::numeric_limits<double>::quiet_NaN();
    res.beta_gxe = beta(1);
    res.se_gxe = (M_inv(1,1) > 0) ? std::sqrt(M_inv(1,1)) : std::numeric_limits<double>::quiet_NaN();
    res.joint_chisq = score.dot(beta);
    res.joint_p = chisq_prob_gxe(res.joint_chisq, 2);
    return res;
}

// (scan_gxe_variants template removed — inline loops used instead)

void GxEModel::scanGwas() {
    std::cout << "Running GxE 2-df joint scan...\n";

    int n = keep_.n_keep;
    Eigen::VectorXd y_sub = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X_, keep_);
    Eigen::VectorXd env = env_covar_.segment(0, n);
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd Xd = X_sub.cast<double>();

    // Environment outer product for K1
    Eigen::MatrixXd eeT = env * env.transpose();

    // Build global GRM
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.algorithm = options_.grm_alg;
    GrmBuilder grm_builder(grm_opts);
    if (is_bgen_)
        grm_builder.compute(*bgen_reader_, keep_.keep_indices);
    else if (is_pgen_)
        grm_builder.compute(*pgen_reader_, keep_.keep_indices);
    else
        grm_builder.compute(*plink_reader_);
    Eigen::MatrixXd G = grm_builder.getGrm();

    // Determine chromosomes
    std::set<std::string> chrom_set;
    if (is_bgen_) {
        bgen_reader_->loadVariantInfo();
        for (const auto& v : bgen_reader_->getVariants()) chrom_set.insert(v.chrom);
    } else if (is_pgen_) {
        for (const auto& v : pgen_reader_->getVariants()) chrom_set.insert(v.chrom);
    } else {
        for (const auto& s : plink_reader_->getBimInfo()) chrom_set.insert(s.chrom);
    }

    bool use_loco = chrom_set.size() > 1 && !options_.no_loco;

    std::vector<std::string> chrom_labels(chrom_set.begin(), chrom_set.end());
    std::sort(chrom_labels.begin(), chrom_labels.end(), [](const std::string& a, const std::string& b) {
        try { return std::stoi(a) < std::stoi(b); }
        catch (...) { return a < b; }
    });

    // Build split GRM for LOCO
    GrmBuilder::SplitGrmResult split;
    double S_global = 1.0;
    Eigen::MatrixXd G_sub;
    if (use_loco) {
        GrmBuilder split_builder(grm_opts);
        if (is_bgen_) {
            bgen_reader_->reset();
            split = split_builder.computeSplit(*bgen_reader_, chrom_set);
        } else if (is_pgen_) {
            pgen_reader_->reset();
            split = split_builder.computeSplit(*pgen_reader_, chrom_set);
        } else {
            split = split_builder.computeSplit(*plink_reader_, chrom_set);
        }
        S_global = (options_.grm_alg == 0) ? (double)split.global_snp_count : split.global_sum_2pq;
        G_sub = split.global;
    } else {
        G_sub = G;
    }

    // Build per-chrom raw K1 matrix: K1_chr = G_chr ⊙ ee' (for LOCO reconstruction)
    // Pre-compute eeT to avoid per-chrom recomputation

    // Output
    std::string out_file = options_.out_prefix + ".gxe.assoc.txt";
    std::ofstream ofs(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA_G\tSE_G\tBETA_GxE\tSE_GxE\tJOINT_CHISQ\tJOINT_P\n";
    ofs.close();

    // Build extraction/exclusion sets
    std::unordered_set<std::string> extract_ids, exclude_ids;
    if (!options_.extract_file.empty()) {
        std::ifstream ef(options_.extract_file);
        std::string id; while (ef >> id) extract_ids.insert(id);
    }
    if (!options_.exclude_file.empty()) {
        std::ifstream xf(options_.exclude_file);
        std::string id; while (xf >> id) exclude_ids.insert(id);
    }

    // Lambda to build V, factor, compute Py, then scan
    auto scan_with_grm = [&](const Eigen::MatrixXd& G_work,
                             const std::string& chrom_filter) {
        // Build K0, K1 from G_work
        Eigen::MatrixXd K0 = G_work;
        Eigen::MatrixXd K1(G_work.rows(), G_work.cols());
        for (int i = 0; i < G_work.rows(); ++i)
            for (int j = 0; j < G_work.cols(); ++j)
                K1(i, j) = G_work(i, j) * eeT(i, j);

        // V = sigma2_g0 * K0 + sigma2_g1 * K1 + sigma2_e * I
        Eigen::MatrixXd V = sigma2_g0_ * K0 + sigma2_g1_ * K1 + sigma2_e_ * I;

        // Factor V and compute P = V⁻¹ - V⁻¹·X·(X'·V⁻¹·X)⁻¹·X'·V⁻¹
        Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
        if (ldlt.info() != Eigen::Success) {
            V.diagonal().array() += 1e-6;
            ldlt.compute(V);
            if (ldlt.info() != Eigen::Success) {
                std::cerr << "  V factorization failed for GxE scan, skipping.\n";
                return;
            }
        }
        Eigen::MatrixXd VinvX = ldlt.solve(Xd);
        Eigen::MatrixXd XtVinvX = Xd.transpose() * VinvX;
        Eigen::LDLT<Eigen::MatrixXd> xldlt(XtVinvX);
        if (xldlt.info() != Eigen::Success) {
            XtVinvX.diagonal().array() += 1e-6;
            xldlt.compute(XtVinvX);
        }
        Eigen::MatrixXd C = xldlt.solve(VinvX.transpose());
        Eigen::VectorXd Vinvy = ldlt.solve(y_sub);
        Eigen::VectorXd Py = Vinvy - VinvX * xldlt.solve(Xd.transpose() * Vinvy);

        // Build LDLT solver of V for per-SNP solves
        Eigen::MatrixXd V_factor = V;  // copy for LDLT
        Eigen::LDLT<Eigen::MatrixXd> ldlt_V(V_factor);
        if (ldlt_V.info() != Eigen::Success) {
            V_factor.diagonal().array() += 1e-6;
            ldlt_V.compute(V_factor);
        }

        // Pre-compute Python residual: Py = P·y where P = V⁻¹ - V⁻¹X·C·X'V⁻¹
        // We already have Py from above

        // Scan variants
        std::vector<GwasResult> results;

        if (is_bgen_) {
            // Reset and sequential read for BGEN
            bgen_reader_->reset();
            size_t n_var = bgen_reader_->getNumVariants();
            std::string chrom, rsid, ref, alt;
            uint32_t pos = 0;
            Eigen::VectorXf dosage_f;
            results.clear();
            for (size_t j = 0; j < n_var; ++j) {
                if (!bgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                if (!extract_ids.empty() && extract_ids.find(rsid) == extract_ids.end()) continue;
                if (!exclude_ids.empty() && exclude_ids.find(rsid) != exclude_ids.end()) continue;

                Eigen::VectorXd g(n);
                double sum_g = 0.0;
                int valid_n = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep_.keep_indices[i]];
                    g[i] = val;
                    if (!std::isnan(val)) { sum_g += val; valid_n++; }
                }
                if (valid_n < 3) continue;
                double freq = sum_g / (2.0 * valid_n);
                double maf = std::min(freq, 1.0 - freq);
                if (maf < (float)options_.min_maf) continue;
                double mean = 2.0 * freq;
                for (int i = 0; i < n; ++i) {
                    if (std::isnan(g[i])) g[i] = 0.0;
                    else g[i] -= mean;
                }

                Eigen::VectorXd gxe(n);
                double sum_gxe = 0.0;
                for (int i = 0; i < n; ++i) { gxe[i] = g[i] * env[i]; sum_gxe += gxe[i]; }
                double gxe_mean = sum_gxe / n;
                for (int i = 0; i < n; ++i) gxe[i] -= gxe_mean;

                Eigen::VectorXd Pg = ldlt_V.solve(g);
                Eigen::VectorXd Pgxe = ldlt_V.solve(gxe);
                if (!Pg.allFinite() || !Pgxe.allFinite()) continue;

                double score_g = g.dot(Py), score_gxe = gxe.dot(Py);
                double m_gg = g.dot(Pg), m_g_gxe = g.dot(Pgxe), m_gxe_gxe = gxe.dot(Pgxe);
                GxeSnpTestResult r = gxe_wald_2df(score_g, score_gxe, m_gg, m_g_gxe, m_gxe_gxe);
                if (!std::isfinite(r.joint_p)) continue;

                GwasResult res;
                res.chrom = chrom; res.snp_id = rsid; res.pos = pos;
                res.a1 = alt; res.a2 = ref; res.freq = freq; res.n_obs = valid_n;
                res.beta = r.beta_g; res.se = r.se_g;
                res.beta_gxe = r.beta_gxe; res.se_gxe = r.se_gxe;
                res.joint_chisq = r.joint_chisq;
                res.p_value = r.joint_p;
                res.log10_p = (res.p_value > 0) ? -std::log10(res.p_value) : 300.0;
                results.push_back(std::move(res));
            }
        } else if (is_pgen_) {
            pgen_reader_->reset();
            size_t n_var = pgen_reader_->getNumVariants();
            std::string chrom, rsid, ref, alt;
            uint32_t pos = 0;
            Eigen::VectorXf dosage_f;
            results.clear();
            for (size_t j = 0; j < n_var; ++j) {
                if (!pgen_reader_->readVariant(chrom, rsid, pos, ref, alt, dosage_f)) break;
                if (!chrom_filter.empty() && chrom != chrom_filter) continue;
                if (!extract_ids.empty() && extract_ids.find(rsid) == extract_ids.end()) continue;
                if (!exclude_ids.empty() && exclude_ids.find(rsid) != exclude_ids.end()) continue;

                Eigen::VectorXd g(n);
                double sum_g = 0.0; int valid_n = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep_.keep_indices[i]];
                    g[i] = val;
                    if (!std::isnan(val)) { sum_g += val; valid_n++; }
                }
                if (valid_n < 3) continue;
                double freq = sum_g / (2.0 * valid_n);
                double maf = std::min(freq, 1.0 - freq);
                if (maf < (float)options_.min_maf) continue;
                double mean = 2.0 * freq;
                for (int i = 0; i < n; ++i) {
                    if (std::isnan(g[i])) g[i] = 0.0; else g[i] -= mean;
                }

                Eigen::VectorXd gxe(n);
                double sum_gxe = 0.0;
                for (int i = 0; i < n; ++i) { gxe[i] = g[i] * env[i]; sum_gxe += gxe[i]; }
                double gxe_mean = sum_gxe / n;
                for (int i = 0; i < n; ++i) gxe[i] -= gxe_mean;

                Eigen::VectorXd Pg = ldlt_V.solve(g);
                Eigen::VectorXd Pgxe = ldlt_V.solve(gxe);
                if (!Pg.allFinite() || !Pgxe.allFinite()) continue;

                double score_g = g.dot(Py), score_gxe = gxe.dot(Py);
                double m_gg = g.dot(Pg), m_g_gxe = g.dot(Pgxe), m_gxe_gxe = gxe.dot(Pgxe);
                GxeSnpTestResult r = gxe_wald_2df(score_g, score_gxe, m_gg, m_g_gxe, m_gxe_gxe);
                if (!std::isfinite(r.joint_p)) continue;

                GwasResult res;
                res.chrom = chrom; res.snp_id = rsid; res.pos = pos;
                res.a1 = alt; res.a2 = ref; res.freq = freq; res.n_obs = valid_n;
                res.beta = r.beta_g; res.se = r.se_g;
                res.beta_gxe = r.beta_gxe; res.se_gxe = r.se_gxe;
                res.joint_chisq = r.joint_chisq;
                res.p_value = r.joint_p;
                res.log10_p = (res.p_value > 0) ? -std::log10(res.p_value) : 300.0;
                results.push_back(std::move(res));
            }
        } else {
            // PLINK: index-based access
            size_t m = plink_reader_->getNumSnps();
            const auto& bim = plink_reader_->getBimInfo();
            results.clear();
            results.reserve(m);
            for (size_t j = 0; j < m; ++j) {
                if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;
                if (!extract_ids.empty() && extract_ids.find(bim[j].id) == extract_ids.end()) continue;
                if (!exclude_ids.empty() && exclude_ids.find(bim[j].id) != exclude_ids.end()) continue;

                Eigen::VectorXf dosage_f;
                plink_reader_->getSnpGenotypes(j, dosage_f);

                Eigen::VectorXd g(n);
                double sum_g = 0.0; int valid_n = 0;
                for (int i = 0; i < n; ++i) {
                    float val = dosage_f[keep_.keep_indices[i]];
                    g[i] = val;
                    if (!std::isnan(val)) { sum_g += val; valid_n++; }
                }
                if (valid_n < 3) continue;
                double freq = sum_g / (2.0 * valid_n);
                double maf = std::min(freq, 1.0 - freq);
                if (maf < (float)options_.min_maf) continue;
                double mean = 2.0 * freq;
                for (int i = 0; i < n; ++i) {
                    if (std::isnan(g[i])) g[i] = 0.0; else g[i] -= mean;
                }

                Eigen::VectorXd gxe(n);
                double sum_gxe = 0.0;
                for (int i = 0; i < n; ++i) { gxe[i] = g[i] * env[i]; sum_gxe += gxe[i]; }
                double gxe_mean = sum_gxe / n;
                for (int i = 0; i < n; ++i) gxe[i] -= gxe_mean;

                Eigen::VectorXd Pg = ldlt_V.solve(g);
                Eigen::VectorXd Pgxe = ldlt_V.solve(gxe);
                if (!Pg.allFinite() || !Pgxe.allFinite()) continue;

                double score_g = g.dot(Py), score_gxe = gxe.dot(Py);
                double m_gg = g.dot(Pg), m_g_gxe = g.dot(Pgxe), m_gxe_gxe = gxe.dot(Pgxe);
                GxeSnpTestResult r = gxe_wald_2df(score_g, score_gxe, m_gg, m_g_gxe, m_gxe_gxe);
                if (!std::isfinite(r.joint_p)) continue;

                GwasResult res;
                res.chrom = bim[j].chrom; res.snp_id = bim[j].id; res.pos = bim[j].bp_pos;
                res.a1 = bim[j].alt; res.a2 = bim[j].ref; res.freq = freq; res.n_obs = valid_n;
                res.beta = r.beta_g; res.se = r.se_g;
                res.beta_gxe = r.beta_gxe; res.se_gxe = r.se_gxe;
                res.joint_chisq = r.joint_chisq;
                res.p_value = r.joint_p;
                res.log10_p = (res.p_value > 0) ? -std::log10(res.p_value) : 300.0;
                results.push_back(std::move(res));
            }
        }

        // Append results
        std::ofstream fout(out_file, std::ios::app);
        for (const auto& r : results) {
            if (std::isnan(r.p_value)) continue;
            fout << r.chrom << "\t" << r.snp_id << "\t" << r.pos << "\t"
                 << r.a1 << "\t" << r.a2 << "\t" << r.freq << "\t"
                 << r.n_obs << "\t" << r.beta << "\t" << r.se << "\t"
                 << r.beta_gxe << "\t" << r.se_gxe << "\t"
                 << r.joint_chisq << "\t" << r.p_value << "\n";
        }
    };

    if (!use_loco) {
        std::cout << "Building global V and scanning (no LOCO)...\n";
        scan_with_grm(G_sub, "");
    } else {
        std::cout << "Running GxE LOCO scan on " << split.chrom_labels.size() << " chromosomes...\n";
        for (size_t ci = 0; ci < split.chrom_labels.size(); ++ci) {
            const std::string& chrom = split.chrom_labels[ci];
            std::cout << "  LOCO chr " << chrom << "...\n";

            double S_chr = (options_.grm_alg == 0)
                ? (double)split.per_chrom_snp_count[ci]
                : split.per_chrom_sum_2pq[ci];
            double S_loco = S_global - S_chr;

            if (S_chr <= 0.0 || S_loco <= 1e-6) {
                scan_with_grm(G_sub, chrom);
                continue;
            }

            // Build G_loco
            const auto& G_chr_raw = split.per_chrom_raw[ci];
            Eigen::MatrixXd G_loco(n, n);
#pragma omp parallel for collapse(2)
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    G_loco(i, j) = G_sub(i, j) * (S_global / S_loco)
                                 - G_chr_raw(i, j) / S_loco;
                }
            }
            scan_with_grm(G_loco, chrom);
        }
    }
    std::cout << "GxE results saved to " << out_file << "\n";
}

} // namespace cosmic
