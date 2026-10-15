#include "models/maternal_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <set>
#include <algorithm>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "grm_builder.h"
#include "pheno_reader.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "eigen_vce.h"

namespace cosmic {

MaternalModel::Options MaternalModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.maternal_id_col = cfg.maternal_id_col;
    opts.maternal_snp_test = cfg.maternal_snp_test;
    return opts;
}

MaternalModel::MaternalModel(const Options& opts) : options_(opts) {}

void MaternalModel::prepare(const std::vector<std::string>& fids,
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

    for (int i = 0; i < (int)iids_.size(); ++i) {
        std::string key = fids_[i] + ":" + iids_[i];
        iid_to_subject_idx_[key] = i;
    }

    loadPhenotype();
    loadCovariates();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " subjects remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid subjects for maternal model.");
}

void MaternalModel::loadPhenotype() {
    std::string file = options_.pheno_file;
    if (file.empty()) throw std::runtime_error("Maternal model requires --pheno file with columns: FID IID DAM_ID Pheno");

    std::ifstream ifs(file);
    if (!ifs.is_open()) throw std::runtime_error("Cannot open phenotype file: " + file);

    std::string line;
    int n_animals = n_samples_;
    std::vector<double> y_vals;
    std::vector<int> animal_idx, dam_idx;
    std::unordered_map<std::string, int> dam_map;

    y_vals.reserve(n_animals * 2);
    animal_idx.reserve(n_animals * 2);
    dam_idx.reserve(n_animals * 2);

    int dam_col = options_.maternal_id_col - 1;
    int pheno_col = 2 + options_.mpheno;
    int min_cols = std::max(options_.maternal_id_col, options_.mpheno + 3);

    int line_no = 0;
    bool header = true;
    while (std::getline(ifs, line)) {
        line_no++;
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> cols;
        while (iss >> token) cols.push_back(token);

        if ((int)cols.size() < min_cols) {
            if (header && line_no == 1) { header = false; continue; }
            std::cerr << "Warning: skipping malformed line " << line_no << "\n";
            continue;
        }
        // Header detection: if first line's pheno column is non-numeric, skip as header.
        // This must happen BEFORE std::stod to avoid "Error: stod" on header rows like
        // "FID IID DAM_ID Pheno".
        if (header && line_no == 1) {
            int pc = pheno_col < (int)cols.size() ? pheno_col : (int)cols.size() - 1;
            bool is_numeric = true;
            try { std::stod(cols[pc]); } catch (...) { is_numeric = false; }
            if (!is_numeric) { header = false; continue; }
        }
        header = false;

        std::string fid = cols[0], iid = cols[1], dam_id = cols[dam_col];
        if (pheno_col >= (int)cols.size()) pheno_col = (int)cols.size() - 1;
        double pheno_val = std::stod(cols[pheno_col]);

        std::string animal_key = fid + ":" + iid;
        std::string dam_key = fid + ":" + dam_id;

        auto it = iid_to_subject_idx_.find(animal_key);
        if (it == iid_to_subject_idx_.end()) continue;

        int ai = it->second;
        y_vals.push_back(pheno_val);
        animal_idx.push_back(ai);

        auto dit = iid_to_subject_idx_.find(dam_key);
        if (dit != iid_to_subject_idx_.end()) {
            dam_idx.push_back(dit->second);
        } else {
            auto dmit = dam_map.find(dam_key);
            if (dmit == dam_map.end()) {
                int new_idx = n_animals + (int)dam_map.size();
                dam_map[dam_key] = new_idx;
                dam_idx.push_back(new_idx);
            } else {
                dam_idx.push_back(dmit->second);
            }
        }
    }
    ifs.close();

    n_animals_ = n_animals;
    n_dams_ = n_animals + (int)dam_map.size();
    record_to_animal_ = std::move(animal_idx);
    record_to_dam_ = std::move(dam_idx);

    int n_records = (int)y_vals.size();
    y_.resize(n_records);
    for (int i = 0; i < n_records; ++i) y_(i) = y_vals[i];

    std::cout << "Loaded " << n_records << " records for " << n_animals
              << " animals and " << dam_map.size() << " unique dams.\n";
    if (n_records < 10) throw std::runtime_error("Too few records for maternal model.");
}

void MaternalModel::loadCovariates() {
    Eigen::MatrixXf X0 = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
    int n_records = (int)y_.size();
    int n_covars = X0.cols();
    X_.resize(n_records, n_covars);
    for (int r = 0; r < n_records; ++r) {
        int ai = record_to_animal_[r];
        for (int c = 0; c < X0.cols(); ++c) {
            X_(r, c) = ai < X0.rows() ? X0(ai, c) : 0.0f;
        }
    }
}

void MaternalModel::loadDamMap() {
    // Dam mapping already done in loadPhenotype
}

void MaternalModel::buildGRMs() {
    std::cout << "Building direct genomic GRM...\n";
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
    G_direct_ = grm_builder.getGrm();
    n_animals_ = G_direct_.rows();

    std::cout << "Maternal GRM uses direct GRM (full implementation pending).\n";
    G_maternal_ = G_direct_;
    n_dams_ = G_maternal_.rows();

    if (options_.direct_maternal_cov) {
        G_dm_ = G_direct_;
    }

    std::cout << "Direct GRM: " << G_direct_.rows() << "x" << G_direct_.cols() << "\n";
    std::cout << "Maternal GRM: " << G_maternal_.rows() << "x" << G_maternal_.cols() << "\n";
}

void MaternalModel::runVCE() {
    std::cout << "Running multi-component VCE for maternal model...\n";

    int n_records = (int)y_.size();
    Eigen::MatrixXd K(n_records, n_records);
    K.setZero();

    for (int i = 0; i < n_records; ++i) {
        for (int j = 0; j < n_records; ++j) {
            double v = G_direct_(record_to_animal_[i], record_to_animal_[j]);
            K(i, j) = v;
        }
    }

    Eigen::MatrixXd K_m(n_records, n_records);
    K_m.setZero();
    for (int i = 0; i < n_records; ++i) {
        for (int j = 0; j < n_records; ++j) {
            double v = G_maternal_(
                record_to_dam_[i] < G_maternal_.rows() ? record_to_dam_[i] : 0,
                record_to_dam_[j] < G_maternal_.rows() ? record_to_dam_[j] : 0);
            K_m(i, j) = v;
        }
    }

    Eigen::MatrixXd K_dm(n_records, n_records);
    K_dm.setZero();
    if (options_.direct_maternal_cov) {
        for (int i = 0; i < n_records; ++i) {
            for (int j = 0; j < n_records; ++j) {
                double v = G_dm_(
                    record_to_animal_[i] < G_dm_.rows() ? record_to_animal_[i] : 0,
                    record_to_dam_[j] < G_dm_.rows() ? record_to_dam_[j] : 0);
                K_dm(i, j) = v;
            }
        }
    }

    // AI-REML for multi-component model
    sigma2_ad_ = 0.3;
    sigma2_am_ = 0.1;
    sigma_adm_ = 0.0;
    sigma2_c_ = 0.0;
    sigma2_e_ = 0.5;

    double vy = (y_.array() - y_.mean()).square().mean();
    sigma2_e_ = vy * 0.4;
    if (sigma2_e_ < 0.01) sigma2_e_ = 0.3;

    // Detect collinear GRMs: if K ≈ K_m, merge into single genetic component
    double grm_diff = (K - K_m).norm() / (K.norm() + 1e-10);
    bool merged_grm = grm_diff < 1e-4;
    if (merged_grm) {
        sigma2_am_ = 0.0;
        sigma2_ad_ = vy * 0.4;
        std::cout << "  Direct and maternal GRMs are near-identical (diff="
                  << grm_diff << "), merging into single genetic component.\n";
    }

    Eigen::MatrixXd X_mat = X_.cast<double>();
    Eigen::MatrixXd I_n = Eigen::MatrixXd::Identity(n_records, n_records);

    double vg_prev = 0.0, ve_prev = 0.0;

    for (int iter = 0; iter < options_.reml_maxit; ++iter) {
        vg_prev = merged_grm ? sigma2_ad_ : (sigma2_ad_ + sigma2_am_);
        ve_prev = sigma2_e_;

        Eigen::MatrixXd V;
        if (merged_grm) {
            V = sigma2_ad_ * K;
        } else {
            V = sigma2_ad_ * K + sigma2_am_ * K_m;
        }
        V.diagonal().array() += sigma2_e_;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::MatrixXd Vinv = ldlt.solve(I_n);
        Eigen::MatrixXd VX = Vinv * X_mat;
        Eigen::MatrixXd XtVX = X_mat.transpose() * VX;
        Eigen::LDLT<Eigen::MatrixXd> xldlt(XtVX);
        if (xldlt.info() != Eigen::Success) break;
        Eigen::MatrixXd C = xldlt.solve(VX.transpose());
        Eigen::MatrixXd P_reml = Vinv - VX * C;

        Eigen::VectorXd Py = P_reml * y_;

        Eigen::VectorXd KPy = K * Py;
        Eigen::VectorXd PKPy = P_reml * KPy;
        Eigen::VectorXd PPy = P_reml * Py;

        // logL_REML = -0.5·[log|V| + log|X'V⁻¹X| + y'Py]
        // Exclude zero/singular D values from X'V⁻¹X logdet (X may be rank-deficient)
        double logdet_V = (ldlt.vectorD().array() > 1e-12).select(ldlt.vectorD().array().log(), 0).sum();
        Eigen::ArrayXd xd = xldlt.vectorD().array();
        double logdet_XVX = (xd > 1e-12).select(xd.log(), 0).sum();
        double old_logL = -0.5 * (logdet_V + logdet_XVX + y_.dot(Py));

        if (merged_grm) {
            // ---- 2-component AI-REML (σ²_g + σ²_e) ----
            Eigen::Matrix2d AI;
            AI(0,0) = 0.5 * KPy.dot(PKPy);
            AI(1,1) = 0.5 * Py.dot(PPy);
            AI(0,1) = AI(1,0) = 0.5 * KPy.dot(PPy);

            double tr_PK = (P_reml * K).trace();
            double tr_P = P_reml.trace();
            Eigen::Vector2d grad;
            grad(0) = -0.5 * tr_PK + 0.5 * KPy.dot(Py);
            grad(1) = -0.5 * tr_P + 0.5 * Py.squaredNorm();

            Eigen::LDLT<Eigen::Matrix2d> ai_ldlt(AI);
            bool ai_ok = (ai_ldlt.info() == Eigen::Success);
            Eigen::Vector2d delta = Eigen::Vector2d::Zero();
            double max_delta_ratio = 0.0;
            if (ai_ok) {
                delta = ai_ldlt.solve(grad);
                max_delta_ratio = std::max(std::abs(delta(0)) / (sigma2_ad_ + 1e-10),
                                          std::abs(delta(1)) / (sigma2_e_ + 1e-10));
                if (!delta.allFinite() || max_delta_ratio > 10.0) ai_ok = false;
            }

            bool accepted = false;
            if (ai_ok) {
                for (int step = 0; step < 12; ++step) {
                    double g_new = sigma2_ad_ + delta(0);
                    double e_new = sigma2_e_ + delta(1);
                    if (g_new > 1e-6 && e_new > 1e-6) {
                        Eigen::MatrixXd V_new = g_new * K;
                        V_new.diagonal().array() += e_new;
                        Eigen::LDLT<Eigen::MatrixXd> new_ldlt(V_new);
                        if (new_ldlt.info() == Eigen::Success) {
                            Eigen::MatrixXd new_Vinv = new_ldlt.solve(I_n);
                            Eigen::MatrixXd new_VX = new_Vinv * X_mat;
                            Eigen::MatrixXd new_XtVX = X_mat.transpose() * new_VX;
                            Eigen::LDLT<Eigen::MatrixXd> new_xldlt(new_XtVX);
                            if (new_xldlt.info() == Eigen::Success) {
                                Eigen::MatrixXd new_C = new_xldlt.solve(new_VX.transpose());
                                Eigen::VectorXd new_Py = (new_Vinv - new_VX * new_C) * y_;
                                double new_logL = -0.5 * (new_ldlt.vectorD().array().log().sum()
                                                  + new_xldlt.vectorD().array().log().sum()
                                                  + y_.dot(new_Py));
                                if (std::isfinite(new_logL) && new_logL >= old_logL - 1e-8) {
                                    sigma2_ad_ = g_new; sigma2_e_ = e_new;
                                    accepted = true;
                                    break;
                                }
                            }
                        }
                    }
                    delta *= 0.5;
                }
            }
            if (!accepted) {
                // EM-REML fallback
                sigma2_ad_ *= KPy.dot(Py) / (tr_PK + 1e-10);
                sigma2_e_ *= Py.squaredNorm() / (tr_P + 1e-10);
                sigma2_e_ = std::max(sigma2_e_, 0.03 * sigma2_ad_);
            }
        } else {
            // ---- 3-component AI-REML (σ²_ad + σ²_am + σ²_e) ----
            Eigen::VectorXd KmPy = K_m * Py;
            Eigen::VectorXd PKmPy = P_reml * KmPy;

            Eigen::Matrix3d AI;
            AI(0,0) = 0.5 * KPy.dot(PKPy);
            AI(1,1) = 0.5 * KmPy.dot(PKmPy);
            AI(2,2) = 0.5 * Py.dot(PPy);
            AI(0,1) = AI(1,0) = 0.5 * KPy.dot(PKmPy);
            AI(0,2) = AI(2,0) = 0.5 * KPy.dot(PPy);
            AI(1,2) = AI(2,1) = 0.5 * KmPy.dot(PPy);

            double tr_PK = (P_reml * K).trace();
            double tr_PKm = (P_reml * K_m).trace();
            double tr_P = P_reml.trace();

            Eigen::Vector3d grad;
            grad(0) = -0.5 * tr_PK + 0.5 * KPy.dot(Py);
            grad(1) = -0.5 * tr_PKm + 0.5 * KmPy.dot(Py);
            grad(2) = -0.5 * tr_P + 0.5 * Py.squaredNorm();

            Eigen::LDLT<Eigen::Matrix3d> ai_ldlt(AI);
            bool ai_ok = (ai_ldlt.info() == Eigen::Success);
            Eigen::Vector3d delta = Eigen::Vector3d::Zero();
            double max_delta_ratio = 0.0;
            if (ai_ok) {
                delta = ai_ldlt.solve(grad);
                max_delta_ratio = std::max({std::abs(delta(0)) / (sigma2_ad_ + 1e-10),
                                           std::abs(delta(1)) / (sigma2_am_ + 1e-10),
                                           std::abs(delta(2)) / (sigma2_e_ + 1e-10)});
                if (!delta.allFinite() || max_delta_ratio > 10.0) ai_ok = false;
            }

            bool accepted = false;
            if (ai_ok) {
                for (int step = 0; step < 12; ++step) {
                    double ad_new = sigma2_ad_ + delta(0);
                    double am_new = sigma2_am_ + delta(1);
                    double e_new = sigma2_e_ + delta(2);
                    if (ad_new > 1e-6 && am_new > 1e-6 && e_new > 1e-6) {
                        Eigen::MatrixXd V_new = ad_new * K + am_new * K_m;
                        V_new.diagonal().array() += e_new;
                        Eigen::LDLT<Eigen::MatrixXd> new_ldlt(V_new);
                        if (new_ldlt.info() == Eigen::Success) {
                            Eigen::MatrixXd new_Vinv = new_ldlt.solve(I_n);
                            Eigen::MatrixXd new_VX = new_Vinv * X_mat;
                            Eigen::MatrixXd new_XtVX = X_mat.transpose() * new_VX;
                            Eigen::LDLT<Eigen::MatrixXd> new_xldlt(new_XtVX);
                            if (new_xldlt.info() == Eigen::Success) {
                                Eigen::MatrixXd new_C = new_xldlt.solve(new_VX.transpose());
                                Eigen::VectorXd new_Py = (new_Vinv - new_VX * new_C) * y_;
                                double new_logL = -0.5 * (new_ldlt.vectorD().array().log().sum()
                                                  + new_xldlt.vectorD().array().log().sum()
                                                  + y_.dot(new_Py));
                                if (std::isfinite(new_logL) && new_logL >= old_logL - 1e-8) {
                                    sigma2_ad_ = ad_new; sigma2_am_ = am_new; sigma2_e_ = e_new;
                                    accepted = true;
                                    break;
                                }
                            }
                        }
                    }
                    delta *= 0.5;
                }
            }
            if (!accepted) {
                // EM-REML fallback
                sigma2_ad_ *= KPy.dot(Py) / (tr_PK + 1e-10);
                sigma2_am_ *= KmPy.dot(Py) / (tr_PKm + 1e-10);
                sigma2_e_ *= Py.squaredNorm() / (tr_P + 1e-10);
                sigma2_e_ = std::max(sigma2_e_, 0.03 * (sigma2_ad_ + sigma2_am_ + sigma2_e_));
            }
        }

        double vg = merged_grm ? sigma2_ad_ : (sigma2_ad_ + sigma2_am_);
        double ve = sigma2_e_;
        {   // clamp
            double v_total = vg + ve;
            if (v_total > 50.0 * vy) { double scale = 50.0 * vy / v_total;
                sigma2_ad_ *= scale; sigma2_am_ *= scale; sigma2_e_ *= scale; }
            if (sigma2_ad_ < 1e-8) sigma2_ad_ = 1e-8;
            if (sigma2_am_ < 1e-8) sigma2_am_ = 1e-8;
            if (sigma2_e_ < 1e-8) sigma2_e_ = 1e-8;
        }

        if (iter % 5 == 0 || iter == options_.reml_maxit - 1) {
            std::cout << "Maternal VCE iter " << iter << ": sig2_ad=" << sigma2_ad_
                      << " sig2_am=" << sigma2_am_ << " sig2_e=" << sigma2_e_
                      << " logL=" << old_logL << "\n";
        }
        if (iter > 1) {
            double vg_new = merged_grm ? sigma2_ad_ : (sigma2_ad_ + sigma2_am_);
            double chg = std::abs(vg_new - vg_prev) / (vg_prev + 1e-10)
                       + std::abs(sigma2_e_ - ve_prev) / (ve_prev + 1e-10);
            if (chg < 1e-8) break;
        }
    }

    std::cout << "Final VCE: Vg(d)=" << sigma2_ad_ << " Vg(m)=" << sigma2_am_
              << " Cov(d,m)=" << sigma_adm_
              << " Vpe=" << sigma2_c_ << " Ve=" << sigma2_e_ << "\n";
}

void MaternalModel::run() {
    buildGRMs();
    runVCE();
    scanGwas();
    if (options_.maternal_snp_test == "maternal" || options_.maternal_snp_test == "joint") {
        scanMaternalSnp();
    }
    writeNullSummary();
}

void MaternalModel::runLocoScan(
    const Eigen::VectorXd& y_agg,
    const Eigen::MatrixXf& X_agg,
    double vg_eff,
    double ve_eff,
    int n_animal,
    const std::string& out_file)
{
    // Determine chromosomes
    std::set<std::string> chrom_labels_set;
    if (is_bgen_) {
        bgen_reader_->loadVariantInfo();
        for (const auto& v : bgen_reader_->getVariants()) chrom_labels_set.insert(v.chrom);
    } else if (is_pgen_) {
        for (const auto& v : pgen_reader_->getVariants()) chrom_labels_set.insert(v.chrom);
    } else {
        for (const auto& snp : plink_reader_->getBimInfo()) chrom_labels_set.insert(snp.chrom);
    }

    bool use_loco = chrom_labels_set.size() > 1 && !options_.no_loco;

    std::vector<std::string> chrom_labels(chrom_labels_set.begin(), chrom_labels_set.end());
    std::sort(chrom_labels.begin(), chrom_labels.end(), [](const std::string& a, const std::string& b) {
        try { return std::stoi(a) < std::stoi(b); }
        catch (...) { return a < b; }
    });

    // Build split GRM for LOCO or use existing G_direct_
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.algorithm = options_.grm_alg;

    Eigen::MatrixXd G_sub;
    GrmBuilder::SplitGrmResult split;
    double S_global = 1.0;
    bool has_split = false;

    if (use_loco) {
        GrmBuilder split_builder(grm_opts);
        if (is_bgen_) {
            bgen_reader_->reset();
            split = split_builder.computeSplit(*bgen_reader_, chrom_labels_set);
        } else if (is_pgen_) {
            pgen_reader_->reset();
            split = split_builder.computeSplit(*pgen_reader_, chrom_labels_set);
        } else {
            split = split_builder.computeSplit(*plink_reader_, chrom_labels_set);
        }
        S_global = (options_.grm_alg == 0) ? (double)split.global_snp_count : split.global_sum_2pq;
        G_sub = DataFilter::subsetSquare(split.global, keep_);
    } else {
        G_sub = DataFilter::subsetSquare(G_direct_, keep_);
        if (G_sub.rows() != n_animal) {
            // If dimension mismatch, use G_direct_ directly
            G_sub = G_direct_;
        }
    }

    // LmmEngine options
    LmmEngine::Options lmm_opts;
    lmm_opts.min_maf = (float)options_.min_maf;
    lmm_opts.min_hwe = (float)options_.min_hwe;
    lmm_opts.fast_mode = options_.fast_lmm;
    lmm_opts.use_spa = options_.use_spa;
    lmm_opts.use_covariate_projection = true;
#ifdef _OPENMP
    lmm_opts.thread_num = options_.n_threads;
#endif

    // Open output file
    std::ofstream ofs(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    ofs.close();

    // Scan lambda: create EigenVce with fixed (effective) VCs → LmmEngine scan
    auto run_scan = [&](const Eigen::MatrixXd& G_work, const std::string& chrom_filter) {
        EigenVce::Options vce_opts;
        vce_opts.use_covariate_projection = true;
        EigenVce vce(vce_opts);
        vce.prepare(y_agg, X_agg, G_work);
        vce.fixedComponents(vg_eff, ve_eff);

        LmmEngine engine(lmm_opts);
        engine.prepareFromVce(vce);

        std::vector<int> all_indices(n_animal);
        for (int i = 0; i < n_animal; ++i) all_indices[i] = i;

        std::vector<GwasResult> results;
        if (is_bgen_) {
            bgen_reader_->reset();
            engine.scan(*bgen_reader_, results, all_indices, chrom_filter);
        } else if (is_pgen_) {
            pgen_reader_->reset();
            engine.scan(*pgen_reader_, results, all_indices, chrom_filter);
        } else {
            engine.scan(*plink_reader_, results, all_indices, chrom_filter);
        }

        std::ofstream fout(out_file, std::ios::app);
        for (const auto& r : results) {
            if (std::isnan(r.p_value)) continue;
            fout << r.chrom << "\t" << r.snp_id << "\t" << r.pos << "\t"
                 << r.a1 << "\t" << r.a2 << "\t" << r.freq << "\t"
                 << r.n_obs << "\t" << r.beta << "\t" << r.se << "\t"
                 << r.p_value << "\t" << r.log10_p << "\n";
        }
    };

    if (!use_loco) {
        run_scan(G_sub, "");
    } else {
        for (size_t ci = 0; ci < split.chrom_labels.size(); ++ci) {
            const std::string& chrom = split.chrom_labels[ci];
            std::cout << "  Maternal LOCO chr " << chrom << "...\n";

            double S_chr = (options_.grm_alg == 0)
                ? (double)split.per_chrom_snp_count[ci]
                : split.per_chrom_sum_2pq[ci];
            double S_loco = S_global - S_chr;

            if (S_chr > 0.0 && S_loco > 1e-6) {
                Eigen::MatrixXd G_chr_sub = DataFilter::subsetSquare(split.per_chrom_raw[ci], keep_);
                Eigen::MatrixXd G_loco(n_animal, n_animal);
#pragma omp parallel for collapse(2)
                for (int i = 0; i < n_animal; ++i) {
                    for (int j = 0; j < n_animal; ++j) {
                        G_loco(i, j) = G_sub(i, j) * (S_global / S_loco)
                                     - G_chr_sub(i, j) / S_loco;
                    }
                }
                run_scan(G_loco, chrom);
            } else {
                run_scan(G_sub, chrom);
            }
        }
    }
}

void MaternalModel::scanGwas() {
    std::cout << "Running maternal model GWAS scan...\n";

    // Compute effective VCs from multi-component model
    double vg_eff = sigma2_ad_ + sigma2_am_;
    if (options_.maternal_pe) vg_eff += sigma2_c_;
    if (options_.direct_maternal_cov) vg_eff += 2.0 * sigma_adm_;
    double ve_eff = sigma2_e_;
    std::cout << "  Effective VCs (from multi-component VCE): Vg=" << vg_eff
              << " Ve=" << ve_eff << "\n";

    int n_records = (int)y_.size();
    int n_animal = (int)keep_.keep_indices.size();

    // Aggregate record-level data to animal level
    std::vector<double> y_sum(n_animal, 0.0);
    std::vector<int> rec_count(n_animal, 0);
    Eigen::MatrixXd X_sum = Eigen::MatrixXd::Zero(n_animal, X_.cols());

    for (int r = 0; r < n_records; ++r) {
        int ai = record_to_animal_[r];
        if (ai < 0 || ai >= n_animal) continue;
        y_sum[ai] += y_(r);
        rec_count[ai]++;
        for (int c = 0; c < X_.cols(); ++c) {
            X_sum(ai, c) += X_(r, c);
        }
    }

    Eigen::VectorXd y_agg(n_animal);
    Eigen::MatrixXf X_agg(n_animal, X_.cols());
    for (int i = 0; i < n_animal; ++i) {
        double cnt = std::max(rec_count[i], 1);
        y_agg(i) = (cnt > 0) ? y_sum[i] / cnt : 0.0;
        for (int c = 0; c < X_.cols(); ++c) {
            X_agg(i, c) = (float)((cnt > 0) ? X_sum(i, c) / cnt : 0.0);
        }
    }

    std::string fname = options_.out_prefix + ".assoc.txt";
    runLocoScan(y_agg, X_agg, vg_eff, ve_eff, n_animal, fname);
    std::cout << "Maternal model results saved to " << fname << "\n";
}

namespace {
// Upper-tail chi-square probability (regularized gamma Q), for the joint 2-df Wald.
double maternal_gamma_q(double a, double x) {
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
double maternal_chisq_p(double stat, int df) {
    if (stat <= 0.0 || df < 1) return 1.0;
    return maternal_gamma_q(df / 2.0, stat / 2.0);
}
// Robust symmetric PD pseudo-inverse via eigenvalue truncation.
// Returns (rank, condition number); zeroed eigenvalues -> dropped directions.
Eigen::MatrixXd maternal_pinv(const Eigen::MatrixXd& A, int* rank_out = nullptr,
                              double* cond_out = nullptr) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    const double tol = 1e-10;
    Eigen::VectorXd ev = es.eigenvalues();
    double e_max = ev.size() > 0 ? std::max(std::abs(ev[ev.size() - 1]), 1e-300) : 1e-300;
    double cutoff = tol * e_max;
    int rank = 0;
    for (int i = 0; i < ev.size(); ++i) if (ev[i] > cutoff) rank++;
    double e_min_kept = e_max, e_max_kept = e_max;
    for (int i = 0; i < ev.size(); ++i) {
        if (ev[i] > cutoff) e_min_kept = std::min(e_min_kept, ev[i]);
    }
    Eigen::VectorXd einv = Eigen::VectorXd::Zero(ev.size());
    for (int i = 0; i < ev.size(); ++i) {
        if (ev[i] > cutoff) einv[i] = 1.0 / ev[i];
    }
    if (rank_out) *rank_out = rank;
    if (cond_out) *cond_out = e_max_kept / e_min_kept;
    return es.eigenvectors() * einv.asDiagonal() * es.eigenvectors().transpose();
}
} // namespace

void MaternalModel::scanMaternalSnp() {
    if (is_bgen_ || is_pgen_) {
        throw std::runtime_error("--maternal-snp-test maternal|joint currently requires PLINK .bed input (BGEN/PGEN path not yet wired).");
    }
    const bool join = (options_.maternal_snp_test == "joint");
    const int n_samples = n_samples_;

    // Effective VCs from the (standard) maternal VCE; fixed, shared across modes.
    double vg_eff = sigma2_ad_ + sigma2_am_;
    if (options_.maternal_pe) vg_eff += sigma2_c_;
    if (options_.direct_maternal_cov) vg_eff += 2.0 * sigma_adm_;
    double ve_eff = sigma2_e_;

    // Subject-indexed animal-level aggregation (independent of record/keep ordering):
    // every GRM in this scan is indexed by genomic subject id (0..n_samples-1), and the
    // eligible set is chosen purely from each subject's dam id.
    int n_records = (int)y_.size();
    std::vector<double> y_sum(n_samples, 0.0);
    std::vector<int> rec_count(n_samples, 0);
    Eigen::MatrixXd X_sum = Eigen::MatrixXd::Zero(n_samples, X_.cols());
    std::vector<int> animal_dam(n_samples, -1);
    for (int r = 0; r < n_records; ++r) {
        int ai = record_to_animal_[r];
        if (ai < 0 || ai >= n_samples) continue;
        y_sum[ai] += y_(r);
        rec_count[ai]++;
        for (int c = 0; c < X_.cols(); ++c) X_sum(ai, c) += X_(r, c);
        int di = record_to_dam_[r];
        if (animal_dam[ai] < 0) animal_dam[ai] = di;
    }

    // Eligible subset: subject has records AND its dam is a genotyped subject (< n_samples).
    std::vector<int> eligible_subj;   // subject index of the animal
    std::vector<int> eligible_dam;    // subject index of the animal's dam
    int n_total = 0, n_with_dam = 0;
    for (int s = 0; s < n_samples; ++s) {
        if (rec_count[s] == 0) continue;
        n_total++;
        int di = animal_dam[s];
        if (di < 0 || di >= n_samples) continue; // missing dam genotype
        eligible_subj.push_back(s);
        eligible_dam.push_back(di);
        n_with_dam++;
    }
    int n_dropped = n_total - n_with_dam;
    int n_el = (int)eligible_subj.size();
    if (n_el < 10) throw std::runtime_error("Too few animals with dam genotypes for maternal/joint SNP test.");

    Eigen::VectorXd y_sub(n_el);
    Eigen::MatrixXf X_sub(n_el, X_.cols());
    for (int e = 0; e < n_el; ++e) {
        int s = eligible_subj[e];
        double cnt = std::max(rec_count[s], 1);
        y_sub(e) = y_sum[s] / cnt;
        for (int c = 0; c < X_.cols(); ++c) X_sub(e, c) = (float)(X_sum(s, c) / cnt);
    }

    // Output files
    std::string out_file = options_.out_prefix + (join ? ".joint.assoc.txt" : ".maternal.assoc.txt");
    std::string header;
    if (join) {
        header = "CHR\tSNP\tBP\tA1\tA2\tFREQ_OFF\tFREQ_DAM\tN"
                 "\tBETA_DIRECT\tSE_DIRECT\tBETA_MATERNAL\tSE_MATERNAL"
                 "\tCOV_BETA_DM\tCORR_BETA_DM"
                 "\tSTAT_JOINT\tDF\tP_JOINT\tLOG10P_JOINT"
                 "\tCOND\tRANK\tSTATUS"
                 "\tN_TOTAL\tN_WITH_DAM_GENOTYPE\tN_DROPPED_MISSING_DAM\n";
    } else {
        header = "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN"
                 "\tBETA_MATERNAL\tSE_MATERNAL\tP\tLOG10P"
                 "\tN_TOTAL\tN_WITH_DAM_GENOTYPE\tN_DROPPED_MISSING_DAM\n";
    }
    {
        std::ofstream ofs(out_file);
        ofs << header;
    }
    std::string meta_file = options_.out_prefix + ".maternal.snp.meta.txt";
    {
        std::ofstream ofs(meta_file);
        ofs << "mode\t" << (join ? "joint" : "maternal") << "\n";
        ofs << "n_total\t" << n_total << "\n";
        ofs << "n_with_dam_genotype\t" << n_with_dam << "\n";
        ofs << "n_dropped_missing_dam\t" << n_dropped << "\n";
        ofs << "n_eligible_animals\t" << n_el << "\n";
        ofs << "null_covariance\tfull_component_maternal\n";
        ofs << "null_components\tsigma2_ad*K_offspring+sigma2_am*K_dam+sigma2_e*I\n";
        ofs << "vg_eff\t" << vg_eff << "\n";
        ofs << "ve_eff\t" << ve_eff << "\n";
        ofs << "sigma2_ad\t" << sigma2_ad_ << "\n";
        ofs << "sigma2_am\t" << sigma2_am_ << "\n";
        ofs << "sigma_adm\t" << sigma_adm_ << "\n";
        ofs << "sigma2_c\t" << sigma2_c_ << "\n";
        ofs << "sigma2_e\t" << sigma2_e_ << "\n";
        ofs << "missing_dam_strategy\tfilter_to_eligible_subset\n";
        ofs << "per_snp_missing\tmean_imputed_in_design_column\n";
    }

    std::cout << "Maternal/joint SNP test: mode=" << (join ? "joint" : "maternal")
              << " N_total=" << n_total << " N_with_dam=" << n_with_dam
              << " N_dropped=" << n_dropped << "\n";

    // Chromosome set + LOCO GRM (same policy as direct scan).
    std::set<std::string> chrom_labels_set;
    for (const auto& snp : plink_reader_->getBimInfo()) chrom_labels_set.insert(snp.chrom);
    bool use_loco = chrom_labels_set.size() > 1 && !options_.no_loco;
    std::vector<std::string> chrom_labels(chrom_labels_set.begin(), chrom_labels_set.end());
    std::sort(chrom_labels.begin(), chrom_labels.end(), [](const std::string& a, const std::string& b) {
        try { return std::stoi(a) < std::stoi(b); }
        catch (...) { return a < b; }
    });

    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.algorithm = options_.grm_alg;

    Eigen::MatrixXd G_sub_all;
    GrmBuilder::SplitGrmResult split;
    double S_global = 1.0;
    if (use_loco) {
        GrmBuilder split_builder(grm_opts);
        split = split_builder.computeSplit(*plink_reader_, chrom_labels_set);
        S_global = (options_.grm_alg == 0) ? (double)split.global_snp_count : split.global_sum_2pq;
        G_sub_all = split.global;   // subject-indexed full GRM (no keep subset)
    } else {
        G_sub_all = G_direct_;
    }

    // Scan lambda: G_work is a subject-indexed GRM (full, or LOCO-corrected for one chromosome).
    // Null covariance mirrors the maternal VCE: σ²ad·K(offspring) + σ²am·K(dams) + σ²e·I.
    // Factored ONCE per chromosome, then every SNP reuses it.
    auto run_scan = [&](const Eigen::MatrixXd& G_work, const std::string& chrom_filter) {
        Eigen::MatrixXd G_dir_el(n_el, n_el);
        Eigen::MatrixXd G_mat_el(n_el, n_el);
        Eigen::MatrixXd G_dm_el;
        if (options_.direct_maternal_cov) G_dm_el.resize(n_el, n_el);
        for (int a = 0; a < n_el; ++a) {
            int sa = eligible_subj[a];
            int da = eligible_dam[a];
            for (int b = 0; b < n_el; ++b) {
                G_dir_el(a, b) = G_work(sa, eligible_subj[b]);
                G_mat_el(a, b) = G_work(da, eligible_dam[b]);
                if (options_.direct_maternal_cov)
                    G_dm_el(a, b) = G_work(sa, eligible_dam[b]);
            }
        }
        Eigen::MatrixXd V = sigma2_ad_ * G_dir_el;
        if (sigma2_am_ > 0.0) V.noalias() += sigma2_am_ * G_mat_el;
        if (options_.direct_maternal_cov && std::abs(sigma_adm_) > 0.0)
            V.noalias() += sigma_adm_ * (G_dm_el + G_dm_el.transpose());
        V.diagonal().array() += sigma2_e_;
        if (options_.maternal_pe && sigma2_c_ > 0.0) V.diagonal().array() += sigma2_c_;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(V);
        if (ldlt.info() != Eigen::Success) {
            std::cerr << "Warning: null covariance factorization failed; skipping chromosome " << chrom_filter << "\n";
            return;
        }
        Eigen::MatrixXd Vinv = ldlt.solve(Eigen::MatrixXd::Identity(n_el, n_el));
        Eigen::MatrixXd X_mat = X_sub.cast<double>();
        const int k2 = (int)X_mat.cols();
        Eigen::VectorXd Py;
        Eigen::MatrixXd VinvX;
        Eigen::MatrixXd XtVx_inv;
        if (k2 > 0) {
            VinvX = Vinv * X_mat;
            Eigen::MatrixXd XtVx = X_mat.transpose() * VinvX;
            XtVx_inv = maternal_pinv(XtVx);
            Eigen::VectorXd XtVinvY = X_mat.transpose() * (Vinv * y_sub);
            Eigen::VectorXd coefX = XtVx_inv * XtVinvY;
            Py = Vinv * y_sub - VinvX * coefX;
        } else {
            Py = Vinv * y_sub;
        }

        const auto& bim = plink_reader_->getBimInfo();
        const uint8_t* raw_data = plink_reader_->getData();
        const size_t bytes_per_snp = plink_reader_->getBytesPerSnp();
        const int m = (int)plink_reader_->getNumSnps();

        std::ofstream fout(out_file, std::ios::app);
        for (int j = 0; j < m; ++j) {
            if (!chrom_filter.empty() && bim[j].chrom != chrom_filter) continue;

            const uint8_t* snp_data = raw_data + 3 + (size_t)j * bytes_per_snp;
            auto dose = [&](int orig_idx) {
                uint8_t byte = snp_data[orig_idx >> 2];
                uint8_t code = (byte >> ((orig_idx & 3) * 2)) & 3;
                if (code == 1) return std::numeric_limits<double>::quiet_NaN();
                return (code == 0) ? 2.0 : (code == 2 ? 1.0 : 0.0);
            };

            // Offspring + dam design columns over the eligible subset.
            Eigen::VectorXd s_off(n_el), s_dam(n_el);
            double sum_off = 0.0, sum_dam = 0.0;
            int v_off = 0, v_dam = 0;
            for (int e = 0; e < n_el; ++e) {
                double vo = dose(eligible_subj[e]);
                double vd = dose(eligible_dam[e]);
                s_off(e) = vo; s_dam(e) = vd;
                if (!std::isnan(vo)) { sum_off += vo; v_off++; }
                if (!std::isnan(vd)) { sum_dam += vd; v_dam++; }
            }
            double freq_off = (v_off > 0) ? sum_off / (2.0 * v_off) : 0.0;
            double freq_dam = (v_dam > 0) ? sum_dam / (2.0 * v_dam) : 0.0;
            double maf_check = join ? freq_off : freq_dam;
            if ((join ? v_off : v_dam) < 3 || std::min(maf_check, 1.0 - maf_check) < options_.min_maf) continue;

            // Center; per-SNP missingness mean-imputed (0 after centering), matching LMM engine.
            double mean_off = 2.0 * freq_off, mean_dam = 2.0 * freq_dam;
            for (int e = 0; e < n_el; ++e) {
                s_off(e) = std::isnan(s_off(e)) ? 0.0 : (s_off(e) - mean_off);
                s_dam(e) = std::isnan(s_dam(e)) ? 0.0 : (s_dam(e) - mean_dam);
            }

            int p = join ? 2 : 1;
            Eigen::MatrixXd S(n_el, p);
            S.col(0) = join ? s_off : s_dam;              // maternal mode tests the DAM genotype
            if (join) S.col(1) = s_dam;

            Eigen::VectorXd sc = S.transpose() * Py;      // p-vector score
            Eigen::MatrixXd VinvS = Vinv * S;
            Eigen::MatrixXd M = S.transpose() * VinvS;
            if (k2 > 0) {
                M -= (S.transpose() * VinvX) * (XtVx_inv * (X_mat.transpose() * VinvS));
            }
            M = 0.5 * (M + M.transpose());

            // Maternal (single-column) mode: mirror the LMM engine and skip non-testable SNPs.
            if (!join && !(M(0, 0) > 1e-12)) continue;

            int rank = 0; double cond = 0.0;
            Eigen::MatrixXd Minv = maternal_pinv(M, &rank, &cond);
            Eigen::VectorXd beta = Minv * sc;
            double stat = sc.dot(beta);
            int df = std::max(1, rank);
            double p_val = maternal_chisq_p(stat, df);
            std::string status = (rank < p) ? "COLLINEAR" : "OK";
            if (cond > 1e6) status = "COLLINEAR";

            if (join) {
                double se_d = std::sqrt(std::max(Minv(0, 0), 0.0));
                double se_m = std::sqrt(std::max(Minv(1, 1), 0.0));
                double cov_dm = Minv(0, 1);
                double corr_dm = (se_d > 1e-300 && se_m > 1e-300) ? cov_dm / (se_d * se_m) : 0.0;
                fout << bim[j].chrom << "\t" << bim[j].id << "\t" << bim[j].bp_pos << "\t"
                     << bim[j].alt << "\t" << bim[j].ref << "\t"
                     << freq_off << "\t" << freq_dam << "\t" << n_el << "\t"
                     << beta(0) << "\t" << se_d << "\t"
                     << beta(1) << "\t" << se_m << "\t"
                     << cov_dm << "\t" << corr_dm << "\t"
                     << stat << "\t" << df << "\t" << p_val << "\t"
                     << (p_val > 0 ? -std::log10(p_val) : 300.0) << "\t"
                     << cond << "\t" << rank << "\t" << status << "\t"
                     << n_total << "\t" << n_with_dam << "\t" << n_dropped << "\n";
            } else {
                double se_m = std::sqrt(std::max(Minv(0, 0), 0.0));
                fout << bim[j].chrom << "\t" << bim[j].id << "\t" << bim[j].bp_pos << "\t"
                     << bim[j].alt << "\t" << bim[j].ref << "\t"
                     << freq_dam << "\t" << n_el << "\t"
                     << beta(0) << "\t" << se_m << "\t" << p_val << "\t"
                     << (p_val > 0 ? -std::log10(p_val) : 300.0) << "\t"
                     << n_total << "\t" << n_with_dam << "\t" << n_dropped << "\n";
            }
        }
    };

    if (!use_loco) {
        run_scan(G_sub_all, "");
    } else {
        for (size_t ci = 0; ci < split.chrom_labels.size(); ++ci) {
            const std::string& chrom = split.chrom_labels[ci];
            double S_chr = (options_.grm_alg == 0) ? (double)split.per_chrom_snp_count[ci]
                                                   : split.per_chrom_sum_2pq[ci];
            double S_loco = S_global - S_chr;
            std::cout << "  Maternal " << (join ? "joint" : "maternal") << " LOCO chr " << chrom << "...\n";
            if (S_chr > 0.0 && S_loco > 1e-6) {
                Eigen::MatrixXd G_chr_sub = split.per_chrom_raw[ci];   // subject-indexed full
                Eigen::MatrixXd G_loco(n_samples, n_samples);
#pragma omp parallel for collapse(2)
                for (int i = 0; i < n_samples; ++i) {
                    for (int j = 0; j < n_samples; ++j) {
                        G_loco(i, j) = G_sub_all(i, j) * (S_global / S_loco) - G_chr_sub(i, j) / S_loco;
                    }
                }
                run_scan(G_loco, chrom);
            } else {
                run_scan(G_sub_all, chrom);
            }
        }
    }
    std::cout << "Maternal " << (join ? "joint" : "maternal")
              << " SNP-test results saved to " << out_file << "\n";
}

void MaternalModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_maternal_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "sigma2_ad\tsigma2_am\tsigma_adm\tsigma2_c\tsigma2_e\n";
    ofs << sigma2_ad_ << "\t" << sigma2_am_ << "\t"
        << sigma_adm_ << "\t" << sigma2_c_ << "\t"
        << sigma2_e_ << "\n";
    ofs.close();
    std::cout << "Maternal null model summary written to " << fname << "\n";
}

} // namespace cosmic
