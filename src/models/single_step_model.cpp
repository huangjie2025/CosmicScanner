#include "models/single_step_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "eigen_vce.h"

namespace {

std::unordered_set<std::string> loadSnpIdSet(const std::string& path) {
    std::unordered_set<std::string> ids;
    if (path.empty()) return ids;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open SNP list: " + path);
    std::string id;
    while (ifs >> id) ids.insert(id);
    if (ids.empty()) throw std::runtime_error("SNP list is empty: " + path);
    return ids;
}

} // anonymous namespace

namespace cosmic {

SingleStepModel::Options SingleStepModel::Options::fromConfig(const GwasConfig& cfg) {
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
    opts.ped_file = cfg.ped_file;
    opts.genotyped_list = cfg.genotyped_list;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.h_tau = cfg.h_tau;
    opts.h_omega = cfg.h_omega;
    opts.blend_g = cfg.blend_g;
    opts.write_ebv = cfg.write_ebv;
    return opts;
}

SingleStepModel::SingleStepModel(const Options& opts) : options_(opts) {}

void SingleStepModel::prepare(const std::vector<std::string>& fids,
                              const std::vector<std::string>& iids,
                              int n_samples,
                              PlinkReader& plink_reader,
                              BgenReader& bgen_reader,
                              PgenReader& pgen_reader,
                              bool is_bgen, bool is_pgen) {
    fids_ = fids;
    iids_ = iids;
    n_samples_ = n_samples;
    plink_reader_ = &plink_reader;
    bgen_reader_ = &bgen_reader;
    pgen_reader_ = &pgen_reader;
    is_bgen_ = is_bgen;
    is_pgen_ = is_pgen;

    loadPhenotype();
    loadCovariates();

    DataFilter filter;
    keep_ = filter.filterSamples(y_, X_);
    std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for ssGBLUP.");
}

void SingleStepModel::loadPhenotype() {
    y_ = DataLoader::loadPhenotype(options_.pheno_file, options_.mpheno, fids_, iids_, n_samples_, *plink_reader_, is_bgen_, is_pgen_);
}

void SingleStepModel::loadCovariates() {
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_, n_samples_);
}

void SingleStepModel::loadPedigree() {
    std::string file = options_.ped_file;
    if (file.empty()) throw std::runtime_error("Single-step model requires --ped file.");

    std::ifstream ifs(file);
    if (!ifs.is_open()) throw std::runtime_error("Cannot open pedigree file: " + file);

    std::string line;
    int idx = 0;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string id, sire, dam;
        if (!(iss >> id >> sire >> dam)) continue;
        if (sire == "NA" || sire == ".") sire = "0";
        if (dam == "NA" || dam == ".") dam = "0";
        pedigree_.push_back({id, sire, dam});
        ped_id_to_idx_[id] = idx++;
    }

    ifs.close();
    std::cout << "Loaded pedigree with " << pedigree_.size() << " individuals.\n";

    if (!options_.genotyped_list.empty()) {
        std::ifstream gf(options_.genotyped_list);
        std::string gid;
        while (gf >> gid) {
            auto it = ped_id_to_idx_.find(gid);
            if (it != ped_id_to_idx_.end()) {
                genotyped_indices_.push_back(it->second);
            }
        }
        n_genotyped_ = (int)genotyped_indices_.size();
        std::cout << n_genotyped_ << " individuals have genotypes.\n";
    }
}

void SingleStepModel::buildAInverse() {
    int n_ped = (int)pedigree_.size();
    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<double> d(n_ped, 0.0);

    for (int i = 0; i < n_ped; ++i) {
        const auto& p = pedigree_[i];
        int si = -1, di = -1;
        if (p.sire != "0") {
            auto sit = ped_id_to_idx_.find(p.sire);
            if (sit != ped_id_to_idx_.end()) si = sit->second;
        }
        if (p.dam != "0") {
            auto dit = ped_id_to_idx_.find(p.dam);
            if (dit != ped_id_to_idx_.end()) di = dit->second;
        }

        if (si < 0 && di < 0) {
            d[i] = 1.0;
        } else if (si >= 0 && di >= 0) {
            d[i] = 2.0;
        } else {
            d[i] = 4.0 / 3.0;
        }
    }

    for (int i = 0; i < n_ped; ++i) {
        const auto& p = pedigree_[i];
        int si = -1, di = -1;
        if (p.sire != "0") {
            auto sit = ped_id_to_idx_.find(p.sire);
            if (sit != ped_id_to_idx_.end()) si = sit->second;
        }
        if (p.dam != "0") {
            auto dit = ped_id_to_idx_.find(p.dam);
            if (dit != ped_id_to_idx_.end()) di = dit->second;
        }

        double alpha = d[i];

        if (si < 0 && di < 0) {
            triplets.emplace_back(i, i, alpha);
        } else if (si >= 0 && di >= 0) {
            double diag_val = 2.0 / alpha;
            triplets.emplace_back(i, i, diag_val);
            triplets.emplace_back(i, si, -1.0 / alpha);
            triplets.emplace_back(i, di, -1.0 / alpha);
            triplets.emplace_back(si, i, -1.0 / alpha);
            triplets.emplace_back(di, i, -1.0 / alpha);
            triplets.emplace_back(si, si, 0.25 / alpha);
            triplets.emplace_back(di, di, 0.25 / alpha);
            triplets.emplace_back(si, di, 0.5 / alpha);
            triplets.emplace_back(di, si, 0.5 / alpha);
        } else {
            int parent = (si >= 0) ? si : di;
            double diag_val = 4.0 / (3.0 * alpha);
            triplets.emplace_back(i, i, diag_val);
            triplets.emplace_back(i, parent, -2.0 / (3.0 * alpha));
            triplets.emplace_back(parent, i, -2.0 / (3.0 * alpha));
            triplets.emplace_back(parent, parent, 1.0 / (3.0 * alpha));
        }
    }

    A_inv_.resize(n_ped, n_ped);
    A_inv_.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "A⁻¹ built with " << triplets.size() << " non-zero entries.\n";
}

void SingleStepModel::buildGRM() {
    std::cout << "Building genomic GRM (G)...\n";
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

    G_ = grm_builder.getGrm();
    n_genotyped_ = G_.rows();
    if (genotyped_indices_.empty()) {
        genotyped_indices_.resize(n_genotyped_);
        for (int i = 0; i < n_genotyped_; ++i) genotyped_indices_[i] = i;
    }
    std::cout << "G matrix: " << n_genotyped_ << "x" << n_genotyped_ << "\n";

    if (options_.blend_g > 0.0) {
        double alpha = options_.blend_g;
        G_ = (1.0 - alpha) * G_ + alpha * Eigen::MatrixXd::Identity(n_genotyped_, n_genotyped_);
        std::cout << "G blended with identity (alpha=" << alpha << ")\n";
    }
}

void SingleStepModel::computeA22() {
    std::cout << "A22 correction will be fully implemented in v1.1.\n";
    std::cout << "MVP: skipping A22 correction, using G⁻¹ + A⁻¹ block structure.\n";
}

void SingleStepModel::computeHInverse() {
    int n_total = (int)pedigree_.size();
    H_inv_ = A_inv_;

    if (n_genotyped_ > 0) {
        Eigen::MatrixXd Ginv = G_.inverse();
        for (int gi = 0; gi < n_genotyped_; ++gi) {
            int idx_i = genotyped_indices_[gi];
            for (int gj = 0; gj < n_genotyped_; ++gj) {
                int idx_j = genotyped_indices_[gj];
                H_inv_.coeffRef(idx_i, idx_j) += Ginv(gi, gj);
            }
        }
    }

    std::cout << "H⁻¹ constructed. Size: " << H_inv_.rows() << "x" << H_inv_.cols() << "\n";
}

void SingleStepModel::runVCE() {
    std::cout << "Running ssGBLUP VCE...\n";

    int n = keep_.n_keep;
    int n_ped = (int)pedigree_.size();
    Eigen::MatrixXd X_mat = X_.cast<double>();
    Eigen::VectorXd y_adj = y_;

    Eigen::MatrixXd M = Eigen::MatrixXd::Identity(n, n) - X_mat * (X_mat.transpose() * X_mat).ldlt().solve(X_mat.transpose());
    Eigen::VectorXd y_proj = M * y_adj;

    double vy = y_proj.squaredNorm() / (n - X_mat.cols());
    sigma2_g_ = vy * 0.5;
    sigma2_e_ = vy * 0.5;

    for (int iter = 0; iter < 20; ++iter) {
        Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_inv_);
        Eigen::MatrixXd H = H_dense.inverse();

        // Using the genotyped subset for VCE
        int ng = n_genotyped_;
        if (ng == 0) ng = n;
        Eigen::MatrixXd G_sub = H.topLeftCorner(ng, ng);

        Eigen::MatrixXd V = sigma2_g_ * G_sub + sigma2_e_ * Eigen::MatrixXd::Identity(ng, ng);
        Eigen::MatrixXd Vinv = V.ldlt().solve(Eigen::MatrixXd::Identity(ng, ng));

        double tr_PG = (Vinv * G_sub).trace();
        Eigen::VectorXd yg = y_adj.head(ng);
        Eigen::VectorXd t = Vinv * yg;
        double yPGy = t.dot(G_sub * t);

        double grad_g = -0.5 * tr_PG + 0.5 * yPGy;
        double grad_e = -0.5 * (double)ng + 0.5 * y_adj.head(ng).dot(Vinv * y_adj.head(ng));

        double step = 0.3;
        sigma2_g_ += step * grad_g;
        sigma2_e_ += step * grad_e;

        if (sigma2_g_ < 1e-6) sigma2_g_ = 1e-6;
        if (sigma2_e_ < 1e-6) sigma2_e_ = 1e-6;

        if (iter % 5 == 0) {
            std::cout << "ssVCE iter " << iter << ": Vg=" << sigma2_g_
                      << " Ve=" << sigma2_e_ << "\n";
        }
    }

    std::cout << "Final ssVCE: Vg=" << sigma2_g_ << " Ve=" << sigma2_e_ << "\n";
}

void SingleStepModel::scanGwas() {
    std::cout << "Running ssGBLUP-informed GWAS scan (using H⁻¹ VCs)...\n";
    std::cout << "  Fixed VCs: Vg=" << sigma2_g_ << " Ve=" << sigma2_e_ << "\n";

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

    // Subset y and X
    Eigen::VectorXd y_sub = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X_, keep_);
    int n_sub = keep_.n_keep;

    // Build split GRMs for LOCO or global GRM
    GrmBuilder::Options grm_opts;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
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
        has_split = true;
    } else {
        G_sub = DataFilter::subsetSquare(G_, keep_);
    }

    // LmmEngine options
    LmmEngine::Options lmm_opts;
    lmm_opts.min_maf = (float)options_.min_maf;
    lmm_opts.min_hwe = (float)options_.min_hwe;
    lmm_opts.fast_mode = options_.fast_lmm;
    lmm_opts.use_spa = options_.use_spa;
    lmm_opts.use_covariate_projection = true;
    lmm_opts.extract_set = loadSnpIdSet(options_.extract_file);
    lmm_opts.exclude_set = loadSnpIdSet(options_.exclude_file);
#ifdef _OPENMP
    lmm_opts.thread_num = options_.n_threads;
#endif

    // Open output file
    std::string out_file = options_.out_prefix + ".assoc.txt";
    std::ofstream ofs(out_file);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    ofs.close();

    // Scan lambda: create EigenVce with fixed VCs → LmmEngine scan → append
    auto run_scan = [&](const Eigen::MatrixXd& G_work, const std::string& chrom_filter) {
        EigenVce::Options vce_opts;
        vce_opts.use_covariate_projection = true;
        EigenVce vce(vce_opts);
        vce.prepare(y_sub, X_sub, G_work);
        vce.fixedComponents(sigma2_g_, sigma2_e_);

        LmmEngine engine(lmm_opts);
        engine.prepareFromVce(vce);

        std::vector<GwasResult> results;
        if (is_bgen_) {
            bgen_reader_->reset();
            engine.scan(*bgen_reader_, results, keep_.keep_indices, chrom_filter);
        } else if (is_pgen_) {
            pgen_reader_->reset();
            engine.scan(*pgen_reader_, results, keep_.keep_indices, chrom_filter);
        } else {
            engine.scan(*plink_reader_, results, keep_.keep_indices, chrom_filter);
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
        std::cout << "Using standard LMM with ssGBLUP VCs (no LOCO)...\n";
        run_scan(G_sub, "");
    } else {
        std::cout << "Running LOCO scan using ssGBLUP VCs on "
                  << split.chrom_labels.size() << " chromosomes...\n";
        for (size_t ci = 0; ci < split.chrom_labels.size(); ++ci) {
            const std::string& chrom = split.chrom_labels[ci];
            std::cout << "  Processing LOCO chr " << chrom << "...\n";

            double S_chr = (options_.grm_alg == 0)
                ? (double)split.per_chrom_snp_count[ci]
                : split.per_chrom_sum_2pq[ci];
            double S_loco = S_global - S_chr;

            if (S_chr > 0.0 && S_loco > 1e-6) {
                Eigen::MatrixXd G_chr_sub = DataFilter::subsetSquare(split.per_chrom_raw[ci], keep_);
                Eigen::MatrixXd G_loco(n_sub, n_sub);
#pragma omp parallel for collapse(2)
                for (int i = 0; i < n_sub; ++i) {
                    for (int j = 0; j < n_sub; ++j) {
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
    std::cout << "ssGBLUP results saved to " << out_file << "\n";
}

void SingleStepModel::writeEbv() {
    if (!options_.write_ebv) return;
    std::cout << "GEBV output will be implemented in v1.1.\n";
}

void SingleStepModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_ssgwas_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "sigma2_g\tsigma2_e\th2\n";
    double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);
    ofs << sigma2_g_ << "\t" << sigma2_e_ << "\t" << h2 << "\n";
    ofs.close();
    std::cout << "ssGBLUP null model summary written to " << fname << "\n";
}

void SingleStepModel::run() {
    loadPedigree();
    buildAInverse();
    buildGRM();
    computeA22();
    computeHInverse();
    runVCE();
    scanGwas();
    writeEbv();
    writeNullSummary();
}

} // namespace cosmic
