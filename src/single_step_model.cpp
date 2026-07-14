#include "single_step_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "eigen_vce.h"
#include "a_matrix.h"
#include "matrix_io.h"

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

Eigen::MatrixXd invertSpdWithRidge(const Eigen::MatrixXd& matrix,
                                   const std::string& label,
                                   double ridge = 1e-8) {
    if (matrix.rows() != matrix.cols()) {
        throw std::runtime_error(label + " must be square.");
    }
    Eigen::MatrixXd sym = 0.5 * (matrix + matrix.transpose());
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(sym.rows(), sym.cols());
    for (int attempt = 0; attempt < 6; ++attempt) {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(sym + ridge * I);
        if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
            return ldlt.solve(I);
        }
        ridge *= 10.0;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(sym + ridge * I);
    if (!lu.isInvertible()) {
        throw std::runtime_error("Cannot invert " + label + " for ssGWAS H-inverse construction.");
    }
    return lu.inverse();
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
    opts.make_agh = cfg.make_agh;
    opts.write_agh = cfg.write_agh;
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

    if (options_.make_agh) {
        keep_.keep_indices.resize(n_samples_);
        std::iota(keep_.keep_indices.begin(), keep_.keep_indices.end(), 0);
        keep_.n_keep = n_samples_;
    } else {
        loadPhenotype();
        loadCovariates();
        DataFilter filter;
        keep_ = filter.filterSamples(y_, X_);
        std::cout << "After filtering, " << keep_.n_keep << " individuals remain.\n";
        if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for ssGWAS.");
    }
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
        if (ped_id_to_idx_.find(id) != ped_id_to_idx_.end()) {
            throw std::runtime_error("Duplicate pedigree ID: " + id);
        }
        pedigree_.push_back({id, sire, dam});
        ped_id_to_idx_[id] = idx++;
    }

    ifs.close();
    std::cout << "Loaded pedigree with " << pedigree_.size() << " individuals.\n";

    if (!options_.genotyped_list.empty()) {
        std::ifstream gf(options_.genotyped_list);
        std::string gid;
        while (gf >> gid) {
            if (!requested_genotyped_ids_.insert(gid).second) {
                throw std::runtime_error("Duplicate ID in --genotyped-list: " + gid);
            }
            if (ped_id_to_idx_.find(gid) == ped_id_to_idx_.end()) {
                throw std::runtime_error("Genotyped ID is absent from pedigree: " + gid);
            }
        }
        std::cout << requested_genotyped_ids_.size() << " requested genotyped IDs loaded.\n";
    }
}

void SingleStepModel::buildAInverse() {
    int n_ped = (int)pedigree_.size();
    std::vector<int> dam(n_ped, -1), sire(n_ped, -1);
    for (int i = 0; i < n_ped; ++i) {
        const auto& p = pedigree_[i];
        if (p.sire != "0") {
            auto sit = ped_id_to_idx_.find(p.sire);
            if (sit == ped_id_to_idx_.end()) throw std::runtime_error("Unknown sire ID in pedigree: " + p.sire);
            sire[i] = sit->second;
        }
        if (p.dam != "0") {
            auto dit = ped_id_to_idx_.find(p.dam);
            if (dit == ped_id_to_idx_.end()) throw std::runtime_error("Unknown dam ID in pedigree: " + p.dam);
            dam[i] = dit->second;
        }
    }
    A_inv_ = genmatrix::invertA_henderson_fixed(dam, sire);
    std::cout << "A⁻¹ built with " << A_inv_.nonZeros() << " non-zero entries.\n";
}

void SingleStepModel::buildGRM() {
    std::cout << "Building genomic GRM (G)...\n";
    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.extract_set = loadSnpIdSet(options_.extract_file);
    grm_opts.exclude_set = loadSnpIdSet(options_.exclude_file);
    grm_opts.thread_num = options_.n_threads;
    GrmBuilder grm_builder(grm_opts);

    if (is_bgen_)
        grm_builder.compute(*bgen_reader_, keep_.keep_indices);
    else if (is_pgen_)
        grm_builder.compute(*pgen_reader_, keep_.keep_indices);
    else
        grm_builder.compute(*plink_reader_);

    G_ = (is_bgen_ || is_pgen_)
        ? grm_builder.getGrm()
        : DataFilter::subsetSquare(grm_builder.getGrm(), keep_);
    n_genotyped_ = G_.rows();
    genotyped_indices_.clear();
    genotyped_ids_.clear();
    std::unordered_set<int> mapped_pedigree_rows;
    for (int sample_index : keep_.keep_indices) {
        const std::string fid_iid = fids_[sample_index] + ":" + iids_[sample_index];
        const std::string& iid = iids_[sample_index];
        auto it = ped_id_to_idx_.find(fid_iid);
        std::string matched_id = fid_iid;
        if (it == ped_id_to_idx_.end()) {
            it = ped_id_to_idx_.find(iid);
            matched_id = iid;
        }
        if (it == ped_id_to_idx_.end()) {
            throw std::runtime_error("Genotype sample is absent from pedigree: " + fid_iid);
        }
        if (!requested_genotyped_ids_.empty() && requested_genotyped_ids_.find(matched_id) == requested_genotyped_ids_.end()) {
            throw std::runtime_error("Genotype sample is absent from --genotyped-list: " + matched_id);
        }
        if (!mapped_pedigree_rows.insert(it->second).second) {
            throw std::runtime_error("Multiple genotype samples map to pedigree ID: " + matched_id);
        }
        genotyped_indices_.push_back(it->second);
        genotyped_ids_.push_back(matched_id);
    }
    if (!requested_genotyped_ids_.empty() && requested_genotyped_ids_.size() != genotyped_ids_.size()) {
        throw std::runtime_error("--genotyped-list must match the filtered genotype sample set exactly");
    }
    if ((int)genotyped_indices_.size() != n_genotyped_) {
        throw std::runtime_error("Genotype-to-pedigree mapping size does not match G matrix dimension");
    }
    std::cout << "G matrix: " << n_genotyped_ << "x" << n_genotyped_ << "\n";

    if (options_.blend_g > 0.0) {
        double alpha = options_.blend_g;
        G_ = (1.0 - alpha) * G_ + alpha * Eigen::MatrixXd::Identity(n_genotyped_, n_genotyped_);
        std::cout << "G blended with identity (alpha=" << alpha << ")\n";
    }
}

void SingleStepModel::computeA22() {
    if (n_genotyped_ <= 0) {
        has_a22_correction_ = false;
        A22_inv_.resize(0, 0);
        return;
    }
    if ((int)genotyped_indices_.size() != n_genotyped_) {
        throw std::runtime_error("ssGWAS A22 correction requires genotyped ID count to match GRM dimension.");
    }

    Eigen::MatrixXd A_dense = invertSpdWithRidge(Eigen::MatrixXd(A_inv_), "pedigree A-inverse");
    Eigen::MatrixXd A22(n_genotyped_, n_genotyped_);
    for (int i = 0; i < n_genotyped_; ++i) {
        int pi = genotyped_indices_[i];
        if (pi < 0 || pi >= A_dense.rows()) {
            throw std::runtime_error("Genotyped individual index is outside the pedigree matrix.");
        }
        for (int j = 0; j < n_genotyped_; ++j) {
            int pj = genotyped_indices_[j];
            A22(i, j) = A_dense(pi, pj);
        }
    }

    A22_inv_ = invertSpdWithRidge(A22, "pedigree A22");
    has_a22_correction_ = true;
    std::cout << "A22 correction enabled for " << n_genotyped_
              << " genotyped individuals (tau=" << options_.h_tau
              << ", omega=" << options_.h_omega << ").\n";
}

void SingleStepModel::computeHInverse() {
    H_inv_ = A_inv_;

    if (n_genotyped_ > 0) {
        Eigen::MatrixXd Ginv = invertSpdWithRidge(G_, "genomic relationship G");
        if (!has_a22_correction_ || A22_inv_.rows() != n_genotyped_) {
            throw std::runtime_error("ssGWAS H-inverse construction requires A22 correction.");
        }
        for (int gi = 0; gi < n_genotyped_; ++gi) {
            int idx_i = genotyped_indices_[gi];
            for (int gj = 0; gj < n_genotyped_; ++gj) {
                int idx_j = genotyped_indices_[gj];
                H_inv_.coeffRef(idx_i, idx_j) += options_.h_tau * Ginv(gi, gj)
                                               - options_.h_omega * A22_inv_(gi, gj);
            }
        }
    }
    H_inv_ = 0.5 * (H_inv_ + Eigen::SparseMatrix<double>(H_inv_.transpose()));
    Eigen::MatrixXd H_dense = Eigen::MatrixXd(H_inv_);
    h_inverse_symmetry_error_ = (H_dense - H_dense.transpose()).cwiseAbs().maxCoeff();

    std::cout << "Standard ssGWAS H⁻¹ constructed: A⁻¹ + [0, tau*G⁻¹ - omega*A22⁻¹]. "
              << "Size: " << H_inv_.rows() << "x" << H_inv_.cols()
              << ", symmetry max error=" << h_inverse_symmetry_error_ << "\n";
}

Eigen::MatrixXd SingleStepModel::buildHSubForG(const Eigen::MatrixXd& G) const {
    if (G.rows() != n_genotyped_ || G.cols() != n_genotyped_) {
        throw std::runtime_error("H construction requires G dimension to match genotyped individuals");
    }
    Eigen::SparseMatrix<double> Hinv = A_inv_;
    Eigen::MatrixXd Ginv = invertSpdWithRidge(G, "genomic relationship G");
    for (int gi = 0; gi < n_genotyped_; ++gi) {
        int idx_i = genotyped_indices_[gi];
        for (int gj = 0; gj < n_genotyped_; ++gj) {
            int idx_j = genotyped_indices_[gj];
            Hinv.coeffRef(idx_i, idx_j) += options_.h_tau * Ginv(gi, gj)
                                          - options_.h_omega * A22_inv_(gi, gj);
        }
    }
    Hinv = 0.5 * (Hinv + Eigen::SparseMatrix<double>(Hinv.transpose()));
    Eigen::MatrixXd H = invertSpdWithRidge(Eigen::MatrixXd(Hinv), "ssGWAS H-inverse");
    Eigen::MatrixXd Hsub(n_genotyped_, n_genotyped_);
    for (int i = 0; i < n_genotyped_; ++i) {
        for (int j = 0; j < n_genotyped_; ++j) {
            Hsub(i, j) = H(genotyped_indices_[i], genotyped_indices_[j]);
        }
    }
    return 0.5 * (Hsub + Hsub.transpose());
}

void SingleStepModel::runVCE() {
    std::cout << "Running ssGWAS VCE with Eigen REML on the genotyped H subset...\n";

    int ng = n_genotyped_;
    if (ng <= 0) throw std::runtime_error("ssGWAS VCE requires at least one genotyped individual.");
    if (ng != keep_.n_keep) {
        throw std::runtime_error("Current ssGWAS VCE requires the genotyped-list order to match the filtered genotype sample set.");
    }

    H_sub_ = buildHSubForG(G_);

    Eigen::VectorXd y_sub = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X_, keep_);

    EigenVce::Options vce_opts;
    vce_opts.use_covariate_projection = true;
    vce_opts.max_iter = options_.reml_maxit;
    vce_opts.tol = options_.reml_tol;
    EigenVce vce(vce_opts);
    vce.prepare(y_sub, X_sub, H_sub_);
    bool ok = vce.runNullModel();
    if (!ok) {
        throw std::runtime_error("ssGWAS Eigen REML null model did not converge.");
    }
    sigma2_g_ = vce.getSigma2g();
    sigma2_e_ = vce.getSigma2e();

    if (!std::isfinite(sigma2_g_) || !std::isfinite(sigma2_e_) ||
        sigma2_g_ < 0.0 || sigma2_e_ < 0.0) {
        throw std::runtime_error("ssGWAS Eigen REML produced invalid variance components.");
    }

    std::cout << "Final ssGWAS VCE: Vg=" << sigma2_g_ << " Ve=" << sigma2_e_
              << " h2=" << sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10) << "\n";
}

void SingleStepModel::scanGwas() {
    std::cout << "Running ssGWAS-informed GWAS scan (using H⁻¹ VCs)...\n";
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
    grm_opts.extract_set = loadSnpIdSet(options_.extract_file);
    grm_opts.exclude_set = loadSnpIdSet(options_.exclude_file);
    grm_opts.thread_num = options_.n_threads;

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
        G_sub = G_;
    }
    auto blend_g = [&](Eigen::MatrixXd G) {
        if (options_.blend_g > 0.0) {
            G = (1.0 - options_.blend_g) * G
              + options_.blend_g * Eigen::MatrixXd::Identity(G.rows(), G.cols());
        }
        return G;
    };
    if (use_loco) G_sub = blend_g(G_sub);

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
        std::cout << "Using H-based LMM scan (no LOCO)...\n";
        run_scan(H_sub_, "");
    } else {
        std::cout << "Running LOCO scan using ssGWAS VCs on "
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
                run_scan(buildHSubForG(blend_g(G_loco)), chrom);
            } else {
                run_scan(buildHSubForG(G_sub), chrom);
            }
        }
    }
    std::cout << "ssGWAS results saved to " << out_file << "\n";
}

void SingleStepModel::writeEbv() {
    if (!options_.write_ebv) return;
    std::cout << "GEBV output will be implemented in v1.1.\n";
}

void SingleStepModel::writeAghMatrices() const {
    if (!options_.write_agh) return;
    std::vector<std::string> pedigree_ids;
    pedigree_ids.reserve(pedigree_.size());
    for (const auto& record : pedigree_) pedigree_ids.push_back(record.id);

    write_matrix_bin(options_.out_prefix + ".Ainv.bin", options_.out_prefix + ".Ainv.id",
                     Eigen::MatrixXd(A_inv_), pedigree_ids);
    write_matrix_bin(options_.out_prefix + ".G.bin", options_.out_prefix + ".G.id",
                     G_, genotyped_ids_);
    write_matrix_bin(options_.out_prefix + ".Hinv.bin", options_.out_prefix + ".Hinv.id",
                     Eigen::MatrixXd(H_inv_), pedigree_ids);
    std::ofstream meta(options_.out_prefix + ".agh.meta.txt");
    meta << "Item\tValue\n"
         << "n_pedigree\t" << pedigree_.size() << "\n"
         << "n_genotyped\t" << n_genotyped_ << "\n"
         << "grm_algorithm\t" << options_.grm_alg << "\n"
         << "h_tau\t" << options_.h_tau << "\n"
         << "h_omega\t" << options_.h_omega << "\n"
         << "blend_g\t" << options_.blend_g << "\n"
         << "id_mapping_verified\t1\n"
         << "hinv_symmetry_max_abs\t" << h_inverse_symmetry_error_ << "\n";
    std::cout << "A/G/H-inverse matrices written with prefix " << options_.out_prefix << "\n";
}

void SingleStepModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + "_ssgwas_null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs << "sigma2_g\tsigma2_e\th2\ta22_correction\th_tau\th_omega\thinv_symmetry_max_abs\n";
    double h2 = sigma2_g_ / (sigma2_g_ + sigma2_e_ + 1e-10);
    ofs << sigma2_g_ << "\t" << sigma2_e_ << "\t" << h2 << "\t"
        << (has_a22_correction_ ? 1 : 0) << "\t"
        << options_.h_tau << "\t" << options_.h_omega << "\t"
        << h_inverse_symmetry_error_ << "\n";
    ofs.close();
    std::cout << "ssGWAS null model summary written to " << fname << "\n";
}

void SingleStepModel::run() {
    loadPedigree();
    if (options_.make_agh && !requested_genotyped_ids_.empty()) {
        keep_.keep_indices.clear();
        std::unordered_set<std::string> matched_ids;
        for (int sample_index = 0; sample_index < n_samples_; ++sample_index) {
            const std::string fid_iid = fids_[sample_index] + ":" + iids_[sample_index];
            const std::string& iid = iids_[sample_index];
            if (requested_genotyped_ids_.count(fid_iid)) {
                keep_.keep_indices.push_back(sample_index);
                matched_ids.insert(fid_iid);
            } else if (requested_genotyped_ids_.count(iid)) {
                keep_.keep_indices.push_back(sample_index);
                matched_ids.insert(iid);
            }
        }
        if (matched_ids.size() != requested_genotyped_ids_.size()) {
            throw std::runtime_error("--genotyped-list contains IDs absent from the genotype samples");
        }
        keep_.n_keep = static_cast<int>(keep_.keep_indices.size());
        std::cout << "Selected " << keep_.n_keep << " genotype samples from --genotyped-list.\n";
    }
    buildAInverse();
    buildGRM();
    computeA22();
    computeHInverse();
    writeAghMatrices();
    if (options_.make_agh) return;
    runVCE();
    scanGwas();
    writeEbv();
    writeNullSummary();
}

} // namespace cosmic
