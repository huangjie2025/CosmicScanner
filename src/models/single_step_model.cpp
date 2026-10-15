#include "models/lanczos_eigensolver.h"
#include "models/single_step_model.h"
#include "models/cross_genomic_operator.h"
#include "models/h_inverse_operator.h"
#include "models/snp_uncertainty_scanner.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <limits>
#include <omp.h>
#include "data_loader.h"
#include "data_filter.h"
#include "pheno_reader.h"
#include "grm_builder.h"
#include "gwas_scanner.h"
#include "lmm_engine.h"
#include "eigen_vce.h"
#include "a_matrix.h"
#include "matrix_io.h"
#include "mme_builder.h"
#include "matrix_adapter.h"
#include "implicit_mme.h"
#include "matrix_free_pcg.h"
#include "vce.h"
#include "design.h"
#include "apy_mme_operator.h"

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

std::string format_dp(double v, int prec) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(prec) << v;
    return os.str();
}

std::string format_sci(double v, int prec) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(prec) << v;
    return os.str();
}

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

    opts.method = (cfg.single_step_method == "backsolve") ? Method::Backsolve : Method::Lmm;
    opts.ssgwas_window_snps = cfg.ssgwas_window_snps;
    opts.ssgwas_tune_g = cfg.ssgwas_tune_g;
    opts.ssgwas_blend_a22 = cfg.ssgwas_blend_a22;
    opts.ssgwas_solver = cfg.ssgwas_solver;
    opts.ssgwas_vce = cfg.ssgwas_vce;
    opts.reml_vg = cfg.reml_vg;
    opts.reml_ve = cfg.reml_ve;
    opts.ssgwas_pcg_tol = cfg.ssgwas_pcg_tol;
    opts.ssgwas_pcg_maxit = cfg.ssgwas_pcg_maxit;
    opts.ssgwas_pvalue_method = cfg.ssgwas_pvalue_method;
    opts.ssgwas_apy_exact_estimator = cfg.ssgwas_apy_exact_estimator;
    opts.ssgwas_apy_exact_estimator_explicit = cfg.ssgwas_apy_exact_estimator_explicit;
    opts.allow_large_exact_ssgwas_pvalue = cfg.allow_large_exact_ssgwas_pvalue;
    opts.ssgwas_grm_mode = cfg.ssgwas_grm_mode;
    opts.ssgwas_gcc_solver = cfg.ssgwas_gcc_solver;
    opts.ssgwas_apy_inverse_mode = cfg.ssgwas_apy_inverse_mode;
    opts.ssgwas_mme_apy_mode = cfg.ssgwas_mme_apy_mode;
    opts.ssgwas_pev_apy_mode = cfg.ssgwas_pev_apy_mode;
    opts.ssgwas_hinv_mode = cfg.ssgwas_hinv_mode;
    opts.ssgwas_pev_scan_mode = cfg.ssgwas_pev_scan_mode;
    opts.ssgwas_pev_scan_mode_explicit = cfg.ssgwas_pev_scan_mode_explicit;
    opts.ssgwas_pev_scan_block_size = cfg.ssgwas_pev_scan_block_size;
    opts.ssgwas_pcg_preconditioner = cfg.ssgwas_pcg_preconditioner;
    opts.ssgwas_pcg_preconditioner_rank = cfg.ssgwas_pcg_preconditioner_rank;
    opts.ssgwas_apy_core_size = cfg.ssgwas_apy_core_size;
    opts.ssgwas_apy_core_file = cfg.ssgwas_apy_core_file;
    opts.ssgwas_apy_seed = cfg.ssgwas_apy_seed;
    opts.ssgwas_apy_core_target_ratio = cfg.ssgwas_apy_core_target_ratio;
    opts.ssgwas_apy_core_variance = cfg.ssgwas_apy_core_variance;
    opts.ssgwas_apy_core_eigensolver = cfg.ssgwas_apy_core_eigensolver;
    return opts;
}

SingleStepModel::SingleStepModel(const Options& opts) : options_(opts) {}

void SingleStepModel::prepare(const std::vector<std::string>& fids,
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

    if (options_.make_agh) {
        keep_.keep_indices.resize(n_samples_);
        std::iota(keep_.keep_indices.begin(), keep_.keep_indices.end(), 0);
        keep_.n_keep = n_samples_;
    } else if (options_.method == Method::Backsolve) {
        // Method 2 (ssGBLUP back-solving) loads phenotypes by PEDIGREE ID
        // inside runBacksolve(); the sample-index phenotype path cannot cover
        // phenotyped animals without genotypes and is intentionally bypassed.
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
    if (options_.ped_file.empty()) throw std::runtime_error("Single-step model requires --ped file.");
    pedigree_.load(options_.ped_file);

    if (!options_.genotyped_list.empty()) {
        std::ifstream gf(options_.genotyped_list);
        std::string gid;
        while (gf >> gid) {
            if (!requested_genotyped_ids_.insert(gid).second) {
                throw std::runtime_error("Duplicate ID in --genotyped-list: " + gid);
            }
            if (pedigree_.idToIndex().find(gid) == pedigree_.idToIndex().end()) {
                throw std::runtime_error("Genotyped ID is absent from pedigree: " + gid);
            }
        }
        std::cout << requested_genotyped_ids_.size() << " requested genotyped IDs loaded.\n";
    }
}

void SingleStepModel::buildAInverse() {
    A_inv_ = pedigree_.buildAInverse();
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
    const auto& ped_id_to_idx = pedigree_.idToIndex();
    std::unordered_set<int> mapped_pedigree_rows;
    for (int sample_index : keep_.keep_indices) {
        const std::string fid_iid = fids_[sample_index] + ":" + iids_[sample_index];
        const std::string& iid = iids_[sample_index];
        auto it = ped_id_to_idx.find(fid_iid);
        std::string matched_id = fid_iid;
        if (it == ped_id_to_idx.end()) {
            it = ped_id_to_idx.find(iid);
            matched_id = iid;
        }
        if (it == ped_id_to_idx.end()) {
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

// ---------------------------------------------------------------------------
// Phase D, Leite et al. 2024 Eq. 9: diagonal effective records after
// fixed-effect absorption,
//   D_animal = diag(W'(I - X(X'X)^-1 X')W)   ->   D_i = sum_{records r of i}
//                                                    (1 - x_r'(X'X)^-1 x_r),
// computed exactly from the record-level design (no off-diagonal dropping
// beyond the paper's own diagonal statement).  Only genotyped animals appear
// (D_c / D_n split for the Eq. 8 core-PEV approximation).
// ---------------------------------------------------------------------------
void SingleStepModel::computeEffectiveRecordDiagonals(const std::vector<GenRecord>& recs,
                                                     const FixedDesignG& fd) {
    const int p = fd.p;
    d_genotyped_.assign(n_genotyped_, 0.0);
    if (p <= 0 || n_genotyped_ <= 0) return;

    Eigen::MatrixXd XtX = Eigen::MatrixXd::Zero(p, p);
    for (const auto& row : fd.rows) {
        for (const auto& a : row) {
            for (const auto& b : row) {
                XtX(a.first, b.first) += a.second * b.second;
            }
        }
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0) {
        throw std::runtime_error(
            "computeEffectiveRecordDiagonals: fixed-effect matrix X'X is not SPD");
    }
    std::unordered_map<int, int> ped_row_to_geno;
    ped_row_to_geno.reserve(n_genotyped_);
    for (int g = 0; g < n_genotyped_; ++g) ped_row_to_geno[genotyped_indices_[g]] = g;
    for (size_t r = 0; r < recs.size() && r < fd.rows.size(); ++r) {
        auto it = ped_row_to_geno.find(recs[r].aid - 1);
        if (it == ped_row_to_geno.end()) continue;
        Eigen::VectorXd xr = Eigen::VectorXd::Zero(p);
        for (const auto& pr : fd.rows[r]) xr(pr.first) = pr.second;
        d_genotyped_[it->second] += 1.0 - xr.dot(ldlt.solve(xr));
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

    Eigen::MatrixXd A22 = pedigree_.buildA22(genotyped_indices_);
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
    pedigree_ids.reserve(pedigree_.nPedigree());
    for (const auto& record : pedigree_.records()) pedigree_ids.push_back(record.id);

    write_matrix_bin(options_.out_prefix + ".Ainv.bin", options_.out_prefix + ".Ainv.id",
                     Eigen::MatrixXd(A_inv_), pedigree_ids);
    write_matrix_bin(options_.out_prefix + ".G.bin", options_.out_prefix + ".G.id",
                     G_, genotyped_ids_);
    write_matrix_bin(options_.out_prefix + ".Hinv.bin", options_.out_prefix + ".Hinv.id",
                     Eigen::MatrixXd(H_inv_), pedigree_ids);
    std::ofstream meta(options_.out_prefix + ".agh.meta.txt");
    meta << "Item\tValue\n"
         << "n_pedigree\t" << pedigree_.nPedigree() << "\n"
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
    if ((options_.ssgwas_apy_exact_estimator_explicit ||
         options_.ssgwas_apy_exact_estimator == "core-v2") && options_.method != Method::Backsolve) {
        throw std::runtime_error("--ssgwas-apy-exact-estimator requires --single-step-method backsolve");
    }
    auto t0_ped = std::chrono::steady_clock::now();
    loadPedigree();
    if (options_.method == Method::Backsolve) {
        t_parse_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_ped).count();
        runBacksolve();
        return;
    }
    if (options_.ssgwas_pvalue_method != "off") {
        throw std::runtime_error(
            "--ssgwas-pvalue-method exact requires --single-step-method backsolve (Method 2); "
            "Method 1 (single-step fixed-SNP LMM association) has its own Wald/Score P and "
            "must not reuse the back-solved SNP-effect P.");
    }
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

// ============================================================================
// Method 2: ssGBLUP genetic evaluation -> GEBV -> SNP back-solving (Phase 6)
// ============================================================================
//
// Statistical convention (docs/SSGWAS_METHOD_AUDIT.md §6/§10):
//   G      = M M' / S,  M = X - mean (dosage),  S = sum_j 2 p_j (1 - p_j)
//   backsolve:  a_hat_g ~= Z a   =>   a = Z' (Z Z')^-1 a_hat_g
//               with (Z Z') = G * S  this equals  a = Z' G^-1 a_hat_g / S
//   implemented as  q = Ginv * a_hat_g ;  zt = Z' q (multiply_Zt_v) ;
//                  effect_allele = zt / (std * S)  (per allele-dosage scale,
//                  matching the M = x - mean centering of the VanRaden G)
//                  effect_std    = effect_allele * sd  (per standardized-genotype
//                  SD; reconstruction-consistent:  Z*effect_std = M*effect_allele
//                  = a_hat_g on the genotyped block)
//   SNP variance   = 2 p_j (1 - p_j) * effect_allele_j^2
//   PCT_SNP        = SNP variance_j / sum_k SNP variance_k
// No SNP P-values / SE / Wald are computed for Method 2 (Phase C).

void SingleStepModel::runBacksolve() {
    std::cout << "=== ssGBLUP back-solving (--single-step-method backsolve) ===\n";
    point_multiply_zt_v_calls_ = 0;
    if (options_.ssgwas_apy_exact_estimator != "core-v2" &&
        options_.ssgwas_apy_exact_estimator != "legacy-hybrid-v1") {
        throw std::runtime_error("--ssgwas-apy-exact-estimator must be core-v2|legacy-hybrid-v1");
    }
    const bool core_v2 = options_.ssgwas_apy_exact_estimator == "core-v2";
    if ((core_v2 || options_.ssgwas_apy_exact_estimator_explicit) &&
        options_.ssgwas_pvalue_method != "apy-exact") {
        throw std::runtime_error("--ssgwas-apy-exact-estimator requires --ssgwas-pvalue-method apy-exact");
    }
    if (options_.ssgwas_pvalue_method == "apy-exact") {
        if (core_v2) {
            // R2A validated scope: actual Hgg equals the final genomic prior,
            // and the same direct Gcc inverse drives point and variance.
            if (options_.h_tau != 1.0 || options_.h_omega != 1.0 ||
                options_.ssgwas_tune_g || options_.ssgwas_blend_a22 ||
                options_.blend_g != 0.05 || options_.ssgwas_gcc_solver != "dense" ||
                options_.ssgwas_apy_core_size == -1 || !options_.ssgwas_vce.empty() ||
                !std::isfinite(options_.reml_vg) || !std::isfinite(options_.reml_ve) ||
                !(options_.reml_vg > 0.0) || !(options_.reml_ve > 0.0)) {
                throw std::runtime_error(
                    "APY-exact core-v2 validated scope requires fixed positive Va/Ve, "
                    "tau=omega=1, tuning disabled, identity blend 0.05, dense Gcc solver "
                    "and fixed core selection; unsupported settings are not silently accepted.");
            }
        } else {
            std::cerr << "WARNING: APY-exact legacy-hybrid-v1 preserves raw/full SNP points "
                         "with core-only uncertainty. This pairing is generally statistically "
                         "inconsistent and is not recommended for calibrated inference. "
                         "Use --ssgwas-apy-exact-estimator core-v2 in its validated scope.\n";
        }
    }
    const bool is_apy = options_.ssgwas_pvalue_method == "apy-exact" ||
                        options_.ssgwas_pvalue_method == "apy-approx";
    if (options_.ssgwas_pvalue_method != "off" && options_.ssgwas_pvalue_method != "exact" &&
        options_.ssgwas_pvalue_method != "apy-exact" && options_.ssgwas_pvalue_method != "apy-approx") {
        throw std::runtime_error("Invalid --ssgwas-pvalue-method value '" + options_.ssgwas_pvalue_method +
                                 "' (expected off|exact|apy-exact|apy-approx)");
    }
    if (options_.ssgwas_pvalue_method == "apy-approx") {
        // Phase D: permitted now; the PCG non-convergence gate refuses
        // formal P-values downstream (single_step_model.cpp solveGeneticEvaluation).
    }
    if ((options_.ssgwas_pvalue_method == "exact" || options_.ssgwas_pvalue_method == "apy-exact") &&
        options_.ssgwas_solver != "ldlt") {
        throw std::runtime_error(
            "--ssgwas-pvalue-method exact/apy-exact currently requires the LDLT MME factorization "
            "because selected elements of the inverse coefficient matrix are needed; "
            "use --ssgwas-solver ldlt (no silent fallback is performed).");
    }
    if (is_apy) {
        if (options_.ssgwas_apy_core_size != -1 && options_.ssgwas_apy_core_size < 1
            && options_.ssgwas_apy_core_file.empty()) {
            throw std::runtime_error(
                "--ssgwas-pvalue-method apy-* requires a core partition: pass "
                "--ssgwas-apy-core-size N (1..n_genotyped; N = n_genotyped reproduces the "
                "exact path), 'auto' for automatic C98/C99 sizing, or "
                "--ssgwas-apy-core-file FILE.");
        }
    }
    if (options_.ssgwas_mme_apy_mode == "operator") {
        if (!is_apy) {
            throw std::runtime_error(
                "--ssgwas-mme-apy-mode operator requires --ssgwas-pvalue-method apy-exact|apy-approx");
        }
        if (options_.ssgwas_solver != "pcg") {
            throw std::runtime_error(
                "--ssgwas-mme-apy-mode operator requires --ssgwas-solver pcg; explicit LDLT MME assembly requires the dense APY backend");
        }
        if (!options_.ssgwas_vce.empty()) {
            throw std::runtime_error(
                "--ssgwas-mme-apy-mode operator currently requires fixed --reml-vg/--reml-ve; operator-backed VCE is not enabled");
        }
    }
    if (options_.ssgwas_hinv_mode == "operator") {
        if (options_.ssgwas_solver != "pcg") {
            throw std::runtime_error(
                "--ssgwas-hinv-mode operator requires --ssgwas-solver pcg; dense/reference LDLT remains available with --ssgwas-hinv-mode dense");
        }
        if (!options_.ssgwas_vce.empty()) {
            throw std::runtime_error(
                "--ssgwas-hinv-mode operator currently requires fixed --reml-vg/--reml-ve; operator-backed VCE is not enabled");
        }
        if (options_.ssgwas_tune_g || options_.ssgwas_blend_a22) {
            std::cout << "WARNING: Hinv operator is active, but tune-G/A22 blending still uses the dense A22 reference path; A22_BLEND_MODE=dense.\n";
        }
    }
    if (options_.ssgwas_pev_apy_mode == "operator" &&
        options_.ssgwas_pvalue_method != "apy-approx") {
        throw std::runtime_error(
            "--ssgwas-pev-apy-mode operator requires --ssgwas-pvalue-method apy-approx");
    }
    if (options_.ssgwas_pev_scan_mode != "scalar" &&
        options_.ssgwas_pev_scan_mode != "batched") {
        throw std::runtime_error(
            "--ssgwas-pev-scan-mode must be scalar|batched");
    }
    if (options_.ssgwas_pev_scan_block_size <= 0) {
        throw std::runtime_error("--ssgwas-pev-scan-block-size must be > 0");
    }
    if (is_bgen_ || is_pgen_) {
        throw std::runtime_error("--single-step-method backsolve requires PLINK .bed/.bim/.fam "
                                 "input in Phase 6 v1 (BGEN/PGEN GEBV + back-solving pending)");
    }
    auto now_t = std::chrono::steady_clock::now;
    auto secs = [](const std::chrono::steady_clock::time_point& t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    loadBacksolvePhenotypes();

    auto t0 = now_t();
    buildAInverse();
    t_ainv_ = secs(t0);

    t0 = now_t();
    setupGenotypedSamples();
    resolveStreamedGrmMode();
    if (grm_streamed_) {
        // Phase E2: block-streamed G (pass 1: allele stats + S + tune scalars;
        // the Gcc/Gcn/diag(Gnn) block pass runs after the APY core selection).
        buildBacksolveGrmStreamed();
    } else {
        buildBacksolveGrm();
    }
    t_g_ = secs(t0);

    t0 = now_t();
    if (is_apy) {
        // Phase C2: no dense n_g^2 G^-1 is built.  The APY block inverse
        // (Leite et al. 2024 Eqs. 5-6) replaces Ginv_bs_ in the MME and the
        // back-solve operator; the uncertainty block is core-sized (Eq. 7).
        // Automatic core-size determination (C98/C99) if requested.
        if (options_.ssgwas_apy_core_size == -1) {
            if (n_genotyped_ < 10) {
                throw std::runtime_error(
                    "Automatic APY core sizing requires at least 10 genotyped animals");
            }
            const double target = options_.ssgwas_apy_core_variance;
            if (target <= 0.0 || target > 1.0) {
                throw std::runtime_error("--ssgwas-apy-core-variance must be in (0, 1]");
            }
            if (target < 0.90) {
                std::cout << "WARNING: --ssgwas-apy-core-variance " << target
                          << " is below 0.90; this may yield an under-sized core.\n";
            }

            // Eigensolver selection: exact (dense) for small n_g, Lanczos for large
            const std::string& solver_pref = options_.ssgwas_apy_core_eigensolver;
            const int LANCZOS_THRESHOLD = 2000;
            bool use_lanczos = grm_streamed_ || (solver_pref == "lanczos") ||
                               (solver_pref == "auto" && n_genotyped_ > LANCZOS_THRESHOLD);

            int K_auto = n_genotyped_;
            int K98 = n_genotyped_;
            int K99 = n_genotyped_;
            double trace_g = 0.0;
            int eig_iterations = 0;

            if (use_lanczos) {
                // Matrix-free Lanczos eigendecomposition
                apy_core_eigensolver_used_ = "lanczos";
                std::cout << "Using matrix-free Lanczos eigensolver for auto core sizing (n_g="
                          << n_genotyped_ << ")\n";

                // Matrix-free operator reproduces the production tuned+blended G:
                //   (1-alpha)*(b*G_raw + a*J) + alpha*T, with T = I (default) or A22
                //   (--ssgwas-blend-a22). tune_a_/tune_b_ are buildBacksolveGrm
                //   results (0/1 defaults when tuneG is off).
                std::unique_ptr<Eigen::MatrixXd> a22_blend;
                if (options_.ssgwas_blend_a22) {
                    a22_blend = std::make_unique<Eigen::MatrixXd>(
                        pedigree_.buildA22(genotyped_indices_));
                }
                MatrixFreeGrmOperator G_op(
                    *plink_reader_, geno_full_idx_, S_2pq_,
                    tune_a_, tune_b_,
                    options_.blend_g > 0.0 ? options_.blend_g : 0.0,
                    options_.ssgwas_blend_a22 ? a22_blend.get() : nullptr,
                    options_.n_threads > 0 ? options_.n_threads : 1);

                LanczosEigenSolver lanczos;
                LanczosEigenSolver::Options lanczos_opts;
                lanczos_opts.max_iterations = std::min(80, n_genotyped_);
                lanczos_opts.tolerance = 1e-6;
                lanczos_opts.seed = static_cast<unsigned long>(options_.ssgwas_apy_seed);
                lanczos_opts.n_threads = options_.n_threads > 0 ? options_.n_threads : 1;

                LanczosResult lres = lanczos.compute(G_op, target, lanczos_opts);

                if (lres.status != "OK") {
                    throw std::runtime_error("Lanczos eigendecomposition failed: " + lres.status);
                }

                K_auto = lres.K_target;
                K98 = lres.K98;
                K99 = lres.K99;
                trace_g = lres.trace_g;
                eig_iterations = lres.iterations;

                std::cout << "Lanczos: " << lres.iterations << " iterations, "
                          << lres.convergedritz << " Ritz values, trace(G)="
                          << std::fixed << std::setprecision(6) << trace_g << "\n";

            } else {
                // Dense exact eigendecomposition (existing code path)
                apy_core_eigensolver_used_ = "dense";
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs_g(G_bs_);
                if (eigs_g.info() != Eigen::Success) {
                    throw std::runtime_error("Eigen decomposition of G_bs_ failed for automatic core sizing");
                }
                Eigen::VectorXd ev = eigs_g.eigenvalues();
                std::sort(ev.data(), ev.data() + ev.size(), std::greater<double>());
                trace_g = ev.sum();
                if (trace_g <= 0.0) {
                    throw std::runtime_error("G_bs_ has non-positive total eigenvalue sum");
                }

                double cum = 0.0;
                for (int k = 0; k < n_genotyped_; ++k) {
                    cum += std::max(0.0, ev(k));
                    if (K_auto == n_genotyped_ && cum >= target * trace_g) K_auto = k + 1;
                    if (K98 == n_genotyped_ && cum >= 0.98 * trace_g) K98 = k + 1;
                    if (K99 == n_genotyped_ && cum >= 0.99 * trace_g) K99 = k + 1;
                }
                eig_iterations = n_genotyped_;
            }

            if (K_auto < 1) K_auto = 1;
            if (K_auto > n_genotyped_) K_auto = n_genotyped_;
            options_.ssgwas_apy_core_size = K_auto;

            std::cout << "Automatic APY core selection:\n"
                      << "  mode                = automatic\n"
                      << "  eigensolver         = " << apy_core_eigensolver_used_ << "\n"
                      << "  variance target     = " << std::fixed << std::setprecision(6) << target << "\n"
                      << "  genomic dimension   = " << K_auto << "\n"
                      << "  genotyped animals   = " << n_genotyped_ << "\n"
                      << "  core fraction       = " << std::fixed << std::setprecision(4)
                      << static_cast<double>(K_auto) / n_genotyped_ << "\n"
                      << "  C98                 = " << K98 << "\n"
                      << "  C99                 = " << K99 << "\n"
                      << "  trace(G)            = " << std::fixed << std::setprecision(6) << trace_g << "\n"
                      << "  eig iterations      = " << eig_iterations << "\n"
                      << "  selection strategy  = seeded-random\n"
                      << "  random seed         = " << options_.ssgwas_apy_seed << "\n";
            // Store for meta file
            apy_core_auto_ = true;
            apy_core_auto_K98_ = K98;
            apy_core_auto_K99_ = K99;
            apy_core_auto_variance_ = target;
            apy_core_trace_g_ = trace_g;
            apy_core_eig_iterations_ = eig_iterations;
        }
        buildApyCoreSelection();
        std::cout << "APY core partition: n_core=" << apy_n_core_
                  << " n_noncore=" << apy_n_noncore_
                  << " selection=" << apy_core_selection_mode_
                  << " seed=" << (apy_core_selection_mode_ == "RANDOM_SEEDED"
                                      ? std::to_string(options_.ssgwas_apy_seed) : "NA") << "\n";
        if (grm_streamed_) {
            buildBacksolveGrmStreamBlocks();
        } else {
            if (options_.ssgwas_mme_apy_mode == "dense") {
                apy_ginv_mme_ = std::make_unique<ApyGinv>();
                apy_ginv_mme_->build(G_bs_, core_indices_, "APY G (MME; final Gw)");
            }
            if (options_.ssgwas_blend_a22) {
                if (options_.ssgwas_apy_inverse_mode == "dense") {
                    apy_ginv_eff_ = std::make_unique<ApyGinv>();
                    apy_ginv_eff_->build(G_bs_tuned_, core_indices_, "APY G* (effect; tuned-unblended)");
                }
                effect_ginv_name_ = "preblend-apy";
            } else {
                if (options_.ssgwas_apy_inverse_mode == "dense") {
                    apy_ginv_eff_ = std::make_unique<ApyGinv>();
                    apy_ginv_eff_->build(G_bs_, core_indices_, "APY G (effect)");
                }
                effect_ginv_name_ = "apy";
            }
        }
        if (options_.ssgwas_apy_inverse_mode == "operator") {
            GcnBlockIteratorFactory factory;
            const Eigen::VectorXd* gnn_diag = nullptr;
            Eigen::VectorXd dense_gnn_diag;
            if (grm_streamed_) {
                const GcnBlockView view = options_.ssgwas_blend_a22
                                              ? GcnBlockView::Tuned
                                              : GcnBlockView::Final;
                factory = [this, view]() {
                    return grm_stream_builder_->makeGcnIterator(-1, view);
                };
                gnn_diag = options_.ssgwas_blend_a22
                               ? &gnn_diag_stream_tuned_
                               : &gnn_diag_stream_final_;
            } else {
                const Eigen::MatrixXd& effect_g = options_.ssgwas_blend_a22 ? G_bs_tuned_ : G_bs_;
                factory = makeDenseGcnBlockIteratorFactory(
                    effect_g, core_indices_, noncore_indices_);
                dense_gnn_diag.resize(noncore_indices_.size());
                for (int j = 0; j < static_cast<int>(noncore_indices_.size()); ++j) {
                    dense_gnn_diag(j) = effect_g(noncore_indices_[j], noncore_indices_[j]);
                }
                gnn_diag = &dense_gnn_diag;
            }
            Eigen::MatrixXd effect_gcc;
            if (grm_streamed_) {
                effect_gcc = options_.ssgwas_blend_a22 ? gcc_stream_tuned_ : gcc_stream_final_;
            } else {
                const Eigen::MatrixXd& effect_g = options_.ssgwas_blend_a22 ? G_bs_tuned_ : G_bs_;
                effect_gcc.resize(core_indices_.size(), core_indices_.size());
                for (int i = 0; i < effect_gcc.rows(); ++i)
                    for (int j = 0; j < effect_gcc.cols(); ++j)
                        effect_gcc(i, j) = effect_g(core_indices_[i], core_indices_[j]);
            }
            apy_inverse_eff_operator_ = std::make_unique<MatrixFreeApyInverseOperator>(
                makeDenseSolver(effect_gcc), std::move(factory), *gnn_diag,
                n_genotyped_, core_indices_, "APY effect operator",
                options_.ssgwas_blend_a22 ? gcn_cross_tuned_.get()
                                          : gcn_cross_final_.get());
        } else {
            apy_inverse_eff_operator_ =
                std::make_unique<DenseApyInverseOperator>(*apy_ginv_eff_);
        }
        if (options_.ssgwas_mme_apy_mode == "operator") {
            GcnBlockIteratorFactory factory;
            const Eigen::VectorXd* gnn_diag = nullptr;
            Eigen::VectorXd dense_gnn_diag;
            if (grm_streamed_) {
                factory = [this]() {
                    return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Final);
                };
                gnn_diag = &gnn_diag_stream_final_;
            } else {
                factory = makeDenseGcnBlockIteratorFactory(
                    G_bs_, core_indices_, noncore_indices_);
                dense_gnn_diag.resize(noncore_indices_.size());
                for (int j = 0; j < static_cast<int>(noncore_indices_.size()); ++j) {
                    dense_gnn_diag(j) = G_bs_(noncore_indices_[j], noncore_indices_[j]);
                }
                gnn_diag = &dense_gnn_diag;
            }
            Eigen::MatrixXd mme_gcc;
            if (grm_streamed_) {
                mme_gcc = gcc_stream_final_;
            } else {
                mme_gcc.resize(core_indices_.size(), core_indices_.size());
                for (int i = 0; i < mme_gcc.rows(); ++i)
                    for (int j = 0; j < mme_gcc.cols(); ++j)
                        mme_gcc(i, j) = G_bs_(core_indices_[i], core_indices_[j]);
            }
            apy_inverse_mme_operator_ = std::make_unique<MatrixFreeApyInverseOperator>(
                makeDenseSolver(mme_gcc), std::move(factory), *gnn_diag,
                n_genotyped_, core_indices_, "APY MME operator",
                gcn_cross_final_.get());
        }
        if (apy_ginv_mme_) {
            apy_gcc_min_diag_d_ = format_sci(apy_ginv_mme_->info().gcc_min_diag_d, 6);
            apy_gcc_max_diag_d_ = format_sci(apy_ginv_mme_->info().gcc_max_diag_d, 6);
            apy_mnn_min_ = apy_ginv_mme_->info().mnn_min;
            apy_mnn_max_ = apy_ginv_mme_->info().mnn_max;
        } else {
            const auto* op = dynamic_cast<const MatrixFreeApyInverseOperator*>(
                apy_inverse_mme_operator_.get());
            apy_gcc_min_diag_d_ = format_sci(op->gccSolverPtr()->minDiagD(), 6);
            apy_gcc_max_diag_d_ = format_sci(op->gccSolverPtr()->maxDiagD(), 6);
            apy_mnn_min_ = op->mnn().size() ? op->mnn().minCoeff() : 0.0;
            apy_mnn_max_ = op->mnn().size() ? op->mnn().maxCoeff() : 0.0;
        }
    } else {
        Ginv_bs_ = invertSpdWithRidge(G_bs_, "backsolve genomic relationship G");
        if (options_.ssgwas_blend_a22) {
            Ginv_bs_tuned_ = invertSpdWithRidge(G_bs_tuned_, "backsolve tuned genomic relationship G*");
        }
    }
    // Phase D5: eigen-spectrum K98 of the final blended G (cumulative 98% of
    // eigenvalue variance); diagnostic only - recorded in meta, the APY core
    // target decision belongs to the benchmark harness.
    k98_g_ = -1;
    if (n_genotyped_ >= 10) {
        if (grm_streamed_) {
            // Block-streamed path has no dense G_bs_: C98 of the final blended
            // G from the matrix-free Lanczos operator (E1).  When automatic
            // core sizing already ran a Lanczos on the same operator its C98
            // is reused so no duplicate factorization occurs.
            if (apy_core_auto_) {
                k98_g_ = apy_core_auto_K98_;
            } else {
                k98_g_ = runStreamedK98();
            }
        } else {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs_g(G_bs_);
            if (eigs_g.info() == Eigen::Success) {
                Eigen::VectorXd ev = eigs_g.eigenvalues();
                std::sort(ev.data(), ev.data() + ev.size(), std::greater<double>());
                const double total = ev.sum();
                if (total > 0.0) {
                    double cum = 0.0;
                    for (int k = 0; k < n_genotyped_; ++k) {
                        cum += std::max(0.0, ev(k));
                        if (cum >= k98_cumfraction_ * total) { k98_g_ = k + 1; break; }
                    }
                    if (k98_g_ < 0) k98_g_ = n_genotyped_;
                }
            }
        }
    }
    t_ginv_ = secs(t0);

    t0 = now_t();
    if (options_.ssgwas_hinv_mode == "operator") {
        a22_inverse_operator_ = std::make_unique<A22InverseOperator>(A_inv_, genotyped_indices_);
        A22_inv_.resize(0, 0);
        has_a22_correction_ = true;
        std::cout << "A22 correction operator enabled for " << n_genotyped_
                  << " genotyped individuals (untyped=" << a22_inverse_operator_->nUntyped()
                  << ", P nnz=" << a22_inverse_operator_->pNonZeros()
                  << ", Q nnz=" << a22_inverse_operator_->qNonZeros()
                  << ", R nnz=" << a22_inverse_operator_->rNonZeros()
                  << ", factor=" << a22_inverse_operator_->factorTimeSeconds() << "s).\n";
    } else {
        computeA22();  // dense reference A22_inv_ for the genotyped block of Hinv
    }
    t_hnv_ = secs(t0);

    t0 = now_t();
    solveGeneticEvaluation();
    t_mme_ = secs(t0);

    t0 = now_t();
    backsolveSnpEffects();
    t_bs_ = secs(t0);

    t0 = now_t();
    if (options_.ssgwas_pvalue_method == "exact") {
        computeSnpPvaluesExact();
        writePvalueFile();
        t_pev_out_ = secs(t0) - (t_pev_inverse_ + t_pev_vhat_ + t_pev_scan_);
    } else if (options_.ssgwas_pvalue_method == "apy-exact") {
        computeSnpPvaluesApyExact();
        writePvalueFile();
        t_pev_out_ = secs(t0) - (t_pev_inverse_ + t_pev_vhat_ + t_pev_scan_);
    } else if (options_.ssgwas_pvalue_method == "apy-approx") {
        if (options_.ssgwas_solver == "pcg" && !solver_converged_) {
            throw std::runtime_error(
                "apy-approx SNP p-values refused: MME_SOLVER_CONVERGED=false "
                "(PCG did not reach --ssgwas-pcg-tol within --ssgwas-pcg-maxit); "
                "no formal P-values are produced from a non-converged MME.");
        }
        if (options_.ssgwas_pev_apy_mode == "dense") {
            materializeApyMmeReference();
        }
        computeSnpPvaluesApyApprox();
        writePvalueFile();
        t_pev_out_ = secs(t0) - (t_pev_inverse_ + t_pev_vhat_ + t_pev_scan_);
    }

    t0 = now_t();
    writeEbvFile();
    writeSnpFile();
    t_snp_ = secs(t0);

    t0 = now_t();
    writeWindowFile();
    t_win_ = secs(t0);

    writeMetaFile();
}

// ---------------------------------------------------------------------------
// Phenotype loading by PEDIGREE ID (Method 2): a phenotyped animal is included
// whether or not it is genotyped; a genotyped animal without phenotype is still
// in the MME through the random-effect equations (Hinv rows).
// ---------------------------------------------------------------------------
void SingleStepModel::loadBacksolvePhenotypes() {
    pheno_ped_rows_.clear();
    pheno_y_.clear();
    PhenoReader pr;
    pr.load(options_.pheno_file, options_.mpheno);
    const auto pheno_records = pr.loadAllRecords(options_.pheno_file, options_.mpheno);
    const auto& ped_id_to_idx = pedigree_.idToIndex();
    std::unordered_set<int> seen_rows;
    int unmatched = 0, dup = 0;
    for (const auto& r : pheno_records) {
        auto it = ped_id_to_idx.find(r.fid + ":" + r.iid);
        if (it == ped_id_to_idx.end()) it = ped_id_to_idx.find(r.iid);
        if (it == ped_id_to_idx.end()) {
            ++unmatched;
            continue;
        }
        if (!seen_rows.insert(it->second).second) {
            ++dup;
            continue;
        }
        pheno_ped_rows_.push_back(it->second);
        pheno_y_.push_back(r.value);
    }
    if (unmatched > 0) {
        std::cout << "NOTE: " << unmatched
                  << " phenotype record(s) had no pedigree ID and were ignored.\n";
    }
    if (dup > 0) {
        std::cout << "NOTE: " << dup
                  << " duplicate phenotype record(s) per animal ignored (first kept).\n";
    }
    if (pheno_ped_rows_.empty()) {
        throw std::runtime_error(
            "No phenotype records matched the pedigree (check FID/IID spelling and --mpheno).");
    }
    std::cout << "Phenotyped animals (by pedigree ID): " << pheno_ped_rows_.size() << "\n";
}

// ---------------------------------------------------------------------------
// Genotyped sample selection + G construction (independent of the Method-1
// keep set: genotype animals may lack phenotype, phenotype animals may lack
// genotype, ancestors may have neither).
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Genotyped sample selection + pedigree mapping (steps 1-2 of the GRM build,
// shared by the dense and the block-streamed (Phase E2) paths).
// ---------------------------------------------------------------------------
void SingleStepModel::setupGenotypedSamples() {
    std::cout << "Building back-solving GRM (G = M M' / S, VanRaden 2p(1-p))...\n";

    // 1. Genotyped sample selection (all samples unless --genotyped-list given)
    geno_full_idx_.clear();
    if (!requested_genotyped_ids_.empty()) {
        for (int i = 0; i < n_samples_; ++i) {
            const std::string fid_iid = fids_[i] + ":" + iids_[i];
            if (requested_genotyped_ids_.count(fid_iid) || requested_genotyped_ids_.count(iids_[i])) {
                geno_full_idx_.push_back(i);
            }
        }
        if ((int)geno_full_idx_.size() != (int)requested_genotyped_ids_.size()) {
            throw std::runtime_error("--genotyped-list contains IDs absent from the genotype samples");
        }
    } else {
        geno_full_idx_.resize(n_samples_);
        std::iota(geno_full_idx_.begin(), geno_full_idx_.end(), 0);
    }
    n_genotyped_ = (int)geno_full_idx_.size();

    // 2. Genotype sample -> pedigree row mapping (explicit ID mapping)
    const auto& ped_id_to_idx = pedigree_.idToIndex();
    const auto& ped_records = pedigree_.records();
    genotyped_indices_.clear();
    genotyped_ids_.clear();
    std::unordered_set<int> mapped_pedigree_rows;
    for (int i = 0; i < n_genotyped_; ++i) {
        int sample_idx = geno_full_idx_[i];
        const std::string fid_iid = fids_[sample_idx] + ":" + iids_[sample_idx];
        const std::string& iid = iids_[sample_idx];
        auto it = ped_id_to_idx.find(fid_iid);
        if (it == ped_id_to_idx.end()) it = ped_id_to_idx.find(iid);
        if (it == ped_id_to_idx.end()) {
            throw std::runtime_error("Genotyped sample is absent from pedigree: " + fid_iid);
        }
        if (!mapped_pedigree_rows.insert(it->second).second) {
            throw std::runtime_error("Multiple genotype samples map to pedigree ID: " + iid);
        }
        genotyped_indices_.push_back(it->second);
        genotyped_ids_.push_back(ped_records[it->second].id);
    }
    std::cout << "Genotyped animals: " << n_genotyped_ << " (pedigree-mapped)\n";
}

// ---------------------------------------------------------------------------
// Phase E2: resolve the block-streamed GRM path from --ssgwas-grm-mode.
// "auto" selects streaming once the genotyped count crosses the dense/Lanczos
// threshold (the same heuristic the auto-core eigensolver uses); "streamed"
// forces it (APY methods only); "dense" always builds the full G.
// ---------------------------------------------------------------------------
void SingleStepModel::resolveStreamedGrmMode() {
    const std::string& mode = options_.ssgwas_grm_mode;
    if (mode != "auto" && mode != "dense" && mode != "streamed") {
        throw std::runtime_error("Invalid --ssgwas-grm-mode value '" + mode +
                                 "' (expected auto|dense|streamed)");
    }
    const bool is_apy = options_.ssgwas_pvalue_method == "apy-exact" ||
                        options_.ssgwas_pvalue_method == "apy-approx";
    grm_streamed_ = false;
    if (mode == "streamed") {
        if (!is_apy) {
            throw std::runtime_error(
                "--ssgwas-grm-mode streamed requires an APY p-value method "
                "(--ssgwas-pvalue-method apy-exact|apy-approx); the exact path "
                "needs the full dense G_bs_.");
        }
        grm_streamed_ = true;
    } else if (mode == "auto") {
        const int STREAM_THRESHOLD = 2000;
        grm_streamed_ = is_apy && n_genotyped_ >= STREAM_THRESHOLD;
    }
    if (grm_streamed_ && options_.grm_alg != 1) {
        throw std::runtime_error(
            "--ssgwas-grm-mode streamed requires --grm-alg 1 (VanRaden) in Phase E2");
    }
    if (grm_streamed_) {
        std::cout << "Block-streamed GRM (Phase E2): Gcc/Gcn/diag(Gnn) derived from the "
                  << "genotype stream without materializing the full n_g^2 G.\n";
    } else {
        std::cout << "Dense GRM: full n_g^2 G built.\n";
    }
}

// ---------------------------------------------------------------------------
// Phase E2 pass 1 (block-streamed): allele statistics, S = sum 2p(1-p), the
// raw trace / off-diagonal sums of the genotyped block, the tuneG scalars and
// the A22 cache.  Pass 2 (computeBlocks) runs only after the APY core
// partition exists (runBacksolve -> buildBacksolveGrmStreamBlocks).
// ---------------------------------------------------------------------------
void SingleStepModel::buildBacksolveGrmStreamed() {
    BlockStreamedGrmBuilder::Options bopts;
    bopts.min_maf = (float)options_.min_maf;
    bopts.min_hwe = (float)options_.min_hwe;
    bopts.extract_set = loadSnpIdSet(options_.extract_file);
    bopts.exclude_set = loadSnpIdSet(options_.exclude_file);
    bopts.n_threads = options_.n_threads > 0 ? options_.n_threads : 1;
    bopts.materialize_gcn = options_.ssgwas_apy_inverse_mode != "operator" &&
                            options_.ssgwas_mme_apy_mode != "operator" &&
                            options_.ssgwas_pev_apy_mode != "operator";
    grm_stream_builder_ =
        std::make_unique<BlockStreamedGrmBuilder>(*plink_reader_, geno_full_idx_, bopts);

    grm_stream_builder_->computeStats();
    S_2pq_ = grm_stream_builder_->S2pq();
    std::cout << "G (streamed): " << n_genotyped_ << " genotyped animals, "
              << "S = sum 2p(1-p) = " << S_2pq_ << "\n";

    const double ngd = (double)n_genotyped_;
    const double n_off_raw = ngd * (ngd - 1.0);
    g_raw_mean_diag_ = grm_stream_builder_->rawDiagSum() / ngd;
    g_raw_mean_off_ = n_off_raw > 0.0 ? grm_stream_builder_->rawOffSum() / n_off_raw : 0.0;

    // A22 of the genotyped animals is materialized only when tuneG or the
    // A22 blend needs it (dense n_g^2 for now; Phase E3 removes that too).
    A22_geno_.resize(0, 0);
    if (options_.ssgwas_tune_g || options_.ssgwas_blend_a22) {
        A22_geno_ = pedigree_.buildA22(genotyped_indices_);
    }

    tune_a_ = 0.0;
    tune_b_ = 1.0;
    if (options_.ssgwas_tune_g) {
        double sum_A_diag = 0.0, sum_A_off = 0.0;
        for (int i = 0; i < n_genotyped_; ++i) {
            sum_A_diag += A22_geno_(i, i);
            for (int j = 0; j < n_genotyped_; ++j) {
                if (i != j) sum_A_off += A22_geno_(i, j);
            }
        }
        const double mean_G_diag = g_raw_mean_diag_;
        const double mean_A_diag = sum_A_diag / ngd;
        const double mean_G_off = g_raw_mean_off_;
        const double mean_A_off = n_off_raw > 0.0 ? sum_A_off / n_off_raw : 0.0;
        tune_b_ = 1.0;
        if (std::abs(mean_G_diag - mean_G_off) > 1e-12) {
            tune_b_ = (mean_A_diag - mean_A_off) / (mean_G_diag - mean_G_off);
        }
        tune_a_ = mean_A_diag - tune_b_ * mean_G_diag;
        a22_mean_diag_ = mean_A_diag;
        a22_mean_off_ = mean_A_off;
        std::cout << "tuneG: aligning G = b*G + a to A22 with a=" << tune_a_
                  << ", b=" << tune_b_ << "\n";
    }
    grm_stream_builder_->setTune(tune_a_, tune_b_);
}

// ---------------------------------------------------------------------------
// Phase E2 pass 2 (block-streamed): accumulate Gcc/Gcn/diag(Gnn) over the
// genotyped rows, apply tuneG + blend block-wise, then build the MME/effect
// APY block inverses from the blocks (ApyGinv::build(blocks)).  Also updates
// the final blended mean diagnostics and prints the matrix footprint.
// ---------------------------------------------------------------------------
void SingleStepModel::buildBacksolveGrmStreamBlocks() {
    if (!grm_stream_builder_ || core_indices_.empty()) {
        throw std::runtime_error("buildBacksolveGrmStreamBlocks: no stream stats / core partition");
    }
    grm_stream_builder_->computeBlocks(core_indices_, noncore_indices_);
    std::cout << "E3-D3 stream audit: " << grm_stream_builder_->streamAuditSummary() << "\n";

    const double alpha = options_.blend_g > 0.0 ? options_.blend_g : 0.0;
    const bool blend_a22 = options_.ssgwas_blend_a22;
    if (blend_a22 && A22_geno_.size() != (long)n_genotyped_ * n_genotyped_) {
        A22_geno_ = pedigree_.buildA22(genotyped_indices_);
    }
    if (blend_a22 && A22_geno_.cols() == n_genotyped_) {
        const double sum_A_diag = A22_geno_.trace();
        double sum_A_off = 0.0;
        for (int i = 0; i < n_genotyped_; ++i) {
            for (int j = 0; j < n_genotyped_; ++j) {
                if (i != j) sum_A_off += A22_geno_(i, j);
            }
        }
        a22_mean_diag_ = sum_A_diag / n_genotyped_;
        const double n_off = (double)n_genotyped_ * (n_genotyped_ - 1);
        a22_mean_off_ = n_off > 0.0 ? sum_A_off / n_off : 0.0;
        blend_target_name_ = "a22";
    } else if (blend_a22) {
        blend_target_name_ = "a22";
    }
    if (!blend_a22) blend_target_name_ = "identity";

    grm_stream_builder_->finalize(alpha, blend_a22 ? &A22_geno_ : nullptr);

    // Hold only Gcc and diag(Gnn). In operator mode Gcn is regenerated from
    // genotypes for each consumer and never persists as K x n_noncore state.
    gcc_stream_final_ = grm_stream_builder_->gccFinal();
    gnn_diag_stream_final_ = grm_stream_builder_->gnnDiagFinal();
    gcc_stream_tuned_ = grm_stream_builder_->gccTuned();
    gnn_diag_stream_tuned_ = grm_stream_builder_->gnnDiagTuned();

    // Phase E3-D2D: with any APY operator backend active the Gcn block is NOT
    // materialized (see buildBacksolveGrmStreamed); replace the per-apply full
    // K x Nn regeneration of the MME / PEV operators with a matrix-free cross
    // operator built on the reader primitives.  The dense / materialized
    // reference paths keep the block-iterator behavior (cross stays null).
    const bool apy_operator_mode =
        options_.ssgwas_apy_inverse_mode == "operator" ||
        options_.ssgwas_mme_apy_mode == "operator" ||
        options_.ssgwas_pev_apy_mode == "operator";
    if (apy_operator_mode && apy_n_noncore_ > 0) {
        std::vector<int> core_keep(apy_n_core_), noncore_keep(apy_n_noncore_);
        for (int a = 0; a < apy_n_core_; ++a) {
            core_keep[a] = geno_full_idx_[core_indices_[a]];
        }
        for (int b = 0; b < apy_n_noncore_; ++b) {
            noncore_keep[b] = geno_full_idx_[noncore_indices_[b]];
        }
        const double alpha = options_.blend_g > 0.0 ? options_.blend_g : 0.0;
        const bool blend_a22 = options_.ssgwas_blend_a22;
        const int threads = options_.n_threads > 0 ? options_.n_threads : 1;
        gcn_cross_final_ = std::make_unique<CrossGenomicOperator>(
            *plink_reader_, core_keep, noncore_keep, S_2pq_, tune_a_, tune_b_,
            alpha, GcnBlockView::Final, blend_a22 ? &A22_geno_ : nullptr,
            core_indices_, noncore_indices_, threads,
            &grm_stream_builder_->snpValidMask(),
            &grm_stream_builder_->snpFrequencies());
        if (blend_a22) {
            gcn_cross_tuned_ = std::make_unique<CrossGenomicOperator>(
                *plink_reader_, core_keep, noncore_keep, S_2pq_, tune_a_,
                tune_b_, 0.0, GcnBlockView::Tuned, nullptr, core_indices_,
                noncore_indices_, threads, &grm_stream_builder_->snpValidMask(),
                &grm_stream_builder_->snpFrequencies());
        }
        gcn_regen_after_construct_ = grm_stream_builder_->gcnRegenCalls();
        std::cout << "Phase E3-D2D: matrix-free Gcn cross operator active (n_core="
                  << apy_n_core_ << " n_noncore=" << apy_n_noncore_
                  << ", view=" << (blend_a22 ? "final+tuned" : "final")
                  << ", gcn_regen_after_construct=" << gcn_regen_after_construct_ << ")\n";
    }

    const double ngd = (double)n_genotyped_;
    const double n_off = ngd * (ngd - 1.0);
    g_final_mean_diag_ = grm_stream_builder_->finalDiagSum() / ngd;
    g_final_mean_off_ = n_off > 0.0 ? grm_stream_builder_->finalOffSum() / n_off : 0.0;

    {
        const size_t bytes = gcc_stream_final_.rows() * gcc_stream_final_.cols() * sizeof(double) +
                             grm_stream_builder_->gcnFinal().size() * sizeof(double) +
                             gnn_diag_stream_final_.size() * sizeof(double) +
                             (blend_a22 ? gcc_stream_tuned_.rows() * gcc_stream_tuned_.cols() * sizeof(double) +
                                              grm_stream_builder_->gcnTuned().size() * sizeof(double) +
                                              gnn_diag_stream_tuned_.size() * sizeof(double)
                                        : 0);
        std::cout << "Block-streamed G: n_core=" << apy_n_core_ << " n_noncore=" << apy_n_noncore_
                  << " held blocks = " << (bytes / (1024.0 * 1024.0)) << " MiB (no full G; full Gcn="
                  << (grm_stream_builder_->fullGcnStored() ? "yes" : "no") << ").\n";
    }

    // APY block inverses: MME on the final blended blocks; effect matrix on the
    // tuned-unblended blocks when the A22-blend reference path is active.
    if (options_.ssgwas_mme_apy_mode == "dense") {
        apy_ginv_mme_ = std::make_unique<ApyGinv>();
        if (grm_stream_builder_->fullGcnStored()) {
            apy_ginv_mme_->build(gcc_stream_final_, grm_stream_builder_->gcnFinal(),
                                 gnn_diag_stream_final_, n_genotyped_, core_indices_,
                                 "APY G (MME; final Gw; materialized Gcn)");
        } else {
            GcnBlockIteratorFactory final_factory = [this]() {
                return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Final);
            };
            apy_ginv_mme_->build(gcc_stream_final_, gnn_diag_stream_final_, n_genotyped_,
                                 core_indices_, "APY G (MME; genotype-streamed Gcn)", final_factory);
        }
    }
    if (blend_a22) {
        if (options_.ssgwas_apy_inverse_mode == "dense") {
            apy_ginv_eff_ = std::make_unique<ApyGinv>();
            if (grm_stream_builder_->fullGcnStored()) {
                apy_ginv_eff_->build(gcc_stream_tuned_, grm_stream_builder_->gcnTuned(),
                                     gnn_diag_stream_tuned_, n_genotyped_, core_indices_,
                                     "APY G* (effect; tuned; materialized Gcn)");
            } else {
                GcnBlockIteratorFactory tuned_factory = [this]() {
                    return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Tuned);
                };
                apy_ginv_eff_->build(gcc_stream_tuned_, gnn_diag_stream_tuned_, n_genotyped_,
                                     core_indices_, "APY G* (effect; genotype-streamed Gcn)",
                                     tuned_factory);
            }
        }
        effect_ginv_name_ = "preblend-apy";
    } else {
        if (options_.ssgwas_apy_inverse_mode == "dense") {
            apy_ginv_eff_ = std::make_unique<ApyGinv>();
            if (grm_stream_builder_->fullGcnStored()) {
                apy_ginv_eff_->build(gcc_stream_final_, grm_stream_builder_->gcnFinal(),
                                     gnn_diag_stream_final_, n_genotyped_, core_indices_,
                                     "APY G (effect; materialized Gcn)");
            } else {
                GcnBlockIteratorFactory final_factory = [this]() {
                    return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Final);
                };
                apy_ginv_eff_->build(gcc_stream_final_, gnn_diag_stream_final_, n_genotyped_,
                                     core_indices_, "APY G (effect; genotype-streamed Gcn)",
                                     final_factory);
            }
        }
        effect_ginv_name_ = "apy";
    }
}

void SingleStepModel::materializeApyMmeReference() {
    if (apy_ginv_mme_) return;
    apy_ginv_mme_ = std::make_unique<ApyGinv>();
    if (!grm_streamed_) {
        apy_ginv_mme_->build(G_bs_, core_indices_, "APY G (PEV reference; final Gw)");
        return;
    }
    if (grm_stream_builder_->fullGcnStored()) {
        apy_ginv_mme_->build(gcc_stream_final_, grm_stream_builder_->gcnFinal(),
                             gnn_diag_stream_final_, n_genotyped_, core_indices_,
                             "APY G (PEV reference; materialized Gcn)");
    } else {
        GcnBlockIteratorFactory final_factory = [this]() {
            return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Final);
        };
        apy_ginv_mme_->build(gcc_stream_final_, gnn_diag_stream_final_, n_genotyped_,
                             core_indices_, "APY G (PEV reference; streamed Gcn)", final_factory);
    }
}

// ---------------------------------------------------------------------------
// Phase E2: eigen-spectrum C98 of the final blended G via the matrix-free
// Lanczos operator (E1).  Approximate spectrum (Lanczos); returns the C98
// cumulative-variance count with the same semantics as the dense diagnostic.
// ---------------------------------------------------------------------------
int SingleStepModel::runStreamedK98() const {
    std::unique_ptr<Eigen::MatrixXd> a22_blend;
    if (options_.ssgwas_blend_a22 && A22_geno_.cols() == n_genotyped_) {
        a22_blend = std::make_unique<Eigen::MatrixXd>(A22_geno_);
    }
    MatrixFreeGrmOperator G_op(*plink_reader_, geno_full_idx_, S_2pq_, tune_a_, tune_b_,
                               options_.blend_g > 0.0 ? options_.blend_g : 0.0,
                               options_.ssgwas_blend_a22 ? a22_blend.get() : nullptr,
                               options_.n_threads > 0 ? options_.n_threads : 1);
    LanczosEigenSolver lanczos;
    LanczosEigenSolver::Options lopts;
    lopts.max_iterations = std::min(80, n_genotyped_);
    lopts.tolerance = 1e-6;
    lopts.seed = static_cast<unsigned long>(options_.ssgwas_apy_seed);
    lopts.n_threads = options_.n_threads > 0 ? options_.n_threads : 1;
    LanczosResult lres = lanczos.compute(G_op, k98_cumfraction_, lopts);
    if (lres.status != "OK") return -1;
    if (k98_cumfraction_ < 0.98) return lres.K_target < 1 ? 1 : lres.K_target;
    return lres.K98 < 1 ? 1 : lres.K98;
}

// ---------------------------------------------------------------------------
// Final hybrid G core-block element in genotyped order (core row / core col):
// the streamed path never materializes G_bs_, so the core block is served from
// the held Gcc block; the dense path reads the precomputed core block of G_bs_.
// ---------------------------------------------------------------------------
double SingleStepModel::gwCc(int a, int b) const {
    if (grm_streamed_) return gcc_stream_final_(a, b);
    return G_bs_(core_indices_[a], core_indices_[b]);
}

// ---------------------------------------------------------------------------
// Dense back-solving GRM (full n_g^2 G; the Phase E2 parity reference).
// ---------------------------------------------------------------------------
void SingleStepModel::buildBacksolveGrm() {
    // 3. G over the FULL genotyped sample set (not the phenotype keep set).
    //    Full-sample computation + subset is used so that G's centering and the
    //    Z used by multiply_Zt_v (standardized with full-file snp_means/snp_inv_stds)
    //    refer to the same allele frequencies, keeping the back-solving identity
    //    a_hat_g = M * a  <->  a = M' G^-1 a_hat_g / S exact even under
    //    --genotyped-list subsets.
    GrmBuilder::Options grm_opts;
    grm_opts.algorithm = options_.grm_alg;
    grm_opts.min_maf = (float)options_.min_maf;
    grm_opts.min_hwe = (float)options_.min_hwe;
    grm_opts.extract_set = loadSnpIdSet(options_.extract_file);
    grm_opts.exclude_set = loadSnpIdSet(options_.exclude_file);
    grm_opts.thread_num = options_.n_threads;
    GrmBuilder grm_builder(grm_opts);
    grm_builder.compute(*plink_reader_);
    DataFilter::KeepInfo geno_keep;
    geno_keep.keep_indices = geno_full_idx_;
    geno_keep.n_keep = (size_t)n_genotyped_;
    G_bs_ = DataFilter::subsetSquare(grm_builder.getGrm(), geno_keep);
    S_2pq_ = grm_builder.getSum2pq();
    if (G_bs_.rows() != n_genotyped_) {
        throw std::runtime_error("Back-solving GRM dimension does not match genotyped sample count");
    }
    std::cout << "G: " << n_genotyped_ << "x" << n_genotyped_
              << ", S = sum 2p(1-p) = " << S_2pq_ << "\n";

    // 4. tuneG: align G to A22 (G = b*G + a), preGSf90-style (CosmicSolver
    //    matrix_builder). Applied to the RAW G (before blending) so that the
    //    tuned matrix G* used by back-solving is the pre-blend alignment,
    //    matching reference implementations (HIBLUP/preGSf90: tune then blend).
    g_raw_mean_diag_ = 0.0;
    g_raw_mean_off_ = 0.0;
    for (int i = 0; i < n_genotyped_; ++i) {
        g_raw_mean_diag_ += G_bs_(i, i);
        for (int j = 0; j < n_genotyped_; ++j) {
            if (i != j) g_raw_mean_off_ += G_bs_(i, j);
        }
    }
    double n_off_raw = (double)n_genotyped_ * (n_genotyped_ - 1);
    g_raw_mean_off_ = n_off_raw > 0 ? g_raw_mean_off_ / n_off_raw : 0.0;
    g_raw_mean_diag_ /= n_genotyped_;
    if (options_.ssgwas_tune_g) {
        Eigen::MatrixXd A22 = pedigree_.buildA22(genotyped_indices_);
        double sum_G_diag = 0.0, sum_A_diag = 0.0, sum_G_off = 0.0, sum_A_off = 0.0;
        for (int i = 0; i < n_genotyped_; ++i) {
            sum_G_diag += G_bs_(i, i); sum_A_diag += A22(i, i);
            for (int j = 0; j < n_genotyped_; ++j) {
                if (i != j) { sum_G_off += G_bs_(i, j); sum_A_off += A22(i, j); }
            }
        }
        double mean_G_diag = sum_G_diag / n_genotyped_;
        double mean_A_diag = sum_A_diag / n_genotyped_;
        double n_off = (double)n_genotyped_ * (n_genotyped_ - 1);
        double mean_G_off = n_off > 0 ? sum_G_off / n_off : 0.0;
        double mean_A_off = n_off > 0 ? sum_A_off / n_off : 0.0;
        tune_b_ = 1.0;
        if (std::abs(mean_G_diag - mean_G_off) > 1e-12) {
            tune_b_ = (mean_A_diag - mean_A_off) / (mean_G_diag - mean_G_off);
        }
        tune_a_ = mean_A_diag - tune_b_ * mean_G_diag;
        a22_mean_diag_ = mean_A_diag;
        a22_mean_off_ = mean_A_off;
        std::cout << "tuneG: aligning G = b*G + a to A22 with a=" << tune_a_
                  << ", b=" << tune_b_ << "\n";
        G_bs_ = tune_b_ * G_bs_;
        G_bs_.array() += tune_a_;
    }

    // 5. Blend (Method-1 convention: (1-alpha) G + alpha I; with
    //    --ssgwas-blend-a22 the reference-compatible (1-alpha) G + alpha A22)
    if (options_.blend_g > 0.0) {
        double alpha = options_.blend_g;
        if (options_.ssgwas_blend_a22) {
            G_bs_tuned_ = G_bs_;  // pre-blend copy used for SNP back-solving
            effect_ginv_name_ = "preblend";
            Eigen::MatrixXd A22 = pedigree_.buildA22(genotyped_indices_);
            double sum_A_diag = 0.0, sum_A_off = 0.0;
            for (int i = 0; i < n_genotyped_; ++i) {
                sum_A_diag += A22(i, i);
                for (int j = 0; j < n_genotyped_; ++j) {
                    if (i != j) sum_A_off += A22(i, j);
                }
            }
            a22_mean_diag_ = sum_A_diag / n_genotyped_;
            double n_off = (double)n_genotyped_ * (n_genotyped_ - 1);
            a22_mean_off_ = n_off > 0 ? sum_A_off / n_off : 0.0;
            G_bs_ = (1.0 - alpha) * G_bs_ + alpha * A22;
            std::cout << "Back-solving G blended with A22 (alpha=" << alpha
                      << "), A22 diag/off mean = " << a22_mean_diag_
                      << " / " << a22_mean_off_ << "\n";
            blend_target_name_ = "a22";
        } else {
            G_bs_ = (1.0 - alpha) * G_bs_ + alpha * Eigen::MatrixXd::Identity(n_genotyped_, n_genotyped_);
            std::cout << "Back-solving G blended with identity (alpha=" << alpha << ")\n";
            blend_target_name_ = "identity";
        }
    }

    // 6. The final G enters H^-1 as the (possibly blended) matrix; the tuned-
    //    but-unblended copy G* is retained for back-solving when the A22-blend
    //    reference path is active (HIBLUP computes SNP effects with G*^-1,
    //    while its MME uses the blended Gw^-1).
    g_final_mean_diag_ = 0.0;
    g_final_mean_off_ = 0.0;
    for (int i = 0; i < n_genotyped_; ++i) {
        g_final_mean_diag_ += G_bs_(i, i);
        for (int j = 0; j < n_genotyped_; ++j) {
            if (i != j) g_final_mean_off_ += G_bs_(i, j);
        }
    }
    double n_off_tot = (double)n_genotyped_ * (n_genotyped_ - 1);
    g_final_mean_off_ = n_off_tot > 0 ? g_final_mean_off_ / n_off_tot : 0.0;
    g_final_mean_diag_ /= n_genotyped_;
}

// ---------------------------------------------------------------------------
// MME:  [X'X   X'Z   ] [b]   [X'y]
//       [Z'X  Z'Z + Hinv*lambda] [a] = [Z'y],  lambda = sigma2_e / sigma2_a
// Hinv = Ainv + [tau*Ginv - omega*A22inv] on the genotyped block
// (SplitMatrixAdapter, same formula as CosmicSolver Hinv construction).
// ---------------------------------------------------------------------------
void SingleStepModel::solveGeneticEvaluation() {
    n_ped_ = pedigree_.nPedigree();
    const auto& ped_records = pedigree_.records();

    // Fixed design: intercept only (y = mu + a + e) in Phase 6 v1.
    FixedDesignG fd;
    fd.p = 1;
    fd.names = {"mu"};

    // Records: one per phenotype observation; aid = 1-based pedigree row index.
    const int n_rec = (int)pheno_ped_rows_.size();
    std::vector<GenRecord> recs;
    recs.reserve(n_rec);
    Eigen::VectorXd y_rec(n_rec);
    for (int k = 0; k < n_rec; ++k) {
        GenRecord g;
        g.aid = pheno_ped_rows_[k] + 1;  // 1-based (MMELHSBuilder uses recs[i].aid - 1)
        g.y = pheno_y_[k];
        g.idstr = ped_records[pheno_ped_rows_[k]].id;
        recs.push_back(g);
        y_rec(k) = pheno_y_[k];
        fd.rows.push_back({{0, 1.0}});
    }
    n_phenotyped_ = n_rec;

    // Phase D: diagonal effective records (Leite et al. 2024 Eq. 9) for the
    // approximate core-PEV route only; harmless (O(n_rec*p^2)) in all modes.
    if (options_.ssgwas_pvalue_method == "apy-approx") {
        computeEffectiveRecordDiagonals(recs, fd);
    }

    // Counts (meta / EBV annotation)
    n_genotyped_with_pheno_ = 0;
    n_ungenotyped_with_pheno_ = 0;
    std::unordered_set<int> geno_rows(genotyped_indices_.begin(), genotyped_indices_.end());
    for (int r : pheno_ped_rows_) {
        if (geno_rows.count(r)) ++n_genotyped_with_pheno_;
        else ++n_ungenotyped_with_pheno_;
    }
    std::cout << "ssGBLUP setup: n_pedigree=" << n_ped_
              << " n_phenotyped=" << n_phenotyped_
              << " (genotyped=" << n_genotyped_with_pheno_
              << ", ungenotyped=" << n_ungenotyped_with_pheno_ << ")\n";

    // Hinv = Ainv + [tau*Ginv - omega*A22inv] on the genotyped block. The
    // legacy path retains dense A22inv; E4 operator mode applies its exact
    // sparse Schur complement without materializing A22/A22inv/Hinv.
    SparseMatrixAdapter ainv_ada(A_inv_);
    Eigen::MatrixXd ginv_scaled;
    DenseMatrixAdapter ginv_ada(ginv_scaled);
    std::unique_ptr<ApyGinvAdapter> apy_ginv_ada;
    std::unique_ptr<APYMMEOperator> apy_mme_operator_ada;
    AbstractMatrix* ginv_ptr = nullptr;
    if (options_.ssgwas_pvalue_method == "apy-exact" ||
        options_.ssgwas_pvalue_method == "apy-approx") {
        // Phase C2/D: the MME genotyped block uses the APY block inverse of Gw
        // (Leite et al. 2024 Eq. 5); no dense n_g^2 G^-1 is materialized.
        if (options_.ssgwas_mme_apy_mode == "operator") {
            apy_mme_operator_ada = std::make_unique<APYMMEOperator>(
                apy_inverse_mme_operator_.get(), options_.h_tau);
            ginv_ptr = apy_mme_operator_ada.get();
        } else {
            apy_ginv_ada = std::make_unique<ApyGinvAdapter>(apy_ginv_mme_.get(), options_.h_tau);
            ginv_ptr = apy_ginv_ada.get();
        }
    } else {
        ginv_scaled = options_.h_tau * Ginv_bs_;
        ginv_ada = DenseMatrixAdapter(ginv_scaled);
        ginv_ptr = &ginv_ada;
    }

    std::unique_ptr<DenseMatrixAdapter> dense_a22inv_ada;
    std::unique_ptr<AbstractMatrix> hinv_ada;
    if (options_.ssgwas_hinv_mode == "operator") {
        if (!a22_inverse_operator_) {
            throw std::runtime_error("Hinv operator mode requires a constructed A22InverseOperator.");
        }
        auto operator_hinv = std::make_unique<HInverseOperator>(
            A_inv_, ginv_ptr, a22_inverse_operator_.get(), genotyped_indices_, options_.h_omega);
        hinv_preconditioner_diagonal_mode_ = operator_hinv->preconditionerDiagonalMode();
        hinv_ada = std::move(operator_hinv);
    } else {
        Eigen::MatrixXd a22inv_scaled = options_.h_omega * A22_inv_;
        dense_a22inv_ada = std::make_unique<DenseMatrixAdapter>(std::move(a22inv_scaled));
        hinv_ada = std::make_unique<SplitMatrixAdapter>(
            &ainv_ada, ginv_ptr, dense_a22inv_ada.get(), genotyped_indices_);
        hinv_preconditioner_diagonal_mode_ = "dense_exact";
    }
    RandomComponent comp;
    comp.Qinv = hinv_ada.get();  // id_map empty -> recs[i].aid - 1

    // Variance components
    if (!options_.ssgwas_vce.empty()) {
        VCEConfig vc_cfg;
        vc_cfg.set_algorithm_from_string(options_.ssgwas_vce);
        std::vector<RandomComponent> comps = {comp};
        auto vce = VCESolver::create(recs, fd, comps, y_rec, vc_cfg);
        vce->solve();
        sigma2_a_backsolve_ = vce->getVarU();
        sigma2_e_backsolve_ = vce->getVarE();
        vce_mode_used_ = options_.ssgwas_vce;
        if (!vce->isConverged()) {
            std::cout << "WARNING: estimates for --ssgwas-vce did not converge (last-iteration used)\n";
        }
        std::cout << "VCE (" << vce_mode_used_ << "): Va=" << sigma2_a_backsolve_
                  << " Ve=" << sigma2_e_backsolve_ << "\n";
    } else {
        sigma2_a_backsolve_ = options_.reml_vg;
        sigma2_e_backsolve_ = options_.reml_ve;
        vce_mode_used_ = "fixed(" + std::to_string(options_.reml_vg) + "," +
                         std::to_string(options_.reml_ve) + ")";
        std::cout << "Fixed VC: Va=" << sigma2_a_backsolve_
                  << " Ve=" << sigma2_e_backsolve_ << "\n";
    }
    if (!(sigma2_a_backsolve_ > 0.0) || !(sigma2_e_backsolve_ > 0.0)) {
        throw std::runtime_error("ssGBLUP requires positive variance components Va > 0 and Ve > 0");
    }
    double lambda = sigma2_e_backsolve_ / sigma2_a_backsolve_;

    // MME
    auto t0 = std::chrono::steady_clock::now();
    std::vector<RandomComponent> comps = {comp};
    Eigen::VectorXd sol;
    if (options_.ssgwas_solver == "pcg") {
        ImplicitMME implicit_mme(recs, fd, comps, {lambda});
        MMELHSBuilder rhs_builder(recs, fd, comps, /*build_matrix=*/false);
        Eigen::VectorXd rhs = rhs_builder.build_rhs(recs, fd, implicit_mme.dim);
        MatrixFreePCG pcg(implicit_mme.dim);
        long long mult_calls = 0;
        pcg.setMultOp([&](const Eigen::VectorXd& v, Eigen::VectorXd& out) {
            ++mult_calls;
            implicit_mme.multiply(v, out);
        });
        pcg.setPrecondOp([&](const Eigen::VectorXd& v, Eigen::VectorXd& out) {
            out = v.cwiseQuotient(implicit_mme.diag);
        });
        sol = pcg.solve(rhs, options_.ssgwas_pcg_tol, options_.ssgwas_pcg_maxit);
        Eigen::VectorXd Ax(implicit_mme.dim);
        implicit_mme.multiply(sol, Ax);
        solver_residual_ = (rhs - Ax).norm() / (1.0 + rhs.norm());
        solver_iterations_ = (int)(mult_calls - 1);
        const bool pcg_converged_early = solver_iterations_ < options_.ssgwas_pcg_maxit;
        solver_converged_ = pcg_converged_early || (solver_residual_ <= options_.ssgwas_pcg_tol);
        std::cout << "MME solved with matrix-free PCG in " << solver_iterations_
                  << " iterations; relative residual=" << solver_residual_
                  << " converged=" << (solver_converged_ ? "yes" : "no") << "\n";
    } else {
        MMELHSBuilder builder(recs, fd, comps, /*build_matrix=*/true);
        Eigen::SparseMatrix<double> lhs = builder.build_lhs({lambda});
        Eigen::VectorXd rhs = builder.build_rhs(recs, fd, builder.get_dim());
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(lhs);
        if (ldlt.info() != Eigen::Success) {
            throw std::runtime_error("ssGBLUP MME SimplicialLDLT factorization failed");
        }
        sol = ldlt.solve(rhs);
        Eigen::VectorXd Ax = lhs * sol;
        solver_residual_ = (rhs - Ax).norm() / (1.0 + rhs.norm());
        solver_converged_ = solver_residual_ <= 1e-6;
        solver_iterations_ = 0;
        std::cout << "MME solved with SimplicialLDLT; dimension=" << sol.size()
                  << " relative residual=" << solver_residual_ << "\n";
        if (options_.ssgwas_pvalue_method == "exact" ||
            options_.ssgwas_pvalue_method == "apy-exact") {
            mme_lhs_ = lhs;  // C_lambda retained for the Phase C1/C2 selected inverse
        }
    }
    t_solve_ = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    beta_hat_ = sol.head(fd.p);
    a_hat_ = sol.tail(sol.size() - fd.p);
    if (a_hat_.size() != n_ped_) {
        throw std::runtime_error("Internal error: breeding-value vector size does not match pedigree size");
    }

    // Genotyped EBV (explicit ID mapping, not row-guessed)
    a_hat_g_.resize(n_genotyped_);
    for (int i = 0; i < n_genotyped_; ++i) a_hat_g_(i) = a_hat_(genotyped_indices_[i]);
}

// ---------------------------------------------------------------------------
// Back-solving:  q = Ginv * a_hat_g;  zt = Z' q;  effect_allele = zt/(std*S)
// ---------------------------------------------------------------------------
void SingleStepModel::backsolveSnpEffects() {
    std::cout << "Back-solving SNP effects (q = G^-1 * a_hat_g -> Z' q)...\n";
    Eigen::VectorXd q;
    if (options_.ssgwas_pvalue_method == "apy-exact" &&
        options_.ssgwas_apy_exact_estimator == "legacy-hybrid-v1") {
        // Phase C2: q = G_APY^-1 a_hat_g (full APY product; exact given the
        // APY partition, Leite et al. 2024 Eq. 5).
        q = apy_inverse_eff_operator_->apply(a_hat_g_);
    } else if (options_.ssgwas_pvalue_method == "apy-approx" ||
               (options_.ssgwas_pvalue_method == "apy-exact" &&
                options_.ssgwas_apy_exact_estimator == "core-v2")) {
        // R2B core-v2 and unchanged Phase D: q = Gcc^-1 (a_hat_g)_c on
        // the core, zeros on the noncore. Exact selected-core covariance is
        // retained for core-v2; only Phase D substitutes approximate PEV.
        // (Bermann 2022a Eq. 21 core-only operator; the same (Gcc^-1)_cc
        // operator drives the Eq. 8/10 variance, so effect and uncertainty
        // share one back-solve operator).
        q = Eigen::VectorXd::Zero(n_genotyped_);
        Eigen::VectorXd a_core(apy_n_core_);
        for (int i = 0; i < apy_n_core_; ++i) {
            a_core(i) = a_hat_g_(core_indices_[i]);
        }
        Eigen::VectorXd q_core = apy_inverse_eff_operator_->applyCoreInverse(a_core);
        for (int i = 0; i < apy_n_core_; ++i) {
            q(core_indices_[i]) = q_core(i);
        }
    } else {
        const Eigen::MatrixXd& ginv_eff =
            options_.ssgwas_blend_a22 ? Ginv_bs_tuned_ : Ginv_bs_;
        q = ginv_eff * a_hat_g_;
    }

    // Z = standardized genotype matrix over the genotyped samples; the subset
    // multiply_Zt_v uses v(k) <-> keep_indices[k], so q is in genotyped order.
    Eigen::VectorXd zt_q;
    plink_reader_->multiply_Zt_v(q, zt_q, geno_full_idx_,
                                 options_.n_threads > 0 ? options_.n_threads : 1);
    ++point_multiply_zt_v_calls_;

    // SNP list: all SNPs of the PLINK files; extract/exclude restrict the output.
    const auto& bim = plink_reader_->getBimInfo();
    const auto& means = plink_reader_->getMeans();
    const auto& inv_stds = plink_reader_->getInvStds();
    const std::unordered_set<std::string> extract_set = loadSnpIdSet(options_.extract_file);
    const std::unordered_set<std::string> exclude_set = loadSnpIdSet(options_.exclude_file);

    const int m = (int)zt_q.size();
    snp_keep_indices_.clear();
    snp_keep_indices_.reserve(m);
    for (int j = 0; j < m; ++j) {
        const std::string& id = bim[j].id;
        if (!extract_set.empty() && !extract_set.count(id)) continue;
        if (exclude_set.count(id)) continue;
        snp_keep_indices_.push_back(j);
    }

    snp_effect_std_.resize(m);
    snp_effect_allele_.resize(m);
    snp_variance_.resize(m);
    total_snp_variance_ = 0.0;
    for (int j = 0; j < m; ++j) {
        double zt = zt_q(j);
        double inv_std = j < (int)inv_stds.size() ? inv_stds[j] : 0.0;
        bool std_ok = (inv_std > 0.0 && std::isfinite(inv_std));
        // Undo the (x-mean)/sd standardization of multiply_Zt_v, and the
        // G = M M' / S scaling: allele effect = (M'q) / S = zt * sd / S.
        double effect_allele = (std_ok && S_2pq_ > 1e-12)
            ? (zt / inv_std) / S_2pq_
            : 0.0;
        // Reconstruction-consistent standardized effect: gamma_j = EA_j * sd,
        // i.e. Z * gamma = M * EA = a_hat_g (sd = 1 / inv_std).
        double effect_std = std_ok ? effect_allele / inv_std : 0.0;
        double mean_dos = j < (int)means.size() ? means[j] : NAN;
        double p = (mean_dos >= 0.0 && std::isfinite(mean_dos)) ? mean_dos / 2.0 : 0.0;
        double snp_var = (2.0 * p * (1.0 - p)) * effect_allele * effect_allele;
        snp_effect_std_(j) = effect_std;
        snp_effect_allele_(j) = effect_allele;
        snp_variance_(j) = snp_var;
        total_snp_variance_ += snp_var;
    }
    std::cout << "Back-solved " << snp_keep_indices_.size() << "/" << m
              << " SNP effects (scale: per allele dosage; variance = 2p(1-p)*effect^2).\n";
}

// ---------------------------------------------------------------------------
// Phase C1: frequentist P-values for back-solved SNP effects (exact mode).
// Reference: Aguilar et al. (2019) GSE 51:28 (PMC6584984); terminology per
// the C0 gate (Garcia et al. 2023 correction, PMC10108542):
//   - PEV(a)      = sigma_e2 * (C_lambda^-1)              breeding-value PEV
//   - Var(u_hat_g)= G_eval*sigma_a2 - sigma_e2*Cinv_gg    variance of the EBV
//                     predictor (G_eval = final blended Gw, the H_gg block)
//   - Var(EA_j)   = B_j * Var(u_hat_g) * B_j'             sampling variance of
//                     the SNP-effect estimate (NOT "SNP PEC"); B_j is the
//                     implemented back-solve operator M_j' G_backsolve^-1 / S.
// z_j = EA_j/SE_j, P = 2*(1-Phi(|z|)) = erfc(|z|/sqrt(2)); standard normal.
// ---------------------------------------------------------------------------
void SingleStepModel::acceptExactScanResult(SnpUncertaintyScanResult&& result) {
    snp_se_allele_ = std::move(result.se);
    snp_z_ = std::move(result.z);
    snp_p_ = std::move(result.p);
    snp_log10p_ = std::move(result.log10p);
    snp_pvalue_status_ = std::move(result.status);
    n_pec_clamped_ = result.clamped;
    n_pec_failure_ = result.failures;
    n_zero_se_ = result.zero_se;
    const auto& a = result.audit;
    exact_scan_actual_mode_ = "BATCHED_EXACT_SELECTED";
    t_pev_scan_ = a.time_total;
    t_pev_scan_setup_ = a.time_setup;
    t_pev_scan_onehot_ = a.time_onehot;
    t_pev_scan_genotype_ = a.time_genotype;
    t_pev_scan_gcc_solve_ = a.time_gcc_solve;
    t_pev_scan_gw_matvec_ = a.time_gw_multiply;
    t_pev_scan_pec_matvec_ = a.time_pec_multiply;
    t_pev_scan_quadratic_ = a.time_quadratic;
    t_pev_scan_statistics_ = a.time_statistics;
    t_pev_scan_other_ = a.time_other;
    pev_scan_multiply_z_v_calls_ = a.multiply_z_v_calls;
    pev_scan_genotype_passes_ = a.genotype_passes;
    pev_scan_bed_bytes_read_estimate_ = a.bed_bytes_read_estimate;
    pev_scan_variants_tested_ = a.variants_tested;
    pev_scan_variants_decoded_ = a.variants_decoded;
    pev_scan_allocations_ = a.heap_allocations;
    pev_scan_blocks_ = a.blocks;
    pev_scan_peak_transient_bytes_ = a.peak_transient_bytes;
    pev_scan_block_reallocations_ = a.block_reallocations;
    pev_scan_gcc_solve_calls_ = a.gcc_solve_calls;
    pev_scan_gcc_rhs_columns_ = a.gcc_rhs_columns;
    pev_scan_gw_gemm_calls_ = a.gw_gemm_calls;
    pev_scan_pec_gemm_calls_ = a.pec_gemm_calls;
    pev_scan_missing_to_zero_ = a.missing_to_zero;
    pev_scan_effective_block_size_ = a.block_size;
}

void SingleStepModel::writeExactScanAudit(std::ostream& out) const {
    out << "SNP_SCAN_REQUESTED_MODE\t" << options_.ssgwas_pev_scan_mode << "\n"
        << "SNP_SCAN_MODE_EXPLICIT\t" << (options_.ssgwas_pev_scan_mode_explicit ? "true" : "false") << "\n"
        << "SNP_SCAN_ACTUAL_MODE\t" << exact_scan_actual_mode_ << "\n"
        << "SNP_POINT_MULTIPLY_ZT_V_CALLS\t" << point_multiply_zt_v_calls_ << "\n"
        << "SNP_SCAN_MULTIPLY_Z_V_CALLS\t" << pev_scan_multiply_z_v_calls_ << "\n"
        << "SNP_SCAN_GENOTYPE_PASSES\t" << pev_scan_genotype_passes_ << "\n"
        << "SNP_SCAN_PASS_DEFINITION\tLOGICAL_VARIANCE_PASS (scalar=full panel per SNP; batched=tested rows/variants once; excludes point/setup)\n"
        << "SNP_SCAN_BED_LOGICAL_BYTES\t" << pev_scan_bed_bytes_read_estimate_ << "\n"
        << "SNP_SCAN_VARIANTS_TESTED\t" << pev_scan_variants_tested_ << "\n"
        << "SNP_SCAN_VARIANTS_DECODED\t" << pev_scan_variants_decoded_ << "\n"
        << "SNP_SCAN_BLOCKS\t" << pev_scan_blocks_ << "\n"
        << "SNP_SCAN_EFFECTIVE_BLOCK_SIZE\t" << pev_scan_effective_block_size_ << "\n"
        << "SNP_SCAN_NAMED_ALLOCATION_EVENTS\t" << pev_scan_allocations_ << "\n"
        << "SNP_SCAN_TRANSIENT_BYTES_BOUND\t" << pev_scan_peak_transient_bytes_ << "\n"
        << "SNP_SCAN_ALLOCATION_SCOPE\tNAMED_WORK_BUFFERS_ONLY; backend scratch excluded; process peak is GNU-time RSS\n"
        << "SNP_SCAN_MISSING_CENTERED_ZERO\t" << pev_scan_missing_to_zero_ << "\n"
        << "SNP_SCAN_INVERSE_BATCH_CALLS\t" << pev_scan_gcc_solve_calls_ << "\n"
        << "SNP_SCAN_INVERSE_RHS_COLUMNS\t" << pev_scan_gcc_rhs_columns_ << "\n"
        << "TIME_EXACT_SCAN_DECODE\t" << t_pev_scan_genotype_ << "\n"
        << "TIME_EXACT_SCAN_INVERSE\t" << t_pev_scan_gcc_solve_ << "\n"
        << "TIME_EXACT_SCAN_PRIOR_PRODUCT\t" << t_pev_scan_gw_matvec_ << "\n"
        << "TIME_EXACT_SCAN_COVARIANCE_PRODUCT\t" << t_pev_scan_pec_matvec_ << "\n"
        << "EXACT_SCAN_COMPONENT_TIME_SCOPE\tBATCHED breakdown; scalar historical loop reports total TIME_PEV_SCAN\n";
}

void SingleStepModel::computeSnpPvaluesExact() {
    auto now_t = std::chrono::steady_clock::now;
    auto secs = [](const std::chrono::steady_clock::time_point& t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    const int ng = n_genotyped_;
    const int dim = (int)mme_lhs_.rows();
    if (ng <= 0 || dim <= 0) {
        throw std::runtime_error("Exact SNP p-values require a genotyped subset and a solved MME");
    }

    // --- Memory estimate + guard (n_genotyped x n_genotyped selected block) ---
    estimated_cinv_gg_bytes_ = 8.0 * ng * ng;
    const double kMaxExactBytes = 8ull * 1024ull * 1024ull * 1024ull;  // 8 GiB
    std::cout << "Exact ssGWAS SNP p-values: n_genotyped=" << ng
              << " selected-inverse block = 8*n_g^2 = "
              << estimated_cinv_gg_bytes_ / (1024.0 * 1024.0) << " MiB\n";
    if (estimated_cinv_gg_bytes_ > kMaxExactBytes && !options_.allow_large_exact_ssgwas_pvalue) {
        throw std::runtime_error(
            "Exact ssGWAS SNP p-values would need > 8 GiB for the selected MME-inverse block "
            "(n_genotyped=" + std::to_string(ng) + "); pass --allow-large-exact-ssgwas-pvalue to "
            "override (not recommended; Phase C2 scalable paths are the production route).");
    }

    auto t0 = now_t();

    // --- Selected inverse: (C_lambda^-1)_{gg} via n_g sparse solves ---
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(mme_lhs_);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("ssGBLUP MME SimplicialLDLT factorization failed (PEV step)");
    }
    ++dense_uncertainty_allocs_;  // debug allocation counter (gate D2)
    cinv_gg_.resize(ng, ng);
    Eigen::VectorXd e = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd x(dim);
    for (int j = 0; j < ng; ++j) {
        e(genotyped_indices_[j] + 1) = 1.0;  // column 0 = intercept
        x = ldlt.solve(e);
        e(genotyped_indices_[j] + 1) = 0.0;
        for (int k = 0; k < ng; ++k) cinv_gg_(k, j) = x(genotyped_indices_[k] + 1);
    }
    double max_asym = 0.0;
    for (int j = 0; j < ng; ++j)
        for (int k = 0; k < ng; ++k) {
            double a = std::fabs(cinv_gg_(k, j) - cinv_gg_(j, k));
            if (a > max_asym) max_asym = a;
        }
    cinv_gg_.noalias() = 0.5 * (cinv_gg_ + cinv_gg_.transpose());  // symmetrize
    cinv_gg_max_asym_ = format_sci(max_asym, 6);
    t_pev_inverse_ = secs(t0);

    // --- Var(u_hat_g) = G_eval * sigma_a2 - sigma_e2 * Cinv_gg ---
    // G_eval = G_bs_ (final blended matrix used in H^-1 / MME; H_gg block).
    // G_backsolve (effect operator inverse) = Ginv_bs_tuned_ when blend-a22,
    // else Ginv_bs_ - captured below via the same selection as backsolveSnpEffects.
    t0 = now_t();
    vhat_g_ = sigma2_a_backsolve_ * G_bs_ - sigma2_e_backsolve_ * cinv_gg_;
    t_pev_vhat_ = secs(t0);

    // --- PSD diagnostics of Var(u_hat_g) (C0 gate; exact mode) ---
    // Tolerance is RELATIVE to the matrix scale: entries originate partly from
    // float32 GRM storage (eps32 ~ 1.2e-7 relative precision), so LDLT-noise
    // eigenvalues of order 1e-6*scale (eps32 * safety factor 10) are clamped,
    // while genuinely negative covariance stops the analysis (INVALID_COVARIANCE).
    double max_abs_eig = vhat_g_.diagonal().cwiseAbs().maxCoeff();
    double min_eig = 0.0, max_eig = 1.0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
    if (ng <= 4096) {  // explicit eigen diagnostics on the small exact path
        eig.compute(vhat_g_);
        if (eig.info() == Eigen::Success) {
            min_eig = eig.eigenvalues()(0);
            max_eig = eig.eigenvalues()(ng - 1);
            max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
            var_ghat_eig_mode_ = "EXPLICIT_EIGEN";
        }
    }
    if (var_ghat_eig_mode_.empty()) {
        // Large exact path: diagonal + random-Rayleigh lower bound (no dense eigen).
        min_eig = vhat_g_.diagonal().minCoeff();
        max_eig = vhat_g_.diagonal().maxCoeff();
        for (int r = 0; r < 24; ++r) {
            Eigen::VectorXd v = Eigen::VectorXd::Random(ng).normalized();
            double ra = v.dot(vhat_g_ * v);
            if (ra < min_eig) min_eig = ra;
            if (ra > max_eig) max_eig = ra;
        }
        max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
        var_ghat_eig_mode_ = "RAYLEIGH_PROBE";
    }
    const double psd_tol = 1e-6 * std::max(1.0, max_abs_eig);
    var_ghat_min_eig_ = format_sci(min_eig, 6);
    var_ghat_max_eig_ = format_sci(max_eig, 6);
    std::cout << "Var(u_hat_g) = G_eval*sigma_a2 - sigma_e2*Cinv_gg: min_eig=" << min_eig
              << " max_eig=" << max_eig << " psd_tol=" << psd_tol
              << " (" << var_ghat_eig_mode_ << ")\n";
    if (min_eig < -psd_tol) {
        throw std::runtime_error(
            "Exact SNP p-values: Var(u_hat_g) residual covariance has min eigenvalue "
            + format_sci(min_eig, 6) + " < -tol " + format_sci(psd_tol, 6) +
            "; INVALID_COVARIANCE - refusing to produce SNP P-values (check variance "
            "components / tuning / blending consistency).");
    }

    if (options_.ssgwas_pev_scan_mode_explicit && options_.ssgwas_pev_scan_mode == "batched") {
        const Eigen::MatrixXd& effect_inverse = options_.ssgwas_blend_a22 ? Ginv_bs_tuned_ : Ginv_bs_;
        ExactSnpUncertaintyScanInput input{
            *plink_reader_, snp_effect_allele_, snp_keep_indices_, geno_full_idx_,
            G_bs_, vhat_g_, &effect_inverse, nullptr, sigma2_a_backsolve_, S_2pq_,
            options_.ssgwas_pev_scan_block_size};
        acceptExactScanResult(BatchedExactSnpUncertaintyScanner().scan(input));
        std::cout << "SNP p-values (exact, batched selected covariance): " << pev_scan_variants_tested_ << " SNPs\n";
        return;
    }

    // --- Per-SNP scan: Var(EA_j) = B_j * Var(u_hat_g) * B_j', B_j = M_j' G_backsolve^-1 / S ---
    t0 = now_t();
    const Eigen::MatrixXd& ginv_eff =
        options_.ssgwas_blend_a22 ? Ginv_bs_tuned_ : Ginv_bs_;
    const auto& bim = plink_reader_->getBimInfo();
    const auto& inv_stds = plink_reader_->getInvStds();
    const int m = static_cast<int>(bim.size());
    const double S = S_2pq_;
    snp_se_allele_.assign(m, 0.0);
    snp_z_.assign(m, 0.0);
    snp_p_.assign(m, 1.0);
    snp_log10p_.assign(m, 0.0);
    snp_pvalue_status_.assign(m, "NA");
    n_pec_clamped_ = 0;
    n_pec_failure_ = 0;
    n_zero_se_ = 0;
    Eigen::VectorXd mcol(ng), w(ng), t1(ng), t2(ng);
    Eigen::VectorXd e_snp = Eigen::VectorXd::Zero(m);
    exact_scan_actual_mode_ = "SCALAR_ONEHOT_REFERENCE";
    pev_scan_effective_block_size_ = 1;
    pev_scan_allocations_ = 5;
    pev_scan_peak_transient_bytes_ = sizeof(double) * (m + (5LL + options_.n_threads) * ng);
    for (int j : snp_keep_indices_) {
        double inv_std = j < (int)inv_stds.size() && inv_stds[j] > 0.0 ? inv_stds[j] : 0.0;
        if (inv_std <= 0.0 || !(S > 1e-12)) {
            snp_pvalue_status_[j] = "NA_MONOMORPHIC";
            continue;
        }
        // Standardized SNP column -> M = x - mean column (same standardization
        // basis as the backsolve Z'q step, so B_j is consistent with EA_j).
        e_snp.setZero();
        e_snp(j) = 1.0;
        Eigen::VectorXd zcol;
        plink_reader_->multiply_Z_v(e_snp, zcol, geno_full_idx_,
                                     options_.n_threads > 0 ? options_.n_threads : 1);
        ++pev_scan_multiply_z_v_calls_;
        ++pev_scan_genotype_passes_;
        ++pev_scan_variants_decoded_;
        ++pev_scan_variants_tested_;
        ++pev_scan_blocks_;
        ++pev_scan_allocations_; // explicit per-SNP zcol; reader/internal scratch excluded
        pev_scan_bed_bytes_read_estimate_ += static_cast<long long>(plink_reader_->getBytesPerSnp()) * m;
        mcol.noalias() = zcol / inv_std;  // a_hat_g-scale M column on genotyped animals
        w.noalias() = ginv_eff * mcol;    // G_backsolve^-1 * M_j (the B_j weight)
        t1.noalias() = vhat_g_ * w;
        t2.noalias() = G_bs_ * w;
        double var_effect = w.dot(t1) / (S * S);
        double var_scale = sigma2_a_backsolve_ * w.dot(t2) / (S * S);  // positive-scale for tolerance
        double ea = snp_effect_allele_(j);
        double tol_j = 1e-10 * std::max(1.0, std::fabs(var_scale));
        if (var_effect < 0.0) {
            if (var_effect >= -tol_j) {
                var_effect = 0.0;
                ++n_pec_clamped_;
                snp_pvalue_status_[j] = "NEGATIVE_VARIANCE_CLAMPED";
            } else {
                ++n_pec_failure_;
                snp_pvalue_status_[j] = "NUMERICAL_FAILURE";
                snp_se_allele_[j] = NAN;
                snp_z_[j] = NAN;
                snp_p_[j] = NAN;
                snp_log10p_[j] = NAN;
                continue;
            }
        }
        double se = std::sqrt(var_effect);
        snp_se_allele_[j] = se;
        double z = ea / se;
        snp_z_[j] = z;
        double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
        if (!(p > 0.0)) p = 0.0;  // underflow to 0 for |z| > ~38
        snp_p_[j] = p;
        snp_log10p_[j] = p > 0.0 ? -std::log10(p) : 320.0;
        if (se <= 1e-15 * std::max(1.0, std::fabs(ea))) {
            ++n_zero_se_;
            snp_pvalue_status_[j] = "ZERO_SE";
        } else if (snp_pvalue_status_[j] == "NA") {
            snp_pvalue_status_[j] = "OK";
        }
    }
    t_pev_scan_ = secs(t0);
    std::cout << "SNP p-values (exact, Aguilar 2019): " << snp_keep_indices_.size()
              << " SNPs, clamped=" << n_pec_clamped_ << " failures=" << n_pec_failure_
              << " zero_se=" << n_zero_se_ << "\n";
}

// ---------------------------------------------------------------------------
// Phase C2: APY core partition (Leite et al. 2024, GSE 56:59).
// Sources: (a) --ssgwas-apy-core-file (pedigree IDs, one per line; the file
// defines the core), (b) deterministic seeded random selection of
// --ssgwas-apy-core-size animals.  core_indices_/noncore_indices_ are stored
// in ascending genotyped order (the APY partition is order-insensitive).
// ---------------------------------------------------------------------------
void SingleStepModel::buildApyCoreSelection() {
    const int ng = n_genotyped_;
    if (ng < 1) {
        throw std::runtime_error("APY core partition requires genotyped animals (n_g >= 1)");
    }
    if (!options_.ssgwas_apy_core_file.empty()) {
        std::ifstream ifs(options_.ssgwas_apy_core_file);
        if (!ifs.is_open()) {
            throw std::runtime_error("Cannot open --ssgwas-apy-core-file '"
                                     + options_.ssgwas_apy_core_file + "'");
        }
        std::unordered_map<std::string, int> id_to_geno;
        for (int i = 0; i < ng; ++i) id_to_geno[genotyped_ids_[i]] = i;
        std::vector<int> core;
        std::unordered_set<int> seen;
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty() || line[0] == '#') continue;
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            if (line.empty()) continue;
            auto it = id_to_geno.find(line);
            if (it == id_to_geno.end()) {
                throw std::runtime_error("--ssgwas-apy-core-file contains an animal without a "
                                         "genotype in this run: '" + line + "'");
            }
            if (!seen.insert(it->second).second) {
                throw std::runtime_error("--ssgwas-apy-core-file contains duplicate IDs: '" + line + "'");
            }
            core.push_back(it->second);
        }
        if (core.empty()) {
            throw std::runtime_error("--ssgwas-apy-core-file defines no core animals");
        }
        if (options_.ssgwas_apy_core_size > 0 && (int)core.size() != options_.ssgwas_apy_core_size) {
            throw std::runtime_error("--ssgwas-apy-core-file has " + std::to_string(core.size()) +
                                     " IDs but --ssgwas-apy-core-size is " +
                                     std::to_string(options_.ssgwas_apy_core_size) +
                                     " (make them consistent)");
        }
        core_indices_ = core;
        std::sort(core_indices_.begin(), core_indices_.end());
        apy_core_selection_mode_ = "FILE";
    } else {
        const int n_core = options_.ssgwas_apy_core_size;
        if (n_core > ng) {
            throw std::runtime_error("--ssgwas-apy-core-size " + std::to_string(n_core) +
                                     " exceeds n_genotyped=" + std::to_string(ng));
        }
        // Deterministic seeded random core (no replacement); ascending order for
        // canonical storage (the APY partition is order-insensitive).
        std::mt19937_64 rng(static_cast<uint64_t>(options_.ssgwas_apy_seed));
        std::vector<int> perm(ng);
        std::iota(perm.begin(), perm.end(), 0);
        for (int i = 0; i < ng; ++i) {
            int j = i + (int)(rng() % (uint64_t)(ng - i));
            std::swap(perm[i], perm[j]);
        }
        core_indices_.assign(perm.begin(), perm.begin() + n_core);
        std::sort(core_indices_.begin(), core_indices_.end());
        apy_core_selection_mode_ = "RANDOM_SEEDED";
    }
    noncore_indices_.clear();
    noncore_indices_.reserve(ng - (int)core_indices_.size());
    {
        size_t ci = 0;
        for (int i = 0; i < ng; ++i) {
            while (ci < core_indices_.size() && core_indices_[ci] < i) ++ci;
            if (ci < core_indices_.size() && core_indices_[ci] == i) continue;
            noncore_indices_.push_back(i);
        }
    }
    apy_n_core_ = (int)core_indices_.size();
    apy_n_noncore_ = (int)noncore_indices_.size();
}

// ---------------------------------------------------------------------------
// Phase C2: exact SNP P-values on the APY route (--ssgwas-pvalue-method
// apy-exact).  Reference: Leite et al. (2024) GSE 56:59 Eqs. 3/7 with the C1
// terminology (docs/SSGWAS_APY_PVALUE_METHOD.md §4.3.1):
//   - MME genotyped block = tau*G_APY^-1 - omega*A22^-1   (ApyGinvAdapter)
//   - Var(u_hat_core) = (Gw)_cc*sigma_a2 - sigma_e2*(C_lambda^-1)_cc
//   - w_c = (G_backsolve^-1)_cc M_cj  (core-only weight; Eq. 7 / Bermann 2022a
//     Eq. 21-22: under the APY model the SNP effect and its sampling variance
//     depend only on the core animals)
//   - z_j = EA_j/SE_j, P = 2(1-Phi(|z|)) (Eq. 3, standard normal; identical
//     formulas to the C1 exact mode)
// Gate C2-G2: n_core = n_genotyped must reproduce the C1 exact path (G_APY = G,
// core slice = full Cinv_gg).
// ---------------------------------------------------------------------------
void SingleStepModel::computeSnpPvaluesApyExact() {
    auto now_t = std::chrono::steady_clock::now;
    auto secs = [](const std::chrono::steady_clock::time_point& t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    const int nc = apy_n_core_;
    const int ng = n_genotyped_;
    const int dim = (int)mme_lhs_.rows();
    const bool core_v2 = options_.ssgwas_apy_exact_estimator == "core-v2";
    if (!apy_ginv_mme_ || !apy_inverse_eff_operator_ || nc <= 0 || dim <= 0) {
        throw std::runtime_error("apy-exact SNP p-values require a built APY partition and a solved MME");
    }

    // --- Memory estimate + guard (n_core x n_core core-slice selected block) ---
    estimated_cinv_gg_bytes_ = 8.0 * nc * nc;
    const double kMaxExactBytes = 8ull * 1024ull * 1024ull * 1024ull;  // 8 GiB
    std::cout << "apy-exact ssGWAS SNP p-values: n_genotyped=" << ng
              << " core=" << nc << " (noncore=" << apy_n_noncore_ << "); "
              << "core-slice selected-inverse block = 8*n_core^2 = "
              << estimated_cinv_gg_bytes_ / (1024.0 * 1024.0) << " MiB\n";
    if (estimated_cinv_gg_bytes_ > kMaxExactBytes && !options_.allow_large_exact_ssgwas_pvalue) {
        throw std::runtime_error(
            "apy-exact ssGWAS SNP p-values would need > 8 GiB for the core-slice "
            "selected-inverse block (n_core=" + std::to_string(nc) + "); pass "
            "--allow-large-exact-ssgwas-pvalue to override (or increase core coverage "
            "carefully; the Phase D approximate route is the scalable production path).");
    }

    auto t0 = now_t();

    // --- Core-slice selected inverse: (C_lambda^-1)_cc via n_core sparse solves ---
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
    ldlt.compute(mme_lhs_);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("ssGBLUP MME SimplicialLDLT factorization failed (PEV step, apy-exact)");
    }
    ++dense_uncertainty_allocs_;  // debug allocation counter (gate D2)
    cinv_cc_.resize(nc, nc);
    Eigen::VectorXd e = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd x(dim);
    for (int j = 0; j < nc; ++j) {
        e(genotyped_indices_[core_indices_[j]] + 1) = 1.0;  // column 0 = intercept
        x = ldlt.solve(e);
        e(genotyped_indices_[core_indices_[j]] + 1) = 0.0;
        for (int k = 0; k < nc; ++k) cinv_cc_(k, j) = x(genotyped_indices_[core_indices_[k]] + 1);
    }
    double max_asym = 0.0;
    for (int j = 0; j < nc; ++j)
        for (int k = 0; k < nc; ++k) {
            double a = std::fabs(cinv_cc_(k, j) - cinv_cc_(j, k));
            if (a > max_asym) max_asym = a;
        }
    cinv_cc_.noalias() = 0.5 * (cinv_cc_ + cinv_cc_.transpose());  // symmetrize
    cinv_cc_max_asym_ = format_sci(max_asym, 6);
    t_pev_inverse_ = secs(t0);

    // --- Var(u_hat_core) = (Gw)_cc * sigma_a2 - sigma_e2 * Cinv_cc ---
    t0 = now_t();
    vhat_c_.resize(nc, nc);
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            vhat_c_(i, j) = sigma2_a_backsolve_ * gwCc(i, j) -
                            sigma2_e_backsolve_ * cinv_cc_(i, j);
        }
    }
    t_pev_vhat_ = secs(t0);

    // --- PSD diagnostics of Var(u_hat_core) (same convention as the C1 exact
    // mode: RELATIVE tolerance; float32 GRM eps32 ~ 1.2e-7 safety factor 10) ---
    double max_abs_eig = vhat_c_.diagonal().cwiseAbs().maxCoeff();
    double min_eig = 0.0, max_eig = 1.0;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
    if (nc <= 4096) {
        eig.compute(vhat_c_);
        if (eig.info() == Eigen::Success) {
            min_eig = eig.eigenvalues()(0);
            max_eig = eig.eigenvalues()(nc - 1);
            max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
            var_ghat_eig_mode_ = "EXPLICIT_EIGEN";
        }
    }
    if (var_ghat_eig_mode_.empty()) {
        min_eig = vhat_c_.diagonal().minCoeff();
        max_eig = vhat_c_.diagonal().maxCoeff();
        for (int r = 0; r < 24; ++r) {
            Eigen::VectorXd v = Eigen::VectorXd::Random(nc).normalized();
            double ra = v.dot(vhat_c_ * v);
            if (ra < min_eig) min_eig = ra;
            if (ra > max_eig) max_eig = ra;
        }
        max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
        var_ghat_eig_mode_ = "RAYLEIGH_PROBE";
    }
    const double psd_tol = 1e-6 * std::max(1.0, max_abs_eig);
    var_ghat_min_eig_ = format_sci(min_eig, 6);
    var_ghat_max_eig_ = format_sci(max_eig, 6);
    std::cout << "Var(u_hat_core) = (Gw)_cc*sigma_a2 - sigma_e2*Cinv_cc: min_eig=" << min_eig
              << " max_eig=" << max_eig << " psd_tol=" << psd_tol
              << " (" << var_ghat_eig_mode_ << ")\n";
    if (min_eig < -psd_tol) {
        throw std::runtime_error(
            "apy-exact SNP p-values: Var(u_hat_core) residual covariance has min eigenvalue "
            + format_sci(min_eig, 6) + " < -tol " + format_sci(psd_tol, 6) +
            "; INVALID_COVARIANCE - refusing to produce SNP P-values (check variance "
            "components / tuning / blending / core-partition consistency).");
    }

    if (core_v2 && options_.ssgwas_pev_scan_mode_explicit && options_.ssgwas_pev_scan_mode == "batched") {
        std::vector<int> core_keep(nc);
        Eigen::MatrixXd prior(nc, nc);
        for (int i = 0; i < nc; ++i) {
            core_keep[i] = geno_full_idx_[core_indices_[i]];
            for (int j = 0; j < nc; ++j) prior(i, j) = gwCc(i, j);
        }
        ExactSnpUncertaintyScanInput input{
            *plink_reader_, snp_effect_allele_, snp_keep_indices_, core_keep,
            prior, vhat_c_, nullptr, apy_inverse_eff_operator_.get(), sigma2_a_backsolve_, S_2pq_,
            options_.ssgwas_pev_scan_block_size};
        acceptExactScanResult(BatchedExactSnpUncertaintyScanner().scan(input));
        std::cout << "SNP p-values (APY-exact core-v2, batched selected-core covariance): " << pev_scan_variants_tested_ << " SNPs\n";
        return;
    }

    // --- Per-SNP scan: Var(EA_j) = w_c' Var(u_hat_core) w_c / S^2,
    //     w_c = (G_backsolve^-1)_cc M_cj (core rows only; Eq. 7) ---
    t0 = now_t();
    const auto& bim = plink_reader_->getBimInfo();
    const auto& inv_stds = plink_reader_->getInvStds();
    const int m = (int)bim.size();
    const double S = S_2pq_;
    snp_se_allele_.assign(m, 0.0);
    snp_z_.assign(m, 0.0);
    snp_p_.assign(m, 1.0);
    snp_log10p_.assign(m, 0.0);
    snp_pvalue_status_.assign(m, "NA");
    n_pec_clamped_ = 0;
    n_pec_failure_ = 0;
    n_zero_se_ = 0;
    Eigen::VectorXd zcol(ng), mcol_c(nc), w_c(nc), t1_c(nc), t2_c(nc);
    Eigen::MatrixXd gw_cc(nc, nc);
    for (int a = 0; a < nc; ++a) {
        for (int b = 0; b < nc; ++b) gw_cc(a, b) = gwCc(a, b);
    }
    Eigen::VectorXd e_snp = Eigen::VectorXd::Zero(m);
    exact_scan_actual_mode_ = "SCALAR_ONEHOT_REFERENCE";
    pev_scan_effective_block_size_ = 1;
    pev_scan_allocations_ = 7;
    pev_scan_peak_transient_bytes_ = sizeof(double) * (m + (1LL + options_.n_threads) * ng + 5LL * nc + 1LL * nc * nc);
    for (int j : snp_keep_indices_) {
        double inv_std = j < (int)inv_stds.size() && inv_stds[j] > 0.0 ? inv_stds[j] : 0.0;
        if (inv_std <= 0.0 || !(S > 1e-12)) {
            snp_pvalue_status_[j] = "NA_MONOMORPHIC";
            continue;
        }
        e_snp.setZero();
        e_snp(j) = 1.0;
        plink_reader_->multiply_Z_v(e_snp, zcol, geno_full_idx_,
                                     options_.n_threads > 0 ? options_.n_threads : 1);
        ++pev_scan_multiply_z_v_calls_;
        ++pev_scan_genotype_passes_;
        ++pev_scan_variants_decoded_;
        ++pev_scan_variants_tested_;
        ++pev_scan_blocks_;
        ++pev_scan_allocations_;
        pev_scan_bed_bytes_read_estimate_ += static_cast<long long>(plink_reader_->getBytesPerSnp()) * m;
        for (int i = 0; i < nc; ++i) {
            mcol_c(i) = zcol(core_indices_[i]) / inv_std;  // core rows of M_j
        }
        w_c.noalias() = apy_inverse_eff_operator_->applyCoreInverse(mcol_c);  // (G_backsolve^-1)_cc M_cj
        t1_c.noalias() = vhat_c_ * w_c;
        t2_c.noalias() = gw_cc * w_c;  // positive-scale tolerance reference (Gw)_cc w_c
        double var_effect = w_c.dot(t1_c) / (S * S);
        double var_scale = sigma2_a_backsolve_ * w_c.dot(t2_c) / (S * S);
        double ea = snp_effect_allele_(j);
        double tol_j = 1e-10 * std::max(1.0, std::fabs(var_scale));
        if (var_effect < 0.0) {
            if (var_effect >= -tol_j) {
                var_effect = 0.0;
                ++n_pec_clamped_;
                snp_pvalue_status_[j] = "NEGATIVE_VARIANCE_CLAMPED";
            } else {
                ++n_pec_failure_;
                snp_pvalue_status_[j] = "NUMERICAL_FAILURE";
                snp_se_allele_[j] = NAN;
                snp_z_[j] = NAN;
                snp_p_[j] = NAN;
                snp_log10p_[j] = NAN;
                continue;
            }
        }
        double se = std::sqrt(var_effect);
        snp_se_allele_[j] = se;
        if (core_v2 && se <= 1e-15 * std::max(1.0, std::fabs(ea))) {
            // A projected locus can be uninformative even when globally
            // polymorphic. Preserve the NA output policy without evaluating
            // an undefined 0/0 statistic or manufacturing a P value.
            ++n_zero_se_;
            snp_pvalue_status_[j] = "ZERO_SE";
            snp_z_[j] = snp_p_[j] = snp_log10p_[j] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        double z = ea / se;
        snp_z_[j] = z;
        double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
        if (!(p > 0.0)) p = 0.0;  // underflow to 0 for |z| > ~38
        snp_p_[j] = p;
        snp_log10p_[j] = p > 0.0 ? -std::log10(p) : 320.0;
        if (se <= 1e-15 * std::max(1.0, std::fabs(ea))) {
            ++n_zero_se_;
            snp_pvalue_status_[j] = "ZERO_SE";
        } else if (snp_pvalue_status_[j] == "NA") {
            snp_pvalue_status_[j] = "OK";
        }
    }
    t_pev_scan_ = secs(t0);
    std::cout << "SNP p-values (apy-exact, Leite 2024 Eq. 7): " << snp_keep_indices_.size()
              << " SNPs, clamped=" << n_pec_clamped_ << " failures=" << n_pec_failure_
              << " zero_se=" << n_zero_se_ << "\n";
}

// ---------------------------------------------------------------------------
// Phase D: frequentist SNP P-values for back-solved effects, Approx_GinvAPY
// route (--ssgwas-pvalue-method apy-approx).  Reference: Leite et al. 2024
// Eqs. 3/8/9; single Schur/Woodbury production form (equivalence gate D1.5):
//
//   C_core_approx^-1 = D_c + lambda*Gcc^-1 + Q*S*Q'
//       D (Eq. 9)     = diag(W'(I-X(X'X)^-1 X')W)   (per genotyped animal)
//       S_j           = lambda*d_j/(lambda + d_j*m_jj), m_jj = Mnn diagonal
//   Var(u_hat_core)   = (Gw)_cc*sigma_a2 - sigma_e2*C_core_approx   (n_c^2)
//   Var(EA_j)         = w_c' Var(u_hat_core) w_c / S^2
//       w_c           = (G_backsolve^-1)_cc M_cj  (same operator as the effect)
//
// Hard guarantees: no dense n_g^2 uncertainty block (DENSE_NG2_ALLOCATED=false,
// dense_uncertainty_allocs_ must be 0 here) and no MME selected inverse.
// ---------------------------------------------------------------------------
void SingleStepModel::computeSnpPvaluesApyApprox() {
    auto now_t = std::chrono::steady_clock::now;
    auto secs = [](const std::chrono::steady_clock::time_point& t0) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    const int ng = n_genotyped_;
    const int nc = apy_n_core_;
    const int nn = apy_n_noncore_;
    if (ng <= 0 || nc <= 0 || !apy_inverse_eff_operator_) {
        throw std::runtime_error(
            "apy-approx SNP p-values require an APY partition and a solved MME");
    }
    if (d_genotyped_.size() != (size_t)ng) {
        throw std::runtime_error(
            "apy-approx SNP p-values: effective-record diagonals (Eq. 9) missing");
    }
    // Hard allocation gate (D2): the approximate route must never touch the
    // exact/apy-exact selected-inverse state.
    if (dense_uncertainty_allocs_ != 0 || cinv_gg_.size() != 0 || vhat_g_.size() != 0 ||
        cinv_cc_.size() != 0 || vhat_c_.size() != 0) {
        throw std::runtime_error(
            "apy-approx internal gate: a dense n_g^2 uncertainty block was "
            "allocated by the exact/apy-exact path before the approximate "
            "route ran; DENSE_NG2_ALLOCATED must stay false.");
    }

    // --- Eq. 8 core PEV approximation (Schur form) ---
    const double lambda = sigma2_e_backsolve_ / sigma2_a_backsolve_;
    Eigen::VectorXd Dc(nc), Dn(nn);
    for (int i = 0; i < nc; ++i) Dc(i) = d_genotyped_[core_indices_[i]];
    for (int j = 0; j < nn; ++j) Dn(j) = d_genotyped_[noncore_indices_[j]];

    auto apply_gw_cc = [this, nc](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
        if (x.size() != nc) {
            throw std::runtime_error("apy-approx Gw_cc application size mismatch");
        }
        if (grm_streamed_) {
            y.noalias() = gcc_stream_final_ * x;
            return;
        }
        y = Eigen::VectorXd::Zero(nc);
        for (int i = 0; i < nc; ++i) {
            for (int j = 0; j < nc; ++j) {
                y(i) += G_bs_(core_indices_[i], core_indices_[j]) * x(j);
            }
        }
    };
    Eigen::MatrixXd gw_cc(nc, nc);
    for (int a = 0; a < nc; ++a) {
        for (int b = 0; b < nc; ++b) gw_cc(a, b) = gwCc(a, b);
    }

    t_pev_inverse_ = 0.0;
    auto t0 = now_t();
    // E3-D2D no-regen gate: capture baseline after all operator construction,
    // before any per-SNP applies. This allows one-time factory use during
    // MME/PEV operator construction but forbids regen during per-SNP applies.
    gcn_regen_after_construct_ = grm_stream_builder_ ? grm_stream_builder_->gcnRegenCalls() : 0;
    std::unique_ptr<APYCorePEVOperator> core_pev;
    pec_c_appx_.resize(0, 0);
    vhat_c_appx_.resize(0, 0);
    apy_pev_max_relative_residual_ = 0.0;
    if (options_.ssgwas_pev_apy_mode == "dense") {
        if (!apy_ginv_mme_) {
            throw std::runtime_error("apy-approx dense PEV requires materialized APY MME coefficients");
        }
        auto dense_pev = std::make_unique<DenseAPYCorePEVOperator>(
            *apy_ginv_mme_, Dc, Dn, lambda);
        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(nc, nc);
        pec_c_appx_ = dense_pev->apply(identity);
        core_pev = std::move(dense_pev);
    } else {
        GccSolverPtr gcc_solver;
        const Eigen::VectorXd* mnn = nullptr;
        if (apy_inverse_mme_operator_) {
            const auto* mme_operator = dynamic_cast<const MatrixFreeApyInverseOperator*>(
                apy_inverse_mme_operator_.get());
            if (!mme_operator) {
                throw std::runtime_error("apy-approx operator PEV requires matrix-free APY MME state");
            }
            gcc_solver = mme_operator->gccSolverPtr();
            mnn = &mme_operator->mnn();
        } else if (apy_ginv_mme_) {
            gcc_solver = apy_ginv_mme_->gccSolverPtr();
            mnn = &apy_ginv_mme_->mnn();
        } else {
            throw std::runtime_error("apy-approx operator PEV has no final-MME Gcc/Mnn source");
        }

        GcnBlockIteratorFactory factory;
        if (nn > 0) {
            if (grm_streamed_) {
                factory = [this]() {
                    return grm_stream_builder_->makeGcnIterator(-1, GcnBlockView::Final);
                };
            } else {
                factory = makeDenseGcnBlockIteratorFactory(
                    G_bs_, core_indices_, noncore_indices_);
            }
        }
        MatrixFreePEVOperator::Options pev_options;
        pev_options.tolerance = std::min(options_.ssgwas_pcg_tol, 1e-10);
        pev_options.max_iterations = options_.ssgwas_pcg_maxit;
        auto operator_pev = std::make_unique<MatrixFreeAPYCorePEVOperator>(
            std::move(gcc_solver), std::move(factory), *mnn, Dc, Dn, lambda,
            pev_options, gcn_cross_final_.get());
        // E3-D3: precompute the nc x nc core PEV inverse ONCE (nc column solves)
        // instead of running a matrix-free PCG solve per SNP.  The dense PEV
        // branch below materializes the identical object (Eq. 8); materializing
        // the small nc x nc operator here keeps every per-SNP p-value a cheap
        // 25 x 25 matvec while preserving the operator equations verbatim.
        Eigen::MatrixXd identity_core = Eigen::MatrixXd::Identity(nc, nc);
        pec_c_appx_ = operator_pev->apply(identity_core);
        for (const auto& info : operator_pev->lastBatchSolveInfo()) {
            if (info.relative_residual > apy_pev_max_relative_residual_) {
                apy_pev_max_relative_residual_ = info.relative_residual;
            }
        }
        core_pev = std::move(operator_pev);
    }
    apy_pev_backend_ = core_pev->backendName();
    t_pev_inverse_ += secs(t0);

    // --- Var(u_hat_core) approximate PSD diagnostics ---
    t0 = now_t();
    double max_abs_eig = 0.0;
    double min_eig = 0.0, max_eig = 0.0;
    var_ghat_eig_mode_.clear();
    if (options_.ssgwas_pev_apy_mode == "dense") {
        vhat_c_appx_ = sigma2_a_backsolve_ * gw_cc -
                       sigma2_e_backsolve_ * pec_c_appx_;
        max_abs_eig = vhat_c_appx_.diagonal().cwiseAbs().maxCoeff();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig;
        if (nc <= 4096) {
            eig.compute(vhat_c_appx_);
            if (eig.info() == Eigen::Success) {
                min_eig = eig.eigenvalues()(0);
                max_eig = eig.eigenvalues()(nc - 1);
                max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
                var_ghat_eig_mode_ = "EXPLICIT_EIGEN";
            }
        }
        if (var_ghat_eig_mode_.empty()) {
            min_eig = vhat_c_appx_.diagonal().minCoeff();
            max_eig = vhat_c_appx_.diagonal().maxCoeff();
            for (int r = 0; r < 24; ++r) {
                Eigen::VectorXd v = Eigen::VectorXd::Random(nc).normalized();
                double ra = v.dot(vhat_c_appx_ * v);
                if (ra < min_eig) min_eig = ra;
                if (ra > max_eig) max_eig = ra;
            }
            max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
            var_ghat_eig_mode_ = "RAYLEIGH_PROBE";
        }
    } else {
        min_eig = std::numeric_limits<double>::infinity();
        max_eig = -std::numeric_limits<double>::infinity();
        for (int r = 0; r < 24; ++r) {
            Eigen::VectorXd v(nc);
            for (int i = 0; i < nc; ++i) {
                v(i) = std::sin((r + 1.0) * (i + 1.0));
            }
            v.normalize();
            Eigen::VectorXd gw_v(nc);
            apply_gw_cc(v, gw_v);
            const double ra = sigma2_a_backsolve_ * v.dot(gw_v) -
                              sigma2_e_backsolve_ * v.dot(pec_c_appx_ * v);
            if (ra < min_eig) min_eig = ra;
            if (ra > max_eig) max_eig = ra;
        }
        max_abs_eig = std::max(std::fabs(min_eig), std::fabs(max_eig));
        var_ghat_eig_mode_ = "RAYLEIGH_PROBE_OPERATOR";
    }
    const double psd_tol = 1e-6 * std::max(1.0, max_abs_eig);
    var_ghat_min_eig_ = format_sci(min_eig, 6);
    var_ghat_max_eig_ = format_sci(max_eig, 6);
    std::cout << "Var(u_hat_core) approx (Eq. 8): min_eig=" << min_eig
              << " max_eig=" << max_eig << " psd_tol=" << psd_tol
              << " (" << var_ghat_eig_mode_ << ")\n";
    if (min_eig < -psd_tol) {
        throw std::runtime_error(
            "apy-approx SNP p-values: approximate Var(u_hat_core) has min eigenvalue "
            + format_sci(min_eig, 6) + " < -tol " + format_sci(psd_tol, 6) +
            "; INVALID_COVARIANCE - refusing to produce SNP P-values.");
    }
    t_pev_vhat_ = secs(t0);

    // --- Estimated peak memory (audit) ---
    {
        double bytes = 0.0;
        if (grm_streamed_) {
            // No dense G_bs_/G_bs_tuned_; held streamed blocks + optional A22 cache.
            if (options_.ssgwas_hinv_mode == "dense")
                bytes += 8.0 * (double)ng * (double)ng;     // A22_inv_ (persistent)
            if (options_.ssgwas_blend_a22 || options_.ssgwas_tune_g)
                bytes += 8.0 * (double)ng * (double)ng;     // A22_geno_ cache (Phase E3 removes)
            bytes += 8.0 * (double)nc * (double)nc;         // Gcc (final)
            if (grm_stream_builder_->fullGcnStored())
                bytes += 8.0 * (double)nc * (double)nn;     // Gcn
            bytes += 8.0 * (double)nn;                      // diag(Gnn)
            if (options_.ssgwas_blend_a22)
                bytes += 8.0 * (double)nc * (double)nc + 8.0 * (double)nn;
            if (options_.ssgwas_blend_a22 && grm_stream_builder_->fullGcnStored())
                bytes += 8.0 * (double)nc * (double)nn;     // tuned Gcn
        } else {
            bytes = 8.0 * (double)ng * (double)ng;  // G_bs_
            if (options_.ssgwas_blend_a22) bytes += 8.0 * (double)ng * (double)ng;  // G_bs_tuned_
            if (options_.ssgwas_hinv_mode == "dense")
                bytes += 8.0 * (double)ng * (double)ng;     // A22_inv_ (persistent)
        }
        if (options_.ssgwas_pev_apy_mode == "dense") {
            bytes += 8.0 * (double)nc * (double)nc;     // C_core approx / Vhat
            bytes += 8.0 * (double)nc * (double)nc;     // Kc / factors
            bytes += 8.0 * (double)nc * (double)nn;     // Q
        }
        bytes += 8.0 * (double)nc;                      // D vectors / misc
        estimated_peak_memory_bytes_ = bytes;
    }

    // --- Per-SNP scan: w_c = (G_backsolve^-1)_cc M_cj; qc = w_c'Vhat_apx w_c/S^2 ---
    const double S = S_2pq_;
    SnpUncertaintyScanInput scan_input{
        *plink_reader_, *apy_inverse_eff_operator_, gw_cc, pec_c_appx_,
        options_.ssgwas_pev_apy_mode == "dense" ? &vhat_c_appx_ : nullptr,
        snp_effect_allele_, snp_keep_indices_, geno_full_idx_, core_indices_,
        sigma2_a_backsolve_, sigma2_e_backsolve_, S,
        std::max(1, options_.n_threads), options_.ssgwas_pev_scan_block_size};
    std::unique_ptr<SnpUncertaintyScanner> scanner;
    if (options_.ssgwas_pev_scan_mode == "batched") {
        scanner = std::make_unique<BatchedOperatorSnpUncertaintyScanner>();
    } else {
        scanner = std::make_unique<ScalarReferenceSnpUncertaintyScanner>();
    }
    SnpUncertaintyScanResult scan_result = scanner->scan(scan_input);
    snp_se_allele_ = std::move(scan_result.se);
    snp_z_ = std::move(scan_result.z);
    snp_p_ = std::move(scan_result.p);
    snp_log10p_ = std::move(scan_result.log10p);
    snp_pvalue_status_ = std::move(scan_result.status);
    n_pec_clamped_ = scan_result.clamped;
    n_pec_failure_ = scan_result.failures;
    n_zero_se_ = scan_result.zero_se;
    const auto& audit = scan_result.audit;
    t_pev_scan_ = audit.time_total;
    t_pev_scan_setup_ = audit.time_setup;
    t_pev_scan_onehot_ = audit.time_onehot;
    t_pev_scan_genotype_ = audit.time_genotype;
    t_pev_scan_core_gather_ = audit.time_core_gather;
    t_pev_scan_gcc_solve_ = audit.time_gcc_solve;
    t_pev_scan_gw_matvec_ = audit.time_gw_multiply;
    t_pev_scan_pec_matvec_ = audit.time_pec_multiply;
    t_pev_scan_quadratic_ = audit.time_quadratic;
    t_pev_scan_statistics_ = audit.time_statistics;
    t_pev_scan_other_ = audit.time_other;
    pev_scan_multiply_z_v_calls_ = audit.multiply_z_v_calls;
    pev_scan_genotype_passes_ = audit.genotype_passes;
    pev_scan_bed_bytes_read_estimate_ = audit.bed_bytes_read_estimate;
    pev_scan_variants_tested_ = audit.variants_tested;
    pev_scan_allocations_ = audit.heap_allocations;
    pev_scan_blocks_ = audit.blocks;
    pev_scan_variants_decoded_ = audit.variants_decoded;
    pev_scan_block_reallocations_ = audit.block_reallocations;
    pev_scan_peak_transient_bytes_ = audit.peak_transient_bytes;
    pev_scan_gcc_solve_calls_ = audit.gcc_solve_calls;
    pev_scan_gcc_rhs_columns_ = audit.gcc_rhs_columns;
    pev_scan_gw_gemm_calls_ = audit.gw_gemm_calls;
    pev_scan_pec_gemm_calls_ = audit.pec_gemm_calls;
    pev_scan_missing_to_zero_ = audit.missing_to_zero;
    pev_scan_effective_block_size_ = audit.block_size;
    estimated_peak_memory_bytes_ += static_cast<double>(audit.peak_transient_bytes);
    std::cout << "SNP p-values (apy-approx, Leite 2024 Eq. 8/10): " << snp_keep_indices_.size()
              << " SNPs, clamped=" << n_pec_clamped_ << " failures=" << n_pec_failure_
              << " zero_se=" << n_zero_se_ << "\n";
    // Phase E3-D2D no-regen gate: with the matrix-free cross operator active the
    // per-SNP PCG applies must never regenerate the full K x Nn Gcn block.
    if (gcn_cross_final_ || gcn_cross_tuned_) {
        const long long regen_after_scan = grm_stream_builder_->gcnRegenCalls();
        std::cout << "Phase E3-D2D no-regen check: gcn block regenerations = "
                  << (regen_after_scan - gcn_regen_after_construct_)
                  << " during operator construction -> per-SNP apply ("
                  << gcn_regen_after_construct_ << " -> " << regen_after_scan
                  << ")\n";
        if (regen_after_scan != gcn_regen_after_construct_) {
            throw std::runtime_error(
                "Phase E3-D2D no-regen gate: Gcn block was regenerated during the "
                "matrix-free operator applies");
        }
    }
}

// ---------------------------------------------------------------------------
// Output files: <out>.ss_backsolve.{ebv,snp,pvalue,window,meta}.txt
// ---------------------------------------------------------------------------
void SingleStepModel::writePvalueFile() const {
    const std::string path = options_.out_prefix + ".ss_backsolve.pvalue.txt";
    std::ofstream ofs(path);
    ofs << std::setprecision(15);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tEFFECT_ALLELE\tSE_ALLELE\tVAR_EFFECT\tZ\tP\tLOG10P\tSTATUS\n";
    const auto& bim = plink_reader_->getBimInfo();
    const auto& means = plink_reader_->getMeans();
    for (int j : snp_keep_indices_) {
        const auto& row = bim[j];
        double mean_dos = j < (int)means.size() ? means[j] : NAN;
        double freq = (mean_dos >= 0.0 && std::isfinite(mean_dos)) ? mean_dos / 2.0 : -1.0;
        const std::string& st = snp_pvalue_status_[j];
        if (st == "NA_MONOMORPHIC" || st == "NUMERICAL_FAILURE" || st == "ZERO_SE") {
            ofs << row.chrom << "\t" << row.id << "\t" << row.bp_pos << "\t"
                << (row.alt.empty() ? "-" : row.alt) << "\t"
                << (row.ref.empty() ? "-" : row.ref) << "\t" << freq << "\t"
                << snp_effect_allele_(j) << "\tNA\tNA\tNA\tNA\tNA\t" << st << "\n";
        } else {
            ofs << row.chrom << "\t" << row.id << "\t" << row.bp_pos << "\t"
                << (row.alt.empty() ? "-" : row.alt) << "\t"
                << (row.ref.empty() ? "-" : row.ref) << "\t" << freq << "\t"
                << snp_effect_allele_(j) << "\t" << snp_se_allele_[j] << "\t"
                << snp_se_allele_[j] * snp_se_allele_[j] << "\t"
                << snp_z_[j] << "\t" << snp_p_[j] << "\t" << snp_log10p_[j] << "\t"
                << st << "\n";
        }
    }
    ofs.close();
    std::cout << "Wrote " << path << "\n";
}

void SingleStepModel::writeEbvFile() const {
    const std::string path = options_.out_prefix + ".ss_backsolve.ebv.txt";
    std::ofstream ofs(path);
    ofs << std::setprecision(15);
    ofs << "ID\tEBV\tGENOTYPED\tPHENOTYPED\n";
    const auto& ped_records = pedigree_.records();
    std::unordered_set<int> geno_rows(genotyped_indices_.begin(), genotyped_indices_.end());
    std::unordered_set<int> pheno_rows(pheno_ped_rows_.begin(), pheno_ped_rows_.end());
    for (int i = 0; i < n_ped_; ++i) {
        ofs << ped_records[i].id << "\t" << a_hat_(i) << "\t"
            << (geno_rows.count(i) ? 1 : 0) << "\t" << (pheno_rows.count(i) ? 1 : 0) << "\n";
    }
    ofs.close();
    std::cout << "Wrote " << path << "\n";
}

void SingleStepModel::writeSnpFile() const {
    const std::string path = options_.out_prefix + ".ss_backsolve.snp.txt";
    std::ofstream ofs(path);
    ofs << std::setprecision(15);
    ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tEFFECT_ALLELE\tSNP_EFFECT_STD\tSNP_VARIANCE\tPCT_GENETIC_VARIANCE\n";
    const auto& bim = plink_reader_->getBimInfo();
    const auto& means = plink_reader_->getMeans();
    for (int j : snp_keep_indices_) {
        const auto& row = bim[j];
        double pct = total_snp_variance_ > 0.0 ? snp_variance_(j) / total_snp_variance_ : 0.0;
        double mean_dos = j < (int)means.size() ? means[j] : NAN;
        double freq = (mean_dos >= 0.0 && std::isfinite(mean_dos)) ? mean_dos / 2.0 : -1.0;
        ofs << row.chrom << "\t" << row.id << "\t" << row.bp_pos << "\t"
            << (row.alt.empty() ? "-" : row.alt) << "\t"
            << (row.ref.empty() ? "-" : row.ref) << "\t"
            << freq << "\t" << snp_effect_allele_(j) << "\t" << snp_effect_std_(j)
            << "\t" << snp_variance_(j) << "\t" << pct << "\n";
    }
    ofs.close();
    std::cout << "Wrote " << path << "\n";
}

void SingleStepModel::writeWindowFile() const {
    const int W = options_.ssgwas_window_snps;
    const std::string path = options_.out_prefix + ".ss_backsolve.window.txt";
    std::ofstream ofs(path);
    ofs << std::setprecision(15);
    ofs << "CHR\tWINDOW_START_BP\tWINDOW_END_BP\tFIRST_SNP\tLAST_SNP\tN_SNP\tWINDOW_VARIANCE\tPCT_GENETIC_VARIANCE\n";
    const auto& bim = plink_reader_->getBimInfo();
    const int n = (int)snp_keep_indices_.size();
    // Sliding SNP-count windows; never crossing chromosome boundaries.
    int s = 0;
    while (s < n) {
        const std::string chrom = bim[snp_keep_indices_[s]].chrom;
        int e = s;
        while (e + 1 < n && bim[snp_keep_indices_[e + 1]].chrom == chrom) ++e;
        const int chrom_len = e - s + 1;
        if (chrom_len <= W) {
            double wvar = 0.0;
            for (int k = s; k <= e; ++k) wvar += snp_variance_(snp_keep_indices_[k]);
            double pct = total_snp_variance_ > 0.0 ? wvar / total_snp_variance_ : 0.0;
            ofs << chrom << "\t" << bim[snp_keep_indices_[s]].bp_pos << "\t"
                << bim[snp_keep_indices_[e]].bp_pos << "\t"
                << bim[snp_keep_indices_[s]].id << "\t" << bim[snp_keep_indices_[e]].id << "\t"
                << chrom_len << "\t" << wvar << "\t" << pct << "\n";
        } else {
            for (int k = s; k + W - 1 <= e; ++k) {
                double wvar = 0.0;
                for (int j = k; j < k + W; ++j) wvar += snp_variance_(snp_keep_indices_[j]);
                double pct = total_snp_variance_ > 0.0 ? wvar / total_snp_variance_ : 0.0;
                ofs << chrom << "\t" << bim[snp_keep_indices_[k]].bp_pos << "\t"
                    << bim[snp_keep_indices_[k + W - 1]].bp_pos << "\t"
                    << bim[snp_keep_indices_[k]].id << "\t" << bim[snp_keep_indices_[k + W - 1]].id << "\t"
                    << W << "\t" << wvar << "\t" << pct << "\n";
            }
        }
        s = e + 1;
    }
    ofs.close();
    std::cout << "Wrote " << path << "\n";
}

void SingleStepModel::writeMetaFile() const {
    std::string path = options_.out_prefix + ".ss_backsolve.meta.txt";
    std::ofstream ofs(path);
    ofs << "MODEL\tssGBLUP-derived back-solving ssGWAS (Method 2)\n"
        << "GENETIC_EVALUATION\tfull-animal MME  y = Xb + Za + e (X = intercept in v1),  a ~ N(0, H*Va)\n"
        << "HINV\tA^-1 + [tau*Ginv - omega*A22inv] on genotyped block ("
        << (options_.ssgwas_hinv_mode == "operator" ? "operator-only" : "dense-reference")
        << " MME path)\n"
        << "INPUT_FORMAT\tPLINK .bed/.bim/.fam only in v1 (BGEN/PGEN rejected)\n"
        << "TAU\t" << options_.h_tau << "\n"
        << "OMEGA\t" << options_.h_omega << "\n"
        << "G_METHOD\tVanRaden (M M' / S)\n"
        << "G_RAW_METHOD\tVanRaden M M' / S, genotype content centered by sample mean\n"
        << "G_RAW_DENOMINATOR\tS = sum_j 2p_j(1-p_j) over the full genotyped sample set (S_SNP)\n"
        << "G_MEAN_DIAG_BEFORE\t" << g_raw_mean_diag_ << "\n"
        << "G_MEAN_OFFDIAG_BEFORE\t" << g_raw_mean_off_ << "\n"
        << "A22_MEAN_DIAG\t" << a22_mean_diag_ << "\n"
        << "A22_MEAN_OFFDIAG\t" << a22_mean_off_ << "\n"
        << "G_TUNED\t" << (options_.ssgwas_tune_g ? "true" : "false") << "\n"
        << "G_SCALING_A\t" << tune_a_ << "\n"
        << "G_SCALING_B\t" << tune_b_ << "\n"
        << "G_MEAN_DIAG_AFTER\t" << g_final_mean_diag_ << "\n"
        << "G_MEAN_OFFDIAG_AFTER\t" << g_final_mean_off_ << "\n"
        << "BLEND_TARGET\t" << blend_target_name_ << "\n"
        << "BLEND_ALPHA\t" << options_.blend_g << "\n"
        << "G_MATERIALIZATION\t"
        << (grm_streamed_ ? (grm_stream_builder_->fullGcnStored()
                                 ? "BLOCK_STREAMED_GCC_GCN_GNN_DIAG"
                                 : "BLOCK_STREAMED_GCC_GNN_DIAG_GENOTYPE_BACKED_GCN")
                          : "DENSE_FULL_NG2") << "\n"
        << "FULL_G_MATERIALIZED\t" << (grm_streamed_ ? "false" : "true") << "\n"
        << "FULL_GNN_MATERIALIZED\t" << (grm_streamed_ ? "false" : "true") << "\n"
        << "FULL_GCN_STORED\t"
        << (grm_streamed_ ? (grm_stream_builder_->fullGcnStored() ? "true" : "false") : "NA") << "\n"
        << "EFFECT_GINV\t" << effect_ginv_name_ << "\n"
        << "RIDGE\t1e-8 geometric growth (invertSpdWithRidge)\n"
        << "ALLELE_FREQUENCY_SOURCE\tgenotyped-sample mean dosage (reader snp_means)\n"
        << "GENOTYPE_CENTERING\tx - mean (M matrix of VanRaden G)\n"
        << "GENOTYPE_SCALING\tnone for G; multiply_Zt_v uses (x-mean)/sd internally, undone by 1/inv_std\n"
        << "BACKSOLVE_FORMULA\teffect_allele = (Z' q)/(std * S); effect_std = effect_allele * std; q = G^-1 a_hat_g; S = sum 2p(1-p)\n"
        << "EFFECT_SCALE\tper allele dosage (M = x - mean centering)\n"
        << "EFFECT_SCALE_STD\tper standardized-genotype SD: EFFECT_ALLELE * sd (SNP_EFFECT_STD); Z * SNP_EFFECT_STD = M * EFFECT_ALLELE = GEBV on genotyped block\n"
        << "SNP_VARIANCE_DEFINITION\t2p_j(1-p_j) * EFFECT_ALLELE_j^2 (per-marker marginal additive variance; the per-SNP-entry interpretation is the sample variance Var(M_j * EA_j) with consistency factor (n-1)/n)\n"
        << "WINDOW_VARIANCE_DEFINITION\tsum over window SNPs of 2p_j(1-p_j)*effect_allele_j^2 (Definition A, marker-marginal variances summed; NOT Var(M_win * alpha_win), which differs under LD)\n"
        << "PCT_VARIANCE_DENOMINATOR\tsum_j 2p_j(1-p_j)*effect_allele_j^2 (total SNP variance); PCT is per-marker share of this sum, not of additive genetic variance sigma2_a\n"
        << "WINDOW_MODE\tsliding SNP-count, never crossing chromosome\n"
        << "WINDOW_SIZE\t" << options_.ssgwas_window_snps << "\n"
        << "VARIANCE_DENOMINATOR\tsum_j 2p_j(1-p_j) * effect_allele_j^2 (total SNP variance)\n"
        << "VCE_MODE\t" << vce_mode_used_ << "\n"
        << "SIGMA2_A\t" << sigma2_a_backsolve_ << "\n"
        << "SIGMA2_E\t" << sigma2_e_backsolve_ << "\n"
        << "SOLVER\t" << options_.ssgwas_solver << "\n"
        << "HINV_MODE\t" << options_.ssgwas_hinv_mode << "\n"
        << "HINV_MATERIALIZED\tfalse\n"
        << "FULL_HINV_STORED\tfalse\n"
        << "A22_MATERIALIZED\t" << (A22_geno_.size() > 0 ? "true" : "false") << "\n"
        << "FULL_A22_STORED\t" << (A22_geno_.size() > 0 ? "true" : "false") << "\n"
        << "A22_INV_MATERIALIZED\t" << (options_.ssgwas_hinv_mode == "dense" ? "true" : "false") << "\n"
        << "FULL_A22_INV_STORED\t" << (options_.ssgwas_hinv_mode == "dense" ? "true" : "false") << "\n"
        << "A22_INV_BACKEND\t" << (a22_inverse_operator_ ? a22_inverse_operator_->backendName() : "dense-reference") << "\n"
        << "A22_BLEND_MODE\t" << (options_.ssgwas_blend_a22 ? "dense" : "not_needed") << "\n"
        << "HINV_PRECONDITIONER_DIAGONAL\t" << hinv_preconditioner_diagonal_mode_ << "\n"
        << "A22_OPERATOR_N_UNTYPED\t" << (a22_inverse_operator_ ? a22_inverse_operator_->nUntyped() : 0) << "\n"
        << "A22_OPERATOR_P_NNZ\t" << (a22_inverse_operator_ ? a22_inverse_operator_->pNonZeros() : 0) << "\n"
        << "A22_OPERATOR_Q_NNZ\t" << (a22_inverse_operator_ ? a22_inverse_operator_->qNonZeros() : 0) << "\n"
        << "A22_OPERATOR_R_NNZ\t" << (a22_inverse_operator_ ? a22_inverse_operator_->rNonZeros() : 0) << "\n"
        << "A22_OPERATOR_FACTOR_NNZ\t" << (a22_inverse_operator_ ? a22_inverse_operator_->factorNonZeros() : 0) << "\n"
        << "A22_OPERATOR_FACTOR_BYTES_ESTIMATE\t" << (a22_inverse_operator_ ? a22_inverse_operator_->factorBytesEstimate() : 0) << "\n"
        << "A22_OPERATOR_STORAGE_BYTES_ESTIMATE\t" << (a22_inverse_operator_ ? a22_inverse_operator_->storageBytesEstimate() : 0) << "\n"
        << "A22_OPERATOR_FACTOR_TIME\t" << (a22_inverse_operator_ ? a22_inverse_operator_->factorTimeSeconds() : 0.0) << "\n"
        << "A22_OPERATOR_SOLVE_COUNT\t" << (a22_inverse_operator_ ? a22_inverse_operator_->solveCount() : 0) << "\n"
        << "A22_OPERATOR_SOLVE_TIME\t" << (a22_inverse_operator_ ? a22_inverse_operator_->solveTimeSeconds() : 0.0) << "\n"
        << "A22_OPERATOR_MAX_RELATIVE_RESIDUAL\t" << (a22_inverse_operator_ ? a22_inverse_operator_->maxSolveRelativeResidual() : 0.0) << "\n"
        << "MME_APY_MODE\t" << options_.ssgwas_mme_apy_mode << "\n"
        << "MME_APY_Q_REQUIRED\t" << (options_.ssgwas_mme_apy_mode == "dense" ? "true" : "false") << "\n"
        << "MME_APY_CORE_BLOCK_REQUIRED\t" << (options_.ssgwas_mme_apy_mode == "dense" ? "true" : "false") << "\n"
        << "PEV_APY_MODE\t" << options_.ssgwas_pev_apy_mode << "\n"
        << "PEV_APY_Q_REQUIRED\t"
        << (options_.ssgwas_pvalue_method == "apy-approx" &&
                    options_.ssgwas_pev_apy_mode == "dense" ? "true" : "false") << "\n"
        << "PEV_APY_CORE_MATRIX_REQUIRED\t"
        << (options_.ssgwas_pvalue_method == "apy-approx" &&
                    options_.ssgwas_pev_apy_mode == "dense" ? "true" : "false") << "\n"
        << "SOLVER_ITERATIONS\t" << solver_iterations_ << "\n"
        << "SOLVER_RESIDUAL\t" << solver_residual_ << "\n"
        << "SOLVER_CONVERGED\t" << (solver_converged_ ? "true" : "false") << "\n"
        << "N_PEDIGREE\t" << n_ped_ << "\n"
        << "N_PHENOTYPED\t" << n_phenotyped_ << "\n"
        << "N_GENOTYPED\t" << n_genotyped_ << "\n"
        << "N_GENOTYPED_WITH_PHENO\t" << n_genotyped_with_pheno_ << "\n"
        << "N_UNGENOTYPED_WITH_PHENO\t" << n_ungenotyped_with_pheno_ << "\n"
        << "S_SNP\t" << S_2pq_ << "\n"
        << "TOTAL_SNP_VARIANCE\t" << total_snp_variance_ << "\n"
        << "TIME_PARSE\t" << t_parse_ << "\n"
        << "TIME_AINV\t" << t_ainv_ << "\n"
        << "TIME_G\t" << t_g_ << "\n"
        << "TIME_GINV\t" << t_ginv_ << "\n"
        << "TIME_HINV\t" << t_hnv_ << "\n"
        << "TIME_MME\t" << t_mme_ << "\n"
        << "TIME_SOLVE\t" << t_solve_ << "\n"
        << "TIME_BACKSOLVE\t" << t_bs_ << "\n"
        << "TIME_SNPVAR\t" << t_snp_ << "\n"
        << "TIME_WINDOW\t" << t_win_ << "\n";
    if (options_.ssgwas_pvalue_method == "exact") {
        ofs << "P_VALUE_METHOD\tAGUILAR2019_EXACT\n"
            << "TEST_STATISTIC\tZ\n"
            << "REFERENCE_DISTRIBUTION\tSTANDARD_NORMAL\n"
            << "SIDES\t2\n"
            << "REFERENCE_FORMULA\tAguilar 2019 GSE 51:28 (PMC6584984); z = effect/sd(effect), P = 2(1-Phi(|z|))\n"
            << "MME_DIMENSION\t" << mme_lhs_.rows() << "\n"
            << "MME_INVERSE_METHOD\tSELECTED_SOLVES_LDLT\n"
            << "MME_SCALE\tLAMBDA_FORM (C_lambda solves [X'X X'Z; Z'X Z'Z+lambda*Hinv], lambda = Ve/Va)\n"
            << "PEV_SCALE\tSIGMA_E2 * C_LAMBDA_INV (PEV(a) = sigma_e2 * (C_lambda^-1); variance of predictor Var(u_hat_g) = G_eval*sigma_a2 - sigma_e2*Cinv_gg)\n"
            << "EFFECT_SCALE\tALLELE_DOSAGE (EFFECT_ALLELE)\n"
            << "BACKSOLVE_OPERATOR\tB_j = M_j' * G_backsolve^-1 / S (same operator for effect and variance)\n"
            << "G_EVALUATION_CONVENTION\tblended Gw (final, as in H^-1 genotyped block; H_gg = Gw)\n"
            << "G_BACKSOLVE_CONVENTION\t" << effect_ginv_name_
            << " (tuned-unblended G*^-1 when --ssgwas-blend-a22, else blended Gw^-1)\n"
            << "TUNE_G\t" << (options_.ssgwas_tune_g ? "true" : "false") << "\n"
            << "TUNE_A\t" << tune_a_ << "\n"
            << "TUNE_B\t" << tune_b_ << "\n"
            << "BLEND_ALPHA\t" << options_.blend_g << "\n"
            << "TAU\t" << options_.h_tau << "\n"
            << "OMEGA\t" << options_.h_omega << "\n"
            << "EXACT_OR_APPROX\tEXACT\n"
            << "DENSE_NG2_ALLOCATED\ttrue\n"
            << "UPSTREAM_DENSE_NG2_PRESENT\ttrue (G/A22 are dense n_g^2; they are the dominant bottleneck, Phase E backlog)\n"
            << "K98_G\t" << k98_g_ << "\n"
            << "K98_G_CUMFRACTION\t" << k98_cumfraction_ << "\n"
            << "APY_CORE_AUTO\t" << (apy_core_auto_ ? "true" : "false") << "\n"
            << "APY_CORE_AUTO_VARIANCE\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_variance_) : "NA") << "\n"
            << "APY_CORE_AUTO_K98\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K98_) : "NA") << "\n"
            << "APY_CORE_AUTO_K99\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K99_) : "NA") << "\n"
            << "APY_CORE_EIGENSOLVER\t" << (apy_core_auto_ ? apy_core_eigensolver_used_ : "NA") << "\n"
            << "APY_CORE_TRACE_G\t" << (apy_core_auto_ ? std::to_string(apy_core_trace_g_) : "NA") << "\n"
            << "APY_CORE_EIG_ITERATIONS\t" << (apy_core_auto_ ? std::to_string(apy_core_eig_iterations_) : "NA") << "\n"
            << "APY_INVERSE_MODE\t" << options_.ssgwas_apy_inverse_mode << "\n"
            << "ESTIMATED_CINV_GG_BYTES\t" << estimated_cinv_gg_bytes_ << "\n"
            << "CINV_GG_MAX_ASYMMETRY\t" << cinv_gg_max_asym_ << "\n"
            << "VAR_GHAT_MIN_EIG\t" << var_ghat_min_eig_ << "\n"
            << "VAR_GHAT_MAX_EIG\t" << var_ghat_max_eig_ << "\n"
            << "VAR_GHAT_EIG_MODE\t" << var_ghat_eig_mode_ << "\n"
            << "PEC_NEGATIVE_CLAMPED\t" << n_pec_clamped_ << "\n"
            << "PEC_NUMERICAL_FAILURE\t" << n_pec_failure_ << "\n"
            << "ZERO_SE_COUNT\t" << n_zero_se_ << "\n"
            << "TIME_PEV_SELECTED_INVERSE\t" << t_pev_inverse_ << "\n"
            << "TIME_PEV_VHAT\t" << t_pev_vhat_ << "\n"
            << "TIME_PEV_SCAN\t" << t_pev_scan_ << "\n"
            << "TIME_PEV_OUTPUT\t" << t_pev_out_ << "\n";
        writeExactScanAudit(ofs);
    } else if (options_.ssgwas_pvalue_method == "apy-exact") {
        const bool core_v2 = options_.ssgwas_apy_exact_estimator == "core-v2";
        const std::string core_path = options_.out_prefix + ".ss_backsolve.core_ids.txt";
        std::ofstream core_ids(core_path);
        if (!core_ids) throw std::runtime_error("Cannot write APY core IDs: " + core_path);
        for (int index : core_indices_) core_ids << genotyped_ids_[index] << "\n";
        const auto previous_precision = ofs.precision();
        ofs << "APY_EXACT_ESTIMATOR\t" << options_.ssgwas_apy_exact_estimator << "\n"
            << "APY_EXACT_ESTIMATOR_EXPLICIT\t" << (options_.ssgwas_apy_exact_estimator_explicit ? "true" : "false") << "\n"
            << "SNP_POINT_PROJECTION\t" << (core_v2 ? "PROJECTED_CORE" : "RAW_FULL_APY") << "\n"
            << "SNP_COVARIANCE_PROJECTION\tCORE\n"
            << "STATISTICAL_CONSISTENCY\t" << (core_v2 ? "MATCHED_POINT_VARIANCE" : "LEGACY_UNMATCHED_POINT_VARIANCE") << "\n"
            << "LEGACY_INFERENCE_RECOMMENDED\tfalse\n"
            << "EXACT_LABEL_SCOPE\t" << (core_v2 ? "CONDITIONAL_FITTED_APY_PROJECTED_CORE" : "SELECTED_CORE_PEV_ONLY_NOT_PAIRED_TEST") << "\n"
            << "CALIBRATION_STATUS\tNOT_VALIDATED_BY_PAIRING_REPAIR\n"
            << "STRICT_TYPE_I_CONTROL\tNOT_GUARANTEED\n"
            << "KAPPA_CONVENTION\tCURRENT_COMPATIBILITY_1_OVER_S_NO_EXTRA_BLEND_FACTOR\n"
            << "KAPPA\t" << std::setprecision(17) << 1.0 / S_2pq_ << "\n"
            << "REPORTING_POLICY\tERFC_CAP320_V1\n"
            << "OUTPUT_SIGNIFICANT_DIGITS\t15\n"
            << "UNCERTAINTY_SCAN_ACTUAL_BACKEND\t" << exact_scan_actual_mode_ << "\n"
            << "APY_CORE_ACTUAL_ID_FILE\t" << core_path << "\n";
        ofs.precision(previous_precision);
        ofs << "P_VALUE_METHOD\t" << (core_v2 ? "BERMANN2022_APY_CORE_EXACT_V2" : "AGUILAR2019_APY_EXACT") << "\n"
            << "TEST_STATISTIC\tZ\n"
            << "REFERENCE_DISTRIBUTION\tSTANDARD_NORMAL\n"
            << "SIDES\t2\n"
            << "REFERENCE_FORMULA\tLeite et al. 2024 GSE 56:59 (PMC11340074) Eqs. 3/7 (APY); z = effect/sd(effect), P = 2(1-Phi(|z|)); Aguilar 2019-style per-SNP statistics\n"
            << "G_INVERSE_METHOD\tGINV_APY\n"
            << "G_INVERSE_REFERENCE\tLeite et al. 2024 Eq. 5-6 (G_APY^-1 = [Gcc^-1 + Q Mnn^-1 Q' -Q Mnn^-1; -Mnn^-1 Q' Mnn^-1], Q = Gcc^-1 Gcn, Mnn = diag(Gnn - Gnc Gcc^-1 Gcn))\n"
            << "APY_CORE_SIZE\t" << apy_n_core_ << "\n"
            << "APY_NONCORE_SIZE\t" << apy_n_noncore_ << "\n"
            << "APY_CORE_SELECTION\t" << apy_core_selection_mode_ << "\n"
            << "APY_SEED\t" << (apy_core_selection_mode_ == "RANDOM_SEEDED"
                                    ? std::to_string(options_.ssgwas_apy_seed) : "NA") << "\n"
            << "APY_CORE_FILE\t" << (apy_core_selection_mode_ == "FILE"
                                         ? options_.ssgwas_apy_core_file : "NA") << "\n"
            << "CORE_TARGET_RATIO\t" << (options_.ssgwas_apy_core_target_ratio > 0.0
                                             ? std::to_string(options_.ssgwas_apy_core_target_ratio) : "NA") << "\n"
            << "K98_G\t" << k98_g_ << "\n"
            << "K98_G_CUMFRACTION\t" << k98_cumfraction_ << "\n"
            << "APY_CORE_AUTO\t" << (apy_core_auto_ ? "true" : "false") << "\n"
            << "APY_CORE_AUTO_VARIANCE\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_variance_) : "NA") << "\n"
            << "APY_CORE_AUTO_K98\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K98_) : "NA") << "\n"
            << "APY_CORE_AUTO_K99\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K99_) : "NA") << "\n"
            << "APY_CORE_EIGENSOLVER\t" << (apy_core_auto_ ? apy_core_eigensolver_used_ : "NA") << "\n"
            << "APY_CORE_TRACE_G\t" << (apy_core_auto_ ? std::to_string(apy_core_trace_g_) : "NA") << "\n"
            << "APY_CORE_EIG_ITERATIONS\t" << (apy_core_auto_ ? std::to_string(apy_core_eig_iterations_) : "NA") << "\n"
            << "APY_GCC_MIN_LDLT_PIVOT\t" << apy_gcc_min_diag_d_ << "\n"
            << "APY_GCC_MAX_LDLT_PIVOT\t" << apy_gcc_max_diag_d_ << "\n"
            << "APY_MNN_MIN\t" << apy_mnn_min_ << "\n"
            << "APY_MNN_MAX\t" << apy_mnn_max_ << "\n"
            << "MME_DIMENSION\t" << mme_lhs_.rows() << "\n"
            << "MME_INVERSE_METHOD\tSELECTED_SOLVES_LDLT_CORE_SLICE\n"
            << "MME_SCALE\tLAMBDA_FORM (C_lambda solves [X'X X'Z; Z'X Z'Z+lambda*Hinv], lambda = Ve/Va; genotyped block = tau*G_APY^-1 - omega*A22^-1)\n"
            << "PEV_SCALE\tSIGMA_E2 * (C_LAMBDA_INV)_CC (core slice; Var(u_hat_core) = (Gw)_cc*sigma_a2 - sigma_e2*Cinv_cc)\n"
            << "EFFECT_SCALE\tALLELE_DOSAGE (EFFECT_ALLELE)\n"
            << "BACKSOLVE_OPERATOR\t"
            << (core_v2 ? "B_j = M_cj' * Gcc_effect^-1 / S; identical projected/core operator for point and exact selected-core covariance"
                        : "B_j = M_cj' * (G_backsolve^-1)_cc / S (core-only weight, Leite 2024 Eq. 7); effect q = G_APY^-1 a_hat_g (full APY product)") << "\n"
            << "G_EVALUATION_CONVENTION\tblended Gw (final, as in H^-1 genotyped block); APY-inverted\n"
            << "G_BACKSOLVE_CONVENTION\t" << effect_ginv_name_
            << " (APY inverse of tuned-unblended G* when --ssgwas-blend-a22, else of blended Gw)\n"
            << "TUNE_G\t" << (options_.ssgwas_tune_g ? "true" : "false") << "\n"
            << "TUNE_A\t" << tune_a_ << "\n"
            << "TUNE_B\t" << tune_b_ << "\n"
            << "BLEND_ALPHA\t" << options_.blend_g << "\n"
            << "TAU\t" << options_.h_tau << "\n"
            << "OMEGA\t" << options_.h_omega << "\n"
            << "GEBV_UNCERTAINTY_METHOD\tEXACT_LHS_CORE_SLICE\n"
            << "EXACT_OR_APPROX\tEXACT\n"
            << "DENSE_NG2_ALLOCATED\tfalse\n"
            << "UPSTREAM_DENSE_NG2_PRESENT\t"
            << (grm_streamed_ && options_.ssgwas_hinv_mode == "operator" && A22_geno_.size() == 0
                    ? "false"
                    : (grm_streamed_ ? "true (A22 is dense for tuning/blending)"
                                     : "true (dense G and/or A22 path)")) << "\n"
            << "ESTIMATED_CINV_CC_BYTES\t" << estimated_cinv_gg_bytes_ << "\n"
            << "CINV_CC_MAX_ASYMMETRY\t" << cinv_cc_max_asym_ << "\n"
            << "VAR_GHAT_MIN_EIG\t" << var_ghat_min_eig_ << "\n"
            << "VAR_GHAT_MAX_EIG\t" << var_ghat_max_eig_ << "\n"
            << "VAR_GHAT_EIG_MODE\t" << var_ghat_eig_mode_ << "\n"
            << "PEC_NEGATIVE_CLAMPED\t" << n_pec_clamped_ << "\n"
            << "PEC_NUMERICAL_FAILURE\t" << n_pec_failure_ << "\n"
            << "ZERO_SE_COUNT\t" << n_zero_se_ << "\n"
            << "TIME_PEV_SELECTED_INVERSE\t" << t_pev_inverse_ << "\n"
            << "TIME_PEV_VHAT\t" << t_pev_vhat_ << "\n"
            << "TIME_PEV_SCAN\t" << t_pev_scan_ << "\n"
            << "TIME_PEV_OUTPUT\t" << t_pev_out_ << "\n";
        writeExactScanAudit(ofs);
    } else if (options_.ssgwas_pvalue_method == "apy-approx") {
        ofs << "P_VALUE_METHOD\tAGUILAR2019_APY_APPROX\n"
            << "APPROXIMATION_STATUS\tEXPERIMENTAL_APPROXIMATE\n"
            << "CALIBRATION_STATUS\tLIMITED / REGIME_DEPENDENT\n"
            << "STRICT_TYPE_I_CONTROL\tNOT_GUARANTEED\n"
            << "GENOMIC_CONTROL_APPLIED\tfalse\n"
            << "SCALABILITY_SCOPE\t"
            << (options_.ssgwas_hinv_mode == "operator"
                    ? "GENOMIC_AND_PEDIGREE_PRECISION_OPERATOR_PATH"
                    : "UNCERTAINTY_LAYER_ONLY (dense A22 inverse remains upstream)") << "\n"
            << "TEST_STATISTIC\tZ\n"
            << "REFERENCE_DISTRIBUTION\tSTANDARD_NORMAL\n"
            << "SIDES\t2\n"
            << "REFERENCE_FORMULA\tLeite et al. 2024 GSE 56:59 (PMC11340074) Eqs. 3/8/9 (Approx_GinvAPY); z = effect/sd(effect), P = 2(1-Phi(|z|))\n"
            << "G_INVERSE_METHOD\tGINV_APY\n"
            << "G_INVERSE_REFERENCE\tLeite et al. 2024 Eq. 5-6 (G_APY^-1 = [Gcc^-1 + Q Mnn^-1 Q' -Q Mnn^-1; -Mnn^-1 Q' Mnn^-1], Q = Gcc^-1 Gcn, Mnn = diag(Gnn - Gnc Gcc^-1 Gcn))\n"
            << "APY_CORE_SIZE\t" << apy_n_core_ << "\n"
            << "APY_NONCORE_SIZE\t" << apy_n_noncore_ << "\n"
            << "APY_CORE_SELECTION\t" << (apy_core_selection_mode_ == "FILE"
                                              ? "USER_FILE" : apy_core_selection_mode_) << "\n"
            << "APY_SEED\t" << (apy_core_selection_mode_ == "RANDOM_SEEDED"
                                    ? std::to_string(options_.ssgwas_apy_seed) : "NA") << "\n"
            << "APY_CORE_FILE\t" << (apy_core_selection_mode_ == "FILE"
                                         ? options_.ssgwas_apy_core_file : "NA") << "\n"
            << "CORE_TARGET_RATIO\t" << (options_.ssgwas_apy_core_target_ratio > 0.0
                                             ? std::to_string(options_.ssgwas_apy_core_target_ratio) : "NA") << "\n"
            << "K98_G\t" << k98_g_ << "\n"
            << "K98_G_CUMFRACTION\t" << k98_cumfraction_ << "\n"
            << "APY_CORE_AUTO\t" << (apy_core_auto_ ? "true" : "false") << "\n"
            << "APY_CORE_AUTO_VARIANCE\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_variance_) : "NA") << "\n"
            << "APY_CORE_AUTO_K98\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K98_) : "NA") << "\n"
            << "APY_CORE_AUTO_K99\t" << (apy_core_auto_ ? std::to_string(apy_core_auto_K99_) : "NA") << "\n"
            << "APY_CORE_EIGENSOLVER\t" << (apy_core_auto_ ? apy_core_eigensolver_used_ : "NA") << "\n"
            << "APY_CORE_TRACE_G\t" << (apy_core_auto_ ? std::to_string(apy_core_trace_g_) : "NA") << "\n"
            << "APY_CORE_EIG_ITERATIONS\t" << (apy_core_auto_ ? std::to_string(apy_core_eig_iterations_) : "NA") << "\n"
            << "APY_GCC_MIN_LDLT_PIVOT\t" << apy_gcc_min_diag_d_ << "\n"
            << "APY_GCC_MAX_LDLT_PIVOT\t" << apy_gcc_max_diag_d_ << "\n"
            << "APY_MNN_MIN\t" << apy_mnn_min_ << "\n"
            << "APY_MNN_MAX\t" << apy_mnn_max_ << "\n"
            << "MME_DIMENSION\t" << (n_ped_ + 1) << "\n"
            << "MME_INVERSE_METHOD\tNONE (no selected inverse; Eq. 8 C_core_approx only)\n"
            << "MME_SCALE\tLAMBDA_FORM (C_lambda solves [X'X X'Z; Z'X Z'Z+lambda*Hinv], lambda = Ve/Va; genotyped block = tau*G_APY^-1 - omega*A22^-1)\n"
            << "PEV_SCALE\tSIGMA_E2 * Kc^-1 (approximate core PEV; Kc is Leite Eq. 8 Schur precision)\n"
            << "EFFECT_SCALE\tALLELE_DOSAGE (EFFECT_ALLELE)\n"
            << "BACKSOLVE_OPERATOR\tB_j = M_cj' * (G_backsolve^-1)_cc / S (core-only weight, Bermann 2022a Eq. 21; identical operator for effect and variance)\n"
            << "G_EVALUATION_CONVENTION\tblended Gw (final, as in H^-1 genotyped block); APY-inverted\n"
            << "G_BACKSOLVE_CONVENTION\t" << effect_ginv_name_
            << " (APY inverse of tuned-unblended G* when --ssgwas-blend-a22, else of blended Gw)\n"
            << "CORE_UNCERTAINTY_METHOD\tLEITE2024_EQ8_BERMANN2022_SCHUR\n"
            << "CORE_UNCERTAINTY_REFERENCE\tC_core^-1 = D_c + lambda*Gcc^-1 + Q*S*Q' with S_j = lambda*d_j/(lambda + d_j*m_jj); Schur/Woodbury form algebraically equivalent to Leite et al. 2024 Eq. 8 written with the APY-inverse blocks (equivalence gate D1.5, 100-case numeric, max|diff| < 1e-12)\n"
            << "D_EFFECTIVE_RECORDS\tLeite et al. 2024 Eq. 9: diag(W'(I - X(X'X)^-1 X')W), exact from record-level design\n"
            << "GEBV_UNCERTAINTY_METHOD\tLEITE2024_EQ8_CORE_PEV_APPROXIMATION\n"
            << "CORE_MATRIX_DIM\t" << apy_n_core_ << "\n"
            << "CORE_MATRIX_BYTES\t"
            << (options_.ssgwas_pev_apy_mode == "dense"
                    ? 8.0 * apy_n_core_ * apy_n_core_ : 0.0) << "\n"
            << "PEV_APY_BACKEND\t" << apy_pev_backend_ << "\n"
            << "PEV_APY_MAX_RELATIVE_RESIDUAL\t" << apy_pev_max_relative_residual_ << "\n"
            << "FULL_Q_REQUIRED\t"
            << (options_.ssgwas_pev_apy_mode == "dense" ||
                        options_.ssgwas_mme_apy_mode == "dense" ||
                        options_.ssgwas_apy_inverse_mode == "dense" ? "true" : "false") << "\n"
            << "CORE_MATRIX_REQUIRED\t"
            << (options_.ssgwas_pev_apy_mode == "dense" ||
                        options_.ssgwas_mme_apy_mode == "dense" ||
                        options_.ssgwas_apy_inverse_mode == "dense" ? "true" : "false") << "\n"
            << "DENSE_NG2_ALLOCATED\tfalse\n"
            << "UPSTREAM_DENSE_NG2_PRESENT\t"
            << (grm_streamed_ && options_.ssgwas_hinv_mode == "operator" && A22_geno_.size() == 0
                    ? "false"
                    : (grm_streamed_ ? "true (A22 is dense for tuning/blending)"
                                     : "true (dense G and/or A22 path)")) << "\n"
            << "DENSE_UNCERTAINTY_ALLOCS\t" << dense_uncertainty_allocs_ << "\n"
            << "ESTIMATED_PEAK_MEMORY_BYTES\t" << estimated_peak_memory_bytes_ << "\n"
            << "SNP_BLOCK_SIZE\t" << pev_scan_effective_block_size_ << "\n"
            << "FULL_MXM_REQUIRED\tfalse\n"
            << "TUNE_G\t" << (options_.ssgwas_tune_g ? "true" : "false") << "\n"
            << "TUNE_A\t" << tune_a_ << "\n"
            << "TUNE_B\t" << tune_b_ << "\n"
            << "BLEND_ALPHA\t" << options_.blend_g << "\n"
            << "TAU\t" << options_.h_tau << "\n"
            << "OMEGA\t" << options_.h_omega << "\n"
            << "EXACT_OR_APPROX\tAPPROXIMATE\n"
            << "VAR_GHAT_MIN_EIG\t" << var_ghat_min_eig_ << "\n"
            << "VAR_GHAT_MAX_EIG\t" << var_ghat_max_eig_ << "\n"
            << "VAR_GHAT_EIG_MODE\t" << var_ghat_eig_mode_ << "\n"
            << "PEC_NEGATIVE_CLAMPED\t" << n_pec_clamped_ << "\n"
            << "PEC_NUMERICAL_FAILURE\t" << n_pec_failure_ << "\n"
            << "ZERO_SE_COUNT\t" << n_zero_se_ << "\n"
            << "TIME_PEV_CORE_INVERSE\t" << t_pev_inverse_ << "\n"
            << "TIME_PEV_VHAT\t" << t_pev_vhat_ << "\n"
            << "TIME_PEV_SCAN\t" << t_pev_scan_ << "\n"
            << "PEV_SCAN_MODE\t" << options_.ssgwas_pev_scan_mode << "\n"
            << "PEV_SCAN_REQUESTED_BLOCK_SIZE\t" << options_.ssgwas_pev_scan_block_size << "\n"
            << "PEV_SCAN_BLOCK_SIZE\t" << pev_scan_effective_block_size_ << "\n"
            << "PEV_SCAN_BLOCKS\t" << pev_scan_blocks_ << "\n"
            << "PEV_SCAN_VARIANTS_DECODED\t" << pev_scan_variants_decoded_ << "\n"
            << "PEV_SCAN_BLOCK_REALLOCATIONS\t" << pev_scan_block_reallocations_ << "\n"
            << "PEV_SCAN_PEAK_TRANSIENT_BYTES\t" << pev_scan_peak_transient_bytes_ << "\n"
            << "PEV_SCAN_MULTIPLY_Z_V_CALLS\t" << pev_scan_multiply_z_v_calls_ << "\n"
            << "PEV_SCAN_GENOTYPE_PASSES\t" << pev_scan_genotype_passes_ << "\n"
            << "PEV_SCAN_BED_BYTES_READ_ESTIMATE\t" << pev_scan_bed_bytes_read_estimate_ << "\n"
            << "PEV_SCAN_VARIANTS_TESTED\t" << pev_scan_variants_tested_ << "\n"
            << "PEV_SCAN_HEAP_ALLOCATIONS\t" << pev_scan_allocations_ << "\n"
            << "PEV_SCAN_GCC_SOLVE_CALLS\t" << pev_scan_gcc_solve_calls_ << "\n"
            << "PEV_SCAN_GCC_RHS_COLUMNS\t" << pev_scan_gcc_rhs_columns_ << "\n"
            << "PEV_SCAN_GW_GEMM_CALLS\t" << pev_scan_gw_gemm_calls_ << "\n"
            << "PEV_SCAN_PEC_GEMM_CALLS\t" << pev_scan_pec_gemm_calls_ << "\n"
            << "PEV_SCAN_MISSING_TO_ZERO\t" << pev_scan_missing_to_zero_ << "\n"
            << "TIME_PEV_SCAN_SETUP\t" << t_pev_scan_setup_ << "\n"
            << "TIME_PEV_SCAN_ONEHOT\t" << t_pev_scan_onehot_ << "\n"
            << "TIME_PEV_SCAN_GENOTYPE_READ_MULTIPLY\t" << t_pev_scan_genotype_ << "\n"
            << "TIME_PEV_SCAN_CORE_GATHER\t" << t_pev_scan_core_gather_ << "\n"
            << "TIME_PEV_SCAN_GCC_SOLVE\t" << t_pev_scan_gcc_solve_ << "\n"
            << "TIME_PEV_SCAN_GW_MATVEC\t" << t_pev_scan_gw_matvec_ << "\n"
            << "TIME_PEV_SCAN_PEC_MATVEC\t" << t_pev_scan_pec_matvec_ << "\n"
            << "TIME_PEV_SCAN_QUADRATIC\t" << t_pev_scan_quadratic_ << "\n"
            << "TIME_PEV_SCAN_SE_Z_P_LOG10P\t" << t_pev_scan_statistics_ << "\n"
            << "TIME_PEV_SCAN_OTHER\t" << t_pev_scan_other_ << "\n"
            << "PEV_SCAN_TIME_PER_SNP\t"
            << (pev_scan_variants_tested_ > 0 ? t_pev_scan_ / pev_scan_variants_tested_ : 0.0) << "\n"
            << "TIME_PEV_OUTPUT\t" << t_pev_out_ << "\n"
            << "SOLVER_CONVERGED\t" << (solver_converged_ ? "true" : "false") << "\n";
    } else {
        ofs << "P_VALUE_METHOD\tOFF\n";
    }
    ofs.close();
    std::cout << "Wrote " << path << "\n";
}

} // namespace cosmic
