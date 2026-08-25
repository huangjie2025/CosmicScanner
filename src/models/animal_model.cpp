#include "models/animal_model.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include "data_loader.h"
#include "eigen_vce.h"
#include "lmm_engine.h"
#include "a_matrix.h"

namespace cosmic {

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

double nowSeconds() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

} // anonymous namespace

AnimalModel::Options AnimalModel::Options::fromConfig(const GwasConfig& cfg) {
    Options opts;
    opts.out_prefix = cfg.out_prefix;
    opts.min_maf = cfg.min_maf;
    opts.min_hwe = cfg.min_hwe;
    opts.reml_maxit = cfg.reml_maxit;
    opts.reml_tol = cfg.reml_tol;
    opts.mpheno = cfg.mpheno;
    opts.pheno_file = cfg.pheno_file;
    opts.qcovar_file = cfg.qcovar_file;
    opts.covar_file = cfg.covar_file;
    opts.ped_file = cfg.ped_file;
    opts.extract_file = cfg.extract_file;
    opts.exclude_file = cfg.exclude_file;
    opts.write_ped_kernel = cfg.write_ped_kernel;
    opts.no_loco = cfg.no_loco;
    opts.n_threads = 1;
#ifdef _OPENMP
    opts.n_threads = omp_get_max_threads();
#endif
    return opts;
}

AnimalModel::AnimalModel(const Options& opts) : options_(opts) {}

void AnimalModel::prepare(const std::vector<std::string>& fids,
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
    if (keep_.n_keep < 10) throw std::runtime_error("Too few valid samples for the animal model.");
}

void AnimalModel::loadPhenotype() {
    y_ = DataLoader::loadPhenotype(options_.pheno_file, options_.mpheno, fids_, iids_,
                                   n_samples_, *plink_reader_, is_bgen_, is_pgen_);
}

void AnimalModel::loadCovariates() {
    X_ = DataLoader::loadCovariates(options_.qcovar_file, options_.covar_file, fids_, iids_,
                                    n_samples_);
}

void AnimalModel::buildPedigreeRelationship() {
    double t0 = nowSeconds();
    pedigree_.load(options_.ped_file);
    t_ped_parse_ = nowSeconds() - t0;

    t0 = nowSeconds();
    const Eigen::SparseMatrix<double>& Ainv = pedigree_.buildAInverse();
    t_ainv_ = nowSeconds() - t0;
    std::cout << "A⁻¹ built with " << Ainv.nonZeros() << " non-zero entries.\n";

    // Map kept (analysis) samples to pedigree indices by ID. Tolerates
    // "FID:IID" and plain "IID" spellings; every analysis individual must be
    // present in the pedigree, and no duplicate mapping is allowed.
    t0 = nowSeconds();
    analysis_pedigree_indices_.clear();
    std::unordered_set<int> seen;
    const auto& id_to_idx = pedigree_.idToIndex();
    for (int s : keep_.keep_indices) {
        const std::string fid_iid = fids_[s] + ":" + iids_[s];
        const std::string& iid = iids_[s];
        auto it = id_to_idx.find(fid_iid);
        if (it == id_to_idx.end()) it = id_to_idx.find(iid);
        if (it == id_to_idx.end()) {
            throw std::runtime_error("Analysis sample is absent from pedigree: " + fid_iid);
        }
        if (!seen.insert(it->second).second) {
            throw std::runtime_error("Multiple analysis samples map to pedigree ID: " + fid_iid);
        }
        analysis_pedigree_indices_.push_back(it->second);
    }

    A22_ = pedigree_.buildA22(analysis_pedigree_indices_);
    t_a22_ = nowSeconds() - t0;
    std::cout << "Pedigree A22 (analysis subset " << A22_.rows() << "x" << A22_.cols()
              << ") built via Schur complement of sparse A⁻¹.\n";
}

void AnimalModel::runVCE() {
    double t0 = nowSeconds();
    if (A22_.rows() != keep_.n_keep) {
        throw std::runtime_error("Animal model A22 dimension does not match the analysis sample set.");
    }

    Eigen::VectorXd y_sub = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X_, keep_);

    EigenVce::Options vce_opts;
    vce_opts.use_covariate_projection = true;
    vce_opts.max_iter = options_.reml_maxit;
    vce_opts.tol = options_.reml_tol;
    EigenVce vce(vce_opts);
    vce.prepare(y_sub, X_sub, A22_);
    bool ok = vce.runNullModel();
    vce_converged_ = ok && vce.converged();
    vce_iterations_ = vce.lastIterations();

    sigma2_a_ = vce.getSigma2g();
    sigma2_e_ = vce.getSigma2e();

    if (!std::isfinite(sigma2_a_) || !std::isfinite(sigma2_e_) ||
        sigma2_a_ < 0.0 || sigma2_e_ < 0.0) {
        throw std::runtime_error("Animal model Eigen REML produced invalid variance components.");
    }
    t_vce_ = nowSeconds() - t0;

    std::cout << "Animal model VCE: sigma2_a=" << sigma2_a_ << " sigma2_e=" << sigma2_e_
              << " h2=" << sigma2_a_ / (sigma2_a_ + sigma2_e_ + 1e-10)
              << " converged=" << (vce_converged_ ? "yes" : "no") << "\n";
}

void AnimalModel::scanGwas() {
    double t0 = nowSeconds();

    // LOCO is not applicable: the pedigree relationship matrix A is identical
    // for every chromosome.
    std::set<std::string> chrom_labels_set;
    if (is_bgen_) {
        bgen_reader_->loadVariantInfo();
        for (const auto& v : bgen_reader_->getVariants()) chrom_labels_set.insert(v.chrom);
    } else if (is_pgen_) {
        for (const auto& v : pgen_reader_->getVariants()) chrom_labels_set.insert(v.chrom);
    } else {
        for (const auto& snp : plink_reader_->getBimInfo()) chrom_labels_set.insert(snp.chrom);
    }
    if (chrom_labels_set.size() > 1 && !options_.no_loco) {
        std::cout << "WARNING: LOCO is not applicable to pedigree-only Animal Model GWAS; "
                     "using the pedigree relationship matrix A for all chromosomes.\n";
    }

    Eigen::VectorXd y_sub = DataFilter::subset(y_, keep_);
    Eigen::MatrixXf X_sub = DataFilter::subsetRows(X_, keep_);

    LmmEngine::Options lmm_opts;
    lmm_opts.min_maf = (float)options_.min_maf;
    lmm_opts.min_hwe = (float)options_.min_hwe;
    lmm_opts.use_covariate_projection = true;
    lmm_opts.extract_set = loadSnpIdSet(options_.extract_file);
    lmm_opts.exclude_set = loadSnpIdSet(options_.exclude_file);
#ifdef _OPENMP
    lmm_opts.thread_num = options_.n_threads;
#endif

    std::string out_file = options_.out_prefix + ".assoc.txt";
    {
        std::ofstream ofs(out_file);
        ofs << "CHR\tSNP\tBP\tA1\tA2\tFREQ\tN\tBETA\tSE\tP\tLOG10P\n";
    }

    EigenVce::Options vce_opts;
    vce_opts.use_covariate_projection = true;
    EigenVce vce(vce_opts);
    vce.prepare(y_sub, X_sub, A22_);
    vce.fixedComponents(sigma2_a_, sigma2_e_);

    LmmEngine engine(lmm_opts);
    engine.prepareFromVce(vce);

    std::vector<GwasResult> results;
    if (is_bgen_) {
        bgen_reader_->reset();
        engine.scan(*bgen_reader_, results, keep_.keep_indices, "");
    } else if (is_pgen_) {
        pgen_reader_->reset();
        engine.scan(*pgen_reader_, results, keep_.keep_indices, "");
    } else {
        engine.scan(*plink_reader_, results, keep_.keep_indices, "");
    }

    std::ofstream fout(out_file, std::ios::app);
    for (const auto& r : results) {
        if (std::isnan(r.p_value)) continue;
        fout << r.chrom << "\t" << r.snp_id << "\t" << r.pos << "\t"
             << r.a1 << "\t" << r.a2 << "\t" << r.freq << "\t"
             << r.n_obs << "\t" << r.beta << "\t" << r.se << "\t"
             << r.p_value << "\t" << r.log10_p << "\n";
    }
    t_scan_ = nowSeconds() - t0;
    std::cout << "Animal model results saved to " << out_file << "\n";
}

void AnimalModel::writeKernel() const {
    if (!options_.write_ped_kernel) return;
    const std::string fname = options_.out_prefix + ".animal_model.k.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    ofs.precision(12);
    for (int i = 0; i < A22_.rows(); ++i) {
        for (int j = 0; j < A22_.cols(); ++j) {
            if (j) ofs << " ";
            ofs << A22_(i, j);
        }
        ofs << "\n";
    }
    ofs.close();
    std::ofstream idf(options_.out_prefix + ".animal_model.k.id");
    for (int s : keep_.keep_indices) {
        idf << fids_[s] << ":" << iids_[s] << "\n";
    }
    idf.close();
    std::cout << "Animal model pedigree kernel A22 written to " << fname << "\n";
}

void AnimalModel::writeNullSummary() const {
    std::string fname = options_.out_prefix + ".animal_model.null.txt";
    std::ofstream ofs(fname);
    if (!ofs.is_open()) {
        std::cerr << "Warning: cannot write " << fname << "\n";
        return;
    }
    const double total = sigma2_a_ + sigma2_e_;
    ofs << "Item\tValue\n"
        << "n_pedigree\t" << pedigree_.nPedigree() << "\n"
        << "n_subjects\t" << keep_.n_keep << "\n"
        << "sigma2_a\t" << sigma2_a_ << "\n"
        << "sigma2_e\t" << sigma2_e_ << "\n"
        << "h2\t" << ((total > 0.0) ? (sigma2_a_ / total) : 0.0) << "\n"
        << "converged\t" << (vce_converged_ ? "yes" : "no") << "\n"
        << "iterations\t" << vce_iterations_ << "\n"
        << "solver\teigen_reml_brent\n"
        << "scan_mode\tanimal_model_gls\n"
        << "timing_pedigree_parse_s\t" << t_ped_parse_ << "\n"
        << "timing_ainv_build_s\t" << t_ainv_ << "\n"
        << "timing_a22_extraction_s\t" << t_a22_ << "\n"
        << "timing_vce_s\t" << t_vce_ << "\n"
        << "timing_scan_s\t" << t_scan_ << "\n";
    ofs.close();
    std::cout << "Animal model null summary written to " << fname << "\n";
}

void AnimalModel::run() {
    std::cout << "--- Pedigree Animal Model GWAS ---\n";
    std::cout << "  y = Xb + x_j*beta_j + Za + e,  a ~ N(0, A*sigma2_a)\n";
    buildPedigreeRelationship();
    runVCE();
    writeKernel();
    scanGwas();
    writeNullSummary();
}

} // namespace cosmic
