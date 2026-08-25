#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "lmm_engine.h"
#include "grm_builder.h"
#include "pheno_reader.h"
#include "ridge_regenie.h"
#include "eigen_vce.h"

#ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
#include "rare_variant_test.h"
#endif

namespace cosmic {

class GwasScanner {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        bool out_ssf = false; // Add switch for EBI GWAS-SSF output
        double min_maf = 0.01;
        double min_hwe = 0.0;
        bool fast_lmm = false;
        std::string scan_tier;
        bool profile_scale = false;
        int lmm_gemma_mode = 0;
        bool lmm_gemma_refit = false;
        bool eb_lmm = false;
        int n_threads = 1;
        bool no_loco = false;
        int grm_alg = 1; // 1: VanRaden, 0: Yang et al
        EigenVce::NullSolver lmm_null_solver = EigenVce::NullSolver::Brent;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        bool biobank_mode = false;
        bool use_covariate_projection = true; // Use exact covariate projection for LMM acceleration
        bool use_spa = false; // Use Saddle Point Approximation for binary traits
        int regenie_block_size = 1000;
        std::string regenie_model_snps;
        std::string regenie_write_pred;
        std::string regenie_use_pred;
        int regenie_step = 0;
        std::string scan_chr;
        long long scan_from_bp = -1;
        long long scan_to_bp = -1;
        bool run_skat = false;
        bool run_burden = false;
        std::string set_file;
        std::string extract_file;
        std::string exclude_file;
        int rrm_curve_grid = 21;
        int rrm_top = 100;
        int rrm_ai_iter = -1;
        double rrm_sig_threshold = 0.05;
        bool rrm_loco_global_null = false;
        std::string rrm_basis = "legendre";
        std::string rrm_gwas_method = "coefficient";

        Eigen::VectorXd obs_weights; // Optional weights for observations (e.g., Repeatability Model D^{-1/2})
    };

    GwasScanner(const Options& opts);

    // Run standard LOCO scan (or single chromosome scan if only 1 chrom)
    void runLocoScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

    // Run Biobank mode scan (Matrix-Free RHE + PCG)
    void runBiobankScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

    #ifdef COSMICSCANNER_ENABLE_EXPERIMENTAL
    // Run Rare Variant Association Test (SKAT / Burden)
    void runRareVariantScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );
#endif

    // Run Random Regression Model (RRM) scan for longitudinal data
    void runRrmScan(
        const std::vector<PhenoReader::PhenoRecord>& records,
        int rrm_order,
        const Eigen::MatrixXf& X_base,
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

    // Run unbalanced longitudinal (RRM) scan with pre-built fixed effects per observation
    void runRrmUnbalancedScan(
        const std::vector<PhenoReader::PhenoRecord>& records,
        int rrm_order,
        const Eigen::MatrixXf& X_base,
        const Eigen::MatrixXf& X_tdep,         // time-dependent covariates (n_obs x n_tdep)
        const std::vector<std::string>& tdep_names,
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        bool is_bgen,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

private:
    Options options;
    double eb_effective_tests_ = 0.0;
    long long eb_tested_markers_ = 0;

    std::unordered_set<std::string> extract_ids_;
    std::unordered_set<std::string> exclude_ids_;

    // Helper to append results to file
    void appendResults(const std::string& filename, const std::vector<GwasResult>& results);

    // Helper to write header
    void writeHeader(const std::string& filename);
    void writeEbSummary(const std::string& filename);
};

} // namespace cosmic
