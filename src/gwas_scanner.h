#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "data_filter.h"
#include "lmm_engine.h"
#include "grm_builder.h"
#include "pheno_reader.h"
#include "rare_variant_test.h"
#include "ridge_regenie.h"
#include "eigen_vce.h"
#include "eigen_vce_ad.h"
#include "ad_engine.h"
#include "dominance_engine.h"
#include "validation.h"

namespace cosmic {

struct BgenAdResearchOptions;
struct BgenAdIntegrationDiagnostics;

class GwasScanner {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        bool out_ssf = false; // Add switch for EBI GWAS-SSF output
        bool out_gz = false; // Compress output with gzip
        bool manhattan_tsv = false; // Generate plot-ready Manhattan TSV
        bool qq_tsv = false; // Generate plot-ready QQ TSV
        double min_maf = 0.01;
        double min_hwe = 0.0;
        double bgen_dom_max_missing = 0.05;
        double bgen_dom_min_info = 0.01;
        double bgen_dom_entropy_warning = 0.50;
        std::string snp_model = "additive";
        std::string gwas_model_md = "additive";
        bool fast_lmm = false;
        bool fastlmm = false;
        double fastlmm_trace_threshold = 0.95;
        int fastlmm_min_rank = 1;
        int fastlmm_max_rank = 256;
        bool fastlmm_matrix_free = false; // Use matrix-free low-rank path
        bool fastlmm_no_fallback = false; // Don't fall back to full rank if trace threshold not met
        int lmm_low_rank = 0;
        int lmm_low_rank_oversample = 10;
        long long lmm_low_rank_seed = 42;
        double lmm_low_rank_error_tol = 1.0;
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
        std::string rare_method;
        bool run_burden = false;
        double rare_max_maf = 0.5;
        int rare_min_mac = 1;
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
    const ValidationReport& validationReport() const { return validation_report_; }

    // Run standard LOCO scan (or single chromosome scan if only 1 chrom)
    void runLocoScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        VcfReader& vcf_reader,
        bool is_bgen,
        bool is_pgen,
        bool is_vcf,
        const DataFilter::KeepInfo& keep,
        SequentialGenotypeMatrix* generic_reader = nullptr
    );

    // Run Additive-Dominance (a+d+A+D) joint scan on PLINK BED or PGEN hardcall input.
    // Builds K_A (algorithm 1) and K_D (algorithm 4) GRMs, fits the AD null
    // per chromosome (LOCO) via EigenVceAD, and scans every SNP with the
    // joint additive + dominance design (see LmmEngineAD).
    // is_pgen selects the PGEN hardcall backend; PGEN supplies exact
    // AA/AB/BB/missing hardcalls, so the design columns match BED exactly.
    void runAdLocoScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        PgenReader& pgen_reader,
        bool is_pgen,
        const DataFilter::KeepInfo& keep
    );

    // Explicit internal research entry only. Normal CLI/univariate dispatch
    // remains guarded for BGEN additive_dominance.
    void runBgenAdResearchScan(const Eigen::VectorXd& y, const Eigen::MatrixXf& X,
                              BgenReader& reader, const DataFilter::KeepInfo& keep,
                              const BgenAdResearchOptions& research,
                              BgenAdIntegrationDiagnostics* diagnostics = nullptr);

    // Run Dominance-only (d+D) scan on BED/PGEN hardcalls or BGEN probabilities.
    void runDominanceLocoScan(
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

    // Run Rare Variant Association Test (SKAT / Burden)
    void runRareVariantScan(
        const Eigen::VectorXd& y,
        const Eigen::MatrixXf& X,
        PlinkReader& plink_reader,
        BgenReader& bgen_reader,
        PgenReader& pgen_reader,
        VcfReader& vcf_reader,
        bool is_bgen,
        bool is_pgen,
        bool is_vcf,
        const DataFilter::KeepInfo& keep
    );

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
        VcfReader& vcf_reader,
        bool is_bgen,
        bool is_pgen,
        bool is_vcf,
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
    ValidationReport validation_report_;

    std::unordered_set<std::string> extract_ids_;
    std::unordered_set<std::string> exclude_ids_;
    bool sumstats_header_written_ = false;
    std::vector<GwasResult> plot_results_;  // accumulated for Manhattan/QQ
    std::vector<GwasResult> assoc_results_; // accumulated for single-pass assoc output

    // Additive-Dominance scan accumulation (--gwas-model additive_dominance)
    std::vector<AdResult> ad_results_;
    long long ad_tested_markers_ = 0;
    long long ad_rank_deficient_markers_ = 0;
    void appendAdResults(const std::vector<AdResult>& results);
    void flushAdResults(const std::string& filename);

    std::vector<DomResult> dom_results_;
    long long dom_tested_markers_ = 0;
    long long dom_rank_deficient_markers_ = 0;
    void appendDomResults(const std::vector<DomResult>& results);
    void flushDomResults(const std::string& filename);

    // Helper to append results to file
    void appendResults(const std::string& filename, const std::vector<GwasResult>& results);

    // Flush accumulated assoc results to file in one pass (header + data)
    void flushAssocResults(const std::string& filename);

    // Helper to write header
    void writeHeader(const std::string& filename);
    void writeEbSummary(const std::string& filename);
    void appendCosmicSumstats(const std::string& filename, const std::vector<GwasResult>& results, const std::string& model_name);
    void writePlotTsv(const std::string& prefix, const std::vector<GwasResult>& results);
};

} // namespace cosmic
