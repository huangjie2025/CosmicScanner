#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"
#include "analysis_dataset.h"
#include "validation.h"

namespace cosmic {

class UnivariateModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        bool out_ssf = false;
        bool out_gz = false;
        bool manhattan_tsv = false;
        bool qq_tsv = false;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        double bgen_dom_max_missing = 0.05;
        double bgen_dom_min_info = 0.01;
        double bgen_dom_entropy_warning = 0.50;
        bool fast_lmm = false;
        bool fastlmm = false;
        double fastlmm_trace_threshold = 0.95;
        int fastlmm_min_rank = 1;
        int fastlmm_max_rank = 256;
        bool fastlmm_matrix_free = false;
        bool fastlmm_no_fallback = false;
        int lmm_low_rank = 0;
        int lmm_low_rank_oversample = 10;
        long long lmm_low_rank_seed = 42;
        double lmm_low_rank_error_tol = 1.0;
        std::string scan_tier;
        bool profile_scale = false;
        int lmm_gemma_mode = 0;
        bool lmm_gemma_refit = false;
        bool linear_model = false;
        int linear_model_mode = 4;
        std::string snp_model = "additive";
        std::string gwas_model = "additive";  // additive | additive_dominance | dominance
        std::string rare_method;
        bool farmcpu = false;
        bool blink = false;
        bool firth_logistic = false;
        bool saige_style = false;
        std::string logistic_test;
        bool logistic_glmm = false;
        bool conditional_joint = false;
        int farmcpu_max_iterations = 5;
        int farmcpu_pseudo_qtn = 10;
        std::string condition_snp_file;
        bool eb_lmm = false;
        int n_threads = 1;
        bool no_loco = false;
        int grm_alg = 1;
        EigenVce::NullSolver lmm_null_solver = EigenVce::NullSolver::Brent;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        bool biobank_mode = false;
        bool use_covariate_projection = true;
        bool use_spa = false;
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
        double rare_max_maf = 0.5;
        int rare_min_mac = 1;
        std::string set_file;
        std::string extract_file;
        std::string exclude_file;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        int mpheno = 1;
        double sparse_cutoff = 0.05;

        static Options fromConfig(const GwasConfig& cfg);
    };

    UnivariateModel(const Options& opts);

    void prepare(const std::vector<std::string>& fids,
                 const std::vector<std::string>& iids,
                 int n_samples,
                 PlinkReader& plink_reader,
                 BgenReader& bgen_reader,
                 PgenReader& pgen_reader,
                 VcfReader& vcf_reader,
                 bool is_bgen, bool is_pgen, bool is_vcf,
                 SequentialGenotypeMatrix* generic_reader = nullptr);

    void run();
    bool hasAnalysisDataset() const { return has_analysis_dataset_; }
    const AnalysisDataset& analysisDataset() const { return analysis_dataset_; }
    const ValidationReport& validationReport() const { return validation_report_; }

private:
    void runLinearModel();
    void runMultiLocusModel();
    void runSaigeStyle();
    void runLogistic();
    bool usesAnalysisDatasetContract() const;
    AnalysisDataset buildAnalysisDatasetForStandardLmm(PlinkReader& plink_reader,
                                                        BgenReader& bgen_reader,
                                                        PgenReader& pgen_reader,
                                                        VcfReader& vcf_reader) const;
    void applyAnalysisDatasetForStandardLmm(const AnalysisDataset& dataset);

    Options options_;
    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false, is_vcf_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    VcfReader* vcf_reader_ = nullptr;
    SequentialGenotypeMatrix* generic_reader_ = nullptr;
    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;
    AnalysisDataset analysis_dataset_;
    ValidationReport validation_report_;
    bool has_analysis_dataset_ = false;
};

} // namespace cosmic
