#pragma once
#include <string>
#include <vector>

namespace cosmic {

struct GwasConfig {
    std::string bfile_prefix;
    std::string bgen_file;
    std::string sample_file;
    std::string pfile_prefix;
    std::string pgen_file;
    std::string pvar_file;
    std::string psam_file;
    std::string out_prefix = "cosmic_gwas";

    std::string pheno_file;
    std::string covar_file;
    std::string qcovar_file;
    int mpheno = 1;
    int tpheno = 2;
    bool tpheno_set = false;
    std::vector<int> mpheno_list;

    bool run_gwas = false;
    std::string scan_tier;
    bool no_loco = false;
    bool fast_lmm = false;
    bool linear_model = false;
    int linear_model_mode = 4;
    int lmm_gemma_mode = 0;
    bool lmm_gemma_refit = false;
    bool repeatability_model = false;
    std::string repeatability_mode = "auto";
    std::string repeatability_grm_cache_prefix;
    std::string repeatability_weight_file;
    bool repeatability_mixed_precision = false;
    bool repeatability_refit_validation = false;
    int repeatability_refit_limit = 0;
    std::string vce_refine_mode = "auto";
    bool multivariate = false;
    std::string multivariate_method = "emmv";
    bool multivariate_method_explicit = false;
    std::string multivariate_compatible;
    std::string multivariate_tests = "wald,lrt,score";
    int multivariate_gemma_mode = 4;
    double multivariate_pnr = 0.001;
    int multivariate_emi = 10000;
    int multivariate_nri = 100;
    double multivariate_emp = 0.0001;
    double multivariate_nrp = 0.0001;
    bool multivariate_crt = false;
    int multivariate_refine_top = 100;
    double multivariate_refine_p = 1e-3;
    int multivariate_refine_discrepancy_top = 0;
    double multivariate_refine_discrepancy = 1.0;
    int multivariate_refine_score_prepass = 0;
    double multivariate_refit_p = 1.0;
    bool multivariate_use_aireml = false;
    int multivariate_em_min_iter = 10;
    int multivariate_ai_max_iter = 50;
    double multivariate_ai_tol = 1e-8;
    bool multivariate_output_vbeta = false;
    bool multivariate_null_only = false;
    int multivariate_gemma_force_iterations = -1;
    std::string k_file;
    std::string multivariate_k_file;
    std::string multivariate_compare_vg;
    std::string multivariate_compare_ve;
    std::string multivariate_init_vg;
    std::string multivariate_init_ve;
    bool multivariate_emmvgwas_output = false;
    bool multivariate_emmvgwas_vc_only = false;
    std::string multivariate_emmvgwas_G_file;
    std::string multivariate_emmvgwas_R_file;
    bool multivariate_eb = false;
    int env_covar_col = -1;
    bool single_step_model = false;
    bool make_agh = false;
    bool write_agh = false;
    std::string ped_file;
    std::string genotyped_list;
    double h_tau = 1.0;
    double h_omega = 1.0;
    double blend_g = 0.05;
    bool write_ebv = false;

    int grm_alg = 1;
    double min_maf = 0.01;
    double min_hwe = 0.0;
    std::string lmm_null_solver = "brent";
    int reml_maxit = 50;
    double reml_tol = 1e-6;
    double reml_vg = -1.0;
    double reml_vpe = -1.0;
    double reml_ve = -1.0;
    bool use_covariate_projection = true;
    bool out_ssf = false;
    int num_threads = 0;
    std::string extract_file;
    std::string exclude_file;
    std::string scan_chr;
    long long scan_from_bp = -1;
    long long scan_to_bp = -1;

    // Unsupported fields kept only so retained implementation files compile.
    bool eb_lmm = false;
    bool use_spa = false;
    bool run_skat = false;
    bool run_burden = false;
    bool biobank_mode = false;
    int block_size = 1000;
    std::string set_file;
    std::string regenie_model_snps;
    std::string regenie_write_pred;
    std::string regenie_use_pred;
    int regenie_step = 0;
    bool make_grm_sparse = false;
    bool force_slow = false;
    double sparse_cutoff = 0.05;
    std::string snp_weight_file;
    std::string pop_class_file;
};

}  // namespace cosmic
