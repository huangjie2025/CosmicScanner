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
    std::string long_pheno_file;
    std::string covar_file;
    std::string qcovar_file;
    std::string input_file;
    std::string input_config;
    std::string write_template_config;
    std::string input_format = "auto";
    std::string input_model = "repeatability";
    std::string prepared_input_dir;
    std::string id_col;
    std::string iid_col;
    std::string time_col;
    std::string trait_cols;
    std::string input_trait;
    std::string record_qcovars;
    std::string record_ccovars;
    std::string subject_qcovars;
    std::string subject_ccovars;
    std::string time_parser = "auto";
    std::string time_unique_policy = "within-time-index";
    std::string time_index_multiplier = "auto";
    bool prepare_input_only = false;
    bool preview_schema = false;
    bool infer_columns = false;
    bool wide_to_long = false;
    int mpheno = 1;
    int tpheno = 2;
    bool tpheno_set = false;
    int rrm_order = 3;
    int rrm_curve_grid = 21;
    int rrm_top = 100;
    int rrm_ai_iter = -1;
    double rrm_sig_threshold = 0.05;
    bool rrm_unbalance = false;
    bool rrm_loco_global_null = false;
    std::string rrm_basis = "legendre";
    std::string rrm_compatible;
    std::string rrm_compatible_varcom;
    std::vector<int> mpheno_list;
    std::string tpoint_col_name;
    std::string tclass_names;
    std::string tcovar_names;
    bool multivariate = false;
    std::string multivariate_method = "emmv"; // emmv, wald, gemma, gemma-full, gemma-exact, gemma-gxe, gemma-gxe-full, hybrid
    bool multivariate_eb = false;
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
    bool multivariate_emmvgwas_output = false;
    bool multivariate_emmvgwas_vc_only = false;
    std::string multivariate_emmvgwas_G_file;
    std::string multivariate_emmvgwas_R_file;
    bool repeatability_model = false;
    bool rrm_model = false;
    bool maternal_model = false;
    bool single_step_model = false;
    bool threshold_model = false;
    bool count_model = false;
    bool survival_model = false;
    bool gxe_model = false;
    std::string ped_file;
    bool make_grm = false;
    bool make_grm_sparse = false;
    int gemma_s_mode = 0; // GEMMA -gs: S matrix for summary-stat VC
    int gemma_q_mode = 0; // GEMMA -gq: q vector for summary-stat VC
    std::string beta_file;
    std::string cor_file;
    std::string category_file;
    std::string multi_category_file;
    std::string category_weight_file;
    std::string snp_weight_summary_file;
    std::string study_prefix;
    std::string ref_prefix;
    std::string multi_study_file;
    std::string multi_ref_file;
    int summary_block_count = 200;
    bool gemma_eigen = false;
    bool gemma_calccor = false;
    long long calccor_window_bp = 1000000;
    double calccor_window_cm = 0.0;
    int calccor_window_ns = 0;
    int vc_mode = 0; // GEMMA -vc: 1 HE, 2 REML, 3 AI/IDUL-style REML
    int ci_mode = 0; // GEMMA -ci: summary-stat confidence interval route
    bool run_gwas = false;
    bool fast_lmm = false;
    int lmm_gemma_mode = 0; // GEMMA -lmm: 1 Wald, 2 LRT, 3 Score, 4 all; 0 disabled
    bool lmm_gemma_refit = false;
    bool linear_model = false;
    int linear_model_mode = 4; // GEMMA -lm: 1 Wald, 2 LRT, 3 Score, 4 all
    bool eb_lmm = false;
    bool use_spa = false;
    bool run_skat = false;
    bool run_burden = false;
    std::string set_file = "";
    bool biobank_mode = false;
    int block_size = 1000;
    std::string regenie_model_snps;
    std::string regenie_write_pred;
    std::string regenie_use_pred;
    int regenie_step = 0; // 0: run Step1+Step2, 1: Step1 only, 2: Step2 only
    std::string predict_effect_file;
    std::string predict_score_file;
    std::string predict_snp_col;
    std::string predict_effect_col;
    int bslmm_mode = 0; // GEMMA -bslmm: 1 BSLMM, 2 ridge/GBLUP, 3 probit BSLMM
    int bslmm_top = 100;
    double bslmm_h2 = 0.5;
    std::string bslmm_bf_file;
    std::string bslmm_hyp_file;
    int bslmm_hgrid = 5;
    int bslmm_rho_grid = 5;
    double bslmm_h_min = 0.05;
    double bslmm_h_max = 0.95;
    double bslmm_rho_min = 0.05;
    double bslmm_rho_max = 0.95;
    double bslmm_p_min = -4.0; // GEMMA -pmin: minimum log10(pi)
    double bslmm_p_max = 0.0;  // GEMMA -pmax: maximum log10(pi)
    bool bslmm_p_min_set = false;
    bool bslmm_p_max_set = false;
    double bslmm_gamma_mean = 0.0;
    double bslmm_h_scale = 1.0;
    double bslmm_rho_scale = 1.0;
    double bslmm_p_scale = 1.0;
    int bslmm_mh_steps = 0;
    int bslmm_s_min = 0;
    int bslmm_rpace = 10;
    int bslmm_wpace = 1000;
    int bslmm_burnin = 100000;
    int bslmm_sampling = 1000000;
    long long bslmm_seed = 0;
    std::string scan_chr;
    long long scan_from_bp = -1;
    long long scan_to_bp = -1;
    bool force_slow = false;
    bool no_loco = false;
    bool use_covariate_projection = true;
    bool repeatability_mixed_precision = false;
    bool repeatability_refit_validation = false;
    int repeatability_refit_limit = 0;
    int grm_alg = 1;
    double min_maf = 0.01;
    double min_hwe = 0.0;
    std::string lmm_null_solver = "brent";
    std::string extract_file;
    std::string exclude_file;
    double sparse_cutoff = 0.05;
    std::string snp_weight_file;
    std::string pop_class_file;
    std::string repeatability_grm_cache_prefix;
    std::string repeatability_weight_file;
    std::string repeatability_mode = "auto";  // auto, exact, fast, stcg
    std::string vce_refine_mode = "auto";     // auto, full-record-exact, subject-level-exact, fast
    bool out_ssf = false;
    int num_threads = 0;                       // 0 = auto (use OMP max)
    int reml_maxit = 50;
    double reml_tol = 1e-6;
    double reml_vg = -1.0;
    double reml_vpe = -1.0;
    double reml_ve = -1.0;

    // Model-specific options (parsed in CLI but forwarded to model Options)
    int maternal_id_col = 3;
    std::string link_type = "probit";
    std::string family_type = "poisson";
    std::string offset_file;
    std::string binomial_file;
    bool ordinal_model = false;
    std::string ordinal_file;
    std::string ordinal_mode = "experimental-cumulative-binary";
    std::string ordinal_stacked_engine = "fixed-gls";
    std::string ordinal_vc_method = "fixed";
    double ordinal_sigma2_g = 0.5;
    double ordinal_vc_grid_min = 0.0;
    double ordinal_vc_grid_max = 2.0;
    int ordinal_vc_grid_points = 21;
    bool ordinal_allow_fallback = false;
    std::string ordinal_test = "wald";
    bool psr_null_only = false;
    int psr_maxit = 50;
    double psr_tol = 1e-6;
    double psr_mu_eps = 1e-6;
    double psr_min_weight = 1e-8;
    double psr_max_weight = 1e8;
    int env_covar_col = -1;
    std::string genotyped_list;
    double h_tau = 1.0;
    double h_omega = 1.0;
    double blend_g = 0.05;
    bool write_ebv = false;
};
}
