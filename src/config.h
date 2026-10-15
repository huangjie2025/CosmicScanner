#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace cosmic {
struct GwasConfig {
    bool pca = false, ld_prune = false, qc_requested = false;
    bool skat_o = false, acat_v = false;
    bool farmcpu = false, blink = false, firth_logistic = false, saige_style = false;
    std::string logistic_test;
    bool logistic_glmm = false;
    bool conditional_joint = false;
    int farmcpu_max_iterations = 5, farmcpu_pseudo_qtn = 10;
    std::string condition_snp_file;
    int pca_k = 10, ld_window = 250000, ld_step = 1;
    double ld_r2 = 0.2, clump_p = 5e-8, geno = 1.0, mind = 1.0;
    std::string ld_clump_file, keep_file, remove_file;
    std::string snp_model = "additive";
    std::string gwas_model = "additive";  // additive | additive_dominance | dominance
    // Phase E2: block-streamed vanRaden G construction (no full n_g^2 G).
    // "auto" (default) = dense for small n_g / streamed for large; "dense" always
    // builds the full G; "streamed" forces the block-streamed path (APY methods only).
    std::string ssgwas_grm_mode = "auto";
    // Phase E3-B1: Gcc solver backend.
    std::string ssgwas_gcc_solver = "dense";  // "dense" (default) | "pcg"
    std::string ssgwas_apy_inverse_mode = "dense";  // "dense" (default) | "operator"
    std::string ssgwas_mme_apy_mode = "dense";  // "dense" (default) | "operator"
    std::string ssgwas_pev_apy_mode = "dense";  // "dense" (default) | "operator"
    std::string ssgwas_hinv_mode = "dense";     // "dense" (default) | "operator"
    // Phase E3-B2: PCG preconditioner.
    std::string ssgwas_pcg_preconditioner = "jacobi";  // "jacobi" | "lowrank"
    int ssgwas_pcg_preconditioner_rank = 50;
    double ssgwas_apy_core_variance = 0.99;    // automatic core-size target: C98=0.98, C99=0.99 (when core_size = -1)
    std::string ssgwas_apy_core_eigensolver = "auto";  // "auto" | "exact" | "lanczos"
    std::string rare_method;
    std::string bfile_prefix;
    std::string bgen_file;
    std::string sample_file;
    std::string pfile_prefix;
    std::string pgen_file;
    std::string pvar_file;
    std::string psam_file;
    std::string vcf_file;
    std::string vcf_tbi_file;
    std::string vcf_sample_file;
    std::string hmp_file;
    std::string geno_matrix_file;
    std::string geno_var_file;
    std::string geno_sample_file;
    std::string geno_orientation = "variants-by-samples";
    std::string npy_file;
    std::string npy_var_file;
    std::string npy_sample_file;
    std::string npy_orientation = "variants-by-samples";
    std::string ped_file;
    std::string map_file;
    std::string convert_from;
    std::string convert_to;
    std::string convert_out;
    std::string batch_models;
    std::string batch_out_prefix;
    std::string out_prefix = "cosmic_gwas";
    std::string pheno_file;
    std::string long_pheno_file;
    std::string common_env_file;
    bool common_env_record_level = false;
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
    std::string rrm_gwas_method = "coefficient"; // record | coefficient (unbalanced RRM-GWAS)
    std::string rrm_compatible;
    std::string rrm_compatible_varcom;
    std::vector<int> mpheno_list;
    std::string tpoint_col_name;
    std::string tclass_names;
    std::string tcovar_names;
    bool multivariate = false;
    std::string multivariate_method = "emmv"; // emmv, wald, gemma, gemma-full, gemma-exact, gemma-gxe, gemma-gxe-full, hybrid
    bool multivariate_method_explicit = false;
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
    bool common_env_model = false;
    bool single_step_model = false;
    bool animal_model = false;
    bool write_ped_kernel = false;
    bool make_agh = false;
    bool write_agh = false;
    bool threshold_model = false;
    bool count_model = false;
    bool survival_model = false;
    bool gxe_model = false;
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
    bool standalone_reml = false;
    std::string reml_method = "ai";
    int ci_mode = 0; // GEMMA -ci: summary-stat confidence interval route
    bool run_gwas = false;
    std::string scan_tier; // approx=fixed null, semi=fixed ratio/profile scale, exact=per-SNP VC refit
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
    int lmm_gemma_mode = 0; // GEMMA -lmm: 1 Wald, 2 LRT, 3 Score, 4 all; 0 disabled
    bool lmm_gemma_refit = false;
    bool linear_model = false;
    int linear_model_mode = 4; // GEMMA -lm: 1 Wald, 2 LRT, 3 Score, 4 all
    bool eb_lmm = false;
    bool use_spa = false;
    bool run_skat = false;
    bool run_burden = false;
    double rare_max_maf = 0.5;
    int rare_min_mac = 1;
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
    int bslmm_thin = 10;
    long long bslmm_seed = 0;
    bool bslmm_mcmc = false;
    bool bvsr = false;
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
    double bgen_dom_max_missing = 0.05;
    double bgen_dom_min_info = 0.01;
    double bgen_dom_entropy_warning = 0.50;
    std::string lmm_null_solver = "brent";
    std::string extract_file;
    std::string exclude_file;
    double sparse_cutoff = 0.05;
    std::string snp_weight_file;
    std::string pop_class_file;

    // COJO-style conditional joint analysis
    bool cojo = false;
    bool cojo_slct = false;
    bool cojo_joint = false;
    std::string cojo_summary;
    std::string cojo_ld_ref;
    double cojo_p = 5e-8;
    double cojo_r2 = 0.9;
    int cojo_window = 10000000;
    int cojo_max_iter = 20;
    std::string cojo_cond_file;

    // Native summary-statistic fine-mapping
    bool fine_map = false;
    std::string finemap_summary;
    std::string finemap_region;
    double credible_level = 0.95;
    double finemap_prior_var = 0.04;
    int finemap_max_signals = 5;
    double finemap_signal_p = 5e-8;

    // Meta-analysis
    bool meta_analysis = false;
    std::string meta_studies_file;
    std::vector<std::string> meta_study_files;
    std::string meta_method = "fixed"; // fixed, random-dl, random-reml
    double meta_het_p = 0.05;
    int meta_min_studies = 2;
    bool meta_gc = false;
    double meta_gc_lambda = 1.0;
    std::string meta_out_format = "metal"; // metal, gwama, simple
    std::string repeatability_grm_cache_prefix;
    std::string repeatability_weight_file;
    std::string repeatability_mode = "auto";  // auto, exact, fast, stcg
    std::string vce_refine_mode = "auto";     // auto, full-record-exact, subject-level-exact, fast
    bool out_ssf = false;
    bool out_gz = false;
    bool manhattan_tsv = false;
    bool qq_tsv = false;
    int num_threads = 0;                       // 0 = auto (use OMP max)
    bool memory_budget_specified = false;      // SolverPlanner diagnostic budget
    std::uint64_t memory_budget_bytes = 0;     // Internal unit: bytes
    std::string memory_budget_source = "UNKNOWN";
    int reml_maxit = 50;
    double reml_tol = 1e-6;
    double reml_vg = -1.0;
    double reml_vpe = -1.0;
    double reml_vc = -1.0;
    double reml_ve = -1.0;

    // Model-specific options (parsed in CLI but forwarded to model Options)
    int maternal_id_col = 3;
    std::string maternal_snp_test = "direct";
    std::string link_type = "probit";
    std::string family_type = "poisson";
    std::string offset_file;
    // Opt-in: allow --threshold-model to dichotomise a non-binary response at 0.5.
    bool threshold_dichotomize = false;
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

    // --- ssGWAS Method 2 (ssGBLUP-derived back-solving), Phase 6 ---
    std::string single_step_method = "lmm";  // "lmm" (default, Method 1) | "backsolve" (Method 2)
    int ssgwas_window_snps = 10;             // sliding SNP-count window for %var aggregation
    bool ssgwas_tune_g = false;              // align G to A22 (G = b*G + a, preGSf90-style)
    bool ssgwas_blend_a22 = false;           // blend G* with A22 instead of I (HIBLUP/preGSf90 reference); backsolve only
    std::string ssgwas_solver = "ldlt";      // "ldlt" (reference) | "pcg" (matrix-free, large n)
    std::string ssgwas_vce;                  // VCE mode ("" = fixed --reml-vg/--reml-ve path; "ai"/"emai"/...)
    double ssgwas_pcg_tol = 1e-8;
    int ssgwas_pcg_maxit = 2000;
    // Phase C1: frequentist P-values for back-solved SNP effects (Aguilar 2019).
    // "off" (default) | "exact" (reference; needs LDLT factorization; exact selected inverse)
    // | "apy-exact" (Phase C2: APY G^-1 in MME + exact core-slice selected inverse;
    //   needs LDLT) | "apy-approx" (Phase D: APY G^-1 + Eq. 8/9 approximate core
    //   PEV, Schur form; no selected inverse / no dense n_g^2 uncertainty block;
    //   LDLT or PCG).
    std::string ssgwas_pvalue_method = "off";
    // Explicit APY-exact SNP estimator; default preserves historical outputs.
    std::string ssgwas_apy_exact_estimator = "legacy-hybrid-v1";
    bool ssgwas_apy_exact_estimator_explicit = false;
    std::string ssgwas_pev_scan_mode = "batched"; // "scalar" | "batched"
    bool ssgwas_pev_scan_mode_explicit = false; // exact/core-v2 opt in; preserve historical omission
    int ssgwas_pev_scan_block_size = 256;
    bool allow_large_exact_ssgwas_pvalue = false;  // bypass the 8 GiB selected-inverse guard
    // Phase C2: APY core partition (Leite et al. 2024).
    int ssgwas_apy_core_size = 0;              // core animals (must be >= 1 for apy-* modes)
    std::string ssgwas_apy_core_file;          // optional core-animal ID file (overrides random)
    int ssgwas_apy_seed = 12345;               // deterministic seed for the random core selection
    double ssgwas_apy_core_target_ratio = 0.0; // harness label: core_size / K98 (0 = unset)
    std::vector<std::string> command_argv;
};
}
