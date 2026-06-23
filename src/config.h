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
    int mpheno = 1;
    int tpheno = 2;
    int rrm_order = 3;
    std::vector<int> mpheno_list;
    bool multivariate = false;
    bool repeatability_model = false;
    bool threshold_model = false;
    bool count_model = false;
    bool rrm_model = false;
    bool make_grm = false;
    bool make_grm_sparse = false;
    bool run_gwas = false;
    bool fast_lmm = false;
    bool use_spa = false;
    bool run_skat = false;
    bool run_burden = false;
    std::string set_file = "";
    bool biobank_mode = false;
    int block_size = 1000;
    bool force_slow = false;
    bool no_loco = false;
    bool use_covariate_projection = true;
    bool repeatability_mixed_precision = false;
    bool repeatability_refit_validation = false;
    int repeatability_refit_limit = 0;
    int grm_alg = 1;
    double min_maf = 0.01;
    double min_hwe = 0.0;
    double sparse_cutoff = 0.05;
    std::string snp_weight_file;
    std::string pop_class_file;
    std::string repeatability_grm_cache_prefix;
    std::string repeatability_mode = "auto";  // auto, exact, fast, stcg
    bool out_ssf = false;
    int reml_maxit = 50;
    double reml_tol = 1e-6;
    double reml_vg = -1.0;
    double reml_vpe = -1.0;
    double reml_ve = -1.0;
    std::string link_type = "probit";
    std::string family_type = "poisson";
    std::string offset_file;
    std::string binomial_file;
    int psr_maxit = 50;
    double psr_tol = 1e-6;
    double psr_mu_eps = 1e-6;
    double psr_min_weight = 1e-8;
    double psr_max_weight = 1e8;
};
}
