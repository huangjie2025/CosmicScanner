#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"

namespace cosmic {

class UnivariateModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_gwas";
        bool out_ssf = false;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        bool fast_lmm = false;
        std::string scan_tier;
        bool profile_scale = false;
        int lmm_gemma_mode = 0;
        bool lmm_gemma_refit = false;
        bool linear_model = false;
        int linear_model_mode = 4;
        bool eb_lmm = false;
        int n_threads = 1;
        bool no_loco = false;
        int grm_alg = 1;
        EigenVce::NullSolver lmm_null_solver = EigenVce::NullSolver::Brent;
        std::string lmm_null_solver_name = "brent";
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
        std::string set_file;
        std::string extract_file;
        std::string exclude_file;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        int mpheno = 1;

        static Options fromConfig(const GwasConfig& cfg);
    };

    UnivariateModel(const Options& opts);

    void prepare(const std::vector<std::string>& fids,
                 const std::vector<std::string>& iids,
                 int n_samples,
                 PlinkReader& plink_reader,
                 BgenReader& bgen_reader,
                 PgenReader& pgen_reader,
                 bool is_bgen,
                 bool is_pgen);

    void run();

private:
    void runLinearModel();

    Options options_;
    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;
    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;
    DataFilter::KeepInfo keep_;
};

} // namespace cosmic
