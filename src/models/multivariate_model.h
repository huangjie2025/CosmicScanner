#pragma once
#include <string>
#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "mvlmm_engine.h"
#include "data_filter.h"

namespace cosmic {

class MultivariateModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_mv";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        std::string scan_tier;
        bool use_covariate_projection = true;
        bool use_emmv = true;
        bool use_gemma_refit = false;
        bool use_gemma_full = false;
        bool use_gemma_exact = false;
        bool use_gemma_gxe = false;
        bool use_hybrid = false;
        bool use_eb = false;
        bool test_wald = true;
        bool test_lrt = true;
        bool test_score = true;
        int gemma_mode = 4;
        double pnr = 0.001;
        int gemma_emi = 10000;
        int gemma_nri = 100;
        double gemma_emp = 0.0001;
        double gemma_nrp = 0.0001;
        bool gemma_crt = false;
        int refine_top = 100;
        double refine_p = 1e-3;
        int refine_discrepancy_top = 0;
        double refine_discrepancy = 1.0;
        int refine_score_prepass = 0;
        double refit_p = 1.0;
        int reml_maxit = 1000;
        double reml_tol = 1e-8;
        bool use_aireml = false;
        bool use_emmvgwas_lambda = false;
        bool output_vbeta = false;
        bool output_emmvgwas = false;
        bool emmvgwas_vc_only = false;
        std::string emmvgwas_G_file;
        std::string emmvgwas_R_file;
        int em_min_iter = 10;
        int ai_max_iter = 50;
        double ai_tol = 1e-8;
        // Data source fields (migrated from GwasConfig)
        std::string pheno_file;
        std::vector<int> mpheno_list;
        std::string qcovar_file;
        std::string covar_file;
        int env_covar_col = -1;
        std::string scan_chrom;
        long long scan_from_bp = -1;
        long long scan_to_bp = -1;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;

        static Options fromConfig(const GwasConfig& cfg);
    };

    MultivariateModel(const Options& opts);

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
    Options options_;

    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;

    Eigen::MatrixXd Y_;
    Eigen::MatrixXd X_;
    Eigen::VectorXd env_;
    DataFilter::KeepInfo keep_;
    int n_keep_ = 0;
    std::vector<char> refine_mask_;
    const Eigen::MatrixXd* refit_grm_override_ = nullptr;
    std::string refit_chrom_filter_;
    bool refit_append_mode_ = false;
    double eb_effective_tests_ = 0.0;
    long long eb_tested_markers_ = 0;

    void loadPhenotype();
    void loadCovariates();
    void runNoLoco();
    void runLoco();
    void runGemmaRefit();
    void runGemmaGxe();
    void runHybrid();
    void writeResults(const std::string& out_file, const std::vector<MvGwasResult>& results);
    void writeEmmvgwasOutputs(const MvEigenVce& vce, const std::vector<MvGwasResult>& results) const;
    void writeEbSummary(const std::string& out_file) const;
};

} // namespace cosmic
