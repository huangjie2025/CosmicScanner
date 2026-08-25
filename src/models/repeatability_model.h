#pragma once
#include <string>
#include <vector>
#include <memory>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "lmm_engine.h"
#include "eigen_vce.h"
#include "design.h"
#include "mme_builder.h"
#include "matrix_adapter.h"

namespace cosmic {

/// Repeatability model strategy (IASBLUP-style)
enum class RepeatabilityMode {
    Auto,   ///< Auto-select based on data size
    Exact,  ///< R1: Record-level exact GLS (medium n, high accuracy)
    Fast,   ///< R2: Aggregated mean approximation + Eigen null
    STCG    ///< R3: Matrix-free STCG (Biobank-scale longitudinal)
};

/// VCE refinement mode for repeatability model
enum class VceRefineMode {
    Auto,              ///< Default: subject-level Woodbury REML; full-record exact only when explicitly requested
    FullRecordExact,   ///< Force full record-level exact REML refinement
    SubjectLevelExact, ///< Force subject-level exact REML refinement (Woodbury-based)
    Fast               ///< Skip exact REML refinement, use EM results only
};

class RepeatabilityModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_repeatability";
        bool out_ssf = false;
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        bool fast_lmm = false;
        std::string scan_tier;
        bool eb_lmm = false;
        bool mixed_precision = false;
        bool refit_validation = false;
        int refit_limit = 0;
        std::string grm_cache_prefix;
        std::string weight_file;
        std::string extract_file;
        std::string exclude_file;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;
        bool use_covariate_projection = true;
        EigenVce::NullSolver lmm_null_solver = EigenVce::NullSolver::Brent;
        int reml_maxit = 50;
        double reml_tol = 1e-6;
        double reml_vg = -1.0;
        double reml_vpe = -1.0;
        double reml_ve = -1.0;
        // Data source fields (migrated from GwasConfig)
        std::string pheno_file;
        int mpheno = 1;
        int tpheno = 0;
        std::string qcovar_file;
        std::string covar_file;
        bool use_spa = false;

        // IASBLUP-style strategy selection
        RepeatabilityMode repeatability_mode = RepeatabilityMode::Auto;

        // VCE refinement mode
        VceRefineMode vce_refine_mode = VceRefineMode::Auto;

        static Options fromConfig(const GwasConfig& cfg);
    };

    RepeatabilityModel(const Options& opts);

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

    std::vector<GenRecord> gen_records_;
    FixedDesignG fd_;
    std::vector<RandomComponent> components_;
    Eigen::VectorXd y_obs_;
    Eigen::VectorXd record_weights_;
    Eigen::MatrixXf X_obs_;
    std::vector<int> obs_to_sample_;
    Eigen::MatrixXd G_subject_;
    Eigen::MatrixXd Ginv_;
    Eigen::SparseMatrix<double> I_mat_;
    std::unique_ptr<DenseMatrixAdapter> Ginv_adapter_;
    std::unique_ptr<SparseMatrixAdapter> I_adapter_;

    double v_g_ = 0.0, v_pe_ = 0.0, v_e_ = 0.0;
    int n_recs_kept_ = 0;
    double grm_time_sec_ = 0.0;
    double vce_time_sec_ = 0.0;
    double scan_vbuild_time_sec_ = 0.0;
    double scan_snp_time_sec_ = 0.0;

    void loadPhenotype();
    void loadCovariates();
    void buildGRM();
    void buildMMEComponents();
    void runVCE();
    void scanFast();
    void scanExact();
    void scanSTCG();
    void scanRefitValidation();

    void writeNullSummary(const std::string& filename, const std::string& scan_mode,
                          const std::vector<std::pair<std::string, std::string>>& extra_fields = {}) const;
};

} // namespace cosmic
