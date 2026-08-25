#pragma once
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>
#include "config.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "data_filter.h"
#include "eigen_vce.h"
#include "multi_random_effect_lmm.h"
#include "design.h"

namespace cosmic {

/// Common-environment GWAS.
///
/// y = Xb + x_j beta_j + Z_a a + Z_c c + e
///
///   a ~ N(0, K sigma_a^2)   direct additive genetic (K = genomic GRM)
///   c ~ N(0, I_g sigma_c^2) common environment / litter / group
///   e ~ N(0, I sigma_e^2)
///
/// Z_a maps records to individuals, Z_c maps records to group levels.
/// This is NOT the permanent-environment effect: Z_c has a completely
/// different incidence structure (group blocks), not record->individual.
///
/// Status: Experimental (internal simulation validation only).
class CommonEnvironmentModel {
public:
    struct Options {
        std::string out_prefix = "cosmic_common_env";
        double min_maf = 0.01;
        double min_hwe = 0.0;
        int grm_alg = 1;
        int n_threads = 1;
        bool no_loco = false;
        int reml_maxit = 100;
        double reml_tol = 1e-6;
        int mpheno = 1;
        int tpheno = 0;          // time column for longitudinal pheno (0 = plain)
        double reml_vg = -1.0;   // debug: fix variance components (skip VCE)
        double reml_vc = -1.0;
        double reml_ve = -1.0;
        std::string pheno_file;
        std::string qcovar_file;
        std::string covar_file;
        std::string common_env_file;   // FID IID GROUP_ID (record or subject level)
        bool group_file_record_level = false; // GROUP column is record-level (FID IID TIME GROUP)
        std::string extract_file;
        std::string exclude_file;

        static Options fromConfig(const GwasConfig& cfg);
    };

    CommonEnvironmentModel(const Options& opts);

    void prepare(const std::vector<std::string>& fids,
                 const std::vector<std::string>& iids,
                 int n_samples,
                 PlinkReader& plink_reader,
                 BgenReader& bgen_reader,
                 PgenReader& pgen_reader,
                 bool is_bgen,
                 bool is_pgen);

    void run();

    // Exposed for tests / diagnostics
    const Options& options() const { return options_; }
    double v_g() const { return v_g_; }
    double v_c() const { return v_c_; }
    double v_e() const { return v_e_; }
    int n_records() const { return n_recs_kept_; }
    int n_subjects_kept() const { return n_subjects_kept_; }
    int n_groups() const { return n_groups_; }

private:
    Options options_;

    int n_samples_ = 0;
    std::vector<std::string> fids_, iids_;
    bool is_bgen_ = false, is_pgen_ = false;
    PlinkReader* plink_reader_ = nullptr;
    BgenReader* bgen_reader_ = nullptr;
    PgenReader* pgen_reader_ = nullptr;

    std::vector<GenRecord> gen_records_;   // kept records
    std::vector<int> obs_to_sample_;       // record -> subject index
    std::vector<int> obs_to_group_;        // record -> group level
    std::vector<std::string> group_labels_; // group level -> label

    Eigen::VectorXd y_obs_;
    Eigen::MatrixXd X_design_;             // n_records x p (double)
    FixedDesignG fd_;
    Eigen::MatrixXd G_subject_;

    int n_recs_kept_ = 0;
    int n_subjects_kept_ = 0;
    int n_groups_ = 0;

    double v_g_ = 0.0, v_c_ = 0.0, v_e_ = 0.0;
    bool vce_converged_ = false;
    int vce_iters_ = 0;
    double vce_best_logL_ = -std::numeric_limits<double>::infinity();

    void loadPhenotype();
    void loadCommonEnvMapping();
    void loadCovariates();
    void buildGRM();

    RandomEffectSpec additiveSpec() const;
    RandomEffectSpec commonEnvSpec() const;

    // REML logL of the record-level 3-component model evaluated at
    // (vg, vc, ve) with the given subject GRM. Returns finite value on success.
    double recordLogL(const Eigen::MatrixXd& G_used, double vg, double vc, double ve);
    bool runVCE(const Eigen::MatrixXd& G_used);

    // Fixed-null scan on one GRM, writing to assoc_file.
    void scanWithGrm(const Eigen::MatrixXd& G_used, const std::string& assoc_file, const std::string& chrom_filter);

    void writeNullSummary(const std::string& filename, const std::string& scan_mode) const;
};

} // namespace cosmic