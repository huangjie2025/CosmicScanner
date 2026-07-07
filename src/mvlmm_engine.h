#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <limits>
#include <unordered_set>
#include <Eigen/Sparse>
#include "lmm_engine.h" // For MvGwasResult struct (if moved) or we redefine
#include "eigen_vce.h"

namespace cosmic {

struct MvGwasResult {
    std::string chrom;
    std::string snp_id;
    long long pos = 0;
    std::string a1;
    std::string a2;
    double freq = std::numeric_limits<double>::quiet_NaN();
    
    // Multivariate stats
    Eigen::VectorXd beta; // d x 1
    Eigen::MatrixXd se;   // d x d covariance matrix of beta? Or just diagonal SE? 
                          // Usually we store diagonal SEs and maybe the full cov matrix if needed.
                          // For summary stats, usually just per-trait beta/se + p-values.
    
    double wald_stat = std::numeric_limits<double>::quiet_NaN();    // Multivariate Wald test statistic
    double p_value = std::numeric_limits<double>::quiet_NaN();      // P-value (Chi-squared df=n_traits)
    double p_wald = std::numeric_limits<double>::quiet_NaN();
    double p_lrt = std::numeric_limits<double>::quiet_NaN();
    double p_score = std::numeric_limits<double>::quiet_NaN();
    double p_wald_gc = std::numeric_limits<double>::quiet_NaN();
    double p_lrt_gc = std::numeric_limits<double>::quiet_NaN();
    double p_score_gc = std::numeric_limits<double>::quiet_NaN();
    double p_wald_crt = std::numeric_limits<double>::quiet_NaN();
    double p_lrt_crt = std::numeric_limits<double>::quiet_NaN();
    double p_score_crt = std::numeric_limits<double>::quiet_NaN();
    double eb_marker_lambda = 0.0;
    double eb_confidence = 0.0;
    Eigen::MatrixXd eb_prior_cov;
    Eigen::MatrixXd eb_post_cov;
};

// Multivariate LMM Engine for GWAS
// Focuses exclusively on genome-wide scanning using pre-computed Variance Components
class MvLmmEngine {
public:
    struct Options {
        int thread_num = 1;
        bool use_projection = true; // Use covariate projection for faster scan
        bool em_mvgwas_test = true; // Use EMmvGWAS semi-exact Wald test
        bool use_general_sigma_inverse = false;
        bool center_genotypes = true;
        bool use_dosage_minus_one = false;
        bool eb_lmm = false;
        std::string scan_chrom;
        long long scan_from_bp = -1;
        long long scan_to_bp = -1;
        std::unordered_set<std::string> extract_set;
        std::unordered_set<std::string> exclude_set;
    };

    MvLmmEngine(const Options& opts);

    // Prepare for Scan using an already solved MvEigenVce Null Model
    void prepareFromVce(const MvEigenVce& vce);

    // Scan (PGEN)
    void scan(class PgenReader& reader, 
              std::vector<MvGwasResult>& results, 
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");

    // Scan (PLINK)
    void scan(class PlinkReader& reader, 
              std::vector<MvGwasResult>& results, 
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");
              
    // Scan (BGEN)
    void scan(class BgenReader& reader, 
              std::vector<MvGwasResult>& results, 
              const std::vector<int>& keep_indices,
              const std::string& chrom_filter = "");

private:
    void processBlock(const Eigen::MatrixXd& S_block, 
                      const std::vector<double>& means, 
                      size_t global_start_idx, 
                      const std::vector<size_t>& global_indices,
                      std::vector<MvGwasResult>& results,
                      const std::vector<struct PlinkBimInfo>& bim,
                      const Eigen::MatrixXd& H_XX_inv,
                      const Eigen::VectorXd& RHS_X,
                      const std::vector<Eigen::MatrixXd>& Sigma_inv,
                      const Eigen::MatrixXd& Cyy_null);

    Options options;

    // Extracted from MvEigenVce
    Eigen::VectorXd D;
    Eigen::MatrixXd W_proj;
    Eigen::MatrixXd Uty;
    Eigen::MatrixXd UtX;
    Eigen::MatrixXd Vg;
    Eigen::MatrixXd Ve;

    int n_samples;
    int n_traits;
    int n_covars;
};

} // namespace cosmic
