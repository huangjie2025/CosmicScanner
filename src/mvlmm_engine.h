#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <Eigen/Sparse>
#include "lmm_engine.h" // For MvGwasResult struct (if moved) or we redefine
#include "eigen_vce.h"

namespace cosmic {

struct MvGwasResult {
    std::string chrom;
    std::string snp_id;
    long long pos;
    std::string a1;
    std::string a2;
    double freq;
    
    // Multivariate stats
    Eigen::VectorXd beta; // d x 1
    Eigen::MatrixXd se;   // d x d covariance matrix of beta? Or just diagonal SE? 
                          // Usually we store diagonal SEs and maybe the full cov matrix if needed.
                          // For summary stats, usually just per-trait beta/se + p-values.
    
    double wald_stat;    // Multivariate Wald test statistic
    double p_value;      // P-value (Chi-squared df=n_traits)
    // double p_value_lrt;
};

// Multivariate LMM Engine for GWAS
// Focuses exclusively on genome-wide scanning using pre-computed Variance Components
class MvLmmEngine {
public:
    struct Options {
        int thread_num = 1;
        bool use_projection = true; // Use covariate projection for faster scan
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
              const std::string& chrom_filter = "") {
        throw std::runtime_error("MvLmmEngine::scan(BgenReader) not implemented");
    }

private:
    void processBlock(const Eigen::MatrixXd& S_block, 
                      const std::vector<double>& means, 
                      size_t global_start_idx, 
                      std::vector<MvGwasResult>& results,
                      const std::vector<struct PlinkBimInfo>& bim,
                      const Eigen::MatrixXd& H_XX_inv,
                      const Eigen::VectorXd& RHS_X,
                      const std::vector<Eigen::MatrixXd>& Sigma_inv);

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
