#pragma once
#include <functional>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <utility> // Added for std::pair
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "eigen_vce.h" // Added for EigenVce dependency
#include "spa_test.h"  // Saddle Point Approximation

namespace cosmic {

class MatrixFreePCG;

// Structure to hold LMM scan results
struct GwasResult {
    std::string chrom;
    std::string snp_id;
    long long pos = 0;
    std::string a1; // Effect allele (Minor)
    std::string a2; // Reference allele (Major)
    double freq = 0.0;    // Frequency of A1
    double beta = 0.0;
    double se = 0.0;
    double p_value = std::numeric_limits<double>::quiet_NaN();
    double log10_p = 0.0;
    double wald_stat = std::numeric_limits<double>::quiet_NaN();
    double lrt_stat = std::numeric_limits<double>::quiet_NaN();
    double score_stat = std::numeric_limits<double>::quiet_NaN();
    double p_wald = std::numeric_limits<double>::quiet_NaN();
    double p_lrt = std::numeric_limits<double>::quiet_NaN();
    double p_score = std::numeric_limits<double>::quiet_NaN();
    int n_obs = 0;
    double eb_marker_lambda = 0.0;
    double eb_prior_var = 0.0;
    double eb_post_var = std::numeric_limits<double>::quiet_NaN();
    double eb_confidence = 0.0;
    double eb_wald = 0.0;
    double eb_lrt = 0.0;
    double eb_lrt_p = 1.0;
    double eb_sigma2 = 0.0;
    int eb_conv = 0;
    // GxE model fields (optional)
    double beta_gxe = 0.0;
    double se_gxe = 0.0;
    double joint_chisq = 0.0;
};

// LMM Engine for GWAS
// Focuses exclusively on genome-wide scanning using pre-computed Variance Components
class LmmEngine {
public:
    struct Options {
        int thread_num = 1;
        double min_maf = 0.01;
        double min_hwe = 0.0;   // Minimum HWE p-value filter (0.0 to disable)
        bool fast_mode = false; // If true, use fixed delta from null model
        bool profile_scale = false; // Fixed null variance ratio, per-SNP residual-scale profile
        bool gemma_refit = false; // If true, profile alternative VCE per SNP for GEMMA LMM Wald/LRT
        bool use_residuals = false; // If true, use residuals from null model for scan (GCTA style)
        bool use_covariate_projection = true; // If true, project out X before Eigen Decomposition (accelerates Exact LMM)
        bool use_spa = false;   // Saddle Point Approximation for binary traits
        Eigen::VectorXd obs_weights; // Optional weights for observations (e.g., D^{-1/2} for repeatability model)
        std::unordered_set<std::string> extract_set; // If not empty, only test these SNP IDs
        std::unordered_set<std::string> exclude_set; // If not empty, exclude these SNP IDs
        int64_t scan_from_bp = -1;
        int64_t scan_to_bp = -1;
    };

    LmmEngine();
    LmmEngine(const Options& opts);

    // 1. Prepare for Scan using an already solved EigenVce Null Model
    void prepareFromVce(const EigenVce& vce);

    // Accessors for Null Model Results
    double getDelta() const { return null_delta; }
    double getSigma2g() const { return sigma2_g; }
    double getSigma2e() const { return sigma2_e; }

    // Load Sparse GRM
    void loadSparseGRM(const std::string& prefix);

    Eigen::SparseMatrix<double> sparse_grm;
    bool has_sparse_grm = false;

    // 2. Scan: Run GWAS on all SNPs in the reader
    // keep_indices: Indices of samples to use (must match prepare's y and X)
    // chrom_filter: If not empty, only scan SNPs on this chromosome
    void scan(const PlinkReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");
    void scanEB(const PlinkReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");

    // Scan BGEN
    void scan(BgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");
    void scanEB(BgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");

    // Scan PGEN
    void scan(PgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");
    void scanEB(PgenReader& reader, std::vector<GwasResult>& results, const std::vector<int>& keep_indices, const std::string& chrom_filter = "");

    // Biobank Mode (Matrix-Free PCG)
    // Prepare with pre-computed VCs (from RHE)
    // Supports optional covariate matrix X (which will be projected out)
    void prepareBiobank(const Eigen::VectorXd& y, 
                        const Eigen::MatrixXf& X, 
                        double s2g, double s2e);
                        
    // Run PCG-based Scan
    void scanBiobank(PlinkReader& reader, 
                     std::vector<GwasResult>& results, 
                     const std::vector<int>& keep_indices,
                     const std::string& chrom_filter = "");
                     
    void scanBiobank(PgenReader& reader, 
                     std::vector<GwasResult>& results, 
                     const std::vector<int>& keep_indices,
                     const std::string& chrom_filter = "");

    void scanBiobank(BgenReader& reader, 
                     std::vector<GwasResult>& results, 
                     const std::vector<int>& keep_indices,
                     const std::string& chrom_filter = "");

private:
    Options options;
    
    // Extracted from EigenVce
    Eigen::VectorXd D;
    Eigen::MatrixXd W_proj; // Rotation matrix (either U or A*U depending on projection)
    Eigen::VectorXd Uty;
    Eigen::MatrixXd UtX;
    
    double sigma2_g = 0.0;
    double sigma2_e = 0.0;
    double null_delta = 1.0;
    
    // Biobank Mode Variables
    Eigen::VectorXd y_raw;
    Eigen::MatrixXf X_raw; // Can also store Q_proj

    double estimateTraceInverse(MatrixFreePCG& pcg, int n, int k);
};

} // namespace cosmic
