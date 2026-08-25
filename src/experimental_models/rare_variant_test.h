#ifndef COSMIC_SCANNER_RARE_VARIANT_TEST_H
#define COSMIC_SCANNER_RARE_VARIANT_TEST_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

namespace cosmic {

// Structure to define a set of SNPs (e.g., a Gene or a regulatory region)
struct VariantSet {
    std::string set_id; // e.g., "GeneA"
    std::vector<std::string> snp_ids; // List of RSIDs belonging to this set
    std::vector<double> weights;      // Optional weights (e.g., from beta(MAF, 1, 25))
};

// Parser for Set Files (Format: SetID \t SNP_ID)
class SetFileReader {
public:
    static std::vector<VariantSet> read(const std::string& filename);
};

// Rare Variant Association Test (SKAT & Burden)
// References:
// 1. Burden Test: Li and Leal (2008), Madsen and Browning (2009)
// 2. SKAT: Wu et al. (2011)
class RareVariantTest {
public:
    enum class Method {
        BURDEN,
        SKAT,
        SKAT_O // Optimal linear combination of Burden and SKAT
    };

    struct Result {
        std::string set_id;
        int n_snps = 0;        // Number of valid SNPs in the set
        double p_value = 1.0;
        double stat = 0.0;     // Test statistic (e.g., Q for SKAT)
    };

    RareVariantTest(Method method = Method::SKAT);
    
    // Set parameters for Beta weight density function (default: a=1, b=25 as in SKAT paper)
    void setBetaWeights(double a = 1.0, double b = 25.0);

    // X_set: N x K matrix of genotypes for K SNPs in the set
    // score_vec: K x 1 vector of marginal scores (Z'P_y)
    // P_matrix: N x N projection matrix from the null model (P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1})
    // MAF: K x 1 vector of Minor Allele Frequencies (for calculating weights if not provided)
    Result runTest(const std::string& set_id, 
                   const Eigen::MatrixXd& X_set, 
                   const Eigen::VectorXd& score_vec, 
                   const Eigen::MatrixXd& P_matrix, 
                   const Eigen::VectorXd& MAF);

private:
    Method method_;
    double beta_a_;
    double beta_b_;

    // Calculate Beta distribution weights based on MAF
    Eigen::VectorXd calculateWeights(const Eigen::VectorXd& MAF) const;

    // Run Burden Test
    Result runBurden(const std::string& set_id, const Eigen::MatrixXd& X_set, const Eigen::VectorXd& score_vec, const Eigen::MatrixXd& P_matrix, const Eigen::VectorXd& W) const;

    // Run standard SKAT
    Result runSKAT(const std::string& set_id, const Eigen::MatrixXd& X_set, const Eigen::VectorXd& score_vec, const Eigen::MatrixXd& P_matrix, const Eigen::VectorXd& W) const;

    // Helper: Davies method or Kuonen method to compute p-value from mixture of Chi-squares
    // For skeleton, we use a Satterthwaite approximation
    double mixtureChiSqPval(double Q, const Eigen::VectorXd& lambda) const;
};

} // namespace cosmic

#endif // COSMIC_SCANNER_RARE_VARIANT_TEST_H