#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <limits>
#include <vector>
#include <functional>
#include "kinship_operator.h"

namespace cosmic {

// Matrix-Free VCE for Single Trait using KinshipOperator
// Implements REML variance component estimation without materializing dense GRM
class MatrixFreeVce {
public:
    enum class NullSolver {
        Brent,
        IDUL,
        IDULPlus,
        AIDUL,
        EM,
        AI
    };

    struct Options {
        bool use_covariate_projection = true;
        NullSolver null_solver = NullSolver::Brent;
        int max_iter = 100;
        double tol = 1e-6;
        // Randomized low-rank decomposition controls
        int rank = 256;
        int oversample = 10;
        long long seed = 42;
        double trace_threshold = 0.0;
        int min_rank = 1;
        int n_threads = 1;
    };

    MatrixFreeVce();
    MatrixFreeVce(const Options& opts);

    // Prepare using matrix-free kinship operator
    // kinship: matrix-free kinship operator (K = Z*Z'/M)
    // y: phenotype vector (n)
    // X: covariate matrix (n x c) as float
    void prepare(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const class KinshipOperator& kinship);

    // Run null model variance component estimation
    bool runNullModel();

    // Skip VCE and use pre-computed variance components directly
    void fixedComponents(double s2g, double s2e);

    bool converged() const { return converged_; }
    int lastIterations() const { return iterations_; }
    double getSigma2g() const { return sigma2_g_; }
    double getSigma2e() const { return sigma2_e_; }
    double getDelta() const { return null_delta_; }
    double getSigma2gSE() const { return sigma2_g_se_; }
    double getSigma2eSE() const { return sigma2_e_se_; }
    double getLogLikelihood() const { return log_likelihood_; }
    double getMaxScore() const { return max_score_; }

    // Accessors for low-rank results
    const Eigen::VectorXd& getEigenvalues() const { return D_; }
    const Eigen::MatrixXd& getEigenvectors() const { return U_; }
    const Eigen::MatrixXd& getProjectedW() const { return W_proj_; }
    const Eigen::VectorXd& getUty() const { return Uty_; }
    const Eigen::MatrixXd& getUtX() const { return UtX_; }
    double getResidualEigenvalue() const { return low_rank_residual_eigenvalue_; }

private:
    Options options_;
    int n_samples_ = 0;
    int n_covars_ = 0;

    // Low-rank decomposition results
    Eigen::VectorXd D_;           // eigenvalues of K
    Eigen::MatrixXd U_;           // eigenvectors of K (n x k)
    Eigen::MatrixXd W_proj_;      // projected W (A*U or U)
    Eigen::VectorXd Uty_;         // U' * y
    Eigen::MatrixXd UtX_;         // U' * X

    const class KinshipOperator* kinship_ptr_ = nullptr;

    bool converged_ = false;
    int iterations_ = 0;
    double null_delta_ = 1.0;
    double sigma2_g_ = 0.0;
    double sigma2_e_ = 0.0;
    double sigma2_g_se_ = std::numeric_limits<double>::quiet_NaN();
    double sigma2_e_se_ = std::numeric_limits<double>::quiet_NaN();
    double log_likelihood_ = std::numeric_limits<double>::quiet_NaN();
    double max_score_ = std::numeric_limits<double>::quiet_NaN();
    double low_rank_residual_eigenvalue_ = 0.0;

    // Stored phenotype and covariates for REML optimization
    Eigen::VectorXd y_;
    Eigen::MatrixXf X_;

    // Matrix-free operations
    void prepareMatrixFree(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const class KinshipOperator& kinship);

    // Matrix-free null model solvers
    bool runIterativeReml(bool average_information);
    bool runIdulNullModel(bool monotone, bool accelerated);
    bool runBrentNullModel();

    // Matrix-free operations using KinshipOperator
    void applyKinship(const Eigen::VectorXd& v, Eigen::VectorXd& out, const class KinshipOperator& kinship, int n_threads) const;
    void applyKinship(const Eigen::MatrixXd& X, Eigen::MatrixXd& out, const class KinshipOperator& kinship, int n_threads) const;

    // Matrix-free trace estimation
    double estimateTraceInverse(int n, int k, const class KinshipOperator& kinship, int n_threads) const;

    void updateRemlDiagnostics();
};

} // namespace cosmic