#pragma once
#include <Eigen/Dense>
#include <vector>
#include <limits>

namespace cosmic {

// Exact VCE for Dominance-Only (d+D) Model using Eigen Decomposition
//
// Model: y = X*beta + u_D + e
//   u_D ~ N(0, sigma_D^2 * K_D)
//   e   ~ N(0, sigma_e^2 * I)
//
// Profile REML on single variance ratio eta_d = sigma_D^2 / sigma_e^2
// Likelihood definition matches EigenVceAD::profileLogL and VCESolver::solveDiagonal exactly.
class EigenVceD {
public:
    struct Options {
        bool use_covariate_projection = true;
        int max_em_iter = 50;
        double tol = 1e-6;
        bool verbose = true;
    };

    EigenVceD();
    EigenVceD(const Options& opts);

    void prepare(const Eigen::VectorXd& y,
                 const Eigen::MatrixXf& X,
                 const Eigen::MatrixXd& K_D);

    bool runNullModel();

    void fixedComponents(double s2d, double s2e);

    bool converged() const { return converged_; }
    int lastIterations() const { return iterations_; }

    double getSigma2d() const { return sigma2_d; }
    double getSigma2e() const { return sigma2_e; }
    double getLogLikelihood() const { return logL_; }

    const Eigen::VectorXd& getEigenvalues() const { return D_V; }
    const Eigen::MatrixXd& getProjectedW() const { return W_proj; }
    const Eigen::VectorXd& getUty() const { return Uty; }
    const Eigen::MatrixXd& getUtX() const { return UtX; }

private:
    Options options;
    int n_samples;
    int n_covars;
    int m;

    Eigen::MatrixXd U_D;
    Eigen::VectorXd D_D;

    Eigen::VectorXd Uty_prepared;
    Eigen::VectorXd Uty;
    Eigen::MatrixXd A_proj;
    Eigen::MatrixXd UtX;      // m x n_covars in the projected space

    Eigen::VectorXd D_V;
    Eigen::MatrixXd W_proj;

    double sigma2_d = 0.0;
    double sigma2_e = 0.0;
    double logL_ = -std::numeric_limits<double>::infinity();
    bool converged_ = false;
    int iterations_ = 0;

    void computeFinalDecomposition();
};

} // namespace cosmic