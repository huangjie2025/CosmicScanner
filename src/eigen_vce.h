#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

namespace cosmic {

// Exact VCE for Single Trait using Eigen Decomposition
class EigenVce {
public:
    enum class NullSolver {
        Brent,
        IDUL,
        IDULPlus,
        AIDUL
    };

    struct Options {
        bool use_covariate_projection = true;
        NullSolver null_solver = NullSolver::Brent;
        int max_iter = 100;
        double tol = 1e-6;
    };

    EigenVce();
    EigenVce(const Options& opts);

    void prepare(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& G);
    void prepareSparse(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::SparseMatrix<double>& G_sparse);
    
    // RRM-specific VCE logic for longitudinal data
    void prepareRRM(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& Phi, const Eigen::MatrixXd& G);
    
    bool runNullModel();

    // Skip VCE and use pre-computed variance components directly
    void fixedComponents(double s2g, double s2e);

    double getSigma2g() const { return sigma2_g; }
    double getSigma2e() const { return sigma2_e; }
    double getDelta() const { return null_delta; }
    const Eigen::VectorXd& getEigenvalues() const { return D; }
    const Eigen::MatrixXd& getEigenvectors() const { return U; }
    const Eigen::MatrixXd& getProjectedW() const { return W_proj; }
    const Eigen::VectorXd& getUty() const { return Uty; }
    const Eigen::MatrixXd& getUtX() const { return UtX; }

private:
    Options options;
    int n_samples;
    int n_covars;
    
    Eigen::VectorXd D;
    Eigen::MatrixXd U;
    Eigen::MatrixXd W_proj;
    Eigen::VectorXd Uty;
    Eigen::MatrixXd UtX;

    double null_delta = 1.0;
    double sigma2_g = 0.0;
    double sigma2_e = 0.0;
    
    double computeNegLogLikelihood(double delta);
    bool runIdulNullModel(bool monotone, bool accelerated);
};

// Exact VCE for Multiple Traits using Eigen Decomposition (EM/REML)
class MvEigenVce {
public:
    enum class VceSolver {
        EM,
        AIREML
    };

    struct Options {
        int max_iter = 100;
        double tol = 1e-5;
        bool use_projection = true;
        bool verbose = true;
        VceSolver solver = VceSolver::EM;
        int em_min_iter = 10;
        int ai_max_iter = 50;
        double ai_tol = 1e-8;
        bool use_emmvgwas_lambda = false;
    };

    MvEigenVce();
    MvEigenVce(const Options& opts);

    void prepare(const Eigen::MatrixXd& Y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& G);
    bool runNullModel();
    void setEmmvgwasFixedComponents(const Eigen::MatrixXd& G_component,
                                    const Eigen::MatrixXd& R_component);

    Eigen::MatrixXd getVg() const { return Vg; }
    Eigen::MatrixXd getVe() const { return Ve; }
    Eigen::MatrixXd getVgSE() const { return Vg_se; }
    Eigen::MatrixXd getVeSE() const { return Ve_se; }
    Eigen::MatrixXd getEmmvgwasG() const { return emmvgwas_G; }
    Eigen::MatrixXd getEmmvgwasR() const { return emmvgwas_R; }
    Eigen::VectorXd getEmmvgwasBeta() const { return emmvgwas_beta; }
    const Eigen::VectorXd& getEigenvalues() const { return D; }
    const Eigen::MatrixXd& getProjectedW() const { return W_proj; }
    const Eigen::MatrixXd& getUty() const { return Uty; }
    const Eigen::MatrixXd& getUtX() const { return UtX; }

private:
    Options options;
    int n_samples;
    int n_traits;
    int n_covars;

    Eigen::VectorXd D;
    Eigen::MatrixXd U;
    Eigen::MatrixXd W_proj;
    Eigen::MatrixXd Uty;
    Eigen::MatrixXd UtX;

    Eigen::MatrixXd Vg;
    Eigen::MatrixXd Ve;
    Eigen::MatrixXd Vg_se;
    Eigen::MatrixXd Ve_se;
    Eigen::MatrixXd B;
    Eigen::MatrixXd emmvgwas_G;
    Eigen::MatrixXd emmvgwas_R;
    Eigen::VectorXd emmvgwas_beta;

    void updateEM(Eigen::MatrixXd& Vg_new, Eigen::MatrixXd& Ve_new);
    bool runAireml(int em_min_iter, int ai_max_iter);
    bool runEmmvgwasLambdaNullModel();
    bool isPositiveDefinite(const Eigen::MatrixXd& M) const;
    double calcLogLikelihood() const;
};

} // namespace cosmic
