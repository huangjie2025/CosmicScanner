#include "eigen_vce_d.h"
#include "vce.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>

namespace cosmic {

EigenVceD::EigenVceD() : options(Options{}) {}
EigenVceD::EigenVceD(const Options& opts) : options(opts) {}

void EigenVceD::prepare(const Eigen::VectorXd& y,
                        const Eigen::MatrixXf& X,
                        const Eigen::MatrixXd& K_D) {
    n_samples = y.size();
    n_covars = X.cols();

    if (X.rows() != n_samples || K_D.rows() != n_samples || K_D.cols() != n_samples) {
        throw std::runtime_error("EigenVceD::prepare: Dimension mismatch");
    }
    if (n_samples == 0 || n_covars >= n_samples) {
        throw std::runtime_error("EigenVceD::prepare: Covariate count must be smaller than sample count");
    }
    if (!y.allFinite() || !X.allFinite() || !K_D.allFinite()) {
        throw std::runtime_error("EigenVceD::prepare: Inputs must be finite");
    }
    if (!options.use_covariate_projection && n_covars > 0) {
        throw std::runtime_error("EigenVceD::prepare: Non-empty covariates require covariate projection");
    }

    Eigen::VectorXd y_work;
    Eigen::MatrixXd K_D_work;

    if (options.use_covariate_projection && n_covars > 0) {
        if (options.verbose)
            std::cout << "EigenVceD: Applying Covariate Projection..." << std::endl;

        Eigen::MatrixXd X_d = X.cast<double>();
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> rank_qr(X_d);
        if (rank_qr.rank() != n_covars) {
            throw std::runtime_error("EigenVceD::prepare: Covariate matrix is rank deficient");
        }
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X_d);
        Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples, n_samples);
        Eigen::MatrixXd A = Q.rightCols(n_samples - n_covars);

        A_proj = A;
        y_work = A.transpose() * y;
        K_D_work = A.transpose() * K_D * A;
        m = n_samples - n_covars;
    } else {
        A_proj = Eigen::MatrixXd::Identity(n_samples, n_samples);
        y_work = y;
        K_D_work = K_D;
        m = n_samples;
    }

    if (options.verbose)
        std::cout << "EigenVceD: Eigendecomposing K_D (" << m << " x " << m << ")..." << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_D(K_D_work);
    if (es_D.info() != Eigen::Success) {
        throw std::runtime_error("EigenVceD::prepare: K_D eigendecomposition failed");
    }
    U_D = es_D.eigenvectors();
    D_D = es_D.eigenvalues();

    for (int i = 0; i < D_D.size(); ++i) {
        if (D_D[i] < 0) D_D[i] = 0;
    }

    if (options.verbose)
        std::cout << "EigenVceD: Precomputing K_D in K_D eigenbasis..." << std::endl;

    Uty_prepared = U_D.transpose() * y_work;
    if (!(Uty_prepared.squaredNorm() > 1e-20)) {
        throw std::runtime_error("EigenVceD::prepare: Phenotype has no residual variation");
    }
    Uty = Uty_prepared;
    UtX = Eigen::MatrixXd::Zero(m, 0);
}

static double profileLogL(const Eigen::VectorXd& D_D,
                          const Eigen::VectorXd& Uty,
                          double eta_d,
                          int m,
                          double& out_s2d, double& out_s2e) {

    Eigen::VectorXd V_shape = eta_d * D_D + Eigen::VectorXd::Ones(D_D.size());

    for (int i = 0; i < V_shape.size(); ++i) {
        if (V_shape[i] < 1e-20) V_shape[i] = 1e-20;
    }

    Eigen::VectorXd V_inv_diag = V_shape.cwiseInverse();
    Eigen::VectorXd Uty_W = Uty;
    Eigen::VectorXd quad_vec = (Uty_W.array().square() * V_inv_diag.array());
    double quad = quad_vec.sum();

    double sigma_e2 = quad / static_cast<double>(m);
    double sigma_d2 = sigma_e2 * eta_d;

    double log_det_V_shape = V_shape.array().cwiseMax(1e-300).log().sum();
    double logL = -0.5 * (log_det_V_shape + m * std::log(sigma_e2) + m);

    out_s2d = sigma_d2;
    out_s2e = sigma_e2;

    return logL;
}

static double goldenSearch1D(const Eigen::VectorXd& D_D,
                             const Eigen::VectorXd& Uty,
                             int m,
                             double& out_s2d, double& out_s2e) {

    auto evalAt = [&](double log_eta) -> double {
        double eta = std::pow(10.0, log_eta);
        return profileLogL(D_D, Uty, eta, m, out_s2d, out_s2e);
    };

    double lo = -5.0, hi = 5.0;
    double gr = (std::sqrt(5.0) + 1.0) / 2.0;

    double c = hi - (hi - lo) / gr;
    double d = lo + (hi - lo) / gr;
    double fc = evalAt(c);
    double fd = evalAt(d);

    while (std::abs(c - d) > 1e-4) {
        if (fc > fd) {
            hi = d; d = c; fd = fc;
            c = hi - (hi - lo) / gr;
            fc = evalAt(c);
        } else {
            lo = c; c = d; fc = fd;
            d = lo + (hi - lo) / gr;
            fd = evalAt(d);
        }
    }

    return (lo + hi) / 2.0;
}

bool EigenVceD::runNullModel() {
    if (m <= 0) {
        throw std::runtime_error("EigenVceD: prepare() must be called before runNullModel()");
    }

    if (options.verbose)
        std::cout << "EigenVceD: Running Profile REML for D model (" << m << " samples)..." << std::endl;

    // A projected zero K_D contains no estimable dominance variance.  Treat it
    // as the residual-only boundary instead of optimizing an unidentifiable ratio.
    if (D_D.size() == 0 || D_D.maxCoeff() <= 1e-12) {
        sigma2_d = 0.0;
        sigma2_e = Uty_prepared.squaredNorm() / static_cast<double>(m);
        if (!(sigma2_e > 0.0) || !std::isfinite(sigma2_e)) {
            throw std::runtime_error("EigenVceD: Residual variance is not positive and finite");
        }
        logL_ = -0.5 * m * (std::log(sigma2_e) + 1.0);
        converged_ = true;
        iterations_ = 0;
        computeFinalDecomposition();
        return true;
    }

    double var_g, var_e, se_g, se_e;
    VCESolver::solveDiagonal(D_D, Uty_prepared, Eigen::MatrixXd::Zero(m, 0), var_g, var_e, se_g, se_e);

    double eta_d_init = var_g / (var_e + 1e-20);

    if (options.verbose)
        std::cout << "EigenVceD: Init from LMM: var_g=" << var_g << ", var_e=" << var_e
                  << " => eta_d=" << eta_d_init << std::endl;

    double log_eta_d = std::log10(std::max(eta_d_init, 1e-5));

    double best_s2d = 0, best_s2e = 0;
    double best_logL = profileLogL(D_D, Uty_prepared,
                                    std::pow(10.0, log_eta_d),
                                    m, best_s2d, best_s2e);

    if (options.verbose)
        std::cout << "EigenVceD: Initial logL=" << best_logL
                  << " (eta_d=" << std::pow(10.0, log_eta_d) << ")" << std::endl;

    int iter = 0;
    bool reached_maxit = true;
    for (; iter < options.max_em_iter; ++iter) {
        double prev_logL = best_logL;
        double prev_eta_d = log_eta_d;

        log_eta_d = goldenSearch1D(D_D, Uty_prepared, m, best_s2d, best_s2e);

        best_logL = profileLogL(D_D, Uty_prepared,
                                std::pow(10.0, log_eta_d),
                                m, best_s2d, best_s2e);

        double rel_diff = std::abs(log_eta_d - prev_eta_d);

        if (options.verbose && (iter + 1) % 5 == 0) {
            std::cout << "  CD iter " << (iter + 1) << ": s2d=" << best_s2d
                      << ", s2e=" << best_s2e
                      << ", logL=" << best_logL << ", delta_param=" << rel_diff << std::endl;
        }

        if (rel_diff < options.tol) {
            if (options.verbose)
                std::cout << "EigenVceD: Converged in " << (iter + 1) << " iterations." << std::endl;
            reached_maxit = false;
            break;
        }
    }

    if (reached_maxit && options.verbose) {
        std::cout << "EigenVceD: Reached max iterations (" << options.max_em_iter << "); "
                  << "tolerance not met." << std::endl;
    }

    iterations_ = reached_maxit ? iter : (iter + 1);
    converged_ = !reached_maxit;
    sigma2_d = best_s2d;
    sigma2_e = best_s2e;
    logL_ = best_logL;

    if (options.verbose) {
        std::cout << "EigenVceD: Final: sigma_D^2=" << sigma2_d
                  << ", sigma_e^2=" << sigma2_e
                  << ", logL=" << logL_
                  << ", converged=" << (converged_ ? "true" : "false")
                  << ", iterations=" << iterations_ << std::endl;
    }

    computeFinalDecomposition();
    return true;
}

void EigenVceD::fixedComponents(double s2d, double s2e) {
    if (m <= 0) {
        throw std::runtime_error("EigenVceD: prepare() must be called before fixedComponents()");
    }
    if (!std::isfinite(s2d) || !std::isfinite(s2e) || s2d < 0.0 || s2e <= 0.0) {
        throw std::runtime_error("EigenVceD::fixedComponents: Require sigma_D^2 >= 0 and sigma_e^2 > 0");
    }
    sigma2_d = s2d;
    sigma2_e = s2e;

    computeFinalDecomposition();
    Eigen::VectorXd V_inv_diag = D_V.cwiseMax(1e-20).cwiseInverse();
    double quad = (Uty.array().square() * V_inv_diag.array()).sum();
    double log_det_V = D_V.array().cwiseMax(1e-300).log().sum();
    logL_ = -0.5 * (log_det_V + quad + m * std::log(2.0 * M_PI));

    converged_ = true;
    iterations_ = 0;
}

void EigenVceD::computeFinalDecomposition() {
    // K_D is diagonal in U_D, so V = sigma_D^2 K_D + sigma_e^2 I
    // has the same eigenvectors.  Keep the eigenvalue order aligned with U_D.
    D_V = sigma2_d * D_D + Eigen::VectorXd::Constant(m, sigma2_e);
    D_V = D_V.cwiseMax(0.0);
    W_proj = A_proj * U_D;
    Uty = Uty_prepared;
}

} // namespace cosmic
