#include "eigen_vce_ad.h"
#include "vce.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <utility>

namespace cosmic {

static Eigen::MatrixXd formVRot(const Eigen::VectorXd& D_A,
                                const Eigen::MatrixXd& K_D_rot,
                                double s2a, double s2d, double s2e) {
    Eigen::MatrixXd V_rot = s2d * K_D_rot;
    V_rot.diagonal() += s2a * D_A;
    V_rot.diagonal().array() += s2e;
    return V_rot;
}

EigenVceAD::EigenVceAD() : options(Options{}) {}
EigenVceAD::EigenVceAD(const Options& opts) : options(opts) {}

void EigenVceAD::prepare(const Eigen::VectorXd& y,
                         const Eigen::MatrixXf& X,
                         const Eigen::MatrixXd& K_A,
                         const Eigen::MatrixXd& K_D) {
    n_samples = y.size();
    n_covars = X.cols();

    if (X.rows() != n_samples || K_A.rows() != n_samples || K_A.cols() != n_samples ||
        K_D.rows() != n_samples || K_D.cols() != n_samples) {
        throw std::runtime_error("EigenVceAD::prepare: Dimension mismatch");
    }

    Eigen::VectorXd y_work;
    Eigen::MatrixXd K_A_work;
    Eigen::MatrixXd K_D_work;

    if (options.use_covariate_projection && n_covars > 0) {
        if (options.verbose)
            std::cout << "EigenVceAD: Applying Covariate Projection..." << std::endl;

        Eigen::MatrixXd X_d = X.cast<double>();
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X_d);
        Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples, n_samples);
        Eigen::MatrixXd A = Q.rightCols(n_samples - n_covars);

        A_proj = A;
        y_work = A.transpose() * y;
        K_A_work = A.transpose() * K_A * A;
        K_D_work = A.transpose() * K_D * A;
        m = n_samples - n_covars;
    } else {
        A_proj = Eigen::MatrixXd::Identity(n_samples, n_samples);
        y_work = y;
        K_A_work = K_A;
        K_D_work = K_D;
        m = n_samples;
    }

    if (options.verbose)
        std::cout << "EigenVceAD: Eigendecomposing K_A (" << m << " x " << m << ")..." << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_A(K_A_work);
    U_A = es_A.eigenvectors();
    D_A = es_A.eigenvalues();

    for (int i = 0; i < D_A.size(); ++i) {
        if (D_A[i] < 0) D_A[i] = 0;
    }

    if (options.verbose)
        std::cout << "EigenVceAD: Precomputing K_D in K_A eigenbasis..." << std::endl;

    K_D_rot = U_A.transpose() * K_D_work * U_A;
    Uty_prepared = U_A.transpose() * y_work;
    Uty = Uty_prepared;
    UtX = Eigen::MatrixXd::Zero(m, 0);
}

// Profile REML log-likelihood matching VCESolver::solveDiagonal exactly.
//
// V_rot = D_A + eta_d*K_D_rot + eta_a*I
// V = sigma_g^2 * V_rot (sigma_g^2 profiled out)
//
// logL = -0.5 * (log|V_shape| + n_eff * log(sigma_g^2) + n_eff)
static double profileLogL(
    const Eigen::VectorXd& D_A,
    const Eigen::MatrixXd& K_D_rot,
    const Eigen::VectorXd& Uty,
    double eta_a, double eta_d,
    int m, int n_covars,
    double& out_s2a, double& out_s2d, double& out_s2e) {

    double n_eff = (double)m;

    // Form V_rot = D_A + eta_d * K_D_rot + eta_a * I
    Eigen::MatrixXd V_rot = eta_d * K_D_rot;
    V_rot.diagonal() += D_A;
    V_rot.diagonal().array() += eta_a;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(V_rot);
    Eigen::VectorXd D_V = es.eigenvalues();
    Eigen::MatrixXd W = es.eigenvectors();

    for (int i = 0; i < D_V.size(); ++i) {
        if (D_V[i] < 1e-20) D_V[i] = 1e-20;
    }

    Eigen::VectorXd Uty_W = W.transpose() * Uty;
    Eigen::VectorXd V_inv_diag = D_V.cwiseInverse();

    double quad = (Uty_W.array().square() * V_inv_diag.array()).sum();
    double sigma_g2 = quad / n_eff;

    double log_det_V_shape = D_V.array().cwiseMax(1e-300).log().sum();
    double logL = -0.5 * (log_det_V_shape + n_eff * std::log(sigma_g2) + n_eff);

    out_s2a = sigma_g2;
    out_s2d = sigma_g2 * eta_d;
    out_s2e = sigma_g2 * eta_a;

    return logL;
}

// Golden section search on log10(eta) for a single coordinate
static double goldenSearch1D(
    const Eigen::VectorXd& D_A,
    const Eigen::MatrixXd& K_D_rot,
    const Eigen::VectorXd& Uty,
    double log_eta_other, bool optimize_a,
    int m, int n_covars,
    double& out_s2a, double& out_s2d, double& out_s2e) {

    auto evalAt = [&](double log_eta) -> double {
        double ea = optimize_a ? std::pow(10.0, log_eta) : std::pow(10.0, log_eta_other);
        double ed = optimize_a ? std::pow(10.0, log_eta_other) : std::pow(10.0, log_eta);
        return profileLogL(D_A, K_D_rot, Uty, ea, ed, m, n_covars, out_s2a, out_s2d, out_s2e);
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

bool EigenVceAD::runNullModel() {
    if (m <= 0) {
        throw std::runtime_error("EigenVceAD: prepare() must be called before runNullModel()");
    }

    if (options.verbose)
        std::cout << "EigenVceAD: Running Profile REML for AD model (" << m << " samples)..." << std::endl;

    // Initialize from 2-component LMM (K_A only)
    double var_g, var_e, se_g, se_e;
    VCESolver::solveDiagonal(D_A, Uty_prepared, UtX, var_g, var_e, se_g, se_e);

    double eta_a_init = var_e / (var_g + 1e-20);

    if (options.verbose)
        std::cout << "EigenVceAD: Init from LMM: var_g=" << var_g << ", var_e=" << var_e
                  << " => eta_a=" << eta_a_init << std::endl;

    // Coordinate descent over (eta_a, eta_d) using golden section search
    double log_eta_a = std::log10(std::max(eta_a_init, 1e-5));
    double log_eta_d = -5.0;

    double best_s2a = 0, best_s2d = 0, best_s2e = 0;
    double best_logL = profileLogL(D_A, K_D_rot, Uty_prepared,
                                   std::pow(10.0, log_eta_a), std::pow(10.0, log_eta_d),
                                   m, n_covars, best_s2a, best_s2d, best_s2e);

    if (options.verbose)
        std::cout << "EigenVceAD: Initial logL=" << best_logL
                  << " (eta_a=" << std::pow(10.0, log_eta_a)
                  << ", eta_d=" << std::pow(10.0, log_eta_d) << ")" << std::endl;

    // Coordinate descent iterations
    int iter = 0;
    bool reached_maxit = true;
    for (; iter < options.max_em_iter; ++iter) {
        double prev_logL = best_logL;
        double prev_eta_a = log_eta_a;
        double prev_eta_d = log_eta_d;

        // Optimize eta_a given eta_d
        log_eta_a = goldenSearch1D(D_A, K_D_rot, Uty_prepared, log_eta_d, true,
                                   m, n_covars, best_s2a, best_s2d, best_s2e);

        // Optimize eta_d given eta_a
        log_eta_d = goldenSearch1D(D_A, K_D_rot, Uty_prepared, log_eta_a, false,
                                   m, n_covars, best_s2a, best_s2d, best_s2e);

        // Recompute with final values
        best_logL = profileLogL(D_A, K_D_rot, Uty_prepared,
                                std::pow(10.0, log_eta_a), std::pow(10.0, log_eta_d),
                                m, n_covars, best_s2a, best_s2d, best_s2e);

        double rel_diff = std::abs(log_eta_a - prev_eta_a) + std::abs(log_eta_d - prev_eta_d);

        if (options.verbose && (iter + 1) % 5 == 0) {
            std::cout << "  CD iter " << (iter + 1) << ": s2a=" << best_s2a
                      << ", s2d=" << best_s2d << ", s2e=" << best_s2e
                      << ", logL=" << best_logL << ", delta_param=" << rel_diff << std::endl;
        }

        if (rel_diff < options.tol) {
            if (options.verbose)
                std::cout << "EigenVceAD: Converged in " << (iter + 1) << " iterations." << std::endl;
            reached_maxit = false;
            break;
        }
    }

    if (reached_maxit && options.verbose) {
        std::cout << "EigenVceAD: Reached max iterations (" << options.max_em_iter << "); "
                  << "tolerance not met." << std::endl;
    }

    // Termination status reflects the actual optimizer outcome, not a hard-coded success.
    // iterations_ counts completed coordinate-descent sweeps (0 when maxit=0 / no sweep ran).
    iterations_ = reached_maxit ? iter : (iter + 1);
    converged_ = !reached_maxit;
    sigma2_a = best_s2a;
    sigma2_d = best_s2d;
    sigma2_e = best_s2e;
    logL_ = best_logL;

    if (options.verbose) {
        std::cout << "EigenVceAD: Final: sigma_A^2=" << sigma2_a
                  << ", sigma_D^2=" << sigma2_d
                  << ", sigma_e^2=" << sigma2_e
                  << ", logL=" << logL_
                  << ", converged=" << (converged_ ? "true" : "false")
                  << ", iterations=" << iterations_ << std::endl;
    }

    computeFinalDecomposition();
    return true;
}

void EigenVceAD::fixedComponents(double s2a, double s2d, double s2e) {
    if (m <= 0) {
        throw std::runtime_error("EigenVceAD: prepare() must be called before fixedComponents()");
    }
    sigma2_a = s2a;
    sigma2_d = s2d;
    sigma2_e = s2e;

    Eigen::MatrixXd V_rot = formVRot(D_A, K_D_rot, s2a, s2d, s2e);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(V_rot);
    D_V = es.eigenvalues();
    Eigen::MatrixXd W = es.eigenvectors();
    for (int i = 0; i < D_V.size(); ++i) {
        if (D_V[i] < 0) D_V[i] = 0;
    }
    Eigen::VectorXd Uty_W = W.transpose() * Uty_prepared;
    Eigen::VectorXd V_inv_diag = D_V.cwiseInverse();
    double quad = (Uty_W.array().square() * V_inv_diag.array()).sum();
    double log_det_V = D_V.array().cwiseMax(1e-300).log().sum();
    logL_ = -0.5 * (log_det_V + quad + m * std::log(2.0 * M_PI));

    computeFinalDecomposition();
    converged_ = true;
    iterations_ = 0;
}

void EigenVceAD::computeFinalDecomposition() {
    Eigen::MatrixXd V_rot = formVRot(D_A, K_D_rot, sigma2_a, sigma2_d, sigma2_e);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(V_rot);
    D_V = es.eigenvalues();
    Eigen::MatrixXd W_local = es.eigenvectors();

    for (int i = 0; i < D_V.size(); ++i) {
        if (D_V[i] < 0) D_V[i] = 0;
    }

    // Compose the full operator back to original n-space: A * U_A * W_local.
    // With W_proj n x m, rotate() of an n-vector g is W_proj' * g, matching the
    // additive EigenVce convention (W_proj = A * U) and giving the P-adjusted
    // (covariate-projected) joint information matrix.
    W_proj = A_proj * U_A * W_local;

    Uty = W_local.transpose() * Uty_prepared;

    UtX = Eigen::MatrixXd::Zero(m, 0);
}

} // namespace cosmic
