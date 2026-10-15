#include "eigen_vce.h"
#include "vce.h" // For VCESolver::solveDiagonal
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

#include <unsupported/Eigen/KroneckerProduct>

namespace cosmic {

// --- EigenVce (Single Trait) ---

EigenVce::EigenVce() : options(Options{}) {}
EigenVce::EigenVce(const Options& opts) : options(opts) {}

void EigenVce::prepare(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& G) {
    n_samples = y.size();
    n_covars = X.cols();

    if (X.rows() != n_samples || G.rows() != n_samples || G.cols() != n_samples) {
        throw std::runtime_error("EigenVce::prepare: Dimension mismatch");
    }

    if (options.use_covariate_projection && n_covars > 0) {
        std::cout << "EigenVce: Applying Covariate Projection to y and G..." << std::endl;
        
        Eigen::MatrixXd X_d = X.cast<double>();
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X_d);
        Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples, n_samples);
        
        // A is the null space of X (n x (n-c))
        Eigen::MatrixXd A = Q.rightCols(n_samples - n_covars);
        
        // Project y and G
        Eigen::VectorXd y_reml = A.transpose() * y;
        Eigen::MatrixXd G_reml = A.transpose() * G * A;
        
        std::cout << "EigenVce: Eigen Decomposition of projected G (" << (n_samples - n_covars) << " x " << (n_samples - n_covars) << ")..." << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G_reml);
        D = es.eigenvalues();
        
        for(int i=0; i<D.size(); ++i) if(D[i] < 0) D[i] = 0;
        
        U = es.eigenvectors();
        Uty = U.transpose() * y_reml;
        W_proj = A * U;
        UtX = Eigen::MatrixXd::Zero(n_samples - n_covars, 0);
    } else {
        std::cout << "EigenVce: Eigen Decomposition of G (" << n_samples << " x " << n_samples << ")..." << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
        D = es.eigenvalues();
        
        for(int i=0; i<D.size(); ++i) if(D[i] < 0) D[i] = 0;
        
        U = es.eigenvectors();
        Uty = U.transpose() * y;
        UtX = U.transpose() * X.cast<double>();
        W_proj = U; // Standard rotation
    }
}

void EigenVce::prepareSparse(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::SparseMatrix<double>& G_sparse) {
    std::cout << "EigenVce: Converting sparse GRM to dense for Eigen Decomposition..." << std::endl;
    Eigen::MatrixXd G_dense = Eigen::MatrixXd(G_sparse);
    prepare(y, X, G_dense);
}

void EigenVce::prepareRRM(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& Phi, const Eigen::MatrixXd& G) {
    std::cout << "EigenVce: Preparing RRM specific VCE logic..." << std::endl;
    // For RRM, the random effect covariance is V = Z * (K \otimes G) * Z' + e
    // If K is scalar (order 0), V = Z G Z' * sigma2_g + I * sigma2_e
    // Let's implement the simplified order 0 case first for exact decomposition
    
    Eigen::MatrixXd ZGZt = Z * G * Z.transpose();
    
    // Fall back to standard dense preparation
    prepare(y, X, ZGZt);
}

bool EigenVce::runNullModel() {
    if (options.null_solver == NullSolver::EM || options.null_solver == NullSolver::AI) {
        return runIterativeReml(options.null_solver == NullSolver::AI);
    }
    if (options.null_solver == NullSolver::IDUL ||
        options.null_solver == NullSolver::IDULPlus ||
        options.null_solver == NullSolver::AIDUL) {
        const bool monotone = (options.null_solver == NullSolver::IDULPlus || options.null_solver == NullSolver::AIDUL);
        const bool accelerated = (options.null_solver == NullSolver::AIDUL);
        return runIdulNullModel(monotone, accelerated);
    }

    std::cout << "EigenVce: Running Null Model (Brent Search)..." << std::endl;
    double var_g, var_e, se_g, se_e;
    Eigen::VectorXd beta, beta_se;
    
    // We reuse the existing VCESolver::solveDiagonal which implements the Brent search
    VCESolver::solveDiagonal(D, Uty, UtX, var_g, var_e, se_g, se_e, &beta, &beta_se);
    bool success = true;
    
    sigma2_g = var_g;
    sigma2_e = var_e;
    sigma2_g_se_ = se_g;
    sigma2_e_se_ = se_e;
    null_delta = sigma2_e / sigma2_g;
    converged_ = success;
    iterations_ = 1;
    updateRemlDiagnostics();
    
    return success;
}

void EigenVce::updateRemlDiagnostics() {
    const int n = static_cast<int>(D.size());
    const int p = static_cast<int>(UtX.cols());
    const int df = n - p;
    const Eigen::VectorXd covariance =
        (sigma2_g * D.array() + sigma2_e).matrix().cwiseMax(1e-12);
    const Eigen::VectorXd inverse = covariance.cwiseInverse();
    Eigen::MatrixXd P = inverse.asDiagonal();
    double log_det_fixed = 0.0;
    if (p > 0) {
        const Eigen::MatrixXd information = UtX.transpose() * inverse.asDiagonal() * UtX;
        Eigen::LDLT<Eigen::MatrixXd> solver(information);
        if (solver.info() != Eigen::Success) return;
        P.noalias() -= inverse.asDiagonal() * UtX *
                       solver.solve(Eigen::MatrixXd::Identity(p, p)) *
                       UtX.transpose() * inverse.asDiagonal();
        log_det_fixed = solver.vectorD().array().abs().max(1e-300).log().sum();
    }
    const Eigen::VectorXd py = P * Uty;
    const Eigen::VectorXd kpy_g = D.array() * py.array();
    Eigen::Vector2d score;
    score[0] = 0.5 * (py.dot(kpy_g) - (P.diagonal().array() * D.array()).sum());
    score[1] = 0.5 * (py.squaredNorm() - P.trace());
    Eigen::Matrix2d ai;
    ai(0, 0) = 0.5 * kpy_g.dot(P * kpy_g);
    ai(0, 1) = ai(1, 0) = 0.5 * kpy_g.dot(P * py);
    ai(1, 1) = 0.5 * py.dot(P * py);
    Eigen::LDLT<Eigen::Matrix2d> ai_solver(ai);
    if (ai_solver.info() == Eigen::Success) {
        const Eigen::Matrix2d parameter_covariance =
            ai_solver.solve(Eigen::Matrix2d::Identity());
        sigma2_g_se_ = std::sqrt(std::max(0.0, parameter_covariance(0, 0)));
        sigma2_e_se_ = std::sqrt(std::max(0.0, parameter_covariance(1, 1)));
    }
    max_score_ = score.cwiseAbs().maxCoeff();
    log_likelihood_ = -0.5 *
        (covariance.array().log().sum() + log_det_fixed + Uty.dot(py) +
         df * std::log(2.0 * 3.14159265358979323846));
}

bool EigenVce::runIterativeReml(bool average_information) {
    const int n = static_cast<int>(D.size());
    const int p = static_cast<int>(UtX.cols());
    const int df = n - p;
    if (df < 2) throw std::runtime_error("EigenVce REML requires residual degrees of freedom");
    const double phenotype_variance = std::max(1e-10, Uty.squaredNorm() / df);
    const double lower = phenotype_variance * 1e-10;
    Eigen::Vector2d theta(0.5 * phenotype_variance, 0.5 * phenotype_variance);

    struct Terms {
        double log_likelihood = -std::numeric_limits<double>::infinity();
        Eigen::Vector2d score = Eigen::Vector2d::Zero();
        Eigen::Matrix2d ai = Eigen::Matrix2d::Zero();
        double trace_g = 0.0;
        double trace_e = 0.0;
        double quad_g = 0.0;
        double quad_e = 0.0;
    };
    auto evaluate = [&](const Eigen::Vector2d& values) {
        Terms terms;
        const Eigen::VectorXd covariance =
            (values[0] * D.array() + values[1]).matrix().cwiseMax(lower);
        const Eigen::VectorXd inverse = covariance.cwiseInverse();
        Eigen::MatrixXd P = inverse.asDiagonal();
        double log_det_fixed = 0.0;
        if (p > 0) {
            const Eigen::MatrixXd information = UtX.transpose() * inverse.asDiagonal() * UtX;
            Eigen::LDLT<Eigen::MatrixXd> solver(information);
            if (solver.info() != Eigen::Success)
                throw std::runtime_error("EigenVce REML fixed-effect information matrix is singular");
            const Eigen::MatrixXd information_inverse = solver.solve(
                Eigen::MatrixXd::Identity(p, p));
            P.noalias() -= inverse.asDiagonal() * UtX * information_inverse *
                           UtX.transpose() * inverse.asDiagonal();
            log_det_fixed = solver.vectorD().array().abs().max(1e-300).log().sum();
        }
        const Eigen::VectorXd py = P * Uty;
        const Eigen::VectorXd kpy_g = D.array() * py.array();
        terms.trace_g = (P.diagonal().array() * D.array()).sum();
        terms.trace_e = P.trace();
        terms.quad_g = py.dot(kpy_g);
        terms.quad_e = py.squaredNorm();
        terms.score[0] = 0.5 * (terms.quad_g - terms.trace_g);
        terms.score[1] = 0.5 * (terms.quad_e - terms.trace_e);
        terms.ai(0, 0) = 0.5 * kpy_g.dot(P * kpy_g);
        terms.ai(0, 1) = terms.ai(1, 0) = 0.5 * kpy_g.dot(P * py);
        terms.ai(1, 1) = 0.5 * py.dot(P * py);
        terms.log_likelihood = -0.5 *
            (covariance.array().log().sum() + log_det_fixed + Uty.dot(py) +
             df * std::log(2.0 * 3.14159265358979323846));
        return terms;
    };

    Terms current = evaluate(theta);
    converged_ = false;
    iterations_ = 0;
    for (int iteration = 0; iteration < options.max_iter; ++iteration) {
        Eigen::Vector2d proposal = theta;
        if (average_information) {
            Eigen::LDLT<Eigen::Matrix2d> solver(current.ai);
            if (solver.info() != Eigen::Success) break;
            const Eigen::Vector2d step = solver.solve(current.score);
            double fraction = 1.0;
            for (int halving = 0; halving < 40; ++halving) {
                proposal = (theta + fraction * step).cwiseMax(lower);
                if (evaluate(proposal).log_likelihood >= current.log_likelihood - 1e-10) break;
                fraction *= 0.5;
            }
        } else {
            proposal[0] = theta[0] * current.quad_g / std::max(1e-12, current.trace_g);
            proposal[1] = theta[1] * current.quad_e / std::max(1e-12, current.trace_e);
            proposal = (0.5 * theta + 0.5 * proposal).cwiseMax(lower);
        }
        const double change = ((proposal - theta).array().abs() /
            (theta.array().abs() + phenotype_variance * 1e-8)).maxCoeff();
        theta = proposal;
        current = evaluate(theta);
        iterations_ = iteration + 1;
        if (change < options.tol) {
            converged_ = true;
            break;
        }
    }
    sigma2_g = std::max(0.0, theta[0] <= 1.01 * lower ? 0.0 : theta[0]);
    sigma2_e = std::max(0.0, theta[1] <= 1.01 * lower ? 0.0 : theta[1]);
    null_delta = sigma2_g > 0.0 ? sigma2_e / sigma2_g
                                : std::numeric_limits<double>::infinity();
    log_likelihood_ = current.log_likelihood;
    max_score_ = current.score.cwiseAbs().maxCoeff();
    Eigen::LDLT<Eigen::Matrix2d> ai_solver(current.ai);
    if (ai_solver.info() == Eigen::Success) {
        const Eigen::Matrix2d covariance = ai_solver.solve(Eigen::Matrix2d::Identity());
        sigma2_g_se_ = std::sqrt(std::max(0.0, covariance(0, 0)));
        sigma2_e_se_ = std::sqrt(std::max(0.0, covariance(1, 1)));
    }
    return converged_ && std::isfinite(sigma2_g) && std::isfinite(sigma2_e) &&
           std::isfinite(log_likelihood_) && std::isfinite(max_score_);
}

void EigenVce::fixedComponents(double s2g, double s2e) {
    if (n_samples <= 0) throw std::runtime_error("EigenVce: prepare() must be called before fixedComponents()");
    sigma2_g = s2g;
    sigma2_e = s2e;
    null_delta = sigma2_e / (sigma2_g + 1e-10);
    converged_ = true;
    iterations_ = 0;
}

bool EigenVce::runIdulNullModel(bool monotone, bool accelerated) {
    const char* name = accelerated ? "AIDUL" : (monotone ? "IDUL_PLUS" : "IDUL");
    std::cout << "EigenVce: Running Null Model (" << name << ")..." << std::endl;

    const int n = (int)D.size();
    const int p = (int)UtX.cols();
    const double n_eff = std::max(1.0, (double)(n - p));
    const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
    const double eta_min = 1e-5;
    const double eta_max = 1e5;

    auto clamp_eta = [&](double eta) {
        if (!std::isfinite(eta)) return eta_min;
        return std::min(eta_max, std::max(eta_min, eta));
    };

    auto fit_residual = [&](double eta,
                            Eigen::VectorXd& H,
                            Eigen::VectorXd& Hinv,
                            Eigen::VectorXd& resid,
                            Eigen::MatrixXd* XtHinvX_out,
                            Eigen::MatrixXd* XtHinvX_inv_out) -> double {
        H = eta * D.array() + 1.0;
        Hinv = H.cwiseInverse();
        Eigen::VectorXd beta;
        if (p > 0) {
            Eigen::MatrixXd Xw = UtX.array().colwise() * Hinv.array().sqrt();
            Eigen::VectorXd yw = Uty.array() * Hinv.array().sqrt();
            Eigen::MatrixXd XtHinvX = Xw.transpose() * Xw;
            XtHinvX.diagonal().array() += 1e-10;
            Eigen::LDLT<Eigen::MatrixXd> solver(XtHinvX);
            if (solver.info() == Eigen::Success) {
                beta = solver.solve(Xw.transpose() * yw);
                if (XtHinvX_inv_out) {
                    *XtHinvX_inv_out = solver.solve(Eigen::MatrixXd::Identity(p, p));
                }
            } else {
                beta = XtHinvX.completeOrthogonalDecomposition().solve(Xw.transpose() * yw);
                if (XtHinvX_inv_out) {
                    *XtHinvX_inv_out = XtHinvX.completeOrthogonalDecomposition().pseudoInverse();
                }
            }
            if (XtHinvX_out) *XtHinvX_out = XtHinvX;
            resid = Uty - UtX * beta;
        } else {
            resid = Uty;
            if (XtHinvX_out) XtHinvX_out->resize(0, 0);
            if (XtHinvX_inv_out) XtHinvX_inv_out->resize(0, 0);
        }
        return (resid.array().square() * Hinv.array()).sum();
    };

    auto reml_loglike = [&](double eta) -> double {
        Eigen::VectorXd H, Hinv, resid;
        Eigen::MatrixXd XtHinvX, unused_inv;
        double quad = fit_residual(eta, H, Hinv, resid, &XtHinvX, &unused_inv);
        if (!(quad > 0.0) || !std::isfinite(quad)) return -std::numeric_limits<double>::infinity();
        double log_det_H = H.array().log().sum();
        double log_det_X = 0.0;
        if (p > 0) {
            Eigen::LDLT<Eigen::MatrixXd> solver(XtHinvX);
            if (solver.info() == Eigen::Success) {
                log_det_X = solver.vectorD().array().abs().log().sum();
            }
        }
        return -0.5 * (log_det_H + log_det_X + n_eff * std::log(quad / n_eff) + n_eff);
    };

    // Fisher-scoring update: g(eta) = eta + fp * eta / n / vh
    auto fisher_update = [&](double eta, Eigen::VectorXd& H, Eigen::VectorXd& Hinv,
                             Eigen::VectorXd& r2, Eigen::MatrixXd& XtHinvX_inv) -> double {
        const double tr_h1 = Hinv.sum();
        const double tr_h2 = Hinv.array().square().sum();
        const double mh = tr_h1 / (double)n;
        const double vh = tr_h2 / (double)n - mh * mh;
        if (vh <= 1e-14) return eta;
        double tr_projected = 0.0;
        if (p > 0) {
            Eigen::MatrixXd Xh2 = UtX.array().colwise() * Hinv.array().square();
            Eigen::MatrixXd XtH2X = UtX.transpose() * Xh2;
            tr_projected = (XtHinvX_inv * XtH2X).trace();
        }
        double rh1r = (r2.array() * Hinv.array()).sum();
        double rh2r = (r2.array() * Hinv.array().square()).sum();
        double fp = tr_h1 - tr_projected - n_eff * rh2r / std::max(rh1r, 1e-20);
        double lazy = 0.5 + 0.5 * mh;
        return eta + fp * eta / (double)n / vh * lazy;
    };

    // Anderson acceleration for scalar fixed-point iteration
    // x = g(x), f = g(x) - x; finds optimal linear combination of recent iterates
    struct AndersonAccel {
        int window;
        std::vector<double> x_hist;
        std::vector<double> gx_hist;

        double apply(double eta, double g_eta) {
            x_hist.push_back(eta);
            gx_hist.push_back(g_eta);
            while ((int)x_hist.size() > window) {
                x_hist.erase(x_hist.begin());
                gx_hist.erase(gx_hist.begin());
            }
            int m = (int)x_hist.size();
            if (m < 2) return g_eta;

            std::vector<double> f(m);
            for (int i = 0; i < m; ++i) f[i] = gx_hist[i] - x_hist[i];

            if (m == 2) {
                double denom = f[0] - f[1];
                if (std::abs(denom) > 1e-15) {
                    double alpha0 = f[1] / denom;
                    return alpha0 * gx_hist[0] + (1.0 - alpha0) * gx_hist[1];
                }
                return g_eta;
            }

            Eigen::MatrixXd A(m-1, m-1);
            Eigen::VectorXd b(m-1);
            double f_m = f[m-1];
            for (int i = 0; i < m-1; ++i) {
                for (int j = 0; j < m-1; ++j) {
                    A(i, j) = (f[i] - f_m) * (f[j] - f_m);
                }
                A(i, i) += 1e-10;
                b(i) = -f_m * (f[i] - f_m);
            }
            Eigen::VectorXd beta = A.ldlt().solve(b);
            double result = 0.0;
            double sum_beta = 0.0;
            for (int i = 0; i < m-1; ++i) {
                result += beta(i) * gx_hist[i];
                sum_beta += beta(i);
            }
            result += (1.0 - sum_beta) * gx_hist[m-1];
            return result;
        }
    };

    AndersonAccel anderson;
    anderson.window = accelerated ? 5 : 1;

    double eta = 0.1;
    double like = reml_loglike(eta);
    bool converged = false;
    int iter = 0;
    for (; iter < options.max_iter; ++iter) {
        Eigen::VectorXd H, Hinv, resid;
        Eigen::MatrixXd XtHinvX, XtHinvX_inv;
        fit_residual(eta, H, Hinv, resid, &XtHinvX, &XtHinvX_inv);
        Eigen::VectorXd r2 = resid.array().square();
        double tau_inv = (r2.array() * Hinv.array()).mean();
        if (!(tau_inv > 0.0) || !std::isfinite(tau_inv)) break;

        double eta_new = eta;
        if (accelerated) {
            double g_eta = fisher_update(eta, H, Hinv, r2, XtHinvX_inv);
            g_eta = clamp_eta(g_eta);
            eta_new = anderson.apply(eta, g_eta);
        } else {
            Eigen::MatrixXd Z(n, 2);
            Z.col(0) = ones;
            Z.col(1) = D;
            Eigen::VectorXd w2 = Hinv.array().square();
            Eigen::MatrixXd Zw = Z.array().colwise() * w2.array().sqrt();
            Eigen::VectorXd r2w = r2.array() * w2.array().sqrt();
            Eigen::Vector2d coef = (Zw.transpose() * Zw).ldlt().solve(Zw.transpose() * r2w);
            double mu = coef(0);
            double gamma = coef(1);
            if (monotone) {
                eta_new = gamma / tau_inv + (1.0 - mu / tau_inv) * eta;
            } else {
                double eta_mu = gamma / std::max(mu, 1e-20);
                double eta_tau = gamma / tau_inv;
                eta_new = 0.5 * (eta_mu + eta_tau);
            }
        }
        eta_new = clamp_eta(eta_new);

        if (monotone) {
            double proposal_like = reml_loglike(eta_new);
            int halves = 0;
            while (proposal_like < like && std::abs(eta_new / (1.0 + eta_new) - eta / (1.0 + eta)) > options.tol && halves < 60) {
                eta_new = 0.5 * (eta_new + eta);
                proposal_like = reml_loglike(eta_new);
                ++halves;
            }
            if (proposal_like >= like) like = proposal_like;
        }

        double diff = std::abs(eta_new / (1.0 + eta_new) - eta / (1.0 + eta));
        eta = eta_new;
        if (diff < options.tol) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        std::cout << "  " << name << " did not converge in " << options.max_iter
                  << " iterations (eta=" << eta << "), falling back to Brent." << std::endl;
        double var_g, var_e, se_g, se_e;
        Eigen::VectorXd beta, beta_se;
        VCESolver::solveDiagonal(D, Uty, UtX, var_g, var_e, se_g, se_e, &beta, &beta_se);
        sigma2_g = var_g;
        sigma2_e = var_e;
        null_delta = sigma2_e / sigma2_g;
        converged_ = true;
        iterations_ = options.max_iter;
        updateRemlDiagnostics();
        return true;
    }

    Eigen::VectorXd H, Hinv, resid;
    Eigen::MatrixXd XtHinvX, XtHinvX_inv;
    double quad = fit_residual(eta, H, Hinv, resid, &XtHinvX, &XtHinvX_inv);
    sigma2_e = quad / n_eff;
    sigma2_g = eta * sigma2_e;
    null_delta = 1.0 / eta;
    std::cout << "  " << name << " eta=" << eta << ", iterations=" << (iter + 1) << std::endl;
    const bool ok = std::isfinite(sigma2_g) && std::isfinite(sigma2_e) && sigma2_g >= 0.0 && sigma2_e >= 0.0;
    converged_ = ok;
    iterations_ = iter + 1;
    updateRemlDiagnostics();
    return ok;
}

// --- MvEigenVce (Multiple Traits) ---

MvEigenVce::MvEigenVce() {}
MvEigenVce::MvEigenVce(const Options& opts) : options(opts) {}

void MvEigenVce::prepare(const Eigen::MatrixXd& Y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& G) {
    n_samples = Y.rows();
    n_traits = Y.cols();
    n_covars = X.cols();

    if (X.rows() != n_samples || G.rows() != n_samples || G.cols() != n_samples) {
        throw std::runtime_error("MvEigenVce::prepare: Dimension mismatch");
    }

    if (options.use_projection && n_covars > 0) {
        if (options.verbose) std::cout << "MvEigenVce: Applying Covariate Projection to Y and G..." << std::endl;
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
        Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples, n_samples);
        
        Eigen::MatrixXd A = Q.rightCols(n_samples - n_covars);
        
        Eigen::MatrixXd Y_reml = A.transpose() * Y;
        Eigen::MatrixXd G_reml = A.transpose() * G * A;
        
        if (options.verbose) std::cout << "MvEigenVce: Eigen Decomposition of projected G (" << (n_samples - n_covars) << " x " << (n_samples - n_covars) << ")..." << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G_reml);
        D = es.eigenvalues();
        for (int i=0; i<D.size(); ++i) if (D[i] < 0) D[i] = 0;
        
        U = es.eigenvectors(); 
        W_proj = A * U; 
        
        Uty = U.transpose() * Y_reml; 
        UtX = Eigen::MatrixXd::Zero(n_samples - n_covars, 0); 
        
        n_samples = n_samples - n_covars;
        n_covars = 0;
    } else {
        if (options.verbose) std::cout << "MvEigenVce: Eigen Decomposition of G (" << G.rows() << " x " << G.cols() << ")..." << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
        D = es.eigenvalues();
        for (int i=0; i<D.size(); ++i) {
            if (D[i] < 0) D[i] = 0;
        }
        
        if (options.verbose) std::cout << "MvEigenVce: Rotating Y and X..." << std::endl;
        U = es.eigenvectors(); 
        Uty = U.transpose() * Y;
        UtX = U.transpose() * X;
        W_proj = U; 
    }
    
    // Initialize Parameters
    Eigen::MatrixXd Y_centered = Y.rowwise() - Y.colwise().mean();
    Eigen::MatrixXd CovY = (Y_centered.transpose() * Y_centered) / (double)(n_samples - 1);
    Vg = CovY * 0.5;
    Ve = CovY * 0.5;
    B = Eigen::MatrixXd::Zero(n_covars, n_traits);
}

void MvEigenVce::updateEM(Eigen::MatrixXd& Vg_new, Eigen::MatrixXd& Ve_new) {
    Eigen::MatrixXd S_g = Eigen::MatrixXd::Zero(n_traits, n_traits);
    Eigen::MatrixXd S_e = Eigen::MatrixXd::Zero(n_traits, n_traits);
    
    int pd = n_covars * n_traits;
    Eigen::MatrixXd LHS_B = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd RHS_B = Eigen::VectorXd::Zero(pd);
    
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd V_i = Vg * D[i] + Ve;
        Eigen::MatrixXd V_i_inv = V_i.inverse();
        
        Eigen::VectorXd y_i = Uty.row(i).transpose();
        Eigen::VectorXd x_i;
        if (n_covars > 0) x_i = UtX.row(i).transpose();
        
        Eigen::VectorXd resid_i = y_i;
        if (n_covars > 0) {
            resid_i -= B.transpose() * x_i;
            
            Eigen::MatrixXd X_i_kron = Eigen::kroneckerProduct(x_i.transpose(), Eigen::MatrixXd::Identity(n_traits, n_traits));
            LHS_B += X_i_kron.transpose() * V_i_inv * X_i_kron;
            RHS_B += X_i_kron.transpose() * V_i_inv * y_i;
        }
        
        Eigen::MatrixXd P_i = V_i_inv; 
        
        Eigen::MatrixXd W_g_i = Vg * D[i] * P_i;
        Eigen::VectorXd u_hat_i = W_g_i * resid_i;

        if (D[i] > 1e-12) {
            S_g += (u_hat_i * u_hat_i.transpose()) / D[i] + Vg - W_g_i * Vg;
        } else {
            S_g += Vg;
        }
        
        Eigen::MatrixXd W_e_i = Ve * P_i;
        Eigen::VectorXd e_hat_i = W_e_i * resid_i;
        
        S_e += e_hat_i * e_hat_i.transpose() + Ve - W_e_i * Ve;
    }
    
    Vg_new = S_g / n_samples;
    Ve_new = S_e / n_samples;
    
    if (n_covars > 0) {
        Eigen::VectorXd b_vec = LHS_B.ldlt().solve(RHS_B);
        for (int c = 0; c < n_covars; ++c) {
            for (int t = 0; t < n_traits; ++t) {
                B(c, t) = b_vec(c * n_traits + t);
            }
        }
    }
}

bool MvEigenVce::isPositiveDefinite(const Eigen::MatrixXd& M) const {
    if (M.rows() != M.cols()) return false;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
    if (ldlt.info() != Eigen::Success) return false;
    return (ldlt.vectorD().array() > 0.0).all();
}

double MvEigenVce::calcLogLikelihood() const {
    // REML log-likelihood after covariate projection: no X remaining
    // l = -0.5 * Σ_i [log|V_i| + y_i' V_i⁻¹ y_i] + const
    double logl = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        Eigen::MatrixXd V_i = D(i) * Vg + Ve;
        V_i.diagonal().array() += 1e-8;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(V_i);
        if (ldlt.info() != Eigen::Success) return -std::numeric_limits<double>::infinity();
        Eigen::VectorXd yi = Uty.row(i).transpose();
        Eigen::VectorXd sol = ldlt.solve(yi);
        double quad = yi.dot(sol);
        double ldet = ldlt.vectorD().array().abs().log().sum();
        logl += -0.5 * (ldet + quad);
    }
    return logl;
}

bool MvEigenVce::runAireml(int em_min_iter, int ai_max_iter) {
    if (options.verbose)
        std::cout << "MvEigenVce: Running Multivariate Null Model (NR/REML)..." << std::endl;

    const int v_size = n_traits * (n_traits + 1) / 2;
    const int n_params = 2 * v_size;

    // Precompute (j,k) pair for each packed index
    std::vector<std::pair<int,int>> pairs(v_size);
    for (int idx = 0, j = 0; j < n_traits; ++j)
        for (int k = j; k < n_traits; ++k, ++idx)
            pairs[idx] = {j, k};

    // --- Phase 1: EM initialization ---
    if (options.verbose)
        std::cout << "  Phase 1: EM initialization (" << em_min_iter << " min iters)..." << std::endl;
    int em_iter = 0;
    for (; em_iter < options.max_iter; ++em_iter) {
        Eigen::MatrixXd Vg_new, Ve_new;
        updateEM(Vg_new, Ve_new);
        double err = ((Vg_new - Vg).squaredNorm() + (Ve_new - Ve).squaredNorm()) /
                     std::max(1.0, 2.0 * n_traits * n_traits);
        Vg = Vg_new;
        Ve = Ve_new;
        if (err < options.tol && em_iter >= em_min_iter) break;
    }
    if (options.verbose)
        std::cout << "  EM done after " << (em_iter + 1) << " iterations." << std::endl;

    // Ensure positive definite starting point via eigenvalue clipping
    auto nearPD = [](Eigen::MatrixXd& M) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
        Eigen::VectorXd ev = es.eigenvalues();
        if (ev.minCoeff() <= 1e-10) {
            ev = ev.cwiseMax(1e-10);
            M = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
        }
    };
    nearPD(Vg);
    nearPD(Ve);

    // --- Phase 2: AI-REML iterations ---
    if (options.verbose)
        std::cout << "  Phase 2: NR/REML refinement..." << std::endl;

    double logl_old = calcLogLikelihood();
    const double logl_initial = logl_old;
    Eigen::MatrixXd Vg_initial = Vg, Ve_initial = Ve;
    if (!std::isfinite(logl_old)) {
        if (options.verbose)
            std::cout << "  Warning: Non-finite log-likelihood at EM exit, falling back to EM." << std::endl;
        return true; // EM result is good enough
    }

    const double reg = 1e-8;
    for (int iter = 0; iter < ai_max_iter; ++iter) {
        // Build score (gradient) and AI matrix
        Eigen::VectorXd score = Eigen::VectorXd::Zero(n_params);
        Eigen::MatrixXd AI = Eigen::MatrixXd::Zero(n_params, n_params);

        for (int i = 0; i < n_samples; ++i) {
            const double d_i = D(i);
            const double d2 = d_i * d_i;

            Eigen::MatrixXd V_i = d_i * Vg + Ve;
            V_i.diagonal().array() += 1e-8;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(V_i);
            if (ldlt.info() != Eigen::Success) return false;

            Eigen::MatrixXd V_inv = ldlt.solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
            Eigen::VectorXd yi = Uty.row(i).transpose();
            Eigen::VectorXd v = V_inv * yi;

            // Precompute w_pq = V_inv * E_pq * v for all (p,q)
            std::vector<Eigen::VectorXd> w(v_size);
            for (int a = 0; a < v_size; ++a) {
                int p = pairs[a].first, q = pairs[a].second;
                if (p == q) {
                    w[a] = v(p) * V_inv.col(p);
                } else {
                    w[a] = v(q) * V_inv.col(p) + v(p) * V_inv.col(q);
                }
            }

            // Score contributions
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                double vjvk = v(j) * v(k);
                double tr  = (j == k) ? V_inv(j, j) : 2.0 * V_inv(j, k);
                double quad = (j == k) ? vjvk : 2.0 * vjvk;
                score(a)          += 0.5 * d_i * (quad - tr);  // Vg
                score(a + v_size) += 0.5 *       (quad - tr);  // Ve
            }

            // AI matrix contributions
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                for (int b = a; b < v_size; ++b) {
                    int p = pairs[b].first, q = pairs[b].second;
                    double qval;
                    if (j == k) {
                        qval = v(j) * w[b](j);
                    } else {
                        qval = v(j) * w[b](k) + v(k) * w[b](j);
                    }
                    qval *= 0.5;

                    // Vg-Vg block
                    AI(a, b)                   += d2 * qval;
                    // Vg-Ve block
                    AI(a, v_size + b)          += d_i * qval;
                    // Ve-Vg block
                    AI(v_size + a, b)          += d_i * qval;
                    // Ve-Ve block
                    AI(v_size + a, v_size + b) +=       qval;
                }
            }
        }

        // Symmetrize AI (some Vg-Ve cross-block entries only accumulated in one triangle)
        AI = 0.5 * (AI + AI.transpose()).eval();

        // Regularize AI and solve: delta = AI⁻¹ * g
        AI.diagonal().array() *= (1.0 + reg);
        Eigen::LDLT<Eigen::MatrixXd> hess_solver(AI);
        if (hess_solver.info() != Eigen::Success) {
            if (options.verbose)
                std::cout << "  Hessian LDLT failed, breaking." << std::endl;
            break;
        }
        Eigen::VectorXd delta = hess_solver.solve(score);

        // AI-REML step-halving: θ_new = θ_old + step * AI⁻¹ * g
        double step = 1.0;
        bool accepted = false;
        Eigen::MatrixXd Vg_test, Ve_test;
        for (int h = 0; h < 12; ++h) {
            Vg_test = Vg;
            Ve_test = Ve;
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                Vg_test(j, k) += step * delta(a);
                Vg_test(k, j) = Vg_test(j, k);
            }
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                Ve_test(j, k) += step * delta(a + v_size);
                Ve_test(k, j) = Ve_test(j, k);
            }
            // Project to PD using eigenvalue clipping
            {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_vg(Vg_test);
                Eigen::VectorXd ev_vg = es_vg.eigenvalues();
                if (ev_vg.minCoeff() <= 1e-10) {
                    ev_vg = ev_vg.cwiseMax(1e-10);
                    Vg_test = es_vg.eigenvectors() * ev_vg.asDiagonal() * es_vg.eigenvectors().transpose();
                }
            }
            {
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_ve(Ve_test);
                Eigen::VectorXd ev_ve = es_ve.eigenvalues();
                if (ev_ve.minCoeff() <= 1e-10) {
                    ev_ve = ev_ve.cwiseMax(1e-10);
                    Ve_test = es_ve.eigenvectors() * ev_ve.asDiagonal() * es_ve.eigenvectors().transpose();
                }
            }

            Eigen::MatrixXd Vg_save = Vg, Ve_save = Ve;
            Vg = Vg_test;
            Ve = Ve_test;
            double logl_new = calcLogLikelihood();
            Vg = Vg_save;
            Ve = Ve_save;

            if (std::isfinite(logl_new) && logl_new >= logl_old - 1e-10) {
                accepted = true;
                break;
            }
            step *= 0.5;
        }

        if (!accepted) {
            if (options.verbose)
                std::cout << "  NR iter " << (iter + 1) << " rejected, stopping." << std::endl;
            break;
        }

        Vg = Vg_test;
        Ve = Ve_test;
        double logl_new = calcLogLikelihood();
        double err = delta.norm() / std::sqrt((double)n_params);
        if (options.verbose)
            std::cout << "  NR iter " << (iter + 1) << ": logl=" << logl_new
                      << " step=" << step << " |Δ|=" << err << std::endl;
        logl_old = logl_new;

        if (err < options.ai_tol) {
            if (options.verbose)
                std::cout << "  NR converged in " << (iter + 1) << " iterations." << std::endl;
            return true;
        }
    }

    // Fall back to EM unless AI-REML converged with clearly better logl
    if (logl_old < logl_initial + 1e-8) {
        Vg = Vg_initial;
        Ve = Ve_initial;
        if (options.verbose && logl_old < logl_initial - 1e-10)
            std::cout << "  Reverting to EM params (logl=" << logl_initial << ")." << std::endl;
    }

    // Compute SE from last AI matrix inverse
    Vg_se = Eigen::MatrixXd::Zero(n_traits, n_traits);
    Ve_se = Eigen::MatrixXd::Zero(n_traits, n_traits);
    {
        // Rebuild AI at final parameters for accurate SE
        Eigen::MatrixXd AI_final = Eigen::MatrixXd::Zero(n_params, n_params);
        for (int i = 0; i < n_samples; ++i) {
            const double d_i = D(i);
            const double d2 = d_i * d_i;
            Eigen::MatrixXd V_i = d_i * Vg + Ve;
            V_i.diagonal().array() += 1e-8;
            Eigen::LDLT<Eigen::MatrixXd> ldlt(V_i);
            if (ldlt.info() != Eigen::Success) break;
            Eigen::MatrixXd V_inv = ldlt.solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
            Eigen::VectorXd yi = Uty.row(i).transpose();
            Eigen::VectorXd v = V_inv * yi;
            std::vector<Eigen::VectorXd> w(v_size);
            for (int a = 0; a < v_size; ++a) {
                int p = pairs[a].first, q = pairs[a].second;
                if (p == q)
                    w[a] = v(p) * V_inv.col(p);
                else
                    w[a] = v(q) * V_inv.col(p) + v(p) * V_inv.col(q);
            }
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                for (int b = a; b < v_size; ++b) {
                    double qval;
                    if (j == k)
                        qval = v(j) * w[b](j);
                    else
                        qval = v(j) * w[b](k) + v(k) * w[b](j);
                    qval *= 0.5;
                    AI_final(a, b)                   += d2 * qval;
                    AI_final(a, v_size + b)          += d_i * qval;
                    AI_final(v_size + a, b)          += d_i * qval;
                    AI_final(v_size + a, v_size + b) +=       qval;
                }
            }
        }
        AI_final = 0.5 * (AI_final + AI_final.transpose()).eval();
        AI_final.diagonal().array() *= (1.0 + 1e-8);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(AI_final);
        if (ldlt.info() == Eigen::Success) {
            Eigen::MatrixXd AI_inv = ldlt.solve(Eigen::MatrixXd::Identity(n_params, n_params));
            for (int a = 0; a < v_size; ++a) {
                int j = pairs[a].first, k = pairs[a].second;
                double se_vg = std::sqrt(std::abs(AI_inv(a, a)));
                double se_ve = std::sqrt(std::abs(AI_inv(v_size + a, v_size + a)));
                Vg_se(j, k) = se_vg;
                Vg_se(k, j) = se_vg;
                Ve_se(j, k) = se_ve;
                Ve_se(k, j) = se_ve;
            }
        }
    }
    if (options.verbose) {
        std::cout << "  Vg:" << std::endl << Vg << std::endl;
        std::cout << "  Ve:" << std::endl << Ve << std::endl;
        if (Vg_se.sum() > 0 || Ve_se.sum() > 0) {
            std::cout << "  SE(Vg):" << std::endl << Vg_se << std::endl;
            std::cout << "  SE(Ve):" << std::endl << Ve_se << std::endl;
        }
    }

    if (options.verbose)
        std::cout << "  AI-REML finished." << std::endl;
    return true;
}

bool MvEigenVce::runEmmvgwasLambdaNullModel() {
    if (options.use_projection || n_covars != 1) {
        throw std::runtime_error("EMmvGWAS-compatible lambda VCE requires explicit single-column covariate mode; use --no-covar-proj and an intercept-only design.");
    }

    if (options.verbose) std::cout << "MvEigenVce: Running EMmvGWAS-compatible lambda VCE..." << std::endl;

    for (int i = 0; i < D.size(); ++i) {
        if (D(i) < 1e-5) D(i) = 1e-5;
    }

    Eigen::MatrixXd aa = Eigen::MatrixXd::Identity(n_traits, n_traits);
    Eigen::MatrixXd ee = Eigen::MatrixXd::Identity(n_traits, n_traits);
    Eigen::VectorXd mu = Eigen::VectorXd::Zero(n_traits);

    for (int iter = 0; iter < options.max_iter; ++iter) {
        Eigen::MatrixXd xx = Eigen::MatrixXd::Zero(n_traits, n_traits);
        Eigen::VectorXd xy = Eigen::VectorXd::Zero(n_traits);
        Eigen::MatrixXd yy = Eigen::MatrixXd::Zero(n_traits, n_traits);
        Eigen::MatrixXd bb = Eigen::MatrixXd::Zero(n_traits, n_traits);

        std::vector<Eigen::MatrixXd> tt_cache(n_samples);
        std::vector<Eigen::VectorXd> yh_cache(n_samples);
        std::vector<double> xh_cache(n_samples, 0.0);

        for (int j = 0; j < n_samples; ++j) {
            double w = D(j);
            double sqrt_w = std::sqrt(w);
            double xh = UtX(j, 0) / sqrt_w;
            Eigen::VectorXd yh = Uty.row(j).transpose() / sqrt_w;
            Eigen::MatrixXd tt = (aa * w + ee + Eigen::MatrixXd::Identity(n_traits, n_traits) * 1e-5).inverse();

            xx += w * xh * xh * tt;
            xy += w * xh * (tt * yh);

            tt_cache[j] = std::move(tt);
            yh_cache[j] = std::move(yh);
            xh_cache[j] = xh;
        }

        Eigen::LDLT<Eigen::MatrixXd> xx_ldlt(xx);
        if (xx_ldlt.info() != Eigen::Success) return false;
        mu = xx_ldlt.solve(xy);

        for (int j = 0; j < n_samples; ++j) {
            double w = D(j);
            const Eigen::MatrixXd& tt = tt_cache[j];
            const Eigen::VectorXd& yh = yh_cache[j];
            double xh = xh_cache[j];
            Eigen::VectorXd resid = yh - xh * mu;
            Eigen::VectorXd ex = w * aa * tt * resid;
            yy += w * (resid - ex) * yh.transpose();
            Eigen::MatrixXd va = aa - w * aa * tt * aa;
            bb += ex * ex.transpose() + va;
        }

        Eigen::MatrixXd aa_new = bb / (double)n_samples;
        Eigen::MatrixXd ee_new = yy / (double)(n_samples - 1);
        aa_new = 0.5 * (aa_new + aa_new.transpose()).eval();
        ee_new = 0.5 * (ee_new + ee_new.transpose()).eval();

        double err = ((aa_new - aa).squaredNorm() + (ee_new - ee).squaredNorm()) /
                     std::max(1.0, 2.0 * n_traits * n_traits);
        aa = aa_new;
        ee = ee_new;

        if (err < options.tol) {
            Eigen::MatrixXd ee_inv = ee.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
            emmvgwas_G = aa;
            emmvgwas_R = ee;
            emmvgwas_beta = mu;
            Vg = aa * ee_inv;
            Ve = Eigen::MatrixXd::Identity(n_traits, n_traits) * (1.0 + 1e-5);
            B = mu.transpose();
            if (options.verbose) {
                std::cout << "MvEigenVce: EMmvGWAS lambda converged in " << iter + 1 << " iterations." << std::endl;
                std::cout << "  lambda=G*R^-1:" << std::endl << Vg << std::endl;
                std::cout << "  scan residual scale:" << std::endl << Ve << std::endl;
            }
            return true;
        }
    }

    Eigen::MatrixXd ee_inv = ee.ldlt().solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
    emmvgwas_G = aa;
    emmvgwas_R = ee;
    emmvgwas_beta = mu;
    Vg = aa * ee_inv;
    Ve = Eigen::MatrixXd::Identity(n_traits, n_traits) * (1.0 + 1e-5);
    B = mu.transpose();
    if (options.verbose) std::cout << "MvEigenVce: Warning: EMmvGWAS lambda VCE did not converge within " << options.max_iter << " iterations." << std::endl;
    return false;
}

void MvEigenVce::setEmmvgwasFixedComponents(const Eigen::MatrixXd& G_component,
                                            const Eigen::MatrixXd& R_component) {
    if (G_component.rows() != n_traits || G_component.cols() != n_traits ||
        R_component.rows() != n_traits || R_component.cols() != n_traits) {
        throw std::runtime_error("MvEigenVce::setEmmvgwasFixedComponents: G/R dimensions must match phenotype trait count");
    }
    Eigen::LDLT<Eigen::MatrixXd> r_ldlt(R_component);
    if (r_ldlt.info() != Eigen::Success) {
        throw std::runtime_error("MvEigenVce::setEmmvgwasFixedComponents: R matrix is not invertible");
    }
    for (int i = 0; i < D.size(); ++i) {
        if (D(i) < 1e-5) D(i) = 1e-5;
    }
    emmvgwas_G = G_component;
    emmvgwas_R = R_component;
    emmvgwas_beta = Eigen::VectorXd::Zero(n_traits);
    Vg = G_component * r_ldlt.solve(Eigen::MatrixXd::Identity(n_traits, n_traits));
    Ve = Eigen::MatrixXd::Identity(n_traits, n_traits) * (1.0 + 1e-5);
    B = Eigen::MatrixXd::Zero(n_covars, n_traits);
    if (options.verbose) {
        std::cout << "MvEigenVce: Using fixed EMmvGWAS lambda=G*R^-1:" << std::endl
                  << Vg << std::endl;
    }
}

bool MvEigenVce::runNullModel() {
    if (options.use_emmvgwas_lambda) {
        return runEmmvgwasLambdaNullModel();
    }

    if (options.solver == VceSolver::AIREML) {
        return runAireml(options.em_min_iter, options.ai_max_iter);
    }

    if (options.verbose) std::cout << "MvEigenVce: Running Multivariate Null Model (EM)..." << std::endl;
    
    Eigen::MatrixXd Vg_new(n_traits, n_traits);
    Eigen::MatrixXd Ve_new(n_traits, n_traits);
    
    for (int iter = 0; iter < options.max_iter; ++iter) {
        updateEM(Vg_new, Ve_new);
        
        double err = ((Vg_new - Vg).squaredNorm() + (Ve_new - Ve).squaredNorm()) /
                     std::max(1.0, 2.0 * n_traits * n_traits);
        
        Vg = Vg_new;
        Ve = Ve_new;
        
        if (err < options.tol) {
            if (options.verbose) {
                std::cout << "MvEigenVce: EM Converged in " << iter + 1 << " iterations." << std::endl;
                std::cout << "  Vg:" << std::endl << Vg << std::endl;
                std::cout << "  Ve:" << std::endl << Ve << std::endl;
            }
            return true;
        }
    }
    
    if (options.verbose) std::cout << "MvEigenVce: Warning: EM did not converge within " << options.max_iter << " iterations." << std::endl;
    return false;
}

} // namespace cosmic
