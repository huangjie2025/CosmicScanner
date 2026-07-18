#include "mv_reml_dense.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>

static void make_data(Eigen::MatrixXd& Y, Eigen::MatrixXd& X, Eigen::MatrixXd& K) {
    Y.resize(5, 2);
    Y << 0.5, -0.7,
         1.2,  0.1,
        -0.4,  0.8,
         0.9, -1.1,
        -0.8,  0.4;
    X.resize(5, 2);
    X << 1.0, -0.8,
         1.0, -0.2,
         1.0,  0.0,
         1.0,  0.4,
         1.0,  1.1;
    Eigen::MatrixXd Z(5, 3);
    Z << 0.9, -0.2, 0.4,
        -0.3,  0.8, 0.1,
         0.2,  0.1, 0.7,
        -0.8, -0.4, 0.2,
         0.0,  0.6, -0.5;
    K = Z * Z.transpose() / 3.0;
}

static Eigen::MatrixXd full_cov(const Eigen::VectorXd& D,
                                const Eigen::MatrixXd& U,
                                const Eigen::MatrixXd& Vg,
                                const Eigen::MatrixXd& Ve) {
    const int n = D.size();
    const int d = Vg.rows();
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n * d, n * d);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Eigen::MatrixXd block = (U.row(i).dot(D.asDiagonal() * U.row(j).transpose())) * Vg;
            if (i == j) block += Ve;
            for (int t = 0; t < d; ++t) {
                for (int s = 0; s < d; ++s) {
                    V(i + t * n, j + s * n) = block(t, s);
                }
            }
        }
    }
    return V;
}

static bool is_spd(const Eigen::MatrixXd& M) {
    Eigen::LLT<Eigen::MatrixXd> llt(M);
    return llt.info() == Eigen::Success;
}

static std::vector<double> pack_spd_cholesky(const Eigen::MatrixXd& M) {
    Eigen::LLT<Eigen::MatrixXd> llt(M + 1e-10 * Eigen::MatrixXd::Identity(M.rows(), M.cols()));
    Eigen::MatrixXd L = llt.matrixL();
    std::vector<double> out;
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j <= i; ++j) {
            out.push_back(i == j ? std::log(std::max(L(i, j), 1e-12)) : L(i, j));
        }
    }
    return out;
}

static Eigen::MatrixXd unpack_spd_cholesky(const std::vector<double>& x, int offset, int d) {
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(d, d);
    int k = offset;
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j <= i; ++j, ++k) {
            L(i, j) = i == j ? std::exp(x[k]) : x[k];
        }
    }
    return L * L.transpose();
}

static double optimize_dense_reml(const Eigen::MatrixXd& UtY,
                                  const Eigen::MatrixXd& UtX,
                                  const Eigen::VectorXd& D,
                                  Eigen::MatrixXd& Vg,
                                  Eigen::MatrixXd& Ve) {
    const int d = Vg.rows();
    const int block = d * (d + 1) / 2;
    std::vector<double> x = pack_spd_cholesky(Vg);
    std::vector<double> xe = pack_spd_cholesky(Ve);
    x.insert(x.end(), xe.begin(), xe.end());
    const int npar = static_cast<int>(x.size());
    auto objective = [&](const std::vector<double>& p) {
        Eigen::MatrixXd Vg2 = unpack_spd_cholesky(p, 0, d);
        Eigen::MatrixXd Ve2 = unpack_spd_cholesky(p, block, d);
        return -cosmic::mv_eval_reml_dense(UtY, UtX, D, Vg2, Ve2).loglik;
    };
    auto gradient = [&](const std::vector<double>& p) {
        Eigen::VectorXd g(npar);
        for (int k = 0; k < npar; ++k) {
            const double h = 1e-5 * std::max(1.0, std::abs(p[k]));
            std::vector<double> xp = p, xm = p;
            xp[k] += h;
            xm[k] -= h;
            g(k) = (objective(xp) - objective(xm)) / (2.0 * h);
        }
        return g;
    };

    double f = objective(x);
    Eigen::VectorXd g = gradient(x);
    Eigen::MatrixXd Hinv = Eigen::MatrixXd::Identity(npar, npar);
    for (int iter = 0; iter < 100; ++iter) {
        if (g.norm() < 1e-6) break;
        Eigen::VectorXd direction = -Hinv * g;
        if (!direction.allFinite() || direction.dot(g) >= 0.0) {
            direction = -g;
            Hinv.setIdentity();
        }
        double alpha = 1.0;
        std::vector<double> xnew(npar);
        double fnew = std::numeric_limits<double>::infinity();
        for (int ls = 0; ls < 30; ++ls) {
            for (int k = 0; k < npar; ++k) xnew[k] = x[k] + alpha * direction(k);
            fnew = objective(xnew);
            if (std::isfinite(fnew) && fnew < f) break;
            alpha *= 0.5;
        }
        if (!std::isfinite(fnew) || fnew >= f) break;
        Eigen::VectorXd gnew = gradient(xnew);
        Eigen::VectorXd s(npar);
        for (int k = 0; k < npar; ++k) s(k) = xnew[k] - x[k];
        Eigen::VectorXd y = gnew - g;
        double ys = y.dot(s);
        if (ys > 1e-12) {
            double rho = 1.0 / ys;
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(npar, npar);
            Hinv = (I - rho * s * y.transpose()) * Hinv * (I - rho * y * s.transpose()) +
                   rho * s * s.transpose();
        } else {
            Hinv.setIdentity();
        }
        x = xnew;
        f = fnew;
        g = gnew;
    }
    Vg = unpack_spd_cholesky(x, 0, d);
    Ve = unpack_spd_cholesky(x, block, d);
    return -f;
}

int main(int argc, char** argv) {
    std::string mode = argc > 1 ? argv[1] : "objective";
    Eigen::MatrixXd Y, X, K;
    make_data(Y, X, K);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
    Eigen::VectorXd D = es.eigenvalues().cwiseMax(0.0);
    Eigen::MatrixXd U = es.eigenvectors();
    Eigen::MatrixXd UtY = U.transpose() * Y;
    Eigen::MatrixXd UtX = U.transpose() * X;
    Eigen::MatrixXd Vg(2, 2), Ve(2, 2);
    Vg << 0.75, 0.12,
          0.12, 0.45;
    Ve << 0.62, -0.08,
         -0.08, 0.81;

    if (mode == "objective") {
        cosmic::MvRemlEval e = cosmic::mv_eval_reml_dense(UtY, UtX, D, Vg, Ve);
        Eigen::MatrixXd V = full_cov(D, U, Vg, Ve);
        Eigen::VectorXd yvec(Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size()));
        Eigen::MatrixXd Xfull = Eigen::MatrixXd::Zero(Y.size(), X.cols() * Y.cols());
        for (int i = 0; i < Y.rows(); ++i) {
            for (int c = 0; c < X.cols(); ++c) {
                for (int t = 0; t < Y.cols(); ++t) Xfull(i + t * Y.rows(), c * Y.cols() + t) = X(i, c);
            }
        }
        Eigen::LDLT<Eigen::MatrixXd> vldlt(V);
        Eigen::MatrixXd Vinv = vldlt.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));
        Eigen::MatrixXd XtVinvX = Xfull.transpose() * Vinv * Xfull;
        Eigen::LDLT<Eigen::MatrixXd> xldlt(XtVinvX);
        Eigen::VectorXd beta = xldlt.solve(Xfull.transpose() * Vinv * yvec);
        Eigen::VectorXd resid = yvec - Xfull * beta;
        double logdetV = vldlt.vectorD().array().log().sum();
        double logdetX = xldlt.vectorD().array().log().sum();
        double quad = resid.dot(Vinv * resid);
        double ref = -0.5 * ((Y.rows() - X.cols()) * Y.cols() * std::log(2.0 * M_PI) +
                             logdetV + logdetX + quad);
        if (std::abs(e.loglik - ref) > 1e-8) {
            std::cerr << "dense REML objective mismatch " << e.loglik << " vs " << ref << "\n";
            return 1;
        }
        return 0;
    }

    if (mode == "gradient") {
        Eigen::VectorXd g1 = cosmic::mv_reml_numeric_gradient(UtY, UtX, D, Vg, Ve, 1e-5);
        Eigen::VectorXd g2 = cosmic::mv_reml_numeric_gradient(UtY, UtX, D, Vg, Ve, 5e-6);
        if ((g1 - g2).cwiseAbs().maxCoeff() > 1e-4) {
            std::cerr << "finite-difference gradient unstable\n";
            return 1;
        }
        return 0;
    }

    if (mode == "monotonicity") {
        Eigen::VectorXd grad = cosmic::mv_reml_numeric_gradient(UtY, UtX, D, Vg, Ve, 1e-5);
        const double base = cosmic::mv_eval_reml_dense(UtY, UtX, D, Vg, Ve).loglik;
        Eigen::MatrixXd Vg2 = Vg, Ve2 = Ve;
        const double step = 1e-4;
        int idx = 0;
        for (int block = 0; block < 2; ++block) {
            for (int i = 0; i < 2; ++i) {
                for (int j = i; j < 2; ++j, ++idx) {
                    Eigen::MatrixXd& M = block == 0 ? Vg2 : Ve2;
                    M(i, j) += step * grad(idx);
                    if (i != j) M(j, i) = M(i, j);
                }
            }
        }
        const double improved = cosmic::mv_eval_reml_dense(UtY, UtX, D, Vg2, Ve2).loglik;
        if (!(improved > base)) {
            std::cerr << "small finite-difference gradient step did not improve REML\n";
            return 1;
        }
        return 0;
    }

    if (mode == "multistart") {
        std::vector<Eigen::MatrixXd> starts_g(3), starts_e(3);
        starts_g[0] = Vg;
        starts_e[0] = Ve;
        starts_g[1].resize(2, 2);
        starts_e[1].resize(2, 2);
        starts_g[1] << 1.1, -0.2,
                      -0.2,  0.3;
        starts_e[1] << 0.8,  0.1,
                       0.1,  1.0;
        starts_g[2] = 0.6 * Vg;
        starts_e[2] = 1.4 * Ve;
        std::vector<double> logliks;
        for (int r = 0; r < 3; ++r) {
            double ll = optimize_dense_reml(UtY, UtX, D, starts_g[r], starts_e[r]);
            Eigen::VectorXd grad = cosmic::mv_reml_numeric_gradient(UtY, UtX, D, starts_g[r], starts_e[r], 1e-5);
            if (grad.norm() > 5e-4) {
                std::cerr << "dense multistart gradient too large: " << grad.norm() << "\n";
                return 1;
            }
            logliks.push_back(ll);
        }
        for (size_t i = 1; i < logliks.size(); ++i) {
            if (std::abs(logliks[i] - logliks[0]) > 1e-4) {
                std::cerr << "dense multistart log-likelihood mismatch\n";
                return 1;
            }
        }
        return 0;
    }

    std::cerr << "unknown mode\n";
    return 1;
}
