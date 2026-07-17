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

    std::cerr << "unknown mode\n";
    return 1;
}
