#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <vector>

namespace cosmic {

struct MvRemlEval {
    double loglik = std::numeric_limits<double>::quiet_NaN();
    double logdet_h = std::numeric_limits<double>::quiet_NaN();
    double logdet_xt_hinv_x = std::numeric_limits<double>::quiet_NaN();
    double quad = std::numeric_limits<double>::quiet_NaN();
    Eigen::VectorXd beta;
};

inline MvRemlEval mv_eval_reml_dense(const Eigen::MatrixXd& UtY,
                                     const Eigen::MatrixXd& UtX,
                                     const Eigen::VectorXd& D,
                                     const Eigen::MatrixXd& Vg,
                                     const Eigen::MatrixXd& Ve) {
    MvRemlEval out;
    const int n = UtY.rows();
    const int d = UtY.cols();
    const int c = UtX.cols();
    const int pd = c * d;
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(pd, pd);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(pd);
    double yhy = 0.0;
    double logdet_h = 0.0;
    for (int i = 0; i < n; ++i) {
        Eigen::MatrixXd S = D(i) * Vg + Ve;
        S.diagonal().array() += 1e-8;
        Eigen::LDLT<Eigen::MatrixXd> ldlt(S);
        if (ldlt.info() != Eigen::Success) return out;
        Eigen::MatrixXd Sinv = ldlt.solve(Eigen::MatrixXd::Identity(d, d));
        for (int k = 0; k < d; ++k) {
            const double val = ldlt.vectorD()(k);
            if (!(val > 0.0)) return out;
            logdet_h += std::log(val);
        }
        Eigen::VectorXd yi = UtY.row(i).transpose();
        yhy += (yi.transpose() * Sinv * yi)(0, 0);
        for (int a = 0; a < c; ++a) {
            rhs.segment(a * d, d) += UtX(i, a) * (Sinv * yi);
            for (int b = 0; b < c; ++b) {
                H.block(a * d, b * d, d, d) += UtX(i, a) * UtX(i, b) * Sinv;
            }
        }
    }
    Eigen::LDLT<Eigen::MatrixXd> hldlt(H);
    if (hldlt.info() != Eigen::Success) return out;
    out.beta = hldlt.solve(rhs);
    Eigen::VectorXd hsolve = hldlt.solve(rhs);
    double logdet_x = 0.0;
    for (int k = 0; k < hldlt.vectorD().size(); ++k) {
        const double val = hldlt.vectorD()(k);
        if (!(val > 0.0)) return out;
        logdet_x += std::log(val);
    }
    out.logdet_h = logdet_h;
    out.logdet_xt_hinv_x = logdet_x;
    out.quad = yhy - rhs.dot(hsolve);
    out.loglik = -0.5 * ((n - c) * d * std::log(2.0 * M_PI) +
                         out.logdet_h + out.logdet_xt_hinv_x + out.quad);
    return out;
}

inline Eigen::VectorXd mv_reml_numeric_gradient(const Eigen::MatrixXd& UtY,
                                                const Eigen::MatrixXd& UtX,
                                                const Eigen::VectorXd& D,
                                                const Eigen::MatrixXd& Vg,
                                                const Eigen::MatrixXd& Ve,
                                                double eps = 1e-5) {
    const int d = Vg.rows();
    const int npar = d * (d + 1);
    Eigen::VectorXd grad(npar);
    int idx = 0;
    auto perturb = [&](Eigen::MatrixXd& M, int i, int j, double h) {
        M(i, j) += h;
        if (i != j) M(j, i) += h;
    };
    for (int block = 0; block < 2; ++block) {
        for (int i = 0; i < d; ++i) {
            for (int j = i; j < d; ++j, ++idx) {
                Eigen::MatrixXd Vg_p = Vg, Vg_m = Vg, Ve_p = Ve, Ve_m = Ve;
                if (block == 0) {
                    perturb(Vg_p, i, j, eps);
                    perturb(Vg_m, i, j, -eps);
                } else {
                    perturb(Ve_p, i, j, eps);
                    perturb(Ve_m, i, j, -eps);
                }
                double fp = mv_eval_reml_dense(UtY, UtX, D, Vg_p, Ve_p).loglik;
                double fm = mv_eval_reml_dense(UtY, UtX, D, Vg_m, Ve_m).loglik;
                grad(idx) = (fp - fm) / (2.0 * eps);
            }
        }
    }
    return grad;
}

} // namespace cosmic
