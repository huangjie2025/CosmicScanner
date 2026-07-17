#include "mv_em_update.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>

static Eigen::MatrixXd old_scaled_contribution(double eigenvalue,
                                               const Eigen::MatrixXd& Vg,
                                               const Eigen::MatrixXd& sigma_inv,
                                               const Eigen::VectorXd& resid) {
    Eigen::MatrixXd Wg = Vg * eigenvalue * sigma_inv;
    Eigen::VectorXd uhat = Wg * resid;
    return uhat * uhat.transpose() + Vg * eigenvalue - Wg * Vg * eigenvalue;
}

int main() {
    Eigen::MatrixXd Vg(2, 2);
    Vg << 1.35, 0.28,
          0.28, 0.74;
    Eigen::MatrixXd Ve(2, 2);
    Ve << 0.62, -0.11,
         -0.11, 0.91;
    Eigen::VectorXd resid(2);
    resid << 0.73, -1.21;
    const double lambda = 2.7;

    Eigen::MatrixXd Sigma = lambda * Vg + Ve;
    Eigen::MatrixXd Sigma_inv = Sigma.ldlt().solve(Eigen::MatrixXd::Identity(2, 2));

    const double sqrt_lambda = std::sqrt(lambda);
    Eigen::VectorXd ghat = sqrt_lambda * Vg * Sigma_inv * resid;
    Eigen::MatrixXd Var_g = Vg - lambda * Vg * Sigma_inv * Vg;
    Eigen::MatrixXd dense_reference = ghat * ghat.transpose() + Var_g;
    dense_reference = 0.5 * (dense_reference + dense_reference.transpose());

    Eigen::MatrixXd fixed = cosmic::mv_em_genetic_contribution(lambda, Vg, Sigma_inv, resid);
    fixed = 0.5 * (fixed + fixed.transpose());
    Eigen::MatrixXd old = old_scaled_contribution(lambda, Vg, Sigma_inv, resid);
    old = 0.5 * (old + old.transpose());

    const double fixed_err = (fixed - dense_reference).cwiseAbs().maxCoeff();
    const double old_err = (old - dense_reference).cwiseAbs().maxCoeff();

    if (!(fixed_err < 1e-12)) {
        std::cerr << "fixed EM update does not match dense reference; max error="
                  << fixed_err << "\n";
        return 1;
    }
    if (!(old_err > 0.1)) {
        std::cerr << "old eigenvalue-scaled update unexpectedly matched reference; max error="
                  << old_err << "\n";
        return 1;
    }
    return 0;
}
