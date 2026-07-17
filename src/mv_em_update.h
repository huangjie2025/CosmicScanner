#pragma once

#include <Eigen/Dense>

namespace cosmic {

inline Eigen::MatrixXd mv_em_genetic_contribution(double eigenvalue,
                                                  const Eigen::MatrixXd& Vg,
                                                  const Eigen::MatrixXd& sigma_inv,
                                                  const Eigen::VectorXd& resid) {
    if (eigenvalue <= 1e-12) return Vg;
    Eigen::MatrixXd Wg = Vg * eigenvalue * sigma_inv;
    Eigen::VectorXd uhat = Wg * resid;
    return (uhat * uhat.transpose()) / eigenvalue + Vg - Wg * Vg;
}

} // namespace cosmic
