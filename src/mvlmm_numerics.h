#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <limits>

namespace cosmic::mvlmm_numerics {

enum class MatrixStatus {
    OK,
    NUMERICAL_ROUNDOFF,
    BOUNDARY_SOLUTION,
    TRUE_SINGULARITY
};

inline const char* statusName(MatrixStatus status) {
    switch (status) {
        case MatrixStatus::OK: return "OK";
        case MatrixStatus::NUMERICAL_ROUNDOFF: return "NUMERICAL_ROUNDOFF";
        case MatrixStatus::BOUNDARY_SOLUTION: return "BOUNDARY_SOLUTION";
        case MatrixStatus::TRUE_SINGULARITY: return "TRUE_SINGULARITY";
    }
    return "TRUE_SINGULARITY";
}

inline MatrixStatus invertInformation(const Eigen::MatrixXd& information,
                                      Eigen::MatrixXd& covariance,
                                      double& min_eigenvalue) {
    min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    if (information.rows() == 0 || information.rows() != information.cols() || !information.allFinite())
        return MatrixStatus::TRUE_SINGULARITY;
    const Eigen::MatrixXd symmetric = 0.5 * (information + information.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetric);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
        return MatrixStatus::TRUE_SINGULARITY;
    min_eigenvalue = eigen.eigenvalues().minCoeff();
    const double scale = std::max(1.0, eigen.eigenvalues().cwiseAbs().maxCoeff());
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale * information.rows();
    if (!(min_eigenvalue > tolerance)) return MatrixStatus::TRUE_SINGULARITY;
    covariance = eigen.eigenvectors() * eigen.eigenvalues().cwiseInverse().asDiagonal() *
                 eigen.eigenvectors().transpose();
    covariance = 0.5 * (covariance + covariance.transpose());
    return covariance.allFinite() ? MatrixStatus::OK : MatrixStatus::TRUE_SINGULARITY;
}

inline MatrixStatus conditionCovariance(Eigen::MatrixXd& covariance,
                                        double& min_eigenvalue) {
    min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
    if (covariance.rows() == 0 || covariance.rows() != covariance.cols() || !covariance.allFinite())
        return MatrixStatus::TRUE_SINGULARITY;
    covariance = 0.5 * (covariance + covariance.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(covariance);
    if (eigen.info() != Eigen::Success || !eigen.eigenvalues().allFinite())
        return MatrixStatus::TRUE_SINGULARITY;
    min_eigenvalue = eigen.eigenvalues().minCoeff();
    const double scale = std::max(1.0, eigen.eigenvalues().cwiseAbs().maxCoeff());
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * scale * covariance.rows();
    if (min_eigenvalue < -tolerance) return MatrixStatus::TRUE_SINGULARITY;
    if (min_eigenvalue < 0.0) {
        covariance = eigen.eigenvectors() * eigen.eigenvalues().cwiseMax(0.0).asDiagonal() *
                     eigen.eigenvectors().transpose();
        covariance = 0.5 * (covariance + covariance.transpose());
        return MatrixStatus::NUMERICAL_ROUNDOFF;
    }
    if (min_eigenvalue <= tolerance) return MatrixStatus::BOUNDARY_SOLUTION;
    return MatrixStatus::OK;
}

} // namespace cosmic::mvlmm_numerics
