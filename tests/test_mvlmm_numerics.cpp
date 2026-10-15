#include "mvlmm_numerics.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using cosmic::mvlmm_numerics::MatrixStatus;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
} // namespace

int main() {
    Eigen::MatrixXd covariance;
    double minimum = 0.0;

    Eigen::Matrix2d well_conditioned;
    well_conditioned << 4.0, 0.5, 0.5, 2.0;
    require(cosmic::mvlmm_numerics::invertInformation(well_conditioned, covariance, minimum) == MatrixStatus::OK,
            "well-conditioned information was rejected");
    require(covariance.allFinite() && (covariance.diagonal().array() > 0.0).all(),
            "well-conditioned standard errors are not identifiable");

    Eigen::Matrix2d near_singular = Eigen::Matrix2d::Zero();
    near_singular.diagonal() << 1.0, 1e-10;
    require(cosmic::mvlmm_numerics::invertInformation(near_singular, covariance, minimum) == MatrixStatus::OK,
            "identifiable near-singular information was rejected");
    require(covariance.allFinite() && (covariance.diagonal().array() > 0.0).all(),
            "near-singular standard errors are not finite");

    Eigen::MatrixXd boundary = Eigen::MatrixXd::Zero(2, 2);
    boundary(0, 0) = 1.0;
    require(cosmic::mvlmm_numerics::conditionCovariance(boundary, minimum) == MatrixStatus::BOUNDARY_SOLUTION,
            "PSD boundary covariance was not classified");

    Eigen::MatrixXd roundoff = Eigen::MatrixXd::Zero(2, 2);
    roundoff.diagonal() << 1.0, -1e-15;
    require(cosmic::mvlmm_numerics::conditionCovariance(roundoff, minimum) == MatrixStatus::NUMERICAL_ROUNDOFF,
            "tiny negative covariance eigenvalue was not classified as roundoff");
    require(roundoff.diagonal().minCoeff() >= 0.0, "roundoff correction did not restore PSD covariance");

    Eigen::Matrix2d singular;
    singular << 1.0, 1.0, 1.0, 1.0;
    require(cosmic::mvlmm_numerics::invertInformation(singular, covariance, minimum) == MatrixStatus::TRUE_SINGULARITY,
            "rank-deficient information was not classified as singular");

    std::cout << "mvLMM covariance numerical tests passed\n";
    return 0;
}
