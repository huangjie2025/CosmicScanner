#include "pca.h"
#include <algorithm>
#include <stdexcept>

namespace cosmic {
namespace genetics {

PCAResult PCA::computeFromGrm(const Eigen::MatrixXd& grm, int k) {
    if (grm.rows() == 0 || grm.rows() != grm.cols() ||
        k < 1 || k > grm.rows() || !grm.allFinite()) {
        throw std::invalid_argument("PCA: expected a finite square GRM and 1 <= k <= n");
    }
    // Scale before symmetry checks and decomposition to avoid overflow.
    const double scale = grm.cwiseAbs().maxCoeff();
    Eigen::MatrixXd matrix = grm;
    if (scale > 0.0) matrix /= scale;
    constexpr double tolerance = 1e-10;
    if ((matrix - matrix.transpose()).cwiseAbs().maxCoeff() > tolerance) {
        throw std::invalid_argument("PCA: GRM must be symmetric");
    }
    matrix = (0.5 * matrix + 0.5 * matrix.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("PCA: eigensolver failed");
    }
    const double spectral_scale = std::max(1.0, solver.eigenvalues().cwiseAbs().maxCoeff());
    if (solver.eigenvalues().minCoeff() < -tolerance * spectral_scale) {
        throw std::invalid_argument("PCA: GRM must be positive semidefinite");
    }
    PCAResult result;
    result.scores.resize(grm.rows(), k);
    result.eigenvalues.resize(k);
    for (int j = 0; j < k; ++j) {
        const Eigen::Index source = grm.rows() - 1 - j;
        result.eigenvalues[j] = std::max(0.0, solver.eigenvalues()[source]) * scale;
        result.scores.col(j) = solver.eigenvectors().col(source);
        Eigen::Index pivot;
        result.scores.col(j).cwiseAbs().maxCoeff(&pivot);
        if (result.scores(pivot, j) < 0.0) result.scores.col(j) *= -1.0;
    }
    if (!result.eigenvalues.allFinite()) {
        throw std::runtime_error("PCA: eigenvalues exceed double precision range");
    }
    return result;
}

} // namespace genetics
} // namespace cosmic
