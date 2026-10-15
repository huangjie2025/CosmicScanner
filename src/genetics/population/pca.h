#pragma once
#include <Eigen/Dense>

namespace cosmic {
namespace genetics {

struct PCAResult {
    Eigen::MatrixXd scores;
    Eigen::VectorXd eigenvalues;
};

class PCA {
public:
    // Nonempty finite symmetric PSD GRM, 1 <= k <= n. Scores are unit
    // eigenvectors (not scaled), descending by eigenvalue. The first largest-
    // magnitude entry in each vector is positive. Relative tolerance: 1e-10.
    static PCAResult computeFromGrm(const Eigen::MatrixXd& grm, int k);
};

} // namespace genetics
} // namespace cosmic
