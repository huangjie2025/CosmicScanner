#include "gcc_preconditioner.h"
#include "gcc_solver_pcg.h"
#include "lanczos_eigensolver.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace cosmic {

// ============================================================
// JacobiPreconditioner implementation
// ============================================================

void JacobiPreconditioner::initialize(const GccOperator& op, int n_probes, unsigned seed) {
    const int n = op.rows();
    if (n <= 0) return;

    Eigen::VectorXd diag_sum = Eigen::VectorXd::Zero(n);
    std::mt19937 rng(seed);
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int p = 0; p < n_probes; ++p) {
        Eigen::VectorXd z(n);
        for (int i = 0; i < n; ++i) z(i) = dist(rng);
        Eigen::VectorXd y(n);
        op.apply(z, y);
        diag_sum += y.cwiseProduct(z);
     }

     Eigen::VectorXd diag = diag_sum / n_probes;

     // Ensure positive diagonal with fallback
    inv_diag_.resize(n);
    for (int i = 0; i < n; ++i) {
        if (diag(i) > 0.0 && std::isfinite(diag(i))) {
            inv_diag_(i) = 1.0 / diag(i);
        } else {
            inv_diag_(i) = 1.0;  // Identity fallback
        }
    }
}

void JacobiPreconditioner::apply(const Eigen::VectorXd& x, Eigen::VectorXd& y) const {
    if (x.size() != inv_diag_.size()) {
        throw std::runtime_error("JacobiPreconditioner: size mismatch");
    }
    y = inv_diag_.cwiseProduct(x);
}

// ============================================================
// LowRankSpectralPreconditioner implementation
// ============================================================

void LowRankSpectralPreconditioner::initialize(const GccOperator& op,
                                                const JacobiPreconditioner& jacobi) {
    const int n = op.rows();
    if (n <= 0) return;

    // Create a wrapper operator for Lanczos that matches the GccOperator interface
    class OperatorWrapper {
    public:
        OperatorWrapper(const GccOperator& op) : op_(op) {}
        void apply(const Eigen::VectorXd& v, Eigen::VectorXd& y) const {
            op_.apply(v, y);
        }
        int rows() const { return op_.rows(); }
    private:
        const GccOperator& op_;
     };

     OperatorWrapper wrapper(op);

     // Run Lanczos to get top-r eigenpairs
    LanczosEigenSolver::Options lopts;
    lopts.max_iterations = opts_.max_lanczos_iter;
    lopts.tolerance = opts_.tol;
    lopts.seed = opts_.seed;
     lopts.n_threads = 1;  // Single-threaded for reproducibility

     LanczosEigenSolver solver;
     LanczosResult result = solver.compute(wrapper, 0.99, lopts);

     if (result.status != "OK") {
        std::cerr << "Warning: Lanczos failed (" << result.status << "), "
                  << "falling back to Jacobi only\n";
        rank_ = 0;
        return;
     }

     int actual_rank = std::min(opts_.rank, result.convergedritz);
    if (actual_rank <= 0) {
        rank_ = 0;
        return;
     }

     // Get the top-r eigenvalues and eigenvectors
    Lambda_ = result.eigenvalues.head(actual_rank);
    U_ = result.eigenvectors.leftCols(actual_rank);
     rank_ = actual_rank;

     // Compute diag(U^T D^-1 U) for the correction term
    const Eigen::VectorXd& D_inv = jacobi.invDiagonal();
     UT_Dinv_U_diag_ = (U_.cwiseProduct(D_inv.asDiagonal() * U_)).rowwise().sum();

     // Store Jacobi diagonal for apply()
     jacobi_diag_ = D_inv;

     if (opts_.verbose) {
        std::cout << "LowRankSpectralPreconditioner: rank=" << rank_
                  << ", min eig=" << Lambda_.minCoeff()
                  << ", max eig=" << Lambda_.maxCoeff() << "\n";
    }
}

void LowRankSpectralPreconditioner::apply(const Eigen::VectorXd& x, Eigen::VectorXd& y) const {
    if (!isValid()) {
        // Should not reach here if properly checked, but fallback to identity
        y = x;
        return;
     }

     // Apply M^-1 = U Lambda^-1 U^T + (I - U U^T) D^-1 (I - U U^T)
     // where D^-1 is the Jacobi diagonal

     // Compute U^T x
     Eigen::VectorXd UT_x = U_.transpose() * x;

     // Low-rank part: U * Lambda^-1 * U^T * x
     Eigen::VectorXd y_lr = U_ * (Lambda_.cwiseInverse().cwiseProduct(UT_x));

     // Diagonal correction: (I - U U^T) D^-1 (I - U U^T) x
    // = D^-1 (x - U UT_x) - U (U^T D^-1 (x - U UT_x))
     // = D^-1 * (x - U * UT_x) - U * (U^T * D^-1 * (x - U * UT_x))

     // Compute residual = x - U * UT_x
     Eigen::VectorXd residual = x - U_ * UT_x;

     // D^-1 * residual
     Eigen::VectorXd D_inv_residual = jacobi_diag_.cwiseProduct(residual);

     // U^T * D^-1 * residual
     Eigen::VectorXd UT_D_inv_residual = U_.transpose() * D_inv_residual;

     // U * (U^T * D^-1 * residual)
     Eigen::VectorXd U_UT_D_inv_residual = U_ * UT_D_inv_residual;

     // Diagonal correction = D^-1 * residual - U * (U^T * D^-1 * residual)
     Eigen::VectorXd y_diag = D_inv_residual - U_UT_D_inv_residual;

     // Total: low-rank part + diagonal correction
    y = y_lr + y_diag;
}

} // namespace cosmic
