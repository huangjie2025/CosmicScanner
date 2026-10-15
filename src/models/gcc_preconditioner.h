#ifndef COSMIC_GCC_PRECONDITIONER_H
#define COSMIC_GCC_PRECONDITIONER_H

#include "gcc_solver.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <random>

namespace cosmic {

// Forward declarations
class GccOperator;

/**
 * Abstract base class for Gcc preconditioners.
 * A preconditioner M approximates Gcc such that M^-1 is easy to apply.
 */
class GccPreconditioner {
public:
    virtual ~GccPreconditioner() = default;

    /**
     * Apply the preconditioner: y = M^-1 x
     */
    virtual void apply(const Eigen::VectorXd& x, Eigen::VectorXd& y) const = 0;

    /**
     * Get the number of rows (dimension)
     */
    virtual int rows() const = 0;

    /**
     * Get the name of the preconditioner
     */
    virtual std::string name() const = 0;

    /**
     * Check if the preconditioner is valid (initialized)
     */
    virtual bool isValid() const = 0;
};

using GccPreconditionerPtr = std::shared_ptr<const GccPreconditioner>;

/**
 * Jacobi (diagonal) preconditioner.
 * M = diag(Gcc), so M^-1 = 1 / diag(Gcc)
 * Diagonal estimated via probe vectors: diag(Gcc) = E[ (Gcc * z) ⊙ z ]
 */
class JacobiPreconditioner : public GccPreconditioner {
public:
    JacobiPreconditioner() = default;

    /**
     * Initialize from a GccOperator by estimating diag(Gcc) using probe vectors.
     * @param op GccOperator to estimate diagonal from
     * @param n_probes Number of probe vectors (default 10)
     * @param seed Random seed for probe vectors
     */
    void initialize(const GccOperator& op, int n_probes = 10, unsigned seed = 12345);

    void apply(const Eigen::VectorXd& x, Eigen::VectorXd& y) const override;
    int rows() const override { return inv_diag_.size(); }
    std::string name() const override { return "Jacobi"; }
    bool isValid() const override { return inv_diag_.size() > 0; }

    const Eigen::VectorXd& invDiagonal() const { return inv_diag_; }

private:
    Eigen::VectorXd inv_diag_;  // 1 / diag(Gcc)
};

/**
  * Low-rank spectral preconditioner using top-r eigenpairs from Lanczos.
 *
  * M^-1 = U_r Lambda_r^-1 U_r^T + (I - U_r U_r^T) D^-1 (I - U_r U_r^T)
 *
 * where:
 * - U_r: top r eigenvectors from Lanczos on GccOperator
 * - Lambda_r: corresponding eigenvalues
  * - D: diagonal of Gcc (Jacobi)
 *
 * This is the "deflation" preconditioner that exactly inverts the top-r
 * eigenmodes and uses Jacobi for the remaining subspace.
 */
class LowRankSpectralPreconditioner : public GccPreconditioner {
public:
    struct Options {
        int rank = 50;              // Number of eigenvectors to compute
        int max_lanczos_iter = 80;  // Max Lanczos iterations
        double tol = 1e-6;          // Lanczos tolerance
        unsigned seed = 42;         // Random seed for Lanczos
        bool verbose = false;       // Print initialization info
    };

    LowRankSpectralPreconditioner() = default;
    explicit LowRankSpectralPreconditioner(const Options& opts) : opts_(opts) {}

    /**
     * Initialize by running Lanczos on the GccOperator to get top-r eigenpairs.
     * @param op GccOperator to compute eigenpairs from
     * @param jacobi Jacobi preconditioner (must be pre-initialized for D^-1)
     */
    void initialize(const GccOperator& op, const JacobiPreconditioner& jacobi);

    void apply(const Eigen::VectorXd& x, Eigen::VectorXd& y) const override;
    int rows() const override { return U_.rows(); }
    std::string name() const override { return "LowRankSpectral(r=" + std::to_string(rank_) + ")"; }
    bool isValid() const override { return rank_ > 0 && U_.cols() > 0; }

    int rank() const { return rank_; }
    const Eigen::MatrixXd& U() const { return U_; }
    const Eigen::VectorXd& Lambda() const { return Lambda_; }

private:
    Options opts_;
    int rank_ = 0;
    Eigen::MatrixXd U_;           // n x r eigenvectors
    Eigen::VectorXd Lambda_;      // r eigenvalues
    Eigen::VectorXd UT_Dinv_U_diag_; // diag(U^T D^-1 U)
    Eigen::VectorXd jacobi_diag_; // Cached Jacobi diagonal for apply()
};

} // namespace cosmic

#endif // COSMIC_GCC_PRECONDITIONER_H
