#ifndef COSMIC_GCC_SOLVER_PCG_H
#define COSMIC_GCC_SOLVER_PCG_H

#include "gcc_solver.h"
#include "gcc_preconditioner.h"
#include "plink_reader.h"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <functional>

namespace cosmic {

// Forward declaration
class PlinkReader;

// Matrix-free operator for Gcc (core block of genomic relationship matrix)
// Gcc = Z_c Z_c' / S
// where Z_c is the standardized genotype matrix for core animals only
class GccOperator {
public:
    GccOperator(PlinkReader& reader,
                const std::vector<int>& core_indices,
                double S_2pq,
                double tune_a = 0.0,
                double tune_b = 1.0,
                double blend_alpha = 0.0,
                const Eigen::MatrixXd* blend_target = nullptr,
                int n_threads = 1);

    void apply(const Eigen::VectorXd& v, Eigen::VectorXd& y) const;
    int rows() const { return n_core_; }
    int cols() const { return n_core_; }

    const std::vector<int>& coreIndices() const { return core_indices_; }
    double S_2pq() const { return S_2pq_; }
    double tune_a() const { return tune_a_; }
    double tune_b() const { return tune_b_; }
    double blend_alpha() const { return blend_alpha_; }

private:
    PlinkReader& reader_;
    std::vector<int> core_indices_;
    double S_2pq_;
    double tune_a_;
    double tune_b_;
    double blend_alpha_;
    const Eigen::MatrixXd* blend_target_;
    int n_threads_;
    int n_core_;
};

// Options for PCG solver
struct GccSolverPCGOptions {
    int max_iter = 2000;
    double tolerance = 1e-9;
    int n_threads = 1;
     bool verbose = false;

    // Preconditioner type: "jacobi" or "lowrank"
    std::string preconditioner = "jacobi";
    int preconditioner_rank = 50;  // For lowrank preconditioner
};

// PCG solver for Gcc x = rhs using the matrix-free operator
class GccSolverPCG : public GccSolver {
public:
     GccSolverPCG() = default;

    // Initialize with operator (lazy - doesn't compute anything yet)
    void compute(const GccOperator& op, const GccSolverPCGOptions& opts = GccSolverPCGOptions());

    // GccSolver interface
    Eigen::VectorXd solve(const Eigen::VectorXd& rhs) const override;
    Eigen::MatrixXd solve(const Eigen::MatrixXd& rhs) const override;
    int rows() const override { return op_ ? op_->rows() : 0; }
    int cols() const override { return op_ ? op_->cols() : 0; }
    bool isValid() const override { return op_ != nullptr; }
    std::string backendName() const override { return "PCG"; }

    // Diagonal of Gcc^-1 (approximated via probe vectors)
     Eigen::VectorXd diagonal() const override;

    // LDLT pivots not available for PCG
    Eigen::VectorXd vectorD() const override { return Eigen::VectorXd(); }
    double minDiagD() const override { return 0.0; }
    double maxDiagD() const override { return 0.0; }

private:
    std::shared_ptr<const GccOperator> op_;
    GccSolverPCGOptions opts_;
     mutable std::shared_ptr<const GccPreconditioner> preconditioner_;

    void computePreconditioner() const;
    Eigen::VectorXd pcgSolve(const Eigen::VectorXd& rhs) const;
};

using GccSolverPCGPtr = std::shared_ptr<GccSolverPCG>;

} // namespace cosmic

#endif // COSMIC_GCC_SOLVER_PCG_H
