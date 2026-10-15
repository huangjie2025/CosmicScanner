#pragma once

#include "matrix_adapter.h"

#include <Eigen/SparseCholesky>
#include <chrono>
#include <string>
#include <vector>

namespace cosmic {

// Exact A22^-1 action from the Schur complement of the untyped block in A^-1:
//   A22^-1 = R - Q' P^-1 Q,  A^-1 = [P Q; Q' R].
// P/Q/R remain sparse. No dense pedigree covariance/precision block is stored.
class A22InverseOperator final : public AbstractMatrix {
public:
    A22InverseOperator(const Eigen::SparseMatrix<double>& a_inverse,
                       const std::vector<int>& genotyped_indices);

    int rows() const override { return n_genotyped_; }
    int cols() const override { return n_genotyped_; }
    size_t nonZeros() const override;

    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const override;
    Eigen::MatrixXd apply(const Eigen::MatrixXd& x) const;
    void add_diagonal_to(Eigen::VectorXd& diag, double alpha) const override;

    Eigen::MatrixXd getBlock(int, int) const override;
    void visit_triplets(std::function<void(int, int, double)>) const override;
    Eigen::MatrixXd toDense() const override;

    // The Schur diagonal is cheap/exact only when P is empty. For a mixed
    // pedigree, diag(R) is exposed solely as a safe preconditioner ingredient.
    bool hasExactDiagonal() const { return n_untyped_ == 0; }
    const Eigen::VectorXd& preconditionerDiagonal() const { return r_diagonal_; }

    int nPedigree() const { return n_pedigree_; }
    int nGenotyped() const { return n_genotyped_; }
    int nUntyped() const { return n_untyped_; }
    long long pNonZeros() const { return p_.nonZeros(); }
    long long qNonZeros() const { return q_.nonZeros(); }
    long long rNonZeros() const { return r_.nonZeros(); }
    long long factorNonZeros() const { return factor_nonzeros_; }
    long long factorBytesEstimate() const { return factor_nonzeros_ * (long long)(sizeof(double) + sizeof(int)); }
    long long storageBytesEstimate() const {
        return (p_.nonZeros() + q_.nonZeros() + r_.nonZeros()) *
               (long long)(sizeof(double) + sizeof(int));
    }
    double factorTimeSeconds() const { return factor_time_s_; }
    long long solveCount() const { return solve_count_; }
    double solveTimeSeconds() const { return solve_time_s_; }
    double maxSolveRelativeResidual() const { return max_solve_relative_residual_; }
    std::string backendName() const { return "sparse-schur-a22-inverse"; }

private:
    int n_pedigree_ = 0;
    int n_genotyped_ = 0;
    int n_untyped_ = 0;
    Eigen::SparseMatrix<double> p_;
    Eigen::SparseMatrix<double> q_;
    Eigen::SparseMatrix<double> r_;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> p_solver_;
    Eigen::VectorXd r_diagonal_;
    long long factor_nonzeros_ = 0;
    double factor_time_s_ = 0.0;
    mutable long long solve_count_ = 0;
    mutable double solve_time_s_ = 0.0;
    mutable double max_solve_relative_residual_ = 0.0;
};

} // namespace cosmic
