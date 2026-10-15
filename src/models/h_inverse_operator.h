#pragma once

#include "a22_inverse_operator.h"

#include <Eigen/Sparse>
#include <string>
#include <vector>

namespace cosmic {

// Exact H^-1 matvec over full pedigree order. The diagonal callback is the
// Jacobi preconditioner contract: exact when A22InverseOperator has an exact
// diagonal (Nu=0), otherwise the negative A22 term is deliberately omitted.
class HInverseOperator final : public AbstractMatrix {
public:
    HInverseOperator(const Eigen::SparseMatrix<double>& a_inverse,
                     AbstractMatrix* scaled_genomic_precision,
                     const A22InverseOperator* a22_inverse,
                     const std::vector<int>& genotyped_indices,
                     double omega);

    int rows() const override { return static_cast<int>(a_inverse_.rows()); }
    int cols() const override { return static_cast<int>(a_inverse_.cols()); }
    size_t nonZeros() const override;
    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const override;
    void add_diagonal_to(Eigen::VectorXd& diag, double alpha) const override;

    Eigen::MatrixXd getBlock(int, int) const override;
    void visit_triplets(std::function<void(int, int, double)>) const override;
    Eigen::MatrixXd toDense() const override;

    bool preconditionerDiagonalExact() const { return preconditioner_diagonal_exact_; }
    std::string preconditionerDiagonalMode() const {
        return preconditioner_diagonal_exact_ ? "exact" : "approx_a22_omitted";
    }

private:
    const Eigen::SparseMatrix<double>& a_inverse_;
    AbstractMatrix* scaled_genomic_precision_ = nullptr;
    const A22InverseOperator* a22_inverse_ = nullptr;
    std::vector<int> genotyped_indices_;
    double omega_ = 1.0;
    Eigen::VectorXd preconditioner_diagonal_;
    bool preconditioner_diagonal_exact_ = false;
};

} // namespace cosmic
