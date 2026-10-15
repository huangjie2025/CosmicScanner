#include "h_inverse_operator.h"

#include <stdexcept>

namespace cosmic {

HInverseOperator::HInverseOperator(
    const Eigen::SparseMatrix<double>& a_inverse,
    AbstractMatrix* scaled_genomic_precision,
    const A22InverseOperator* a22_inverse,
    const std::vector<int>& genotyped_indices,
    double omega)
    : a_inverse_(a_inverse),
      scaled_genomic_precision_(scaled_genomic_precision),
      a22_inverse_(a22_inverse),
      genotyped_indices_(genotyped_indices),
      omega_(omega) {
    if (a_inverse_.rows() != a_inverse_.cols() || !scaled_genomic_precision_ || !a22_inverse_) {
        throw std::runtime_error("HInverseOperator requires square A-inverse and non-null correction operators.");
    }
    if (scaled_genomic_precision_->rows() != static_cast<int>(genotyped_indices_.size()) ||
        a22_inverse_->rows() != static_cast<int>(genotyped_indices_.size())) {
        throw std::runtime_error("HInverseOperator genotyped correction dimension mismatch.");
    }

    std::vector<char> seen(a_inverse_.rows(), 0);
    for (int idx : genotyped_indices_) {
        if (idx < 0 || idx >= a_inverse_.rows() || seen[idx]) {
            throw std::runtime_error("HInverseOperator mapping must contain unique, in-range indices.");
        }
        seen[idx] = 1;
    }

    preconditioner_diagonal_ = a_inverse_.diagonal();
    Eigen::VectorXd genomic_diag = Eigen::VectorXd::Zero(genotyped_indices_.size());
    scaled_genomic_precision_->add_diagonal_to(genomic_diag, 1.0);
    if (a22_inverse_->hasExactDiagonal()) {
        genomic_diag.noalias() -= omega_ * a22_inverse_->preconditionerDiagonal();
        preconditioner_diagonal_exact_ = true;
    }
    for (int g = 0; g < static_cast<int>(genotyped_indices_.size()); ++g) {
        preconditioner_diagonal_(genotyped_indices_[g]) += genomic_diag(g);
    }
}

size_t HInverseOperator::nonZeros() const {
    return static_cast<size_t>(a_inverse_.nonZeros()) +
           scaled_genomic_precision_->nonZeros() + a22_inverse_->nonZeros();
}

Eigen::VectorXd HInverseOperator::operator*(const Eigen::VectorXd& x) const {
    if (x.size() != a_inverse_.rows()) {
        throw std::runtime_error("HInverseOperator vector dimension mismatch.");
    }
    Eigen::VectorXd y = a_inverse_ * x;
    Eigen::VectorXd xg(genotyped_indices_.size());
    for (int g = 0; g < xg.size(); ++g) xg(g) = x(genotyped_indices_[g]);
    Eigen::VectorXd correction = scaled_genomic_precision_->operator*(xg);
    correction.noalias() -= omega_ * a22_inverse_->operator*(xg);
    for (int g = 0; g < correction.size(); ++g) y(genotyped_indices_[g]) += correction(g);
    return y;
}

void HInverseOperator::add_diagonal_to(Eigen::VectorXd& diag, double alpha) const {
    if (diag.size() != preconditioner_diagonal_.size()) {
        throw std::runtime_error("HInverseOperator diagonal dimension mismatch.");
    }
    diag.noalias() += alpha * preconditioner_diagonal_;
}

Eigen::MatrixXd HInverseOperator::getBlock(int, int) const {
    throw std::runtime_error("HInverseOperator is operator-only; dense block extraction is disabled.");
}

void HInverseOperator::visit_triplets(std::function<void(int, int, double)>) const {
    throw std::runtime_error("HInverseOperator is operator-only; triplet materialization is disabled.");
}

Eigen::MatrixXd HInverseOperator::toDense() const {
    throw std::runtime_error("HInverseOperator is operator-only; dense materialization is disabled.");
}

} // namespace cosmic
