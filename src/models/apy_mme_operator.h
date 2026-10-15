#pragma once

#include <stdexcept>

#include "apy_inverse_operator.h"
#include "matrix_adapter.h"

namespace cosmic {

// Matrix-free, scaled APY precision adapter for ImplicitMME. Explicit sparse
// assembly remains the responsibility of the dense ApyGinvAdapter path.
class APYMMEOperator final : public AbstractMatrix {
public:
    APYMMEOperator(const ApyInverseOperator* inverse, double scale)
        : inverse_(inverse), scale_(scale) {
        if (!inverse_) throw std::runtime_error("APYMMEOperator requires an APY inverse operator");
    }

    int rows() const override { return inverse_->rows(); }
    int cols() const override { return inverse_->rows(); }
    size_t nonZeros() const override { return 0; }

    Eigen::VectorXd operator*(const Eigen::VectorXd& x) const override {
        return scale_ * inverse_->apply(x);
    }

    void add_diagonal_to(Eigen::VectorXd& diag, double alpha) const override {
        if (diag.size() < rows()) {
            throw std::runtime_error("APYMMEOperator diagonal destination is too small");
        }
        diag.head(rows()) += alpha * scale_ * inverse_->diagonal();
    }

    Eigen::MatrixXd getBlock(int, int) const override {
        throw std::runtime_error("APYMMEOperator::getBlock is unavailable in matrix-free mode");
    }

    void visit_triplets(std::function<void(int, int, double)>) const override {
        throw std::runtime_error("APYMMEOperator cannot be used for explicit sparse MME assembly");
    }

    Eigen::MatrixXd toDense() const override {
        throw std::runtime_error("APYMMEOperator::toDense is unavailable in matrix-free mode");
    }

private:
    const ApyInverseOperator* inverse_;
    double scale_;
};

} // namespace cosmic
