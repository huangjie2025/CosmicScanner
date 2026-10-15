#include "apy_core_pev_operator.h"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "cross_genomic_operator.h"

namespace cosmic {

namespace {

Eigen::VectorXd makeSchurWeights(const Eigen::VectorXd& mnn,
                                 const Eigen::VectorXd& dn,
                                 double lambda) {
    if (!(lambda > 0.0) || !std::isfinite(lambda)) {
        throw std::runtime_error("APY core PEV requires finite positive lambda");
    }
    if (mnn.size() != dn.size()) {
        throw std::runtime_error("APY core PEV Mnn/Dn size mismatch");
    }
    if (!mnn.allFinite() || !dn.allFinite() ||
        (mnn.array() <= 0.0).any() || (dn.array() < 0.0).any()) {
        throw std::runtime_error("APY core PEV requires positive Mnn and nonnegative Dn");
    }
    const Eigen::VectorXd denominator =
        (lambda + dn.array() * mnn.array()).matrix();
    if ((denominator.array() <= 0.0).any()) {
        throw std::runtime_error("APY core PEV Schur denominator is not positive");
    }
    return (lambda * dn.array() / denominator.array()).matrix();
}

void validateDc(const Eigen::VectorXd& dc) {
    if (dc.size() == 0 || !dc.allFinite() || (dc.array() < 0.0).any()) {
        throw std::runtime_error("APY core PEV requires nonempty nonnegative Dc");
    }
}

} // namespace

DenseAPYCorePEVOperator::DenseAPYCorePEVOperator(
    const ApyGinv& apy, const Eigen::VectorXd& dc,
    const Eigen::VectorXd& dn, double lambda) {
    validateDc(dc);
    if (apy.nCore() != dc.size() || apy.nNoncore() != dn.size()) {
        throw std::runtime_error("DenseAPYCorePEVOperator partition size mismatch");
    }
    schur_weights_ = makeSchurWeights(apy.mnn(), dn, lambda);
    coefficient_ = dc.asDiagonal().toDenseMatrix();
    coefficient_.noalias() += lambda * apy.gccInvFull();
    if (apy.nNoncore() > 0) {
        coefficient_.noalias() +=
            apy.q() * schur_weights_.asDiagonal() * apy.q().transpose();
    }
    coefficient_ = 0.5 * (coefficient_ + coefficient_.transpose()).eval();
    Eigen::LDLT<Eigen::MatrixXd> spd_check(coefficient_);
    if (spd_check.info() != Eigen::Success ||
        spd_check.vectorD().minCoeff() <= 0.0) {
        throw std::runtime_error(
            "DenseAPYCorePEVOperator Eq. 8 coefficient is not SPD");
    }
    inverse_ = std::make_unique<DensePEVOperator>(coefficient_);
}

Eigen::VectorXd DenseAPYCorePEVOperator::apply(const Eigen::VectorXd& x) const {
    return inverse_->apply(x);
}

Eigen::MatrixXd DenseAPYCorePEVOperator::apply(const Eigen::MatrixXd& x) const {
    return inverse_->apply(x);
}

void DenseAPYCorePEVOperator::applyCoefficient(
    const Eigen::VectorXd& x, Eigen::VectorXd& y) const {
    if (x.size() != rows()) {
        throw std::runtime_error("DenseAPYCorePEVOperator coefficient size mismatch");
    }
    y.noalias() = coefficient_ * x;
}

MatrixFreeAPYCorePEVOperator::MatrixFreeAPYCorePEVOperator(
    GccSolverPtr gcc_solver, GcnBlockIteratorFactory gcn_factory,
    const Eigen::VectorXd& mnn, const Eigen::VectorXd& dc,
    const Eigen::VectorXd& dn, double lambda,
    const MatrixFreePEVOperator::Options& options,
    const CrossGenomicOperator* cross_operator)
    : gcc_solver_(std::move(gcc_solver)),
      gcn_factory_(std::move(gcn_factory)), dc_(dc), lambda_(lambda),
      n_noncore_(static_cast<int>(dn.size())),
      cross_operator_(cross_operator) {
    validateDc(dc_);
    if (!gcc_solver_ || !gcc_solver_->isValid() || gcc_solver_->rows() != dc_.size()) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator invalid Gcc solver");
    }
    schur_weights_ = makeSchurWeights(mnn, dn, lambda_);
    const bool use_cross = cross_operator_ != nullptr && n_noncore_ > 0 &&
                           cross_operator_->coreRows() == rows() &&
                           cross_operator_->noncoreRows() == n_noncore_;
    if (!use_cross && n_noncore_ > 0 && !gcn_factory_) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator missing Gcn stream");
    }
    if (!use_cross) {
        validateStream();
    }

    diagonal_ = dc_ + lambda_ * gcc_solver_->diagonal();
    if (n_noncore_ > 0) {
        if (use_cross) {
            // Compute diagonal using cross operator: diag(Gcn * diag(s) * Gcn')
            // For each core i: sum_j Gcn_ij^2 * s_j
            for (int i = 0; i < rows(); ++i) {
                Eigen::VectorXd basis = Eigen::VectorXd::Zero(rows());
                basis(i) = 1.0;
                Eigen::VectorXd u = cross_operator_->applyNoncoreFromCore(basis);
                diagonal_(i) += (u.array().square() * schur_weights_.array()).sum();
            }
        } else {
            auto iterator = gcn_factory_();
            if (!iterator) {
                throw std::runtime_error(
                    "MatrixFreeAPYCorePEVOperator Gcn factory returned null during diagonal");
            }
            int expected_start = 0;
            GcnBlock block;
            while (iterator->next(block)) {
                const int cols = block.end_col - block.start_col;
                if (block.start_col != expected_start || cols <= 0 ||
                    block.end_col > n_noncore_ || block.block.rows() != rows() ||
                    block.block.cols() != cols ||
                    static_cast<int>(block.noncore_indices.size()) != cols ||
                    !block.block.allFinite()) {
                    throw std::runtime_error(
                        "MatrixFreeAPYCorePEVOperator invalid Gcn block during diagonal");
                }
                const Eigen::MatrixXd q = gcc_solver_->solve(block.block);
                diagonal_.noalias() += q.array().square().matrix() *
                                       schur_weights_.segment(block.start_col,
                                                              cols);
                expected_start = block.end_col;
            }
            if (expected_start != n_noncore_) {
                throw std::runtime_error(
                    "MatrixFreeAPYCorePEVOperator incomplete Gcn stream during diagonal");
            }
        }
    }
    if (!diagonal_.allFinite() || (diagonal_.array() <= 0.0).any()) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator invalid Jacobi diagonal");
    }

    inverse_ = std::make_unique<MatrixFreePEVOperator>(
        rows(),
        [this](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            applyCoefficient(x, y);
        },
        [this](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            y = x.cwiseQuotient(diagonal_);
        },
        options);
}

void MatrixFreeAPYCorePEVOperator::validateStream() const {
    if (n_noncore_ == 0) return;
    auto iterator = gcn_factory_();
    if (!iterator) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator Gcn factory returned null");
    }
    int expected_start = 0;
    GcnBlock block;
    while (iterator->next(block)) {
        const int cols = block.end_col - block.start_col;
        if (block.start_col != expected_start || cols <= 0 ||
            block.end_col > n_noncore_ || block.block.rows() != rows() ||
            block.block.cols() != cols ||
            static_cast<int>(block.noncore_indices.size()) != cols ||
            !block.block.allFinite()) {
            throw std::runtime_error("MatrixFreeAPYCorePEVOperator invalid Gcn block sequence");
        }
        expected_start = block.end_col;
    }
    if (expected_start != n_noncore_) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator incomplete Gcn stream");
    }
}

void MatrixFreeAPYCorePEVOperator::applyCoefficient(
    const Eigen::VectorXd& x, Eigen::VectorXd& y) const {
    if (x.size() != rows() || !x.allFinite()) {
        throw std::runtime_error("MatrixFreeAPYCorePEVOperator coefficient input mismatch");
    }
    const Eigen::VectorXd a = gcc_solver_->solve(x);
    const bool use_cross =
        cross_operator_ != nullptr && n_noncore_ > 0 &&
        cross_operator_->coreRows() == rows() &&
        cross_operator_->noncoreRows() == n_noncore_;
    if (use_cross) {
        // Phase E3-D2D:  streamed term = Gcn diag(s) Gcn' a, computed as two
        // O(M*(K+Nn)) matrix-free products instead of regenerating the full
        // K x Nn Gcn block per PCG iteration:
        //   u = Gcn' a ; w = s o u ; streamed = Gcn w.
        const Eigen::VectorXd u = cross_operator_->applyNoncoreFromCore(a);
        const Eigen::VectorXd w = schur_weights_.cwiseProduct(u);
        const Eigen::VectorXd streamed = cross_operator_->applyCoreFromNoncore(w);
        y = dc_.cwiseProduct(x) + lambda_ * a + gcc_solver_->solve(streamed);
        return;
    }
    Eigen::VectorXd streamed = Eigen::VectorXd::Zero(rows());
    if (n_noncore_ > 0) {
        auto iterator = gcn_factory_();
        if (!iterator) {
            throw std::runtime_error("MatrixFreeAPYCorePEVOperator Gcn factory returned null during apply");
        }
        int expected_start = 0;
        GcnBlock block;
        while (iterator->next(block)) {
            const int cols = block.end_col - block.start_col;
            if (block.start_col != expected_start || cols <= 0 ||
                block.end_col > n_noncore_ || block.block.rows() != rows() ||
                block.block.cols() != cols ||
                static_cast<int>(block.noncore_indices.size()) != cols ||
                !block.block.allFinite()) {
                throw std::runtime_error("MatrixFreeAPYCorePEVOperator invalid Gcn block during apply");
            }
            const Eigen::VectorXd projected = block.block.transpose() * a;
            streamed.noalias() += block.block *
                schur_weights_.segment(block.start_col, cols).cwiseProduct(projected);
            expected_start = block.end_col;
        }
        if (expected_start != n_noncore_) {
            throw std::runtime_error("MatrixFreeAPYCorePEVOperator incomplete Gcn stream during apply");
        }
    }
    y = dc_.cwiseProduct(x) + lambda_ * a + gcc_solver_->solve(streamed);
}

Eigen::VectorXd MatrixFreeAPYCorePEVOperator::apply(const Eigen::VectorXd& x) const {
    return inverse_->apply(x);
}

Eigen::MatrixXd MatrixFreeAPYCorePEVOperator::apply(const Eigen::MatrixXd& x) const {
    return inverse_->apply(x);
}

const PEVSolveInfo& MatrixFreeAPYCorePEVOperator::lastSolveInfo() const {
    return inverse_->lastSolveInfo();
}

const std::vector<PEVSolveInfo>&
MatrixFreeAPYCorePEVOperator::lastBatchSolveInfo() const {
    return inverse_->lastBatchSolveInfo();
}

} // namespace cosmic
