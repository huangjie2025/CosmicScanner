#include "apy_inverse_operator.h"

#include <algorithm>
#include <stdexcept>

#include "cross_genomic_operator.h"

namespace cosmic {

namespace {

class DenseGcnBlockIterator final : public IGcnBlockIterator {
public:
    DenseGcnBlockIterator(const Eigen::MatrixXd& g, std::vector<int> core_indices,
                          std::vector<int> noncore_indices, int block_size)
        : g_(g), core_indices_(std::move(core_indices)),
          noncore_indices_(std::move(noncore_indices)), block_size_(block_size) {}

    bool next(GcnBlock& out) override {
        if (!hasNext()) return false;
        const int end = std::min(start_ + block_size_, static_cast<int>(noncore_indices_.size()));
        out.start_col = start_;
        out.end_col = end;
        out.block.resize(core_indices_.size(), end - start_);
        out.noncore_indices.assign(noncore_indices_.begin() + start_, noncore_indices_.begin() + end);
        for (int i = 0; i < static_cast<int>(core_indices_.size()); ++i) {
            for (int j = start_; j < end; ++j) {
                out.block(i, j - start_) = g_(core_indices_[i], noncore_indices_[j]);
            }
        }
        start_ = end;
        return true;
    }
    bool hasNext() const override { return start_ < static_cast<int>(noncore_indices_.size()); }
    int totalBlocks() const override {
        return (static_cast<int>(noncore_indices_.size()) + block_size_ - 1) / block_size_;
    }
    int blockSize() const override { return block_size_; }

private:
    const Eigen::MatrixXd& g_;
    std::vector<int> core_indices_;
    std::vector<int> noncore_indices_;
    int block_size_;
    int start_ = 0;
};

} // namespace

GcnBlockIteratorFactory makeDenseGcnBlockIteratorFactory(
    const Eigen::MatrixXd& g, const std::vector<int>& core_indices,
    const std::vector<int>& noncore_indices, int block_size) {
    if (block_size <= 0) throw std::runtime_error("Gcn block size must be positive");
    return [&g, core_indices, noncore_indices, block_size]() {
        return std::make_unique<DenseGcnBlockIterator>(
            g, core_indices, noncore_indices, block_size);
    };
}

Eigen::VectorXd DenseApyInverseOperator::apply(const Eigen::VectorXd& x) const {
    return inverse_.applyInverse(x);
}

Eigen::MatrixXd DenseApyInverseOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != rows()) {
        throw std::runtime_error("DenseApyInverseOperator::apply row mismatch");
    }
    Eigen::MatrixXd y(rows(), x.cols());
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
        y.col(j) = inverse_.applyInverse(x.col(j));
    }
    return y;
}

Eigen::VectorXd DenseApyInverseOperator::applyCoreInverse(const Eigen::VectorXd& x) const {
    return inverse_.applyCoreInverse(x);
}

Eigen::MatrixXd DenseApyInverseOperator::applyCoreInverse(const Eigen::MatrixXd& x) const {
    if (x.rows() != inverse_.nCore()) {
        throw std::runtime_error("DenseApyInverseOperator::applyCoreInverse row mismatch");
    }
    return inverse_.gccSolver().solve(x);
}

const Eigen::VectorXd& DenseApyInverseOperator::diagonal() const {
    if (diagonal_.size() != rows()) {
        diagonal_ = Eigen::VectorXd::Zero(rows());
        inverse_.addDiagonalTo(diagonal_, 1.0);
    }
    return diagonal_;
}

MatrixFreeApyInverseOperator::MatrixFreeApyInverseOperator(
    GccSolverPtr gcc_solver,
    GcnBlockIteratorFactory gcn_factory,
    const Eigen::VectorXd& gnn_diag,
    int ng,
    const std::vector<int>& core_indices,
    const std::string& tag,
    const CrossGenomicOperator* cross_operator)
    : gcc_solver_(std::move(gcc_solver)),
      gcn_factory_(std::move(gcn_factory)),
      core_indices_(core_indices),
      tag_(tag),
      cross_operator_(cross_operator) {
    if (!gcc_solver_ || !gcc_solver_->isValid()) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: invalid Gcc solver (" + tag_ + ")");
    }
    std::sort(core_indices_.begin(), core_indices_.end());
    if (core_indices_.empty() || core_indices_.front() < 0 || core_indices_.back() >= ng ||
        std::adjacent_find(core_indices_.begin(), core_indices_.end()) != core_indices_.end()) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: invalid core partition (" + tag_ + ")");
    }
    n_core_ = static_cast<int>(core_indices_.size());
    n_noncore_ = ng - n_core_;
    if (gcc_solver_->rows() != n_core_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: Gcc solver size mismatch (" + tag_ + ")");
    }
    noncore_indices_.reserve(n_noncore_);
    size_t c = 0;
    for (int i = 0; i < ng; ++i) {
        if (c < core_indices_.size() && core_indices_[c] == i) {
            ++c;
        } else {
            noncore_indices_.push_back(i);
        }
    }
    validateAndBuildMnn(gnn_diag);
}

void MatrixFreeApyInverseOperator::validateAndBuildMnn(const Eigen::VectorXd& gnn_diag) {
    if (gnn_diag.size() != n_noncore_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: Gnn diagonal size mismatch (" + tag_ + ")");
    }
    mnn_ = gnn_diag;
    if (n_noncore_ == 0) {
        mnn_inv_.resize(0);
        diagonal_ = gcc_solver_->diagonal();
        return;
    }
    if (!gcn_factory_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: missing Gcn iterator factory (" + tag_ + ")");
    }

    auto iterator = gcn_factory_();
    if (!iterator) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: Gcn iterator factory returned null (" + tag_ + ")");
    }
    int expected_start = 0;
    GcnBlock block;
    while (iterator->next(block)) {
        const int cols = block.end_col - block.start_col;
        if (block.start_col != expected_start || cols <= 0 || block.end_col > n_noncore_ ||
            block.block.rows() != n_core_ || block.block.cols() != cols ||
            static_cast<int>(block.noncore_indices.size()) != cols) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: invalid Gcn block sequence (" + tag_ + ")");
        }
        for (int j = 0; j < cols; ++j) {
            if (block.noncore_indices[j] != noncore_indices_[block.start_col + j]) {
                throw std::runtime_error("MatrixFreeApyInverseOperator: Gcn mapping mismatch (" + tag_ + ")");
            }
        }
        Eigen::MatrixXd q_block = gcc_solver_->solve(block.block);
        for (int j = 0; j < cols; ++j) {
            mnn_(block.start_col + j) -= block.block.col(j).dot(q_block.col(j));
        }
        expected_start = block.end_col;
    }
    if (expected_start != n_noncore_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: incomplete Gcn stream (" + tag_ + ")");
    }
    if ((mnn_.array() <= 0.0).any()) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: Mnn contains non-positive values (" + tag_ + ")");
    }
    mnn_inv_ = mnn_.cwiseInverse();

    diagonal_ = Eigen::VectorXd::Zero(rows());
    const Eigen::VectorXd gcc_inv_diag = gcc_solver_->diagonal();
    for (int i = 0; i < n_core_; ++i) diagonal_(core_indices_[i]) = gcc_inv_diag(i);
    for (int j = 0; j < n_noncore_; ++j) {
        diagonal_(noncore_indices_[j]) = mnn_inv_(j);
    }

    iterator = gcn_factory_();
    expected_start = 0;
    while (iterator && iterator->next(block)) {
        const int cols = block.end_col - block.start_col;
        if (block.start_col != expected_start || block.block.rows() != n_core_ ||
            block.block.cols() != cols) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: invalid Gcn block during diagonal (" + tag_ + ")");
        }
        const Eigen::MatrixXd q_block = gcc_solver_->solve(block.block);
        const Eigen::VectorXd contribution =
            q_block.array().square().matrix() * mnn_inv_.segment(block.start_col, cols);
        for (int i = 0; i < n_core_; ++i) diagonal_(core_indices_[i]) += contribution(i);
        expected_start = block.end_col;
    }
    if (expected_start != n_noncore_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator: incomplete Gcn stream during diagonal (" + tag_ + ")");
    }
}

Eigen::VectorXd MatrixFreeApyInverseOperator::apply(const Eigen::VectorXd& x) const {
    if (x.size() != rows()) {
        throw std::runtime_error("MatrixFreeApyInverseOperator::apply size mismatch (" + tag_ + ")");
    }
    Eigen::VectorXd x_core(n_core_);
    for (int i = 0; i < n_core_; ++i) x_core(i) = x(core_indices_[i]);
    Eigen::VectorXd y_core = gcc_solver_->solve(x_core);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(rows());

    // Phase E3-D2D: replace the per-apply full K x Nn Gcn regeneration with the
    // matrix-free cross operator (O(M*(K+Nn)) per apply).  The exact identity
    //   y = D * x + Q diag(mnn^-1) (x_n - Q' x_c)   (D = block-diag(Gcc^-1, 0))
    // with Q_b = Gcc^-1 Gcn_b and t = mnn^-1 o (x_n - Gcn' Gcc^-1 x_c) satisfies
    //   u = Gcn' a, t = mnn^-1 o (x_n - u), v = Gcn t,
    //   y_core = a - Gcc^-1 v ;  y_n = t   (a = Gcc^-1 x_core).
    const bool use_cross =
        cross_operator_ != nullptr && n_noncore_ > 0 &&
        cross_operator_->coreRows() == n_core_ &&
        cross_operator_->noncoreRows() == n_noncore_;
    if (use_cross) {
        const Eigen::VectorXd u = cross_operator_->applyNoncoreFromCore(y_core);
        Eigen::VectorXd x_noncore(n_noncore_);
        for (int j = 0; j < n_noncore_; ++j) x_noncore(j) = x(noncore_indices_[j]);
        const Eigen::VectorXd t = mnn_inv_.cwiseProduct(x_noncore - u);
        const Eigen::VectorXd v = cross_operator_->applyCoreFromNoncore(t);
        const Eigen::VectorXd y_core2 = y_core - gcc_solver_->solve(v);
        for (int i = 0; i < n_core_; ++i) y(core_indices_[i]) = y_core2(i);
        for (int j = 0; j < n_noncore_; ++j) y(noncore_indices_[j]) = t(j);
        return y;
    }

    if (n_noncore_ > 0) {
        auto iterator = gcn_factory_();
        if (!iterator) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: null Gcn iterator during apply (" + tag_ + ")");
        }
        GcnBlock block;
        int expected_start = 0;
        while (iterator->next(block)) {
            const int cols = block.end_col - block.start_col;
            if (block.start_col != expected_start || block.block.rows() != n_core_ ||
                block.block.cols() != cols) {
                throw std::runtime_error("MatrixFreeApyInverseOperator: invalid Gcn block during apply (" + tag_ + ")");
            }
            Eigen::MatrixXd q_block = gcc_solver_->solve(block.block);
            Eigen::VectorXd x_noncore(cols);
            for (int j = 0; j < cols; ++j) {
                x_noncore(j) = x(noncore_indices_[block.start_col + j]);
            }
            Eigen::VectorXd t = mnn_inv_.segment(block.start_col, cols).cwiseProduct(
                x_noncore - q_block.transpose() * x_core);
            y_core.noalias() -= q_block * t;
            for (int j = 0; j < cols; ++j) y(noncore_indices_[block.start_col + j]) = t(j);
            expected_start = block.end_col;
        }
        if (expected_start != n_noncore_) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: incomplete Gcn stream during apply (" + tag_ + ")");
        }
    }
    for (int i = 0; i < n_core_; ++i) y(core_indices_[i]) = y_core(i);
    return y;
}

Eigen::MatrixXd MatrixFreeApyInverseOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != rows()) {
        throw std::runtime_error("MatrixFreeApyInverseOperator::apply row mismatch (" + tag_ + ")");
    }
    // Phase E3-D2D: matrix RHS via the vector path when the cross operator is
    // active (parallelism then comes from the inner reader primitives).
    const bool use_cross =
        cross_operator_ != nullptr && n_noncore_ > 0 &&
        cross_operator_->coreRows() == n_core_ &&
        cross_operator_->noncoreRows() == n_noncore_;
    if (use_cross) {
        Eigen::MatrixXd y(rows(), x.cols());
        for (Eigen::Index j = 0; j < x.cols(); ++j) {
            y.col(j) = apply(Eigen::VectorXd(x.col(j)));
        }
        return y;
    }
    Eigen::MatrixXd x_core(n_core_, x.cols());
    for (int i = 0; i < n_core_; ++i) x_core.row(i) = x.row(core_indices_[i]);
    Eigen::MatrixXd y_core = gcc_solver_->solve(x_core);
    Eigen::MatrixXd y = Eigen::MatrixXd::Zero(rows(), x.cols());

    if (n_noncore_ > 0) {
        auto iterator = gcn_factory_();
        if (!iterator) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: null Gcn iterator during apply (" + tag_ + ")");
        }
        GcnBlock block;
        int expected_start = 0;
        while (iterator->next(block)) {
            const int cols = block.end_col - block.start_col;
            if (block.start_col != expected_start || block.block.rows() != n_core_ ||
                block.block.cols() != cols) {
                throw std::runtime_error("MatrixFreeApyInverseOperator: invalid Gcn block during apply (" + tag_ + ")");
            }
            Eigen::MatrixXd q_block = gcc_solver_->solve(block.block);
            Eigen::MatrixXd x_noncore(cols, x.cols());
            for (int j = 0; j < cols; ++j) {
                x_noncore.row(j) = x.row(noncore_indices_[block.start_col + j]);
            }
            Eigen::MatrixXd t = x_noncore - q_block.transpose() * x_core;
            t.array().colwise() *= mnn_inv_.segment(block.start_col, cols).array();
            y_core.noalias() -= q_block * t;
            for (int j = 0; j < cols; ++j) {
                y.row(noncore_indices_[block.start_col + j]) = t.row(j);
            }
            expected_start = block.end_col;
        }
        if (expected_start != n_noncore_) {
            throw std::runtime_error("MatrixFreeApyInverseOperator: incomplete Gcn stream during apply (" + tag_ + ")");
        }
    }
    for (int i = 0; i < n_core_; ++i) y.row(core_indices_[i]) = y_core.row(i);
    return y;
}

Eigen::VectorXd MatrixFreeApyInverseOperator::applyCoreInverse(const Eigen::VectorXd& x) const {
    if (x.size() != n_core_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator::applyCoreInverse size mismatch (" + tag_ + ")");
    }
    return gcc_solver_->solve(x);
}

Eigen::MatrixXd MatrixFreeApyInverseOperator::applyCoreInverse(const Eigen::MatrixXd& x) const {
    if (x.rows() != n_core_) {
        throw std::runtime_error("MatrixFreeApyInverseOperator::applyCoreInverse row mismatch (" + tag_ + ")");
    }
    return gcc_solver_->solve(x);
}

std::string MatrixFreeApyInverseOperator::backendName() const {
    return "matrix-free-apy/" + gcc_solver_->backendName();
}

size_t MatrixFreeApyInverseOperator::persistentBytes() const {
    return static_cast<size_t>(mnn_.size() + mnn_inv_.size() + diagonal_.size()) * sizeof(double) +
           (core_indices_.size() + noncore_indices_.size()) * sizeof(int);
}

} // namespace cosmic
