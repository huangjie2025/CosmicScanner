#include "apy_ginv.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cosmic {

void ApyGinv::build(const Eigen::MatrixXd& G, const std::vector<int>& core_idx,
                    const std::string& tag) {
    tag_ = tag;
    const int ng = (int)G.rows();
    if (G.cols() != ng) {
        throw std::runtime_error("ApyGinv: G must be square (" + tag_ + ")");
    }
    initPartition(core_idx, ng);

    // Gcc, Gcn in genotyped-core / genotyped-noncore order; Gnn diagonal.
    Eigen::MatrixXd Gcc(n_core_, n_core_);
    Eigen::MatrixXd Gcn(n_core_, n_noncore_);
    Eigen::VectorXd Gnn_diag(n_noncore_);
    for (int i = 0; i < n_core_; ++i) {
        const int ci = core_idx_[i];
        for (int j = 0; j < n_core_; ++j) Gcc(i, j) = G(ci, core_idx_[j]);
        for (int j = 0; j < n_noncore_; ++j) Gcn(i, j) = G(ci, noncore_idx_[j]);
    }
    for (int j = 0; j < n_noncore_; ++j) Gnn_diag(j) = G(noncore_idx_[j], noncore_idx_[j]);
    buildBlocksImpl(Gcc, Gcn, Gnn_diag);
}

void ApyGinv::build(const Eigen::MatrixXd& Gcc, const Eigen::MatrixXd& Gcn,
                    const Eigen::VectorXd& Gnn_diag, int ng,
                    const std::vector<int>& core_idx, const std::string& tag) {
    tag_ = tag;
    initPartition(core_idx, ng);
    if (Gcc.rows() != n_core_ || Gcc.cols() != n_core_) {
        throw std::runtime_error(
            "ApyGinv: Gcc must be n_core x n_core (" + tag_ + ")");
    }
    if (Gcn.rows() != n_core_ || Gcn.cols() != n_noncore_) {
        throw std::runtime_error(
            "ApyGinv: Gcn must be n_core x n_noncore (" + tag_ + ")");
    }
    if (Gnn_diag.size() != n_noncore_) {
        throw std::runtime_error(
            "ApyGinv: Gnn diagonal must have n_noncore entries (" + tag_ + ")");
    }
    buildBlocksImpl(Gcc, Gcn, Gnn_diag);
}

void ApyGinv::build(const Eigen::MatrixXd& Gcc, const Eigen::VectorXd& Gnn_diag,
                    int ng, const std::vector<int>& core_idx, const std::string& tag,
                    const GcnBlockIteratorFactory& gcn_factory) {
    tag_ = tag;
    initPartition(core_idx, ng);
    if (Gcc.rows() != n_core_ || Gcc.cols() != n_core_) {
        throw std::runtime_error(
            "ApyGinv: Gcc must be n_core x n_core (" + tag_ + ")");
    }
    if (Gnn_diag.size() != n_noncore_) {
        throw std::runtime_error(
            "ApyGinv: Gnn diagonal must have n_noncore entries (" + tag_ + ")");
    }
    factorGcc(Gcc);
    q_.resize(n_core_, n_noncore_);
    Eigen::VectorXd mnn = Gnn_diag;
    int expected_start = 0;
    if (n_noncore_ > 0) {
        if (!gcn_factory) {
            throw std::runtime_error("ApyGinv: missing Gcn iterator factory (" + tag_ + ")");
        }
        auto iterator = gcn_factory();
        if (!iterator) {
            throw std::runtime_error("ApyGinv: null Gcn iterator (" + tag_ + ")");
        }
        GcnBlock block;
        while (iterator->next(block)) {
            const int cols = block.end_col - block.start_col;
            if (block.start_col != expected_start || cols <= 0 || block.end_col > n_noncore_ ||
                block.block.rows() != n_core_ || block.block.cols() != cols ||
                static_cast<int>(block.noncore_indices.size()) != cols) {
                throw std::runtime_error("ApyGinv: invalid Gcn block sequence (" + tag_ + ")");
            }
            for (int j = 0; j < cols; ++j) {
                if (block.noncore_indices[j] != noncore_idx_[block.start_col + j]) {
                    throw std::runtime_error("ApyGinv: Gcn noncore mapping mismatch (" + tag_ + ")");
                }
            }
            Eigen::MatrixXd q_block = gcc_solver_->solve(block.block);
            q_.middleCols(block.start_col, cols) = q_block;
            for (int j = 0; j < cols; ++j) {
                for (int i = 0; i < n_core_; ++i) {
                    mnn(block.start_col + j) -= block.block(i, j) * q_block(i, j);
                }
            }
            expected_start = block.end_col;
        }
    }
    if (expected_start != n_noncore_) {
        throw std::runtime_error("ApyGinv: incomplete Gcn stream (" + tag_ + ")");
    }
    int n_nonpos = 0;
    for (int j = 0; j < n_noncore_; ++j) {
        if (!(mnn(j) > 0.0)) ++n_nonpos;
    }
    finishBlocks(mnn, n_nonpos);
}

// Phase E3-B1: Matrix-free PCG backend. Builds using a GccOperator
// (matrix-free apply) and PCG options instead of dense Gcc.
void ApyGinv::build(const GccOperator& op, int ng,
                    const std::vector<int>& core_idx, const std::string& tag,
                    const GccSolverPCGOptions& pcg_opts) {
    tag_ = tag;
    initPartition(core_idx, ng);

    // Create GccSolverPCG with the matrix-free operator
    auto pcg_solver = makePCGSolver();
    pcg_solver->compute(op, pcg_opts);

    if (!pcg_solver->isValid()) {
        throw std::runtime_error("ApyGinv: PCG solver initialization failed (" + tag_ + ")");
    }

    // Store the PCG solver
    gcc_solver_ = pcg_solver;

    // Now compute Q, Mnn, core_block using the PCG solver
    const int n_core = n_core_;
    const int n_noncore = n_noncore_;

    // We need Gcn and Gnn_diag blocks for the PCG path.
    // Since we only have the operator, we need to extract these blocks.
    // For now, we'll need to compute them by applying the operator to basis vectors.
    // This is a simplified implementation - a full implementation would need
    // to extract Gcn and Gnn_diag from the genotype data.

    // For now, throw an error indicating this path needs full implementation
    throw std::runtime_error("ApyGinv: Matrix-free PCG build not fully implemented yet. "
                             "Need to extract Gcn and Gnn_diag from genotype data. (" + tag_ + ")");
}

void ApyGinv::initPartition(const std::vector<int>& core_idx, int ng) {
    // Validate the partition: unique, ascending, in range.
    core_idx_ = core_idx;
    std::sort(core_idx_.begin(), core_idx_.end());
    core_idx_.erase(std::unique(core_idx_.begin(), core_idx_.end()), core_idx_.end());
    if (core_idx_.empty() || core_idx_.front() < 0 || core_idx_.back() >= ng) {
        throw std::runtime_error("ApyGinv: core indices out of range [0, n_g) (" + tag_ + ")");
    }
    if ((int)core_idx_.size() != (int)core_idx.size()) {
        throw std::runtime_error("ApyGinv: duplicate core indices not allowed (" + tag_ + ")");
    }
    n_core_ = (int)core_idx_.size();
    info_.core_is_all = (n_core_ == ng);
    n_noncore_ = ng - n_core_;
    noncore_idx_.clear();
    noncore_idx_.reserve(n_noncore_);
    {
        int ci = 0;
        for (int i = 0; i < ng; ++i) {
            while (ci < n_core_ && core_idx_[ci] < i) ++ci;
            if (ci < n_core_ && core_idx_[ci] == i) continue;
            noncore_idx_.push_back(i);
        }
    }
    info_.n_core = n_core_;
    info_.n_noncore = n_noncore_;
    info_.core_is_all = (n_core_ == ng);
}

void ApyGinv::buildBlocksImpl(const Eigen::MatrixXd& Gcc,
                               const Eigen::MatrixXd& Gcn,
                               const Eigen::VectorXd& Gnn_diag) {
    factorGcc(Gcc);

    // Q = Gcc^-1 Gcn  (n_core x n_noncore).
    q_.noalias() = gcc_solver_->solve(Gcn);

    // Noncore conditional variances  m_jj = Gnn(j,j) - Gnc(j,:) Gcc^-1 Gcn(:,j)
    // = Gnn_diag(j) - sum_i Gcn(i, j) Q(i, j)  (Gcn core rows x noncore cols).
    Eigen::VectorXd mnn(n_noncore_);
    int n_nonpos = 0;
    for (int j = 0; j < n_noncore_; ++j) {
        double m = Gnn_diag(j);
        for (int i = 0; i < n_core_; ++i) {
            m -= Gcn(i, j) * q_(i, j);
        }
        mnn(j) = m;
        if (!(m > 0.0)) ++n_nonpos;
    }
    finishBlocks(mnn, n_nonpos);
}

void ApyGinv::factorGcc(const Eigen::MatrixXd& Gcc) {
    // Factor Gcc (dense LDLT) and store the solver (n_core^2 LDLT factors).
    // We do NOT store the explicit Gcc^-1; it is computed on-demand via solve().
    gcc_solver_ = makeDenseSolver(Gcc);
    const GccSolverDense& dense_solver = dynamic_cast<const GccSolverDense&>(*gcc_solver_);
    const Eigen::LDLT<Eigen::MatrixXd>& ldlt = dense_solver.ldlt();
    const Eigen::VectorXd& D = ldlt.vectorD();
    if (D.size() > 0) {
        info_.gcc_min_diag_d = D.minCoeff();
        info_.gcc_max_diag_d = D.maxCoeff();
    }
    // Relative near-singularity gate (Eigen's LDLT does not signal
    // semidefinite inputs): a pivot <= 1e-12 of the largest pivot means
    // Gcc is (numerically) singular and the APY inverse is ill-defined.
    if (D.size() > 0) {
        const double scale = std::max(std::fabs(D.maxCoeff()), 1.0);
        if (D.minCoeff() <= 1e-12 * scale) {
            throw std::runtime_error(
                "ApyGinv: Gcc is not positive definite (min/max LDLT pivot " +
                std::to_string(D.minCoeff()) + " / " + std::to_string(D.maxCoeff()) +
                "); the core block is (near-)singular (" + tag_ + ")");
        }
    }
}

void ApyGinv::finishBlocks(const Eigen::VectorXd& mnn, int n_nonpos) {
    info_.n_mnn_nonpositive = n_nonpos;
    if (n_noncore_ > 0) {
        info_.mnn_min = mnn.minCoeff();
        info_.mnn_max = mnn.maxCoeff();
    }
    if (n_nonpos > 0) {
        throw std::runtime_error(
            "ApyGinv: " + std::to_string(n_nonpos) +
            " noncore conditional variance(s) Mnn <= 0 (invalid APY partition; "
            "increase --ssgwas-apy-core-size or seed) (" + tag_ + ")");
    }
    mnn_ = mnn;
    mnn_inv_ = mnn.cwiseInverse();

    // Core block Gcc^-1 + Q Mnn^-1 Q' (dense n_core^2; the exact core block of
    // G_APY^-1). For the dense LDLT backend, we compute the full Gcc^-1 on-demand
    // during build (not stored as member) to maintain exact numerical parity.
    // Future PCG backend will use a different approach.
    const Eigen::MatrixXd gcc_inv_full = gcc_solver_->solve(Eigen::MatrixXd::Identity(n_core_, n_core_).eval());
    core_block_.noalias() = gcc_inv_full + q_ * mnn_inv_.asDiagonal() * q_.transpose();

    // Pattern size of the full G_APY^-1 (both orientations, symmetric).
    size_t nnz = 0;
    nnz += (size_t)n_core_ * n_core_;
    nnz += 2 * (size_t)n_core_ * n_noncore_;
    nnz += (size_t)n_noncore_;
    info_.dense_nnz = nnz;
}

Eigen::VectorXd ApyGinv::applyInverse(const Eigen::VectorXd& x) const {
    if (x.size() != n()) {
        throw std::runtime_error("ApyGinv::applyInverse size mismatch (" + tag_ + ")");
    }
    Eigen::VectorXd y(n());
    if (n_noncore_ == 0) {
        y.noalias() = gcc_solver_->solve(x);
        return y;
    }
    // t = Mnn^-1 (x_n - Q' x_c);  y_c = Gcc^-1 x_c - Q t;  y_n = t
    Eigen::VectorXd x_c(n_core_), x_n(n_noncore_);
    for (int i = 0; i < n_core_; ++i) x_c(i) = x(core_idx_[i]);
    for (int j = 0; j < n_noncore_; ++j) x_n(j) = x(noncore_idx_[j]);
    Eigen::VectorXd t = mnn_inv_.cwiseProduct(x_n - q_.transpose() * x_c);
    Eigen::VectorXd y_c = gcc_solver_->solve(x_c) - q_ * t;
    for (int i = 0; i < n_core_; ++i) y(core_idx_[i]) = y_c(i);
    for (int j = 0; j < n_noncore_; ++j) y(noncore_idx_[j]) = t(j);
    return y;
}

Eigen::VectorXd ApyGinv::applyCoreInverse(const Eigen::VectorXd& x_core) const {
    if (x_core.size() != n_core_) {
        throw std::runtime_error("ApyGinv::applyCoreInverse size mismatch (" + tag_ + ")");
    }
    Eigen::VectorXd y = gcc_solver_->solve(x_core);
    return y;
}

void ApyGinv::visitTriplets(const std::function<void(int, int, double)>& cb) const {
    for (int i = 0; i < n_core_; ++i) {
        const int ri = core_idx_[i];
        for (int j = 0; j < n_core_; ++j) {
            cb(ri, core_idx_[j], core_block_(i, j));
        }
    }
    for (int j = 0; j < n_noncore_; ++j) {
        const int rn = noncore_idx_[j];
        const double dj = mnn_inv_(j);
        for (int i = 0; i < n_core_; ++i) {
            const double v = -dj * q_(i, j);
            cb(core_idx_[i], rn, v);
            cb(rn, core_idx_[i], v);
        }
        cb(rn, rn, dj);
    }
}

void ApyGinv::addDiagonalTo(Eigen::VectorXd& diag, double alpha) const {
    for (int i = 0; i < n_core_; ++i) {
        const int r = core_idx_[i];
        if (r < diag.size()) diag(r) += core_block_(i, i) * alpha;
    }
    for (int j = 0; j < n_noncore_; ++j) {
        const int r = noncore_idx_[j];
        if (r < diag.size()) diag(r) += mnn_inv_(j) * alpha;
    }
}

Eigen::VectorXd ApyGinv::gccInvDiagonal() const {
    return gcc_solver_->diagonal();
}

Eigen::MatrixXd ApyGinv::gccInvFull() const {
    const int n = n_core_;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    return gcc_solver_->solve(I);
}

} // namespace cosmic
