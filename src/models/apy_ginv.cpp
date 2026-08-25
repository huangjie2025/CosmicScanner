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

    // Gcc, Gcn in genotyped-core / genotyped-noncore order.
    Eigen::MatrixXd Gcc(n_core_, n_core_);
    Eigen::MatrixXd Gcn(n_core_, n_noncore_);
    for (int i = 0; i < n_core_; ++i) {
        const int ci = core_idx_[i];
        for (int j = 0; j < n_core_; ++j) Gcc(i, j) = G(ci, core_idx_[j]);
        for (int j = 0; j < n_noncore_; ++j) Gcn(i, j) = G(ci, noncore_idx_[j]);
    }

    // Factor Gcc (dense LDLT) and store its explicit inverse (n_core^2; the
    // paper's Gcc^-1 is likewise stored explicitly).
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);
    if (ldlt.info() != Eigen::Success) {
        throw std::runtime_error("ApyGinv: Gcc factorization failed (not SPD?) (" + tag_ + ")");
    }
    {
        const Eigen::MatrixXd& D = ldlt.vectorD();
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
    gcc_inv_.noalias() = ldlt.solve(Eigen::MatrixXd::Identity(n_core_, n_core_));

    // Q = Gcc^-1 Gcn  (n_core x n_noncore).
    q_.noalias() = ldlt.solve(Gcn);

    // Noncore conditional variances  m_jj = Gnn(j,j) - Gnc(j,:) Gcc^-1 Gcn(:,j)
    // = G(ci+...,...) - sum_i G(ci, c_i) Q(i,j).
    Eigen::VectorXd mnn(n_noncore_);
    int n_nonpos = 0;
    for (int j = 0; j < n_noncore_; ++j) {
        const int nj = noncore_idx_[j];
        double m = G(nj, nj);
        for (int i = 0; i < n_core_; ++i) {
            m -= G(nj, core_idx_[i]) * q_(i, j);
        }
        mnn(j) = m;
        if (!(m > 0.0)) ++n_nonpos;
    }
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
    // G_APY^-1; O(n_core^2 * n_noncore) build for the exact APY route).
    core_block_.noalias() = q_ * mnn_inv_.asDiagonal() * q_.transpose();
    core_block_.noalias() += gcc_inv_;

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
        y.noalias() = gcc_inv_ * x;
        return y;
    }
    // t = Mnn^-1 (x_n - Q' x_c);  y_c = Gcc^-1 x_c - Q t;  y_n = t
    Eigen::VectorXd x_c(n_core_), x_n(n_noncore_);
    for (int i = 0; i < n_core_; ++i) x_c(i) = x(core_idx_[i]);
    for (int j = 0; j < n_noncore_; ++j) x_n(j) = x(noncore_idx_[j]);
    Eigen::VectorXd t = mnn_inv_.cwiseProduct(x_n - q_.transpose() * x_c);
    Eigen::VectorXd y_c = gcc_inv_ * x_c - q_ * t;
    for (int i = 0; i < n_core_; ++i) y(core_idx_[i]) = y_c(i);
    for (int j = 0; j < n_noncore_; ++j) y(noncore_idx_[j]) = t(j);
    return y;
}

Eigen::VectorXd ApyGinv::applyCoreInverse(const Eigen::VectorXd& x_core) const {
    if (x_core.size() != n_core_) {
        throw std::runtime_error("ApyGinv::applyCoreInverse size mismatch (" + tag_ + ")");
    }
    Eigen::VectorXd y(n_core_);
    y.noalias() = gcc_inv_ * x_core;
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

} // namespace cosmic