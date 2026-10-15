#include "a22_inverse_operator.h"

#include <algorithm>
#include <stdexcept>

namespace cosmic {

namespace {
using Clock = std::chrono::steady_clock;
double secondsSince(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
}

A22InverseOperator::A22InverseOperator(
    const Eigen::SparseMatrix<double>& a_inverse,
    const std::vector<int>& genotyped_indices)
    : n_pedigree_(static_cast<int>(a_inverse.rows())),
      n_genotyped_(static_cast<int>(genotyped_indices.size())),
      n_untyped_(n_pedigree_ - n_genotyped_) {
    if (a_inverse.rows() != a_inverse.cols() || n_pedigree_ <= 0) {
        throw std::runtime_error("A22InverseOperator requires a nonempty square A-inverse.");
    }
    if (n_genotyped_ <= 0 || n_untyped_ < 0) {
        throw std::runtime_error("A22InverseOperator requires 1..N unique genotyped indices.");
    }

    std::vector<int> full_to_g(n_pedigree_, -1);
    for (int g = 0; g < n_genotyped_; ++g) {
        const int idx = genotyped_indices[g];
        if (idx < 0 || idx >= n_pedigree_) {
            throw std::runtime_error("A22InverseOperator genotyped index is out of range.");
        }
        if (full_to_g[idx] >= 0) {
            throw std::runtime_error("A22InverseOperator genotyped indices must be unique.");
        }
        full_to_g[idx] = g;
    }
    std::vector<int> full_to_u(n_pedigree_, -1);
    int u_next = 0;
    for (int i = 0; i < n_pedigree_; ++i) {
        if (full_to_g[i] < 0) full_to_u[i] = u_next++;
    }

    std::vector<Eigen::Triplet<double>> p_trips, q_trips, r_trips;
    p_trips.reserve(a_inverse.nonZeros());
    q_trips.reserve(a_inverse.nonZeros() / 4 + 1);
    r_trips.reserve(a_inverse.nonZeros());
    for (int outer = 0; outer < a_inverse.outerSize(); ++outer) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(a_inverse, outer); it; ++it) {
            const int gr = full_to_g[it.row()];
            const int gc = full_to_g[it.col()];
            if (gr >= 0 && gc >= 0) {
                r_trips.emplace_back(gr, gc, it.value());
            } else if (gr < 0 && gc < 0) {
                p_trips.emplace_back(full_to_u[it.row()], full_to_u[it.col()], it.value());
            } else if (gr < 0) {
                q_trips.emplace_back(full_to_u[it.row()], gc, it.value());
            }
            // The transposed mixed block is represented by Q' and is not copied.
        }
    }

    p_.resize(n_untyped_, n_untyped_);
    q_.resize(n_untyped_, n_genotyped_);
    r_.resize(n_genotyped_, n_genotyped_);
    p_.setFromTriplets(p_trips.begin(), p_trips.end());
    q_.setFromTriplets(q_trips.begin(), q_trips.end());
    r_.setFromTriplets(r_trips.begin(), r_trips.end());
    p_.makeCompressed();
    q_.makeCompressed();
    r_.makeCompressed();

    r_diagonal_ = r_.diagonal();
    if (n_untyped_ > 0) {
        const auto start = Clock::now();
        p_solver_.compute(p_);
        factor_time_s_ = secondsSince(start);
        if (p_solver_.info() != Eigen::Success) {
            throw std::runtime_error("A22InverseOperator sparse P factorization failed; no dense fallback is allowed.");
        }
        // Eigen stores the unit diagonal of L implicitly; include it in the
        // diagnostic so the reported factor size is not misleading.
        factor_nonzeros_ = n_untyped_ + p_solver_.matrixL().nestedExpression().nonZeros();
    }
}

size_t A22InverseOperator::nonZeros() const {
    return static_cast<size_t>(p_.nonZeros() + q_.nonZeros() + r_.nonZeros());
}

Eigen::VectorXd A22InverseOperator::operator*(const Eigen::VectorXd& x) const {
    if (x.size() != n_genotyped_) {
        throw std::runtime_error("A22InverseOperator vector dimension mismatch.");
    }
    Eigen::VectorXd y = r_ * x;
    if (n_untyped_ == 0) return y;

    const Eigen::VectorXd rhs = q_ * x;
    const auto start = Clock::now();
    const Eigen::VectorXd z = p_solver_.solve(rhs);
    solve_time_s_ += secondsSince(start);
    ++solve_count_;
    if (p_solver_.info() != Eigen::Success || !z.allFinite()) {
        throw std::runtime_error("A22InverseOperator sparse P solve failed.");
    }
    const double rel = (p_ * z - rhs).norm() / (1.0 + rhs.norm());
    max_solve_relative_residual_ = std::max(max_solve_relative_residual_, rel);
    y.noalias() -= q_.transpose() * z;
    return y;
}

Eigen::MatrixXd A22InverseOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != n_genotyped_) {
        throw std::runtime_error("A22InverseOperator matrix-RHS dimension mismatch.");
    }
    Eigen::MatrixXd y = r_ * x;
    if (n_untyped_ == 0) return y;

    const Eigen::MatrixXd rhs = q_ * x;
    const auto start = Clock::now();
    const Eigen::MatrixXd z = p_solver_.solve(rhs);
    solve_time_s_ += secondsSince(start);
    solve_count_ += x.cols();
    if (p_solver_.info() != Eigen::Success || !z.allFinite()) {
        throw std::runtime_error("A22InverseOperator sparse P matrix solve failed.");
    }
    const double rel = (p_ * z - rhs).norm() / (1.0 + rhs.norm());
    max_solve_relative_residual_ = std::max(max_solve_relative_residual_, rel);
    y.noalias() -= q_.transpose() * z;
    return y;
}

void A22InverseOperator::add_diagonal_to(Eigen::VectorXd& diag, double alpha) const {
    if (!hasExactDiagonal()) {
        throw std::runtime_error(
            "A22InverseOperator exact Schur diagonal is unavailable for mixed pedigrees; use the labeled preconditioner diagonal.");
    }
    if (diag.size() != n_genotyped_) {
        throw std::runtime_error("A22InverseOperator diagonal dimension mismatch.");
    }
    diag.noalias() += alpha * r_diagonal_;
}

Eigen::MatrixXd A22InverseOperator::getBlock(int, int) const {
    throw std::runtime_error("A22InverseOperator is operator-only; dense block extraction is disabled.");
}

void A22InverseOperator::visit_triplets(std::function<void(int, int, double)>) const {
    throw std::runtime_error("A22InverseOperator is operator-only; triplet materialization is disabled.");
}

Eigen::MatrixXd A22InverseOperator::toDense() const {
    throw std::runtime_error("A22InverseOperator is operator-only; dense materialization is disabled.");
}

} // namespace cosmic
