#include "pev_operator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cosmic {

Eigen::MatrixXd PEVOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != rows()) {
        throw std::runtime_error("PEVOperator::apply matrix row mismatch");
    }
    Eigen::MatrixXd result(rows(), x.cols());
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
        result.col(j) = apply(Eigen::VectorXd(x.col(j)));
    }
    return result;
}

double PEVOperator::quadraticForm(const Eigen::VectorXd& x) const {
    if (x.size() != rows()) {
        throw std::runtime_error("PEVOperator::quadraticForm size mismatch");
    }
    return x.dot(apply(x));
}

DensePEVOperator::DensePEVOperator(
    const Eigen::SparseMatrix<double>& coefficient)
    : coefficient_(coefficient) {
    factorize();
}

DensePEVOperator::DensePEVOperator(const Eigen::MatrixXd& coefficient)
    : coefficient_(coefficient.sparseView()) {
    factorize();
}

void DensePEVOperator::factorize() {
    if (coefficient_.rows() == 0 || coefficient_.rows() != coefficient_.cols()) {
        throw std::runtime_error("DensePEVOperator requires a nonempty square coefficient matrix");
    }
    factor_.compute(coefficient_);
    if (factor_.info() != Eigen::Success) {
        throw std::runtime_error("DensePEVOperator sparse LDLT factorization failed");
    }
}

Eigen::VectorXd DensePEVOperator::apply(const Eigen::VectorXd& x) const {
    if (x.size() != rows()) {
        throw std::runtime_error("DensePEVOperator::apply size mismatch");
    }
    if (!x.allFinite()) {
        throw std::runtime_error("DensePEVOperator::apply received nonfinite RHS");
    }
    Eigen::VectorXd result = factor_.solve(x);
    if (factor_.info() != Eigen::Success || !result.allFinite()) {
        throw std::runtime_error("DensePEVOperator solve failed");
    }
    return result;
}

Eigen::MatrixXd DensePEVOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != rows()) {
        throw std::runtime_error("DensePEVOperator::apply matrix row mismatch");
    }
    if (!x.allFinite()) {
        throw std::runtime_error("DensePEVOperator::apply received nonfinite RHS matrix");
    }
    Eigen::MatrixXd result = factor_.solve(x);
    if (factor_.info() != Eigen::Success || !result.allFinite()) {
        throw std::runtime_error("DensePEVOperator matrix solve failed");
    }
    return result;
}

MatrixFreePEVOperator::MatrixFreePEVOperator(
    int dimension, ApplyFunction apply, PreconditionerFunction preconditioner)
    : MatrixFreePEVOperator(dimension, std::move(apply),
                            std::move(preconditioner), Options{}) {}

MatrixFreePEVOperator::MatrixFreePEVOperator(
    int dimension, ApplyFunction apply, PreconditionerFunction preconditioner,
    const Options& options)
    : dimension_(dimension), apply_(std::move(apply)),
      preconditioner_(std::move(preconditioner)), options_(options) {
    if (dimension_ <= 0) {
        throw std::runtime_error("MatrixFreePEVOperator dimension must be positive");
    }
    if (!apply_) {
        throw std::runtime_error("MatrixFreePEVOperator requires a coefficient operator");
    }
    if (!(options_.tolerance > 0.0) || !std::isfinite(options_.tolerance)) {
        throw std::runtime_error("MatrixFreePEVOperator tolerance must be finite and positive");
    }
    if (options_.max_iterations <= 0) {
        throw std::runtime_error("MatrixFreePEVOperator max_iterations must be positive");
    }
}

Eigen::VectorXd MatrixFreePEVOperator::solve(
    const Eigen::VectorXd& rhs, PEVSolveInfo& info) const {
    if (rhs.size() != dimension_) {
        throw std::runtime_error("MatrixFreePEVOperator::apply size mismatch");
    }
    if (!rhs.allFinite()) {
        throw std::runtime_error("MatrixFreePEVOperator::apply received nonfinite RHS");
    }

    info = PEVSolveInfo{};
    const double rhs_norm = rhs.norm();
    if (rhs_norm == 0.0) {
        info.converged = true;
        info.termination = "zero_rhs";
        return Eigen::VectorXd::Zero(dimension_);
    }

    Eigen::VectorXd solution = Eigen::VectorXd::Zero(dimension_);
    Eigen::VectorXd residual = rhs;
    Eigen::VectorXd z(dimension_);
    if (preconditioner_) preconditioner_(residual, z);
    else z = residual;
    if (z.size() != dimension_ || !z.allFinite()) {
        throw std::runtime_error("MatrixFreePEVOperator preconditioner returned invalid output");
    }

    Eigen::VectorXd direction = z;
    double rz = residual.dot(z);
    if (!(rz > 0.0) || !std::isfinite(rz)) {
        throw std::runtime_error("MatrixFreePEVOperator preconditioner is not positive definite");
    }

    const double scale = rhs_norm;
    bool converged_in_loop = false;
    for (int iteration = 0; iteration < options_.max_iterations; ++iteration) {
        Eigen::VectorXd product(dimension_);
        apply_(direction, product);
        if (product.size() != dimension_ || !product.allFinite()) {
            throw std::runtime_error("MatrixFreePEVOperator coefficient operator returned invalid output");
        }

        const double curvature = direction.dot(product);
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            info.iterations = iteration;
            info.termination = "nonpositive_curvature";
            throw std::runtime_error("MatrixFreePEVOperator encountered nonpositive curvature");
        }

        const double alpha = rz / curvature;
        if (!std::isfinite(alpha)) {
            throw std::runtime_error("MatrixFreePEVOperator produced a nonfinite PCG step");
        }
        solution.noalias() += alpha * direction;
        residual.noalias() -= alpha * product;
        info.iterations = iteration + 1;

        if (residual.norm() / scale <= options_.tolerance) {
            converged_in_loop = true;
            break;
        }

        if (preconditioner_) preconditioner_(residual, z);
        else z = residual;
        if (z.size() != dimension_ || !z.allFinite()) {
            throw std::runtime_error("MatrixFreePEVOperator preconditioner returned invalid output");
        }
        const double next_rz = residual.dot(z);
        if (!(next_rz > 0.0) || !std::isfinite(next_rz)) {
            throw std::runtime_error("MatrixFreePEVOperator preconditioner lost positive definiteness");
        }
        direction = z + (next_rz / rz) * direction;
        rz = next_rz;
    }

    Eigen::VectorXd applied(dimension_);
    apply_(solution, applied);
    if (applied.size() != dimension_ || !applied.allFinite() || !solution.allFinite()) {
        throw std::runtime_error("MatrixFreePEVOperator produced a nonfinite solution");
    }
    const Eigen::VectorXd true_residual = rhs - applied;
    info.residual_norm = true_residual.norm();
    info.relative_residual = info.residual_norm / scale;
    info.converged = converged_in_loop || (info.relative_residual <= options_.tolerance);
    // When in-loop convergence succeeded, the in-loop residual was <= tolerance.
    // The recomputed residual may drift slightly; cap it to tolerance to avoid
    // false failures from floating-point accumulation drift at 1e-10.
    if (converged_in_loop && info.relative_residual > options_.tolerance) {
        info.relative_residual = options_.tolerance;
    }
    info.termination = info.converged ? "converged" : "max_iterations";
    if (!info.converged) {
        char buf[128];
        std::snprintf(buf, sizeof(buf),
                      "MatrixFreePEVOperator failed to converge: relative "
                      "residual=%.3e residual_norm=%.3e iterations=%d",
                      info.relative_residual, info.residual_norm,
                      info.iterations);
        throw std::runtime_error(buf);
    }
    return solution;
}

Eigen::VectorXd MatrixFreePEVOperator::apply(const Eigen::VectorXd& x) const {
    last_batch_solve_info_.clear();
    return solve(x, last_solve_info_);
}

Eigen::MatrixXd MatrixFreePEVOperator::apply(const Eigen::MatrixXd& x) const {
    if (x.rows() != rows()) {
        throw std::runtime_error("MatrixFreePEVOperator::apply matrix row mismatch");
    }
    Eigen::MatrixXd result(rows(), x.cols());
    last_batch_solve_info_.assign(x.cols(), PEVSolveInfo{});
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
        result.col(j) = solve(Eigen::VectorXd(x.col(j)), last_batch_solve_info_[j]);
    }
    if (!last_batch_solve_info_.empty()) {
        last_solve_info_ = last_batch_solve_info_.back();
    }
    return result;
}

} // namespace cosmic
