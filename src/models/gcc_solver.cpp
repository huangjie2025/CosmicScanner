#include "gcc_solver.h"
#include "gcc_solver_pcg.h"
#include <stdexcept>
#include <iostream>

namespace cosmic {

GccSolverDense::GccSolverDense(const Eigen::MatrixXd& Gcc) {
    compute(Gcc);
}

void GccSolverDense::compute(const Eigen::MatrixXd& Gcc) {
    ldlt_.compute(Gcc);
    if (ldlt_.info() != Eigen::Success) {
        throw std::runtime_error("GccSolverDense: LDLT factorization failed (Gcc not SPD?)");
    }
    computed_ = true;
    const Eigen::VectorXd& D = ldlt_.vectorD();
    if (D.size() > 0) {
        min_diag_d_ = D.minCoeff();
        max_diag_d_ = D.maxCoeff();
    }
}

Eigen::VectorXd GccSolverDense::solve(const Eigen::VectorXd& rhs) const {
    if (!computed_) {
        throw std::runtime_error("GccSolverDense: not computed");
    }
    if (rhs.size() != rows()) {
        throw std::runtime_error("GccSolverDense: solve size mismatch");
    }
    return ldlt_.solve(rhs);
}

Eigen::MatrixXd GccSolverDense::solve(const Eigen::MatrixXd& rhs) const {
    if (!computed_) {
        throw std::runtime_error("GccSolverDense: not computed");
    }
    if (rhs.rows() != rows()) {
        throw std::runtime_error("GccSolverDense: solve rows mismatch");
    }
    return ldlt_.solve(rhs);
}

Eigen::VectorXd GccSolverDense::diagonal() const {
    if (!computed_) {
        throw std::runtime_error("GccSolverDense: not computed");
    }
    int n = rows();
    Eigen::VectorXd diag(n);
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
        e(i) = 1.0;
        Eigen::VectorXd col = solve(e);
        diag(i) = col(i);
    }
    return diag;
}

GccSolverPtr makeDenseSolver() {
    return std::make_shared<GccSolverDense>();
}

GccSolverPtr makeDenseSolver(const Eigen::MatrixXd& Gcc) {
    return std::make_shared<GccSolverDense>(Gcc);
}

GccSolverPCGPtr makePCGSolver() {
    return std::make_shared<GccSolverPCG>();
}

} // namespace cosmic
