#include "gcc_solver_pcg.h"
#include "gcc_preconditioner.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>
#include <memory>

namespace cosmic {

// ============================================================
// GccOperator implementation
// ============================================================

GccOperator::GccOperator(PlinkReader& reader,
                         const std::vector<int>& core_indices,
                         double S_2pq,
                         double tune_a,
                         double tune_b,
                         double blend_alpha,
                         const Eigen::MatrixXd* blend_target,
                         int n_threads)
    : reader_(reader)
    , core_indices_(core_indices)
    , S_2pq_(S_2pq)
    , tune_a_(tune_a)
    , tune_b_(tune_b)
    , blend_alpha_(blend_alpha)
    , blend_target_(blend_target)
    , n_threads_(n_threads)
    , n_core_(static_cast<int>(core_indices.size()))
{
}

void GccOperator::apply(const Eigen::VectorXd& v, Eigen::VectorXd& y) const {
    if (v.size() != n_core_) {
        throw std::runtime_error("GccOperator::apply size mismatch");
    }
    if (y.size() != n_core_) y.resize(n_core_);

    // Step 1: tmp = Z_c' * v  (M-vector)
    // Z_c is the standardized genotype matrix for core animals
    Eigen::VectorXd zt_std;
    reader_.multiply_Zt_v(v, zt_std, core_indices_, n_threads_);

    // Undo standardization: tmp = Z_centered' * v
    const auto& inv_stds = reader_.getInvStds();
    const int m = reader_.getNumSnps();
    Eigen::VectorXd tmp(m);
    for (int j = 0; j < m; ++j) {
        double std_j = (inv_stds[j] > 0.0) ? (1.0 / inv_stds[j]) : 1.0;
        tmp(j) = std_j * zt_std(j);
    }

    // Step 2: y = Z_c * tmp / S  (raw Gcc * v)
    Eigen::VectorXd weighted(m);
    for (int j = 0; j < m; ++j) {
        double std_j = (inv_stds[j] > 0.0) ? (1.0 / inv_stds[j]) : 1.0;
        weighted(j) = std_j * tmp(j);
    }
    reader_.multiply_Z_v(weighted, y, core_indices_, n_threads_);

    if (S_2pq_ > 0.0) {
        y /= S_2pq_;
    }

    // Step 3: tuneG  (G_tuned = b*G_raw + a*J, J*v = sum(v)*ones)
    if (tune_b_ != 1.0 || tune_a_ != 0.0) {
        const double sum_v = v.sum();
        y = tune_b_ * y;
        y.array() += tune_a_ * sum_v;
    }

    // Step 4: blend  (nullptr target == identity blend; else dense target T)
    if (blend_alpha_ > 0.0) {
        if (blend_target_ == nullptr) {
            y = (1.0 - blend_alpha_) * y + blend_alpha_ * v;
        } else {
            Eigen::VectorXd blended = (*blend_target_) * v;
            y = (1.0 - blend_alpha_) * y + blend_alpha_ * blended;
        }
    }
}

// ============================================================
// GccSolverPCG implementation
// ============================================================

void GccSolverPCG::compute(const GccOperator& op, const GccSolverPCGOptions& opts) {
    op_ = std::make_shared<const GccOperator>(op);
    opts_ = opts;
    preconditioner_.reset(); // Will be computed lazily on first solve
}

void GccSolverPCG::computePreconditioner() const {
    if (!op_) return;
    if (preconditioner_) return; // Already computed

    const int n = op_->rows();
    if (n <= 0) return;

    if (opts_.preconditioner == "lowrank") {
        // Low-rank spectral preconditioner
        auto jacobi = std::make_shared<JacobiPreconditioner>();
         jacobi->initialize(*op_, 10, 12345);

        auto lr = std::make_shared<LowRankSpectralPreconditioner>();
        LowRankSpectralPreconditioner::Options lr_opts;
        lr_opts.rank = opts_.preconditioner_rank;
        lr_opts.max_lanczos_iter = 80;
        lr_opts.tol = 1e-6;
        lr_opts.seed = 42;
         lr_opts.verbose = opts_.verbose;

         lr->initialize(*op_, *jacobi);

        if (lr->isValid()) {
            preconditioner_ = lr;
            if (opts_.verbose) {
                std::cout << "PCG: Using " << lr->name() << "\n";
            }
        } else {
            // Fallback to Jacobi
            auto jacobi_pc = std::make_shared<JacobiPreconditioner>();
            jacobi_pc->initialize(*op_, 10, 12345);
            preconditioner_ = jacobi_pc;
            if (opts_.verbose) {
                std::cout << "PCG: Low-rank failed, falling back to Jacobi\n";
            }
        }
    } else {
        // Jacobi (diagonal) preconditioner
        auto jacobi_pc = std::make_shared<JacobiPreconditioner>();
        jacobi_pc->initialize(*op_, 10, 12345);
        preconditioner_ = jacobi_pc;
        if (opts_.verbose) {
            std::cout << "PCG: Using Jacobi preconditioner\n";
        }
    }
}

Eigen::VectorXd GccSolverPCG::pcgSolve(const Eigen::VectorXd& rhs) const {
    if (!op_) {
        throw std::runtime_error("GccSolverPCG: not initialized");
    }
    if (rhs.size() != op_->rows()) {
        throw std::runtime_error("GccSolverPCG: rhs size mismatch");
    }

    const int n = op_->rows();
    computePreconditioner();

    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
     Eigen::VectorXd r = rhs;

    Eigen::VectorXd z(n);
    preconditioner_->apply(r, z);
    Eigen::VectorXd p = z;
    double rz_old = r.dot(z);

    int max_iter = std::min(opts_.max_iter, n * 2);
    double tol = opts_.tolerance * rhs.norm();
    if (tol < opts_.tolerance) tol = opts_.tolerance; // Absolute floor

    int iter = 0;
    double initial_residual = r.norm();

    for (iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd Ap(n);
        op_->apply(p, Ap);
        double pAp = p.dot(Ap);
        if (pAp <= 0.0) {
            if (opts_.verbose) {
                std::cout << "PCG: pAp <= 0 at iter " << iter << ", breaking\n";
            }
            break;
        }
        double alpha = rz_old / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        double residual_norm = r.norm();
        if (residual_norm <= tol) {
            if (opts_.verbose) {
                std::cout << "PCG: Converged at iter " << iter
                          << " (residual=" << residual_norm << ", tol=" << tol << ")\n";
            }
            break;
        }

        preconditioner_->apply(r, z);
        double rz_new = r.dot(z);
        double beta = rz_new / rz_old;
        p = z + beta * p;
        rz_old = rz_new;
    }

    if (opts_.verbose) {
        std::cout << "PCG: " << iter + 1 << " iterations, "
                  << "initial_res=" << initial_residual
                  << ", final_res=" << r.norm() << "\n";
    }

    return x;
}

Eigen::VectorXd GccSolverPCG::solve(const Eigen::VectorXd& rhs) const {
    return pcgSolve(rhs);
}

Eigen::MatrixXd GccSolverPCG::solve(const Eigen::MatrixXd& rhs) const {
    if (!op_) {
        throw std::runtime_error("GccSolverPCG: not initialized");
    }
    if (rhs.rows() != op_->rows()) {
        throw std::runtime_error("GccSolverPCG: rhs rows mismatch");
     }

     int n_rhs = rhs.cols();
     Eigen::MatrixXd result(op_->rows(), n_rhs);

    if (opts_.verbose) {
        std::cout << "PCG: Solving for " << n_rhs << " RHS vectors\n";
     }

     for (int j = 0; j < n_rhs; ++j) {
        if (opts_.verbose && n_rhs > 1) {
            std::cout << "  RHS " << j+1 << "/" << n_rhs << "\n";
        }
        result.col(j) = pcgSolve(rhs.col(j));
    }
    return result;
}

Eigen::VectorXd GccSolverPCG::diagonal() const {
    if (!op_) {
        return Eigen::VectorXd();
    }
    // Approximate diagonal of Gcc^-1 using probe vectors
    // diag(A^-1) = E[ x ⊙ (A^-1 x) ] for random x with E[x_i^2] = 1
    const int n = op_->rows();
    const int n_probes = 20;
    Eigen::VectorXd diag_sum = Eigen::VectorXd::Zero(n);
    std::mt19937 rng(54321);
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int p = 0; p < n_probes; ++p) {
        Eigen::VectorXd x(n);
        for (int i = 0; i < n; ++i) x(i) = dist(rng);
        Eigen::VectorXd y = pcgSolve(x);
        diag_sum += x.cwiseProduct(y);
    }
    return diag_sum / n_probes;
}

} // namespace cosmic
