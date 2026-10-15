#include "matrix_free_vce.h"
#include "plink_reader.h"
#include "bgen_reader.h"
#include "pgen_reader.h"
#include "vcf/vcf_reader.h"
#include "kinship_operator.h"
#include "matrix_free_pcg.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <vector>
#include <random>
#include <chrono>

namespace cosmic {

// ============================================================================
// MatrixFreeVce Implementation
// ============================================================================

MatrixFreeVce::MatrixFreeVce() : options_(Options{}), n_samples_(0), n_covars_(0), converged_(false), iterations_(0) {}

MatrixFreeVce::MatrixFreeVce(const Options& opts) : options_(opts), n_samples_(0), n_covars_(0), converged_(false), iterations_(0) {}

void MatrixFreeVce::prepare(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const class KinshipOperator& kinship) {
    prepareMatrixFree(y, X, kinship);
}

void MatrixFreeVce::prepareMatrixFree(const Eigen::VectorXd& y, const Eigen::MatrixXf& X, const class KinshipOperator& kinship) {
    n_samples_ = y.size();
    n_covars_ = X.cols();
    kinship_ptr_ = &kinship;
    y_ = y;
    X_ = X;

    if (X.rows() != n_samples_ || kinship.rows() != n_samples_ || kinship.cols() != n_samples_) {
        throw std::runtime_error("MatrixFreeVce::prepare: Dimension mismatch");
    }

    // Step 1: Perform randomized SVD on the matrix-free kinship operator
    const int n = kinship.rows();
    const int k = std::max(1, std::min(options_.rank, n - 1));
    const int oversample = std::max(0, options_.oversample);
    const int sketch_size = std::min(n, k + oversample);

    std::mt19937_64 generator(static_cast<uint64_t>(options_.seed));
    std::normal_distribution<double> normal(0.0, 1.0);

    // Generate random test matrix Omega (n x sketch_size)
    Eigen::MatrixXd omega(n_samples_, sketch_size);
    for (int j = 0; j < sketch_size; ++j)
        for (int i = 0; i < n_samples_; ++i) omega(i, j) = normal(generator);

    // Step 2: Compute Y = K * Omega using matrix-free operator
    Eigen::MatrixXd sample(n_samples_, sketch_size);
    kinship.apply(omega, sample, options_.n_threads);

    // Power iterations for better accuracy
    for (int iteration = 0; iteration < 2; ++iteration) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(sample);
        Eigen::MatrixXd q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples_, sketch_size);
        kinship.apply(q, sample, options_.n_threads);
    }

    // Final QR to get orthonormal basis Q
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(sample);
    Eigen::MatrixXd q = qr.householderQ() * Eigen::MatrixXd::Identity(n_samples_, sketch_size);

    // Step 3: Compute small projected matrix B = Q^T * K * Q
    Eigen::MatrixXd kq(n_samples_, sketch_size);
    kinship.apply(q, kq, options_.n_threads);
    Eigen::MatrixXd small = q.transpose() * kq;

    // Step 4: Eigendecomposition of small matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(small);
    if (eig.info() != Eigen::Success) throw std::runtime_error("Low-rank GRM eigensolver failed");

    Eigen::VectorXd leading_eigenvalues(k);
    Eigen::MatrixXd leading_u(n_samples_, k);
    for (int j = 0; j < k; ++j) {
        const int source = sketch_size - 1 - j;
        leading_eigenvalues[j] = std::max(0.0, eig.eigenvalues()[source]);
        leading_u.col(j) = q * eig.eigenvectors().col(source);
    }

    const double total_trace = kinship.trace();
    int selected_rank = k;
    if (options_.trace_threshold > 0.0 && total_trace > 0.0) {
        const int first_allowed = std::min(k, std::max(1, options_.min_rank));
        double retained = 0.0;
        for (int j = 0; j < k; ++j) {
            retained += leading_eigenvalues[j];
            if (j + 1 >= first_allowed &&
                retained / total_trace >= options_.trace_threshold) {
                selected_rank = j + 1;
                break;
            }
        }
    }

    D_ = leading_eigenvalues.head(selected_rank);
    U_ = leading_u.leftCols(selected_rank);
    int rank = selected_rank;

    low_rank_residual_eigenvalue_ = std::max(
        0.0, (total_trace - D_.head(selected_rank).sum()) / static_cast<double>(n_samples_ - rank));

    // Project covariates and phenotype onto the low-rank basis
    if (n_covars_ > 0) {
        UtX_ = U_.transpose() * X.cast<double>();
    } else {
        UtX_ = Eigen::MatrixXd::Zero(U_.cols(), 0);
    }
    Uty_ = U_.transpose() * y.cast<double>();

    W_proj_ = U_;

    converged_ = true;
    iterations_ = 1;
}

bool MatrixFreeVce::runNullModel() {
    return runBrentNullModel();
}

bool MatrixFreeVce::runIterativeReml(bool average_information) {
    return runBrentNullModel();
}

bool MatrixFreeVce::runIdulNullModel(bool monotone, bool accelerated) {
    return runBrentNullModel();
}

bool MatrixFreeVce::runBrentNullModel() {
    // Residual-augmented low-rank REML, matching LmmEngine::prepareLowRank:
    //   R = U diag(lambda) U' + r0 (I - U U')
    //   V = sigma2_e * R,   eta = sigma2_g / sigma2_e
    //   R^{-1} = base * I + U (1/(1+eta*lambda) - base) U'
    //   base   = 1 / (1 + eta * r0)
    const int n = n_samples_;
    const int p = n_covars_;
    const int rank = static_cast<int>(D_.size());
    const int df = n - p;
    if (df <= 0) throw std::invalid_argument("MatrixFreeVce: no residual degrees of freedom");
    if (rank <= 0) throw std::runtime_error("MatrixFreeVce: empty low-rank basis");

    const double r0 = low_rank_residual_eigenvalue_;
    const Eigen::MatrixXd Xd = X_.cast<double>();
    const Eigen::VectorXd yd = y_.cast<double>();

    auto apply_ratio_inverse = [&](double eta, const Eigen::MatrixXd& rhs) {
        const double base = 1.0 / (1.0 + eta * r0);
        const Eigen::VectorXd correction =
            (1.0 / (1.0 + eta * D_.array()) - base).matrix();
        return (base * rhs + U_ *
                (correction.asDiagonal() * (U_.transpose() * rhs))).eval();
    };

    auto objective = [&](double log_eta, double* quadratic) {
        const double eta = std::exp(log_eta);
        Eigen::MatrixXd rinv_x = apply_ratio_inverse(eta, Xd);
        Eigen::MatrixXd information = Xd.transpose() * rinv_x;
        Eigen::LDLT<Eigen::MatrixXd> solver(information);
        if (solver.info() != Eigen::Success || !solver.isPositive()) {
            return std::numeric_limits<double>::infinity();
        }
        const Eigen::VectorXd rinv_y = apply_ratio_inverse(eta, Eigen::MatrixXd(yd)).col(0);
        const Eigen::VectorXd beta = solver.solve(Xd.transpose() * rinv_y);
        const Eigen::VectorXd residual = yd - Xd * beta;
        const double qf = residual.dot(
            apply_ratio_inverse(eta, Eigen::MatrixXd(residual)).col(0));
        if (!(qf > 0.0) || !std::isfinite(qf)) {
            return std::numeric_limits<double>::infinity();
        }
        if (quadratic) *quadratic = qf;
        const double logdet_r = (n - rank) * std::log1p(eta * r0) +
            (1.0 + eta * D_.array()).log().sum();
        const double logdet_x = solver.vectorD().array().abs().max(1e-300).log().sum();
        return logdet_r + logdet_x + df * std::log(qf / static_cast<double>(df));
    };

    // Coarse grid over log(eta) followed by a bounded golden-section search,
    // matching the reference low-rank solver.
    const int grid_size = 81;
    int best_index = 0;
    double best_log_eta = -12.0;
    double best_objective = objective(best_log_eta, nullptr);
    for (int i = 1; i < grid_size; ++i) {
        const double candidate =
            -12.0 + 24.0 * i / static_cast<double>(grid_size - 1);
        const double value = objective(candidate, nullptr);
        if (value < best_objective) {
            best_objective = value;
            best_log_eta = candidate;
            best_index = i;
        }
    }
    double left = -12.0 + 24.0 * std::max(0, best_index - 1) /
        static_cast<double>(grid_size - 1);
    double right = -12.0 + 24.0 * std::min(grid_size - 1, best_index + 1) /
        static_cast<double>(grid_size - 1);

    const double golden = (std::sqrt(5.0) - 1.0) / 2.0;
    double c = right - golden * (right - left);
    double d = left + golden * (right - left);
    double fc = objective(c, nullptr);
    double fd = objective(d, nullptr);
    for (int iteration = 0; iteration < 80 && right - left > 1e-8; ++iteration) {
        if (fc < fd) {
            right = d; d = c; fd = fc; c = right - golden * (right - left);
            fc = objective(c, nullptr);
        } else {
            left = c; c = d; fc = fd; d = left + golden * (right - left);
            fd = objective(d, nullptr);
        }
    }

    best_log_eta = 0.5 * (left + right);
    double quadratic = 0.0;
    best_objective = objective(best_log_eta, &quadratic);
    if (!std::isfinite(best_objective)) {
        throw std::runtime_error("MatrixFreeVce: REML optimization failed");
    }

    null_delta_ = std::exp(-best_log_eta);
    sigma2_e_ = quadratic / static_cast<double>(df);
    sigma2_g_ = std::exp(best_log_eta) * sigma2_e_;
    log_likelihood_ = -0.5 * best_objective;
    converged_ = true;
    iterations_ = 1;
    return true;
}


void MatrixFreeVce::fixedComponents(double s2g, double s2e) {
    sigma2_g_ = s2g;
    sigma2_e_ = s2e;
    converged_ = true;
    iterations_ = 1;
}

void MatrixFreeVce::updateRemlDiagnostics() {
    // Placeholder
}

void MatrixFreeVce::applyKinship(const Eigen::VectorXd& v, Eigen::VectorXd& out, const class KinshipOperator& kinship, int n_threads) const {
    kinship.apply(v, const_cast<Eigen::VectorXd&>(out), 1);
}

void MatrixFreeVce::applyKinship(const Eigen::MatrixXd& X, Eigen::MatrixXd& out, const class KinshipOperator& kinship, int n_threads) const {
    kinship.apply(X, const_cast<Eigen::MatrixXd&>(out), 1);
}

double MatrixFreeVce::estimateTraceInverse(int n, int k, const class KinshipOperator& kinship, int n_threads) const {
    return 0.0;
}

} // namespace cosmic