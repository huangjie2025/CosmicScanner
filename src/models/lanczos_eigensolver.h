#pragma once
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <random>
#include "plink_reader.h"

namespace cosmic {

struct LanczosResult {
    Eigen::VectorXd eigenvalues;   // descending order
    Eigen::MatrixXd eigenvectors;  // n x m (matching eigenvalues)
    int iterations = 0;             // actual Lanczos steps
    int convergedritz = 0;          // converged Ritz values
    double trace_g = 0.0;           // trace(G)
    double target_variance = 0.0;   // requested target
    int K_target = -1;              // K for target variance
    int K98 = -1;
    int K99 = -1;
    std::string status = "OK";
};

class MatrixFreeGrmOperator {
public:
    MatrixFreeGrmOperator(PlinkReader& reader,
                          const std::vector<int>& keep_indices,
                          double S_2pq,
                          double tune_a = 0.0,
                          double tune_b = 1.0,
                          double blend_alpha = 0.0,
                          const Eigen::MatrixXd* blend_target = nullptr,
                          int n_threads = 1)
        : reader_(reader), keep_indices_(keep_indices),
          S_2pq_(S_2pq), tune_a_(tune_a), tune_b_(tune_b),
          blend_alpha_(blend_alpha), blend_target_(blend_target),
          n_threads_(n_threads) {}

    // Production dense pipeline (single_step_model.cpp buildBacksolveGrm):
    //   G_tuned  = tune_b_ * G_raw + tune_a_ * J        (J = all-ones, preGSf90 tuneG)
    //   G_blend  = (1-alpha) * G_tuned + alpha * T      (T = I when blend_target_==nullptr,
    //                                                     else T = *blend_target_ = A22)
    // Matrix-free apply derives to:
    //   y = (1-alpha) * ( tune_b_ * G_raw*v + tune_a_ * sum(v) * ones ) + alpha * T*v
    void apply(const Eigen::VectorXd& v, Eigen::VectorXd& y) const {
        const int n = static_cast<int>(keep_indices_.size());
        const int m = reader_.getNumSnps();
        if (y.size() != n) y.resize(n);

        // Step 1: tmp = Z_centered' * v (M-vector)
        //   Z_standardized = (X - 2p) * inv_std ; so scaled zt recovers Z_centered' v
        Eigen::VectorXd zt_std;
        reader_.multiply_Zt_v(v, zt_std, keep_indices_, n_threads_);

        const auto& inv_stds = reader_.getInvStds();
        Eigen::VectorXd tmp(m);
        for (int j = 0; j < m; ++j) {
            double std_j = (inv_stds[j] > 0.0) ? (1.0 / inv_stds[j]) : 1.0;
            tmp(j) = std_j * zt_std(j);
        }
        //   tmp = Z_centered' * v

        // Step 2: y = Z_centered * tmp / S (raw VanRaden GRM application)
        Eigen::VectorXd weighted(m);
        for (int j = 0; j < m; ++j) {
            double std_j = (inv_stds[j] > 0.0) ? (1.0 / inv_stds[j]) : 1.0;
            weighted(j) = std_j * tmp(j);
        }
        reader_.multiply_Z_v(weighted, y, keep_indices_, n_threads_);

        if (S_2pq_ > 0.0) {
            y /= S_2pq_;
        }
        //   y = G_raw * v

        // Step 3: tuneG  (G_tuned = b*G_raw + a*J, J = ones;  J*v = sum(v)*ones)
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

    int rows() const { return static_cast<int>(keep_indices_.size()); }

private:
    PlinkReader& reader_;
    std::vector<int> keep_indices_;
    double S_2pq_;
    double tune_a_;
    double tune_b_;
    double blend_alpha_;
    const Eigen::MatrixXd* blend_target_;
    int n_threads_;
};

class LanczosEigenSolver {
public:
    struct Options {
        int max_iterations = 80;
        double tolerance = 1e-6;
        unsigned long seed = 42;
        int n_threads = 1;
    };

    LanczosEigenSolver() = default;

template <typename Operator>
LanczosResult compute(
    const Operator& G_op,
    double target_variance,
    const Options& opts) const {

        const int n = G_op.rows();
        LanczosResult result;
        result.target_variance = target_variance;

        if (n < 1) {
            result.status = "ERROR: empty operator";
            return result;
        }

        const int m = std::min(opts.max_iterations, n);

        Eigen::MatrixXd V(n, m + 1);
        Eigen::VectorXd alpha = Eigen::VectorXd::Zero(m);
        Eigen::VectorXd beta = Eigen::VectorXd::Zero(m);

        Eigen::VectorXd z(n);
        {
            std::mt19937_64 rng(opts.seed);
            std::normal_distribution<double> dist(0.0, 1.0);
            for (int i = 0; i < n; ++i) z(i) = dist(rng);
        }
        double norm_z = z.norm();
        if (norm_z < 1e-15) {
            z = Eigen::VectorXd::Ones(n) / std::sqrt(static_cast<double>(n));
            norm_z = 1.0;
        }
        V.col(0) = z / norm_z;

        int actual_m = m;
        for (int j = 0; j < m; ++j) {
            Eigen::VectorXd w(n);
            G_op.apply(V.col(j), w);

            if (j > 0) {
                w -= beta(j - 1) * V.col(j - 1);
            }

            alpha(j) = w.dot(V.col(j));
            w -= alpha(j) * V.col(j);

            for (int k = 0; k <= j; ++k) {
                double proj = w.dot(V.col(k));
                w -= proj * V.col(k);
            }

            beta(j) = w.norm();

            if (beta(j) < 1e-12) {
                actual_m = j + 1;
                break;
            }

            V.col(j + 1) = w / beta(j);
        }

        result.iterations = actual_m;

        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(actual_m, actual_m);
        for (int j = 0; j < actual_m; ++j) {
            T(j, j) = alpha(j);
            if (j < actual_m - 1) {
                T(j + 1, j) = beta(j);
                T(j, j + 1) = beta(j);
            }
        }

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
        if (es.info() != Eigen::Success) {
            result.status = "ERROR: tridiagonal eigendecomposition failed";
            return result;
        }

        Eigen::VectorXd evals = es.eigenvalues();
         Eigen::MatrixXd evecs_T = es.eigenvectors();  // m x m

        // Sort eigenvalues and eigenvectors in descending order
        std::vector<std::pair<double, int>> idx_val;
        for (int i = 0; i < evals.size(); ++i) {
            idx_val.emplace_back(evals(i), i);
        }
        std::sort(idx_val.begin(), idx_val.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

         Eigen::VectorXd sorted_evals(actual_m);
        Eigen::MatrixXd sorted_evecs_T(actual_m, actual_m);
        for (int i = 0; i < actual_m; ++i) {
            sorted_evals(i) = idx_val[i].first;
            sorted_evecs_T.col(i) = evecs_T.col(idx_val[i].second);
        }

         // Recover eigenvectors of original matrix: V_m * Q_T
        // where V_m is n x m Lanczos vectors, Q_T is m x m eigenvectors of T
        Eigen::MatrixXd V_m = V.leftCols(actual_m);  // n x m
        result.eigenvectors = V_m * sorted_evecs_T;  // n x m

        result.eigenvalues = sorted_evals;

        // trace(G) = sum of ALL Ritz eigenvalues. Full sum (including any tiny
        // negatives) mirrors the dense path (single_step_model.cpp uses
        // ev.sum()) so that the K98/K99 variance denominator is identical.
        result.trace_g = evals.sum();  // Use original unsorted sum for trace
        result.convergedritz = actual_m;

        // Use the SAME denominator as dense: full trace including negatives.
        // Only positive eigenvalues contribute to cumulative sum.
        // Sorted eigenvalues (descending) accumulate largest first, matching dense.
        if (result.trace_g <= 0.0) {
            result.status = "ERROR: non-positive trace(G)";
            result.K_target = n;
            result.K98 = n;
            result.K99 = n;
            return result;
        }

        auto findK = [&](double target) -> int {
            double cum = 0.0;
            for (int i = 0; i < sorted_evals.size(); ++i) {
                double val = std::max(0.0, sorted_evals(i));
                cum += val;
                if (cum >= target * result.trace_g) return i + 1;
            }
            return n;
        };

        result.K98 = findK(0.98);
        result.K99 = findK(0.99);
        result.K_target = findK(target_variance);
        result.status = "OK";

        return result;
    }
};

} // namespace cosmic
