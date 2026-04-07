#include "spa_test.h"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace cosmic {

SPATest::SPATest(const Options& opts) : options_(opts) {}

void SPATest::set_null_model(const Eigen::VectorXd& mu, const Eigen::VectorXd& variance) {
    mu_ = mu;
    V_ = variance;
    std::cout << "[CosmicScanner] SPATest initialized with null model. N = " << mu_.size() << std::endl;
}

double SPATest::compute_pvalue(const Eigen::VectorXd& X, double score_stat, double score_var) {
    // Calculate normal approximation p-value first
    double z = score_stat / std::sqrt(score_var);
    double pval_norm = 2.0 * pnorm(-std::abs(z));

    // If the p-value is significant (e.g. < 0.05), we run SPA to correct it
    if (pval_norm > options_.pval_threshold) {
        return pval_norm;
    }

    // Run Saddle Point Approximation (SPA)
    // 1. Find Saddle point t*
    double t_star = find_saddle_point(score_stat, X);

    // 2. Compute K(t*), K''(t*)
    double k_t = K(t_star, X);
    double k_dp = K_double_prime(t_star, X);

    // 3. Compute w and v (Lugannani-Rice formula)
    double sign_t = t_star > 0 ? 1.0 : -1.0;
    double w = sign_t * std::sqrt(2.0 * (t_star * score_stat - k_t));
    double v = t_star * std::sqrt(k_dp);

    // Prevent division by zero
    if (std::abs(w) < 1e-8 || std::abs(v) < 1e-8) {
        return pval_norm;
    }

    // 4. SPA p-value: p = 1 - Phi(w) + phi(w)(1/w - 1/v)
    double abs_w = std::abs(w);
    double abs_v = std::abs(v);
    double pval_spa_one_sided = pnorm(-abs_w) + dnorm(abs_w) * (1.0 / abs_w - 1.0 / abs_v);
    double pval_spa = 2.0 * pval_spa_one_sided;
    
    // Ensure p-value is valid
    return std::min(1.0, std::max(0.0, pval_spa));
}

double SPATest::K(double t, const Eigen::VectorXd& X) const {
    // Skeleton: Cumulant Generating Function for Score Test
    // K(t) = \sum log(1 - mu_i + mu_i * e^{t X_i}) - t \sum X_i mu_i
    double k = 0.0;
    for (int i = 0; i < X.size(); ++i) {
        k += std::log(1.0 - mu_(i) + mu_(i) * std::exp(t * X(i))) - t * X(i) * mu_(i);
    }
    return k;
}

double SPATest::K_prime(double t, const Eigen::VectorXd& X) const {
    // First derivative of CGF
    double k_p = 0.0;
    for (int i = 0; i < X.size(); ++i) {
        double exp_tx = std::exp(t * X(i));
        double denom = 1.0 - mu_(i) + mu_(i) * exp_tx;
        k_p += (X(i) * mu_(i) * exp_tx) / denom - X(i) * mu_(i);
    }
    return k_p;
}

double SPATest::K_double_prime(double t, const Eigen::VectorXd& X) const {
    // Second derivative of CGF (Variance at saddle point)
    double k_dp = 0.0;
    for (int i = 0; i < X.size(); ++i) {
        double exp_tx = std::exp(t * X(i));
        double denom = 1.0 - mu_(i) + mu_(i) * exp_tx;
        double term1 = (X(i) * X(i) * mu_(i) * exp_tx) / denom;
        double term2 = std::pow((X(i) * mu_(i) * exp_tx) / denom, 2);
        k_dp += term1 - term2;
    }
    return k_dp;
}

double SPATest::find_saddle_point(double score, const Eigen::VectorXd& X) const {
    // Skeleton: Newton-Raphson to find root of K'(t) - score = 0
    double t = 0.0; // Initial guess
    for (int iter = 0; iter < options_.max_iter; ++iter) {
        double f = K_prime(t, X) - score;
        double f_prime = K_double_prime(t, X);
        
        if (std::abs(f_prime) < 1e-12) break; // Avoid division by zero
        
        double dt = f / f_prime;
        t -= dt;
        
        if (std::abs(dt) < options_.tol) break; // Converged
    }
    return t;
}

double SPATest::pnorm(double x) const {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double SPATest::dnorm(double x) const {
    const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

} // namespace cosmic