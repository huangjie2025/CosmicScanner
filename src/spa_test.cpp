#include "spa_test.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>

namespace cosmic {
namespace {

double spa_logistic(double value) {
    if (value >= 0.0) {
        const double z = std::exp(-value);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(value);
    return z / (1.0 + z);
}

double spa_softplus(double value) {
    return value > 0.0 ? value + std::log1p(std::exp(-value))
                       : std::log1p(std::exp(value));
}

} // namespace

SPATest::SPATest(const Options& opts) : options_(opts) {}

void SPATest::set_null_model(const Eigen::VectorXd& mu, const Eigen::VectorXd& variance) {
    mu_ = mu;
    V_ = variance;
    valid_binary_null_ = (mu_.size() > 0 && mu_.size() == V_.size());
    for (int i = 0; i < mu_.size() && valid_binary_null_; ++i) {
        if (!std::isfinite(mu_(i)) || mu_(i) <= 1e-8 || mu_(i) >= 1.0 - 1e-8) valid_binary_null_ = false;
        if (!std::isfinite(V_(i)) || V_(i) <= 0.0) valid_binary_null_ = false;
    }
    if (valid_binary_null_) {
        std::cout << "[CosmicScanner] SPATest initialized with binary null model. N = " << mu_.size() << std::endl;
    }
}

double SPATest::compute_pvalue(const Eigen::VectorXd& X, double score_stat, double score_var) {
    // Calculate normal approximation p-value first
    double z = score_stat / std::sqrt(score_var);
    double pval_norm = 2.0 * pnorm(-std::abs(z));

    if (!valid_binary_null_) {
        return pval_norm;
    }

    // If the p-value is significant (e.g. < 0.05), we run SPA to correct it
    if (pval_norm > options_.pval_threshold) {
        return pval_norm;
    }

    // Run Saddle Point Approximation (SPA)
    // 1. Find Saddle point t*
    double t_star = find_saddle_point(score_stat, X);
    if (!std::isfinite(t_star)) return pval_norm;

    // 2. Compute K(t*), K''(t*)
    double k_t = K(t_star, X);
    double k_dp = K_double_prime(t_star, X);

    // 3. Compute w and v (Lugannani-Rice formula)
    double sign_t = t_star > 0 ? 1.0 : -1.0;
    const double radicand = 2.0 * (t_star * score_stat - k_t);
    if (!(radicand > 0.0) || !(k_dp > 0.0) || !std::isfinite(radicand) || !std::isfinite(k_dp))
        return pval_norm;
    double w = sign_t * std::sqrt(radicand);
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
        const double logit_mu = std::log(mu_(i)) - std::log1p(-mu_(i));
        k += std::log1p(-mu_(i)) + spa_softplus(logit_mu + t * X(i)) -
             t * X(i) * mu_(i);
    }
    return k;
}

double SPATest::K_prime(double t, const Eigen::VectorXd& X) const {
    // First derivative of CGF
    double k_p = 0.0;
    for (int i = 0; i < X.size(); ++i) {
        const double logit_mu = std::log(mu_(i)) - std::log1p(-mu_(i));
        const double tilted = spa_logistic(logit_mu + t * X(i));
        k_p += X(i) * (tilted - mu_(i));
    }
    return k_p;
}

double SPATest::K_double_prime(double t, const Eigen::VectorXd& X) const {
    // Second derivative of CGF (Variance at saddle point)
    double k_dp = 0.0;
    for (int i = 0; i < X.size(); ++i) {
        const double logit_mu = std::log(mu_(i)) - std::log1p(-mu_(i));
        const double tilted = spa_logistic(logit_mu + t * X(i));
        k_dp += X(i) * X(i) * tilted * (1.0 - tilted);
    }
    return k_dp;
}

double SPATest::find_saddle_point(double score, const Eigen::VectorXd& X) const {
    if (std::abs(score) < options_.tol) return 0.0;
    double lower = 0.0, upper = 0.0;
    if (score > 0.0) {
        upper = 1.0;
        while (upper < 1e6 && K_prime(upper, X) < score) upper *= 2.0;
        if (K_prime(upper, X) < score) return std::numeric_limits<double>::quiet_NaN();
    } else {
        lower = -1.0;
        while (lower > -1e6 && K_prime(lower, X) > score) lower *= 2.0;
        if (K_prime(lower, X) > score) return std::numeric_limits<double>::quiet_NaN();
    }
    double t = 0.5 * (lower + upper);
    for (int iter = 0; iter < options_.max_iter; ++iter) {
        double f = K_prime(t, X) - score;
        double f_prime = K_double_prime(t, X);
        if (std::abs(f) < options_.tol) return t;
        if (f > 0.0) upper = t;
        else lower = t;
        double candidate = t;
        if (f_prime > 1e-14) candidate = t - f / f_prime;
        if (!std::isfinite(candidate) || candidate <= lower || candidate >= upper)
            candidate = 0.5 * (lower + upper);
        if (std::abs(candidate - t) < options_.tol) return candidate;
        t = candidate;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double SPATest::pnorm(double x) const {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double SPATest::dnorm(double x) const {
    const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

} // namespace cosmic
