#ifndef COSMIC_SCANNER_SPA_TEST_H
#define COSMIC_SCANNER_SPA_TEST_H

#include <Eigen/Dense>
#include <vector>

namespace cosmic {

/**
 * @brief Saddle Point Approximation (SPA) Test Skeleton
 * Reference: SAIGE (Zhou et al., 2018), SPAGRM (Ma et al., 2020)
 * 
 * Core Logic:
 * 1. For highly unbalanced Case-Control studies, the score statistic S = X'(Y - mu) is not normally distributed.
 * 2. SPA uses the Cumulant Generating Function (CGF) of the Score statistic.
 * 3. Find the saddle point t* such that K'(t*) = S using Newton-Raphson or Secant method.
 * 4. Calculate the SPA p-value using the normal and standard normal distribution CDF/PDF around the saddle point.
 * 5. This prevents Type I error inflation for rare variants.
 */
class SPATest {
public:
    struct Options {
        double pval_threshold = 0.05; // Below this p-value from normal approx, we run SPA
        int max_iter = 100;
        double tol = 1e-5;
    };

    SPATest(const Options& opts);
    ~SPATest() = default;

    // Set the null model predictions (mu) and variance (V)
    void set_null_model(const Eigen::VectorXd& mu, const Eigen::VectorXd& variance);

    // Compute exact P-value for a given SNP (X) using SPA
    double compute_pvalue(const Eigen::VectorXd& X, double score_stat, double score_var);

private:
    Options options_;
    Eigen::VectorXd mu_;
    Eigen::VectorXd V_;
    bool valid_binary_null_ = false;

    // Calculate Cumulant Generating Function (K(t))
    double K(double t, const Eigen::VectorXd& X) const;

    // First derivative of CGF (K'(t))
    double K_prime(double t, const Eigen::VectorXd& X) const;

    // Second derivative of CGF (K''(t))
    double K_double_prime(double t, const Eigen::VectorXd& X) const;

    // Find the saddle point t* where K'(t) = S
    double find_saddle_point(double score, const Eigen::VectorXd& X) const;

    // Helper: Standard Normal CDF
    double pnorm(double x) const;
    // Helper: Standard Normal PDF
    double dnorm(double x) const;
};

} // namespace cosmic

#endif // COSMIC_SCANNER_SPA_TEST_H
