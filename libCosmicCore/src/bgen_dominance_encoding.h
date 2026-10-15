#ifndef COSMIC_BGEN_DOMINANCE_ENCODING_H
#define COSMIC_BGEN_DOMINANCE_ENCODING_H

// Probabilistic BGEN dominance encoding shared by the dominance-only GRM and
// marker scan. The additive-dominance scanner does not use this path.
//
// Math frozen in audit/BGEN_PROBABILISTIC_DOMINANCE_DESIGN.md:
//   dominance marker  Xd = P(AB) - mean(P(AB))
//   dominance GRM     K_D = sum_j w_j mu_j mu_j' / sum_j w_j 4 p_j^2 q_j^2
//   posterior-mean Vitezica code (counted allele frequency p, q = 1-p):
//        G=0 -> -2 p^2
//        G=1 ->  2 p q
//        G=2 -> -2 q^2
//        mu = -2 p^2 P(G=0) + 2 p q P(G=1) - 2 q^2 P(G=2)
//   No posterior-variance diagonal correction.

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace cosmic {
namespace bgen {

inline constexpr const char* kDominanceQcPolicyVersion = "BGEN_DOM_V1";
inline constexpr double kDominanceMaxMissing = 0.05;
inline constexpr double kDominanceMissingWarning = 0.02;
inline constexpr double kDominanceMinInfo = 0.01;
inline constexpr double kDominanceInfoWarning = 0.10;
inline constexpr double kDominanceEntropyWarning = 0.50;

struct DominanceQcMetrics {
    long long total = 0;
    long long observed = 0;
    double frequency = 0.0;
    double missing_rate = 1.0;
    double info_d = 0.0;
    double mean_entropy = 0.0;
    double mean_max_probability = 0.0;
};

// Hardcall Vitezica code, identical to the production BED/PGEN algorithm-4
// lookup table in GrmBuilder (grm_builder.cpp algorithm 4 pop_val_map).
// G is the number of copies of the counted allele.
inline double vitezicaCodeHardcall(int G, double p) {
    const double q = 1.0 - p;
    if (G == 0) return -2.0 * p * p;
    if (G == 1) return 2.0 * p * q;
    if (G == 2) return -2.0 * q * q;
    throw std::invalid_argument("Vitezica code requires G in {0,1,2}");
}

// Per-marker Vitezica normalization contribution sum_j 4 p^2 q^2.
inline double vitezicaDenominator(double p) {
    const double q = 1.0 - p;
    return 4.0 * p * p * q * q;
}

// Posterior-mean Vitezica code for one sample.
//
// p0/p1/p2 are P(G=0), P(G=1), P(G=2) for the counted allele; p is the
// counted-allele frequency. One-hot probabilities must reproduce
// vitezicaCodeHardcall() exactly.
inline double posteriorMeanVitezica(double p0, double p1, double p2, double p) {
    const double q = 1.0 - p;
    return -2.0 * p * p * p0 + 2.0 * p * q * p1 - 2.0 * q * q * p2;
}

// Probability QC is evaluated on the same analysis samples used by Xd and K_D.
// INFO_D is the observed-sample variance of posterior-mean Vitezica codes,
// divided by the hardcall HWE variance 4p^2q^2. Entropy is normalized to [0,1]
// by log(3) and is diagnostic only in policy BGEN_DOM_V1.
inline DominanceQcMetrics dominanceQcMetrics(const Eigen::VectorXf& p0,
                                             const Eigen::VectorXf& p1,
                                             const Eigen::VectorXf& p2,
                                             const Eigen::VectorXf& missing,
                                             const std::vector<int>& keep_indices = {}) {
    if (p0.size() != p1.size() || p1.size() != p2.size() || p2.size() != missing.size()) {
        throw std::invalid_argument("dominanceQcMetrics: size mismatch");
    }
    DominanceQcMetrics result;
    result.total = keep_indices.empty() ? p0.size() : static_cast<long long>(keep_indices.size());
    double dosage_sum = 0.0;
    double entropy_sum = 0.0;
    double max_probability_sum = 0.0;
    const auto accumulate = [&](int source) {
        if (source < 0 || source >= p0.size()) {
            throw std::out_of_range("dominanceQcMetrics: sample index out of range");
        }
        if (missing[source] != 0.0f) return;
        const double v0 = p0[source];
        const double v1 = p1[source];
        const double v2 = p2[source];
        dosage_sum += v1 + 2.0 * v2;
        double entropy = 0.0;
        if (v0 > 0.0) entropy -= v0 * std::log(v0);
        if (v1 > 0.0) entropy -= v1 * std::log(v1);
        if (v2 > 0.0) entropy -= v2 * std::log(v2);
        entropy_sum += entropy / std::log(3.0);
        max_probability_sum += std::max(v0, std::max(v1, v2));
        ++result.observed;
    };
    if (keep_indices.empty()) {
        for (int i = 0; i < p0.size(); ++i) accumulate(i);
    } else {
        for (int source : keep_indices) accumulate(source);
    }
    if (result.total <= 0 || result.observed <= 0) return result;

    result.missing_rate = static_cast<double>(result.total - result.observed) /
        static_cast<double>(result.total);
    result.frequency = dosage_sum / (2.0 * static_cast<double>(result.observed));
    result.mean_entropy = entropy_sum / static_cast<double>(result.observed);
    result.mean_max_probability = max_probability_sum / static_cast<double>(result.observed);

    const double denominator = vitezicaDenominator(result.frequency);
    if (denominator <= 0.0) return result;
    double mu_sum = 0.0;
    double mu_sq_sum = 0.0;
    const auto accumulate_mu = [&](int source) {
        if (missing[source] != 0.0f) return;
        const double mu = posteriorMeanVitezica(
            p0[source], p1[source], p2[source], result.frequency);
        mu_sum += mu;
        mu_sq_sum += mu * mu;
    };
    if (keep_indices.empty()) {
        for (int i = 0; i < p0.size(); ++i) accumulate_mu(i);
    } else {
        for (int source : keep_indices) accumulate_mu(source);
    }
    const double mean_mu = mu_sum / static_cast<double>(result.observed);
    const double variance_mu = std::max(0.0,
        mu_sq_sum / static_cast<double>(result.observed) - mean_mu * mean_mu);
    result.info_d = variance_mu / denominator;
    return result;
}

// Vector form. `missing` is 1.0 where the genotype is hard-missing; those
// entries are set to 0, matching the production hardcall missing code in
// GrmBuilder (missing hardcalls contribute zero to the GRM numerator).
inline Eigen::VectorXd posteriorMeanVitezicaVector(const Eigen::VectorXf& p0,
                                                   const Eigen::VectorXf& p1,
                                                   const Eigen::VectorXf& p2,
                                                   const Eigen::VectorXf& missing,
                                                   double p) {
    if (p0.size() != p1.size() || p1.size() != p2.size() || p2.size() != missing.size()) {
        throw std::invalid_argument("posteriorMeanVitezicaVector: size mismatch");
    }
    Eigen::VectorXd mu(p0.size());
    for (int i = 0; i < p0.size(); ++i) {
        mu[i] = (missing[i] != 0.0f)
                    ? 0.0
                    : posteriorMeanVitezica(p0[i], p1[i], p2[i], p);
    }
    return mu;
}

// Centered heterozygosity-probability dominance marker Xd = P(AB) - mean(P(AB)).
// Missing entries are mean-imputed (assigned the observed mean) and therefore
// contribute 0 after centering, matching prepareAssociationMarker() semantics.
inline Eigen::VectorXd centeredHeterozygosity(const Eigen::VectorXf& p1,
                                              const Eigen::VectorXf& missing) {
    if (p1.size() != missing.size()) {
        throw std::invalid_argument("centeredHeterozygosity: size mismatch");
    }
    double sum = 0.0;
    long long n = 0;
    for (int i = 0; i < p1.size(); ++i) {
        if (missing[i] != 0.0f) continue;
        sum += static_cast<double>(p1[i]);
        ++n;
    }
    if (n == 0) throw std::runtime_error("centeredHeterozygosity: no observed samples");
    const double mean = sum / static_cast<double>(n);
    Eigen::VectorXd xd(p1.size());
    for (int i = 0; i < p1.size(); ++i) {
        const double h = (missing[i] != 0.0f) ? mean : static_cast<double>(p1[i]);
        xd[i] = h - mean;
    }
    return xd;
}

} // namespace bgen
} // namespace cosmic

#endif // COSMIC_BGEN_DOMINANCE_ENCODING_H
