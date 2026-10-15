#include "rare_variant_extensions.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cosmic {
namespace rare_variant {
namespace {

constexpr long double kPi = 3.141592653589793238462643383279502884L;

double inverseNormal(double p) {
    // Acklam's rational approximation; Halton points never include 0 or 1.
    static constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01,
        2.506628277459239e+00};
    static constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
    static constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00,
        2.938163982698783e+00};
    static constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00};
    if (p < 0.02425) {
        const double q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    if (p > 0.97575) return -inverseNormal(1.0 - p);
    const double q = p - 0.5, r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
}

double radicalInverse(unsigned long long index, unsigned base) {
    double value = 0.0, factor = 1.0 / base;
    while (index) {
        value += factor * (index % base);
        index /= base;
        factor /= base;
    }
    return value;
}

std::vector<unsigned> primes(int count) {
    std::vector<unsigned> result;
    for (unsigned candidate = 2; static_cast<int>(result.size()) < count; ++candidate) {
        bool prime = true;
        for (unsigned divisor = 2; divisor * divisor <= candidate; ++divisor)
            if (candidate % divisor == 0) { prime = false; break; }
        if (prime) result.push_back(candidate);
    }
    return result;
}

double probability(long double p) {
    // Positive mathematical tails must not become reported zeros on underflow.
    return static_cast<double>(std::min(1.0L, std::max(
        static_cast<long double>(std::numeric_limits<double>::denorm_min()), p)));
}

void validateWeights(const Eigen::VectorXd& weights, Eigen::Index n) {
    if (weights.size() != n || !weights.allFinite() || (weights.array() < 0.0).any())
        throw std::invalid_argument("Rare-variant weights must be finite, nonnegative and match input size");
    if (n && weights.maxCoeff() == 0.0)
        throw std::invalid_argument("Rare-variant weights must have positive total weight");
}

struct ScoreInput {
    Eigen::VectorXd scores;
    Eigen::VectorXd coordinates;
    Eigen::MatrixXd covariance;
    Eigen::MatrixXd root;
    double covariance_scale = 1.0;
};

ScoreInput validateScores(const Eigen::VectorXd& scores, const Eigen::MatrixXd& covariance,
                          const Eigen::VectorXd& weights) {
    const Eigen::Index n = scores.size();
    validateWeights(weights, n);
    if (covariance.rows() != n || covariance.cols() != n ||
        !scores.allFinite() || !covariance.allFinite())
        throw std::invalid_argument("Rare-variant scores/covariance have invalid dimensions or nonfinite values");
    ScoreInput input;
    input.scores = scores;
    input.coordinates = Eigen::VectorXd::Zero(n);
    input.covariance = covariance;
    input.root = Eigen::MatrixXd::Zero(n, n);
    if (!n) return input;
    const double scale = covariance.cwiseAbs().maxCoeff();
    if (scale == 0.0) {
        if ((scores.array() != 0.0).any())
            throw std::invalid_argument("Nonzero scores with zero covariance");
        return input;
    }
    input.covariance_scale = scale;
    input.covariance /= scale;
    input.scores /= std::sqrt(scale);
    if (!input.scores.allFinite()) throw std::overflow_error("Scaled rare-variant scores overflow");
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * n;
    if ((input.covariance - input.covariance.transpose()).cwiseAbs().maxCoeff() > tolerance)
        throw std::invalid_argument("Score covariance must be symmetric");
    input.covariance = (0.5 * (input.covariance + input.covariance.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(input.covariance);
    if (es.info() != Eigen::Success) throw std::runtime_error("Score covariance eigendecomposition failed");
    const double cutoff = tolerance * std::max(1.0, es.eigenvalues().maxCoeff());
    if (es.eigenvalues().minCoeff() < -cutoff)
        throw std::invalid_argument("Score covariance must be positive semidefinite");
    Eigen::VectorXd eigenvalues = es.eigenvalues().cwiseMax(0.0);
    const double score_scale = input.scores.cwiseAbs().maxCoeff();
    const Eigen::VectorXd scaled_scores = input.scores / (score_scale > 0.0 ? score_scale : 1.0);
    const Eigen::VectorXd rotated = es.eigenvectors().transpose() * scaled_scores;
    for (Eigen::Index i = 0; i < n; ++i) {
        // Retain every positive mode: even a tiny variance can dominate after
        // weighting. Only zero/roundoff-negative modes have no null support.
        if (eigenvalues[i] == 0.0) {
            if (std::abs(rotated[i]) > 1e-7 * scaled_scores.norm())
                throw std::invalid_argument("Scores lie outside the numerical covariance support");
        } else {
            input.coordinates[i] = static_cast<double>(static_cast<long double>(rotated[i]) * score_scale /
                                                       std::sqrt(static_cast<long double>(eigenvalues[i])));
        }
        if (covariance(i, i) == 0.0 && scores[i] != 0.0)
            throw std::invalid_argument("Nonzero score at a zero-variance variant");
    }
    if (!input.coordinates.allFinite()) throw std::overflow_error("Whitened rare-variant scores overflow");
    input.root = es.eigenvectors() * eigenvalues.cwiseSqrt().asDiagonal();
    return input;
}

// Integer/half-integer gamma survival, using only positive recurrence terms.
long double chiSquareTail(long double q, Eigen::Index df) {
    if (q <= 0.0L) return 1.0L;
    const long double x = q / 2.0L;
    long double a = (df % 2) ? 0.5L : 1.0L;
    long double tail = (df % 2) ? std::erfc(std::sqrt(x)) : std::exp(-x);
    for (; a < static_cast<long double>(df) / 2.0L; a += 1.0L)
        tail += std::exp(a * std::log(x) - x - std::lgamma(a + 1.0L));
    return std::min(1.0L, tail);
}

} // namespace

MixtureTailResult mixtureChiSquareUpperTail(double q, const Eigen::VectorXd& lambda) {
    if (!std::isfinite(q) || q < 0.0 || !lambda.allFinite() || (lambda.array() < 0.0).any())
        throw std::invalid_argument("Mixture statistic/eigenvalues must be finite and nonnegative");
    if (q == 0.0) return {};
    std::vector<long double> values;
    for (Eigen::Index i = 0; i < lambda.size(); ++i)
        if (lambda[i] > 0.0) values.push_back(lambda[i]);
    if (values.empty()) return {0.0, false};
    const long double maximum = *std::max_element(values.begin(), values.end());
    const long double minimum = *std::min_element(values.begin(), values.end());
    const long double observed = static_cast<long double>(q) / maximum;
    if (minimum == maximum)
        return {probability(chiSquareTail(observed, values.size())), false};

    // Any t in (0,1/(2*lambda_max)) gives a rigorous MGF/Markov upper bound.
    // Optimize in lambda_max units; bisection accuracy affects power, not validity.
    long double mean = 0.0L;
    for (auto& value : values) { value /= maximum; mean += value; }
    long double bound = 1.0L;
    if (observed > mean) {
        long double lo = 0.0L, hi = 0.5L;
        for (int iteration = 0; iteration < 100; ++iteration) {
            const long double t = (lo + hi) / 2.0L;
            long double derivative = -observed;
            for (const auto value : values) derivative += value / (1.0L - 2.0L * t * value);
            if (derivative > 0.0L) hi = t; else lo = t;
        }
        long double log_bound = -lo * observed;
        for (const auto value : values) log_bound -= 0.5L * std::log1p(-2.0L * lo * value);
        bound = std::exp(log_bound);
    }

    // Ruben's positive series: sum c_j * P(chi-square_(k+2j) >= q/beta),
    // beta=min(lambda). The c_j are nonnegative and sum to one. Bounding all
    // omitted survival terms by one makes truncation conservative, even when
    // an ill-conditioned spectrum prevents useful numerical convergence.
    constexpr int kTerms = 512;
    const long double beta = minimum / maximum;
    long double log_c0 = 0.0L;
    std::vector<long double> ratios, powers(values.size(), 1.0L);
    for (const auto value : values) {
        log_c0 += 0.5L * std::log(beta / value);
        ratios.push_back(1.0L - beta / value);
    }
    std::vector<long double> coefficients(kTerms), recurrence(kTerms);
    coefficients[0] = std::exp(log_c0);
    if (coefficients[0] > 0.0L) {
        const long double x = observed / (2.0L * beta);
        const long double a = static_cast<long double>(values.size()) / 2.0L;
        long double tail = chiSquareTail(observed / beta, values.size());
        long double mass = coefficients[0], partial = coefficients[0] * tail;
        for (int j = 1; j < kTerms; ++j) {
            for (size_t i = 0; i < values.size(); ++i) {
                powers[i] *= ratios[i];
                recurrence[j] += 0.5L * powers[i];
            }
            for (int k = 1; k <= j; ++k)
                coefficients[j] += recurrence[k] * coefficients[j - k];
            coefficients[j] /= j;
            tail = std::min(1.0L, tail + std::exp((a + j - 1.0L) * std::log(x) - x - std::lgamma(a + j)));
            mass += coefficients[j];
            partial += coefficients[j] * tail;
            if (1.0L - mass < 1e-16L) break;
        }
        const long double roundoff = 64.0L * std::numeric_limits<long double>::epsilon() *
                                     (kTerms + values.size());
        bound = std::min(bound, partial + std::max(0.0L, 1.0L - mass) + roundoff);
    }
    return {probability(bound * (1.0L + 64.0L * std::numeric_limits<double>::epsilon())), true};
}

SkatOGridResult skatOConservativeGrid(const Eigen::VectorXd& scores,
                                     const Eigen::MatrixXd& covariance,
                                     const Eigen::VectorXd& weights,
                                     const std::vector<double>& rho_grid) {
    if (rho_grid.empty()) throw std::invalid_argument("SKAT-O rho grid must not be empty");
    for (const auto rho : rho_grid)
        if (!std::isfinite(rho) || rho < 0.0 || rho > 1.0)
            throw std::invalid_argument("SKAT-O rho must be finite and in [0,1]");
    const auto input = validateScores(scores, covariance, weights);
    SkatOGridResult result;
    result.rho_grid = rho_grid;
    std::sort(result.rho_grid.begin(), result.rho_grid.end());
    result.rho_grid.erase(std::unique(result.rho_grid.begin(), result.rho_grid.end()), result.rho_grid.end());
    result.selected_rho = result.rho_grid.front();
    if (!scores.size()) {
        result.rho_statistics.assign(result.rho_grid.size(), 0.0);
        result.rho_p_values.assign(result.rho_grid.size(), 1.0);
        return result;
    }
    const Eigen::VectorXd w = weights / weights.maxCoeff();
    const Eigen::MatrixXd h = w.asDiagonal() * input.root;
    // Use the same covariance support for Q and its null kernel. Otherwise a
    // tolerated nullspace rounding residual can give Q>0 with exactly zero variance.
    const Eigen::VectorXd weighted_scores = h * input.coordinates;
    const Eigen::MatrixXd skat_kernel = h.transpose() * h;
    const Eigen::VectorXd burden_root = h.colwise().sum().transpose();
    const double skat_q = weighted_scores.squaredNorm();
    const double burden_q = std::pow(burden_root.dot(input.coordinates), 2.0);
    if (!std::isfinite(skat_q) || !std::isfinite(burden_q))
        throw std::overflow_error("Rare-variant quadratic statistic overflow");
    for (Eigen::Index i = 0; i < scores.size(); ++i)
        if (weights[i] > 0.0 && covariance(i, i) > 0.0) ++result.n_variants;
    for (const auto rho : result.rho_grid) {
        const double q = (1.0 - rho) * skat_q + rho * burden_q;
        Eigen::VectorXd lambda;
        if (rho == 1.0) {
            lambda = Eigen::VectorXd::Constant(1, burden_root.squaredNorm());
        } else {
            const Eigen::MatrixXd kernel = (1.0 - rho) * skat_kernel +
                                          rho * burden_root * burden_root.transpose();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kernel, Eigen::EigenvaluesOnly);
            if (es.info() != Eigen::Success) throw std::runtime_error("SKAT-O kernel eigendecomposition failed");
            lambda = es.eigenvalues().cwiseMax(0.0);
        }
        const double p = mixtureChiSquareUpperTail(q, lambda).p_value;
        const long double original_q = static_cast<long double>(q) * input.covariance_scale;
        if (original_q > std::numeric_limits<double>::max())
            throw std::overflow_error("Rare-variant output statistic overflow");
        result.rho_statistics.push_back(static_cast<double>(original_q));
        result.rho_p_values.push_back(p);
        if (result.rho_p_values.size() == 1 || p < result.min_p_value) {
            result.min_p_value = p;
            result.selected_rho = rho;
            result.statistic = result.rho_statistics.back();
        }
    }
    result.p_value = std::min(1.0, result.min_p_value * result.rho_grid.size());
    return result;
}

SkatOResult skatOCanonical(const Eigen::VectorXd& scores,
                           const Eigen::MatrixXd& covariance,
                           const Eigen::VectorXd& weights,
                           const std::vector<double>& rho_grid,
                           int null_draws) {
    if (null_draws < 1024) throw std::invalid_argument("Canonical SKAT-O requires at least 1024 null draws");
    const auto marginal = skatOConservativeGrid(scores, covariance, weights, rho_grid);
    SkatOResult result;
    static_cast<SkatOGridResult&>(result) = marginal;
    result.p_value = marginal.min_p_value;
    result.null_draws = null_draws;
    if (scores.size() <= 1 || result.rho_grid.size() == 1 || result.n_variants <= 1) return result;

    const auto input = validateScores(scores, covariance, weights);
    const Eigen::VectorXd w = weights / weights.maxCoeff();
    const Eigen::MatrixXd h = w.asDiagonal() * input.root;
    const int rank = input.root.cols();
    const auto bases = primes(rank);
    const int rhos = static_cast<int>(result.rho_grid.size());
    Eigen::MatrixXd simulated(null_draws, rhos);
    Eigen::VectorXd z(rank);
    for (int draw = 0; draw < null_draws; ++draw) {
        for (int j = 0; j < rank; ++j) {
            double u = radicalInverse(static_cast<unsigned long long>(draw) + 1, bases[j]);
            u = std::clamp(u, std::numeric_limits<double>::epsilon(),
                           1.0 - std::numeric_limits<double>::epsilon());
            z[j] = inverseNormal(u);
        }
        const Eigen::VectorXd wu = h * z;
        const double skat = wu.squaredNorm();
        const double burden = std::pow(wu.sum(), 2.0);
        for (int r = 0; r < rhos; ++r)
            simulated(draw, r) = (1.0 - result.rho_grid[r]) * skat + result.rho_grid[r] * burden;
    }

    // Marginal empirical tails preserve each rho's scale; the second pass
    // integrates the union event under the same shared Gaussian score vector.
    std::vector<double> observed_p(rhos);
    double min_p = 1.0;
    int best = 0;
    int best_exceed = null_draws;
    for (int r = 0; r < rhos; ++r) {
        int exceed = 0;
        const double observed = result.rho_statistics[r] / input.covariance_scale;
        for (int draw = 0; draw < null_draws; ++draw) exceed += simulated(draw, r) >= observed;
        observed_p[r] = (exceed + 0.5) / (null_draws + 1.0);
        if (observed_p[r] < min_p) { min_p = observed_p[r]; best = r; best_exceed = exceed; }
    }
    if (best_exceed == 0) {
        result.p_value = std::min(1.0, marginal.min_p_value * rhos);
        result.p_method = "BONFERRONI_FALLBACK";
        result.fallback_used = true;
        return result;
    }
    std::vector<double> thresholds(rhos);
    for (int r = 0; r < rhos; ++r) {
        std::vector<double> values(simulated.col(r).data(), simulated.col(r).data() + null_draws);
        const int index = std::clamp(static_cast<int>(std::floor((1.0 - min_p) * null_draws)), 0, null_draws - 1);
        std::nth_element(values.begin(), values.begin() + index, values.end());
        thresholds[r] = values[index];
    }
    int union_exceed = 0;
    for (int draw = 0; draw < null_draws; ++draw) {
        bool exceed = false;
        for (int r = 0; r < rhos; ++r) exceed = exceed || simulated(draw, r) >= thresholds[r];
        union_exceed += exceed;
    }
    result.rho_p_values = observed_p;
    result.min_p_value = min_p;
    result.selected_rho = result.rho_grid[best];
    result.statistic = result.rho_statistics[best];
    result.p_value = std::min(1.0, (union_exceed + 0.5) / (null_draws + 1.0));
    if (!std::isfinite(result.p_value)) {
        result.p_value = std::min(1.0, min_p * rhos);
        result.p_method = "BONFERRONI_FALLBACK";
        result.fallback_used = true;
    }
    return result;
}

AcatVResult acatV(const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights) {
    validateWeights(weights, p_values.size());
    if (!p_values.allFinite() || (p_values.array() < 0.0).any() || (p_values.array() > 1.0).any())
        throw std::invalid_argument("ACAT-V p-values must be finite and in [0,1]");
    AcatVResult result;
    if (!p_values.size()) return result;
    bool zero = false, one = false;
    long double weight_sum = 0.0L, scale = 0.5L;
    const long double max_weight = weights.maxCoeff();
    for (Eigen::Index i = 0; i < p_values.size(); ++i) {
        if (weights[i] == 0.0) continue;
        ++result.n_variants;
        zero = zero || p_values[i] == 0.0;
        one = one || p_values[i] == 1.0;
        weight_sum += weights[i] / max_weight;
        scale = std::min(scale, std::min(static_cast<long double>(p_values[i]), 1.0L - p_values[i]));
    }
    if (zero && one) throw std::invalid_argument("ACAT-V cannot combine both exact p=0 and p=1 with positive weights");
    if (zero || one) {
        result.p_value = zero ? 0.0 : 1.0;
        result.statistic = zero ? std::numeric_limits<long double>::infinity() :
                                 -std::numeric_limits<long double>::infinity();
        return result;
    }
    long double sum = 0.0L, correction = 0.0L;
    for (Eigen::Index i = 0; i < p_values.size(); ++i) {
        if (weights[i] == 0.0 || p_values[i] == 0.5) continue;
        const long double p = p_values[i];
        // Scale cotangents BEFORE division to avoid 1/p overflow, including
        // subnormal p. Reflect the upper half instead of evaluating tan near pi.
        const long double cot = p < 0.5L ? scale / std::tan(kPi * p) :
                                         -scale / std::tan(kPi * (1.0L - p));
        const long double term = (weights[i] / max_weight) / weight_sum * cot - correction;
        const long double next = sum + term;
        correction = (next - sum) - term;
        sum = next;
    }
    result.statistic = sum / scale;
    result.p_value = probability(std::atan2(scale, sum) / kPi);
    return result;
}

AcatVResult acatVFromScores(const Eigen::VectorXd& scores, const Eigen::MatrixXd& covariance,
                           const Eigen::VectorXd& weights) {
    const auto input = validateScores(scores, covariance, weights);
    std::vector<double> p_values, active_weights;
    for (Eigen::Index i = 0; i < scores.size(); ++i) {
        if (weights[i] == 0.0 || covariance(i, i) <= 0.0) continue;
        const long double z = static_cast<long double>(input.scores[i]) /
                              std::sqrt(static_cast<long double>(input.covariance(i, i)));
        p_values.push_back(probability(std::erfc(std::abs(z) / std::sqrt(2.0L))));
        active_weights.push_back(weights[i]);
    }
    return acatV(Eigen::Map<const Eigen::VectorXd>(p_values.data(), p_values.size()),
                 Eigen::Map<const Eigen::VectorXd>(active_weights.data(), active_weights.size()));
}

} // namespace rare_variant
} // namespace cosmic
