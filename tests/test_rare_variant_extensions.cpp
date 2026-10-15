#include "models/rare_variant/rare_variant_extensions.h"
#include "rare_variant_test.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace cosmic::rare_variant;

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance * std::max(std::abs(expected), 1e-300),
          message + ": actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected));
}

template <typename F> void invalid(F operation, const std::string& message) {
    bool caught = false;
    try { operation(); } catch (const std::invalid_argument&) { caught = true; }
    check(caught, message);
}

VectorXd vector(std::initializer_list<double> entries) {
    VectorXd result(entries.size());
    Eigen::Index i = 0;
    for (double entry : entries) result[i++] = entry;
    return result;
}

void testMixture() {
    close(mixtureChiSquareUpperTail(4.0, vector({1.0})).p_value,
          std::erfc(std::sqrt(2.0)), 1e-14, "Rank-one chi-square tail");
    close(mixtureChiSquareUpperTail(10.0, vector({2.0, 2.0})).p_value,
          std::exp(-2.5), 1e-14, "Scaled chi-square_2 tail");
    close(mixtureChiSquareUpperTail(10.0, vector({1.0, 1.0, 1.0, 1.0})).p_value,
          std::exp(-5.0) * 6.0, 1e-14, "Chi-square_4 recurrence");
    check(!mixtureChiSquareUpperTail(4.0, vector({0.0, 1.0})).conservative_bound,
          "Zero eigenvalues do not change rank-one calibration");
    double previous = 1.0;
    for (double q : {0.1, 1.0, 5.0, 20.0, 50.0, 1000.0}) {
        const auto result = mixtureChiSquareUpperTail(q, vector({1.0, 1.0, 2.0, 2.0}));
        const double exact = 2.0 * std::exp(-q / 4.0) - std::exp(-q / 2.0);
        check(result.conservative_bound && result.p_value >= exact * (1.0 - 1e-13),
              "Unequal-spectrum bound must dominate exact convolution tail");
        check(result.p_value <= previous, "Mixture tail monotone in Q");
        if (q <= 20.0) close(result.p_value, exact, 1e-11, "Positive-series calibration");
        previous = result.p_value;
    }
    close(mixtureChiSquareUpperTail(5e-200, vector({1e-200, 1e-200, 2e-200, 2e-200})).p_value,
          mixtureChiSquareUpperTail(5.0, vector({1.0, 1.0, 2.0, 2.0})).p_value, 1e-13,
          "Mixture scale invariance");
    check(mixtureChiSquareUpperTail(2000.0, vector({1.0})).p_value > 0.0, "Tail underflow is not zero");
    check(mixtureChiSquareUpperTail(0.0, VectorXd()).p_value == 1.0, "Empty point mass at zero");
    check(mixtureChiSquareUpperTail(1.0, VectorXd()).p_value == 0.0, "Point mass positive threshold");
    const auto ill = mixtureChiSquareUpperTail(20.0, vector({1e-100, 1.0}));
    check(ill.p_value >= std::erfc(std::sqrt(10.0)) && ill.p_value < 0.001,
          "Ill-conditioned mixture has a useful conservative bound");
    invalid([] { mixtureChiSquareUpperTail(-1.0, vector({1.0})); }, "Reject negative Q");
    invalid([] { mixtureChiSquareUpperTail(1.0, vector({-1.0})); }, "Reject negative eigenvalues");
    invalid([] { mixtureChiSquareUpperTail(INFINITY, vector({1.0})); }, "Reject infinite Q");
    invalid([] { mixtureChiSquareUpperTail(1.0, vector({NAN})); }, "Reject NaN eigenvalues");
}

void testGrid() {
    const VectorXd scores = vector({3.0, -1.0});
    const VectorXd weights = vector({1.0, 1.0});
    const MatrixXd covariance = MatrixXd::Identity(2, 2);
    const auto skat = skatOConservativeGrid(scores, covariance, weights, {0.0});
    close(skat.statistic, 10.0, 1e-14, "rho=0 SKAT statistic");
    close(skat.p_value, std::exp(-5.0), 1e-13, "rho=0 exact equal-spectrum SKAT");
    const auto burden = skatOConservativeGrid(scores, covariance, weights, {1.0});
    close(burden.statistic, 4.0, 1e-14, "rho=1 burden Q");
    close(burden.p_value, std::erfc(1.0), 1e-13, "rho=1 exact burden tail");
    const auto middle = skatOConservativeGrid(scores, covariance, weights, {0.5});
    close(middle.statistic, 7.0, 1e-14, "Intermediate-rho quadratic statistic");
    close(middle.p_value, mixtureChiSquareUpperTail(7.0, vector({0.5, 1.5})).p_value,
          1e-13, "Intermediate-rho kernel eigenvalues");
    const auto grid = skatOConservativeGrid(scores, covariance, weights, {1.0, 0.0, 0.0});
    check(grid.rho_grid.size() == 2 && grid.selected_rho == 0.0, "Deduplicate/sort grid and select rho");
    close(grid.p_value, 2.0 * skat.p_value, 1e-13, "Bonferroni correction, not minimum p");
    check(grid.p_value > grid.min_p_value, "Never expose unadjusted minimum as omnibus p");
    const auto one = skatOConservativeGrid(vector({3.0}), MatrixXd::Identity(1, 1), vector({7.0}));
    for (double p : one.rho_p_values)
        close(p, std::erfc(3.0 / std::sqrt(2.0)), 1e-13, "Single-variant rho invariance");
    close(one.p_value, one.rho_grid.size() * one.min_p_value, 1e-13, "Single-variant conservative grid penalty");
    const auto weighted = skatOConservativeGrid(scores, covariance, vector({1.0, 2.0}));
    const auto rescaled = skatOConservativeGrid(scores * 1e-100, covariance * 1e-200, vector({1e200, 2e200}));
    close(rescaled.p_value, weighted.p_value, 1e-12, "Score/covariance and common weight scaling invariance");
    MatrixXd tiny_variance = MatrixXd::Identity(2, 2);
    tiny_variance(0, 0) = 1e-200;
    close(skatOConservativeGrid(vector({3e-100, 0.0}), tiny_variance, vector({1.0, 0.0}), {0.0}).p_value,
          std::erfc(3.0 / std::sqrt(2.0)), 1e-13, "Do not discard small positive covariance modes before weighting");
    MatrixXd ld(2, 2);
    ld << 1.0, 0.6, 0.6, 2.0;
    const auto correlated = skatOConservativeGrid(scores, ld, vector({1.0, 2.0}), {1.0});
    close(correlated.p_value, std::erfc(0.5 / std::sqrt(2.0 * 2.85)), 1e-13, "Burden uses off-diagonal covariance");
    const auto singular = skatOConservativeGrid(vector({3.0, 3.0}), MatrixXd::Ones(2, 2), weights, {1.0});
    close(singular.p_value, std::erfc(3.0 / std::sqrt(2.0)), 1e-13, "Perfect LD burden limit");
    MatrixXd opposing_ld(2, 2);
    opposing_ld << 1.0, -1.0, -1.0, 1.0;
    check(skatOConservativeGrid(vector({1.0, -1.0 + 1e-10}), opposing_ld, weights, {1.0}).p_value == 1.0,
          "Roundoff score residual must not create a significant zero-variance burden");
    check(skatOConservativeGrid(VectorXd(), MatrixXd(0, 0), VectorXd()).p_value == 1.0, "Empty grid set");
    check(skatOConservativeGrid(VectorXd::Zero(2), MatrixXd::Zero(2, 2), weights).n_variants == 0,
          "Zero covariance produces neutral result");
    invalid([&] { skatOConservativeGrid(scores, covariance, weights, {}); }, "Reject empty rho grid");
    invalid([&] { skatOConservativeGrid(scores, covariance, weights, {NAN}); }, "Reject nonfinite rho");
    invalid([&] { skatOConservativeGrid(scores, covariance, weights, {1.01}); }, "Reject rho out of range");
    invalid([&] { skatOConservativeGrid(scores, covariance, VectorXd::Zero(2)); }, "Reject zero weights");
    invalid([&] { skatOConservativeGrid(scores, covariance, vector({1.0})); }, "Reject dimensions");
    invalid([&] { skatOConservativeGrid(vector({NAN, 0.0}), covariance, weights); }, "Reject NaN scores");
    invalid([&] { skatOConservativeGrid(scores, -covariance, weights); }, "Reject indefinite covariance");
    invalid([&] { skatOConservativeGrid(scores, MatrixXd::Zero(2, 2), weights); }, "Reject inconsistent zero covariance");
    invalid([&] { skatOConservativeGrid(scores, MatrixXd::Ones(2, 2), weights); }, "Reject score outside covariance support");
    ld(0, 1) = 0.7;
    invalid([&] { skatOConservativeGrid(scores, ld, weights); }, "Reject asymmetric covariance");
    ld(0, 1) = INFINITY;
    invalid([&] { skatOConservativeGrid(scores, ld, weights); }, "Reject infinite covariance");
}

void testAcat() {
    for (double p : {1e-300, 1e-100, 1e-12, 0.01, 0.2, 0.5, 0.8, std::nextafter(1.0, 0.0)}) {
        close(acatV(vector({p}), vector({7.0})).p_value, p, 1e-14, "Single-variant ACAT identity");
        close(acatV(vector({p, p}), vector({1.0, 3.0})).p_value, p, 1e-14, "Identical p-value identity");
    }
    const double subnormal = std::numeric_limits<double>::denorm_min();
    check(acatV(vector({subnormal}), vector({1.0})).p_value == subnormal, "Subnormal ACAT identity");
    close(acatV(vector({1e-300, 0.5}), vector({1.0, 3.0})).p_value, 4e-300, 1e-14, "Tiny weighted ACAT tail");
    close(acatV(vector({0.25, 0.75}), vector({1.0, 1.0})).p_value, 0.5, 1e-14, "Cauchy cancellation");
    const double expected = std::atan(2.0) / std::acos(-1.0);
    close(acatV(vector({0.25, 0.75}), vector({3.0, 1.0})).p_value, expected, 1e-14, "Weighted analytic Cauchy sum");
    close(acatV(vector({0.25, 0.75}), vector({3e300, 1e300})).p_value, expected, 1e-14, "Huge weight normalization");
    close(acatV(vector({0.0, 1.0, 0.25}), vector({0.0, 0.0, 1.0})).p_value, 0.25, 1e-14, "Ignore zero-weight endpoints");
    check(acatV(vector({0.0, 0.2}), vector({1.0, 1.0})).p_value == 0.0, "Exact p=0 boundary");
    check(acatV(vector({1.0, 0.2}), vector({1.0, 1.0})).p_value == 1.0, "Exact p=1 boundary");
    check(acatV(VectorXd(), VectorXd()).p_value == 1.0, "Empty ACAT set");
    check(acatVFromScores(VectorXd(), MatrixXd(0, 0), VectorXd()).p_value == 1.0, "Empty ACAT score set");
    check(acatVFromScores(VectorXd::Zero(2), MatrixXd::Zero(2, 2), VectorXd::Ones(2)).p_value == 1.0,
          "Neutral all-uninformative ACAT score set");
    invalid([] { acatV(vector({0.0, 1.0}), vector({1.0, 1.0})); }, "Ambiguous exact endpoints rejected");
    invalid([] { acatV(vector({NAN}), vector({1.0})); }, "Reject NaN p");
    invalid([] { acatV(vector({INFINITY}), vector({1.0})); }, "Reject infinite p");
    invalid([] { acatV(vector({-0.1}), vector({1.0})); }, "Reject negative p");
    invalid([] { acatV(vector({1.1}), vector({1.0})); }, "Reject p>1");
    invalid([] { acatV(vector({0.1}), vector({-1.0})); }, "Reject negative weights");
    invalid([] { acatV(vector({0.1}), vector({INFINITY})); }, "Reject infinite weights");
    invalid([] { acatV(vector({0.1}), vector({NAN})); }, "Reject NaN weights");
    invalid([] { acatV(vector({0.1}), vector({0.0})); }, "Reject zero total weight");
    invalid([] { acatV(vector({0.1}), VectorXd()); }, "Reject mismatched weights");
    MatrixXd covariance = MatrixXd::Zero(2, 2);
    covariance(0, 0) = 4.0;
    const auto from_scores = acatVFromScores(vector({6.0, 0.0}), covariance, vector({1.0, 1.0}));
    check(from_scores.n_variants == 1, "Omit zero-variance variants");
    close(from_scores.p_value, std::erfc(3.0 / std::sqrt(2.0)), 1e-13, "ACAT score covariance scaling");
    check(acatVFromScores(vector({50.0}), MatrixXd::Identity(1, 1), vector({1.0})).p_value > 0.0,
          "Score p-value underflow is not an exact endpoint");
    invalid([&] { acatVFromScores(vector({6.0, 1.0}), covariance, vector({1.0, 1.0})); },
            "Reject nonzero zero-variance score");
}

void testFrameworkAndNull() {
    const MatrixXd x = MatrixXd::Identity(2, 2);
    const VectorXd scores = vector({3.0, -1.0});
    const VectorXd maf = vector({0.01, 0.1});
    const VectorXd beta_weights = (1.0 - maf.array()).pow(24.0);
    const VectorXd equal_maf = vector({0.1, 0.1});
    cosmic::RareVariantTest burden_test(cosmic::RareVariantTest::Method::BURDEN);
    const auto burden_wrapped = burden_test.runTest("gene", x, scores, x, equal_maf);
    close(burden_wrapped.stat, 2.0, 1e-14, "Framework burden chi-square oracle");
    close(burden_wrapped.p_value, std::erfc(1.0), 1e-14,
          "Framework burden normal-tail oracle");
    cosmic::RareVariantTest skat_test(cosmic::RareVariantTest::Method::SKAT);
    const auto skat_wrapped = skat_test.runTest("gene", x, scores, x, equal_maf);
    const double equal_weight = std::pow(0.9, 24.0);
    const double expected_skat_q = 10.0 * equal_weight * equal_weight;
    close(skat_wrapped.stat, expected_skat_q, 1e-14, "Framework SKAT statistic oracle");
    close(skat_wrapped.p_value, std::exp(-5.0), 1e-12,
          "Framework SKAT equal-spectrum chi-square oracle");
    cosmic::RareVariantTest grid_test(cosmic::RareVariantTest::Method::SKAT_O);
    const auto wrapped = grid_test.runTest("gene", x, scores, x, maf);
    const auto direct = skatOCanonical(scores, x, beta_weights);
    close(wrapped.p_value, direct.p_value, 1e-14, "Framework dispatches canonical grid omnibus");
    check(wrapped.method_name == kSkatOCanonicalMethod && wrapped.set_id == "gene" &&
          wrapped.p_method == "DEPENDENCE_ADJUSTED_GAUSSIAN_QMC" && !wrapped.fallback_used,
          "Canonical SKAT-O method metadata");
    const auto one = skatOCanonical(vector({3.0}), MatrixXd::Identity(1, 1), vector({1.0}));
    close(one.p_value, std::erfc(3.0 / std::sqrt(2.0)), 1e-13,
          "Single-variant canonical SKAT-O equals its score test");
    const auto one_rho = skatOCanonical(scores, x, beta_weights, {0.0});
    close(one_rho.p_value, one_rho.rho_p_values[0], 1e-14,
          "Single-rho canonical SKAT-O equals its marginal test");
    const auto analytic_two_rho = skatOCanonical(scores, x, VectorXd::Ones(2), {0.0, 1.0});
    // Independent polar-normal quadrature gives 0.010658284823866517 for
    // P(X^2+Y^2 >= 10 or |(X+Y)/sqrt(2)| >= 2.709525773141792).
    check(std::abs(analytic_two_rho.p_value - 0.010658284823866517) < 3e-4,
          "Canonical SKAT-O differs from independent two-dimensional quadrature oracle");
    std::cout << std::setprecision(15)
              << "BURDEN_ORACLE statistic=" << burden_wrapped.stat
              << " expected=2 p=" << burden_wrapped.p_value
              << " expected_p=" << std::erfc(1.0) << '\n'
              << "SKAT_ORACLE statistic=" << skat_wrapped.stat
              << " expected=" << expected_skat_q
              << " p=" << skat_wrapped.p_value
              << " expected_p=" << std::exp(-5.0) << '\n'
              << "SKAT_O_ORACLE p=" << analytic_two_rho.p_value
              << " expected_p=0.0106582848238665 tolerance=0.0003\n";
    check(direct.p_value >= direct.min_p_value &&
          direct.p_value <= std::min(1.0, direct.min_p_value * direct.rho_grid.size()),
          "Dependence adjustment lies between min-p and Bonferroni");
    cosmic::RareVariantTest acat_test(cosmic::RareVariantTest::Method::ACAT_V);
    const auto acat = acat_test.runTest("gene", x, scores, x, maf);
    close(acat.p_value, acatVFromScores(scores, x,
          beta_weights.array().square() * maf.array() * (1.0 - maf.array())).p_value,
          1e-14, "Framework ACAT weighting convention");
    check(acat.method_name == kAcatVMethod, "ACAT metadata");
    check(grid_test.runTest("empty", MatrixXd(2, 0), VectorXd(), x, VectorXd()).p_value == 1.0,
          "Framework empty set");
    invalid([&] { grid_test.runTest("bad", x, scores, x, vector({0.1, NAN})); }, "Framework rejects invalid MAF");
    // Regression guard only, not a claim of empirical validation of canonical
    // SKAT-O: the {0,1} grid has exact marginal tails and a union-bound correction.
    std::mt19937_64 rng(20260908);
    std::normal_distribution<double> normal;
    int rejected = 0;
    for (int i = 0; i < 4000; ++i) {
        const VectorXd u = vector({normal(rng), normal(rng)});
        rejected += skatOConservativeGrid(u, x, VectorXd::Ones(2), {0.0, 1.0}).p_value <= 0.05;
    }
    check(rejected <= 250, "Gaussian null grid must not resemble unadjusted min-p calibration");
}

} // namespace

int main() {
    try {
        testMixture();
        testGrid();
        testAcat();
        testFrameworkAndNull();
        std::cout << "Rare-variant extension tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
