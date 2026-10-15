#include "association/generalized/logistic.h"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cosmic;
using namespace cosmic::association;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    constexpr int n = 80;
    Eigen::MatrixXd phenotype(n, 1), covariates = Eigen::MatrixXd::Ones(n, 1);
    Eigen::MatrixXd genotypes(n, 2);
    std::vector<Eigen::Index> identity(n);
    for (int i = 0; i < n; ++i) {
        phenotype(i, 0) = i % 2;
        genotypes(i, 0) = (i * 7 + i / 5) % 3;
        genotypes(i, 1) = phenotype(i, 0);
        identity[i] = i;
    }
    const Eigen::VectorXd marker = genotypes.col(0);
    const Eigen::VectorXd adjusted = marker.array() - marker.mean();
    const double score = adjusted.dot((phenotype.col(0).array() - 0.5).matrix());
    const double variance = 0.25 * adjusted.squaredNorm();
    const double expected_p = std::erfc(std::abs(score) / std::sqrt(2.0 * variance));

    LogisticTest score_test;
    score_test.prepare(ScanContext{phenotype, covariates, identity});
    auto score_results = score_test.scanBlock(GenotypeBlock{genotypes.leftCols(1), {"marker"}});
    check(score_results[0].status == AssociationStatus::SUCCESS, "Logistic score test failed");
    check(std::abs(score_results[0].p_value - expected_p) < 1e-12,
          "Logistic score differs from analytic reference");

    LogisticOptions wald_options;
    wald_options.test = LogisticTestType::WALD;
    LogisticTest wald_test(wald_options);
    wald_test.prepare(ScanContext{phenotype, covariates, identity});
    auto wald_results = wald_test.scanBlock(GenotypeBlock{genotypes, {"marker", "separated"}});
    check(wald_results[0].status == AssociationStatus::SUCCESS &&
          std::isfinite(wald_results[0].p_value), "Logistic Wald test failed");
    check(wald_results[1].status == AssociationStatus::NUMERICAL_FAILURE &&
          wald_test.results()[1].separated, "Logistic Wald did not explicitly detect separation");

    phenotype(0, 0) = 2.0;
    bool rejected = false;
    try { score_test.prepare(ScanContext{phenotype, covariates, identity}); }
    catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Logistic accepted a non-binary phenotype");
}
