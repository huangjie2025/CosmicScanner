#include "association/generalized/firth_logistic.h"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cosmic;
using namespace cosmic::association;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
} // namespace

int main() {
    constexpr int n = 60;
    Eigen::MatrixXd phenotype(n, 1);
    Eigen::MatrixXd covariates = Eigen::MatrixXd::Ones(n, 1);
    Eigen::MatrixXd genotypes(n, 2);
    std::vector<Eigen::Index> identity(n);
    for (int i = 0; i < n; ++i) {
        identity[i] = i;
        phenotype(i, 0) = i >= n / 2 ? 1.0 : 0.0;
        genotypes(i, 0) = phenotype(i, 0) * 2.0;
        genotypes(i, 1) = (i * 7 + i / 4) % 3;
    }

    const ScanContext context{phenotype, covariates, identity};
    const GenotypeBlock block{genotypes, {"separated", "ordinary"}};
    FirthLogisticTest test;
    test.prepare(context);
    const auto generic = test.scanBlock(block);
    check(generic.size() == 2 && test.results().size() == 2, "Missing Firth results");
    const auto& separated = test.results()[0];
    check(separated.converged && separated.beta > 0.0 && std::isfinite(separated.beta) &&
          std::isfinite(separated.se) && separated.se > 0.0 &&
          std::isfinite(separated.p_value) && separated.p_value >= 0.0 && separated.p_value <= 1.0,
          "Firth failed to produce finite inference under complete separation");
    check(test.getModelType() == AssociationModel::FIRTH_LOGISTIC &&
          test.capabilities().block_scan, "Firth adapter contract mismatch");

    Eigen::MatrixXd invalid = phenotype;
    invalid(0, 0) = 2.0;
    try {
        FirthLogisticTest invalid_test;
        invalid_test.prepare(ScanContext{invalid, covariates, identity});
    } catch (const std::invalid_argument&) {
        return 0;
    }
    throw std::runtime_error("Firth accepted a non-binary phenotype");
}
