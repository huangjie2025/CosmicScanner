#include "association/generalized/saige_style.h"
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
} // namespace

int main() {
    constexpr int n = 80;
    Eigen::VectorXd y(n), marker(n);
    Eigen::MatrixXd X = Eigen::MatrixXd::Ones(n, 1);
    Eigen::SparseMatrix<double> kinship(n, n);
    kinship.setIdentity();
    for (int i = 0; i < n; ++i) {
        y[i] = i % 2;
        marker[i] = (i * 7 + i / 5) % 3;
    }

    SaigeStyleNullModel null_model;
    null_model.fit(y, X, kinship);
    null_model.calibrateVarianceRatio({marker});
    const SaigeStyleResult result = null_model.testMarker(marker);
    const Eigen::VectorXd adjusted = marker.array() - marker.mean();
    const double expected_score = adjusted.dot((y.array() - 0.5).matrix());
    const double expected_variance = 0.25 * adjusted.squaredNorm();
    const double expected_p = std::erfc(std::abs(expected_score) /
                                       std::sqrt(2.0 * expected_variance));
    check(result.tested && std::isfinite(result.beta) && std::isfinite(result.se) &&
          result.se > 0.0 && std::isfinite(result.p_value),
          "SAIGE-style score test did not produce finite inference");
    check(std::abs(result.score - expected_score) < 1e-8 &&
          std::abs(result.score_variance - expected_variance) < 1e-8 &&
          std::abs(result.p_value - expected_p) < 1e-10,
          "SAIGE-style unrelated score test differs from analytic logistic score reference");
    check(null_model.metadata().converged && null_model.metadata().sparse_nnz == n &&
          std::abs(null_model.tau()) < 1e-12 &&
          std::abs(null_model.varianceRatio() - 1.0) < 1e-10,
          "SAIGE-style null-model metadata or variance ratio is invalid");

    const auto exact = null_model.testMarkerExact(marker);
    LogisticTest logistic;
    Eigen::MatrixXd phenotype = y;
    std::vector<Eigen::Index> identity(n);
    for (int i = 0; i < n; ++i) identity[i] = i;
    logistic.prepare(ScanContext{phenotype, X, identity});
    Eigen::MatrixXd genotypes(n, 1);
    genotypes.col(0) = marker;
    const auto logistic_result = logistic.scanBlock(GenotypeBlock{genotypes, {"marker"}})[0];
    check(exact.tested && std::abs(null_model.tau()) < 1e-12 &&
          std::abs(exact.p_value - logistic_result.p_value) < 1e-10 &&
          std::abs(exact.beta - logistic_result.effects[0]) < 1e-10,
          "Tau-zero logistic GLMM does not match ordinary logistic score inference");

    SaigeStyleTest adapter(kinship);
    adapter.prepare(ScanContext{phenotype, X, identity});
    GenotypeBlock block{genotypes, {"marker"}};
    const auto generic = adapter.scanBlock(block);
    check(generic.size() == 1 && generic[0].status == AssociationStatus::SUCCESS,
          "SAIGE-style AssociationTest adapter failed");

    // Deterministic 1:9 case-control simulation with a removable relatedness component.
    constexpr int ni = 120;
    Eigen::VectorXd yi = Eigen::VectorXd::Zero(ni);
    Eigen::MatrixXd Xi = Eigen::MatrixXd::Ones(ni, 1);
    Eigen::SparseMatrix<double> full(ni, ni), loco(ni, ni);
    full.setIdentity();
    loco.setIdentity();
    for (int pair = 0; pair < 6; ++pair) {
        const int a = 2 * pair, b = a + 1;
        yi[a] = yi[b] = 1.0;
        full.coeffRef(a, b) = full.coeffRef(b, a) = 0.25;
    }
    Eigen::VectorXd rare = Eigen::VectorXd::Zero(ni);
    rare[0] = rare[2] = rare[4] = rare[6] = 1.0;
    std::vector<Eigen::VectorXd> calibration{rare};
    for (int j = 0; j < 12; ++j) {
        Eigen::VectorXd candidate(ni);
        for (int i = 0; i < ni; ++i) candidate[i] = (i * (j + 3) + j * j) % 3;
        calibration.push_back(std::move(candidate));
    }
    SaigeStyleNullModel imbalanced({100, 1e-5, 10.0});
    imbalanced.fit(yi, Xi, full);
    imbalanced.calibrateVarianceRatio(calibration);
    const auto approximate_rare = imbalanced.testMarker(rare);
    const auto exact_rare = imbalanced.testMarkerExact(rare);
    check(imbalanced.metadata().converged && imbalanced.tau() > 0.0 &&
          approximate_rare.tested && exact_rare.tested &&
          std::abs(approximate_rare.p_value - exact_rare.p_value) < 1e-10,
          "Imbalanced SAIGE-style MAC calibration differs from exact mixed score variance");

    SaigeStyleNullModel loco_null({100, 1e-5, 10.0});
    loco_null.fit(yi, Xi, loco);
    loco_null.calibrateVarianceRatio(calibration);
    const auto loco_rare = loco_null.testMarker(rare);
    check(loco_null.metadata().converged && loco_rare.tested &&
          std::isfinite(loco_rare.p_value) && std::abs(loco_null.tau()) < 1e-12,
          "SAIGE-style LOCO null did not remove the chromosome-specific relatedness component");

    // 1:10 imbalance test
    constexpr int n10 = 110;
    Eigen::VectorXd y10 = Eigen::VectorXd::Zero(n10);
    Eigen::MatrixXd X10 = Eigen::MatrixXd::Ones(n10, 1);
    Eigen::SparseMatrix<double> kinship10(n10, n10);
    kinship10.setIdentity();
    for (int i = 0; i < 10; ++i) y10[i] = 1.0;
    Eigen::VectorXd marker10(n10);
    for (int i = 0; i < n10; ++i) marker10[i] = (i * 7 + i / 5) % 3;
    SaigeStyleNullModel null10({100, 1e-5, 10.0});
    null10.fit(y10, X10, kinship10);
    null10.calibrateVarianceRatio({marker10});
    const auto r10 = null10.testMarker(marker10);
    check(null10.metadata().converged && std::isfinite(null10.tau()) && null10.tau() >= 0.0 &&
          null10.varianceRatio() > 0.0 && std::isfinite(null10.varianceRatio()) &&
          r10.tested && std::isfinite(r10.beta) && std::isfinite(r10.se) && r10.se > 0.0 &&
          std::isfinite(r10.p_value) && r10.p_value >= 0.0 && r10.p_value <= 1.0 &&
          (!r10.spa_used || (r10.p_value >= 0.0 && r10.p_value <= 1.0)),
          "SAIGE 1:10 imbalance failed");

    // 1:50 imbalance test
    constexpr int n50 = 102;
    Eigen::VectorXd y50 = Eigen::VectorXd::Zero(n50);
    Eigen::MatrixXd X50 = Eigen::MatrixXd::Ones(n50, 1);
    Eigen::SparseMatrix<double> kinship50(n50, n50);
    kinship50.setIdentity();
    y50[0] = 1.0;
    y50[1] = 1.0;
    Eigen::VectorXd marker50(n50);
    for (int i = 0; i < n50; ++i) marker50[i] = (i * 7 + i / 5) % 3;
    SaigeStyleNullModel null50({100, 1e-5, 10.0});
    null50.fit(y50, X50, kinship50);
    null50.calibrateVarianceRatio({marker50});
    const auto r50 = null50.testMarker(marker50);
    check(null50.metadata().converged && std::isfinite(null50.tau()) && null50.tau() >= 0.0 &&
          null50.varianceRatio() > 0.0 && std::isfinite(null50.varianceRatio()) &&
          r50.tested && std::isfinite(r50.beta) && std::isfinite(r50.se) && r50.se > 0.0 &&
          std::isfinite(r50.p_value) && r50.p_value >= 0.0 && r50.p_value <= 1.0 &&
          (!r50.spa_used || (r50.p_value >= 0.0 && r50.p_value <= 1.0)),
          "SAIGE 1:50 imbalance failed");
}
