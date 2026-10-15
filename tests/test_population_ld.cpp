#include "../src/genetics/population/pca.h"
#include "../src/genetics/ld/ld_calculator.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cosmic::genetics;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void rejects(Function function, const char* message) {
    try {
        function();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

void testPca() {
    MatrixXd grm(3, 3);
    grm << 2, -1, 0, -1, 2, 0, 0, 0, 0;
    const MatrixXd original = grm;
    const auto result = PCA::computeFromGrm(grm, 3);
    VectorXd expected(3);
    expected << 3, 1, 0;
    check(result.eigenvalues.isApprox(expected, 1e-12), "PCA descending eigenvalues");
    check((result.scores.transpose() * result.scores).isApprox(MatrixXd::Identity(3, 3)),
          "PCA scores must be unit eigenvectors");
    check((grm * result.scores).isApprox(result.scores * result.eigenvalues.asDiagonal()),
          "PCA eigen equation");
    for (int j = 0; j < 3; ++j) {
        Eigen::Index pivot;
        result.scores.col(j).cwiseAbs().maxCoeff(&pivot);
        check(result.scores(pivot, j) > 0, "PCA deterministic sign");
    }
    check(result.scores.isApprox(PCA::computeFromGrm(grm, 3).scores), "PCA repeatability");
    check(PCA::computeFromGrm(grm, 1).scores.cols() == 1, "PCA k truncation");
    check(grm == original, "PCA input mutation");
    check(PCA::computeFromGrm(MatrixXd::Zero(2, 2), 2).eigenvalues.isZero(), "zero GRM");
    check(PCA::computeFromGrm(MatrixXd::Identity(2, 2), 2).scores.allFinite(), "tied eigenvalues");
    rejects([&] { PCA::computeFromGrm(grm, 0); }, "PCA k=0 accepted");
    rejects([&] { PCA::computeFromGrm(grm, 4); }, "PCA k>n accepted");
    rejects([] { PCA::computeFromGrm(MatrixXd(0, 0), 1); }, "empty PCA accepted");
    rejects([] { PCA::computeFromGrm(MatrixXd::Zero(2, 3), 1); }, "nonsquare PCA accepted");
    grm(0, 1) = 0;
    rejects([&] { PCA::computeFromGrm(grm, 1); }, "asymmetric PCA accepted");
    grm = original;
    grm(2, 2) = -0.1;
    rejects([&] { PCA::computeFromGrm(grm, 1); }, "negative discarded eigenvalue accepted");
    grm *= 1e-100;
    rejects([&] { PCA::computeFromGrm(grm, 1); }, "tiny indefinite GRM accepted");
    grm = original;
    grm(2, 2) = -1e-12;
    check(PCA::computeFromGrm(grm, 3).eigenvalues[2] == 0, "PSD roundoff not clamped");
    grm(0, 0) = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { PCA::computeFromGrm(grm, 1); }, "NaN PCA accepted");
    grm(0, 0) = std::numeric_limits<double>::infinity();
    rejects([&] { PCA::computeFromGrm(grm, 1); }, "infinite PCA accepted");
}

void testLd() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    MatrixXd g(4, 3);
    g << 0, 2, 0, 1, 1, 2, 2, 0, 0, nan, -9, 2;
    MatrixXd expected(3, 3);
    expected << 1, -1, 0, -1, 1, 0, 0, 0, 1;
    check(LDCalculator::computePairwiseLD(g).isApprox(expected, 1e-12), "signed imputed LD");
    check(std::isnan(g(3, 0)) && g(3, 1) == -9, "LD input mutation");
    g(3, 0) = 1;
    g(3, 1) = 1;
    check(LDCalculator::computePairwiseLD(g).isApprox(expected, 1e-12), "explicit mean imputation");
    g.col(1) << 0, 0.5, 1, 1.5;
    check(std::abs(LDCalculator::computePairwiseLD(g)(0, 1) - std::sqrt(0.4)) < 1e-12,
          "fractional dosage correlation");
    check(LDCalculator::computePairwiseLD(MatrixXd(2, 0)).size() == 0, "empty LD columns");
    rejects([] { LDCalculator::computePairwiseLD(MatrixXd::Zero(1, 2)); }, "one sample accepted");
    for (double value : {0.1, nan, -9.0, 3.0, std::numeric_limits<double>::infinity(),
                         -std::numeric_limits<double>::infinity()}) {
        rejects([&] { LDCalculator::computePairwiseLD(MatrixXd::Constant(4, 1, value)); },
                "invalid LD column accepted");
    }
}

void testWindows() {
    MatrixXd g(4, 5);
    g.col(0) << 0, 1, 2, 1;
    g.col(1) = 2.0 - g.col(0).array();
    g.col(2) << 0, 2, 0, 2;
    g.col(3) = g.col(0);
    g.col(4) = g.col(0);
    std::vector<std::string> chrom{"1", "1", "1", "1", "2"};
    std::vector<long long> pos{100, 110, 120, 1000, 100};
    VectorXd p(5);
    p << 0.02, 0.01, 0.9, 0.8, 0.001;
    check(LDPruner::prune(g, chrom, pos, 20, 1, 1.0) == std::vector<int>({0, 2, 3, 4}),
          "pruning boundaries, sign, chromosome isolation");
    check(LDPruner::prune(g, chrom, pos, 9, 1, 0.5) == std::vector<int>({0, 1, 2, 3, 4}),
          "pruning physical window");
    check(LDClumper::clump(g, chrom, pos, 20, 1, 1.0, p, 0.01) ==
          std::vector<int>({1, 1, -1, -1, 4}), "clump lead index, threshold, local LD");
    check(LDClumper::clump(g, chrom, pos, 20, 3, 1.0, p, 0.01) ==
          std::vector<int>({1, 1, -1, -1, 4}), "clump batch invariance");
    p[0] = p[1];
    check(LDClumper::clump(g, chrom, pos, 20, 1, 1.0, p, 0.01)[1] == 0, "clump p tie");
    p.setConstant(0.9);
    check(LDClumper::clump(g, chrom, pos, 20, 1, 0.5, p, 0.01) ==
          std::vector<int>(5, -1), "clump no significant leads");
    MatrixXd chain(4, 3);
    chain.col(0) = g.col(0);
    chain.col(1) = g.col(2);
    chain.col(2) = g.col(2);
    check(LDPruner::prune(chain, {"1", "1", "1"}, {0, 10, 20}, 10, 1, 0.5) ==
          std::vector<int>({0, 1}), "overlapping pruning windows");
    check(LDPruner::prune(chain, {"1", "1", "1"}, {0, 10, 20}, 10, 2, 0.5) ==
          std::vector<int>({0, 1, 2}), "pruning step ignored");
    const long long last = std::numeric_limits<long long>::max();
    check(LDPruner::prune(g.leftCols(2), {"1", "1"}, {last - 10, last}, 10, 1, 1) ==
          std::vector<int>({0}), "large position overflow");
    check(LDPruner::prune(g.leftCols(2), {"1", "1"}, {10, 10}, 0, 1, 1) ==
          std::vector<int>({0}), "duplicate position zero window");
    check(LDPruner::prune(MatrixXd(2, 0), {}, {}, 0, 1, 0.5).empty(), "empty prune");
    check(LDClumper::clump(MatrixXd(2, 0), {}, {}, 0, 1, 0.5, VectorXd(0), 0.05).empty(),
          "empty clump");
    rejects([&] { LDPruner::prune(g, {}, pos, 20, 1, 0.5); }, "metadata mismatch accepted");
    rejects([&] { LDPruner::prune(g, chrom, pos, -1, 1, 0.5); }, "negative window accepted");
    rejects([&] { LDPruner::prune(g, chrom, pos, 20, 0, 0.5); }, "zero step accepted");
    for (double r2 : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
        rejects([&] { LDPruner::prune(g, chrom, pos, 20, 1, r2); }, "invalid r2 accepted");
    }
    pos[1] = 99;
    rejects([&] { LDPruner::prune(g, chrom, pos, 20, 1, 0.5); }, "unsorted positions accepted");
    pos[1] = 110;
    chrom[1] = "2";
    rejects([&] { LDClumper::clump(g, chrom, pos, 20, 1, 0.5, p, 0.05); }, "unsorted chromosomes accepted");
    chrom[1] = "1";
    rejects([&] { LDClumper::clump(g, chrom, pos, 20, 1, 0.5, VectorXd(0), 0.05); }, "p mismatch accepted");
    for (double value : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
        rejects([&] { LDClumper::clump(g, chrom, pos, 20, 1, 0.5, p, value); }, "invalid p threshold accepted");
        p[0] = value;
        rejects([&] { LDClumper::clump(g, chrom, pos, 20, 1, 0.5, p, 0.05); }, "invalid p accepted");
        p[0] = 0.9;
    }
    g.col(0).setConstant(1);
    rejects([&] { LDPruner::prune(g, chrom, pos, 20, 1, 0.5); }, "prune monomorphic accepted");
    rejects([&] { LDClumper::clump(g, chrom, pos, 20, 1, 0.5, p, 0.05); }, "clump monomorphic accepted");

    // Widely spaced variants would require >3 GB for an accidental dense LD matrix.
    const int m = 20000;
    MatrixXd sparse(2, m);
    sparse.row(0).setZero();
    sparse.row(1).setConstant(2);
    chrom.assign(m, "1");
    pos.resize(m);
    for (int j = 0; j < m; ++j) pos[j] = 1000LL * j;
    check(LDPruner::prune(sparse, chrom, pos, 10, 1, 0.5).size() == m, "bounded prune");
    const auto leads = LDClumper::clump(sparse, chrom, pos, 10, 1, 0.5,
                                      VectorXd::Constant(m, 0.01), 0.05);
    for (int j = 0; j < m; ++j) check(leads[j] == j, "bounded clump");
}
} // namespace

int main() {
    try {
        testPca();
        testLd();
        testWindows();
        std::cout << "population/LD tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "population/LD test failure: " << error.what() << '\n';
        return 1;
    }
}
