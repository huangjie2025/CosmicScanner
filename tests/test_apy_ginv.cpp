// Phase C2 unit tests: ApyGinv (Leite et al. 2024 GSE 56:59, Eqs. 5-6).
// Verifies the APY block-inverse operator against explicit dense reference
// computations on small random SPD matrices.
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>
#include "models/apy_ginv.h"

using namespace cosmic;

static int g_failures = 0;

#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) {                                                            \
            ++g_failures;                                                         \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, msg);             \
        }                                                                         \
    } while (0)

static Eigen::MatrixXd randomSpdG(int ng, int m, unsigned seed) {
    std::mt19937 rng(seed);
    Eigen::MatrixXd M(ng, m);
    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < m; ++j) M(i, j) = (rng() % 3) - 1.0;  // -1/0/1 dosages
    Eigen::MatrixXd G = (M * M.transpose()) / (double)m;
    // Stabilize: dominant diagonal (like a blended/tuned GRM).
    double mean_diag = G.diagonal().mean();
    G.diagonal().array() += 0.2 * mean_diag;
    return G;
}

static void testInverseIdentity(int ng, const std::vector<int>& core) {
    const int m = 200;
    Eigen::MatrixXd G = randomSpdG(ng, m, 42);
    ApyGinv apy;
    apy.build(G, core, "test");

    // Materialize G_APY^-1 from the triplets.
    Eigen::MatrixXd Ginv_full = Eigen::MatrixXd::Zero(ng, ng);
    apy.visitTriplets([&](int r, int c, double v) { Ginv_full(r, c) = v; });

    // Symmetry.
    double max_asym = (Ginv_full - Ginv_full.transpose()).cwiseAbs().maxCoeff();
    CHECK(max_asym < 1e-12, "G_APY^-1 must be symmetric");

    // Build the matrix whose EXACT inverse is the APY operator (Leite Eq. 5):
    //   G_apy = G with the noncore block replaced by Mnn + Q' Gcn
    //         = G with noncore off-diagonals replaced by the conditional
    //           prediction given the core (Q = Gcc^-1 Gcn, Mnn = diag(m_jj)).
    Eigen::MatrixXd Gcc = Eigen::MatrixXd::Zero((int)core.size(), (int)core.size());
    for (int a = 0; a < (int)core.size(); ++a)
        for (int b = 0; b < (int)core.size(); ++b) Gcc(a, b) = G(core[a], core[b]);
    auto gcc_ldlt = Gcc.ldlt();
    std::vector<int> noncore;
    for (int i = 0; i < ng; ++i)
        if (std::find(core.begin(), core.end(), i) == core.end()) noncore.push_back(i);
    const int nn = (int)noncore.size();
    Eigen::MatrixXd Gcn((int)core.size(), nn);
    for (int a = 0; a < (int)core.size(); ++a)
        for (int j = 0; j < nn; ++j) Gcn(a, j) = G(core[a], noncore[j]);
    Eigen::MatrixXd Q = gcc_ldlt.solve(Gcn);
    Eigen::MatrixXd Mnn = Eigen::MatrixXd::Zero(nn, nn);
    Eigen::MatrixXd G_apy = G;
    for (int j = 0; j < nn; ++j) {
        // m_jj = Gnn(j,j) - Gnc(j,:) Gcc^-1 Gcn(:,j)  (noncore conditional variance)
        Eigen::VectorXd g_cj((int)core.size());
        for (int a = 0; a < (int)core.size(); ++a) g_cj(a) = G(core[a], noncore[j]);
        double m_jj = G(noncore[j], noncore[j]) - g_cj.dot(gcc_ldlt.solve(g_cj));
        // m_jj must equal 1 / G_APY^-1(noncore_j, noncore_j) (the Mnn^-1 block).
        CHECK(std::fabs(m_jj - 1.0 / Ginv_full(noncore[j], noncore[j])) < 1e-10,
              "Mnn diagonal consistent");
        Mnn(j, j) = m_jj;
    }
    for (int a = 0; a < nn; ++a)
        for (int b = 0; b < nn; ++b) G_apy(noncore[a], noncore[b]) = Mnn(a, b) + (Q.transpose() * Gcn)(a, b);

    // G_apy * G_APY^-1 ~= I  (the APY operator is the EXACT inverse of G_apy).
    Eigen::MatrixXd prod = G_apy * Ginv_full;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(ng, ng);
    CHECK((prod - I).cwiseAbs().maxCoeff() < 1e-8, "G_apy * G_APY^-1 ~= I");

    // applyInverse must agree with the materialized matrix.
    std::mt19937 rng(7);
    for (int rep = 0; rep < 5; ++rep) {
        Eigen::VectorXd x(ng);
        for (int i = 0; i < ng; ++i) x(i) = (double)(rng() % 1000) / 100.0 - 5.0;
        Eigen::VectorXd y1 = apy.applyInverse(x);
        Eigen::VectorXd y2 = Ginv_full * x;
        CHECK((y1 - y2).cwiseAbs().maxCoeff() < 1e-10, "applyInverse matches materialized inverse");
    }

    // coreBlock() must equal the core rows/cols of G_APY^-1.
    Eigen::MatrixXd cb = apy.coreBlock();
    for (int a = 0; a < (int)core.size(); ++a)
        for (int b = 0; b < (int)core.size(); ++b)
            CHECK(std::fabs(cb(a, b) - Ginv_full(core[a], core[b])) < 1e-12,
                  "coreBlock equals core slice of G_APY^-1");

    // applyCoreInverse = Gcc^-1 on core vectors.
    for (int rep = 0; rep < 3; ++rep) {
        Eigen::VectorXd xc((int)core.size());
        for (int i = 0; i < (int)core.size(); ++i) xc(i) = (double)(rng() % 1000) / 100.0 - 5.0;
        Eigen::VectorXd yc = apy.applyCoreInverse(xc);
        // dense Gcc^-1 from the same Gcc above (recompute here for independence).
        Eigen::MatrixXd Gcc = Eigen::MatrixXd::Zero((int)core.size(), (int)core.size());
        for (int a = 0; a < (int)core.size(); ++a)
            for (int b = 0; b < (int)core.size(); ++b) Gcc(a, b) = G(core[a], core[b]);
        Eigen::VectorXd yref = Gcc.ldlt().solve(xc);
        CHECK((yc - yref).cwiseAbs().maxCoeff() < 1e-10, "applyCoreInverse = Gcc^-1 x_c");
    }
}

static void testCoreIsAll(int ng) {
    Eigen::MatrixXd G = randomSpdG(ng, 200, 3);
    std::vector<int> core(ng);
    for (int i = 0; i < ng; ++i) core[i] = i;

    ApyGinv apy;
    apy.build(G, core, "core=all");
    CHECK(apy.info().core_is_all, "core_is_all flag");
    CHECK(apy.nNoncore() == 0, "no noncore when core=n_g");

    Eigen::MatrixXd Ginv_full = Eigen::MatrixXd::Zero(ng, ng);
    apy.visitTriplets([&](int r, int c, double v) { Ginv_full(r, c) = v; });
    Eigen::MatrixXd Ginv_ref = G.inverse();
    CHECK((Ginv_full - Ginv_ref).cwiseAbs().maxCoeff() < 1e-9,
          "core=n_g degenerates to the exact G^-1 (C2 gate G2)");

    std::mt19937 rng(11);
    Eigen::VectorXd x(ng);
    for (int i = 0; i < ng; ++i) x(i) = (double)(rng() % 1000) / 100.0 - 5.0;
    Eigen::VectorXd y1 = apy.applyInverse(x);
    Eigen::VectorXd y2 = Ginv_ref * x;
    CHECK((y1 - y2).cwiseAbs().maxCoeff() < 1e-9, "applyInverse = G^-1 at core=n_g");
}

static void testInvalidInputs() {
    Eigen::MatrixXd G = randomSpdG(20, 100, 5);
    ApyGinv apy;
    // Empty core.
    bool threw = false;
    try { apy.build(G, {}, "empty"); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "empty core throws");
    // Out-of-range core index.
    threw = false;
    try { apy.build(G, {0, 25}, "oor"); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "out-of-range core throws");
    // Duplicate core indices.
    threw = false;
    try { apy.build(G, {3, 3, 5}, "dup"); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "duplicate core throws");
    // Non-SPD Gcc (singular Gcc: identical core rows AND columns -> symmetric).
    threw = false;
    try {
        Eigen::MatrixXd Gbad = randomSpdG(10, 50, 9);
        Gbad.row(1) = Gbad.row(2);
        Gbad.col(1) = Gbad.col(2);  // symmetric singular core block {1,2}
        ApyGinv bad;
        bad.build(Gbad, {1, 2}, "singular");
    } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "singular Gcc throws");
    // Non-square G.
    threw = false;
    try {
        Eigen::MatrixXd Grect = Eigen::MatrixXd::Random(10, 11);
        ApyGinv bad;
        bad.build(Grect, {0, 1}, "rect");
    } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw, "non-square G throws");
}

static void testOrderInvariance() {
    Eigen::MatrixXd G = randomSpdG(30, 150, 13);
    ApyGinv a, b;
    a.build(G, {2, 7, 11, 19}, "a");
    b.build(G, {7, 19, 2, 11}, "b");  // same set, different order
    std::mt19937 rng(21);
    Eigen::VectorXd x(30);
    for (int i = 0; i < 30; ++i) x(i) = (double)(rng() % 1000) / 100.0 - 5.0;
    Eigen::VectorXd ya = a.applyInverse(x);
    Eigen::VectorXd yb = b.applyInverse(x);
    CHECK((ya - yb).cwiseAbs().maxCoeff() < 1e-12, "order-insensitive partition");
}

static void testDiagonal() {
    Eigen::MatrixXd G = randomSpdG(20, 100, 8);
    ApyGinv apy;
    apy.build(G, {1, 4, 9, 15}, "diag");
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(20);
    apy.addDiagonalTo(diag, 1.0);
    Eigen::MatrixXd Ginv_full = Eigen::MatrixXd::Zero(20, 20);
    apy.visitTriplets([&](int r, int c, double v) { Ginv_full(r, c) = v; });
    CHECK((diag - Ginv_full.diagonal()).cwiseAbs().maxCoeff() < 1e-14,
          "addDiagonalTo matches materialized diagonal");
}

int main() {
    testInverseIdentity(40, {0, 3, 7, 11, 19, 22, 31, 38});
    testInverseIdentity(25, {2, 5, 8, 13, 21});
    testInverseIdentity(20, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    testCoreIsAll(35);
    testCoreIsAll(8);
    testInvalidInputs();
    testOrderInvariance();
    testDiagonal();
    if (g_failures == 0) {
        std::printf("test_apy_ginv: ALL PASS\n");
        return 0;
    }
    std::printf("test_apy_ginv: %d FAILURE(S)\n", g_failures);
    return 1;
}