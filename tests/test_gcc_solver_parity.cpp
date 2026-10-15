// Phase E3-A unit test: GccSolver parity against dense inverse.
//
// Verifies that the new GccSolver abstraction produces identical results
// to the old explicit inverse approach.
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#include <string>

#include "gcc_solver.h"
#include "apy_ginv.h"

using namespace cosmic;

static int g_failures = 0;

#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) {                                                            \
            ++g_failures;                                                         \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, msg);             \
        }                                                                         \
    } while (0)

static Eigen::MatrixXd makeSpdMatrix(int n, unsigned seed) {
    std::mt19937 rng(seed);
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
    return A * A.transpose() + Eigen::MatrixXd::Identity(n, n) * 0.1;
}

int main() {
    try {
        std::printf("=== GccSolver Parity Test ===\n");

        // Test 1: solve(I) == inverse
        {
            const int n = 50;
            Eigen::MatrixXd Gcc = makeSpdMatrix(n, 12345);

            GccSolverDense solver(Gcc);
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
            Eigen::MatrixXd inv_solver = solver.solve(I);

            Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);
            Eigen::MatrixXd inv_ldlt = ldlt.solve(I);

            double diff = (inv_solver - inv_ldlt).cwiseAbs().maxCoeff();
            std::printf("Test 1: solve(I) vs LDLT.solve(I) max|diff| = %.3e\n", diff);
            CHECK(diff < 1e-12, "solve(I) must match LDLT.solve(I)");
        }

        // Test 2: solve(random vectors)
        {
            const int n = 50;
            Eigen::MatrixXd Gcc = makeSpdMatrix(n, 12345);
            GccSolverDense solver(Gcc);

            Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);

            std::mt19937 rng(67890);
            std::normal_distribution<> dist(0.0, 1.0);

            double max_diff = 0.0;
            for (int i = 0; i < 100; ++i) {
                Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
                for (int j = 0; j < n; ++j) b(j) = dist(rng);

                Eigen::VectorXd x_solver = solver.solve(b);
                Eigen::VectorXd x_ldlt = ldlt.solve(b);

                double diff = (x_solver - x_ldlt).cwiseAbs().maxCoeff();
                max_diff = std::max(max_diff, diff);
            }
            std::printf("Test 2: solve(random) vs LDLT.solve(random) max|diff| = %.3e\n", max_diff);
            CHECK(max_diff < 1e-12, "solve(random) must match LDLT.solve(random)");
        }

        // Test 3: solve(Gcn) - multiple RHS
        {
            const int n_core = 30;
            const int n_noncore = 15;
            Eigen::MatrixXd Gcc = makeSpdMatrix(n_core, 111);
            Eigen::MatrixXd Gcn = Eigen::MatrixXd::Random(n_core, n_noncore);

            GccSolverDense solver(Gcc);
            Eigen::MatrixXd Q_solver = solver.solve(Gcn);

            Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);
            Eigen::MatrixXd Q_ldlt = ldlt.solve(Gcn);

            double diff = (Q_solver - Q_ldlt).cwiseAbs().maxCoeff();
            std::printf("Test 3: solve(Gcn) max|diff| = %.3e\n", diff);
            CHECK(diff < 1e-12, "solve(Gcn) must match LDLT.solve(Gcn)");
        }

        // Test 4: ApyGinv build + applyInverse parity (full core, n_core = n)
        {
            const int ng = 20;
            Eigen::MatrixXd G = makeSpdMatrix(ng, 222);
            std::vector<int> core_idx(ng);
            for (int i = 0; i < ng; ++i) core_idx[i] = i;

            ApyGinv g_old;
            g_old.build(G, core_idx, "old");

            ApyGinv g_new;
            g_new.build(G, core_idx, "new");

            std::mt19937 rng(333);
            std::normal_distribution<> dist(0.0, 1.0);
            Eigen::VectorXd x(ng);
            for (int i = 0; i < ng; ++i) x(i) = dist(rng);

            Eigen::VectorXd y_old = g_old.applyInverse(x);
            Eigen::VectorXd y_new = g_new.applyInverse(x);

            double diff = (y_old - y_new).cwiseAbs().maxCoeff();
            std::printf("Test 4a: applyInverse (core=n) max|diff| = %.3e\n", diff);
            CHECK(diff < 1e-12, "applyInverse must match between old/new");

            Eigen::VectorXd x_core(ng);
            for (int i = 0; i < ng; ++i) x_core(i) = dist(rng);
            Eigen::VectorXd yc_old = g_old.applyCoreInverse(x_core);
            Eigen::VectorXd yc_new = g_new.applyCoreInverse(x_core);

            double diff_c = (yc_old - yc_new).cwiseAbs().maxCoeff();
            std::printf("Test 4b: applyCoreInverse (core=n) max|diff| = %.3e\n", diff_c);
            CHECK(diff_c < 1e-12, "applyCoreInverse must match between old/new");
        }

        // Test 5: ApyGinv build + applyInverse parity (core < n)
        {
            const int ng = 20;
            const int n_core = 12;
            Eigen::MatrixXd G = makeSpdMatrix(ng, 444);
            std::vector<int> core_idx(n_core);
            for (int i = 0; i < n_core; ++i) core_idx[i] = i;

            ApyGinv g_old;
            g_old.build(G, core_idx, "old");

            ApyGinv g_new;
            g_new.build(G, core_idx, "new");

            std::mt19937 rng(555);
            std::normal_distribution<> dist(0.0, 1.0);
            Eigen::VectorXd x(ng);
            for (int i = 0; i < ng; ++i) x(i) = dist(rng);

            Eigen::VectorXd y_old = g_old.applyInverse(x);
            Eigen::VectorXd y_new = g_new.applyInverse(x);

            double diff = (y_old - y_new).cwiseAbs().maxCoeff();
            std::printf("Test 5a: applyInverse (core<n) max|diff| = %.3e\n", diff);
            CHECK(diff < 1e-12, "applyInverse (core<n) must match");

            Eigen::VectorXd x_core(n_core);
            for (int i = 0; i < n_core; ++i) x_core(i) = dist(rng);
            Eigen::VectorXd yc_old = g_old.applyCoreInverse(x_core);
            Eigen::VectorXd yc_new = g_new.applyCoreInverse(x_core);

            double diff_c = (yc_old - yc_new).cwiseAbs().maxCoeff();
            std::printf("Test 5b: applyCoreInverse (core<n) max|diff| = %.3e\n", diff_c);
            CHECK(diff_c < 1e-12, "applyCoreInverse (core<n) must match");

            // Check Q matrix
            double diff_q = (g_old.q() - g_new.q()).cwiseAbs().maxCoeff();
            std::printf("Test 5c: Q matrix max|diff| = %.3e\n", diff_q);
            CHECK(diff_q < 1e-12, "Q matrix must match");

            // Check mnn
            double diff_mnn = (g_old.mnn() - g_new.mnn()).cwiseAbs().maxCoeff();
            std::printf("Test 5d: mnn vector max|diff| = %.3e\n", diff_mnn);
            CHECK(diff_mnn < 1e-12, "mnn vector must match");

            // Check core_block
            double diff_cb = (g_old.coreBlock() - g_new.coreBlock()).cwiseAbs().maxCoeff();
            std::printf("Test 5e: core_block max|diff| = %.3e\n", diff_cb);
            CHECK(diff_cb < 1e-12, "core_block must match");
        }

        // Test 6: gccInvDiagonal and gccInvFull
        {
            const int n_core = 15;
            Eigen::MatrixXd Gcc = makeSpdMatrix(n_core, 777);

            GccSolverDense solver(Gcc);
            Eigen::VectorXd diag_solver = solver.diagonal();
            Eigen::MatrixXd inv_full = solver.solve(Eigen::MatrixXd::Identity(n_core, n_core).eval());
            Eigen::VectorXd diag_full = inv_full.diagonal();

            double diff = (diag_solver - diag_full).cwiseAbs().maxCoeff();
            std::printf("Test 6: diagonal() vs inv.diagonal() max|diff| = %.3e\n", diff);
            CHECK(diff < 1e-12, "diagonal() must match full inverse diagonal");

            // Test gccInvFull() through ApyGinv
            std::vector<int> core_idx(n_core);
            for (int i = 0; i < n_core; ++i) core_idx[i] = i;
            ApyGinv g;
            g.build(Gcc, core_idx, "diag_test");

            Eigen::VectorXd diag_ap = g.gccInvDiagonal();
            Eigen::MatrixXd inv_ap = g.gccInvFull();

            double diff_diag = (diag_solver - diag_ap).cwiseAbs().maxCoeff();
            double diff_inv = (inv_full - inv_ap).cwiseAbs().maxCoeff();
            std::printf("Test 6b: ApyGinv gccInvDiagonal() max|diff| = %.3e\n", diff_diag);
            std::printf("Test 6c: ApyGinv gccInvFull() max|diff| = %.3e\n", diff_inv);
            CHECK(diff_diag < 1e-12, "gccInvDiagonal() must match solver diagonal");
            CHECK(diff_inv < 1e-12, "gccInvFull() must match solver full inverse");
        }

        // Test 7: LDLT pivot diagnostics
        {
            const int n_core = 10;
            Eigen::MatrixXd Gcc = makeSpdMatrix(n_core, 888);
            GccSolverDense solver(Gcc);
            ApyGinv g;
            std::vector<int> core_idx(n_core);
            for (int i = 0; i < n_core; ++i) core_idx[i] = i;
            g.build(Gcc, core_idx, "pivot_test");

            CHECK(solver.minDiagD() > 0.0, "minDiagD must be positive");
            CHECK(solver.maxDiagD() > 0.0, "maxDiagD must be positive");
            CHECK(std::fabs(solver.minDiagD() - g.info().gcc_min_diag_d) < 1e-12,
                  "minDiagD must match ApyGinv info");
            CHECK(std::fabs(solver.maxDiagD() - g.info().gcc_max_diag_d) < 1e-12,
                  "maxDiagD must match ApyGinv info");
            std::printf("Test 7: LDLT pivot diagnostics match\n");
        }

        // Test 8: ApplyInverse with core < n uses correct formula
        // G_APY^-1 = [Gcc^-1 + Q Mnn^-1 Q'  -Q Mnn^-1; -Mnn^-1 Q'  Mnn^-1]
        {
            const int ng = 16;
            const int n_core = 10;
            Eigen::MatrixXd G = makeSpdMatrix(ng, 999);
            std::vector<int> core_idx(n_core);
            for (int i = 0; i < n_core; ++i) core_idx[i] = i;

            ApyGinv g;
            g.build(G, core_idx, "formula_test");

            // Extract submatrices
            Eigen::MatrixXd Gcc(n_core, n_core);
            for (int i = 0; i < n_core; ++i)
                for (int j = 0; j < n_core; ++j)
                    Gcc(i, j) = G(core_idx[i], core_idx[j]);

            Eigen::MatrixXd Gcn(n_core, ng - n_core);
            for (int i = 0; i < n_core; ++i)
                for (int j = 0; j < ng - n_core; ++j)
                    Gcn(i, j) = G(core_idx[i], j + n_core);

            std::vector<int> noncore_idx;
            for (int i = n_core; i < ng; ++i) noncore_idx.push_back(i);

            Eigen::MatrixXd Gnn(ng - n_core, ng - n_core);
            for (int i = 0; i < ng - n_core; ++i)
                for (int j = 0; j < ng - n_core; ++j)
                    Gnn(i, j) = G(noncore_idx[i], noncore_idx[j]);

            // Manual APY inverse construction
            Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);
            Eigen::MatrixXd Gcc_inv = ldlt.solve(Eigen::MatrixXd::Identity(n_core, n_core));
            Eigen::MatrixXd Q = ldlt.solve(Gcn);
            Eigen::VectorXd mnn_diag(ng - n_core);
            for (int j = 0; j < ng - n_core; ++j) {
                double m = Gnn(j, j) - Gcn.col(j).dot(Q.col(j));
                mnn_diag(j) = m;
            }
            Eigen::VectorXd mnn_inv = mnn_diag.cwiseInverse();
            Eigen::MatrixXd core_block_manual = Gcc_inv + Q * mnn_inv.asDiagonal() * Q.transpose();

            // Compare with ApyGinv core_block_
            double diff_cb = (g.coreBlock() - core_block_manual).cwiseAbs().maxCoeff();
            std::printf("Test 8: core_block vs manual APY formula max|diff| = %.3e\n", diff_cb);
            CHECK(diff_cb < 1e-12, "core_block must match manual APY formula");
        }

    } catch (const std::exception& e) {
        ++g_failures;
        std::printf("FATAL: %s\n", e.what());
    }

    if (g_failures == 0) {
        std::printf("\nALL GccSolver PARITY TESTS PASSED\n");
    } else {
        std::printf("\nFAILURES: %d\n", g_failures);
    }
    return g_failures == 0 ? 0 : 1;
}