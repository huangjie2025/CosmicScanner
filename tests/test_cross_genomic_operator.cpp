// Phase E3-D2D unit test: matrix-free CrossGenomicOperator (Gcn = Zc Zn' / S)
// parity with the block-streamed (materialized) Gcn across Raw / Tuned / Final
// (identity and A22) views, transposition consistency, the APY MME / Eq. 8 PEV
// operator applies with the cross operator active, and the no-regen gate
// (the per-apply PCG applies never regenerate the full K x Nn Gcn block).
//
// Self-contained: generates a deterministic synthetic PLINK dataset.
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

#include "plink_reader.h"
#include "models/block_streamed_grm.h"
#include "models/cross_genomic_operator.h"
#include "models/apy_ginv.h"
#include "models/apy_inverse_operator.h"
#include "models/apy_core_pev_operator.h"
#include "models/gcc_solver.h"

using namespace cosmic;

static int g_failures = 0;

#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) {                                                            \
            ++g_failures;                                                         \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, msg);             \
        }                                                                         \
    } while (0)

// --- deterministic synthetic PLINK writer (mirrors test_e1_operator) ---------
static void writeBed(FILE* f, const std::vector<uint8_t>& dosages, int n, int n_snps) {
    const unsigned char magic[3] = {0x6c, 0x1b, 0x01};
    fwrite(magic, 1, 3, f);
    for (int j = 0; j < n_snps; ++j) {
        std::vector<unsigned char> bytes((n + 3) / 4, 0);
        for (int i = 0; i < n; ++i) {
            uint8_t code = 1;  // missing
            if (dosages[(size_t)j * n + i] < 250) {
                switch (dosages[(size_t)j * n + i]) {
                    case 0: code = 3; break;  // 11 -> 0.0
                    case 1: code = 2; break;  // 10 -> 1.0
                    case 2: code = 0; break;  // 00 -> 2.0
                }
            }
            bytes[(size_t)i / 4] |= (code << (2 * (i % 4)));
        }
        fwrite(bytes.data(), 1, bytes.size(), f);
    }
}

static std::string writeSyntheticData(const std::filesystem::path& dir, int n, int m,
                                      unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<uint8_t> dosages((size_t)m * n, 250);  // 250 = missing marker
    for (int j = 0; j < m; ++j) {
        const double p = 0.10 + 0.40 * ((rng() % 1000) / 1000.0);
        const double p2 = p * p;
        const double h = 2.0 * p * (1.0 - p);
        for (int i = 0; i < n; ++i) {
            const double uu = u(rng);
            if (uu < 0.03) continue;  // missing
            if (uu - 0.03 < p2) dosages[(size_t)j * n + i] = 2;
            else if (uu - 0.03 < p2 + h) dosages[(size_t)j * n + i] = 1;
            else dosages[(size_t)j * n + i] = 0;
        }
    }
    std::string prefix = (dir / "e3d2syn").string();
    {
        FILE* f = std::fopen((prefix + ".bed").c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write bed");
        writeBed(f, dosages, n, m);
        std::fclose(f);
    }
    {
        std::ofstream f(prefix + ".bim");
        for (int j = 0; j < m; ++j) f << "1\tsnp_" << j << "\t0\t" << (j + 1) << "\tA\tC\n";
    }
    {
        std::ofstream f(prefix + ".fam");
        for (int i = 0; i < n; ++i) f << "0\tind_" << i << "\t0\t0\t0\t-9\n";
    }
    return prefix;
}

// Simple SPD "pedigree-like" A22 target (as in test_e1_operator).

static Eigen::MatrixXd reconstruct(const CrossGenomicOperator& op, bool cols_to_core) {
    // Reconstruct the K x Nn matrix Gcn by applying to basis vectors.
    const int nc = op.coreRows();
    const int nn = op.noncoreRows();
    Eigen::MatrixXd G(nc, nn);
    if (cols_to_core) {
        for (int b = 0; b < nn; ++b) {
            Eigen::VectorXd e = Eigen::VectorXd::Zero(nn); e(b) = 1.0;
            G.col(b) = op.applyCoreFromNoncore(e);
        }
    } else {
        for (int a = 0; a < nc; ++a) {
            Eigen::VectorXd e = Eigen::VectorXd::Zero(nc); e(a) = 1.0;
            G.row(a) = op.applyNoncoreFromCore(e).transpose();
        }
    }
    return G;
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path data_root =
            std::filesystem::temp_directory_path() / "cosmic_cross_operator_data";
        std::filesystem::remove_all(data_root);
        std::filesystem::create_directories(data_root);

        const int N = 64, M = 200;   // M > N keeps the APY conditional variances positive
        const int NC = 22;  // core size
        const int NN = N - NC;
        std::string prefix =
            argc > 1 ? std::string(argv[1]) : writeSyntheticData(data_root, N, M, 20260925u);

        PlinkReader reader;
        reader.load(prefix);
        reader.computeStats(1);
        if ((int)reader.getNumSamples() != N || (int)reader.getNumSnps() != M) {
            throw std::runtime_error("test data dims mismatch");
        }

        std::vector<int> keep(N);
        for (int i = 0; i < N; ++i) keep[i] = i;
        std::vector<int> core_idx(NC), noncore_idx(NN);
        for (int a = 0; a < NC; ++a) core_idx[a] = a;
        for (int b = 0; b < NN; ++b) noncore_idx[b] = NC + b;

        BlockStreamedGrmBuilder::Options bopts;
        bopts.min_maf = 0.01f;
        bopts.n_threads = 1;
        bopts.materialize_gcn = true;

        const double tune_a = 0.03, tune_b = 0.91;
        const double alpha = 0.05;
        // A22 target of the streamed-fixture shape: identity plus a small
        // constant relationship baseline (proven-positive-mnn configuration).
        const Eigen::MatrixXd A22 =
            Eigen::MatrixXd::Identity(N, N).array() + 0.01;

        // Reference builders: identity-blend and A22-blend materialized Gcn.
        BlockStreamedGrmBuilder refI(reader, keep, bopts);
        refI.computeStats();
        const double S = refI.S2pq();
        refI.setTune(tune_a, tune_b);
        refI.computeBlocks(core_idx, noncore_idx);
        refI.finalize(alpha, nullptr);
        const Eigen::MatrixXd gcn_raw_ref = refI.gcnRaw();
        const Eigen::MatrixXd gcn_tuned_ref = refI.gcnTuned();
        const Eigen::MatrixXd gcn_finalI_ref = refI.gcnFinal();
        const Eigen::MatrixXd gcc_final_ref = refI.gccFinal();
        const Eigen::VectorXd gnn_diag_final_ref = refI.gnnDiagFinal();
        {
            Eigen::LDLT<Eigen::MatrixXd> ldlt(gcc_final_ref);
            const Eigen::MatrixXd Gccinv = ldlt.solve(Eigen::MatrixXd::Identity(NC, NC));
            const Eigen::MatrixXd Q = Gccinv * gcn_finalI_ref;
            const Eigen::VectorXd mnn_probe =
                gnn_diag_final_ref -
                ((Q.array() * gcn_finalI_ref.array()).colwise().sum().transpose()).matrix();
            std::printf("  [probe] mnn min=%.6e max=%.6e mean=%.6e\n",
                        mnn_probe.minCoeff(), mnn_probe.maxCoeff(), mnn_probe.mean());
        }

        BlockStreamedGrmBuilder refA(reader, keep, bopts);
        refA.computeStats();
        refA.setTune(tune_a, tune_b);
        refA.computeBlocks(core_idx, noncore_idx);
        refA.finalize(alpha, &A22);
        const Eigen::MatrixXd gcn_finalA_ref = refA.gcnFinal();

        // Non-materialized builders (streamed, genotype-backed Gcn) for the
        // operator reference path and the no-regen gate.
        BlockStreamedGrmBuilder::Options nopts;
        nopts.min_maf = 0.01f;
        nopts.n_threads = 1;
        nopts.materialize_gcn = false;

        BlockStreamedGrmBuilder streamI(reader, keep, nopts);
        streamI.computeStats();
        streamI.setTune(tune_a, tune_b);
        streamI.computeBlocks(core_idx, noncore_idx);
        streamI.finalize(alpha, nullptr);

        BlockStreamedGrmBuilder nbuilder(reader, keep, nopts);
        nbuilder.computeStats();
        nbuilder.setTune(tune_a, tune_b);
        nbuilder.computeBlocks(core_idx, noncore_idx);
        nbuilder.finalize(alpha, &A22);

        // keep arrays are SAMPLE indices (identity here); pass the genotyped-order
        // index arrays for the A22 slice.
        const std::vector<int> core_keep(core_idx.begin(), core_idx.end());
        const std::vector<int> noncore_keep(noncore_idx.begin(), noncore_idx.end());

        CrossGenomicOperator op_raw(reader, core_keep, noncore_keep, S, 0.0, 1.0,
                                    0.0, GcnBlockView::Raw, nullptr,
                                    core_idx, noncore_idx, 1);
        CrossGenomicOperator op_tuned(reader, core_keep, noncore_keep, S, tune_a, tune_b,
                                      0.0, GcnBlockView::Tuned, nullptr,
                                      core_idx, noncore_idx, 1);
        CrossGenomicOperator op_finalI(reader, core_keep, noncore_keep, S, tune_a, tune_b,
                                       alpha, GcnBlockView::Final, nullptr,
                                       core_idx, noncore_idx, 1);
        CrossGenomicOperator op_finalA(reader, core_keep, noncore_keep, S, tune_a, tune_b,
                                       alpha, GcnBlockView::Final, &A22,
                                       core_idx, noncore_idx, 1);
        CrossGenomicOperator op_effect_tuned(reader, core_keep, noncore_keep, S, tune_a,
                                             tune_b, 0.0, GcnBlockView::Tuned, nullptr,
                                             core_idx, noncore_idx, 1);

        std::printf(
            "E3-D2D cross operator parity (n=%d, m=%d, nc=%d, nn=%d, S=%.6f, a=%.6f, b=%.6f, alpha=%.4f)\n",
            N, M, NC, NN, S, tune_a, tune_b, alpha);

        auto report = [](const char* label, const Eigen::MatrixXd& opG,
                         const Eigen::MatrixXd& ref) {
            const Eigen::MatrixXd D = opG - ref;
            const double max_abs = D.cwiseAbs().maxCoeff();
            const double rel_frob = ref.norm() > 0 ? D.norm() / ref.norm() : 0.0;
            std::printf("  [%s] max|D|=%.3e frob_rel=%.3e\n", label, max_abs, rel_frob);
            return std::make_pair(max_abs, rel_frob);
        };

        // --- Cross-block parity (Raw / Tuned / Final-identity / Final-A22) ---
        {
            const double tol = 1e-5;
            auto [e1, r1] = report("raw", reconstruct(op_raw, true), gcn_raw_ref);
            CHECK(e1 < tol, "raw Gcn cross-operator parity exceeded (1e-5)");
            CHECK(r1 < 1e-4, "raw Gcn frob-rel exceeded (1e-4)");
            auto [e2, r2] = report("tuned", reconstruct(op_tuned, true), gcn_tuned_ref);
            CHECK(e2 < tol, "tuned Gcn cross-operator parity exceeded");
            auto [e3, r3] =
                report("final_I", reconstruct(op_finalI, true), gcn_finalI_ref);
            CHECK(e3 < tol, "final(identity-blend) Gcn cross-operator parity exceeded");
            auto [e4, r4] =
                report("final_A22", reconstruct(op_finalA, true), gcn_finalA_ref);
            CHECK(e4 < tol, "final(A22-blend) Gcn cross-operator parity exceeded");
            // Transpose consistency: the Gcn' path (applyNoncoreFromCore on
            // core basis vectors) must reproduce the same Gcn as the forward
            // path up to floating-point associativity.
            auto [e5, r5] = report("Gcn'_vs_Gcn", reconstruct(op_raw, true),
                                   reconstruct(op_raw, false));
            CHECK(e5 < 1e-12, "raw Gcn'/Gcn reconstruction paths must agree (1e-12)");
            // Effect operator (Tuned view) parity with the tuned reference.
            auto [e6, r6] = report("effect_tuned", reconstruct(op_effect_tuned, true),
                                   gcn_tuned_ref);
            CHECK(e6 < tol, "tuned (effect) Gcn cross-operator parity exceeded");
        }

        // --- Operator linearity + symmetry (raw) ---
        {
            std::mt19937 rng(99);
            Eigen::VectorXd xn(NN), xc(NC);
            for (int i = 0; i < NN; ++i) xn(i) = ((rng() % 1000) / 1000.0) - 0.3;
            for (int i = 0; i < NC; ++i) xc(i) = ((rng() % 1000) / 1000.0) + 0.2;
            const Eigen::VectorXd g1 = op_raw.applyCoreFromNoncore(xn + xn);
            const Eigen::VectorXd g2 = op_raw.applyCoreFromNoncore(2.0 * xn);
            CHECK((g1 - g2).cwiseAbs().maxCoeff() < 1e-6,
                  "applyCoreFromNoncore must be homogeneous");
            const double ym = xc.dot(op_raw.applyCoreFromNoncore(xn));
            const double mt = xn.dot(op_raw.applyNoncoreFromCore(xc));
            std::printf("  symmetry: x_c'(Gcn x_n)=%.10f vs x_n'(Gcn'x_c)=%.10f (diff %.3e)\n",
                        ym, mt, std::abs(ym - mt));
            CHECK(std::abs(ym - mt) < 1e-6,
                  "applyCoreFromNoncore/applyNoncoreFromCore must be adjoints");
        }

        // --- generateGcnBlock (non-materialized, float conventions) parity ---
        {
            GcnBlock blk;
            nbuilder.makeGcnIterator(-1, GcnBlockView::Tuned)->next(blk);
            const Eigen::MatrixXd cross_tuned = reconstruct(op_tuned, true);
            const Eigen::MatrixXd D(blk.block - cross_tuned);
            const double max_abs = D.cwiseAbs().maxCoeff();
            std::printf("  [generateGcnBlock(tuned) vs cross] max|D|=%.3e\n", max_abs);
            CHECK(max_abs < 1e-5,
                  "generateGcnBlock(tuned) vs cross operator parity exceeded (1e-5)");
        }

        // --- APY MME operator apply: cross vs streamed block iterator --------
        // Production pattern (single_step_model.cpp): the reference factory is
        // the genotype-backed Gcn iterator of a non-materialized builder.
        const Eigen::MatrixXd mme_gcc = streamI.gccFinal();
        const Eigen::VectorXd& mme_gnn_diag = streamI.gnnDiagFinal();
        GcnBlockIteratorFactory stream_I_factory = [&streamI]() {
            return streamI.makeGcnIterator(-1, GcnBlockView::Final);
        };
        MatrixFreeApyInverseOperator mme_ref(makeDenseSolver(mme_gcc), stream_I_factory,
                                             mme_gnn_diag, N, core_idx,
                                             "E3D2D MME ref");
        MatrixFreeApyInverseOperator mme_cross(makeDenseSolver(mme_gcc), stream_I_factory,
                                               mme_gnn_diag, N, core_idx,
                                               "E3D2D MME cross", &op_finalI);
        ApyGinv mme_dense;
        mme_dense.build(refI.gccFinal(), refI.gcnFinal(), refI.gnnDiagFinal(), N,
                        core_idx, "E3D2D MME dense");
        DenseApyInverseOperator mme_adapter(mme_dense);
        {
            std::mt19937 rng(20260925u);
            Eigen::VectorXd x(N);
            for (int i = 0; i < N; ++i) x(i) = ((rng() % 10000) / 10000.0) - 0.5;
            const Eigen::VectorXd y_ref = mme_ref.apply(x);
            const Eigen::VectorXd y_cross = mme_cross.apply(x);
            const Eigen::VectorXd y_dense = mme_adapter.apply(x);
            const double rel_cross =
                y_ref.norm() > 0 ? (y_cross - y_ref).norm() / y_ref.norm()
                                 : (y_cross - y_ref).norm();
            const double rel_dense =
                y_ref.norm() > 0 ? (y_dense - y_ref).norm() / y_ref.norm()
                                 : (y_dense - y_ref).norm();
            std::printf("  MME apply: cross-vs-ref rel=%.3e ; blockiter-vs-dense rel=%.3e\n",
                        rel_cross, rel_dense);
            CHECK(rel_dense < 1e-6, "block-iterator MME apply must match dense APY inverse");
            CHECK(rel_cross < 1e-4, "cross-operator MME apply parity exceeded (1e-4)");
        }

        // --- Eq. 8 PEV operator applyCoefficient: cross vs block iterator ----
        MatrixFreePEVOperator::Options pev_opts;
        pev_opts.tolerance = 1e-8;
        pev_opts.max_iterations = 2000;
        const double lambda = 0.5;
        Eigen::VectorXd Dc = Eigen::VectorXd::Ones(NC);
        Eigen::VectorXd Dn = Eigen::VectorXd::Ones(NN);
        const Eigen::VectorXd& mnn = mme_dense.mnn();
        GccSolverPtr gcc_solver = makeDenseSolver(mme_gcc);
        MatrixFreeAPYCorePEVOperator pev_ref(gcc_solver, stream_I_factory, mnn, Dc, Dn,
                                             lambda, pev_opts, nullptr);
        MatrixFreeAPYCorePEVOperator pev_cross(gcc_solver, stream_I_factory, mnn, Dc, Dn,
                                               lambda, pev_opts, &op_finalI);
        {
            std::mt19937 rng(11);
            for (int trial = 0; trial < 5; ++trial) {
                Eigen::VectorXd x(NC);
                for (int i = 0; i < NC; ++i) x(i) = ((rng() % 10000) / 10000.0) - 0.5;
                Eigen::VectorXd y_ref(NC), y_cross(NC);
                pev_ref.applyCoefficient(x, y_ref);
                pev_cross.applyCoefficient(x, y_cross);
                const double rel =
                    y_ref.norm() > 0 ? (y_cross - y_ref).norm() / y_ref.norm()
                                     : (y_cross - y_ref).norm();
                std::printf("  PEV applyCoefficient trial %d: cross-vs-ref rel=%.3e\n",
                            trial, rel);
                CHECK(rel < 1e-4, "PEV applyCoefficient cross parity exceeded (1e-4)");
            }
            // Full PCG solve parity (each does its own PCG iteration; results
            // differ only by convergence rounding).
            Eigen::VectorXd x(NC), b(NC);
            std::mt19937 rng2(13);
            for (int i = 0; i < NC; ++i) { x(i) = ((rng2() % 10000) / 10000.0) - 0.5; b(i) = x(i); }
            Eigen::VectorXd q_ref = pev_ref.apply(b);
            Eigen::VectorXd q_cross = pev_cross.apply(b);
            const double rel =
                q_ref.norm() > 0 ? (q_cross - q_ref).norm() / q_ref.norm()
                                 : (q_cross - q_ref).norm();
            std::printf("  PEV full apply: cross-vs-ref rel=%.3e (ref iters=%d, cross iters=%d)\n",
                        rel, pev_ref.lastSolveInfo().iterations,
                        pev_cross.lastSolveInfo().iterations);
            CHECK(rel < 1e-3, "PEV full apply cross parity exceeded (1e-3)");
        }

        // --- No-regen gate: applies never regenerate the K x Nn Gcn block ---
        {
            GcnBlockIteratorFactory stream_factory = [&nbuilder]() {
                return nbuilder.makeGcnIterator(-1, GcnBlockView::Final);
            };
            const Eigen::VectorXd gnn_diag_stream = nbuilder.gnnDiagFinal();
            MatrixFreeApyInverseOperator mme_stream(
                makeDenseSolver(nbuilder.gccFinal()), stream_factory,
                gnn_diag_stream, N, core_idx, "E3D2D MME stream", &op_finalA);
            const long long regen_mme_after_construct = nbuilder.gcnRegenCalls();
            std::mt19937 rng(17);
            for (int t = 0; t < 8; ++t) {
                Eigen::VectorXd x(N);
                for (int i = 0; i < N; ++i) x(i) = ((rng() % 10000) / 10000.0) - 0.5;
                (void)mme_stream.apply(x);
            }
            CHECK(nbuilder.gcnRegenCalls() == regen_mme_after_construct,
                  "MME cross applies must not regenerate the Gcn block");

            MatrixFreeAPYCorePEVOperator pev_stream(
                makeDenseSolver(nbuilder.gccFinal()), stream_factory, mnn, Dc, Dn,
                lambda, pev_opts, &op_finalA);
            const long long regen_pev_after_construct = nbuilder.gcnRegenCalls();
            for (int t = 0; t < 8; ++t) {
                Eigen::VectorXd x(NC);
                for (int i = 0; i < NC; ++i) x(i) = ((rng() % 10000) / 10000.0) - 0.5;
                (void)pev_stream.apply(x);
            }
            std::printf(
                "  no-regen gate: regen calls MME=%lld PEV=%lld (baseline %lld)\n",
                regen_mme_after_construct, regen_pev_after_construct,
                nbuilder.gcnRegenCalls());
            CHECK(nbuilder.gcnRegenCalls() == regen_pev_after_construct,
                  "PEV cross applies must not regenerate the Gcn block");
        }

    } catch (const std::exception& e) {
        ++g_failures;
        std::printf("FATAL: %s\n", e.what());
    }

    if (g_failures == 0) {
        std::printf("E3-D2D CROSS OPERATOR PARITY: PASS\n");
    } else {
        std::printf("E3-D2D CROSS OPERATOR PARITY: FAIL (%d failures)\n", g_failures);
    }
    return g_failures == 0 ? 0 : 1;
}