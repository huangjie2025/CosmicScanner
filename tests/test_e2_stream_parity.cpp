// Phase E2 unit test: block-streamed VanRaden G construction parity.
//
// The streamed path (BlockStreamedGrmBuilder) must reproduce exactly the APY-
// relevant blocks of the production dense GRM (GrmBuilder, algorithm 1):
//   GCC_PARITY, GCN_PARITY, GNN_DIAG_PARITY  (raw, tune, identity and A22 blends)
//   S_PARITY                                 (sum 2p(1-p))
//   SCALAR_SUM_PARITY                        (raw/final diag & off-diag means)
//   APY_OPERATOR_PARITY                      (ApyGinv::build(G) == build(blocks):
//                                             q_, mnn_, core_block_, applyInverse*)
//   CORE_N_EQUALS_NG_DEGENERACY              (core = n_g reproduces the full
//                                             inverse build)
//
// Self-contained: generates a deterministic synthetic PLINK dataset. When a
// dataset prefix is given as argv[1] (e.g. benchmarks/.../sim), it is used
// instead (freeze-dataset parity check).
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>

#include "plink_reader.h"
#include "grm_builder.h"
#include "data_filter.h"
#include "models/block_streamed_grm.h"
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

// --- deterministic synthetic PLINK writer (identical conventions to E1) ------
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

static std::string writeSyntheticData(const std::filesystem::path& dir, int n, int m, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);

    std::vector<uint8_t> dosages((size_t)m * n, 250);  // 250 = missing marker
    std::vector<double> p(m);
    for (int j = 0; j < m; ++j) p[j] = 0.10 + 0.40 * ((rng() % 1000) / 1000.0);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            double uu = u(rng);
            if (uu < 0.03) continue;  // missing
            double p2 = p[j] * p[j];
            double h = 2.0 * p[j] * (1.0 - p[j]);
            if (uu - 0.03 < p2) dosages[(size_t)j * n + i] = 2;
            else if (uu - 0.03 < p2 + h) dosages[(size_t)j * n + i] = 1;
            else dosages[(size_t)j * n + i] = 0;
        }
    }

    std::string prefix = (dir / "e2syn").string();
    {
        FILE* f = std::fopen((prefix + ".bed").c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write bed");
        writeBed(f, dosages, n, m);
        std::fclose(f);
    }
    {
        std::ofstream f(prefix + ".bim");
        for (int j = 0; j < m; ++j)
            f << "1\tsnp_" << j << "\t0\t" << (j + 1) << "\tA\tC\n";
    }
    {
        std::ofstream f(prefix + ".fam");
        for (int i = 0; i < n; ++i)
            f << "0\tind_" << i << "\t0\t0\t0\t-9\n";
    }
    return prefix;
}

// --- dense production-mirror transforms --------------------------------------
static Eigen::MatrixXd tuneG(const Eigen::MatrixXd& G, double a, double b) {
    Eigen::MatrixXd J = Eigen::MatrixXd::Ones(G.rows(), G.cols());
    return b * G + a * J;
}
static Eigen::MatrixXd blendI(const Eigen::MatrixXd& Gt, double alpha) {
    return (1.0 - alpha) * Gt + alpha * Eigen::MatrixXd::Identity(Gt.rows(), Gt.cols());
}
static Eigen::MatrixXd blendT(const Eigen::MatrixXd& Gt, double alpha, const Eigen::MatrixXd& T) {
    return (1.0 - alpha) * Gt + alpha * T;
}

static void computeTune(const Eigen::MatrixXd& G, const Eigen::MatrixXd& A,
                        double& a, double& b) {
    const int n = (int)G.rows();
    double sum_Gd = 0, sum_Ad = 0, sum_Go = 0, sum_Ao = 0;
    for (int i = 0; i < n; ++i) {
        sum_Gd += G(i, i); sum_Ad += A(i, i);
        for (int j = 0; j < n; ++j) if (i != j) { sum_Go += G(i, j); sum_Ao += A(i, j); }
    }
    double n_off = (double)n * (n - 1);
    double mean_Gd = sum_Gd / n, mean_Ad = sum_Ad / n;
    double mean_Go = n_off > 0 ? sum_Go / n_off : 0.0;
    double mean_Ao = n_off > 0 ? sum_Ao / n_off : 0.0;
    b = 1.0;
    if (std::abs(mean_Gd - mean_Go) > 1e-12) b = (mean_Ad - mean_Ao) / (mean_Gd - mean_Go);
    a = mean_Ad - b * mean_Gd;
}

static double maxAbsRel(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    double worst = 0.0;
    for (Eigen::Index j = 0; j < A.cols(); ++j)
        for (Eigen::Index i = 0; i < A.rows(); ++i) {
            const double ref = std::max(1.0, std::fabs(B(i, j)));
            worst = std::max(worst, std::fabs(A(i, j) - B(i, j)) / ref);
        }
    return worst;
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path data_root =
            std::filesystem::temp_directory_path() / "cosmic_e2_stream_parity_data";
        std::filesystem::remove_all(data_root);
        std::filesystem::create_directories(data_root);

        const int N = 48, M = 80;
        std::string prefix =
            argc > 1 ? std::string(argv[1])
                     : writeSyntheticData(data_root, N, M, 20260602u);

        PlinkReader reader;
        reader.load(prefix);
        reader.computeStats(1);
        const int n = (int)reader.getNumSamples();
        const int m = (int)reader.getNumSnps();

        std::vector<int> keep(n);
        for (int i = 0; i < n; ++i) keep[i] = i;

        // Dense raw G (production VanRaden path, algorithm 1).
        GrmBuilder::Options gopts;
        gopts.algorithm = 1;
        gopts.min_maf = 0.05f;
        gopts.min_hwe = 0.0f;
        gopts.thread_num = 1;
        GrmBuilder builder(gopts);
        builder.compute(reader);
        DataFilter::KeepInfo ki;
        ki.keep_indices = keep;
        ki.n_keep = (size_t)n;
        const Eigen::MatrixXd G_raw = DataFilter::subsetSquare(builder.getGrm(), ki);
        const double S_dense = builder.getSum2pq();

        // Pedigree-like A22 target (independent SPD matrix).
        Eigen::MatrixXd A22(n, n);
        {
            Eigen::MatrixXd B(n, 8);
            std::mt19937 rng(7);
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 8; ++k) B(i, k) = ((int)(rng() % 3) - 1) * 0.5;
            A22 = B * B.transpose();
            A22.diagonal().array() += 0.5;
            A22 /= 1.2;
        }
        double tune_a = 0.0, tune_b = 1.0;
        computeTune(G_raw, A22, tune_a, tune_b);
        const double alpha = 0.05;

        const Eigen::MatrixXd G_tuned = tuneG(G_raw, tune_a, tune_b);
        const Eigen::MatrixXd G_fid = blendI(G_tuned, alpha);      // identity blend
        const Eigen::MatrixXd G_fa22 = blendT(G_tuned, alpha, A22);  // A22 blend

        // APY core partition: first K genotyped animals as the core.
        const int K = n / 2;
        std::vector<int> core(K), noncore;
        for (int a = 0; a < K; ++a) core[a] = a;
        for (int i = K; i < n; ++i) noncore.push_back(i);

        std::printf("E2 stream parity (n=%d, m=%d, K=%d, S=%.6f, a=%.6f, b=%.6f, alpha=%.4f)\n",
                    n, m, K, S_dense, tune_a, tune_b, alpha);

        // --- Block-streamed builder (same statistics as the dense build) ---
        BlockStreamedGrmBuilder::Options bopts;
        bopts.min_maf = 0.05f;
        bopts.min_hwe = 0.0f;
        bopts.n_threads = 1;
        BlockStreamedGrmBuilder sb(reader, keep, bopts);
        sb.computeStats();

        // S parity (identical per-SNP float divisor accumulation order).
        CHECK(sb.S2pq() == S_dense, "S = sum 2p(1-p) must match dense bitwise");
        std::printf("  S: dense=%.12f streamed=%.12f\n", S_dense, sb.S2pq());

        // Raw trace / off-diag sums vs dense elementwise sums.
        {
            double dd = 0.0, doff = 0.0;
            for (int i = 0; i < n; ++i) {
                dd += G_raw(i, i);
                for (int j = 0; j < n; ++j) if (i != j) doff += G_raw(i, j);
            }
            std::printf("  raw diag sum: dense=%.12f streamed=%.12f ; raw off sum: dense=%.12f streamed=%.12f\n",
                        dd, sb.rawDiagSum(), doff, sb.rawOffSum());
            CHECK(std::fabs(sb.rawDiagSum() - dd) < 1e-9 * std::max(1.0, std::fabs(dd)),
                  "raw diagonal sum must match dense");
            CHECK(std::fabs(sb.rawOffSum() - doff) < 1e-9 * std::max(1.0, std::fabs(doff)),
                  "raw off-diagonal sum must match dense");
        }

        sb.setTune(tune_a, tune_b);
        sb.computeBlocks(core, noncore);

        // --- RAW block parity -------------------------------------------------
        {
            Eigen::MatrixXd cc(K, K), cn(K, (int)noncore.size());
            Eigen::VectorXd diag_n(noncore.size());
            for (int a = 0; a < K; ++a) {
                for (int b = 0; b < K; ++b) cc(a, b) = G_raw(core[a], core[b]);
                for (int b = 0; b < (int)noncore.size(); ++b) cn(a, b) = G_raw(core[a], noncore[b]);
            }
            for (int b = 0; b < (int)noncore.size(); ++b) diag_n(b) = G_raw(noncore[b], noncore[b]);
            const double t_cc = maxAbsRel(sb.gccRaw(), cc);
            const double t_cn = maxAbsRel(sb.gcnRaw(), cn);
            const double t_dg = (sb.gnnDiagRaw() - diag_n).cwiseAbs().maxCoeff();
            std::printf("  RAW: Gcc rel=%.3e  Gcn rel=%.3e  diag(Gnn) abs=%.3e\n", t_cc, t_cn, t_dg);
            CHECK(t_cc < 1e-9, "raw Gcc must match dense");
            CHECK(t_cn < 1e-9, "raw Gcn must match dense");
            CHECK(t_dg < 1e-9, "raw diag(Gnn) must match dense");
        }

        // --- FINALIZED blocks, identity blend --------------------------------
        sb.finalize(alpha, nullptr);
        {
            Eigen::MatrixXd cc(K, K), cn(K, (int)noncore.size());
            Eigen::VectorXd diag_n(noncore.size());
            for (int a = 0; a < K; ++a) {
                for (int b = 0; b < K; ++b) cc(a, b) = G_fid(core[a], core[b]);
                for (int b = 0; b < (int)noncore.size(); ++b) cn(a, b) = G_fid(core[a], noncore[b]);
            }
            for (int b = 0; b < (int)noncore.size(); ++b) diag_n(b) = G_fid(noncore[b], noncore[b]);
            const double t_cc = maxAbsRel(sb.gccFinal(), cc);
            const double t_cn = maxAbsRel(sb.gcnFinal(), cn);
            const double t_dg = (sb.gnnDiagFinal() - diag_n).cwiseAbs().maxCoeff();
            std::printf("  BLEND_I: Gcc rel=%.3e  Gcn rel=%.3e  diag(Gnn) abs=%.3e\n", t_cc, t_cn, t_dg);
            CHECK(t_cc < 1e-9, "blend-I Gcc must match dense");
            CHECK(t_cn < 1e-9, "blend-I Gcn must match dense");
            CHECK(t_dg < 1e-9, "blend-I diag(Gnn) must match dense");
            // scalar final sums vs dense trace / off sum
            double dd = 0.0, doff = 0.0;
            for (int i = 0; i < n; ++i) {
                dd += G_fid(i, i);
                for (int j = 0; j < n; ++j) if (i != j) doff += G_fid(i, j);
            }
            CHECK(std::fabs(sb.finalDiagSum() - dd) < 1e-9 * std::max(1.0, std::fabs(dd)),
                  "final diagonal sum (identity blend) must match dense");
            CHECK(std::fabs(sb.finalOffSum() - doff) < 1e-9 * std::max(1.0, std::fabs(doff)),
                  "final off-diagonal sum (identity blend) must match dense");
        }

        // --- FINALIZED blocks, A22 blend (also check tuned-unblended copies) --
        sb.finalize(alpha, &A22);
        {
            Eigen::MatrixXd cc(K, K), cn(K, (int)noncore.size());
            Eigen::VectorXd diag_n(noncore.size());
            Eigen::MatrixXd tc(K, K), tcn(K, (int)noncore.size());
            Eigen::VectorXd tdiag_n(noncore.size());
            for (int a = 0; a < K; ++a) {
                for (int b = 0; b < K; ++b) {
                    cc(a, b) = G_fa22(core[a], core[b]);
                    tc(a, b) = G_tuned(core[a], core[b]);
                }
                for (int b = 0; b < (int)noncore.size(); ++b) {
                    cn(a, b) = G_fa22(core[a], noncore[b]);
                    tcn(a, b) = G_tuned(core[a], noncore[b]);
                }
            }
            for (int b = 0; b < (int)noncore.size(); ++b) {
                diag_n(b) = G_fa22(noncore[b], noncore[b]);
                tdiag_n(b) = G_tuned(noncore[b], noncore[b]);
            }
            const double t_cc = maxAbsRel(sb.gccFinal(), cc);
            const double t_cn = maxAbsRel(sb.gcnFinal(), cn);
            const double t_dg = (sb.gnnDiagFinal() - diag_n).cwiseAbs().maxCoeff();
            const double t_tc = maxAbsRel(sb.gccTuned(), tc);
            const double t_tcn = maxAbsRel(sb.gcnTuned(), tcn);
            const double t_tdg = (sb.gnnDiagTuned() - tdiag_n).cwiseAbs().maxCoeff();
            std::printf("  BLEND_A22: final Gcc rel=%.3e  Gcn rel=%.3e  diag rel=%.3e ; "
                        "tuned Gcc rel=%.3e  Gcn rel=%.3e  diag rel=%.3e\n",
                        t_cc, t_cn, t_dg, t_tc, t_tcn, t_tdg);
            CHECK(t_cc < 1e-9, "A22-blend Gcc must match dense");
            CHECK(t_cn < 1e-9, "A22-blend Gcn must match dense");
            CHECK(t_dg < 1e-9, "A22-blend diag(Gnn) must match dense");
            CHECK(t_tc < 1e-9, "tuned Gcc must match dense");
            CHECK(t_tcn < 1e-9, "tuned Gcn must match dense");
            CHECK(t_tdg < 1e-9, "tuned diag(Gnn) must match dense");
            // scalar final sums vs dense
            double dd = 0.0, doff = 0.0;
            for (int i = 0; i < n; ++i) {
                dd += G_fa22(i, i);
                for (int j = 0; j < n; ++j) if (i != j) doff += G_fa22(i, j);
            }
            std::printf("  A22 scalar sums: dense diag=%.12f streamed=%.12f ; dense off=%.12f streamed=%.12f\n",
                        dd, sb.finalDiagSum(), doff, sb.finalOffSum());
            CHECK(std::fabs(sb.finalDiagSum() - dd) < 1e-9 * std::max(1.0, std::fabs(dd)),
                  "final diagonal sum (A22 blend) must match dense");
            CHECK(std::fabs(sb.finalOffSum() - doff) < 1e-9 * std::max(1.0, std::fabs(doff)),
                  "final off-diagonal sum (A22 blend) must match dense");
        }

        // --- APY_OPERATOR_PARITY: ApyGinv::build(G) vs ApyGinv::build(blocks) --
        {
            // finalize with A22 blend to also exercise the tuned block copies.
            sb.finalize(alpha, &A22);

            ApyGinv g_dense;
            g_dense.build(G_fa22, core, "dense-G");  // A22-blend final (matches the finalized blocks)
            std::vector<int> nc = noncore;
            ApyGinv g_blk;
            g_blk.build(sb.gccFinal(), sb.gcnFinal(), sb.gnnDiagFinal(), n, core, "blocks");

            const double dq = maxAbsRel(g_blk.q(), g_dense.q());
            const double dmnn = (g_blk.mnn() - g_dense.mnn()).cwiseAbs().maxCoeff();
            const double dgcc_inv = maxAbsRel(g_blk.gccInvFull(), g_dense.gccInvFull());
            const double dcore = maxAbsRel(g_blk.coreBlock(), g_dense.coreBlock());
            std::printf("  APY(g) vs APY(blocks): q rel=%.3e  mnn abs=%.3e  gcc_inv rel=%.3e  "
                        "core_block rel=%.3e\n", dq, dmnn, dgcc_inv, dcore);
            CHECK(dq < 1e-9, "APY build(blocks).q must match build(G).q");
            CHECK(dmnn < 1e-9, "APY build(blocks).mnn must match build(G).mnn");
            CHECK(dgcc_inv < 1e-9, "APY build(blocks).gccInv must match build(G).gccInv");
            CHECK(dcore < 1e-9, "APY build(blocks).coreBlock must match build(G).coreBlock");

            // applyInverse / applyCoreInverse equivalence.
            std::mt19937 rng(1234);
            Eigen::VectorXd v(n), r1(n), r2(n);
            for (int i = 0; i < n; ++i) v(i) = ((rng() % 1000) / 1000.0) - 0.5;
            r1 = g_dense.applyInverse(v);
            r2 = g_blk.applyInverse(v);
            const double d_ai = (r1 - r2).cwiseAbs().maxCoeff();
            Eigen::VectorXd vc(K), c1(K), c2(K);
            for (int i = 0; i < K; ++i) vc(i) = ((rng() % 1000) / 1000.0) - 0.5;
            c1 = g_dense.applyCoreInverse(vc);
            c2 = g_blk.applyCoreInverse(vc);
            const double d_ci = (c1 - c2).cwiseAbs().maxCoeff();
            std::printf("  APY applyInverse diff=%.3e  applyCoreInverse diff=%.3e\n", d_ai, d_ci);
            CHECK(d_ai < 1e-9, "APY applyInverse must match between build paths");
            CHECK(d_ci < 1e-9, "APY applyCoreInverse must match between build paths");

            // Structural info.
            CHECK((int)g_blk.nCore() == K, "block build n_core must match");
            CHECK((int)g_blk.nNoncore() == (int)noncore.size(), "block build n_noncore must match");
        }

        // --- CORE_N_EQUALS_NG_DEGENERACY ---------------------------------------
        {
            std::vector<int> all(n);
            for (int i = 0; i < n; ++i) all[i] = i;
            std::vector<int> empty;
            sb.computeBlocks(all, empty);
            sb.finalize(alpha, nullptr);  // identity blend -> G_fid

            ApyGinv g_blk_all;
            g_blk_all.build(sb.gccFinal(), sb.gcnFinal(), sb.gnnDiagFinal(), n, all, "blocks-all");
            CHECK((int)g_blk_all.nNoncore() == 0, "core=n_g must give empty noncore");
            CHECK(g_blk_all.mnn().size() == 0, "core=n_g must give empty mnn");
            CHECK(g_blk_all.q().size() == 0, "core=n_g must give empty Q");

            Eigen::MatrixXd G_inv = G_fid.inverse();
            Eigen::VectorXd v(n), r(n);
            std::mt19937 rng(55);
            for (int i = 0; i < n; ++i) v(i) = ((rng() % 1000) / 1000.0) - 0.5;
            r = g_blk_all.applyInverse(v);
            const double d = (r - G_inv * v).cwiseAbs().maxCoeff();
            std::printf("  CORE=n: applyInverse vs G^-1 diff=%.3e\n", d);
            CHECK(d < 1e-8, "core=n_g block build must reproduce the full inverse");
        }

    } catch (const std::exception& e) {
        ++g_failures;
        std::printf("FATAL: %s\n", e.what());
    }

    if (g_failures == 0) {
        std::printf("E2 STREAM PARITY: PASS\n");
    } else {
        std::printf("E2 STREAM PARITY: FAIL (%d failures)\n", g_failures);
    }
    return g_failures == 0 ? 0 : 1;
}
