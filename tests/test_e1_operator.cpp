// Phase E1 unit test: matrix-free genomic operator parity with the production
// dense GRM pipeline (tuneG + blend), operator linearity, trace parity, full
// matrix reconstruction, and Lanczos-vs-dense K98/K99 agreement.
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
#include "models/lanczos_eigensolver.h"

using namespace cosmic;

static int g_failures = 0;

#define CHECK(cond, msg)                                                          \
    do {                                                                          \
        if (!(cond)) {                                                            \
            ++g_failures;                                                         \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, msg);             \
        }                                                                         \
    } while (0)

// --- deterministic synthetic PLINK writer ------------------------------------
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
            double h  = 2.0 * p[j] * (1.0 - p[j]);
            if (uu - 0.03 < p2) dosages[(size_t)j * n + i] = 2;
            else if (uu - 0.03 < p2 + h) dosages[(size_t)j * n + i] = 1;
            else dosages[(size_t)j * n + i] = 0;
        }
    }

    std::string prefix = (dir / "e1syn").string();
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

// Compute production-style tune factors (preGSf90: align G to target A22).
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

struct StageResult {
    const char* label;
    double max_abs;      // max |Op(e_i) - G_dense e_i|
    double frob_rel;     // ||G_op_recon - G_dense||_F / ||G_dense||_F
    double trace_dense;
    double trace_basis;  // sum_i e_i' Op(e_i)
    double trace_analytic;
};

// Reconstruct the full matrix from Op applied to basis vectors and compare.
static StageResult parityStage(const MatrixFreeGrmOperator& op,
                               const Eigen::MatrixXd& G_dense,
                               const char* label,
                               double analytic_trace) {
    StageResult r;
    r.label = label;
    const int n = (int)G_dense.rows();

    // Column-by-column reconstruction.
    Eigen::MatrixXd G_op_recon(n, n);
    for (int i = 0; i < n; ++i) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(n);
        e(i) = 1.0;
        Eigen::VectorXd col(n);
        op.apply(e, col);
        G_op_recon.col(i) = col;
    }

    const Eigen::MatrixXd diff = G_op_recon - G_dense;
    r.max_abs = diff.cwiseAbs().maxCoeff();
    const double normF = G_dense.norm();
    r.frob_rel = normF > 0 ? diff.norm() / normF : 0.0;

    r.trace_dense = G_dense.trace();
    r.trace_basis = 0.0;
    for (int i = 0; i < n; ++i) r.trace_basis += G_op_recon(i, i);
    r.trace_analytic = analytic_trace;

    std::printf("  [%s] max|diff|=%.3e frob_rel=%.3e trace(dense)=%.12f "
                "trace(basis)=%.12f trace(analytic)=%.12f\n",
                label, r.max_abs, r.frob_rel, r.trace_dense, r.trace_basis,
                r.trace_analytic);
    return r;
}

// Analytic raw-G trace from genotype columns: sum_j sum_{i obs} (x_ij-mean_j)^2 / S,
// with missing imputed to the mean (Z=0 contribution), mirroring PlinkReader/GrmBuilder.
static double analyticRawTrace(PlinkReader& reader, double S) {
    const int n = (int)reader.getNumSamples();
    const int m = (int)reader.getNumSnps();
    double sum = 0.0;
    for (int j = 0; j < m; ++j) {
        Eigen::VectorXf x;
        reader.getSnpGenotypes(j, x);
        double mean = 0.0; int cnt = 0;
        for (int i = 0; i < n; ++i) if (!std::isnan(x(i))) { mean += x(i); ++cnt; }
        mean = cnt > 0 ? mean / cnt : 0.0;
        double ss = 0.0;
        for (int i = 0; i < n; ++i)
            if (!std::isnan(x(i))) { double d = x(i) - mean; ss += d * d; }
        sum += ss;
    }
    return S > 0 ? sum / S : 0.0;
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path data_root =
            std::filesystem::temp_directory_path() / "cosmic_e1_operator_data";
        std::filesystem::remove_all(data_root);
        std::filesystem::create_directories(data_root);

        const int N = 48, M = 80;
        std::string prefix =
            argc > 1 ? std::string(argv[1])
                     : writeSyntheticData(data_root, N, M, 20260601u);

        PlinkReader reader;
        reader.load(prefix);
        reader.computeStats(1);
        const int n = (int)reader.getNumSamples();
        const int m = (int)reader.getNumSnps();

        std::vector<int> keep(n);
        for (int i = 0; i < n; ++i) keep[i] = i;

        // Dense raw G (production VanRaden path).
        GrmBuilder::Options gopts;
        gopts.algorithm = 1;             // VanRaden
        gopts.min_maf = 0.05f;
        gopts.min_hwe = 0.0f;
        gopts.thread_num = 1;
        GrmBuilder builder(gopts);
        builder.compute(reader);
        DataFilter::KeepInfo ki;
        ki.keep_indices = keep;
        ki.n_keep = (size_t)n;
        const Eigen::MatrixXd G_raw = DataFilter::subsetSquare(builder.getGrm(), ki);
        const double S = builder.getSum2pq();
        const double analytic_raw = analyticRawTrace(reader, S);

        const double alpha = 0.05;

        // Tune factors (align raw G to a stand-in A22): use the raw G itself's
        // block structure via a diagonal-dominant target - simplest is to reuse
        // the raw G shifted, but use an independent pedigree-like SPD target.
        Eigen::MatrixXd A22(n, n);
        {
            // Simple SPD "pedigree-like" block: A = 0.5 I + low-rank kinship.
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

        const Eigen::MatrixXd G_tuned = tuneG(G_raw, tune_a, tune_b);
        const Eigen::MatrixXd G_blendI = blendI(G_tuned, alpha);
        const Eigen::MatrixXd G_blendA = blendT(G_tuned, alpha, A22);

        // Operator instances.
        MatrixFreeGrmOperator op_raw(reader, keep, S, 0.0, 1.0, 0.0, nullptr, 1);
        MatrixFreeGrmOperator op_tuned(reader, keep, S, tune_a, tune_b, 0.0, nullptr, 1);
        MatrixFreeGrmOperator op_blendI(reader, keep, S, tune_a, tune_b, alpha, nullptr, 1);
        MatrixFreeGrmOperator op_blendA(reader, keep, S, tune_a, tune_b, alpha, &A22, 1);

        std::printf("E1 operator parity (n=%d, m=%d, S=%.6f, a=%.6f, b=%.6f, alpha=%.4f)\n",
                    n, m, S, tune_a, tune_b, alpha);

        // --- Operator linearity (raw) ---
        {
            Eigen::VectorXd x(n), y(n), r1(n), r2(n), r3(n);
            std::mt19937 rng(99);
            for (int i = 0; i < n; ++i) { x(i) = ((rng() % 1000) / 1000.0) - 0.3; y(i) = ((rng() % 1000) / 1000.0) + 0.2; }
            Eigen::VectorXd xy = x + y;
            op_raw.apply(xy, r1);
            op_raw.apply(x, r2); op_raw.apply(y, r3);
            double add_diff = (r1 - (r2 + r3)).cwiseAbs().maxCoeff();
            Eigen::VectorXd cx = 1.5 * x, cr(n);
            op_raw.apply(cx, r1);
            op_raw.apply(x, cr);
            double scale_diff = (r1 - 1.5 * cr).cwiseAbs().maxCoeff();
            std::printf("  linearity: Op(x+y)-Op(x)-Op(y) max abs = %.3e ; Op(cx)-c*Op(x) = %.3e\n",
                        add_diff, scale_diff);
            // The operator's per-SNP genotype math is float32 inside PlinkReader
            // (val_map computed in the reader); residual nonlinearity is bounded
            // by that rounding. Gate at the float32 path precision (~1e-7).
            CHECK(add_diff < 1e-6, "operator must be additive (Op(x+y)==Op(x)+Op(y))");
            CHECK(scale_diff < 1e-6, "operator must be homogeneous (Op(cx)==c*Op(x))");
            // Symmetric conjugate: y' Op x = x' Op y (in exact arithmetic Op is symmetric)
            double yox = y.dot(r2);  // y' (Op x)
            Eigen::VectorXd oy(n); op_raw.apply(y, oy);
            double xoy = x.dot(oy);  // x' (Op y)
            std::printf("  symmetry: y'Op x = %.10f vs x'Op y = %.10f (diff %.3e)\n",
                        yox, xoy, std::abs(yox - xoy));
            CHECK(std::abs(yox - xoy) < 1e-6, "operator must be symmetric in exact arithmetic");
        }

        // --- Parity stages ---
        std::printf("  Stage results:\n");
        StageResult r_raw = parityStage(op_raw, G_raw, "raw", analytic_raw);
        // Tuned/blended analytic traces derived from analytic raw trace.
        const double n_d = (double)n;
        StageResult r_tuned = parityStage(op_tuned, G_tuned, "tuned",
                                          tune_b * analytic_raw + tune_a * n_d);
        StageResult r_blendI = parityStage(op_blendI, G_blendI, "blend_I",
                                           (1.0 - alpha) * (tune_b * analytic_raw + tune_a * n_d) + alpha * n_d);
        StageResult r_blendA = parityStage(op_blendA, G_blendA, "blend_A22",
                                           (1.0 - alpha) * (tune_b * analytic_raw + tune_a * n_d) + alpha * A22.trace());

        auto gateStage = [&](const StageResult& rsp, double tol_abs) {
            CHECK(rsp.max_abs < tol_abs, "operator column max-abs parity exceeded");
            CHECK(rsp.trace_basis == rsp.trace_dense ||
                  std::abs(rsp.trace_basis - rsp.trace_dense) < 1e-6,
                  "operator traced from basis must match dense trace");
        };
        gateStage(r_raw, 1e-6);
        gateStage(r_tuned, 1e-6);
        gateStage(r_blendI, 1e-6);
        gateStage(r_blendA, 1e-6);

        // --- K98 / K99: dense final-blend vs Lanczos ---
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G_blendI);
            Eigen::VectorXd ev = es.eigenvalues();
            std::sort(ev.data(), ev.data() + ev.size(), std::greater<double>());
            const double tr = G_blendI.trace();
            auto findK = [&](double t) { double c = 0; for (int k = 0; k < n; ++k) { c += std::max(0.0, ev(k)); if (c >= t * tr) return k + 1; } return n; };
            const int k98_d = findK(0.98), k99_d = findK(0.99);

            LanczosEigenSolver lz;
            LanczosEigenSolver::Options lo;
            lo.max_iterations = n;   // full basis
            lo.tolerance = 1e-10;
            lo.seed = 20260601u;
            lo.n_threads = 1;
            LanczosResult lr = lz.compute(op_blendI, 0.99, lo);

            std::printf("  K98 dense=%d lanczos=%d ; K99 dense=%d lanczos=%d ; "
                        "trace dense=%.12f lanczos(ritz)=%.12f\n",
                        k98_d, lr.K98, k99_d, lr.K99, tr, lr.trace_g);
            CHECK(lr.status == "OK", "lanczos status must be OK");
            CHECK(lr.K98 == k98_d, "K98 must match dense exactly");
            CHECK(lr.K99 == k99_d, "K99 must match dense exactly");
            CHECK(std::abs(lr.trace_g - tr) < 1e-6, "lanczos trace must match dense trace");
        }

        // --- Reported full reconstruction on the final production matrix ---
        {
            const Eigen::MatrixXd& G = G_blendI;
            Eigen::MatrixXd Rec(n, n);
            for (int i = 0; i < n; ++i) {
                Eigen::VectorXd e = Eigen::VectorXd::Zero(n); e(i) = 1.0;
                Eigen::VectorXd c(n); op_blendI.apply(e, c); Rec.col(i) = c;
            }
            const Eigen::MatrixXd D = Rec - G;
            double maxdiff = D.cwiseAbs().maxCoeff();
            double rmse = D.norm() / std::sqrt((double)n * n);
            double frob = D.norm() / G.norm();
            // largest differing cell
            Eigen::Index gi = 0, gj = 0;
            D.cwiseAbs().maxCoeff(&gi, &gj);
            std::printf("  recon final blended: max|D|=%.3e RMSE=%.3e frob_rel=%.3e "
                        "largest cell (%d,%d)=%.6f vs %.6f\n",
                        maxdiff, rmse, frob, (int)gi, (int)gj, Rec((int)gi, (int)gj), G((int)gi, (int)gj));
            CHECK(maxdiff < 1e-6, "full reconstructed blended matrix must match dense");
        }

    } catch (const std::exception& e) {
        ++g_failures;
        std::printf("FATAL: %s\n", e.what());
    }

    if (g_failures == 0) {
        std::printf("E1 OPERATOR PARITY: PASS\n");
    } else {
        std::printf("E1 OPERATOR PARITY: FAIL (%d failures)\n", g_failures);
    }
    return g_failures == 0 ? 0 : 1;
}