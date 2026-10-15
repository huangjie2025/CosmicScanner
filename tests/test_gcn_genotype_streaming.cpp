#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "models/apy_ginv.h"
#include "models/apy_inverse_operator.h"
#include "models/block_streamed_grm.h"
#include "models/gcc_solver_pcg.h"
#include "plink_reader.h"

using namespace cosmic;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { ++failures; std::printf("FAIL: %s\n", msg); } } while (0)

static void writeBed(FILE* file, const std::vector<unsigned char>& dosages, int n, int m) {
    const unsigned char magic[3] = {0x6c, 0x1b, 0x01};
    std::fwrite(magic, 1, 3, file);
    for (int j = 0; j < m; ++j) {
        std::vector<unsigned char> bytes((n + 3) / 4, 0);
        for (int i = 0; i < n; ++i) {
            const unsigned char dosage = dosages[static_cast<size_t>(j) * n + i];
            unsigned char code = 1;
            if (dosage < 3) code = dosage == 0 ? 3 : (dosage == 1 ? 2 : 0);
            bytes[i / 4] |= code << (2 * (i % 4));
        }
        std::fwrite(bytes.data(), 1, bytes.size(), file);
    }
}

static std::string writeFixture() {
    const int n = 36;
    const int m = 90;
    std::mt19937 rng(20260920);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<unsigned char> dosages(static_cast<size_t>(n) * m, 3);
    for (int j = 0; j < m; ++j) {
        const double p = 0.12 + 0.36 * ((rng() % 1000) / 1000.0);
        for (int i = 0; i < n; ++i) {
            const double u = uniform(rng);
            if (u < 0.025) continue;
            const double q = (u - 0.025) / 0.975;
            if (q < p * p) dosages[static_cast<size_t>(j) * n + i] = 2;
            else if (q < p * p + 2.0 * p * (1.0 - p)) dosages[static_cast<size_t>(j) * n + i] = 1;
            else dosages[static_cast<size_t>(j) * n + i] = 0;
        }
    }
    const auto dir = std::filesystem::temp_directory_path() / "cosmic_gcn_genotype_streaming";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const std::string prefix = (dir / "fixture").string();
    FILE* bed = std::fopen((prefix + ".bed").c_str(), "wb");
    if (!bed) throw std::runtime_error("cannot create Gcn streaming BED fixture");
    writeBed(bed, dosages, n, m);
    std::fclose(bed);
    {
        std::ofstream out(prefix + ".bim");
        for (int j = 0; j < m; ++j) out << "1 rs" << j << " 0 " << j + 1 << " A C\n";
    }
    {
        std::ofstream out(prefix + ".fam");
        for (int i = 0; i < n; ++i) out << "0 id" << i << " 0 0 0 -9\n";
    }
    return prefix;
}

static Eigen::MatrixXd assemble(BlockStreamedGrmBuilder& builder, int block_size,
                                GcnBlockView view) {
    Eigen::MatrixXd result(builder.nCore(), builder.nNoncore());
    auto iterator = builder.makeGcnIterator(block_size, view);
    GcnBlock block;
    int expected = 0;
    while (iterator->next(block)) {
        CHECK(block.start_col == expected, "contiguous Gcn block sequence");
        result.middleCols(block.start_col, block.end_col - block.start_col) = block.block;
        expected = block.end_col;
    }
    CHECK(expected == builder.nNoncore(), "complete Gcn block sequence");
    return result;
}

static void reportParity(const char* label, const Eigen::MatrixXd& streamed,
                         const Eigen::MatrixXd& reference) {
    const Eigen::MatrixXd delta = streamed - reference;
    const double max_abs = delta.size() ? delta.cwiseAbs().maxCoeff() : 0.0;
    const double rmse = delta.size() ? std::sqrt(delta.squaredNorm() / delta.size()) : 0.0;
    const double trace_contribution =
        std::fabs((streamed * streamed.transpose()).trace() -
                  (reference * reference.transpose()).trace());
    std::printf("%s max_abs=%.3e rmse=%.3e trace_contribution_diff=%.3e\n",
                label, max_abs, rmse, trace_contribution);
    CHECK(max_abs < 1e-12, "Gcn block max-absolute parity");
    CHECK(rmse < 1e-13, "Gcn block RMSE parity");
    CHECK(trace_contribution < 1e-11, "Gcn trace-contribution parity");
}

int main(int argc, char** argv) {
    try {
        PlinkReader reader;
        reader.load(argc > 1 ? std::string(argv[1]) : writeFixture());
        reader.computeStats(1);
        const int n = static_cast<int>(reader.getNumSamples());
        std::vector<int> keep(n);
        for (int i = 0; i < n; ++i) keep[i] = i;
        const std::vector<int> core{0, 2, 4, 7, 9, 12, 15, 18, 22, 26, 30, 34};
        std::vector<int> noncore;
        for (int i = 0; i < n; ++i) {
            if (!std::binary_search(core.begin(), core.end(), i)) noncore.push_back(i);
        }

        BlockStreamedGrmBuilder::Options reference_options;
        reference_options.min_maf = 0.01f;
        reference_options.n_threads = 1;
        reference_options.materialize_gcn = true;
        BlockStreamedGrmBuilder reference(reader, keep, reference_options);
        reference.computeStats();
        reference.setTune(0.03, 0.91);
        reference.computeBlocks(core, noncore);

        BlockStreamedGrmBuilder::Options streaming_options = reference_options;
        streaming_options.materialize_gcn = false;
        BlockStreamedGrmBuilder streaming(reader, keep, streaming_options);
        streaming.computeStats();
        streaming.setTune(0.03, 0.91);
        streaming.computeBlocks(core, noncore);
        CHECK(!streaming.fullGcnStored(), "genotype-backed builder must not store full Gcn");

        Eigen::MatrixXd a22 = Eigen::MatrixXd::Identity(n, n);
        a22.array() += 0.01;
        reference.finalize(0.05, &a22);
        streaming.finalize(0.05, &a22);

        for (int block_size : {1, 3, 7, 64}) {
            const std::string raw_label = "raw block=" + std::to_string(block_size);
            const std::string tuned_label = "tuned block=" + std::to_string(block_size);
            const std::string final_label = "final block=" + std::to_string(block_size);
            reportParity(raw_label.c_str(), assemble(streaming, block_size, GcnBlockView::Raw),
                         reference.gcnRaw());
            reportParity(tuned_label.c_str(), assemble(streaming, block_size, GcnBlockView::Tuned),
                         reference.gcnTuned());
            reportParity(final_label.c_str(), assemble(streaming, block_size, GcnBlockView::Final),
                         reference.gcnFinal());
        }

        ApyGinv dense_apy;
        dense_apy.build(reference.gccFinal(), reference.gcnFinal(), reference.gnnDiagFinal(),
                        n, core, "materialized-Gcn-reference");
        GcnBlockIteratorFactory final_factory = [&streaming]() {
            return streaming.makeGcnIterator(5, GcnBlockView::Final);
        };
        ApyGinv streamed_apy;
        streamed_apy.build(streaming.gccFinal(), streaming.gnnDiagFinal(), n, core,
                           "genotype-backed-Gcn", final_factory);
        const double q_error = (streamed_apy.q() - dense_apy.q()).cwiseAbs().maxCoeff();
        const double mnn_error = (streamed_apy.mnn() - dense_apy.mnn()).cwiseAbs().maxCoeff();
        std::printf("Q max_abs=%.3e MNN max_abs=%.3e\n", q_error, mnn_error);
        CHECK(q_error < 1e-12, "Q parity");
        CHECK(mnn_error < 1e-12, "MNN parity");

        MatrixFreeApyInverseOperator apy_operator(
            streamed_apy.gccSolverPtr(), final_factory, streaming.gnnDiagFinal(),
            n, core, "genotype-backed-operator");
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Random(n, 3);
        Eigen::MatrixXd dense_result(n, rhs.cols());
        for (int j = 0; j < rhs.cols(); ++j) {
            dense_result.col(j) = dense_apy.applyInverse(rhs.col(j));
        }
        const double apy_error = (apy_operator.apply(rhs) - dense_result).cwiseAbs().maxCoeff();
        std::printf("APY operator max_abs=%.3e\n", apy_error);
        CHECK(apy_error < 1e-11, "APY operator parity");

        // PCG backend over the same genotype-backed Gcn source. The synthetic
        // fixture has no marker-filter mismatch with GccOperator.
        if (argc == 1) {
            Eigen::MatrixXd a22_cc(core.size(), core.size());
            for (int i = 0; i < static_cast<int>(core.size()); ++i) {
                for (int j = 0; j < static_cast<int>(core.size()); ++j) {
                    a22_cc(i, j) = a22(core[i], core[j]);
                }
            }
            GccOperator gcc_operator(reader, core, streaming.S2pq(), 0.03, 0.91,
                                     0.05, &a22_cc, 1);
            auto pcg = makePCGSolver();
            GccSolverPCGOptions pcg_options;
            pcg_options.tolerance = 1e-11;
            pcg_options.max_iter = 200;
            pcg_options.preconditioner = "jacobi";
            pcg->compute(gcc_operator, pcg_options);
            MatrixFreeApyInverseOperator pcg_operator(
                pcg, final_factory, streaming.gnnDiagFinal(), n, core,
                "genotype-backed-PCG-operator");
            const double pcg_error = (pcg_operator.apply(rhs) - dense_result).cwiseAbs().maxCoeff();
            std::printf("PCG APY operator max_abs=%.3e\n", pcg_error);
            CHECK(pcg_error < 1e-6, "PCG genotype-backed APY parity");
        }

        // core=n_g: no Gcn allocation and APY reduces to the exact core solve.
        std::vector<int> all(n);
        for (int i = 0; i < n; ++i) all[i] = i;
        BlockStreamedGrmBuilder all_core(reader, keep, streaming_options);
        all_core.computeStats();
        all_core.computeBlocks(all, {});
        all_core.finalize(0.05, nullptr);
        CHECK(!all_core.fullGcnStored(), "core=n_g must not allocate Gcn");
        GcnBlockIteratorFactory empty_factory = [&all_core]() {
            return all_core.makeGcnIterator(3, GcnBlockView::Final);
        };
        ApyGinv all_core_apy;
        all_core_apy.build(all_core.gccFinal(), all_core.gnnDiagFinal(), n, all,
                           "all-core-genotype-stream", empty_factory);
        const Eigen::VectorXd x = Eigen::VectorXd::Random(n);
        CHECK((all_core_apy.applyInverse(x) - all_core.gccFinal().ldlt().solve(x))
                  .cwiseAbs().maxCoeff() < 1e-10,
              "core=n_g degeneracy");
    } catch (const std::exception& error) {
        ++failures;
        std::printf("FATAL: %s\n", error.what());
    }

    std::printf("GCN GENOTYPE STREAMING E3-C3: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
