#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "models/apy_ginv.h"
#include "models/apy_inverse_operator.h"
#include "models/apy_mme_operator.h"
#include "models/gcc_solver_pcg.h"
#include "plink_reader.h"

using namespace cosmic;

static int failures = 0;

#define CHECK(cond, msg) do { if (!(cond)) { ++failures; std::printf("FAIL: %s\n", msg); } } while (0)

class MatrixBlockIterator final : public IGcnBlockIterator {
public:
    MatrixBlockIterator(const Eigen::MatrixXd& gcn, const std::vector<int>& noncore, int block_size)
        : gcn_(gcn), noncore_(noncore), block_size_(block_size) {}

    bool next(GcnBlock& out) override {
        if (!hasNext()) return false;
        const int end = std::min(start_ + block_size_, static_cast<int>(gcn_.cols()));
        out.start_col = start_;
        out.end_col = end;
        out.block = gcn_.middleCols(start_, end - start_);
        out.noncore_indices.assign(noncore_.begin() + start_, noncore_.begin() + end);
        start_ = end;
        return true;
    }
    bool hasNext() const override { return start_ < gcn_.cols(); }
    int totalBlocks() const override { return (gcn_.cols() + block_size_ - 1) / block_size_; }
    int blockSize() const override { return block_size_; }

private:
    const Eigen::MatrixXd& gcn_;
    const std::vector<int>& noncore_;
    int block_size_;
    int start_ = 0;
};

static std::vector<int> complement(int n, const std::vector<int>& core) {
    std::vector<int> result;
    for (int i = 0; i < n; ++i) {
        if (!std::binary_search(core.begin(), core.end(), i)) result.push_back(i);
    }
    return result;
}

static GcnBlockIteratorFactory matrixFactory(const Eigen::MatrixXd& gcn,
                                             const std::vector<int>& noncore,
                                             int block_size) {
    return [&gcn, &noncore, block_size]() {
        return std::make_unique<MatrixBlockIterator>(gcn, noncore, block_size);
    };
}

static Eigen::MatrixXd randomSpd(int n, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    Eigen::MatrixXd a(n, n + 12);
    for (Eigen::Index j = 0; j < a.cols(); ++j)
        for (Eigen::Index i = 0; i < a.rows(); ++i) a(i, j) = normal(rng);
    Eigen::MatrixXd g = a * a.transpose() / static_cast<double>(a.cols());
    g.diagonal().array() += 0.25;
    return g;
}

static void extractBlocks(const Eigen::MatrixXd& g, const std::vector<int>& core,
                          const std::vector<int>& noncore, Eigen::MatrixXd& gcc,
                          Eigen::MatrixXd& gcn, Eigen::VectorXd& gnn_diag) {
    gcc.resize(core.size(), core.size());
    gcn.resize(core.size(), noncore.size());
    gnn_diag.resize(noncore.size());
    for (int i = 0; i < static_cast<int>(core.size()); ++i) {
        for (int j = 0; j < static_cast<int>(core.size()); ++j) gcc(i, j) = g(core[i], core[j]);
        for (int j = 0; j < static_cast<int>(noncore.size()); ++j) gcn(i, j) = g(core[i], noncore[j]);
    }
    for (int j = 0; j < static_cast<int>(noncore.size()); ++j) gnn_diag(j) = g(noncore[j], noncore[j]);
}

static void testDenseAndMatrixFreeParity() {
    const int n = 31;
    const std::vector<int> core{0, 2, 5, 9, 14, 20, 27};
    const std::vector<int> noncore = complement(n, core);
    const Eigen::MatrixXd g = randomSpd(n, 7123);
    Eigen::MatrixXd gcc, gcn;
    Eigen::VectorXd gnn_diag;
    extractBlocks(g, core, noncore, gcc, gcn, gnn_diag);

    ApyGinv reference;
    reference.build(gcc, gcn, gnn_diag, n, core, "operator-reference");
    DenseApyInverseOperator dense(reference);
    MatrixFreeApyInverseOperator matrix_free(
        reference.gccSolverPtr(), matrixFactory(gcn, noncore, 5), gnn_diag, n, core, "operator-dense");

    std::mt19937 rng(99);
    std::normal_distribution<double> normal(0.0, 1.0);
    for (int rep = 0; rep < 3; ++rep) {
        Eigen::VectorXd x(n);
        for (int i = 0; i < n; ++i) x(i) = normal(rng);
        const Eigen::VectorXd expected = reference.applyInverse(x);
        const Eigen::VectorXd wrapped = dense.apply(x);
        const Eigen::VectorXd actual = matrix_free.apply(x);
        const double abs_error = (actual - expected).cwiseAbs().maxCoeff();
        const double rel_error = (actual - expected).norm() / std::max(1.0, expected.norm());
        std::printf("random vector %d: max_abs=%.3e relative=%.3e\n", rep + 1, abs_error, rel_error);
        CHECK((wrapped - expected).cwiseAbs().maxCoeff() < 1e-13, "dense wrapper vector parity");
        CHECK(abs_error < 1e-11, "matrix-free random-vector absolute parity");
        CHECK(rel_error < 1e-12, "matrix-free random-vector relative parity");
    }

    Eigen::MatrixXd x(n, 4);
    for (Eigen::Index j = 0; j < x.cols(); ++j)
        for (Eigen::Index i = 0; i < x.rows(); ++i) x(i, j) = normal(rng);
    const Eigen::MatrixXd yd = dense.apply(x);
    const Eigen::MatrixXd yo = matrix_free.apply(x);
    CHECK((yo - yd).cwiseAbs().maxCoeff() < 1e-11, "multiple-RHS parity");
    CHECK((matrix_free.mnn() - reference.mnn()).cwiseAbs().maxCoeff() < 1e-12, "streamed Mnn parity");
}

static void testCoreIsAll() {
    const int n = 17;
    std::vector<int> core(n);
    for (int i = 0; i < n; ++i) core[i] = i;
    const std::vector<int> noncore;
    const Eigen::MatrixXd g = randomSpd(n, 41);
    const Eigen::MatrixXd gcn(n, 0);
    const Eigen::VectorXd gnn_diag(0);
    ApyGinv reference;
    reference.build(g, core, "operator-core-all");
    MatrixFreeApyInverseOperator matrix_free(
        reference.gccSolverPtr(), matrixFactory(gcn, noncore, 3), gnn_diag, n, core, "operator-core-all");
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, 3);
    CHECK((matrix_free.apply(x) - g.ldlt().solve(x)).cwiseAbs().maxCoeff() < 1e-11,
          "core=n_g degeneracy");
}

static std::string writePcgFixture() {
    const int n = 12;
    const int m = 80;
    const std::filesystem::path dir =
        std::filesystem::temp_directory_path() / "cosmic_apy_operator_pcg";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    const std::string prefix = (dir / "fixture").string();
    {
        FILE* f = std::fopen((prefix + ".bed").c_str(), "wb");
        const unsigned char magic[3] = {0x6c, 0x1b, 0x01};
        std::fwrite(magic, 1, 3, f);
        for (int j = 0; j < m; ++j) {
            std::vector<unsigned char> bytes((n + 3) / 4, 0);
            for (int i = 0; i < n; ++i) {
                const int dosage = (i * 7 + j * 5 + (i * j) % 3) % 3;
                const unsigned char code = dosage == 0 ? 3 : (dosage == 1 ? 2 : 0);
                bytes[i / 4] |= code << (2 * (i % 4));
            }
            std::fwrite(bytes.data(), 1, bytes.size(), f);
        }
        std::fclose(f);
    }
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

static void testPcgBackend() {
    PlinkReader reader;
    reader.load(writePcgFixture());
    reader.computeStats(1);
    const std::vector<int> core{0, 1, 2, 3, 4, 5};
    const int n = 12;
    const std::vector<int> noncore = complement(n, core);
    GccOperator gcc_op(reader, core, 1.0, 0.0, 1.0, 0.25, nullptr, 1);
    Eigen::MatrixXd gcc(core.size(), core.size());
    for (int j = 0; j < static_cast<int>(core.size()); ++j) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(core.size());
        Eigen::VectorXd col;
        e(j) = 1.0;
        gcc_op.apply(e, col);
        gcc.col(j) = col;
    }
    gcc = 0.5 * (gcc + gcc.transpose());
    Eigen::MatrixXd gcn = 0.01 * Eigen::MatrixXd::Random(core.size(), noncore.size());
    Eigen::VectorXd gnn_diag = Eigen::VectorXd::Ones(noncore.size());

    ApyGinv reference;
    reference.build(gcc, gcn, gnn_diag, n, core, "operator-pcg-reference");
    auto pcg = makePCGSolver();
    GccSolverPCGOptions options;
    options.tolerance = 1e-11;
    options.max_iter = 100;
    options.preconditioner = "jacobi";
    pcg->compute(gcc_op, options);
    MatrixFreeApyInverseOperator matrix_free(
        pcg, matrixFactory(gcn, noncore, 2), gnn_diag, n, core, "operator-pcg");
    const Eigen::MatrixXd x = Eigen::MatrixXd::Random(n, 3);
    Eigen::MatrixXd expected(n, 3);
    for (int j = 0; j < x.cols(); ++j) expected.col(j) = reference.applyInverse(x.col(j));
    const double error = (matrix_free.apply(x) - expected).cwiseAbs().maxCoeff();
    std::printf("PCG backend max_abs=%.3e\n", error);
    CHECK(error < 1e-7, "PCG backend operator parity");
    APYMMEOperator mme_operator(&matrix_free, 1.25);
    const double mme_error =
        (mme_operator * x.col(0) - 1.25 * expected.col(0)).cwiseAbs().maxCoeff();
    std::printf("PCG backend MME adapter max_abs=%.3e\n", mme_error);
    CHECK(mme_error < 1e-7, "PCG backend MME adapter parity");
}

int main() {
    try {
        testDenseAndMatrixFreeParity();
        testCoreIsAll();
        testPcgBackend();
    } catch (const std::exception& e) {
        ++failures;
        std::printf("FATAL: %s\n", e.what());
    }
    std::printf("APY OPERATOR E3-C2: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
