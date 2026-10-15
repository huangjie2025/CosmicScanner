#include <Eigen/Dense>

#include <cstdio>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "apy_ginv.h"
#include "apy_inverse_operator.h"
#include "apy_mme_operator.h"
#include "design.h"
#include "mme_builder.h"
#include "implicit_mme.h"
#include "matrix_adapter.h"
#include "matrix_free_pcg.h"

using namespace cosmic;

namespace {

int failures = 0;
#define CHECK(c, m) do { if (!(c)) { std::printf("FAIL: %s\n", m); ++failures; } } while (0)

Eigen::MatrixXd randomSpd(int n) {
    Eigen::MatrixXd a = Eigen::MatrixXd::Random(n, n);
    return a * a.transpose() + 2.0 * Eigen::MatrixXd::Identity(n, n);
}

std::vector<int> complement(int n, const std::vector<int>& core) {
    std::vector<int> out;
    for (int i = 0; i < n; ++i) {
        if (std::find(core.begin(), core.end(), i) == core.end()) out.push_back(i);
    }
    return out;
}

void runTest() {
    constexpr int ng = 9;
    constexpr int nped = 12;
    const std::vector<int> core{0, 2, 5, 7};
    const std::vector<int> noncore = complement(ng, core);
    const std::vector<int> map{0, 1, 3, 4, 6, 7, 8, 10, 11};
    const Eigen::MatrixXd g = randomSpd(ng);

    ApyGinv reference;
    reference.build(g, core, "E3-D1 reference");

    Eigen::MatrixXd gcc(core.size(), core.size());
    Eigen::VectorXd gnn_diag(noncore.size());
    for (int i = 0; i < gcc.rows(); ++i)
        for (int j = 0; j < gcc.cols(); ++j) gcc(i, j) = g(core[i], core[j]);
    for (int j = 0; j < gnn_diag.size(); ++j) gnn_diag(j) = g(noncore[j], noncore[j]);

    MatrixFreeApyInverseOperator inverse(
        makeDenseSolver(gcc), makeDenseGcnBlockIteratorFactory(g, core, noncore, 2),
        gnn_diag, ng, core, "E3-D1 operator");
    APYMMEOperator operator_matrix(&inverse, 1.3);

    const Eigen::VectorXd x = Eigen::VectorXd::Random(ng);
    const double apply_error =
        (operator_matrix * x - 1.3 * reference.applyInverse(x)).cwiseAbs().maxCoeff();
    Eigen::VectorXd expected_diag = Eigen::VectorXd::Zero(ng);
    Eigen::VectorXd actual_diag = Eigen::VectorXd::Zero(ng);
    reference.addDiagonalTo(expected_diag, 1.3);
    operator_matrix.add_diagonal_to(actual_diag, 1.0);
    const double diag_error = (actual_diag - expected_diag).cwiseAbs().maxCoeff();
    std::printf("APY MME apply max_abs=%.3e diagonal max_abs=%.3e\n", apply_error, diag_error);
    CHECK(apply_error < 1e-11, "APY MME operator matvec parity");
    CHECK(diag_error < 1e-11, "APY MME operator diagonal parity");

    Eigen::SparseMatrix<double> identity(nped, nped);
    identity.setIdentity();
    SparseMatrixAdapter ainv(identity);
    DenseMatrixAdapter a22(0.2 * Eigen::MatrixXd::Identity(ng, ng));
    SplitMatrixAdapter hinv(&ainv, &operator_matrix, &a22, map);

    Eigen::MatrixXd h_dense = Eigen::MatrixXd::Identity(nped, nped);
    for (int j = 0; j < ng; ++j) {
        Eigen::VectorXd e = Eigen::VectorXd::Zero(ng);
        e(j) = 1.0;
        const Eigen::VectorXd col = 1.3 * reference.applyInverse(e);
        for (int i = 0; i < ng; ++i) h_dense(map[i], map[j]) += col(i);
        h_dense(map[j], map[j]) -= 0.2;
    }
    const Eigen::VectorXd ped_x = Eigen::VectorXd::Random(nped);
    const double split_error = (hinv * ped_x - h_dense * ped_x).cwiseAbs().maxCoeff();
    Eigen::VectorXd split_diag = Eigen::VectorXd::Zero(nped);
    hinv.add_diagonal_to(split_diag, 1.0);
    CHECK(split_error < 1e-11, "mapped split H-inverse matvec parity");
    CHECK((split_diag - h_dense.diagonal()).cwiseAbs().maxCoeff() < 1e-11,
          "mapped split H-inverse diagonal parity");

    FixedDesignG fd;
    fd.p = 1;
    fd.names = {"mu"};
    std::vector<GenRecord> records(nped);
    for (int i = 0; i < nped; ++i) {
        records[i].aid = i + 1;
        records[i].y = 0.2 * i - 0.4;
        fd.rows.push_back({{0, 1.0}});
    }
    RandomComponent component;
    component.Qinv = &hinv;
    constexpr double lambda = 0.75;
    ImplicitMME implicit(records, fd, {component}, {lambda});

    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(nped + 1, nped + 1);
    lhs(0, 0) = nped;
    lhs.block(0, 1, 1, nped).setOnes();
    lhs.block(1, 0, nped, 1).setOnes();
    lhs.block(1, 1, nped, nped) = Eigen::MatrixXd::Identity(nped, nped) + lambda * h_dense;
    Eigen::VectorXd rhs(nped + 1);
    rhs(0) = 0.0;
    for (const auto& record : records) rhs(0) += record.y;
    for (int i = 0; i < nped; ++i) rhs(i + 1) = records[i].y;
    const Eigen::VectorXd expected_solution = lhs.ldlt().solve(rhs);

    MatrixFreePCG pcg(nped + 1);
    pcg.setMultOp([&](const Eigen::VectorXd& v, Eigen::VectorXd& out) { implicit.multiply(v, out); });
    pcg.setPrecondOp([&](const Eigen::VectorXd& v, Eigen::VectorXd& out) {
        out = v.cwiseQuotient(implicit.diag);
    });
    const Eigen::VectorXd actual_solution = pcg.solve(rhs, 1e-12, 2000);
    const double solve_error = (actual_solution - expected_solution).cwiseAbs().maxCoeff();
    std::printf("Synthetic implicit MME solve max_abs=%.3e\n", solve_error);
    CHECK(solve_error < 1e-9, "synthetic MME solve parity");
}

} // namespace

int main() {
    try {
        runTest();
    } catch (const std::exception& e) {
        std::printf("FATAL: %s\n", e.what());
        ++failures;
    }
    std::printf("APY MME OPERATOR E3-D1: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
