#include <Eigen/Dense>

#include <cstdio>
#include <stdexcept>
#include <vector>

#include "a22_inverse_operator.h"
#include "design.h"
#include "h_inverse_operator.h"
#include "mme_builder.h"
#include "implicit_mme.h"
#include "matrix_adapter.h"
#include "matrix_free_pcg.h"
#include "pedigree_relationship.h"

using namespace cosmic;

namespace {

int failures = 0;
#define CHECK(c, m) do { if (!(c)) { std::printf("FAIL: %s\n", m); ++failures; } } while (0)

double maxAbs(const Eigen::MatrixXd& x) {
    return x.size() == 0 ? 0.0 : x.cwiseAbs().maxCoeff();
}

std::vector<PedigreeRecord> testPedigree() {
    // Deliberately not topologically ordered.
    return {
        {"7", "3", "5"}, {"2", "0", "0"}, {"5", "1", "2"},
        {"1", "0", "0"}, {"8", "4", "6"}, {"3", "1", "2"},
        {"6", "3", "4"}, {"4", "1", "0"}
    };
}

Eigen::MatrixXd denseA22Inverse(const Eigen::SparseMatrix<double>& ainv,
                                const std::vector<int>& typed) {
    const int n = ainv.rows();
    std::vector<int> typed_pos(n, -1), untyped;
    for (int i = 0; i < static_cast<int>(typed.size()); ++i) typed_pos[typed[i]] = i;
    for (int i = 0; i < n; ++i) if (typed_pos[i] < 0) untyped.push_back(i);

    const Eigen::MatrixXd dense = Eigen::MatrixXd(ainv);
    Eigen::MatrixXd r(typed.size(), typed.size());
    for (int i = 0; i < r.rows(); ++i)
        for (int j = 0; j < r.cols(); ++j) r(i, j) = dense(typed[i], typed[j]);
    if (untyped.empty()) return r;

    Eigen::MatrixXd p(untyped.size(), untyped.size());
    Eigen::MatrixXd q(untyped.size(), typed.size());
    for (int i = 0; i < p.rows(); ++i) {
        for (int j = 0; j < p.cols(); ++j) p(i, j) = dense(untyped[i], untyped[j]);
        for (int j = 0; j < q.cols(); ++j) q(i, j) = dense(untyped[i], typed[j]);
    }
    return r - q.transpose() * p.ldlt().solve(q);
}

void testMixedPedigree() {
    PedigreeRelationship pedigree;
    pedigree.load(testPedigree());
    const Eigen::SparseMatrix<double>& ainv = pedigree.buildAInverse();
    const std::vector<int> typed{0, 2, 4, 6};
    const Eigen::MatrixXd reference = denseA22Inverse(ainv, typed);
    A22InverseOperator a22(ainv, typed);

    Eigen::MatrixXd rhs(typed.size(), 3);
    rhs << 0.2, -0.4, 0.7,
           1.1,  0.3, 0.2,
          -0.5,  0.8, 0.1,
           0.9, -0.6, 1.3;
    const Eigen::MatrixXd matrix_actual = a22.apply(rhs);
    const Eigen::VectorXd vector_actual = a22.apply(rhs.col(1));
    const double matrix_error = maxAbs(matrix_actual - reference * rhs);
    const double vector_error = maxAbs(vector_actual - reference * rhs.col(1));
    std::printf("A22 mixed matrix max_abs=%.3e vector max_abs=%.3e residual=%.3e\n",
                matrix_error, vector_error, a22.maxSolveRelativeResidual());
    CHECK(matrix_error < 1e-12, "mixed-pedigree A22 inverse matrix RHS parity");
    CHECK(vector_error < 1e-12, "mixed-pedigree A22 inverse vector RHS parity");
    CHECK(a22.backendName() == "sparse-schur-a22-inverse", "mixed-pedigree backend metadata");
    CHECK(a22.pNonZeros() > 0 && a22.factorNonZeros() > 0,
          "mixed-pedigree sparse factor diagnostics");
    CHECK(a22.maxSolveRelativeResidual() < 1e-12, "mixed-pedigree solve residual");

    constexpr double tau = 1.25;
    constexpr double omega = 0.65;
    Eigen::MatrixXd b(typed.size(), typed.size());
    b << 2.0, 0.2, 0.1, 0.0,
         0.2, 1.7, 0.3, 0.1,
         0.1, 0.3, 1.9, 0.2,
         0.0, 0.1, 0.2, 1.6;
    const Eigen::MatrixXd ginv = b.ldlt().solve(Eigen::MatrixXd::Identity(b.rows(), b.cols()));
    DenseMatrixAdapter scaled_ginv(tau * ginv);
    HInverseOperator hinv(ainv, &scaled_ginv, &a22, typed, omega);

    Eigen::MatrixXd dense_h = Eigen::MatrixXd(ainv);
    const Eigen::MatrixXd correction = tau * ginv - omega * reference;
    for (int i = 0; i < static_cast<int>(typed.size()); ++i)
        for (int j = 0; j < static_cast<int>(typed.size()); ++j)
            dense_h(typed[i], typed[j]) += correction(i, j);
    const Eigen::VectorXd ped_rhs = Eigen::VectorXd::LinSpaced(ainv.rows(), -0.7, 1.1);
    const double h_error = maxAbs(hinv * ped_rhs - dense_h * ped_rhs);
    std::printf("H-inverse mixed matvec max_abs=%.3e diagonal=%s\n",
                h_error, hinv.preconditionerDiagonalMode().c_str());
    CHECK(h_error < 1e-12, "mixed-pedigree H inverse matvec parity");
    CHECK(!hinv.preconditionerDiagonalExact(), "mixed-pedigree diagonal is explicitly approximate");
    CHECK(hinv.preconditionerDiagonalMode() == "approx_a22_omitted",
          "mixed-pedigree preconditioner metadata");

    FixedDesignG design;
    design.p = 1;
    design.names = {"mu"};
    std::vector<GenRecord> records(ainv.rows());
    for (int i = 0; i < ainv.rows(); ++i) {
        records[i].aid = i + 1;
        records[i].y = 0.15 * i - 0.35;
        design.rows.push_back({{0, 1.0}});
    }
    RandomComponent component;
    component.Qinv = &hinv;
    constexpr double lambda = 0.55;
    ImplicitMME implicit(records, design, {component}, {lambda});

    Eigen::MatrixXd lhs = Eigen::MatrixXd::Zero(ainv.rows() + 1, ainv.rows() + 1);
    lhs(0, 0) = ainv.rows();
    lhs.block(0, 1, 1, ainv.rows()).setOnes();
    lhs.block(1, 0, ainv.rows(), 1).setOnes();
    lhs.block(1, 1, ainv.rows(), ainv.rows()) =
        Eigen::MatrixXd::Identity(ainv.rows(), ainv.rows()) + lambda * dense_h;
    Eigen::VectorXd mme_rhs(ainv.rows() + 1);
    mme_rhs(0) = 0.0;
    for (int i = 0; i < ainv.rows(); ++i) {
        mme_rhs(0) += records[i].y;
        mme_rhs(i + 1) = records[i].y;
    }
    const Eigen::VectorXd expected = lhs.ldlt().solve(mme_rhs);
    MatrixFreePCG pcg(mme_rhs.size());
    pcg.setQuiet(true);
    pcg.setMultOp([&](const Eigen::VectorXd& x, Eigen::VectorXd& out) { implicit.multiply(x, out); });
    pcg.setPrecondOp([&](const Eigen::VectorXd& x, Eigen::VectorXd& out) {
        out = x.cwiseQuotient(implicit.diag);
    });
    const Eigen::VectorXd actual = pcg.solve(mme_rhs, 1e-13, 2000);
    const double mme_error = maxAbs(actual - expected);
    std::printf("Implicit MME mixed solution max_abs=%.3e\n", mme_error);
    CHECK(mme_error < 1e-10, "operator H inverse implicit MME solution parity");
}

void testAllTyped() {
    PedigreeRelationship pedigree;
    pedigree.load(testPedigree());
    const Eigen::SparseMatrix<double>& ainv = pedigree.buildAInverse();
    std::vector<int> typed(ainv.rows());
    for (int i = 0; i < ainv.rows(); ++i) typed[i] = i;

    A22InverseOperator a22(ainv, typed);
    CHECK(a22.nUntyped() == 0, "all-typed A22 branch");
    CHECK(a22.backendName() == "sparse-schur-a22-inverse", "all-typed backend metadata");
    const Eigen::VectorXd rhs = Eigen::VectorXd::LinSpaced(ainv.rows(), -1.0, 0.8);
    CHECK(maxAbs(a22.apply(rhs) - ainv * rhs) < 1e-13,
          "all-typed A22 inverse equals A inverse");

    Eigen::MatrixXd genomic = 0.4 * Eigen::MatrixXd::Identity(ainv.rows(), ainv.rows());
    DenseMatrixAdapter ginv(genomic);
    HInverseOperator hinv(ainv, &ginv, &a22, typed, 0.3);
    Eigen::VectorXd diag = Eigen::VectorXd::Zero(ainv.rows());
    hinv.add_diagonal_to(diag, 1.0);
    const Eigen::VectorXd expected_diag = Eigen::MatrixXd(ainv).diagonal() +
        genomic.diagonal() - 0.3 * Eigen::MatrixXd(ainv).diagonal();
    CHECK(hinv.preconditionerDiagonalExact(), "all-typed H diagonal is exact");
    CHECK(hinv.preconditionerDiagonalMode() == "exact", "all-typed diagonal metadata");
    CHECK(maxAbs(diag - expected_diag) < 1e-13, "all-typed H diagonal parity");
}

} // namespace

int main() {
    try {
        testMixedPedigree();
        testAllTyped();
    } catch (const std::exception& e) {
        std::printf("FATAL: %s\n", e.what());
        ++failures;
    }
    std::printf("A22/HINV OPERATOR E4-C: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
