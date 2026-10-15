#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "models/apy_ginv.h"
#include "models/apy_core_pev_operator.h"
#include "models/apy_inverse_operator.h"
#include "models/apy_mme_operator.h"
#include "models/pev_operator.h"
#include "design.h"
#include "grm_builder.h"
#include "mme_builder.h"
#include "implicit_mme.h"
#include "plink_reader.h"

using namespace cosmic;

namespace {

int failures = 0;
#define CHECK(c, m) do { if (!(c)) { std::printf("FAIL: %s\n", m); ++failures; } } while (0)

double maxAbs(const Eigen::MatrixXd& value) {
    return value.size() ? value.cwiseAbs().maxCoeff() : 0.0;
}

void testSynthetic() {
    constexpr int n = 24;
    Eigen::MatrixXd random = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd coefficient = random * random.transpose();
    coefficient.diagonal().array() += 3.0;

    DensePEVOperator dense(coefficient);
    MatrixFreePEVOperator::Options options;
    options.tolerance = 1e-12;
    options.max_iterations = 500;
    MatrixFreePEVOperator matrix_free(
        n,
        [&coefficient](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            y.noalias() = coefficient * x;
        },
        [&coefficient](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            y = x.cwiseQuotient(coefficient.diagonal());
        },
        options);

    const Eigen::VectorXd rhs = Eigen::VectorXd::Random(n);
    const Eigen::VectorXd legacy_vector = coefficient.ldlt().solve(rhs);
    const Eigen::VectorXd dense_vector = dense.apply(rhs);
    const Eigen::VectorXd operator_vector = matrix_free.apply(rhs);
    const double vector_error = maxAbs(operator_vector - dense_vector);
    const double quadratic_error = std::fabs(
        matrix_free.quadraticForm(rhs) - dense.quadraticForm(rhs));

    const Eigen::MatrixXd rhs_matrix = Eigen::MatrixXd::Random(n, 5);
    const Eigen::MatrixXd legacy_matrix = coefficient.ldlt().solve(rhs_matrix);
    const Eigen::MatrixXd dense_matrix = dense.apply(rhs_matrix);
    const Eigen::MatrixXd operator_matrix = matrix_free.apply(rhs_matrix);
    const double matrix_error = maxAbs(operator_matrix - dense_matrix);
    const double dense_error = std::max(
        maxAbs(dense_vector - legacy_vector), maxAbs(dense_matrix - legacy_matrix));

    std::printf(
        "Synthetic PEV dense=%.3e apply=%.3e matrix=%.3e quadratic=%.3e residual=%.3e\n",
        dense_error, vector_error, matrix_error, quadratic_error,
        matrix_free.lastSolveInfo().relative_residual);
    CHECK(dense_error < 1e-12, "dense interface parity with legacy LDLT solve");
    CHECK(vector_error < 1e-10, "synthetic vector inverse-apply parity");
    CHECK(matrix_error < 1e-10, "synthetic matrix inverse-apply parity");
    CHECK(quadratic_error < 1e-10, "synthetic quadratic parity");
    CHECK(matrix_free.lastSolveInfo().converged, "synthetic PCG convergence");
    CHECK(matrix_free.lastBatchSolveInfo().size() == 5, "synthetic batch diagnostics");

    const Eigen::VectorXd small_rhs = 1e-12 * rhs;
    const Eigen::VectorXd small_solution = matrix_free.apply(small_rhs);
    const double small_relative_residual =
        (small_rhs - coefficient * small_solution).norm() / small_rhs.norm();
    CHECK(small_relative_residual <= options.tolerance,
          "small RHS uses a true relative residual");

    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(n);
    CHECK(matrix_free.apply(zero).norm() == 0.0, "zero RHS solution");
    CHECK(matrix_free.lastSolveInfo().termination == "zero_rhs", "zero RHS diagnostics");

    bool rejected_indefinite = false;
    MatrixFreePEVOperator indefinite(
        2,
        [](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            y.resize(2);
            y << x(0), -x(1);
        });
    try {
        Eigen::VectorXd bad_rhs(2);
        bad_rhs << 0.0, 1.0;
        (void)indefinite.apply(bad_rhs);
    } catch (const std::runtime_error&) {
        rejected_indefinite = true;
    }
    CHECK(rejected_indefinite, "matrix-free PEV rejects non-SPD curvature");
}

void testSyntheticApyCore() {
    constexpr int n = 18;
    constexpr int nc = 7;
    Eigen::MatrixXd random = Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd g = random * random.transpose();
    g.diagonal().array() += 4.0;
    std::vector<int> core(nc);
    std::vector<int> noncore(n - nc);
    for (int i = 0; i < nc; ++i) core[i] = i;
    for (int i = nc; i < n; ++i) noncore[i - nc] = i;

    ApyGinv apy;
    apy.build(g, core, "E3-D2C synthetic APY core PEV");
    Eigen::VectorXd dc = Eigen::VectorXd::LinSpaced(nc, 0.0, 2.0);
    Eigen::VectorXd dn = Eigen::VectorXd::LinSpaced(n - nc, 0.0, 1.5);
    constexpr double lambda = 1.7;
    DenseAPYCorePEVOperator dense(apy, dc, dn, lambda);
    MatrixFreePEVOperator::Options options;
    options.tolerance = 1e-12;
    options.max_iterations = 1000;
    MatrixFreeAPYCorePEVOperator matrix_free(
        apy.gccSolverPtr(), makeDenseGcnBlockIteratorFactory(g, core, noncore, 3),
        apy.mnn(), dc, dn, lambda, options);

    const Eigen::VectorXd probe = Eigen::VectorXd::Random(nc);
    Eigen::VectorXd dense_product;
    Eigen::VectorXd operator_product;
    dense.applyCoefficient(probe, dense_product);
    matrix_free.applyCoefficient(probe, operator_product);
    const double coefficient_error = maxAbs(operator_product - dense_product);
    const double apply_error = maxAbs(matrix_free.apply(probe) - dense.apply(probe));
    const double quadratic_error = std::fabs(
        matrix_free.quadraticForm(probe) - dense.quadraticForm(probe));
    std::printf(
        "Synthetic APY core PEV coefficient=%.3e apply=%.3e quadratic=%.3e residual=%.3e\n",
        coefficient_error, apply_error, quadratic_error,
        matrix_free.lastSolveInfo().relative_residual);
    CHECK(coefficient_error < 1e-11, "synthetic APY core coefficient parity");
    CHECK(apply_error < 1e-10, "synthetic APY core inverse-apply parity");
    CHECK(quadratic_error < 1e-10, "synthetic APY core quadratic parity");
    CHECK(matrix_free.lastSolveInfo().converged, "synthetic APY core PCG convergence");

    bool rejected_nonrestartable = false;
    const auto dense_factory = makeDenseGcnBlockIteratorFactory(g, core, noncore, 3);
    const auto factory_calls = std::make_shared<int>(0);
    GcnBlockIteratorFactory nonrestartable = [dense_factory, factory_calls]() {
        if ((*factory_calls)++ == 0) return dense_factory();
        return GcnBlockIteratorPtr{};
    };
    try {
        MatrixFreeAPYCorePEVOperator invalid_stream(
            apy.gccSolverPtr(), nonrestartable, apy.mnn(), dc, dn, lambda,
            options);
    } catch (const std::runtime_error&) {
        rejected_nonrestartable = true;
    }
    CHECK(rejected_nonrestartable, "APY core PEV rejects non-restartable Gcn streams");

    bool rejected_nonconvergence = false;
    MatrixFreePEVOperator::Options insufficient_options;
    insufficient_options.tolerance = 1e-16;
    insufficient_options.max_iterations = 1;
    MatrixFreeAPYCorePEVOperator insufficient(
        apy.gccSolverPtr(), makeDenseGcnBlockIteratorFactory(g, core, noncore, 3),
        apy.mnn(), dc, dn, lambda, insufficient_options);
    try {
        (void)insufficient.apply(probe);
    } catch (const std::runtime_error&) {
        rejected_nonconvergence = true;
    }
    CHECK(rejected_nonconvergence, "APY core PEV rejects PCG nonconvergence");

    std::vector<int> full_core(n);
    for (int i = 0; i < n; ++i) full_core[i] = i;
    ApyGinv full_apy;
    full_apy.build(g, full_core, "E3-D2C full-core APY PEV");
    const Eigen::VectorXd full_dc = Eigen::VectorXd::Zero(n);
    const Eigen::VectorXd empty;
    DenseAPYCorePEVOperator full_dense(full_apy, full_dc, empty, lambda);
    MatrixFreeAPYCorePEVOperator full_operator(
        full_apy.gccSolverPtr(), {}, full_apy.mnn(), full_dc, empty, lambda,
        options);
    const Eigen::VectorXd full_probe = Eigen::VectorXd::Random(n);
    CHECK(maxAbs(full_operator.apply(full_probe) - g * full_probe / lambda) < 1e-10,
          "full-core no-data limit Kc^-1 = Gcc/lambda");
    CHECK(maxAbs(full_operator.apply(full_probe) - full_dense.apply(full_probe)) < 1e-10,
          "full-core APY core PEV parity");
}

void testFrozen44(const std::string& prefix) {
    PlinkReader reader;
    reader.load(prefix);
    reader.computeStats(1);
    const int n = static_cast<int>(reader.getNumSamples());
    CHECK(n == 44, "frozen fixture must contain 44 animals");
    if (n != 44) return;

    GrmBuilder::Options grm_options;
    grm_options.algorithm = 1;
    grm_options.min_maf = 0.01f;
    grm_options.min_hwe = 0.0f;
    grm_options.thread_num = 1;
    GrmBuilder grm_builder(grm_options);
    grm_builder.compute(reader);
    Eigen::MatrixXd g = 0.95 * grm_builder.getGrm();
    g.diagonal().array() += 0.05;

    std::vector<int> core(28);
    for (int i = 0; i < static_cast<int>(core.size()); ++i) core[i] = i;
    std::vector<int> noncore;
    for (int i = static_cast<int>(core.size()); i < n; ++i) noncore.push_back(i);

    ApyGinv dense_apy;
    dense_apy.build(g, core, "E3-D2B frozen dense APY");
    DenseApyInverseOperator dense_inverse(dense_apy);

    Eigen::MatrixXd gcc(core.size(), core.size());
    Eigen::VectorXd gnn_diag(noncore.size());
    for (int i = 0; i < gcc.rows(); ++i)
        for (int j = 0; j < gcc.cols(); ++j) gcc(i, j) = g(core[i], core[j]);
    for (int j = 0; j < gnn_diag.size(); ++j) gnn_diag(j) = g(noncore[j], noncore[j]);
    MatrixFreeApyInverseOperator streamed_inverse(
        makeDenseSolver(gcc), makeDenseGcnBlockIteratorFactory(g, core, noncore, 7),
        gnn_diag, n, core, "E3-D2B frozen APY operator");

    Eigen::VectorXd dc(core.size());
    Eigen::VectorXd dn(noncore.size());
    for (int i = 0; i < dc.size(); ++i) dc(i) = (i % 5 == 0) ? 0.0 : 0.4 + 0.03 * i;
    for (int i = 0; i < dn.size(); ++i) dn(i) = (i % 4 == 0) ? 0.0 : 0.2 + 0.02 * i;
    MatrixFreePEVOperator::Options core_options;
    core_options.tolerance = 1e-12;
    core_options.max_iterations = 2000;
    DenseAPYCorePEVOperator dense_core_pev(dense_apy, dc, dn, 2.0);
    MatrixFreeAPYCorePEVOperator operator_core_pev(
        streamed_inverse.gccSolverPtr(),
        makeDenseGcnBlockIteratorFactory(g, core, noncore, 7),
        streamed_inverse.mnn(), dc, dn, 2.0, core_options);
    const Eigen::MatrixXd core_probes = Eigen::MatrixXd::Random(core.size(), 4);
    const double core_pev_error =
        maxAbs(operator_core_pev.apply(core_probes) - dense_core_pev.apply(core_probes));
    double core_quadratic_error = 0.0;
    for (Eigen::Index j = 0; j < core_probes.cols(); ++j) {
        core_quadratic_error = std::max(
            core_quadratic_error,
            std::fabs(operator_core_pev.quadraticForm(core_probes.col(j)) -
                      dense_core_pev.quadraticForm(core_probes.col(j))));
    }
    std::printf("Frozen44 APY core PEV apply=%.3e quadratic=%.3e\n",
                core_pev_error, core_quadratic_error);
    CHECK(core_pev_error < 1e-9, "frozen APY core inverse-apply parity");
    CHECK(core_quadratic_error < 1e-9, "frozen APY core quadratic parity");

    ApyGinvAdapter dense_apy_matrix(&dense_apy, 1.0);
    APYMMEOperator apy_mme_matrix(&streamed_inverse, 1.0);

    FixedDesignG fixed_design;
    fixed_design.p = 1;
    fixed_design.names = {"mu"};
    std::vector<GenRecord> records(n);
    for (int i = 0; i < n; ++i) {
        records[i].aid = i + 1;
        records[i].y = std::sin(0.17 * i) + 0.01 * i;
        fixed_design.rows.push_back({{0, 1.0}});
    }

    RandomComponent dense_component;
    dense_component.Qinv = &dense_apy_matrix;
    constexpr double lambda = 2.0;
    MMELHSBuilder builder(records, fixed_design, {dense_component}, true);
    const Eigen::SparseMatrix<double> coefficient = builder.build_lhs({lambda});
    DensePEVOperator dense_pev(coefficient);

    RandomComponent operator_component;
    operator_component.Qinv = &apy_mme_matrix;
    ImplicitMME implicit(records, fixed_design, {operator_component}, {lambda});
    MatrixFreePEVOperator::Options options;
    options.tolerance = 1e-11;
    options.max_iterations = 2000;
    MatrixFreePEVOperator operator_pev(
        implicit.dim,
        [&implicit](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            implicit.multiply(x, y);
        },
        [&implicit](const Eigen::VectorXd& x, Eigen::VectorXd& y) {
            y = x.cwiseQuotient(implicit.diag);
        },
        options);

    Eigen::MatrixXd probes = Eigen::MatrixXd::Random(implicit.dim, 4);
    const double apply_error = maxAbs(operator_pev.apply(probes) - dense_pev.apply(probes));
    double quadratic_error = 0.0;
    for (Eigen::Index j = 0; j < probes.cols(); ++j) {
        quadratic_error = std::max(
            quadratic_error,
            std::fabs(operator_pev.quadraticForm(probes.col(j)) -
                      dense_pev.quadraticForm(probes.col(j))));
    }

    Eigen::MatrixXd selected_rhs = Eigen::MatrixXd::Zero(implicit.dim, n);
    for (int i = 0; i < n; ++i) selected_rhs(1 + i, i) = 1.0;
    const Eigen::MatrixXd dense_selected = dense_pev.apply(selected_rhs);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> legacy_factor(coefficient);
    CHECK(legacy_factor.info() == Eigen::Success, "frozen legacy LDLT factorization");
    const Eigen::MatrixXd legacy_selected = legacy_factor.solve(selected_rhs);
    const Eigen::MatrixXd operator_selected = operator_pev.apply(selected_rhs);
    constexpr double sigma_e2 = 2.0;
    Eigen::VectorXd dense_pev_values(n);
    Eigen::VectorXd operator_pev_values(n);
    Eigen::VectorXd dense_prediction_error(n);
    Eigen::VectorXd operator_prediction_error(n);
    for (int i = 0; i < n; ++i) {
        dense_pev_values(i) = sigma_e2 * dense_selected(1 + i, i);
        operator_pev_values(i) = sigma_e2 * operator_selected(1 + i, i);
        dense_prediction_error(i) = g(i, i) - dense_pev_values(i);
        operator_prediction_error(i) = g(i, i) - operator_pev_values(i);
    }
    const double pev_error = maxAbs(operator_pev_values - dense_pev_values);
    const double dense_error = maxAbs(dense_selected - legacy_selected);
    const double prediction_error =
        maxAbs(operator_prediction_error - dense_prediction_error);
    double max_residual = 0.0;
    for (const auto& info : operator_pev.lastBatchSolveInfo()) {
        max_residual = std::max(max_residual, info.relative_residual);
        CHECK(info.converged, "frozen selected PEV solve convergence");
    }

    std::printf(
        "Frozen44 PEV dense=%.3e apply=%.3e quadratic=%.3e PEV=%.3e prediction_error=%.3e residual=%.3e\n",
        dense_error, apply_error, quadratic_error, pev_error, prediction_error, max_residual);
    CHECK(dense_error < 1e-13, "frozen dense interface parity with legacy selected solves");
    CHECK(apply_error < 1e-8, "frozen inverse-apply parity");
    CHECK(quadratic_error < 1e-8, "frozen quadratic parity");
    CHECK(pev_error < 1e-8, "frozen animal PEV parity");
    CHECK(prediction_error < 1e-8, "frozen prediction-error parity");
}

} // namespace

int main(int argc, char** argv) {
    try {
        testSynthetic();
        testSyntheticApyCore();
        if (argc > 1) testFrozen44(argv[1]);
    } catch (const std::exception& error) {
        std::printf("FATAL: %s\n", error.what());
        ++failures;
    }
    std::printf("PEV OPERATORS E3-D2B/D2C: %s\n", failures == 0 ? "PASS" : "FAIL");
    return failures == 0 ? 0 : 1;
}
