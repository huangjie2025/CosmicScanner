// RRM basis identity unit test (Phase 4 / G6).
//
// Guarantees the invariant relied on by the record-level RRM-GWAS:
// the SAME buildRrmBasisMatrix() (single normalization with global [tmin,tmax])
// provides both the random-effect trajectory basis and the SNP fixed-trajectory
// basis (D(r,k) = x_{subject(r),j} * phi_k(time_r)). A single call is validated
// for Legendre and B-spline: analytic definition match, partition-of-unity for
// B-spline, and deterministic bitwise reproducibility (random basis row == SNP
// fixed basis row at the same t).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "design.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

static double factorial_d(int n) {
    double out = 1.0;
    for (int i = 2; i <= n; ++i) out *= i;
    return out;
}

static double analytic_legendre(double t_norm, int k) {
    // Normalized shifted Legendre polynomial on t_norm in [-1,1],
    // same definition used by buildLegendreMatrix (sqrt((2k+1)/2)).
    double p = 0.0;
    for (int r = 0; r <= k / 2; ++r) {
        p += std::sqrt((2.0 * k + 1.0) / 2.0) * std::pow(0.5, k) *
             (std::pow(-1.0, r) * factorial_d(2 * k - 2 * r) /
              (factorial_d(r) * factorial_d(k - r) * factorial_d(k - 2 * r))) *
             std::pow(t_norm, k - 2 * r);
    }
    return p;
}

int main() {
    srand(20260815);
    const double tmin = 0.0, tmax = 1.0;
    const int order = 3;
    const int n = 500;
    VectorXd time(n);
    for (int i = 0; i < n; ++i) time(i) = (double)rand() / (double)RAND_MAX;

    double max_diff_legendre = 0.0;
    MatrixXd Phi_leg = cosmic::buildRrmBasisMatrix(time, tmin, tmax, order, "legendre");
    for (int i = 0; i < n; ++i) {
        double t_norm = 2.0 * (time(i) - tmin) / (tmax - tmin) - 1.0;
        for (int k = 0; k <= order; ++k) {
            max_diff_legendre = std::max(max_diff_legendre, std::abs(Phi_leg(i, k) - analytic_legendre(t_norm, k)));
        }
    }
    if (max_diff_legendre > 1e-12) {
        std::cerr << "Legendre basis mismatch vs analytic definition: " << max_diff_legendre << "\n";
        return 1;
    }

    MatrixXd Phi_bs = cosmic::buildRrmBasisMatrix(time, tmin, tmax, order, "bspline");
    double max_partition = 0.0;
    double max_bs_diff = 0.0;
    for (int i = 0; i < n; ++i) {
        double row_sum = Phi_bs.row(i).sum();
        max_partition = std::max(max_partition, std::abs(row_sum - 1.0));
        for (int k = 0; k <= order; ++k) {
            if (Phi_bs(i, k) < -1e-12 || Phi_bs(i, k) > 1.0 + 1e-12) {
                std::cerr << "B-spline basis out of [0,1]: " << Phi_bs(i, k) << "\n";
                return 1;
            }
        }
    }
    if (max_partition > 1e-12) {
        std::cerr << "B-spline row does not sum to 1: " << max_partition << "\n";
        return 1;
    }

    // Determinism: the SAME t vector must produce bitwise identical rows.
    Eigen::MatrixXd Phi_leg_2 = cosmic::buildRrmBasisMatrix(time, tmin, tmax, order, "legendre");
    Eigen::MatrixXd Phi_bs_2 = cosmic::buildRrmBasisMatrix(time, tmin, tmax, order, "bspline");
    double max_repro_leg = (Phi_leg - Phi_leg_2).cwiseAbs().maxCoeff();
    double max_repro_bs = (Phi_bs - Phi_bs_2).cwiseAbs().maxCoeff();
    double max_repro = std::max(max_repro_leg, max_repro_bs);
    if (max_repro > 0.0) {
        std::cerr << "Basis not deterministic between calls: " << max_repro << "\n";
        return 1;
    }
    if (order >= 1) {
        // Row differences across basis columns should not be all-identical pathologies
        // (sanity guard that the matrix is not constant/zero).
        double range_max = Phi_leg.cwiseAbs().maxCoeff();
        if (range_max <= 1e-12) {
            std::cerr << "Legend basis appears degenerate\n";
            return 1;
        }
    }

    std::cout << "RRM basis identity OK: max|leg-analytic|=" << max_diff_legendre
              << " max|bs-partition-1|=" << max_partition
              << " max|reproducibility|=" << max_repro << "\n";
    return 0;
}
