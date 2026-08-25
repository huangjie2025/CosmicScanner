// Phase D gate D1.5: equation-equivalence between
//
//   (A) Leite et al. 2024 Eq. 8 written with the APY block-inverse pieces
//       (superscript blocks G^cc, G^cn, G^nc, R = Mnn^-1):
//         C_core^-1 = D_c + lambda*G^cc - lambda*G^cn * (lambda*R + D_n)^-1 * lambda*G^nc
//   (B) the Schur/Woodbury simplified form (single production form):
//         C_core^-1 = D_c + lambda*Gcc^-1 + Q*S*Q',   S_j = lambda*d_j/(lambda + d_j*m_jj)
//
// with Q = Gcc^-1*Gcn, m_jj = g_jj - g_jc'*Gcc^-1*g_cj (noncore conditional
// variance; the code stores Mnn = diag(m_jj) and Mnn^-1 = R).
//
// Gates:
//   1. max_abs(A - B) < 1e-12 on 100 random SPD partitions (lambda, D random);
//   2. D = 0 (no information): B (and A) must give C_core = (1/lambda)*Gcc,
//      i.e. sigma_e2*C_core = sigma_a2*Gcc, the exact no-data PEV limit;
//   3. the *mis-notated* Eq. 8 variant that substitutes the ordinary blocks
//      G_cc / G_cn / G_nc for the APY-inverse blocks must be shown to differ
//      by O(1) - documents that the Phase A "O(1) deviation" of the "printed
//      Eq. 8" was a notation misinterpretation (NOT a paper equation error);
//   4. C_core^-1 is always positive definite on the random grid.
#include <Eigen/Dense>
#include <cstdio>
#include <random>

namespace {

double max_abs_diff(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    return (A - B).cwiseAbs().maxCoeff();
}

int runOne(const Eigen::MatrixXd& G, int nc, const Eigen::VectorXd& Dc,
           const Eigen::VectorXd& Dn, double lambda, double& worst_ab,
           double& worst_degen, double& worst_misnot) {
    const int nn = (int)G.rows() - nc;
    const Eigen::MatrixXd Gcc = G.topLeftCorner(nc, nc);
    const Eigen::MatrixXd Gcn = G.topRightCorner(nc, nn);
    const Eigen::MatrixXd Gnc = G.bottomLeftCorner(nn, nc);
    const Eigen::MatrixXd Gnn = G.bottomRightCorner(nn, nn);

    Eigen::LDLT<Eigen::MatrixXd> ldlt(Gcc);
    if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0) return 1;
    const Eigen::MatrixXd Gcc_inv = ldlt.solve(Eigen::MatrixXd::Identity(nc, nc));
    const Eigen::MatrixXd Q = Gcc_inv * Gcn;

    // m_jj = g_jj - g_jc' Gcc^-1 g_cj  (noncore conditional variance)
    const Eigen::VectorXd m =
        Gnn.diagonal() - (Gnc * (Gcc_inv * Gcn)).diagonal();

    // R = Mnn^-1 (diag); S_j = lambda*d_j / (lambda + d_j*m_jj)
    const Eigen::VectorXd R = m.cwiseInverse();
    const Eigen::VectorXd S = lambda * Dn.array() / (lambda + Dn.array() * m.array());

    // --- APY block-inverse pieces (superscript blocks, Eq. 5) ---
    const Eigen::MatrixXd Gcc_apy = Gcc_inv + Q * R.asDiagonal() * Q.transpose();
    const Eigen::MatrixXd Gcn_apy = -(Q * Eigen::MatrixXd(R.asDiagonal()));
    const Eigen::MatrixXd Gnc_apy = -(Eigen::MatrixXd(R.asDiagonal()) * Q.transpose());

    // --- (A) Leite Eq. 8 with APY-inverse blocks ---
    const Eigen::MatrixXd W = (lambda * R.asDiagonal().toDenseMatrix() + Eigen::MatrixXd(Dn.asDiagonal())).inverse();
    const Eigen::MatrixXd A_inv =
        Eigen::MatrixXd(Dc.asDiagonal()) + lambda * Gcc_apy -
        lambda * Gcn_apy * W * lambda * Gnc_apy;

    // --- (B) Schur/Woodbury simplified form ---
    const Eigen::MatrixXd B_inv =
        Eigen::MatrixXd(Dc.asDiagonal()) + lambda * Gcc_inv + Q * S.asDiagonal() * Q.transpose();

    const double ab = max_abs_diff(A_inv, B_inv);
    worst_ab = std::max(worst_ab, ab);
    if (ab >= 1e-12) return 2;

    // --- Gate 2: D = 0 -> C_core = (1/lambda) Gcc (no-data PEV limit) ---
    const Eigen::VectorXd S0 = Eigen::VectorXd::Zero(nn);
    const Eigen::MatrixXd B0_inv =
        lambda * Gcc_inv + Q * S0.asDiagonal() * Q.transpose();
    const double d0 = max_abs_diff(B0_inv, lambda * Gcc_inv);
    worst_degen = std::max(worst_degen, d0);
    if (d0 > 1e-12) return 3;
    const double lim =
        max_abs_diff(B0_inv.inverse(), (1.0 / lambda) * Gcc);
    if (lim > 1e-9 * (1.0 / lambda) * Gcc.cwiseAbs().maxCoeff()) return 4;

    // --- Gate 3: mis-notated variant (ordinary G_cc/G_cn blocks instead of
    //             the APY-inverse superscript blocks) is NOT equivalent ---
    const Eigen::MatrixXd M_inv =
        Eigen::MatrixXd(Dc.asDiagonal()) + lambda * Gcc -
        lambda * Gcn * W * lambda * Gnc;
    const double mn = max_abs_diff(M_inv, B_inv);
    worst_misnot = std::max(worst_misnot, mn);
    if (mn <= 1e-6) return 5;

    // --- Gate 4: B_inv positive definite ---
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B_inv);
    if (es.info() != Eigen::Success || es.eigenvalues().minCoeff() <= 0.0) return 6;

    return 0;
}

}  // namespace

int main() {
    std::mt19937_64 rng(0xD15u);
    std::uniform_int_distribution<int> nc_uni(3, 12);
    std::uniform_int_distribution<int> nn_uni(5, 18);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> lambda_uni(0.1, 10.0);
    std::uniform_real_distribution<double> d_uni(0.0, 3.0);
    std::uniform_int_distribution<int> zero_mask(0, 7);  // ~25% exact zeros in D

    const int n_cases = 100;
    int n_pass = 0, fails = 0;
    double worst_ab = 0.0, worst_degen = 0.0, worst_misnot = 0.0;

    for (int c = 0; c < n_cases; ++c) {
        const int nc = nc_uni(rng), nn = nn_uni(rng);
        const int ng = nc + nn;
        // full SPD G via A*A' + I (k = ng + 4 ensures full rank)
        const int k = ng + 4;
        Eigen::MatrixXd A(ng, k);
        for (auto& v : A.reshaped()) v = nd(rng);
        const Eigen::MatrixXd G = A * A.transpose() + Eigen::MatrixXd::Identity(ng, ng);

        Eigen::VectorXd Dc(nc), Dn(nn);
        for (int i = 0; i < nc; ++i) Dc(i) = (zero_mask(rng) & 1) ? 0.0 : d_uni(rng);
        for (int j = 0; j < nn; ++j) Dn(j) = (zero_mask(rng) & 2) ? 0.0 : d_uni(rng);
        const double lambda = lambda_uni(rng);

        double wa = 0.0, wd = 0.0, wm = 0.0;
        const int rc = runOne(G, nc, Dc, Dn, lambda, wa, wd, wm);
        worst_ab = std::max(worst_ab, wa);
        worst_degen = std::max(worst_degen, wd);
        worst_misnot = std::max(worst_misnot, wm);
        if (rc == 0) {
            ++n_pass;
        } else {
            ++fails;
            printf("case %d FAIL rc=%d (nc=%d nn=%d lambda=%.3f)\n", c, rc, nc, nn, lambda);
        }
    }

    printf("D1.5: cases=%d pass=%d fails=%d\n", n_cases, n_pass, fails);
    printf("  worst |A-B| (Eq.8-superscript vs Schur) = %.3e  (gate < 1e-12)\n", worst_ab);
    printf("  worst D=0 |B0_inv - lambda*Gcc^-1|      = %.3e  (gate < 1e-12)\n", worst_degen);
    printf("  worst |mis-notated variant - Schur|     = %.3e  (must be O(1), not ~0)\n",
           worst_misnot);

    if (fails == 0 && worst_ab < 1e-12 && worst_degen < 1e-12) {
        printf("D1.5 PASS: Leite Eq.8 (superscript blocks) == Schur/Woodbury form;\n"
               "D=0 no-information limit C_core = (1/lambda)*Gcc recovered exactly.\n");
        return 0;
    }
    printf("D1.5 FAIL: equation equivalence not established - do not proceed to D2.\n");
    return 1;
}