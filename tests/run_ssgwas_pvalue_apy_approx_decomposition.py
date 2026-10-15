#!/usr/bin/env python3
"""Phase D gate D4 (tiny): four-layer uncertainty-error decomposition.

Reference tower (diagnostic-only, NOT exposed via CLI):
  R0  = exact C1: dense G^-1 in the lambda-form MME, full (C_lambda^-1)_gg
        PEV, full operator w = G_op^-1 M_j            -> C++ "exact"
  R1  = apy-exact: APY G^-1 in the lambda-form MME, exact core slice
        (C_lambda^-1)_cc of the ACTUAL MME, core operator w_c
        -> C++ "apy-exact" (this IS the exact core Schur complement of the
        actual MME for the core block; R1 == R2-reference by construction,
        noted in the report)
  R2  = Eq. 8 of Leite et al. 2024 with the FULL absorption matrix
        D_full = W'(I - X(X'X)^-1 X')W (full blocks; block Schur inversion) -- Python
  R3  = production apy-approx: Eq. 8 with Eq. 9 diagonal D (diag(D_full)),
        S_j = lambda d_j/(lambda + d_j m_jj)          -> C++ "apy-approx"

Error decomposition:
  R0->R1 : APY G inverse + core-representation error (per core size)
  R1->R2 : structural model error (actual MME vs (D + lambda G_APY^-1):
           Ainv coupling, tau/omega structure, ungenotyped records)
  R2->R3 : diagonal-D approximation (off-diagonal absorption loss)
  B==C gate: production Eq. 9 diagonal D == diag(D_full) (single-record data)
  operator_full_vs_core: same PEV, full w vs core w_c (Bermann Eq. 21)

Facets:
  * intercept-only N dependence: n_ped in {20,50,100,500,1000}, core=n_g,
    C++ exact vs apy-approx SE/Z/LOG10P slopes
  * fixed-effect designs (Python-only, EE fixed from C++ exact): intercept /
    +binary / +continuous / +categorical(4)
  * tau/omega cases (C++ exact + apy-approx core=n_g): (identity,w=0),
    (a22,w=0), (a22,w=1=production), (identity,w=1); tau=1 always; production
    values never modified.

Output: d4_summary.json
Usage:
  python3 tests/run_ssgwas_pvalue_apy_approx_decomposition.py \
      --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_apy_approx_decomposition --seed 123
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

# ---------------------------------------------------------------------------
# shared utilities
# ---------------------------------------------------------------------------


def meuwissen_luo(pedigree):
    n = len(pedigree)
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    A = np.zeros((n, n))
    for i, r in enumerate(pedigree):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        if s >= 0 and d >= 0:
            A[i, i] = 1.0 + 0.5 * A[s, d]
        elif s >= 0 or d >= 0:
            A[i, i] = 1.0
        else:
            A[i, i] = 1.0
        for j in range(i):
            A[i, j] = 0.5 * ((A[s, j] if s >= 0 else 0.0) +
                             (A[d, j] if d >= 0 else 0.0))
            A[j, i] = A[i, j]
    return A


def sparse_ainv_ml(pedigree):
    n = len(pedigree)
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    A = meuwissen_luo(pedigree)
    Ainv = np.zeros((n, n))
    for i, r in enumerate(pedigree):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        if s >= 0 and d >= 0:
            delta = 1.0 - 0.25 * A[s, s] - 0.25 * A[d, d]
        elif s >= 0 or d >= 0:
            delta = 1.0 - 0.25 * A[max(s, d), max(s, d)]
        else:
            delta = 1.0
        if delta <= 0.0:
            raise RuntimeError("Ainv delta <= 0")
        Ainv[i, i] += 1.0 / delta
        for p in (s, d):
            if p < 0:
                continue
            Ainv[i, p] -= 0.5 / delta
            Ainv[p, i] -= 0.5 / delta
            Ainv[p, p] += 0.25 / delta
        if s >= 0 and d >= 0:
            Ainv[s, d] += 0.25 / delta
            Ainv[d, s] += 0.25 / delta
    return Ainv


def write_plink_trio(bed_prefix, samples, m, p, rng):
    with open(bed_prefix + ".fam", "w") as f:
        for iid in samples:
            f.write(f"F{iid} {iid} 0 0 0 -9\n")
    with open(bed_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"1 rs{j + 1} 0 {j + 1} A G\n")
    n = len(samples)
    dos = np.array([rng.binomial(2, p[j], size=n) for j in range(m)]).T.astype(float)
    bit = {0: 0b11, 1: 0b10, 2: 0b00}
    bytes_per = (n + 3) // 4
    raw = bytearray()
    for j in range(m):
        for b in range(bytes_per):
            byte = 0
            for k in range(4):
                i = b * 4 + k
                if i < n:
                    byte |= bit[int(dos[i, j])] << (2 * k)
            raw.append(byte)
    with open(bed_prefix + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01" + bytes(raw))
    return dos


def run_cs(bin_, work, out, extra, seed=7):
    os.makedirs(out, exist_ok=True)
    cmd = [bin_, "--bfile", f"{work}/sim", "--pheno", f"{work}/sim_y.phe",
           "--gwas", "--single-step", "--single-step-method", "backsolve",
           "--ped", f"{work}/sim.cs.ped", "--reml-vg", "1.0000000000",
           "--reml-ve", "3.0000000000", "--no-loco", "--ssgwas-window-snps",
           "10", "--ssgwas-apy-seed", str(seed), "--out", f"{out}/m2"] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"run failed: {cmd}\n{r.stderr[-2000:]}")
    return r


def load_pvalue(path):
    rows = sorted([l.split("\t") for l in
                   open(path).read().splitlines()[1:]],
                  key=lambda r: int(r[1][2:]))
    rows = [r for r in rows if all(v != "NA" for v in
                                   (r[6], r[7], r[8], r[9], r[11]))]
    return (np.array([float(r[6]) for r in rows]),   # EA
            np.array([float(r[7]) for r in rows]),   # SE
            np.array([float(r[8]) for r in rows]),   # VAR_EFFECT
            np.array([float(r[9]) for r in rows]),   # Z
            np.array([float(r[11]) for r in rows]))  # LOG10P


def load_meta(path):
    return {l.split("\t")[0]: l.split("\t")[1] for l in
            open(path).read().splitlines() if "\t" in l}


def pearson(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return 1.0 if np.allclose(a, b) else float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def metrics(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    slope = 1.0
    if np.std(a) > 0.0:
        slope = float(np.cov(a, b, ddof=0)[0, 1] / np.var(a))
    return {"pearson": pearson(a, b), "slope": slope,
            "rmse": float(np.sqrt(np.mean((a - b) ** 2))),
            "max_abs": float(np.max(np.abs(a - b)))}


def log10p_from_z(z):
    p = np.array([float(math.erfc(abs(v) / np.sqrt(2.0))) for v in z])
    p = np.maximum(p, 0.0)
    return np.array([-np.log10(v) if v > 0.0 else 320.0 for v in p])


# ---------------------------------------------------------------------------
# mt19937_64 replica + seeded APY core (must match C++ exactly)
# ---------------------------------------------------------------------------

MT64_N, MT64_M = 312, 156
MT64_A = 0xB5026F5AA96619E9
MT64_UM = 0xFFFFFFFF80000000
MT64_LM = 0x7FFFFFFF
MT64_F = 6364136223846793005
MT64_B = 0x71D67FFFEDA60000
MT64_C = 0xFFF7EEE000000000
MT64_D = 0x5555555555555555


class Mt19937_64:
    def __init__(self, seed):
        self.mt = [seed & 0xFFFFFFFFFFFFFFFF]
        for i in range(1, MT64_N):
            x = self.mt[i - 1] ^ (self.mt[i - 1] >> 62)
            self.mt.append((MT64_F * x + i) & 0xFFFFFFFFFFFFFFFF)
        self.mti = MT64_N

    def _twist(self):
        for i in range(MT64_N):
            x = (self.mt[i] & MT64_UM) | (self.mt[(i + 1) % MT64_N] & MT64_LM)
            xa = x >> 1
            if x & 1:
                xa ^= MT64_A
            self.mt[i] = self.mt[(i + MT64_M) % MT64_N] ^ xa
        self.mti = 0

    def next(self):
        if self.mti >= MT64_N:
            self._twist()
        x = self.mt[self.mti]
        self.mti += 1
        x ^= (x >> 29) & MT64_D
        x ^= (x << 17) & MT64_B
        x ^= (x << 37) & MT64_C
        x ^= x >> 43
        return x & 0xFFFFFFFFFFFFFFFF


def seeded_core_indices(ng, n_core, seed):
    rng = Mt19937_64(seed)
    perm = list(range(ng))
    for i in range(ng):
        j = i + (rng.next() % (ng - i))
        perm[i], perm[j] = perm[j], perm[i]
    return sorted(perm[:n_core])


# ---------------------------------------------------------------------------
# matrix machinery
# ---------------------------------------------------------------------------


def build_g_bs(dos, alpha=0.05):
    M = dos - dos.mean(axis=0)
    p = dos.mean(axis=0) / 2.0
    S = float(np.sum(2.0 * p * (1.0 - p)))
    G = M @ M.T / S
    return (1.0 - alpha) * G + alpha * np.eye(G.shape[0])


def tune_a22(G, A22):
    ng = G.shape[0]
    dG, oG = np.trace(G) / ng, (G.sum() - np.trace(G)) / (ng * (ng - 1))
    dA, oA = np.trace(A22) / ng, (A22.sum() - np.trace(A22)) / (ng * (ng - 1))
    b = (dA - oA) / (dG - oG)
    a = dA - b * dG
    return b * G + a, b, a


def apy_inverse(G, core):
    ng = G.shape[0]
    nset = set(core)
    noncore = [i for i in range(ng) if i not in nset]
    Gcc = G[np.ix_(core, core)]
    Gcn = G[np.ix_(core, noncore)]
    Gnn = G[np.ix_(noncore, noncore)]
    Gcc_inv = np.linalg.inv(Gcc)
    Q = Gcc_inv @ Gcn
    Mnn_inv = np.diag(1.0 / np.diag(Gnn - Gcn.T @ Q))
    out = np.zeros((ng, ng))
    c = np.array(core, dtype=int)
    o = np.array(noncore, dtype=int)
    out[np.ix_(c, c)] = Gcc_inv + Q @ Mnn_inv @ Q.T
    out[np.ix_(c, o)] = -Q @ Mnn_inv
    out[np.ix_(o, c)] = -Mnn_inv @ Q.T
    out[np.ix_(o, o)] = Mnn_inv
    return out


def build_c_lambda(X, W, Ainv, gblock, geno_idx_ped, lam, h_omega, A22_inv):
    p = X.shape[1]
    n_ped = Ainv.shape[0]
    XtX = X.T @ X
    XtW = X.T @ W
    big = np.zeros((p + n_ped, p + n_ped))
    big[:p, :p] = XtX
    big[:p, p:] = XtW
    big[p:, :p] = XtW.T
    embed = np.zeros((n_ped, n_ped))
    embed[np.ix_(geno_idx_ped, geno_idx_ped)] = gblock - h_omega * A22_inv
    big[p:, p:] = W.T @ W + lam * (Ainv + embed)
    return big


def scan_se(vhat, op_G, core, M, S2):
    """per-SNP SE with core operator w_c = (op_G_cc)^-1 M_cj (Bermann Eq. 21)."""
    Gcc_inv = np.linalg.inv(op_G[np.ix_(core, core)])
    m = M.shape[1]
    out = np.zeros(m)
    for j in range(m):
        w_c = Gcc_inv @ M[np.ix_(core, [j])][:, 0]
        v = w_c @ (vhat @ w_c)
        out[j] = np.sqrt(max(v, 0.0)) / S2
    return out


def d_full_from_X(X, W):
    n_rec = X.shape[0]
    Px = np.eye(n_rec) - X @ np.linalg.inv(X.T @ X) @ X.T
    return W.T @ Px @ W


# ---------------------------------------------------------------------------
# data + runs
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    work = f"{args.outdir}/data"
    os.makedirs(work, exist_ok=True)

    # --- tiny dataset (identical recipe to D3: n_ped=30, n_g=20, m=30) ---
    n_ped, n_geno, m = 30, 20, 30
    n_f = n_ped // 3
    pedigree = []
    for i in range(1, n_f + 1):
        pedigree.append([i, 0, 0])
    for i in range(n_f + 1, n_ped + 1):
        s = int(rng.integers(1, n_f + 1))
        d = int(rng.integers(1, n_f + 1))
        while d == s:
            d = int(rng.integers(1, n_f + 1))
        pedigree.append([i, s, d])
    geno_ids = list(range(n_f + 1, n_ped + 1))
    p_maf = rng.uniform(0.15, 0.45, size=m)
    dos = write_plink_trio(f"{work}/sim", geno_ids, m, p_maf, rng)
    pheno_ids = sorted(set(rng.choice(list(range(1, n_ped + 1)), size=21,
                                      replace=False)))
    with open(f"{work}/sim.phe", "w") as f:
        f.write("FID IID Pheno_1\n")
        for i in pheno_ids:
            f.write(f"F{i} {i} 0.0\n")
    with open(f"{work}/sim.cs.ped", "w") as f:
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")

    idx = {r[0]: i for i, r in enumerate(pedigree)}
    geno_idx_ped = np.array([idx[i] for i in geno_ids])
    A = meuwissen_luo(pedigree)
    A22 = A[np.ix_(geno_idx_ped, geno_idx_ped)]
    A22_inv = np.linalg.inv(A22)
    Ainv_ml = sparse_ainv_ml(pedigree)
    ml_gate = float(np.max(np.abs(Ainv_ml - np.linalg.inv(A))))
    if ml_gate > 1e-8:
        raise RuntimeError(f"Ainv ML gate failed: {ml_gate}")

    M = dos - dos.mean(axis=0)
    pj = dos.mean(axis=0) / 2.0
    S_2pq = float(np.sum(2.0 * pj * (1.0 - pj)))

    G_raw = M @ M.T / S_2pq
    G_star_draw, _, _ = tune_a22(G_raw.copy(), A22)
    Gw_draw = (1.0 - 0.05) * G_star_draw + 0.05 * A22
    H = A.copy()
    H[np.ix_(geno_idx_ped, geno_idx_ped)] = Gw_draw
    w_eig, V_eig = np.linalg.eigh(H)
    w_eig = np.maximum(w_eig, 0.0)
    pheno_rows = np.array([idx[i] for i in pheno_ids])
    Va, Ve = 1.0, 3.0
    a_true = (V_eig @ (np.sqrt(w_eig) * rng.normal(size=n_ped)))
    y = 1.0 + a_true[pheno_rows] + rng.normal(scale=np.sqrt(Ve),
                                              size=len(pheno_rows))
    with open(f"{work}/sim_y.phe", "w") as f:
        f.write("FID IID Pheno_1\n")
        for k, r in enumerate(pheno_rows):
            f.write(f"F{pedigree[r][0]} {pedigree[r][0]} {y[k]:.10f}\n")

    G_raw = M @ M.T / S_2pq
    G_bs_id = (1.0 - 0.05) * G_raw + 0.05 * np.eye(n_geno)
    G_star_a22, tune_b, tune_a = tune_a22(G_raw.copy(), A22)
    G_bs_a22 = (1.0 - 0.05) * G_star_a22 + 0.05 * A22

    Va, Ve = 1.0, 3.0
    lam = Ve / Va

    pheno_rows = np.array([idx[i] for i in pheno_ids])
    n_rec = len(pheno_rows)
    W = np.zeros((n_rec, n_ped))
    for k, r in enumerate(pheno_rows):
        W[k, r] = 1.0
    X0 = np.ones((n_rec, 1))

    # --- C++ runs: modes x (exact, apy-exact, apy-approx) x cores ---
    runs = {}
    for mode, flags, key in [("identity", [], "identity"),
                             ("a22", ["--ssgwas-tune-g",
                                      "--ssgwas-blend-a22"], "a22")]:
        for method, cores in [("exact", [None]),
                              ("apy-exact", [n_geno, n_geno - 5]),
                              ("apy-approx", [n_geno, n_geno - 5])]:
            for core in cores:
                tag = f"{key}_{method}" + ("" if core is None
                                           else f"_core{core}")
                extra = ["--ssgwas-pvalue-method", method] + (
                    [] if core is None else
                    ["--ssgwas-apy-core-size", str(core)])
                run_cs(args.bin, work, f"{args.outdir}/{tag}", flags + extra)
                runs[tag + "_p"] = load_pvalue(
                    f"{args.outdir}/{tag}/m2.ss_backsolve.pvalue.txt")
                if method == "apy-approx":
                    runs[tag + "_m"] = load_meta(
                        f"{args.outdir}/{tag}/m2.ss_backsolve.meta.txt")

    for tag, flags in [("case_identity_omega0", ["--h2-omega", "0"]),
                       ("case_a22_omega0", ["--h2-omega", "0",
                                            "--ssgwas-tune-g",
                                            "--ssgwas-blend-a22"])]:
        for method, core in [("exact", None), ("apy-approx", n_geno)]:
            t = f"{tag}_{method}"
            extra = ["--ssgwas-pvalue-method", method] + (
                [] if core is None else ["--ssgwas-apy-core-size", str(core)])
            run_cs(args.bin, work, f"{args.outdir}/{t}", flags + extra)
            runs[t + "_p"] = load_pvalue(
                f"{args.outdir}/{t}/m2.ss_backsolve.pvalue.txt")

    ea = runs["identity_exact_p"][0]

    # --- Python layer SE models ---
    se = {}
    for mode_key, G_bs, op_G in [("identity", G_bs_id, G_bs_id),
                                 ("a22", G_bs_a22, G_star_a22)]:
        Ginv_dense = np.linalg.inv(G_bs)
        for n_core in [n_geno, n_geno - 5]:
            core = seeded_core_indices(n_geno, n_core, 7)
            nset = set(core)
            noncore = [i for i in range(n_geno) if i not in nset]
            Ginv_apy = apy_inverse(G_bs, core)
            Gcc = G_bs[np.ix_(core, core)]
            # R0: dense-Ginv MME (only meaningful at core = n_g)
            if n_core == n_geno:
                C = build_c_lambda(X0, W, Ainv_ml, Ginv_dense, geno_idx_ped,
                                   lam, 1.0, A22_inv)
                Cinv_gg = np.linalg.inv(C)[1:, 1:][
                    np.ix_(geno_idx_ped, geno_idx_ped)]
                vhat_r0 = Va * G_bs - Ve * Cinv_gg
                se[f"R0_{mode_key}_core{n_core}"] = scan_se(
                    vhat_r0, op_G, core, M, S_2pq)
                se[f"R0full_{mode_key}"] = scan_se(
                    vhat_r0, op_G, list(range(n_geno)), M, S_2pq)
            # R1: APY-Ginv MME, exact core slice of the actual MME inverse
            C = build_c_lambda(X0, W, Ainv_ml, Ginv_apy, geno_idx_ped,
                               lam, 1.0, A22_inv)
            Cinv = np.linalg.inv(C)
            Cinv_cc = Cinv[1:, 1:][np.ix_(geno_idx_ped, geno_idx_ped)][
                np.ix_(core, core)]
            vhat_r1 = Va * Gcc - Ve * Cinv_cc
            se[f"R1_{mode_key}_core{n_core}"] = scan_se(
                vhat_r1, op_G, core, M, S_2pq)
            # R2: Eq. 8 with FULL D; R3: Eq. 8 with diag(D) (production)
            D_full_gg = d_full_from_X(X0, W)[np.ix_(geno_idx_ped,
                                                    geno_idx_ped)]
            D_c = D_full_gg[np.ix_(core, core)]
            D_o = D_full_gg[np.ix_(noncore, noncore)]
            D_co = D_full_gg[np.ix_(core, noncore)]
            G11, G12, G22 = (Ginv_apy[np.ix_(core, core)],
                             Ginv_apy[np.ix_(core, noncore)],
                             Ginv_apy[np.ix_(noncore, noncore)])
            M11 = D_c + lam * G11
            M12 = D_co + lam * G12
            M22 = D_o + lam * G22
            P2 = np.linalg.inv(M11 - M12 @ np.linalg.inv(M22) @ M12.T)
            M11d = np.diag(np.diag(D_c)) + lam * G11
            M12d = lam * G12
            M22d = np.diag(np.diag(D_o)) + lam * G22
            P3 = np.linalg.inv(M11d - M12d @ np.linalg.inv(M22d) @ M12d.T)
            v2 = Va * Gcc - Ve * P2
            v3 = Va * Gcc - Ve * P3
            se[f"R2_{mode_key}_core{n_core}"] = scan_se(v2, op_G, core, M,
                                                        S_2pq)
            se[f"R3_{mode_key}_core{n_core}"] = scan_se(v3, op_G, core, M,
                                                        S_2pq)

    # operator_full_vs_core: same exact PEV, full vs core operator
    op_blocks = {}
    for mode_key, G_bs, op_G in [("identity", G_bs_id, G_bs_id)]:
        C = build_c_lambda(X0, W, Ainv_ml, np.linalg.inv(G_bs), geno_idx_ped,
                           lam, 1.0, A22_inv)
        Cinv_gg = np.linalg.inv(C)[1:, 1:][np.ix_(geno_idx_ped,
                                                  geno_idx_ped)]
        vhat_full = Va * G_bs - Ve * Cinv_gg
        Ginv_dense = np.linalg.inv(op_G)
        full_se = []
        for j in range(m):
            w_f = Ginv_dense @ M[:, j]
            full_se.append(np.sqrt(max(w_f @ (vhat_full @ w_f), 0.0)) / S_2pq)
        for n_core in [n_geno, n_geno - 5]:
            core = seeded_core_indices(n_geno, n_core, 7)
            Gcc_inv = np.linalg.inv(op_G[np.ix_(core, core)])
            vhat_cc = vhat_full[np.ix_(core, core)]
            core_se = []
            for j in range(m):
                w_c = Gcc_inv @ M[core, j]
                core_se.append(np.sqrt(max(w_c @ (vhat_cc @ w_c), 0.0)) / S_2pq)
            op_blocks[f"identity_core{n_core}"] = (
                np.array(full_se), np.array(core_se))

    # --- gates ---
    fails = []
    gtol = 1e-6
    for mode_key, key in [("identity", "identity"), ("a22", "a22")]:
        cpp_ex = runs[f"{key}_exact_p"][1]
        if np.max(np.abs(se[f"R0_{mode_key}_core{n_geno}"] - cpp_ex)) > gtol:
            fails.append(f"{key} exact SE Python-vs-C++")
        for n_core in [n_geno, n_geno - 5]:
            ae = f"{key}_apy-exact_core{n_core}"
            aa = f"{key}_apy-approx_core{n_core}"
            if np.max(np.abs(se[f"R1_{mode_key}_core{n_core}"] -
                             runs[ae + "_p"][1])) > gtol:
                fails.append(f"{key} apy-exact core{n_core} SE match")
            if np.max(np.abs(se[f"R3_{mode_key}_core{n_core}"] -
                             runs[aa + "_p"][1])) > gtol:
                fails.append(f"{key} apy-approx core{n_core} SE match")

    # B==C gate: production Eq.9 diagonal D == diag(D_full) per animal
    # (genotyped with a record: 1 - x'(X'X)^-1 x; without: 0)
    D_full_run = d_full_from_X(X0, W)
    D_full_gg = D_full_run[np.ix_(geno_idx_ped, geno_idx_ped)]
    has_rec = np.isin(geno_idx_ped, pheno_rows)
    prod_d = np.where(has_rec, 1.0 - 1.0 / n_rec, 0.0)
    bc = float(np.max(np.abs(prod_d - np.diag(D_full_gg))))
    if bc > 1e-10:
        fails.append(f"B==C gate failed: {bc}")

    # --- JSON blocks ---
    d4 = {}

    def zlp(se_arr):
        z = ea / np.asarray(se_arr)
        return z, log10p_from_z(z)

    sup = " (R1 == R2-reference by construction: the exact core slice "
    sup += "of the actual MME inverse IS the core Schur complement)"
    for mode_key, key in [("identity", "identity"), ("a22", "a22")]:
        for n_core in [n_geno, n_geno - 5]:
            tag = f"{mode_key}_core{n_core}"
            ex = runs[f"{key}_exact_p"]
            ae = runs[f"{key}_apy-exact_core{n_core}_p"]
            aa = runs[f"{key}_apy-approx_core{n_core}_p"]
            z0, lp0 = zlp(ex[1])
            z1, lp1 = zlp(ae[1])
            z3, lp3 = zlp(aa[1])
            # R0 vs R1 (APY-G + core representation)
            d4[f"exact_vs_apy_exact_{tag}"] = {
                "EFFECT": metrics(ex[0], ae[0]),
                "SE": metrics(ex[1], ae[1]),
                "Z": metrics(z0, z1),
                "LOG10P": metrics(lp0, lp1)}
            # R1(core=c) vs R1(core=n_g)=R0 (core representation with exact PEV)
            d4[f"apy_exact_vs_fullcore_{tag}"] = {
                "EFFECT": metrics(ex[0], ae[0]),
                "SE": metrics(ex[1], ae[1]),
                "Z": metrics(z0, z1),
                "LOG10P": metrics(lp0, lp1)}
            # R1 vs R3 (whole Eq8/Eq9 approximation, same core)
            d4[f"fullcore_vs_eq9_{tag}"] = {
                "EFFECT": metrics(ae[0], aa[0]),
                "SE": metrics(ae[1], aa[1]),
                "Z": metrics(z1, z3),
                "LOG10P": metrics(lp1, lp3)}
            # R2 vs R3 (full-D vs diag-D, Python)
            z2, lp2 = zlp(se[f"R2_{mode_key}_core{n_core}"])
            z3b, lp3b = zlp(se[f"R3_{mode_key}_core{n_core}"])
            d4[f"fullD_vs_diagD_{tag}"] = {
                "SE": metrics(se[f"R2_{mode_key}_core{n_core}"],
                              se[f"R3_{mode_key}_core{n_core}"]),
                "Z": metrics(z2, z3b),
                "LOG10P": metrics(lp2, lp3b)}
            # R1 vs R2 (structural: actual MME vs D-model, Python)
            z1b, lp1b = zlp(se[f"R1_{mode_key}_core{n_core}"])
            d4[f"MME_vs_Dmodel_{tag}"] = {
                "SE": metrics(se[f"R1_{mode_key}_core{n_core}"],
                              se[f"R2_{mode_key}_core{n_core}"]),
                "Z": metrics(z1b, z2),
                "LOG10P": metrics(lp1b, lp2)}
    for n_core in [n_geno, n_geno - 5]:
        fs, cs = op_blocks[f"identity_core{n_core}"]
        z_f, lp_f = zlp(fs)
        z_c, lp_c = zlp(cs)
        d4[f"operator_full_vs_core_identity_core{n_core}"] = {
            "SE": metrics(fs, cs), "Z": metrics(z_f, z_c),
            "LOG10P": metrics(lp_f, lp_c)}
    d4["block_notes"] = {
        "R1_is_R2reference": sup,
        "operator_block": "same exact PEV; full operator w vs core operator "
                          "w_c = (op_G_cc)^-1 M_cj (Bermann Eq. 21)"}

    # --- N dependence ---
    n_dep = []
    for npd in [20, 50, 100, 500, 1000]:
        ng_ = 2 * npd // 3
        nf = npd // 3
        r2 = np.random.default_rng(args.seed * 1000 + npd)
        ped = [[i, 0, 0] for i in range(1, nf + 1)]
        for i in range(nf + 1, npd + 1):
            s = int(r2.integers(1, nf + 1))
            d = int(r2.integers(1, nf + 1))
            while d == s:
                d = int(r2.integers(1, nf + 1))
            ped.append([i, s, d])
        gids = list(range(nf + 1, npd + 1))
        pm = r2.uniform(0.15, 0.45, size=60)
        d2 = write_plink_trio(f"{work}/sim_n{npd}", gids, 60, pm, r2)
        pids = sorted(set(r2.choice(list(range(1, npd + 1)),
                                    size=int(0.7 * npd), replace=False)))
        pidx = {r[0]: i for i, r in enumerate(ped)}
        gpidx = np.array([pidx[i] for i in gids])
        A_n = meuwissen_luo(ped)
        A22_n = A_n[np.ix_(gpidx, gpidx)]
        M_n = d2 - d2.mean(axis=0)
        pj_n = d2.mean(axis=0) / 2.0
        S_n = float(np.sum(2.0 * pj_n * (1.0 - pj_n)))
        G_n = M_n @ M_n.T / S_n
        Gs_n, _, _ = tune_a22(G_n.copy(), A22_n)
        Gw_n = (1.0 - 0.05) * Gs_n + 0.05 * A22_n
        H_n = A_n.copy()
        H_n[np.ix_(gpidx, gpidx)] = Gw_n
        w_e, V_e = np.linalg.eigh(H_n)
        w_e = np.maximum(w_e, 0.0)
        p_rows = np.array([pidx[i] for i in pids])
        a_n = V_e @ (np.sqrt(w_e) * r2.normal(size=npd))
        y_n = 1.0 + a_n[p_rows] + r2.normal(scale=np.sqrt(Ve),
                                            size=len(p_rows))
        with open(f"{work}/sim_n{npd}.phe", "w") as f:
            f.write("FID IID Pheno_1\n")
            for k, r in enumerate(p_rows):
                f.write(f"F{ped[r][0]} {ped[r][0]} {y_n[k]:.10f}\n")
        with open(f"{work}/sim_n{npd}.cs.ped", "w") as f:
            for r in ped:
                f.write(f"{r[0]} {r[1]} {r[2]}\n")
        pvs = {}
        for meth, extra in [("exact", ["--ssgwas-pvalue-method", "exact"]),
                            ("approx", ["--ssgwas-pvalue-method",
                                        "apy-approx",
                                        "--ssgwas-apy-core-size", str(ng_)])]:
            out = f"{args.outdir}/n{npd}_{meth}"
            os.makedirs(out, exist_ok=True)
            cmd = [args.bin, "--bfile", f"{work}/sim_n{npd}", "--pheno",
                   f"{work}/sim_n{npd}.phe", "--gwas", "--single-step",
                   "--single-step-method", "backsolve", "--ped",
                   f"{work}/sim_n{npd}.cs.ped", "--reml-vg", "1.0",
                   "--reml-ve", "3.0", "--no-loco", "--ssgwas-window-snps",
                   "10", "--ssgwas-apy-seed", "7"] + extra + [
                   "--out", f"{out}/m2"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"N-run {npd}/{meth} failed: "
                                   f"{r.stderr[-1000:]}")
            pvs[meth] = load_pvalue(f"{out}/m2.ss_backsolve.pvalue.txt")
        n_dep.append({"n_ped": npd, "n_g": ng_,
                      "SE": metrics(pvs["exact"][1], pvs["approx"][1]),
                      "Z": metrics(pvs["exact"][3], pvs["approx"][3]),
                      "LOG10P": metrics(pvs["exact"][4], pvs["approx"][4])})
    d4["N_dependence"] = n_dep

    # --- fixed-effect designs (Python-only; EA fixed from C++ exact) ---
    fec = []
    designs = [("intercept", X0),
               ("+binary", np.hstack([X0, rng.integers(0, 2,
                (n_rec, 1)).astype(float)])),
               ("+continuous", np.hstack([X0, rng.normal(size=(n_rec, 1))])),
               ("+categorical4", np.hstack(
                   [X0] + [(np.tile(np.arange(4), n_rec // 4 + 1)[:n_rec]
                            == (k + 1)).astype(float)[:, None]
                           for k in range(3)])) ]
    for des, Xd in designs:
        Df_gg = d_full_from_X(Xd, W)[np.ix_(geno_idx_ped, geno_idx_ped)]
        C = build_c_lambda(Xd, W, Ainv_ml, np.linalg.inv(G_bs_id),
                           geno_idx_ped, lam, 1.0, A22_inv)
        Cinv_gg = np.linalg.inv(C)[Xd.shape[1]:, Xd.shape[1]:][
            np.ix_(geno_idx_ped, geno_idx_ped)]
        vhat_ex = Va * G_bs_id - Ve * Cinv_gg
        Ginv_apy = apy_inverse(G_bs_id, list(range(n_geno)))
        P2b = np.linalg.inv(Df_gg + lam * Ginv_apy)
        P3b = np.linalg.inv(np.diag(np.diag(Df_gg)) + lam * Ginv_apy)
        v2b = Va * G_bs_id - Ve * P2b
        v3b = Va * G_bs_id - Ve * P3b
        Ginvd = np.linalg.inv(G_bs_id)
        def se_of(vhat):
            out = np.zeros(m)
            for j in range(m):
                w = Ginvd @ M[:, j]
                out[j] = np.sqrt(max(w @ (vhat @ w), 0.0)) / S_2pq
            return out
        se_ex_f, se_fl_f, se_di_f = map(se_of, (vhat_ex, v2b, v3b))
        zE, zD = zlp(se_ex_f)[0], zlp(se_di_f)[0]
        zFl = zlp(se_fl_f)[0]
        fec.append({"design": des,
                    "SE": metrics(se_ex_f, se_di_f),
                    "Z": metrics(zE, zD),
                    "LOG10P": metrics(log10p_from_z(zE),
                                      log10p_from_z(zD)),
                    "full_vs_diag_SE": metrics(se_fl_f, se_di_f),
                    "full_vs_diag_Z": metrics(zFl, zD),
                    "D_offdiag_max_abs": float(np.max(np.abs(
                        Df_gg - np.diag(np.diag(Df_gg))))),
                    "D_offdiag_diag_ratio": float(np.max(np.abs(
                        Df_gg - np.diag(np.diag(Df_gg)))) /
                        np.max(np.diag(Df_gg)))})
    d4["fixed_effect_designs"] = fec

    # --- tau/omega four cases ---
    tw = []
    for tag, label in [("case_identity_omega0", "identity, omega=0"),
                       ("case_a22_omega0", "a22, omega=0"),
                       ("identity", "identity, omega=1 (production)"),
                       ("a22", "a22, omega=1 (production)")]:
        pe = runs[f"{tag}_exact_p"]
        pa = runs[f"{tag}_apy-approx_p"] if tag.startswith("case") else runs[
            f"{tag}_apy-approx_core{n_geno}_p"]
        zE_, lpE_ = zlp(pe[1])
        zA_, lpA_ = zlp(pa[1])
        tw.append({"case": label,
                   "SE": metrics(pe[1], pa[1]),
                   "Z": metrics(zE_, zA_),
                   "LOG10P": metrics(lpE_, lpA_)})
    d4["tau_omega_cases"] = tw

    d4["gates"] = {"Ainv_ML_maxdiff": ml_gate,
                   "B_eq_C_diag_maxdiff": bc,
                   "python_vs_cpp_SE_maxdiff_ok": len(fails) == 0}
    d4["D_full_offdiag"] = {
        "intercept_only_prediction": -1.0 / n_rec,
        "observed_max": float(np.max(np.abs(
            D_full_gg - np.diag(np.diag(D_full_gg)))))}
    d4["gc_applied_any_mode"] = any(
        runs.get(t + "_m", {}).get("GENOMIC_CONTROL_APPLIED") == "true"
        for t in ["identity_apy-approx_core20",
                  "identity_apy-approx_core" + str(n_geno)])

    with open(f"{args.outdir}/d4_summary.json", "w") as f:
        json.dump({"fails": fails, "d4": d4}, f, indent=2, default=str)

    # --- console report ---
    k = "identity_core20"
    b0 = d4[f"exact_vs_apy_exact_{k}"]
    b1 = d4[f"MME_vs_Dmodel_{k}"]
    b2 = d4[f"fullcore_vs_eq9_{k}"]
    b3 = d4[f"fullD_vs_diagD_{k}"]
    print("\n=== D4 identity core=n_g slopes (SE / Z / LOG10P) ===")
    def line(name, blk):
        print(f"{name:28s} SE {blk['SE']['slope']:.4f} | "
              f"Z {blk['Z']['slope']:.4f} | LOG10P {blk['LOG10P']['slope']:.4f}")
    line("R0->R1 (APY-G+core)", b0)
    line("R1->R2 (MME-vs-Dmodel)", b1)
    line("R2->R3 (fullD-vs-diagD)", b3)
    line("R1->R3 (Eq8/Eq9 total)", b2)
    b4 = d4["operator_full_vs_core_identity_core20"]
    print(f"{'operator full-vs-core':28s} SE {b4['SE']['slope']:.4f} | "
          f"Z {b4['Z']['slope']:.4f} | LOG10P {b4['LOG10P']['slope']:.4f}")
    print("B==C maxdiff:", d4["gates"]["B_eq_C_diag_maxdiff"],
          " offdiag:", d4["D_full_offdiag"])
    if fails:
        print("FAILURES:\n" + "\n".join(fails))
        sys.exit(1)
    print("PASS: D4 decomposition gates")


if __name__ == "__main__":
    main()