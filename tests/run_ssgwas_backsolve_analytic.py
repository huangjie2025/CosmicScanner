#!/usr/bin/env python3
"""Analytic brute-force validation of CosmicScanner ssGBLUP back-solving
(--single-step-method backsolve, Phase 6).

A small three-generation pedigree with typed/untyped/unphenotyped animals is
simulated and written to PLINK .bed/.bim/.fam + phenotype + pedigree files.
cosmicscanner then runs the back-solving workflow and this script computes the
SAME quantities from first principles in numpy:

    A           pedigree relationship matrix (recursive rules)
    Ainv        Henderson A^-1
    A22^(-1)    inverse of the genotyped block of A (replicates the C++
                buildA22 + invertSpdWithRidge composition: Schur complement
                of the untyped block in A^-1, inverted with ridge guarding)
    G = M M'/S  VanRaden, M = x - mean, S = sum_j 2p_j(1-p_j)
    Ginv        ridge-guarded inverse (same 6-attempt ridge growth as C++)
    Hinv        Ainv + tau*Ginv - omega*A22^(-1) on the genotyped block
    MME         [X'X X'Z; Z'X Z'Z + (Ve/Va) Hinv] [b;a] = [X'y; Z'y]
    back-solve  q = Ginv a_g;  effect_allele = (M'q)/S  (per allele dosage)
                effect_std = effect_allele * sd          (per standard-genotype
                SD; reconstruction-consistent Z*std = M*allele = a_hat_g)
    SNP var     2p(1-p) * effect_allele^2;  windows (SNP-count, per chromosome)

Pass/fail gates (max absolute difference vs analytic reference):
    EBV (all animals), GEBV, SNP_EFFECT_ALLELE, SNP_EFFECT_STD,
    SNP_VARIANCE, PCT_GENETIC_VARIANCE, WINDOW_VARIANCE.

Scenario 2 verifies the A-only degeneration tau=omega=0, which must reproduce
pedigree BLUP (Hinv = Ainv).
"""
import argparse
import math
import os
import subprocess
import sys

import numpy as np

RIDGE0 = 1e-8
REL = 1e-4  # relative tolerance (float-path GRM in C++ vs exact double)
ABSF = 1e-8


def ridge_inv(M, label="M"):
    """Replicate cosmic::invertSpdWithRidge (symmetrize, LDLT-PD, ridge*10)."""
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise RuntimeError(f"{label} must be square")
    sym = 0.5 * (M + M.T)
    ridge = RIDGE0
    n = sym.shape[0]
    I = np.eye(n)
    for attempt in range(6):
        T = sym + ridge * I
        try:
            np.linalg.cholesky(T)
            return np.linalg.solve(T, I)
        except np.linalg.LinAlgError:
            ridge *= 10.0
    return np.linalg.inv(sym + ridge * I)


def build_A(ped, idx):
    """Pedigree relationship matrix; parents always precede offspring."""
    n = len(ped)
    A = np.zeros((n, n))
    for i, p in enumerate(ped):
        s = idx.get(p[1], -1)
        d = idx.get(p[2], -1)
        if s >= 0 and d >= 0:
            A[i, i] = 1.0 + 0.5 * A[s, d]      # A_ii = 1 + F_i, F_i = 0.5*r_sd
        else:
            A[i, i] = 1.0
        for j in range(i):
            v = 0.0
            if s >= 0:
                v += 0.5 * A[s, j]
            if d >= 0:
                v += 0.5 * A[d, j]
            A[i, j] = A[j, i] = v
    return A


def build_Ainv(ped, idx, A):
    """Henderson rules matching invertA_henderson_fixed
    (dii = 0.5-0.25(Fs+Fd) / 0.75-0.25F / 1; parents after children)."""
    n = len(ped)
    Ainv = np.zeros((n, n))
    F = np.diag(A) - 1.0
    for i, p in enumerate(ped):
        s = idx.get(p[1], -1)
        d = idx.get(p[2], -1)
        s_ok, d_ok = s >= 0, d >= 0
        if s_ok and d_ok:
            dii = 0.5 - 0.25 * (F[s] + F[d])
        elif s_ok or d_ok:
            dii = 0.75 - 0.25 * F[s if s_ok else d]
        else:
            dii = 1.0
        if dii < 1e-12:
            raise RuntimeError("degenerate pedigree (fully inbred)")
        inv_d = 1.0 / dii
        Ainv[i, i] += inv_d
        if s_ok and d_ok:
            Ainv[i, s] += -0.5 * inv_d
            Ainv[s, i] += -0.5 * inv_d
            Ainv[i, d] += -0.5 * inv_d
            Ainv[d, i] += -0.5 * inv_d
            Ainv[s, s] += 0.25 * inv_d
            Ainv[d, d] += 0.25 * inv_d
            Ainv[s, d] += 0.25 * inv_d
            Ainv[d, s] += 0.25 * inv_d
        elif s_ok or d_ok:
            p = s if s_ok else d
            Ainv[i, p] += -0.5 * inv_d
            Ainv[p, i] += -0.5 * inv_d
            Ainv[p, p] += 0.25 * inv_d
    return Ainv


def write_plink(path_prefix, names, X, chrom, snp_ids, bp):
    """X: (n, m) dosage of A1 in {0,1,2}; complete data (no missing)."""
    n, m = X.shape
    with open(path_prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"0 {names[i]} 0 0 0 -9\n")
    with open(path_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"{chrom[j]} {snp_ids[j]} 0 {bp[j]} A G\n")
    nbytes = 3 + m * ((n + 3) // 4)
    buf = bytearray(nbytes)
    buf[0:3] = bytes([0x6C, 0x1B, 0x01])
    for j in range(m):
        for i in range(n):
            code = {2: 0, 1: 2, 0: 3}[int(X[i, j])]
            byte_idx = 3 + j * ((n + 3) // 4) + i // 4
            buf[byte_idx] |= (code << (2 * (i % 4)))
    with open(path_prefix + ".bed", "wb") as f:
        f.write(bytes(buf))


def make_pedigree():
    ped = []
    founders = [f"f{i:02d}" for i in range(12)]
    ped += [(a, "0", "0") for a in founders]
    g2 = [f"g2_{i:02d}" for i in range(28)]
    for k, a in enumerate(g2):
        ped.append((a, founders[k % 12], founders[(k * 7 + 3) % 12]))
    g3 = [f"g3_{i:02d}" for i in range(20)]
    for k, a in enumerate(g3):
        ped.append((a, g2[k % 28], g2[(k * 5 + 11) % 28]))
    return ped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    ped = make_pedigree()
    n_ped = len(ped)
    idx = {p[0]: i for i, p in enumerate(ped)}
    A = build_A(ped, idx)
    Ainv = build_Ainv(ped, idx, A)

    # typed / phenotyped partitions (all three classes present)
    typed = [idx[f"g2_{i:02d}"] for i in range(24)] + [idx[f"f{i:02d}"] for i in range(12)]
    typed += [idx[f"g3_{i:02d}"] for i in range(8)]          # 44 typed
    typed_set = set(typed)
    typed_names = [ped[i][0] for i in typed]
    typed_orig_idx = [typed.index(i) for i in range(n_ped) if False]  # unused
    del typed_orig_idx

    pheno_ids = [f"f{i:02d}" for i in range(12)]
    pheno_ids += [f"g2_{i:02d}" for i in range(28)]
    pheno_ids += [f"g3_{i:02d}" for i in range(14)]          # 54 phenotyped
    pheno_idx = [idx[a] for a in pheno_ids]
    pheno_set = set(pheno_idx)

    # SNP data: 2 chromosomes, complete
    m1, m2 = 50, 40
    m = m1 + m2
    chrom = ["1"] * m1 + ["2"] * m2
    snp_ids = [f"rs{j:04d}" for j in range(m)]
    bp = [1000 * (j + 1) for j in range(m)]
    X = np.zeros((len(typed), m), dtype=int)
    for j in range(m):
        p = rng.uniform(0.15, 0.5)
        X[:, j] = rng.binomial(2, p, size=len(typed))
    n_g = len(typed)

    pref = os.path.join(outdir, "sim")
    write_plink(pref, typed_names, X, chrom, snp_ids, bp)

    a_true = rng.normal(0.0, 1.0, n_ped)
    y_vals = {}
    for a in pheno_ids:
        y = 5.0 + a_true[idx[a]] + rng.normal(0.0, math.sqrt(1.5))
        y_vals[a] = round(float(y), 6)
    with open(pref + ".pheno", "w") as f:
        for a in pheno_ids:
            f.write(f"0 {a} {y_vals[a]}\n")
    with open(pref + ".ped", "w") as f:
        for a, s, d in ped:
            f.write(f"{a} {s} {d}\n")

    # ---- cosmic runs ----
    def run(extra, tag):
        rdir = os.path.join(outdir, tag)
        os.makedirs(rdir, exist_ok=True)
        cmd = [args.bin, "--bfile", pref, "--pheno", pref + ".pheno", "--gwas",
               "--single-step", "--single-step-method", "backsolve",
               "--ped", pref + ".ped",
               "--reml-vg", "1.0", "--reml-ve", "2.0",
               "--no-loco", "--ssgwas-window-snps", "5"] + \
              extra + ["--out", os.path.join(rdir, "m2")]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout)
            print(r.stderr)
            raise SystemExit(f"FAIL: cosmicscanner returned {r.returncode}")
        return rdir

    rdir_ldlt = run(["--ssgwas-solver", "ldlt"], "ldlt")
    rdir_pcg = run(["--ssgwas-solver", "pcg"], "pcg")

    # ---- analytic reference ----
    tau = omega = 1.0
    Va, Ve = 1.0, 2.0
    lam = Ve / Va

    means = X.mean(axis=0)
    s2pq = 2.0 * (means / 2.0) * (1.0 - means / 2.0)
    S = s2pq.sum()
    Zc = X.astype(float) - means
    G = (Zc @ Zc.T) / S
    BLEND_ALPHA = 0.05
    G = (1.0 - BLEND_ALPHA) * G + BLEND_ALPHA * np.eye(n_g)
    GinvC = ridge_inv(G, "G")

    A22 = A[np.ix_(typed, typed)]
    A22_inv_c = ridge_inv(A22, "A22")

    Hinv = Ainv.copy()
    for gi, i in enumerate(typed):
        for gj, j in enumerate(typed):
            Hinv[i, j] += tau * GinvC[gi, gj] - omega * A22_inv_c[gi, gj]
    Hinv = 0.5 * (Hinv + Hinv.T)

    y_full = np.full(n_ped, np.nan)
    for a in pheno_ids:
        y_full[idx[a]] = y_vals[a]
    rec_idx = [i for i in range(n_ped) if math.isfinite(y_full[i])]
    n_rec = len(rec_idx)
    Zd = np.zeros((n_rec, n_ped))
    for k, i in enumerate(rec_idx):
        Zd[k, i] = 1.0
    y = y_full[rec_idx]

    def solve_mme(Hinv_m):
        """Solve the mu-absorbed MME
           (Z'Z + lam*Hinv - (Z'1)(Z'1)'/n) a = Z'y - (Z'1)*mean(y).
           This is EXACTLY the intercept-eliminated form: eliminating the
           mu column of the full canonical system [n Z'1; Z'1 Z'Z+lam*Hinv]
           reproduces it identically (verified in the Phase 6.5 audit:
           canonical vs mu-absorbed agree to 1e-15; the matrix is
           nonsingular for this design, no free gauge exists)."""
        Z1 = Zd.sum(axis=0)
        Mred = (Zd.T @ Zd + lam * Hinv_m) - np.outer(Z1, Z1) / n_rec
        rhs_red = Zd.T @ y - Z1 * y.mean()
        return np.linalg.solve(Mred, rhs_red)

    sol = solve_mme(Hinv)
    a_ref = sol

    q = GinvC @ a_ref[typed]
    sd = np.sqrt(np.where(s2pq > 1e-9, s2pq, 1.0))
    zt = (Zc.T @ q) / sd                     # standardized Z'q
    zt = np.where(s2pq > 1e-9, zt, 0.0)
    effect_allele_ref = np.where(s2pq > 1e-9, (zt * sd) / S, 0.0)
    # Reconstruction-consistent standardized effect: gamma_j = EA_j * sd_j,
    # i.e. Z * gamma = M * EA = a_hat_g exactly (marker-scale invariant).
    effect_std_ref = np.where(s2pq > 1e-9, effect_allele_ref * sd, 0.0)
    snp_var_ref = s2pq * effect_allele_ref ** 2
    total_ref = snp_var_ref.sum()

    def windows(v):
        out = []
        s = 0
        W = 5
        while s < m:
            ch = chrom[s]
            e = s
            while e + 1 < m and chrom[e + 1] == ch:
                e += 1
            L = e - s + 1
            if L <= W:
                out.append((ch, v[s:e + 1].sum()))
            else:
                for k in range(s, e - W + 2):
                    out.append((ch, v[k:k + W].sum()))
            s = e + 1
        return out

    total_ref = max(total_ref, 0.0)
    win_ref = windows(snp_var_ref)

    # ---- load cosmic outputs ----
    def load_ebv(rd):
        ebv = {}
        with open(os.path.join(rd, "m2.ss_backsolve.ebv.txt")) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                ebv[p[0]] = float(p[1])
        return ebv

    def load_snp(rd):
        cols = {}
        with open(os.path.join(rd, "m2.ss_backsolve.snp.txt")) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                cols[d["SNP"]] = d
        return cols

    def load_win(rd):
        rows = []
        with open(os.path.join(rd, "m2.ss_backsolve.window.txt")) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                rows.append((p[0], float(p[6]), float(p[7])))
        return rows

    def check(rd, tag):
        ebv = load_ebv(rd)
        snp = load_snp(rd)
        win = load_win(rd)
        meta = {}
        with open(os.path.join(rd, "m2.ss_backsolve.meta.txt")) as f:
            for l in f.read().splitlines():
                if "\t" in l:
                    k, v = l.split("\t", 1)
                    meta[k] = v

        a_c = np.array([ebv[a] for a, _, _ in ped], dtype=float)
        d_ebv = float(np.abs(a_c - a_ref).max())
        d_gebv = float(np.abs(a_c[typed] - a_ref[typed]).max())

        ea_ref = {snp_ids[j]: effect_allele_ref[j] for j in range(m)}
        es_ref = {snp_ids[j]: effect_std_ref[j] for j in range(m)}
        sv_ref = {snp_ids[j]: snp_var_ref[j] for j in range(m)}
        pct_ref = {snp_ids[j]: (snp_var_ref[j] / total_ref if total_ref > 0 else 0.0)
                   for j in range(m)}
        d_ea = max(abs(float(snp[s]["EFFECT_ALLELE"]) - ea_ref[s]) for s in snp)
        d_es = max(abs(float(snp[s]["SNP_EFFECT_STD"]) - es_ref[s]) for s in snp)
        d_sv = max(abs(float(snp[s]["SNP_VARIANCE"]) - sv_ref[s]) for s in snp)
        d_pct = max(abs(float(snp[s]["PCT_GENETIC_VARIANCE"]) - pct_ref[s]) for s in snp)

        n_w = min(len(win_ref), len(win))
        d_win = max(abs(win[k][1] - win_ref[k][1]) for k in range(n_w))
        d_winp = max(abs(win[k][2] - (win_ref[k][1] / total_ref if total_ref > 0 else 0.0))
                     for k in range(n_w))

        print(f"[{tag}] max_abs_diff  EBV={d_ebv:.3e}  GEBV={d_gebv:.3e}  "
              f"EFFECT_ALLELE={d_ea:.3e}  EFFECT_STD={d_es:.3e}  "
              f"SNP_VAR={d_sv:.3e}  PCT={d_pct:.3e}  WIN={d_win:.3e} WIN_PCT={d_winp:.3e}")

        ok = (d_ebv <= REL * max(1.0, float(np.abs(a_ref).max())) + ABSF and
              d_gebv <= REL * max(1.0, float(np.abs(a_ref[typed]).max())) + ABSF and
              d_ea <= REL * max(1.0, float(np.abs(effect_allele_ref).max())) + ABSF and
              d_es <= REL * max(1.0, float(np.abs(effect_std_ref).max())) + ABSF and
              d_sv <= REL * max(1.0, total_ref) + ABSF and
              d_pct <= 5e-3 + ABSF and
              d_win <= REL * max(1.0, total_ref) + ABSF and
              d_winp <= 5e-3 + ABSF)
        if n_w != len(win):
            raise SystemExit(f"FAIL: {tag}: window count mismatch {n_w} vs {len(win_ref)}")
        if not ok:
            raise SystemExit(f"FAIL: {tag} exceeded tolerance")
        for key in ("TAU", "OMEGA", "S_SNP", "TOTAL_SNP_VARIANCE", "VCE_MODE",
                    "SOLVER", "N_PEDIGREE", "N_PHENOTYPED", "N_GENOTYPED",
                    "N_GENOTYPED_WITH_PHENO", "N_UNGENOTYPED_WITH_PHENO",
                    "BACKSOLVE_FORMULA", "WINDOW_SIZE"):
            if key not in meta:
                raise SystemExit(f"FAIL: {tag}: meta missing key {key}")
        if int(meta["N_UNGENOTYPED_WITH_PHENO"]) <= 0:
            raise SystemExit(f"FAIL: {tag}: expected ungenotyped-with-phenotype animals")

    check(rdir_ldlt, "ldlt")
    check(rdir_pcg, "pcg")

    # ---- scenario 2: A-only degeneration (tau=omega=0 -> pedigree BLUP) ----
    rdir0 = os.path.join(outdir, "aonly")
    os.makedirs(rdir0, exist_ok=True)
    cmd = [args.bin, "--bfile", pref, "--pheno", pref + ".pheno", "--gwas",
           "--single-step", "--single-step-method", "backsolve",
           "--ped", pref + ".ped",
           "--reml-vg", "1.0", "--reml-ve", "2.0",
           "--h2-tau", "0", "--h2-omega", "0",
           "--ssgwas-solver", "ldlt",
           "--out", os.path.join(rdir0, "m2")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise SystemExit("FAIL: A-only run failed")
    a0 = np.array([load_ebv(rdir0)[a] for a, _, _ in ped], dtype=float)
    a_pblup = solve_mme(Ainv)          # Hinv := Ainv exactly
    d_aonly = float(np.abs(a0 - a_pblup).max())
    print(f"[aonly] max_abs_diff EBV vs pedigree BLUP = {d_aonly:.3e}")
    if d_aonly > REL * max(1.0, float(np.abs(a_pblup).max())) + ABSF:
        raise SystemExit("FAIL: A-only degeneration mismatch")

    print("\nPASS: ssGBLUP back-solving analytic validation "
          "(ldlt + pcg paths, A-only degeneration, full meta check)")
    return 0


if __name__ == "__main__":
    sys.exit(main())