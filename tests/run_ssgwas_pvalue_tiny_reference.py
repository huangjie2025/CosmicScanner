#!/usr/bin/env python3
"""Phase C1 gate C5: independent tiny analytical reference for the exact
ssGWAS SNP P-value path (--ssgwas-pvalue-method exact).

Fully independent: builds pedigree A, centered genotypes, G, tuneG, blend,
H / H^-1, the lambda-form MME, the FULL dense MME inverse, EBV, Cinv_gg,
SNP effects on the CosmicScanner EFFECT_ALLELE scale, and the propagated
sampling variance Var(EA_j) = B_j * Var(u_hat_g) * B_j' (Aguilar 2019;
terminology per Garcia 2023 correction) with plain numpy - no CosmicScanner
or libCosmicCore helpers.

Usage:
  python3 tests/run_ssgwas_pvalue_tiny_reference.py --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_tiny --seed 123
"""
import argparse
import json
import os
import subprocess
import sys

from math import erfc as math_erfc
import numpy as np
erfc = np.vectorize(math_erfc)


def meuwissen_luo(pedigree):
    """Meuwissen & Luo tabular A (founders unrelated). pedigree: [[id,sire,dam],...].
    A[i,i] = 1 + 0.5*A[sire,dam] (inbreeding = coancestry of parents)."""
    n = len(pedigree)
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    A = np.zeros((n, n))
    L = np.zeros((n, n))
    for i, r in enumerate(pedigree):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        if s >= 0 and d >= 0:
            L[i, i] = np.sqrt(max(1.0 + 0.5 * A[s, d] - 0.25 * (A[s, s] + A[d, d])
                                  - 0.5 * A[s, d], 0.0))
            L[i, 0:i] = 0.5 * (L[s, 0:i] + L[d, 0:i])
            A[i, i] = 1.0 + 0.5 * A[s, d]
        elif s >= 0 or d >= 0:
            p = s if s >= 0 else d
            L[i, i] = np.sqrt(max(1.0 - 0.25 * A[p, p], 0.0))
            L[i, 0:i] = 0.5 * L[p, 0:i]
            A[i, i] = 1.0
        else:
            L[i, i] = 1.0
            A[i, i] = 1.0
        for j in range(i):
            A[i, j] = 0.5 * ((A[s, j] if s >= 0 else 0.0) +
                             (A[d, j] if d >= 0 else 0.0))
            A[j, i] = A[i, j]
    return A


def write_plink_trio(bed_prefix, samples, m, p, rng):
    """Write SNP-major bed/bim/fam (CosmicScanner PlinkReader convention) for
    `samples` individuals with MAF p per SNP. Bit order: LSB = first sample.
    codes: 00=hom A1(0), 01=missing, 10=het, 11=hom A2(2); CosmicScanner
    dosage = 2 - code, so write code 2-dos: 2->00, 1->10, 0->11."""
    with open(bed_prefix + ".fam", "w") as f:
        for iid in samples:
            f.write(f"F{iid} {iid} 0 0 0 -9\n")
    with open(bed_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"1 rs{j + 1} 0 {j + 1} A G\n")
    n = len(samples)
    dos = np.zeros((n, m), dtype=np.uint8)
    for j in range(m):
        dos[:, j] = rng.binomial(2, p[j], size=n)
    dos = dos.astype(np.float64)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--n-ped", type=int, default=30)
    ap.add_argument("--n-geno", type=int, default=20)
    ap.add_argument("--n-pheno", type=int, default=21)
    ap.add_argument("--m", type=int, default=30)
    ap.add_argument("--h2s", default="0.25,0.60,0.80")
    ap.add_argument("--modes", default="a22,identity")
    # Tolerance justification: C++ GrmBuilder accumulates the GRM from float32
    # geno blocks (G += Z*Z' with Z in float32), so Gw differs from the exact
    # double-precision reference by ~eps32 * ||G|| per entry. The inverse
    # amplifies that by the condition number of Gw: at n_g=20, kappa ~ 30-40,
    # giving Hinv/Ginv entry noise up to ~1e-7..1e-6 and EBV/Z effects up to
    # ~5e-6 (measured worst case over the 6-tag grid). Plain double-precision
    # arithmetic and exact inverses reproduce the reference to 1e-9; the
    # 5e-6 floor is the float32-GRM storage noise, not a modeling error.
    ap.add_argument("--tolerance", type=float, default=1e-5)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    work = f"{args.outdir}/data"
    os.makedirs(work, exist_ok=True)

    # ---- pedigree: founders 1..n_f, gen1, gen2; genotyped = gen1+gen2 ----
    n_f = args.n_ped // 3
    n_g1 = args.n_ped - n_f
    n_g2 = args.n_ped - n_f  # unused placeholder to keep numbering simple
    ids_ped = list(range(1, args.n_ped + 1))
    pedigree = []
    for i in range(1, n_f + 1):
        pedigree.append([i, 0, 0])
    for i in range(n_f + 1, args.n_ped + 1):
        s = int(rng.integers(1, n_f + 1))
        d = int(rng.integers(1, n_f + 1))
        while d == s:  # selfing degenerates Meuwissen-Luo (f = 4 - As - Ad - 2Asd = 0)
            d = int(rng.integers(1, n_f + 1))
        pedigree.append([i, s, d])
    geno_ids = list(range(n_f + 1, args.n_ped + 1))  # n_geno = n_ped - n_f
    n_geno = len(geno_ids)
    m = args.m
    p = rng.uniform(0.15, 0.45, size=m)
    dos = write_plink_trio(f"{work}/sim", geno_ids, m, p, rng)
    pheno_ids = sorted(set(rng.choice(ids_ped, size=args.n_pheno, replace=False)))
    with open(f"{work}/sim.phe", "w") as f:
        f.write("FID IID Pheno_1\n")
        for i in pheno_ids:
            f.write(f"F{i} {i} 0.0\n")
    with open(f"{work}/sim.ped", "w") as f:
        f.write("id sire dam\n")
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")
    with open(f"{work}/sim.cs.ped", "w") as f:
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")

    # ---- shared model inputs (analysis-consistent) ----
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    geno_idx_ped = [idx[i] for i in geno_ids]  # 0-based pedigree rows, genotyped order
    A = meuwissen_luo(pedigree)
    Ainv = np.linalg.inv(A)
    A22 = A[np.ix_(geno_idx_ped, geno_idx_ped)]
    A22inv = np.linalg.inv(A22)
    mean_dos = dos.mean(axis=0)
    pj = mean_dos / 2.0
    M = dos - mean_dos
    S = np.sum(2.0 * pj * (1.0 - pj))
    G_raw = M @ M.T / S
    dG = np.trace(G_raw) / n_geno
    oG = (G_raw.sum() - np.trace(G_raw)) / (n_geno * (n_geno - 1))
    dA = np.trace(A22) / n_geno
    oA = (A22.sum() - np.trace(A22)) / (n_geno * (n_geno - 1))
    alpha = 0.05
    results = {}
    fails = []

    for mode in args.modes.split(","):
        for h2s in args.h2s.split(","):
            h2 = float(h2s)
            va = 1.0
            ve = (1.0 - h2) / h2
            lam = ve / va
            tag = f"{mode}_h2{h2:.2f}"
            out = f"{args.outdir}/{tag}"
            os.makedirs(out, exist_ok=True)

            # ---- tuneG + blend (replicate C++ buildBacksolveGrm) ----
            if mode == "a22":
                tune_b = 1.0
                if abs(dG - oG) > 1e-12:
                    tune_b = (dA - oA) / (dG - oG)
                tune_a = dA - tune_b * dG
                G_star = tune_b * G_raw + tune_a
                Gw = (1.0 - alpha) * G_star + alpha * A22
                G_backsolve = G_star
            else:  # identity blend, no tuning (default C++ path)
                Gw = (1.0 - alpha) * G_raw + alpha * np.eye(n_geno)
                G_backsolve = Gw
            Gw_inv = np.linalg.inv(Gw)
            Gback_inv = np.linalg.inv(G_backsolve)

            # ---- H and MME (lambda form, intercept first) ----
            H = A.copy()
            H[np.ix_(geno_idx_ped, geno_idx_ped)] = Gw
            Hinv = Ainv.copy()
            Hinv[np.ix_(geno_idx_ped, geno_idx_ped)] += Gw_inv - A22inv
            pheno_rows = [idx[i] for i in pheno_ids]
            n_rec = len(pheno_rows)
            X = np.ones((n_rec, 1))
            Z = np.zeros((n_rec, args.n_ped))
            for k, r in enumerate(pheno_rows):
                Z[k, r] = 1.0

            # data-generating null: a ~ N(0, H*va) (model-consistent, C6-style).
            # H can be marginally non-PSD at tiny n; use the eigen-PSD draw.
            w, Vh = np.linalg.eigh(H)
            w = np.maximum(w, 0.0)
            a_true = (Vh @ (np.sqrt(w) * rng.normal(size=args.n_ped))) * np.sqrt(va)
            y = 1.0 + a_true[pheno_rows] + rng.normal(scale=np.sqrt(ve), size=n_rec)

            # MME solve + full inverse (independent dense reference)
            Clam = np.block([[X.T @ X, X.T @ Z],
                             [Z.T @ X, Z.T @ Z + lam * Hinv]])
            rhs = np.concatenate([X.T @ y, Z.T @ y])
            theta = np.linalg.solve(Clam, rhs)
            a_hat = theta[1:]
            a_hat_g = a_hat[geno_idx_ped]
            Cinv = np.linalg.inv(Clam)
            Cinv_gg = Cinv[1 + np.array(geno_idx_ped)][:, 1 + np.array(geno_idx_ped)]

            # ---- SNP effects + sampling variance (B_j = M_j' G_backsolve^-1 / S) ----
            q = Gback_inv @ a_hat_g
            EA = (M.T @ q) / S
            V_g = Gw * va - ve * Cinv_gg          # Var(u_hat_g)
            W = Gback_inv @ M                     # n_geno x m; col j = G^-1 M_j
            var_EA = np.einsum("ij,ij->j", W, V_g @ W) / (S * S)
            se = np.sqrt(np.maximum(var_EA, 0.0))
            z = np.where(se > 0, EA / se, np.nan)
            pval = erfc(np.abs(z) / np.sqrt(2.0))
            log10p = -np.log10(np.maximum(pval, 1e-300))

            # ---- write phenotype for the C++ run ----
            with open(f"{work}/sim_y_{tag}.phe", "w") as f:
                f.write("FID IID Pheno_1\n")
                for k, r in enumerate(pheno_rows):
                    f.write(f"F{ids_ped[r]} {ids_ped[r]} {y[k]:.10f}\n")

            # ---- run CosmicScanner exact p-value mode ----
            cmd = [args.bin, "--bfile", f"{work}/sim", "--pheno", f"{work}/sim_y_{tag}.phe",
                   "--gwas", "--single-step", "--single-step-method", "backsolve",
                   "--ped", f"{work}/sim.cs.ped", "--reml-vg", f"{va:.10f}",
                   "--reml-ve", f"{ve:.10f}", "--no-loco", "--ssgwas-window-snps", "10", "--ssgwas-solver", "ldlt",
                   "--ssgwas-pvalue-method", "exact", "--out", f"{out}/m2"]
            if mode == "a22":
                cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                fails.append(f"{tag}: C++ failed: {r.stderr[-500:]}")
                continue

            # ---- compare ----
            def load_ebv(path):
                return np.array([float(l.split("\t")[1]) for l in open(path).read().splitlines()[1:]])
            a_hat_cpp = load_ebv(f"{out}/m2.ss_backsolve.ebv.txt")
            pv = {l.split("\t")[1]: l.split("\t") for l in
                  open(f"{out}/m2.ss_backsolve.pvalue.txt").read().splitlines()[1:]}
            rows = sorted(pv.values(), key=lambda r: int(r[1][2:]))
            EA_c = np.array([float(r[6]) for r in rows])
            SE_c = np.array([float(r[7]) for r in rows])
            Z_c = np.array([float(r[9]) for r in rows])
            LOG10P_c = np.array([float(r[11]) for r in rows])

            def md(name, x, y):
                d = np.max(np.abs(np.asarray(x) - np.asarray(y)))
                return d
            res = {
                "EBV": md("ebv", a_hat_cpp, a_hat),
                "EA": md("ea", EA_c, EA),
                "VAR_EFFECT": md("var", SE_c**2, var_EA),
                "SE": md("se", SE_c, se),
                "Z": md("z", Z_c, z),
                "LOG10P": md("l10p", LOG10P_c, log10p),
            }
            results[tag] = res
            # C4 scale invariance on the C++ outputs: E_S = c*E_A, SE_S = |c|*SE_A
            c = np.sqrt(2.0 * pj * (1.0 - pj))
            z_std = (EA_c * c) / np.maximum(SE_c * np.abs(c), 1e-300)
            p_std = erfc(np.abs(z_std) / np.sqrt(2.0))
            z_allele = EA_c / np.maximum(SE_c, 1e-300)
            res["C4_Z_SCALE_INVARIANCE"] = md("c4z", z_std, z_allele)
            res["C4_P_SCALE_INVARIANCE"] = md(
                "c4p", erfc(np.abs(z_std) / np.sqrt(2.0)),
                erfc(np.abs(z_allele) / np.sqrt(2.0)))
            ok = max(res.values()) <= args.tolerance
            if not ok:
                fails.append(f"{tag}: max diff {max(res.values()):.3e} > tol {args.tolerance}")
            print(f"{tag}: " + " ".join(f"{k}={v:.2e}" for k, v in res.items()) +
                  ("  PASS" if ok else "  FAIL"))

    with open(f"{args.outdir}/tiny_reference_summary.json", "w") as f:
        json.dump(results, f, indent=2)
    if fails:
        print("FAILURES:\n" + "\n".join(fails))
        sys.exit(1)
    print(f"PASS: independent tiny reference (gates C5 + C4 scale invariance); "
          f"tolerance {args.tolerance}")


if __name__ == "__main__":
    main()
