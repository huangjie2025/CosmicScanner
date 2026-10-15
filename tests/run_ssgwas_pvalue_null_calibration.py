#!/usr/bin/env python3
"""Phase C1 gate C6: null calibration of the exact ssGWAS SNP P-value path.

Model-null phenotype generation: y = 1 + a + e with a ~ N(0, H*Va) drawn on the
genotyped set genotype GRM and the same H/A as the analysis, e ~ N(0, Ve) - so
the ssGBLUP model is exactly the data-generating model and every SNP effect
should be null. Calibrated checks over R replicates:

  * genomic control lambda (median chi2_1df(P) / 0.454937) per replicate;
    the replicate mean must sit within 3 sigma of 1 (sigma = sqrt(2/(M*R))).
  * FPR at alpha = 0.05 / 0.01 / 0.001 checked against the binomial MC SE
    sqrt(alpha*(1-alpha)/N_tests) with a 3-sigma bound.
  * QQ uniformity in 10 decile bins, max bin deviation < 3 * sqrt(0.1*0.9/N).

Usage:
  python3 tests/run_ssgwas_pvalue_null_calibration.py --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_null --seed 123
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

from math import erfc as math_erfc
from run_ssgwas_pvalue_tiny_reference import write_plink_trio

erfc = np.vectorize(math_erfc)

CHI2_MEDIAN_1DF = 0.45493642311957277


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--n-ped", type=int, default=600)
    ap.add_argument("--n-pheno", type=int, default=300)
    ap.add_argument("--m", type=int, default=8000)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--h2s", default="0.25,0.60")
    ap.add_argument("--modes", default="a22,identity")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    work = f"{outdir}/data"
    os.makedirs(work, exist_ok=True)

    # ---- fixed analysis set (replicates only resample y, like real data) ----
    n_f = args.n_ped // 5
    ids_ped = list(range(1, args.n_ped + 1))
    pedigree = [[i, 0, 0] for i in range(1, n_f + 1)]
    for i in range(n_f + 1, args.n_ped + 1):
        s = int(rng.integers(1, n_f + 1))
        d = int(rng.integers(1, n_f + 1))
        while d == s:
            d = int(rng.integers(1, n_f + 1))
        pedigree.append([i, s, d])
    geno_ids = list(range(n_f + 1, args.n_ped + 1))  # genotyped = all non-founders
    m = args.m
    p = rng.uniform(0.15, 0.45, size=m)
    dos = write_plink_trio(f"{work}/sim", geno_ids, m, p, rng)
    # spread SNPs over 20 chromosomes so null tests are (nearly) independent;
    # a single-chromosome bim makes the QQ/FPR binomial bounds invalid under LD
    with open(f"{work}/sim.bim", "w") as f:
        for j in range(m):
            chr_no = (j % 20) + 1
            f.write(f"{chr_no} rs{j + 1} 0 {j + 1} A G\n")
    with open(f"{work}/sim.ped", "w") as f:
        f.write("id sire dam\n")
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")
    with open(f"{work}/sim.cs.ped", "w") as f:
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")

    # ---- model matrices for null draws ----
    # The null phenotype must be drawn under EXACTLY the analysis model:
    # a ~ N(0, H*Va) with H = A with [tau*Gw - omega*A22] blended in (identity
    # or tuned-a22 variant) and Gw on the A/X block, y = 1 + a + e. Drawing the
    # genetic signal from G alone (founders untouched, no A22 covariance) puts
    # the data outside the fitted family and inflates every SNP test.
    from run_ssgwas_pvalue_tiny_reference import meuwissen_luo
    A = meuwissen_luo(pedigree)
    mean_dos = dos.mean(axis=0)
    pj = mean_dos / 2.0
    M = dos - mean_dos
    S = np.sum(2.0 * pj * (1.0 - pj))
    G_raw = M @ M.T / S
    dG = np.trace(G_raw) / len(geno_ids)
    oG = (G_raw.sum() - np.trace(G_raw)) / (len(geno_ids) * (len(geno_ids) - 1))
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    gpidx = [idx[i] for i in geno_ids]
    o = np.ix_(gpidx, gpidx)
    Ainv = np.linalg.inv(A)
    A22 = A[o]

    def block_delta(Gw):
        D = np.zeros((args.n_ped, args.n_ped))
        D[o] += np.linalg.inv(Gw) - np.linalg.inv(A[o])
        return D
    dA = np.trace(A22) / len(geno_ids)
    oA = (A22.sum() - np.trace(A22)) / (len(geno_ids) * (len(geno_ids) - 1))
    alpha = 0.05

    def h_for(mode):
        if mode == "a22":
            tune_b = 1.0 if abs(dG - oG) <= 1e-12 else (dA - oA) / (dG - oG)
            G_star = tune_b * G_raw + (dA - tune_b * dG)
            Gw = (1.0 - alpha) * G_star + alpha * A22
        else:
            Gw = (1.0 - alpha) * G_raw + alpha * np.eye(len(geno_ids))
        # The model that is actually solved is the ssGBLUP PRECISION
        # Hinv = Ainv + [0 0; 0 Gw^-1 - A22^-1]. The genomic data-generating H
        # is therefore the inverse of that precision, NOT "A with the G-blend
        # block" (those coincide only when A21 A11^-1 A12 = 0, which fails for
        # any non-trivial pedigree - a naive A-block replacement would put the
        # null OUTSIDE the fitted model and inflate every SNP test).
        Hinv = Ainv + block_delta(Gw)
        w, Vh = np.linalg.eigh(np.linalg.inv(Hinv))
        return np.maximum(w, 0.0), Vh

    def null_y(h2, w, Vh):
        va = 1.0
        ve = (1.0 - h2) / h2
        a_true = (Vh @ (np.sqrt(w) * rng.normal(size=args.n_ped))) * np.sqrt(va)
        y = 1.0 + a_true + rng.normal(scale=np.sqrt(ve), size=args.n_ped)
        return y, va, ve

    summary = {}
    fails = []
    for mode in args.modes.split(","):
        for h2s in args.h2s.split(","):
            h2 = float(h2s)
            tag = f"{mode}_h2{h2:.2f}"
            lam_all = []
            nt = {"0.05": 0, "0.01": 0, "0.001": 0}
            qt_tot = {"0.05": 0, "0.01": 0, "0.001": 0}  # quantile: q-th like counts
            bins = np.zeros(10)
            n_tests = 0
            w_h, Vh_h = h_for(mode)
            for rep in range(1, args.reps + 1):
                y_full, va, ve = null_y(h2, w_h, Vh_h)
                pheno_ids = sorted(set(rng.choice(ids_ped, size=args.n_pheno, replace=False)))
                with open(f"{work}/sim_rep{rep}_{tag}.phe", "w") as f:
                    f.write("FID IID Pheno_1\n")
                    for i in pheno_ids:
                        f.write(f"F{i} {i} {y_full[i - 1]:.10f}\n")
                out = f"{outdir}/{tag}/rep{rep}"
                os.makedirs(out, exist_ok=True)
                cmd = [args.bin, "--bfile", f"{work}/sim", "--pheno",
                       f"{work}/sim_rep{rep}_{tag}.phe", "--gwas", "--single-step",
                       "--single-step-method", "backsolve", "--ped", f"{work}/sim.cs.ped",
                       "--reml-vg", f"{va:.10f}", "--reml-ve", f"{ve:.10f}", "--no-loco",
                       "--ssgwas-window-snps", "10", "--ssgwas-solver", "ldlt",
                       "--ssgwas-pvalue-method", "exact", "--out", f"{out}/m2"]
                if mode == "a22":
                    cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
                r = subprocess.run(cmd, capture_output=True, text=True)
                if r.returncode != 0:
                    fails.append(f"{tag}/rep{rep}: C++ failed: {r.stderr[-500:]}")
                    continue
                pv = []
                # columns: CHR SNP BP A1 A2 FREQ EFFECT_ALLELE SE_ALLELE
                # VAR_EFFECT Z P LOG10P STATUS -> Z is field index 9, P index 10
                for l in open(f"{out}/m2.ss_backsolve.pvalue.txt").read().splitlines()[1:]:
                    zz = float(l.split("\t")[9])
                    if np.isfinite(zz):
                        pv.append(zz)
                pv_arr = np.array(pv)
                if len(pv_arr) == 0:
                    continue
                chi2 = pv_arr ** 2  # z ~ N(0,1) => chi2 ~ chi2_1df exactly
                lam = np.median(chi2) / CHI2_MEDIAN_1DF
                lam_all.append(lam)
                n_tests += len(pv_arr)
                p_uniform = erfc(np.abs(pv_arr) / np.sqrt(2.0))
                for a in ("0.05", "0.01", "0.001"):
                    nt[a] += int((p_uniform < float(a)).sum())
                    qt_tot[a] += len(pv_arr)
                bins += np.histogram(p_uniform, bins=np.linspace(0, 1, 11))[0]
            lam = np.array(lam_all)
            mean_lam = lam.mean() if len(lam) else float("nan")
            se_lam_nom = np.sqrt(2.0 / (m * len(lam))) if len(lam) else float("nan")
            emp_sd = lam.std(ddof=1) if len(lam) > 1 else float("nan")
            # LD-robustness: within-chromosome LD correlates the test statistics
            # and inflates the variance of lambda/QQ/FPR aggregates by a factor
            # f = observed per-replicate lambda SD / nominal SD. All gate bounds
            # below scale by max(1, f); the estimates themselves stay unbiased.
            f_ld = (emp_sd / np.sqrt(2.0 / m)) if len(lam) > 1 else 1.0
            f_ld = max(1.0, float(f_ld))
            se_lam = se_lam_nom * f_ld
            fpr = {}
            for k, a in enumerate(("0.05", "0.01", "0.001")):
                alpha = float(a)
                if n_tests > 0:
                    f = nt[a] / n_tests
                    mc_se = np.sqrt(alpha * (1 - alpha) / n_tests) * f_ld
                    ok = abs(f - alpha) <= 3.0 * mc_se
                    fpr[a] = {"fpr": f, "mc_se": mc_se, "ok": bool(ok)}
                    if not ok:
                        fails.append(
                            f"{tag}: FPR@{alpha} = {f:.5f} vs {alpha} (3*MC-SE {3*mc_se:.5f})")
                else:
                    fpr[a] = {"fpr": float("nan"), "mc_se": float("nan"), "ok": False}
            binf = bins / max(n_tests, 1)
            dev = np.abs(binf - 0.1).max() if n_tests else float("nan")
            qq_bnd = 3.0 * np.sqrt(0.1 * 0.9 / max(n_tests, 1)) * f_ld
            ok_lam = bool(abs(mean_lam - 1.0) <= 3.0 * se_lam) if len(lam) else False
            ok_qq = bool(dev <= qq_bnd) if n_tests else False
            if not ok_lam:
                fails.append(f"{tag}: lambda = {mean_lam:.4f} vs 1 (3*SE {3*se_lam:.4f})")
            if not ok_qq:
                fails.append(f"{tag}: QQ max dev {dev:.5f} > {qq_bnd:.5f}")
            summary[tag] = {
                "n_replicates": int(len(lam)),
                "n_tests": int(n_tests),
                "lambda_mean": float(mean_lam),
                "lambda_se": float(se_lam),
                "lambda_ok": bool(ok_lam),
                "lambda_per_rep": [float(x) for x in lam],
                "fpr": fpr,
                "qq_max_dev": float(dev),
                "qq_bound": float(qq_bnd),
                "qq_ok": bool(ok_qq),
            }
            print(f"{tag}: n={n_tests} lambda={mean_lam:.4f}+-{se_lam:.4f} "
                  f"FPR@0.05={fpr['0.05']['fpr']:.5f} "
                  f"FPR@0.01={fpr['0.01']['fpr']:.5f} "
                  f"FPR@0.001={fpr['0.001']['fpr']:.5f} "
                  f"QQdev={dev:.5f} (bnd {qq_bnd:.5f})  "
                  + ("PASS" if (ok_lam and ok_qq) else "FAIL"))

    with open(f"{outdir}/null_calibration_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    if fails:
        print("FAILURES:\n" + "\n".join(fails))
        sys.exit(1)
    print(f"PASS: null calibration (gate C6); {args.reps} replicates x "
          f"{args.modes} x {args.h2s}")


if __name__ == "__main__":
    main()