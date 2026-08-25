#!/usr/bin/env python3
"""Phase C1 gate C7: causal simulation for the ssGBLUP SNP P-value path.

Data: y = 1 + a + Xc*b + e, a ~ N(0, H*Va) with H the inverse of the fitted
ssGBLUP precision (A^-1 + [0 0; 0 Gw^-1 - A22^-1]); 10 causal SNPs carry
FIXED effects b_j = +-b (alternating sign; deterministic per-SNP size, the
"QTL" scenario; random per-SNP draws of N(0, sigma^2) are too weak for a
30-test power gate at this sample size); total SNP variance = 15-50% of Va;
e white.

Inference on the same data (identity blend + a22 tuning on both):
  M1: --single-step-method lmm          (fixed-SNP LMM: BETA/SE/Wald-P)
  M2: --single-step-method backsolve --ssgwas-pvalue-method exact
      (back-solved effects with propagated Var(u_hat_g); Aguilar 2019)

Checks over `--reps` replicates:
  * causal recovery: hits in top-10 / top-50, empirical power at
    P < 1e-4 / 1e-3 / 0.01 across all causal tests.
  * M1-vs-M2 agreement: Spearman rho of -log10(P), top-50 overlap.

Usage:
  python3 tests/run_ssgwas_pvalue_causal_sim.py --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_causal --seed 123
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
from scipy import stats

from run_ssgwas_pvalue_tiny_reference import meuwissen_luo, write_plink_trio


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--n-ped", type=int, default=600)
    ap.add_argument("--n-pheno", type=int, default=500)
    ap.add_argument("--m", type=int, default=8000)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--n-causal", type=int, default=10)
    ap.add_argument("--snp-var-frac", type=float, default=0.5)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    work = f"{outdir}/data"
    os.makedirs(work, exist_ok=True)

    n_f = args.n_ped // 5
    ped = [[i, 0, 0] for i in range(1, n_f + 1)]
    for i in range(n_f + 1, args.n_ped + 1):
        s = int(rng.integers(1, n_f + 1))
        d = int(rng.integers(1, n_f + 1))
        while d == s:
            d = int(rng.integers(1, n_f + 1))
        ped.append([i, s, d])
    geno_ids = list(range(n_f + 1, args.n_ped + 1))
    n_g = len(geno_ids)
    m = args.m
    p = rng.uniform(0.15, 0.45, size=m)
    dos = write_plink_trio(f"{work}/sim", geno_ids, m, p, rng)
    with open(f"{work}/sim.bim", "w") as f:
        for j in range(m):
            chr_no = (j % 20) + 1
            f.write(f"{chr_no} rs{j + 1} 0 {j + 1} A G\n")
    with open(f"{work}/sim.cs.ped", "w") as f:
        for r in ped:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")

    A = meuwissen_luo(ped)
    idx = {r[0]: i for i, r in enumerate(ped)}
    gpidx = [idx[i] for i in geno_ids]
    gp = np.ix_(gpidx, gpidx)
    mean_dos = dos.mean(axis=0)
    pj = mean_dos / 2.0
    M = dos - mean_dos
    S = np.sum(2 * pj * (1 - pj))
    G = M @ M.T / S
    Gw = 0.95 * G + 0.05 * np.eye(n_g)
    A22 = A[gp]
    dG = np.trace(G) / n_g
    oG = (G.sum() - np.trace(G)) / (n_g * (n_g - 1))
    dA = np.trace(A22) / n_g
    oA = (A22.sum() - np.trace(A22)) / (n_g * (n_g - 1))
    tune_b = 1.0 if abs(dG - oG) <= 1e-12 else (dA - oA) / (dG - oG)
    G_star = tune_b * G + (dA - tune_b * dG)
    Hinv = np.linalg.inv(A)
    Hinv[gp] += np.linalg.inv(Gw) - np.linalg.inv(A22)
    wa, Vha = np.linalg.eigh(np.linalg.inv(Hinv))
    wa = np.maximum(wa, 0.0)
    Lh = Vha @ np.diag(np.sqrt(wa))
    var_x = 2.0 * pj * (1.0 - pj)
    b_eff = np.sqrt(args.snp_var_frac / (args.n_causal * var_x.mean()))

    causal = sorted(rng.choice(m, size=args.n_causal, replace=False))
    summ = {"causal_snps": [f"rs{j + 1}" for j in causal],
            "beta_fixed": float(b_eff), "reps": {}}
    fails = []
    hit10 = {"M1": 0, "M2": 0}
    hit50 = {"M1": 0, "M2": 0}
    n_pow = {1e-4: 0, 1e-3: 0, 0.01: 0}
    n_pow_m1 = {1e-4: 0, 1e-3: 0, 0.01: 0}
    rhos = []
    ovs = []
    for rep in range(1, args.reps + 1):
        va, ve = 1.0, 1.0
        a_true = Lh @ rng.normal(size=args.n_ped)
        # fixed per-SNP effects, alternating sign: b^2 * var_x summed = frac*Va
        b = np.sqrt(args.snp_var_frac / (args.n_causal * var_x.mean()))
        beta = np.array([b if j % 2 == 0 else -b for j in range(args.n_causal)])
        sge = (dos[:, causal] - 2.0 * pj[causal]) @ beta
        y = 1.0 + a_true + np.concatenate([np.zeros(n_f), sge]) + \
            rng.normal(scale=np.sqrt(ve), size=args.n_ped)
        pheno_ids = sorted(set(rng.choice(list(range(1, args.n_ped + 1)),
                                          size=args.n_pheno, replace=False)))
        with open(f"{work}/sim_y_rep{rep}.phe", "w") as f:
            f.write("FID IID Pheno_1\n")
            for i in pheno_ids:
                f.write(f"F{i} {i} {y[i - 1]:.10f}\n")
        out = f"{outdir}/rep{rep}"
        os.makedirs(out, exist_ok=True)
        base = [args.bin, "--bfile", f"{work}/sim",
                "--pheno", f"{work}/sim_y_rep{rep}.phe", "--gwas",
                "--single-step", "--ped", f"{work}/sim.cs.ped",
                "--reml-vg", f"{va:.10f}", "--reml-ve", f"{ve:.10f}",
                "--no-loco", "--ssgwas-window-snps", "10"]
        cmd2 = base + ["--single-step-method", "backsolve",
                       "--ssgwas-solver", "ldlt", "--ssgwas-pvalue-method", "exact",
                       "--ssgwas-tune-g", "--ssgwas-blend-a22", "--out", f"{out}/m2"]
        r = subprocess.run(cmd2, capture_output=True, text=True)
        if r.returncode != 0:
            fails.append(f"rep{rep} M2 failed: {r.stderr[-500:]}")
            continue
        pv2 = {l.split("\t")[1]: float(l.split("\t")[10])
               for l in open(f"{out}/m2.ss_backsolve.pvalue.txt").read().splitlines()[1:]}
        cmd1 = base + ["--single-step-method", "lmm", "--out", f"{out}/m1"]
        r = subprocess.run(cmd1, capture_output=True, text=True)
        if r.returncode != 0:
            fails.append(f"rep{rep} M1 failed: {r.stderr[-500:]}")
            continue
        pv1 = {l.split("\t")[1]: float(l.split("\t")[9])
               for l in open(f"{out}/m1.assoc.txt").read().splitlines()[1:]}
        all_snps = [f"rs{j + 1}" for j in range(m)]
        l1 = np.array([-np.log10(pv1[s]) for s in all_snps])
        l2 = np.array([-np.log10(pv2[s]) for s in all_snps])
        caus = np.array([all_snps.index(s) for s in summ["causal_snps"]])
        t1 = np.argsort(-l1)
        t2 = np.argsort(-l2)

        def hits(t, k):
            return sum(1 for c in caus if np.where(t == c)[0][0] < k)

        h2_10 = hits(t2, 10)
        h2_50 = hits(t2, 50)
        h1_10 = hits(t1, 10)
        h1_50 = hits(t1, 50)
        hit10["M1"] += h1_10
        hit10["M2"] += h2_10
        hit50["M1"] += h1_50
        hit50["M2"] += h2_50
        for a_t in n_pow:
            n_pow[a_t] += sum(1 for c in caus if pv2[all_snps[c]] < a_t)
            n_pow_m1[a_t] += sum(1 for c in caus if pv1[all_snps[c]] < a_t)
        rho = stats.spearmanr(l1, l2).statistic
        ov = len(set(t1[:50]) & set(t2[:50])) / 50.0
        rhos.append(float(rho))
        ovs.append(float(ov))
        summ["reps"][str(rep)] = {
            "M2_top10": int(h2_10), "M2_top50": int(h2_50),
            "M1_top10": int(h1_10), "M1_top50": int(h1_50),
            "rho_m1_m2": round(float(rho), 6), "top50_overlap": round(ov, 3)}
        print(f"rep{rep}: M2 top10={h2_10} top50={h2_50} | "
              f"M1 top10={h1_10} top50={h1_50} | "
              f"rho={rho:.4f} ov={ov:.3f}")

    total_c = args.n_causal * args.reps
    summ["totals"] = {
        "M2_top10": int(hit10["M2"]), "M2_top50": int(hit50["M2"]),
        "M1_top10": int(hit10["M1"]), "M1_top50": int(hit50["M1"]),
        "power_M2": {str(k): round(v / total_c, 4) for k, v in n_pow.items()},
        "power_M1": {str(k): round(v / total_c, 4) for k, v in n_pow_m1.items()},
        "rho_m1_m2_mean": round(float(np.mean(rhos)), 6),
        "top50_overlap_mean": round(float(np.mean(ovs)), 3)}
    # Gates (UDP): expect the back-solved path to enrich: >= half the causal
    # SNPs in top-50 and >= half detected at P<0.01; M1/M2 -log10(P) ordering
    # should be strongly concordant (rho >= 0.8) with top-50 overlap >= 0.4.
    ok = True
    if hit50["M2"] < 0.5 * total_c:
        fails.append(f"M2 top-50 hits {hit50['M2']}/{total_c} < 50%")
        ok = False
    if n_pow[0.01] < 0.5 * total_c:
        fails.append(f"M2 power@0.01 {n_pow[0.01]}/{total_c} < 50%")
        ok = False
    if float(np.mean(rhos)) < 0.8:
        fails.append(f"rho(M1,M2) {np.mean(rhos):.4f} < 0.8")
        ok = False
    if float(np.mean(ovs)) < 0.4:
        fails.append(f"top50 overlap {np.mean(ovs):.3f} < 0.4")
        ok = False
    with open(f"{outdir}/causal_sim_summary.json", "w") as f:
        json.dump(summ, f, indent=2)
    if fails:
        print("FAILURES:\n" + "\n".join(fails))
        sys.exit(1)
    print(f"PASS: causal simulation (gate C7); {args.reps} reps x "
          f"{args.n_causal} causal SNPs")


if __name__ == "__main__":
    main()