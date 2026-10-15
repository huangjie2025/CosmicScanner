#!/usr/bin/env python3
"""Phase D6.5: paired null-calibration audit for ssGWAS APY P-values.

This confirmatory audit reuses the D5 generators but changes the sampling plan:
independent genotype/pedigree replicates AND independent phenotype replicates.
For each paired replicate it runs exactly four methods on the same data:

  R0  exact
  R1  apy-exact  @ K98
  R3K apy-approx @ K98
  R3F apy-approx @ full core

0.5*K98 is not recomputed here; D6 keeps it as UNDER_CORE_STRESS reference.
The null phenotype is strictly model-consistent:

  Hinv = A^-1 + [0 0; 0 Gw^-1 - A22^-1],  a ~ N(0, Va H),  y = 1 + a + e.

No genomic control, rescale, empirical correction, or result-dependent core
selection is performed.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
from math import erfc as math_erfc

from run_ssgwas_pvalue_apy_medium_benchmark import (
    ALPHA, TIERS, VE,
    build_g, build_pedigree, build_pedigree_livestocklike, k98_of,
    load_meta, load_pvalue, relationship_matrix, run_cs,
    simulate_dosages, simulate_dosages_livestock, write_plink_trio)

erfc = np.vectorize(math_erfc)
CHI2_MEDIAN_1DF = 0.45493642311957277


def parse_counts(s):
    vals = [int(x) for x in s.split(",")]
    assert len(vals) == 3
    return vals


def per_chrom_counts(m, n_chr):
    base = m // n_chr
    per = [base] * n_chr
    for k in range(m - base * n_chr):
        per[k % n_chr] += 1
    return per


def overwrite_bim_chromosomes(prefix, per_ch):
    with open(prefix + ".bim", "w") as f:
        j = 0
        for chrom, n in enumerate(per_ch, start=1):
            for _ in range(n):
                j += 1
                f.write(f"{chrom} rs{j} 0 {j} A G\n")


def null_chol(ped, typed_idx_ped, dos):
    A, _ = relationship_matrix(ped)
    o = np.ix_(typed_idx_ped, typed_idx_ped)
    G, _, _ = build_g(dos)
    Gw = (1.0 - ALPHA) * G + ALPHA * np.eye(G.shape[0])
    Hinv = np.linalg.inv(A)
    Hinv[o] += np.linalg.inv(Gw) - np.linalg.inv(A[o])
    Hinv = 0.5 * (Hinv + Hinv.T)
    return np.linalg.cholesky(Hinv)


def draw_null(L, p_rows, rng):
    z = rng.normal(size=L.shape[0])
    try:
        from scipy.linalg import solve_triangular
        a = solve_triangular(L.T, z, lower=False)
    except ImportError:
        a = np.linalg.solve(L.T, z)
    return 1.0 + a[p_rows] + rng.normal(scale=math.sqrt(VE), size=len(p_rows))


def load_zp(path):
    ids, ea, se, z, lp = load_pvalue(path)
    return {sid: float(zz) for sid, zz in zip(ids, z) if math.isfinite(zz)}


def metrics_from_z(z):
    z = np.asarray(z, float)
    p = erfc(np.abs(z) / math.sqrt(2.0))
    return {
        "n": int(len(z)),
        "lambda": float(np.median(z * z) / CHI2_MEDIAN_1DF),
        "fpr05": float(np.mean(p < 0.05)),
        "fpr01": float(np.mean(p < 0.01)),
        "fpr001": float(np.mean(p < 0.001)),
        "qq_dev": float(np.max(np.abs(np.histogram(p, bins=np.linspace(0, 1, 11))[0] / len(z) - 0.1))),
        "median_z": float(np.median(z)),
        "var_z": float(np.var(z)),
    }


def paired_quant(z0, z1, z3k, z3f):
    m0, m1, m3k, m3f = [metrics_from_z(x) for x in (z0, z1, z3k, z3f)]
    return {
        "lambda_R0": m0["lambda"], "lambda_R1": m1["lambda"],
        "lambda_R3K": m3k["lambda"], "lambda_R3F": m3f["lambda"],
        "ratio_L1": m1["lambda"] / m0["lambda"],
        "ratio_L2": m3k["lambda"] / m1["lambda"],
        "ratio_TOTAL": m3k["lambda"] / m0["lambda"],
        "ratio_FULL": m3f["lambda"] / m0["lambda"],
        "delta_fpr05_total": m3k["fpr05"] - m0["fpr05"],
        "delta_fpr01_total": m3k["fpr01"] - m0["fpr01"],
        "delta_fpr001_total": m3k["fpr001"] - m0["fpr001"],
        "delta_fpr05_L1": m1["fpr05"] - m0["fpr05"],
        "delta_fpr05_L2": m3k["fpr05"] - m1["fpr05"],
        "delta_fpr05_FULL": m3f["fpr05"] - m0["fpr05"],
        "qq_diff_TOTAL": m3k["qq_dev"] - m0["qq_dev"],
        "qq_diff_FULL": m3f["qq_dev"] - m0["qq_dev"],
    }


def cancellation(q):
    l1 = math.log(max(q["ratio_L1"], 1e-300))
    l2 = math.log(max(q["ratio_L2"], 1e-300))
    tot = math.log(max(q["ratio_TOTAL"], 1e-300))
    return bool(abs(tot) < 0.5 * max(abs(l1), abs(l2)) and
                abs(l1) > 0.03 and abs(l2) > 0.03 and l1 * l2 < 0.0)


def verdict(ratio, delta_fpr):
    if ratio >= 1.10 or delta_fpr >= 0.015:
        return "SUBSTANTIAL_INFLATION"
    if ratio >= 1.03 or delta_fpr >= 0.005:
        return "MILD_INFLATION"
    if ratio <= 0.90 or delta_fpr <= -0.015:
        return "MILD_DEFLATION"
    if ratio <= 0.97 or delta_fpr <= -0.005:
        return "MILD_DEFLATION"
    return "CALIBRATED_WITHIN_SIMULATION_ERROR"


def summarise(vals):
    vals = np.asarray(vals, float)
    return {"mean": float(np.mean(vals)), "median": float(np.median(vals)),
            "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
            "ci95": [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]}


def bootstrap(blocks, n_boot, seed):
    rng = np.random.default_rng(seed)
    out = {"ratio_L1": [], "ratio_L2": [], "ratio_TOTAL": [], "ratio_FULL": [],
           "delta_fpr05_total": [], "delta_fpr01_total": [], "delta_fpr001_total": [],
           "qq_diff_TOTAL": [], "qq_diff_FULL": []}
    n = len(blocks)
    for _ in range(n_boot):
        pick = rng.integers(0, n, size=n)
        z0 = np.concatenate([blocks[i]["z0"] for i in pick])
        z1 = np.concatenate([blocks[i]["z1"] for i in pick])
        z3k = np.concatenate([blocks[i]["z3k"] for i in pick])
        z3f = np.concatenate([blocks[i]["z3f"] for i in pick])
        q = paired_quant(z0, z1, z3k, z3f)
        for k in out:
            out[k].append(q[k])
    return {k: summarise(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--regime", choices=["stress", "livestock"], default="stress")
    ap.add_argument("--tiers", default="all")
    ap.add_argument("--gen-reps", default="3,3,2")
    ap.add_argument("--pheno-reps", default="3,3,3")
    ap.add_argument("--n-boot", type=int, default=500)
    ap.add_argument("--omp", type=int, default=8)
    args = ap.parse_args()

    tiers = list(range(3)) if args.tiers == "all" else [int(x) for x in args.tiers.split(",")]
    gen_reps = parse_counts(args.gen_reps)
    pheno_reps = parse_counts(args.pheno_reps)
    regime_name = "HIGH_DIMENSIONAL_STRESS" if args.regime == "stress" else "LIVESTOCK_LIKE_LD"
    os.makedirs(args.outdir, exist_ok=True)
    summary = {"regime": regime_name, "tiers": [], "cells": {}, "replicates": []}

    for ti in tiers:
        n_g, m, n_chr = TIERS[ti]
        n_ped = 3 * n_g // 2
        blocks = []
        reps = []
        for gi in range(1, gen_reps[ti] + 1):
            work = f"{args.outdir}/data/t{ti}_g{gi}"
            os.makedirs(work, exist_ok=True)
            rng_t = np.random.default_rng(args.seed * 7919 + ti * 104729 + gi * 1000003)
            if args.regime == "livestock":
                ped, typed_ids, info = build_pedigree_livestocklike(n_ped, n_g, m, rng_t)
                idx = {r[0]: i for i, r in enumerate(ped)}
                typed_idx = np.array([idx[i] for i in typed_ids])
                dos, per_ch = simulate_dosages_livestock(ped, idx, typed_idx, m, rng_t, info)
            else:
                ped, typed_ids, *_ = build_pedigree(n_ped, rng_t)
                idx = {r[0]: i for i, r in enumerate(ped)}
                typed_idx = np.array([idx[i] for i in typed_ids])
                dos = simulate_dosages(ped, idx, typed_idx, m, n_chr, rng_t)
                per_ch = per_chrom_counts(m, n_chr)
            G, _, _ = build_g(dos)
            Gbs = (1.0 - ALPHA) * G + ALPHA * np.eye(n_g)
            k98 = k98_of(Gbs)
            L = null_chol(ped, typed_idx, dos)
            write_plink_trio(f"{work}/sim", typed_ids, dos, rng_t)
            overwrite_bim_chromosomes(f"{work}/sim", per_ch)
            with open(f"{work}/sim.cs.ped", "w") as f:
                for r in ped:
                    f.write(f"{r[0]} {r[1]} {r[2]}\n")
            n_rec = int(0.7 * n_ped)
            pids = sorted(set(rng_t.choice(list(range(1, n_ped + 1)), size=n_rec, replace=False)))
            p_rows = np.array([idx[i] for i in pids])

            for pi in range(1, pheno_reps[ti] + 1):
                rep_tag = f"t{ti}_g{gi}_p{pi}"
                rng_y = np.random.default_rng(args.seed * 104729 + ti * 101 + gi * 1009 + pi)
                y = draw_null(L, p_rows, rng_y)
                with open(f"{work}/sim_y.phe", "w") as f:
                    f.write("FID IID Pheno_1\n")
                    for k, r in enumerate(p_rows):
                        f.write(f"F{ped[r][0]} {ped[r][0]} {y[k]:.10f}\n")
                out = f"{args.outdir}/t{ti}_ng{n_g}/{rep_tag}"
                specs = [("R0", "exact", None, None),
                         ("R1K", "apy-exact", k98, 1.0),
                         ("R3K", "apy-approx", k98, 1.0),
                         ("R3F", "apy-approx", n_g, round(n_g / k98, 4))]
                zs = {}
                for label, method, core, ratio in specs:
                    d = f"{out}/{label}"
                    if not os.path.exists(f"{d}/m2.ss_backsolve.pvalue.txt"):
                        run_cs(args.bin, work, d, method, core=core, ratio=ratio, omp=args.omp)
                    meta = load_meta(f"{d}/m2.ss_backsolve.meta.txt")
                    if method == "apy-approx":
                        assert meta.get("GENOMIC_CONTROL_APPLIED") == "false"
                    zs[label] = load_zp(f"{d}/m2.ss_backsolve.pvalue.txt")
                common = sorted(set(zs["R0"]) & set(zs["R1K"]) & set(zs["R3K"]) & set(zs["R3F"]),
                                key=lambda s: int(s[2:]))
                chrom = np.array([(int(s[2:]) - 1) for s in common])
                starts = np.cumsum([0] + per_ch)
                chrom = np.searchsorted(starts[1:], chrom, side="right") + 1
                z0 = np.array([zs["R0"][s] for s in common])
                z1 = np.array([zs["R1K"][s] for s in common])
                z3k = np.array([zs["R3K"][s] for s in common])
                z3f = np.array([zs["R3F"][s] for s in common])
                q = paired_quant(z0, z1, z3k, z3f)
                q["replicate"] = rep_tag
                q["n_snp"] = int(len(common))
                q["error_cancellation"] = cancellation(q)
                reps.append(q)
                for c in sorted(set(chrom)):
                    mask = chrom == c
                    blocks.append({"rep": rep_tag, "chrom": int(c), "z0": z0[mask],
                                   "z1": z1[mask], "z3k": z3k[mask], "z3f": z3f[mask]})
                print(f"{regime_name} t{ti} g{gi} p{pi}: K98={k98} n_snp={len(common)} "
                      f"ratio total={q['ratio_TOTAL']:.3f} full={q['ratio_FULL']:.3f}")

        z0 = np.concatenate([b["z0"] for b in blocks])
        z1 = np.concatenate([b["z1"] for b in blocks])
        z3k = np.concatenate([b["z3k"] for b in blocks])
        z3f = np.concatenate([b["z3f"] for b in blocks])
        pooled = paired_quant(z0, z1, z3k, z3f)
        boot = bootstrap(blocks, args.n_boot, args.seed + ti * 17)
        canc = [r["error_cancellation"] for r in reps]
        cell = {"regime": regime_name, "tier": f"t{ti}_ng{n_g}", "n_g": n_g, "M": m,
                "K98": k98, "K98_FRAC": k98 / n_g,
                "replicates": len(reps), "genotype_reps": gen_reps[ti],
                "phenotype_reps": pheno_reps[ti], "pooled": pooled,
                "bootstrap": boot, "replicate_summary": {k: summarise([r[k] for r in reps])
                    for k in ("ratio_L1", "ratio_L2", "ratio_TOTAL", "ratio_FULL",
                              "delta_fpr05_total", "delta_fpr01_total", "delta_fpr001_total",
                              "qq_diff_TOTAL", "qq_diff_FULL")},
                "error_cancellation_fraction": float(np.mean(canc)),
                "error_cancellation_any": bool(any(canc)),
                "verdict_total": verdict(pooled["ratio_TOTAL"], pooled["delta_fpr05_total"]),
                "verdict_full": verdict(pooled["ratio_FULL"], pooled["delta_fpr05_FULL"])}
        summary["tiers"].append(cell)
        summary["replicates"].extend(reps)

    with open(f"{args.outdir}/d65_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    with open(f"{args.outdir}/d65_paired_table.tsv", "w") as f:
        f.write("REGIME\tN_G\tK98_FRAC\tCORE\tREPLICATES\tLAMBDA_R0\tLAMBDA_R1\tLAMBDA_R3"
                "\tLAMBDA_RATIO_L1\tLAMBDA_RATIO_L2\tLAMBDA_RATIO_TOTAL\tLAMBDA_RATIO_FULL"
                "\tDELTA_FPR05_TOTAL\tDELTA_FPR01_TOTAL\tDELTA_FPR001_TOTAL"
                "\tERROR_CANCELLATION\t95CI_TOTAL\t95CI_FULL\tVERDICT\n")
        for c in summary["tiers"]:
            p = c["pooled"]
            bt = c["bootstrap"]["ratio_TOTAL"]["ci95"]
            bf = c["bootstrap"]["ratio_FULL"]["ci95"]
            for core, v in (("K98", c["verdict_total"]), ("full", c["verdict_full"])):
                f.write(f"{c['regime']}\t{c['n_g']}\t{c['K98_FRAC']:.4f}\t{core}\t{c['replicates']}"
                        f"\t{p['lambda_R0']:.4f}\t{p['lambda_R1']:.4f}\t"
                        f"{(p['lambda_R3K'] if core == 'K98' else p['lambda_R3F']):.4f}"
                        f"\t{p['ratio_L1']:.4f}\t{p['ratio_L2']:.4f}\t{p['ratio_TOTAL']:.4f}\t{p['ratio_FULL']:.4f}"
                        f"\t{p['delta_fpr05_total']:.5f}\t{p['delta_fpr01_total']:.5f}\t{p['delta_fpr001_total']:.6f}"
                        f"\t{c['error_cancellation_any']}\t[{bt[0]:.4f},{bt[1]:.4f}]\t[{bf[0]:.4f},{bf[1]:.4f}]\t{v}\n")
    print(f"wrote {args.outdir}/d65_summary.json")


if __name__ == "__main__":
    main()
