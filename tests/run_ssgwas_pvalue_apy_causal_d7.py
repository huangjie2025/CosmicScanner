#!/usr/bin/env python3
"""Phase D7: APY-approx causal ranking and power characterization.

D7 is NOT calibration validation. It evaluates whether APY-approx preserves
causal SNP prioritization under known causal SNP effects, given the D6/D6.5
raw-P calibration limitations. No production formulas, genomic control, P-value
rescale, empirical correction, or result-dependent core selection are used.
"""
import argparse
import json
import math
import os

import numpy as np
from math import erfc as math_erfc

from run_ssgwas_pvalue_apy_medium_benchmark import (
    ALPHA, TIERS, VE,
    build_g, build_pedigree, build_pedigree_livestocklike, full_metrics,
    k98_of, load_meta, load_pvalue, relationship_matrix, run_cs,
    simulate_dosages, simulate_dosages_livestock, spearman, topk_overlap,
    write_plink_trio)
from run_ssgwas_pvalue_apy_calibration_d65 import (
    draw_null, null_chol, overwrite_bim_chromosomes, per_chrom_counts)

erfc = np.vectorize(math_erfc)


def choose_causal(m, per_ch, n_causal, rng, maf):
    candidates = np.where(maf > 0.05)[0].tolist()
    rng.shuffle(candidates)
    chosen = []
    for j in candidates:
        if all(abs(j - k) > 10 for k in chosen):
            chosen.append(j)
        if len(chosen) == n_causal:
            break
    if len(chosen) < n_causal:
        chosen = candidates[:n_causal]
    chosen = sorted(chosen)
    starts = np.cumsum([0] + per_ch)
    chrom = np.searchsorted(starts[1:], chosen, side="right") + 1
    return np.array(chosen, int), np.array(chrom, int)


def causal_effects(M, causal_idx, rng):
    q = len(causal_idx)
    n_large = max(1, int(round(0.2 * q)))
    n_med = max(1, int(round(0.3 * q)))
    n_small = q - n_large - n_med
    mags = np.array([0.55] * n_large + [0.30] * n_med + [0.15] * n_small)
    rng.shuffle(mags)
    signs = rng.choice([-1.0, 1.0], size=q)
    beta = mags * signs
    score = M[:, causal_idx] @ beta
    v = float(np.var(score))
    if v > 0:
        beta *= math.sqrt(0.50 / v)  # fixed causal component variance target
    return beta


def load_stats(path):
    ids, ea, se, z, lp = load_pvalue(path)
    p = erfc(np.abs(z) / math.sqrt(2.0))
    return {"ids": ids, "EA": ea, "SE": se, "Z": z, "LOG10P": lp, "P": p}


def align2(a, b):
    common = set(a["ids"]) & set(b["ids"])
    ma = np.array([x in common for x in a["ids"]])
    mb = np.array([x in common for x in b["ids"]])
    return ma, mb


def ranking_metrics(ex, method, causal_idx, per_ch, topk=100):
    ma, mb = align2(ex, method)
    e0, e1 = ex["EA"][ma], method["EA"][mb]
    l0, l1 = ex["LOG10P"][ma], method["LOG10P"][mb]
    p1 = method["P"][mb]
    ids = np.array(ex["ids"])[ma]
    pos = np.array([int(s[2:]) - 1 for s in ids])
    causal_set = set(causal_idx.tolist())
    top = {k: set(pos[np.argsort(-l1)[:k]].tolist()) for k in (10, 50, 100)}
    causal_recall = len(causal_set & top[100]) / max(len(causal_set), 1)
    precision = len(causal_set & top[100]) / 100.0
    starts = np.cumsum([0] + per_ch)
    chrom = np.searchsorted(starts[1:], pos, side="right") + 1
    top100_idx = np.argsort(-l1)[:100]
    recovered = 0
    for cj in causal_idx:
        cc = np.searchsorted(starts[1:], cj, side="right") + 1
        ok = False
        for ii in top100_idx:
            if chrom[ii] == cc and abs(int(pos[ii]) - int(cj)) <= 10:
                ok = True
                break
        recovered += int(ok)
    causal_mask = np.array([p in causal_set for p in pos])
    return {
        "EFFECT_R": full_metrics(e0, e1)["pearson"],
        "EFFECT_SPEARMAN": spearman(e0, e1),
        "EFFECT_SLOPE": full_metrics(e0, e1)["slope"],
        "EFFECT_RMSE": full_metrics(e0, e1)["rmse"],
        "TOP10": topk_overlap(l0, l1, 10),
        "TOP50": topk_overlap(l0, l1, 50),
        "TOP100": topk_overlap(l0, l1, 100),
        "CAUSAL_RECALL": causal_recall,
        "PRECISION": precision,
        "WINDOW_RECALL": recovered / max(len(causal_idx), 1),
        "POWER_0.05": float(np.mean(p1[causal_mask] < 0.05)),
        "POWER_0.01": float(np.mean(p1[causal_mask] < 0.01)),
        "POWER_0.001": float(np.mean(p1[causal_mask] < 0.001)),
    }


def status_relative(label, m, r0):
    if label == "R0":
        return "REFERENCE"
    cr_base = max(r0["CAUSAL_RECALL"], 1e-12)
    wr_base = max(r0["WINDOW_RECALL"], 1e-12)
    cr_rel = m["CAUSAL_RECALL"] / cr_base
    wr_rel = m["WINDOW_RECALL"] / wr_base
    if m["TOP100"] >= 0.90 and cr_rel >= 0.80 and wr_rel >= 0.80:
        return "RANKING_PRESERVED"
    if m["TOP100"] >= 0.75 and cr_rel >= 0.60 and wr_rel >= 0.60:
        return "MILD_RANKING_LOSS"
    return "SUBSTANTIAL_RANKING_LOSS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--regime", choices=["stress", "livestock"], default="stress")
    ap.add_argument("--tiers", default="all")
    ap.add_argument("--causal-counts", default="10,50")
    ap.add_argument("--omp", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    regime = "HIGH_DIMENSIONAL_STRESS" if args.regime == "stress" else "LIVESTOCK_LIKE_LD"
    tiers = list(range(3)) if args.tiers == "all" else [int(x) for x in args.tiers.split(",")]
    causal_counts = [int(x) for x in args.causal_counts.split(",")]
    rows, decomp, causal_records = [], [], []

    for ti in tiers:
        n_g, m, n_chr = TIERS[ti]
        n_ped = 3 * n_g // 2
        rng_t = np.random.default_rng(args.seed * 7919 + ti * 104729)
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
        M = dos - dos.mean(axis=0)
        maf = np.minimum(dos.mean(axis=0) / 2.0, 1.0 - dos.mean(axis=0) / 2.0)
        G, _, _ = build_g(dos)
        k98 = k98_of((1.0 - ALPHA) * G + ALPHA * np.eye(n_g))
        work = f"{args.outdir}/data/t{ti}"
        os.makedirs(work, exist_ok=True)
        write_plink_trio(f"{work}/sim", typed_ids, dos, rng_t)
        overwrite_bim_chromosomes(f"{work}/sim", per_ch)
        with open(f"{work}/sim.cs.ped", "w") as f:
            for r in ped:
                f.write(f"{r[0]} {r[1]} {r[2]}\n")
        L = null_chol(ped, typed_idx, dos)
        p_rows = typed_idx  # all observed animals are genotyped; causal genotype defined for every y

        for q in causal_counts:
            rng_c = np.random.default_rng(args.seed * 65537 + ti * 977 + q)
            causal_idx, causal_chr = choose_causal(m, per_ch, q, rng_c, maf)
            beta = causal_effects(M, causal_idx, rng_c)
            y_null = draw_null(L, p_rows, rng_c)
            signal = M[:, causal_idx] @ beta
            y = y_null + signal
            tag = f"t{ti}_ng{n_g}_q{q}"
            with open(f"{work}/sim_y.phe", "w") as f:
                f.write("FID IID Pheno_1\n")
                for k, r in enumerate(p_rows):
                    f.write(f"F{ped[r][0]} {ped[r][0]} {y[k]:.10f}\n")
            for j, b, c in zip(causal_idx, beta, causal_chr):
                causal_records.append({"REGIME": regime, "N": n_g, "SCENARIO": q,
                                       "SNP": f"rs{j + 1}", "EFFECT": float(b),
                                       "CHR": int(c), "WINDOW": 10})
            specs = [("R0", "exact", None, None),
                     ("R1_K98", "apy-exact", k98, 1.0),
                     ("R3_K98", "apy-approx", k98, 1.0),
                     ("R3_full", "apy-approx", n_g, round(n_g / k98, 4))]
            stats = {}
            for label, method, core, ratio in specs:
                out = f"{args.outdir}/{tag}/{label}"
                if not os.path.exists(f"{out}/m2.ss_backsolve.pvalue.txt"):
                    run_cs(args.bin, work, out, method, core=core, ratio=ratio, omp=args.omp)
                meta = load_meta(f"{out}/m2.ss_backsolve.meta.txt")
                if method == "apy-approx":
                    assert meta.get("GENOMIC_CONTROL_APPLIED") == "false"
                stats[label] = load_stats(f"{out}/m2.ss_backsolve.pvalue.txt")
            mets = {label: ranking_metrics(stats["R0"], stats[label], causal_idx, per_ch)
                    for label in ("R0", "R1_K98", "R3_K98", "R3_full")}
            for label in ("R0", "R1_K98", "R3_K98", "R3_full"):
                met = mets[label]
                core = "NA" if label == "R0" else ("K98" if "K98" in label else "full")
                method = {"R0": "Exact_Ginv", "R1_K98": "APY-exact",
                          "R3_K98": "APY-approx", "R3_full": "APY-approx"}[label]
                rows.append({"REGIME": regime, "N": n_g, "SCENARIO": q,
                             "CORE": core, "METHOD": method, "K98": k98,
                             "K98_FRAC": k98 / n_g, **met,
                             "STATUS": status_relative(label, met, mets["R0"])})
            for pair, a, b in (("exact_vs_apy_exact", "R0", "R1_K98"),
                               ("apy_exact_vs_apy_approx", "R1_K98", "R3_K98"),
                               ("exact_vs_apy_approx_full", "R0", "R3_full")):
                ma, mb = align2(stats[a], stats[b])
                decomp.append({"REGIME": regime, "N": n_g, "SCENARIO": q,
                               "PAIR": pair,
                               "EFFECT_R": full_metrics(stats[a]["EA"][ma], stats[b]["EA"][mb])["pearson"],
                               "EFFECT_SLOPE": full_metrics(stats[a]["EA"][ma], stats[b]["EA"][mb])["slope"],
                               "LOG10P_R": full_metrics(stats[a]["LOG10P"][ma], stats[b]["LOG10P"][mb])["pearson"],
                               "LOG10P_SLOPE": full_metrics(stats[a]["LOG10P"][ma], stats[b]["LOG10P"][mb])["slope"]})
            print(f"{regime} {tag}: K98={k98} causal={q} done")

    with open(f"{args.outdir}/d7_causal_summary.tsv", "w") as f:
        f.write("REGIME\tN\tSCENARIO\tCORE\tMETHOD\tEFFECT_R\tEFFECT_SLOPE\tTOP10\tTOP50\tTOP100\tCAUSAL_RECALL\tWINDOW_RECALL\tPOWER_0.05\tPOWER_0.01\tPOWER_0.001\tSTATUS\n")
        for r in rows:
            f.write(f"{r['REGIME']}\t{r['N']}\t{r['SCENARIO']}\t{r['CORE']}\t{r['METHOD']}"
                    f"\t{r['EFFECT_R']:.6f}\t{r['EFFECT_SLOPE']:.6f}\t{r['TOP10']:.3f}\t{r['TOP50']:.3f}\t{r['TOP100']:.3f}\t{r['CAUSAL_RECALL']:.3f}\t{r['WINDOW_RECALL']:.3f}\t{r['POWER_0.05']:.3f}\t{r['POWER_0.01']:.3f}\t{r['POWER_0.001']:.3f}\t{r['STATUS']}\n")
    with open(f"{args.outdir}/d7_error_decomposition.tsv", "w") as f:
        f.write("REGIME\tN\tSCENARIO\tPAIR\tEFFECT_R\tEFFECT_SLOPE\tLOG10P_R\tLOG10P_SLOPE\n")
        for r in decomp:
            f.write(f"{r['REGIME']}\t{r['N']}\t{r['SCENARIO']}\t{r['PAIR']}\t{r['EFFECT_R']:.6f}\t{r['EFFECT_SLOPE']:.6f}\t{r['LOG10P_R']:.6f}\t{r['LOG10P_SLOPE']:.6f}\n")
    with open(f"{args.outdir}/d7_causal_snps.tsv", "w") as f:
        f.write("REGIME\tN\tSCENARIO\tSNP\tEFFECT\tCHR\tWINDOW\n")
        for r in causal_records:
            f.write(f"{r['REGIME']}\t{r['N']}\t{r['SCENARIO']}\t{r['SNP']}\t{r['EFFECT']:.8f}\t{r['CHR']}\t{r['WINDOW']}\n")
    with open(f"{args.outdir}/d7_summary.json", "w") as f:
        json.dump({"rows": rows, "decomposition": decomp,
                   "causal_snps": causal_records}, f, indent=2, default=str)
    print(f"wrote {args.outdir}/d7_causal_summary.tsv")


if __name__ == "__main__":
    main()
