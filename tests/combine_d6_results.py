#!/usr/bin/env python3
"""Combine the two per-regime D6 calibration summaries into one final report:
cells table + pairwise attribution + Q1-Q6 answers + overall verdicts."""
import json
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RES = f"{BASE}/results"


def load(regime_name, outdir):
    with open(f"{RES}/{outdir}/d6_summary.json") as f:
        return json.load(f)


STRESS = load("stress", "ssgwas_pvalue_apy_d6")
LIVESTOCK = load("livestock", "ssgwas_pvalue_apy_d6_livestock")


def cal_row(s, key):
    c = s["cells"][key]
    return (c["lambda_pooled"], c["fpr05"], c.get("lambda_se", float("nan")),
            c.get("fpr05_se", float("nan")), c["n_tests"], c["verdict"])


final = {"Q1": {}, "Q2": {}, "Q3": {}, "Q4": {}, "Q5": {}, "Q6": {},
         "verdicts": {"stress": {}, "livestock": {}}, "cells": {},
         "pairs": {}}

REG_KEY = {"HIGH_DIMENSIONAL_STRESS": "stress", "LIVESTOCK_LIKE_LD": "livestock"}


def regkey(regime):
    return REG_KEY[regime]


for s in (STRESS, LIVESTOCK):
    reg = regkey(s["regime"])
    for ti, t in enumerate(s["tiers"]):
        tag = t["tier"]
        k = t["K98"]
        h = int(round(0.5 * k))
        n_g = t["n_g"]
        r0 = cal_row(s, f"{tag}_R0_exact")
        r1k = cal_row(s, f"{tag}_R1_apy-exact_core{k}")
        r3k = cal_row(s, f"{tag}_R3_apy-approx_core{k}")
        r1half = cal_row(s, f"{tag}_R1_apy-exact_core{h}")
        r3half = cal_row(s, f"{tag}_R3_apy-approx_core{h}")
        r1full = cal_row(s, f"{tag}_R1_apy-exact_core{n_g}")
        r3full = cal_row(s, f"{tag}_R3_apy-approx_core{n_g}")
        row = {"n_g": n_g, "M": t["M"], "K98": k, "K98_N_G": t["K98_N_G"],
               "R0_exact": {"lambda": r0[0], "FPR05": r0[1], "verdict": r0[5]},
               "R1_K98": {"lambda": r1k[0], "FPR05": r1k[1],
                          "verdict": r1k[5]},
               "R3_K98": {"lambda": r3k[0], "FPR05": r3k[1],
                          "lambda_se": r3k[2], "fpr05_se": r3k[3],
                          "n_tests": r3k[4], "verdict": r3k[5]},
               "R1_halfK98": {"lambda": r1half[0], "FPR05": r1half[1],
                              "verdict": r1half[5]},
               "R3_halfK98": {"lambda": r3half[0], "FPR05": r3half[1],
                              "verdict": r3half[5]},
               "R1_full": {"lambda": r1full[0], "FPR05": r1full[1],
                           "verdict": r1full[5]},
               "R3_full": {"lambda": r3full[0], "FPR05": r3full[1],
                           "verdict": r3full[5]}}
        final["cells"][f"{reg.split('_')[0]}_t{ti}"] = row
        final["verdicts"][reg][f"t{ti}"] = row

# ---------- Q answers ----------
def lam(rowkey, cell):
    return final["cells"][rowkey][cell]["lambda"]


def fpr(rowkey, cell):
    return final["cells"][rowkey][cell]["FPR05"]


# Q1: APY-approx null calibration at K98 (primary) and full.
q1 = {"stress": {}, "livestock": {}}
for reg in ("stress", "livestock"):
    for ti in (0, 1, 2):
        k = f"{reg}_t{ti}"
        q1[reg][f"t{ti}"] = {
            "R3_lambda_K98": lam(k, "R3_K98"),
            "R3_lambda_full": lam(k, "R3_full"),
            "R3_FPR05_K98": fpr(k, "R3_K98"),
            "R3_FPR05_full": fpr(k, "R3_full"),
            "verdict_K98": final["verdicts"][reg][f"t{ti}"]["R3_K98"][
                "verdict"],
            "verdict_full": final["verdicts"][reg][f"t{ti}"]["R3_full"][
                "verdict"]}
final["Q1"] = q1

# Q2: K98 vs 0.5*K98.
q2 = {}
for reg in ("stress", "livestock"):
    q2[reg] = {}
    for ti in (0, 1, 2):
        k = f"{reg}_t{ti}"
        q2[reg][f"t{ti}"] = {
            "delta_lambda_K98_minus_half":
                lam(k, "R3_K98") - lam(k, "R3_halfK98"),
            "R1_K98": lam(k, "R1_K98"), "R1_halfK98": lam(k, "R1_halfK98"),
            "R3_K98": lam(k, "R3_K98"),
            "R3_halfK98": lam(k, "R3_halfK98")}
final["Q2"] = q2

# Q3: stress high K98/N -> inflation? (R0 reference plus R3@K98/full, both
# regimes side by side; per-t0/t1 the LD-robust SE is available.)
q3 = {}
for ti in (0, 1, 2):
    q3[f"t{ti}"] = {
        "stress_K98_N": final["cells"][f"stress_t{ti}"]["K98_N_G"],
        "livestock_K98_N": final["cells"][f"livestock_t{ti}"]["K98_N_G"],
        "stress_R3_lambda_K98": lam("stress_t" + str(ti), "R3_K98"),
        "livestock_R3_lambda_K98": lam("livestock_t" + str(ti), "R3_K98"),
        "stress_R0_lambda": lam("stress_t" + str(ti), "R0_exact"),
        "livestock_R0_lambda": lam("livestock_t" + str(ti), "R0_exact"),
        "stress_R3_lambda_full": lam("stress_t" + str(ti), "R3_full"),
        "livestock_R3_lambda_full": lam("livestock_t" + str(ti), "R3_full")}
final["Q3"] = q3

# Q4: livestock-like closer to practical APY?
q4 = {}
for ti in (0, 1, 2):
    q4[f"t{ti}"] = {
        "stress_lambda_3SE_ok_K98": bool(
            abs(lam("stress_t" + str(ti), "R3_K98") - 1.0) <= 3 *
            final["cells"][f"stress_t{ti}"]["R3_K98"]["lambda_se"]),
        "livestock_lambda_3SE_ok_K98": bool(
            abs(lam("livestock_t" + str(ti), "R3_K98") - 1.0) <= 3 *
            final["cells"][f"livestock_t{ti}"]["R3_K98"]["lambda_se"])}
final["Q4"] = q4

# Q5: L1 vs L2 contribution (lambda attribution at K98, from pairs).
q5 = {}
for reg in ("stress", "livestock"):
    s = STRESS if reg == "stress" else LIVESTOCK
    q5[reg] = {}
    for ti, t in enumerate(s["tiers"]):
        tag = t["tier"]
        k = t["K98"]
        p = s["pairs"][f"{tag}_core{k}"]
        q5[reg][f"t{ti}"] = {"lambda_share_L1": p["lambda_share_L1"],
                             "lambda_share_L2": p["lambda_share_L2"],
                             "lambda_R0": p["lambda_R0"],
                             "lambda_R1": p["lambda_R1"],
                             "lambda_R3": p["lambda_R3"]}
final["Q5"] = q5

# Q6: localization verdict.
r3k_stress = [final["verdicts"]["stress"][f"t{ti}"]["R3_K98"]["verdict"]
              for ti in (0, 1, 2)]
r3k_livestock = [final["verdicts"]["livestock"][f"t{ti}"]["R3_K98"][
    "verdict"] for ti in (0, 1, 2)]
r1k_stress = [final["verdicts"]["stress"][f"t{ti}"]["R1_K98"]["verdict"]
              for ti in (0, 1, 2)]
r1k_livestock = [final["verdicts"]["livestock"][f"t{ti}"]["R1_K98"][
    "verdict"] for ti in (0, 1, 2)]
final["Q6"] = {
    "R3_K98_verdicts": {"stress": r3k_stress, "livestock": r3k_livestock},
    "R1_K98_verdicts": {"stress": r1k_stress, "livestock": r1k_livestock},
    "localization": ("APY-approx = EXPERIMENTAL approximate uncertainty "
                     "layer: (i) L1 (APY G^-1 + core-only operator) inflates "
                     "severely at 0.5*K98 (lambda 1.4-1.9) and mildly at K98;"
                     " (ii) L2 (Eq.8/9) persists at full core with SE "
                     "underestimation -> lambda 1.15-1.21 (stress), "
                     "1.03-1.05 (livestock); (iii) R0 exact within +-8% "
                     "drift. No single calibration formula is defendable "
                     "from this data; no GC/rescale applied."),
    "overall": ("CONCERN: approximation bias is real and regime/core-"
                "dependent; NOT a numerical failure; no production formula "
                "changed.")}

with open(f"{RES}/ssgwas_pvalue_apy_d6/d6_final_combined.json", "w") as f:
    json.dump(final, f, indent=2, default=str)

# ---------- TSV table ----------
with open(f"{RES}/ssgwas_pvalue_apy_d6/d6_calibration_table.tsv", "w") as f:
    f.write("REGIME\tTIER\tN_G\tM\tK98\tK98_FRAC\tMETHOD\tCORE\tCORE_FRAC"
            "\tlambda_GC\tFPR@0.05\tFPR@0.01\tFPR@0.001\tmedian_Z\tVar_Z"
            "\tQQ_dev\tKS_D\tVERDICT\n")
    for s in (STRESS, LIVESTOCK):
        reg = regkey(s["regime"])
        for t in s["tiers"]:
            tag = t["tier"]
            n_g = t["n_g"]
            k = t["K98"]
            for key in s["cells"]:
                if not key.startswith(tag):
                    continue
                c = s["cells"][key]
                meth = key.split("_")[2:]
                core = "NA"
                if "core" in key:
                    core = key.split("core")[1]
                mname = ("R0_exact" if "R0" in key else
                         "R1_apy-exact" if "R1" in key else "R3_apy-approx")
                frac = ("" if core == "NA"
                        else f"\t{int(core) / n_g:.3f}")
                f.write(f"{reg}\t{tag}\t{n_g}\t{t['M']}\t{k}\t{k / n_g:.4f}"
                        f"\t{mname}\t{core}{frac}"
                        f"\t{c['lambda_pooled']:.4f}\t{c['fpr05']:.5f}"
                        f"\t{c['fpr01']:.5f}\t{c['fpr001']:.6f}"
                        f"\t{c['median_Z']:+.4f}\t{c['Var_Z']:.4f}"
                        f"\t{c['QQ_dev']:.4f}\t{c['KS_D']:.4f}"
                        f"\t{c['verdict']}\n")

print("wrote d6_final_combined.json + d6_calibration_table.tsv")
print(json.dumps(final["Q1"], indent=1))
print(json.dumps(final["Q3"], indent=1))
print(json.dumps(final["Q6"], indent=1))