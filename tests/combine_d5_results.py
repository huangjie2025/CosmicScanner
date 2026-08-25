#!/usr/bin/env python3
"""D5 post-processing: combine stress + livestock summaries into FINAL D5
table (per user spec) + core/full-core curves + Q1-Q6 answers."""
import json
import glob
import math
import os
import sys


def load(path):
    with open(path) as f:
        return json.load(f)


def main():
    base = "/disk228/huangjie/PanOmics/CosmicScanner/tests/results"
    dirs = [f"{base}/ssgwas_pvalue_apy_medium",
            f"{base}/ssgwas_pvalue_apy_medium_livestock"]
    summaries = []
    for d in dirs:
        p = f"{d}/medium_summary.json"
        if os.path.exists(p):
            s = load(p)
            if s["tiers"]:
                summaries.append(s)
    if not summaries:
        print("no summaries yet")
        sys.exit(1)

    out = {"D5_final_table": [], "core_curve": [], "full_core_curves": [],
           "Q1_Q6": {}}
    for s in summaries:
        for t in s["tiers"]:
            row = {"REGIME": s["regime"], "N_G": t["n_g"], "M": t["M"],
                   "K90": t["K90"], "K95": t["K95"], "K98": t["K98"],
                   "K99": t["K99"], "K98_FRAC": t["K98_N_G"],
                   "EFFECTIVE_RANK": round(t["effective_rank"], 1),
                   "MAF_MEDIAN": t["MAF"]["median"],
                   "MAF_MONO": t["MAF"]["monomorphic"],
                   "LD_R2_25BP": t["LD_DECAY_MEAN_R_BY_BP"].get("25"),
                   "LD_R2_MAX": max(t["LD_DECAY_MEAN_R_BY_BP"].values())}
            out["D5_final_table"].append(row)
        for k, v in s["classifications"].items():
            row = {"REGIME": s["regime"], "CASE": k, "CORE": v["core"],
                   "VERDICT": v["verdict"]}
            out["D5_final_table"].append(row)

    # per-core curves from TSV for the two regimes (merge rows)
    import csv
    for d in dirs:
        p = f"{d}/core_curve.tsv"
        if not os.path.exists(p):
            continue
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                out["core_curve"].append(r)

    # Q1: full-core Z/LP slopes by regime over tiers
    for s in summaries:
        reg = s["regime"]
        fc = [t["N_dep_full_core"] for t in s["tiers"]]
        out["Q1_Q6"].setdefault(reg, {})["full_core_curve"] = [
            {"N_G": t["n_g"], "SE_slope": t["N_dep_full_core"]["SE_slope"],
             "Z_slope": t["N_dep_full_core"]["Z_slope"],
             "P_slope": t["N_dep_full_core"]["P_slope"],
             "SE_median_ratio": t["N_dep_full_core"]["SE_median_ratio"]}
            for t in s["tiers"]]
        out["Q1_Q6"][reg]["Q"] = s["q"]

    with open("/disk228/huangjie/PanOmics/CosmicScanner/tests/results/d5_final_combined.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print("wrote d5_final_combined.json")
    for reg, v in out["Q1_Q6"].items():
        print(reg, "full-core:", v["full_core_curve"])
        print("  Q1:", v["Q"]["Q1"])
        print("  Q6:", v["Q"]["Q6"])


if __name__ == "__main__":
    main()