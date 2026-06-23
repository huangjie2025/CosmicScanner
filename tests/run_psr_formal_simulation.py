#!/usr/bin/env python3
"""Configurable PSR-GWAS formal simulation driver.

Defaults are intentionally small so the script is safe to smoke-test. Set
PSR_FORMAL_REPS, PSR_FORMAL_N, and PSR_FORMAL_M for paper-level runs.
"""

import csv
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_psr_simulation_summary as smoke


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "tests" / "results" / "psr_formal_simulation"
WORK_DIR = OUT_DIR / "work"
ENGINE = "cosmicscanner_psr"
BASE_SEED = int(os.environ.get("PSR_FORMAL_SEED", "20260624"))


def norm_ppf(p):
    # Peter J. Acklam's rational approximation; sufficient for lambda GC summaries.
    if p <= 0.0 or p >= 1.0:
        return math.nan
    a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.3577518672690, -30.66479806614716, 2.506628277459239]
    b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572]
    c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783]
    d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416]
    plow = 0.02425
    phigh = 1.0 - plow
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


def lambda_gc(pvals):
    chisq = []
    for p in pvals:
        p = min(max(float(p), 1e-300), 1.0)
        z = norm_ppf(1.0 - p / 2.0)
        if math.isfinite(z):
            chisq.append(z * z)
    if not chisq:
        return math.nan
    return float(np.median(chisq) / 0.454936423119572)


def scenario_effect(scenario):
    return 0.0 if scenario.startswith("null") else 0.35 if scenario.endswith("low") else 0.7


def run_one(seed, trait_type, scenario, n, m, rep):
    rng = np.random.default_rng(seed)
    effect = scenario_effect(scenario)
    prefix = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}"
    out = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}.cosmic"
    if trait_type == "ordinal":
        G, causal, payload = simulate_ordinal_data(rng, n, m, effect)
    else:
        G, causal, payload = smoke.simulate_data(rng, n, m, trait_type, effect)
    if trait_type == "binary" and scenario == "null_unbalanced":
        payload["pheno"] = rng.binomial(1, 0.1, size=n)
        if payload["pheno"].min() == payload["pheno"].max():
            payload["pheno"][0] = 0
            payload["pheno"][-1] = 1
    if trait_type == "binary" and scenario == "null_balanced":
        payload["pheno"] = np.r_[np.zeros(n // 2, dtype=int), np.ones(n - n // 2, dtype=int)]
        rng.shuffle(payload["pheno"])
    if trait_type == "poisson" and scenario == "null_low_mean":
        payload["pheno"] = rng.poisson(0.2, size=n)
    if trait_type == "poisson" and scenario == "null_high_mean":
        payload["pheno"] = rng.poisson(10.0, size=n)
    if trait_type == "binomial" and scenario == "null_variable_trial":
        trials = rng.integers(3, 31, size=n)
        payload = {"events": rng.binomial(trials, 0.25), "trials": trials}
    smoke.write_plink(prefix, G)
    start = time.time()
    if trait_type == "ordinal":
        pvals, runtime, conv = run_cosmic_ordinal(prefix, payload, out)
    else:
        pvals, runtime, conv = smoke.run_cosmic(prefix, trait_type, payload, out)
    runtime = max(runtime, time.time() - start)
    converged, warning_flag, status_class = conv
    top = np.argsort(pvals)[:10]
    return {
        "engine": ENGINE,
        "seed": seed,
        "trait_type": trait_type,
        "scenario": scenario,
        "n": n,
        "m": m,
        "rep": rep,
        "runtime": runtime,
        "converged": converged,
        "warning_flag": warning_flag,
        "status_class": status_class,
        "effect_size": effect,
        "case_ratio_or_mean": scenario_case_ratio_or_mean(trait_type, payload),
        "trial_scheme": scenario_trial_scheme(trait_type, payload),
        "ordinal_levels": len(set(payload.get("pheno", []))) if trait_type == "ordinal" else "NA",
        "causal_snp": f"sim{causal+1:06d}",
        "causal_p": float(pvals[causal]),
        "causal_rank": int(np.where(np.argsort(pvals) == causal)[0][0] + 1),
        "causal_top10": int(causal in top),
        "top1_hit": int(int(np.argsort(pvals)[0]) == int(causal)),
        "top10_hit": int(causal in top),
        "window_hit": int(any(abs(int(x) - int(causal)) <= 5 for x in top)),
        "type1_1e_2": float(np.mean(pvals < 1e-2)),
        "type1_1e_3": float(np.mean(pvals < 1e-3)),
        "type1_rate": float(np.mean(pvals < 0.05)) if scenario.startswith("null") else "NA",
        "power": int(causal in top) if scenario.startswith("power") else "NA",
        "lambda_gc": lambda_gc(pvals),
        "min_p": float(np.min(pvals)),
        "median_neglog10p": float(np.median(-np.log10(np.clip(pvals, 1e-300, 1.0)))),
        "pvals": pvals,
    }


def simulate_ordinal_data(rng, n, m, effect):
    freqs = rng.uniform(0.08, 0.45, size=m)
    G = rng.binomial(2, freqs, size=(n, m)).astype(int)
    causal = m // 3
    pop = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
    eta = 0.25 * pop + 0.25 * rng.normal(size=n) + effect * (G[:, causal] - G[:, causal].mean())
    latent = eta + rng.normal(size=n)
    y = np.ones(n, dtype=int)
    y[latent > -0.7] = 2
    y[latent > 0.1] = 3
    y[latent > 0.9] = 4
    if len(set(y.tolist())) < 3:
        y[: n // 3] = 1
        y[n // 3 : 2 * n // 3] = 2
        y[2 * n // 3 :] = 3
    return G, causal, {"pheno": y}


def run_cosmic_ordinal(prefix, payload, out_prefix):
    pheno = smoke.add_ext(prefix, ".ordinal.pheno")
    smoke.write_one_col_pheno(pheno, payload["pheno"])
    cmd = [str(smoke.BIN), "--bfile", str(prefix), "--gwas", "--no-loco", "--maf", "0.0",
           "--psr-maxit", "20", "--ordinal-model", "--ordinal-pheno", str(pheno),
           "--link", "logit", "--out", str(out_prefix)]
    env = os.environ.copy()
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    start = time.time()
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    runtime = time.time() - start
    pvals = []
    with smoke.add_ext(out_prefix, ".gwas.tsv").open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pvals.append(float(row.get("p_value", "1")))
            except ValueError:
                pvals.append(1.0)
    return np.array(pvals), runtime, smoke.convergence_status(out_prefix)


def scenario_case_ratio_or_mean(trait_type, payload):
    if trait_type in ("binary", "ordinal"):
        return float(np.mean(payload["pheno"]))
    if trait_type == "poisson":
        return float(np.mean(payload["pheno"]))
    return float(np.mean(payload["events"] / np.maximum(payload["trials"], 1)))


def scenario_trial_scheme(trait_type, payload):
    if trait_type != "binomial":
        return "NA"
    trials = np.asarray(payload["trials"])
    return "fixed" if np.min(trials) == np.max(trials) else "variable"


def write(path, rows, exclude=()):
    clean = [{k: v for k, v in row.items() if k not in exclude} for row in rows]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(clean[0].keys()))
        writer.writeheader()
        writer.writerows(clean)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    n = int(os.environ.get("PSR_FORMAL_N", "80"))
    m = int(os.environ.get("PSR_FORMAL_M", "80"))
    reps = int(os.environ.get("PSR_FORMAL_REPS", "2"))
    traits = [x.strip() for x in os.environ.get("PSR_FORMAL_TRAITS", "binary,binomial,poisson").split(",") if x.strip()]
    scenarios = {
        "binary": ["null_balanced", "null_unbalanced", "power_low", "power_high"],
        "poisson": ["null_low_mean", "null_medium_mean", "null_high_mean", "power_low", "power_high"],
        "binomial": ["null_fixed_trial", "null_variable_trial", "power_low", "power_high"],
        "ordinal": ["null_balanced", "power_low", "power_high"],
    }
    rows, qq_rows = [], []
    trait_order = ["binary", "binomial", "poisson", "ordinal"]
    for trait_type in traits:
        if trait_type not in scenarios:
            raise SystemExit(f"Unknown PSR_FORMAL_TRAITS entry: {trait_type}")
        scenario_list = scenarios[trait_type]
        for sidx, scenario in enumerate(scenario_list):
            for rep in range(reps):
                seed = BASE_SEED + rep + 100 * sidx + 1000 * trait_order.index(trait_type)
                row = run_one(seed, trait_type, scenario, n, m, rep)
                rows.append(row)
                if scenario.startswith("null"):
                    obs = np.sort(-np.log10(np.clip(row["pvals"], 1e-300, 1.0)))
                    exp = -np.log10((np.arange(1, len(obs) + 1) - 0.5) / len(obs))[::-1]
                    for k in np.linspace(0, len(obs) - 1, 50, dtype=int):
                        qq_rows.append({k2: row[k2] for k2 in ("engine", "seed", "trait_type", "scenario", "n", "m", "rep", "runtime", "converged", "status_class")} | {"qq_index": int(k), "expected_neglog10p": float(exp[k]), "observed_neglog10p": float(obs[k])})
    null_rows = [r for r in rows if r["scenario"].startswith("null")]
    power_rows = [r for r in rows if r["scenario"].startswith("power")]
    write(OUT_DIR / "formal_type1_summary.tsv", null_rows, exclude=("pvals", "causal_p", "causal_rank", "causal_top10"))
    write(OUT_DIR / "formal_power_summary.tsv", power_rows, exclude=("pvals", "type1_1e_2", "type1_1e_3", "lambda_gc"))
    write(OUT_DIR / "formal_runtime_summary.tsv", rows, exclude=("pvals",))
    write(OUT_DIR / "formal_convergence_summary.tsv", rows, exclude=("pvals",))
    write(OUT_DIR / "formal_qq_points.tsv", qq_rows)
    print(f"[PSR-GWAS] Formal simulation outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
