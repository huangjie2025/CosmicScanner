#!/usr/bin/env python3
"""PSR-GWAS v0.5 large-N calibration runner.

Systematically tests whether stacked PSR method conservativeness (lambda GC ~0)
resolves at larger sample sizes. Simulates ordinal data across a grid of sample
sizes, runs all ordinal methods, and produces summary tables + markdown report.

Environment variables:
  PSR_CALIB_REPS       Number of replicates (default: 3)
  PSR_CALIB_N_LIST     Comma-separated list of sample sizes (default: "100,200")
  PSR_CALIB_M          Number of SNPs (default: 200)
  PSR_CALIB_LEVELS     Number of ordinal levels (default: 3)
  PSR_CALIB_METHODS    Comma-separated method list (default: "fixed-gls,dense-zkz-fixed,dense-zkz-grid-reml")
  PSR_CALIB_SCENARIOS  Comma-separated scenario list (default: "null,power")
  PSR_CALIB_OUTDIR     Output directory override (default: tests/results/psr_v05_large_n_calibration)
  PSR_CALIB_SEED       Base seed (default: 20260601)
  PSR_CALIB_EFFECT_SIZE Effect size for power scenario (default: 0.5)
  PSR_CALIB_TESTS    Test statistics to compute: wald|score|robust-score|all (default: wald)
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
CALIB_REPS = int(os.environ.get("PSR_CALIB_REPS", "3"))
RAW_N_LIST = os.environ.get("PSR_CALIB_N_LIST", "100,200")
CALIB_N_LIST = sorted(int(x.strip()) for x in RAW_N_LIST.split(",") if x.strip())
CALIB_M = int(os.environ.get("PSR_CALIB_M", "200"))
CALIB_LEVELS = int(os.environ.get("PSR_CALIB_LEVELS", "3"))
RAW_METHODS = os.environ.get("PSR_CALIB_METHODS", "fixed-gls,dense-zkz-fixed,dense-zkz-grid-reml")
RAW_SCENARIOS = os.environ.get("PSR_CALIB_SCENARIOS", "null,power")
CALIB_SCENARIOS = [x.strip() for x in RAW_SCENARIOS.split(",") if x.strip()]
CALIB_SEED = int(os.environ.get("PSR_CALIB_SEED", "20260601"))
CALIB_EFFECT_SIZE = float(os.environ.get("PSR_CALIB_EFFECT_SIZE", "0.5"))
RAW_TESTS = os.environ.get("PSR_CALIB_TESTS", "wald")
CALIB_TESTS = RAW_TESTS.strip()
OUT_DIR = Path(os.environ.get("PSR_CALIB_OUTDIR", ROOT / "tests" / "results" / "psr_v05_large_n_calibration"))
WORK_DIR = OUT_DIR / "work"

DENSE_LIMIT = 2000

ORDINAL_METHODS_MAP = {
    "cumulative-binary": {
        "mode": "experimental-cumulative-binary",
        "stacked_engine": None,
        "vc_method": None,
        "extra_args": [],
    },
    "fixed-gls": {
        "mode": "stacked-psr",
        "stacked_engine": "fixed-gls",
        "vc_method": "fixed",
        "extra_args": [],
    },
    "dense-zkz-fixed": {
        "mode": "stacked-psr",
        "stacked_engine": "dense-zkz",
        "vc_method": "fixed",
        "extra_args": ["--ordinal-sigma2-g", "0.5"],
    },
    "dense-zkz-grid-reml": {
        "mode": "stacked-psr",
        "stacked_engine": "dense-zkz",
        "vc_method": "grid-reml",
        "extra_args": [
            "--ordinal-vc-grid-min", "0",
            "--ordinal-vc-grid-max", "2",
            "--ordinal-vc-grid-points", "11",
        ],
    },
}


def resolve_methods():
    if RAW_METHODS.strip() == "all":
        return list(ORDINAL_METHODS_MAP.keys())
    return [x.strip() for x in RAW_METHODS.split(",") if x.strip()]


def norm_ppf(p):
    if p <= 0.0 or p >= 1.0:
        return math.nan
    a = [-39.69683028665376, 220.9460984245205, -275.9285104469687,
         138.3577518672690, -30.66479806614716, 2.506628277459239]
    b = [-54.47609879822406, 161.5858368580409, -155.6989798598866,
         66.80131188771972, -13.28068155288572]
    c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
         -2.549732539343734, 4.374664141464968, 2.938163982698783]
    d = [0.007784695709041462, 0.3224671290700398,
         2.445134137142996, 3.754408661907416]
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


def simulate_ordinal_data(rng, n, m, effect, levels):
    freqs = rng.uniform(0.08, 0.45, size=m)
    G = rng.binomial(2, freqs, size=(n, m)).astype(int)
    causal = m // 3
    pop = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
    eta = 0.25 * pop + 0.25 * rng.normal(size=n) + effect * (G[:, causal] - G[:, causal].mean())
    latent = eta + rng.normal(size=n)
    spacing = 1.2
    thresholds = [(-(levels - 2) / 2.0 + k) * spacing for k in range(levels - 1)]
    y = np.ones(n, dtype=int)
    for k, th in enumerate(thresholds):
        y[latent > th] = k + 2
    unique_levels = len(set(y.tolist()))
    if unique_levels < min(3, levels):
        bins = np.array_split(np.arange(n), levels)
        for k, idx in enumerate(bins):
            y[idx] = k + 1
    return G, causal, {"pheno": y}


def run_ordinal_method(prefix, payload, out_prefix, method_name, config, m):
    pheno = smoke.add_ext(prefix, ".ordinal.pheno")
    smoke.write_one_col_pheno(pheno, payload["pheno"])
    cmd = [str(smoke.BIN), "--bfile", str(prefix), "--gwas", "--no-loco",
           "--maf", "0.0", "--psr-maxit", "20", "--ordinal-model",
           "--ordinal-pheno", str(pheno), "--link", "logit",
           "--out", str(out_prefix)]
    if config["mode"] == "experimental-cumulative-binary":
        cmd += ["--ordinal-mode", "experimental-cumulative-binary"]
    else:
        cmd += ["--ordinal-mode", "stacked-psr",
                "--ordinal-stacked-engine", config["stacked_engine"],
                "--ordinal-vc-method", config["vc_method"]]
    cmd += config.get("extra_args", [])
    if CALIB_TESTS != "wald" and config.get("mode") == "stacked-psr":
        cmd += ["--ordinal-test", CALIB_TESTS]
    env = os.environ.copy()
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, env=env)
    pvals = []
    gwas_path = smoke.add_ext(out_prefix, ".gwas.tsv")
    with gwas_path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                pvals.append(float(row.get("p_value", "1")))
            except ValueError:
                pvals.append(1.0)
    pvals_arr = np.array(pvals)
    score_data = {}
    score_path = smoke.add_ext(out_prefix, ".ordinal.score.tsv")
    if os.path.exists(score_path):
        with open(score_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                snp = row.get("snp", "")
                score_data[snp] = {
                    "wald_stat": row.get("wald_stat", "NA"),
                    "wald_p": row.get("wald_p", "NA"),
                    "wald_log10p": row.get("wald_log10p", "NA"),
                    "psr_score_stat": row.get("psr_score_stat", "NA"),
                    "psr_score_p": row.get("psr_score_p", "NA"),
                    "psr_score_log10p": row.get("psr_score_log10p", "NA"),
                    "raw_score_stat": row.get("raw_score_stat", "NA"),
                    "raw_score_p": row.get("raw_score_p", "NA"),
                    "raw_score_log10p": row.get("raw_score_log10p", "NA"),
                    "cluster_score_stat": row.get("cluster_score_stat", "NA"),
                    "cluster_score_p": row.get("cluster_score_p", "NA"),
                    "cluster_score_log10p": row.get("cluster_score_log10p", "NA"),
                    "wald_psr_score_equivalent": row.get("wald_psr_score_equivalent", "NA"),
                    "raw_score_projection": row.get("raw_score_projection", "NA"),
                    "n_clusters": row.get("n_clusters", "NA"),
                }
    conv = smoke.convergence_status(out_prefix)
    sigma2_g = "NA"
    sigma2_g_boundary = "NA"
    reml_loglik = "NA"
    n_grid_valid = "NA"
    vc_path = smoke.add_ext(out_prefix, ".ordinal.vc.tsv")
    if vc_path.exists():
        with vc_path.open("r", encoding="utf-8") as f:
            vc_rows = list(csv.DictReader(f, delimiter="\t"))
            if vc_rows:
                last = vc_rows[-1]
                sg_raw = last.get("sigma2_g", "NA")
                try:
                    sigma2_g = float(sg_raw) if sg_raw not in ("", "NA", "nan") else "NA"
                except (ValueError, TypeError):
                    sigma2_g = "NA"
                bv = last.get("boundary", "").strip().upper()
                if bv in ("TRUE", "1", "YES"):
                    sigma2_g_boundary = 1
                elif bv in ("FALSE", "0", "NO", ""):
                    sigma2_g_boundary = 0
                ll_raw = last.get("reml_loglik", "").strip()
                if ll_raw and ll_raw not in ("", "NA", "nan"):
                    try:
                        reml_loglik = float(ll_raw)
                    except ValueError:
                        pass
                ngv_raw = last.get("n_grid_valid", "").strip()
                if ngv_raw and ngv_raw not in ("", "NA", "nan"):
                    try:
                        n_grid_valid = int(float(ngv_raw))
                    except ValueError:
                        pass
    has_blup = 1 if smoke.add_ext(out_prefix, ".ordinal.blup.tsv").exists() else 0
    n_stack = "NA"
    conv_path = smoke.add_ext(out_prefix, ".psr.convergence.tsv")
    if conv_path.exists():
        with conv_path.open("r", encoding="utf-8") as f:
            conv_rows = list(csv.DictReader(f, delimiter="\t"))
            if conv_rows:
                ns_raw = conv_rows[-1].get("n_stack", "").strip()
                if ns_raw and ns_raw not in ("", "NA", "nan"):
                    try:
                        n_stack = int(float(ns_raw))
                    except ValueError:
                        pass
    ordinal_extra = {
        "ordinal_method": method_name,
        "ordinal_stacked_engine": config.get("stacked_engine", "NA") or "NA",
        "ordinal_vc_method": config.get("vc_method", "NA") or "NA",
        "sigma2_g": sigma2_g,
        "sigma2_g_boundary": sigma2_g_boundary,
        "has_blup": has_blup,
        "n_stack": n_stack,
        "reml_loglik": reml_loglik,
        "n_grid_valid": n_grid_valid,
        "vc_update": 1 if config.get("vc_method") == "grid-reml" else 0,
    }
    return pvals_arr, conv, ordinal_extra, score_data


def build_row(seed, scenario, n, m, levels, rep, effect, causal, pvals, runtime, conv, extra, method_status, error_message="", score_data=None):
    converged, warning_flag, status_class = conv
    if len(pvals) < m:
        padded = np.ones(m)
        padded[:len(pvals)] = pvals[:]
        pvals = padded
    top10 = np.argsort(pvals)[:10]
    causal_present = causal is not None and causal < len(pvals)
    top1_hit = int(int(np.argsort(pvals)[0]) == int(causal)) if causal_present else "NA"
    top10_hit = 1 if causal_present and causal in top10 else (0 if causal_present else "NA")
    window_hit = int(any(abs(int(x) - int(causal)) <= 5 for x in top10)) if causal_present else "NA"
    max_log10p = float(-np.log10(np.clip(np.min(pvals), 1e-300, 1.0)))
    mean_log10p = float(np.mean(-np.log10(np.clip(pvals, 1e-300, 1.0))))
    row = {
        "rep": rep,
        "seed": seed,
        "n": n,
        "m": m,
        "levels": levels,
        "n_stack": extra.get("n_stack", "NA"),
        "scenario": scenario,
        "method": extra.get("ordinal_method", "NA"),
        "ordinal_mode": extra.get("ordinal_mode", "NA"),
        "ordinal_implementation": extra.get("ordinal_implementation", "NA"),
        "ordinal_stacked_engine": extra.get("ordinal_stacked_engine", "NA"),
        "ordinal_vc_method": extra.get("ordinal_vc_method", "NA"),
        "runtime_sec": runtime,
        "converged": converged,
        "status_class": status_class,
        "method_status": method_status,
        "lambda_gc": lambda_gc(pvals),
        "top1_hit": top1_hit,
        "top10_hit": top10_hit,
        "window_hit": window_hit,
        "valid_snp_count": len(pvals),
        "max_log10p": max_log10p,
        "mean_log10p": mean_log10p,
        "sigma2_g": extra.get("sigma2_g", "NA"),
        "sigma2_g_boundary": extra.get("sigma2_g_boundary", "NA"),
        "has_blup": extra.get("has_blup", 0),
        "error_message": error_message,
    }
    if score_data:
        psr_p = []
        raw_p = []
        cluster_p = []
        wald_psr_eq = "NA"
        raw_proj = "NA"
        for snp, sd in score_data.items():
            if sd.get("psr_score_p", "NA") != "NA":
                try: psr_p.append(float(sd["psr_score_p"]))
                except ValueError: pass
            if sd.get("raw_score_p", "NA") != "NA":
                try: raw_p.append(float(sd["raw_score_p"]))
                except ValueError: pass
            if sd.get("cluster_score_p", "NA") != "NA":
                try: cluster_p.append(float(sd["cluster_score_p"]))
                except ValueError: pass
            if wald_psr_eq == "NA" and sd.get("wald_psr_score_equivalent", "NA") != "NA":
                wald_psr_eq = sd["wald_psr_score_equivalent"]
            if raw_proj == "NA" and sd.get("raw_score_projection", "NA") != "NA":
                raw_proj = sd["raw_score_projection"]
        row["psr_score_lambda_gc"] = lambda_gc(psr_p) if psr_p else "NA"
        row["raw_score_lambda_gc"] = lambda_gc(raw_p) if raw_p else "NA"
        row["cluster_score_lambda_gc"] = lambda_gc(cluster_p) if cluster_p else "NA"
        row["wald_psr_score_equivalent"] = wald_psr_eq
        row["raw_score_projection"] = raw_proj
        psr_arr = np.array(psr_p) if psr_p else np.array([])
        raw_arr = np.array(raw_p) if raw_p else np.array([])
        cluster_arr = np.array(cluster_p) if cluster_p else np.array([])
        for name, spvals in [("psr_score", psr_arr), ("raw_score", raw_arr), ("cluster_score", cluster_arr)]:
            if len(spvals) == 0:
                row[f"{name}_top1_hit"] = "NA"
                row[f"{name}_top10_hit"] = "NA"
                row[f"{name}_window_hit"] = "NA"
                row[f"{name}_max_log10p"] = "NA"
            else:
                causal_present = causal is not None and causal < len(spvals)
                stop10 = np.argsort(spvals)[:10]
                row[f"{name}_top1_hit"] = int(int(np.argsort(spvals)[0]) == int(causal)) if causal_present else "NA"
                row[f"{name}_top10_hit"] = 1 if causal_present and causal in stop10 else (0 if causal_present else "NA")
                row[f"{name}_window_hit"] = int(any(abs(int(x) - int(causal)) <= 5 for x in stop10)) if causal_present else "NA"
                row[f"{name}_max_log10p"] = float(-np.log10(np.clip(np.min(spvals), 1e-300, 1.0)))
    return row


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    fieldnames = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compute_lambda_summary(rows):
    groups = {}
    for r in rows:
        key = (r["n"], r["m"], r["levels"], r["scenario"], r["method"])
        groups.setdefault(key, []).append(r)
    out = []
    for (n, m, levels, scenario, method), grp in groups.items():
        lams = [r["lambda_gc"] for r in grp if isinstance(r["lambda_gc"], (int, float))]
        convs = [r["converged"] for r in grp if isinstance(r["converged"], (int, float))]
        n_reps = len(grp)
        if lams:
            lam_arr = np.array(lams)
            near1 = float(np.mean((lam_arr >= 0.8) & (lam_arr <= 1.2)))
            lt05 = float(np.mean(lam_arr < 0.5))
        else:
            near1 = "NA"
            lt05 = "NA"
        out.append({
            "n": n,
            "m": m,
            "levels": levels,
            "scenario": scenario,
            "method": method,
            "n_reps": n_reps,
            "mean_lambda_gc": float(np.mean(lams)) if lams else "NA",
            "sd_lambda_gc": float(np.std(lams, ddof=1)) if len(lams) > 1 else "NA",
            "min_lambda_gc": float(np.min(lams)) if lams else "NA",
            "max_lambda_gc": float(np.max(lams)) if lams else "NA",
            "lambda_lt_0_5_rate": lt05,
            "lambda_near_1_rate": near1,
            "convergence_rate": float(np.mean(convs)) if convs else "NA",
        })
    return out


def compute_power_summary(rows):
    power_rows = [r for r in rows if r["scenario"] == "power"]
    if not power_rows:
        return []
    groups = {}
    for r in power_rows:
        key = (r["n"], r["m"], r["levels"], r["method"])
        groups.setdefault(key, []).append(r)
    out = []
    for (n, m, levels, method), grp in groups.items():
        t1 = [r["top1_hit"] for r in grp if isinstance(r["top1_hit"], int)]
        t10 = [r["top10_hit"] for r in grp if isinstance(r["top10_hit"], int)]
        wh = [r["window_hit"] for r in grp if isinstance(r["window_hit"], int)]
        mx = [r["max_log10p"] for r in grp if isinstance(r["max_log10p"], (int, float))]
        out.append({
            "n": n,
            "m": m,
            "levels": levels,
            "method": method,
            "n_reps": len(grp),
            "top1_rate": float(np.mean(t1)) if t1 else "NA",
            "top10_rate": float(np.mean(t10)) if t10 else "NA",
            "window_hit_rate": float(np.mean(wh)) if wh else "NA",
            "mean_max_log10p": float(np.mean(mx)) if mx else "NA",
        })
    return out


def compute_score_lambda_summary(rows):
    test_types = ["psr_score", "raw_score", "cluster_score"]
    groups = {}
    for r in rows:
        for tt in test_types:
            val = r.get(f"{tt}_lambda_gc")
            if val not in ("NA", None):
                key = (r["n"], r["m"], r["levels"], r["scenario"], r["method"], tt)
                groups.setdefault(key, []).append(r)
    out = []
    for (n, m, levels, scenario, method, test_type), grp in groups.items():
        lam_key = f"{test_type}_lambda_gc"
        lams = [r[lam_key] for r in grp if isinstance(r.get(lam_key), (int, float))]
        convs = [r["converged"] for r in grp if isinstance(r["converged"], (int, float))]
        n_reps = len(grp)
        if lams:
            lam_arr = np.array(lams)
            near1 = float(np.mean((lam_arr >= 0.8) & (lam_arr <= 1.2)))
            lt05 = float(np.mean(lam_arr < 0.5))
        else:
            near1 = "NA"
            lt05 = "NA"
        out.append({
            "n": n,
            "m": m,
            "levels": levels,
            "scenario": scenario,
            "method": method,
            "test_type": test_type,
            "n_reps": n_reps,
            "mean_lambda_gc": float(np.mean(lams)) if lams else "NA",
            "sd_lambda_gc": float(np.std(lams, ddof=1)) if len(lams) > 1 else "NA",
            "min_lambda_gc": float(np.min(lams)) if lams else "NA",
            "max_lambda_gc": float(np.max(lams)) if lams else "NA",
            "lambda_lt_0_5_rate": lt05,
            "lambda_near_1_rate": near1,
            "convergence_rate": float(np.mean(convs)) if convs else "NA",
        })
    return out


def compute_score_power_summary(rows):
    power_rows = [r for r in rows if r["scenario"] == "power"]
    if not power_rows:
        return []
    test_types = ["psr_score", "raw_score", "cluster_score"]
    groups = {}
    for r in power_rows:
        for tt in test_types:
            val = r.get(f"{tt}_max_log10p")
            if val not in ("NA", None):
                key = (r["n"], r["m"], r["levels"], r["method"], tt)
                groups.setdefault(key, []).append(r)
    out = []
    for (n, m, levels, method, test_type), grp in groups.items():
        t1 = [r[f"{test_type}_top1_hit"] for r in grp if isinstance(r.get(f"{test_type}_top1_hit"), int)]
        t10 = [r[f"{test_type}_top10_hit"] for r in grp if isinstance(r.get(f"{test_type}_top10_hit"), int)]
        wh = [r[f"{test_type}_window_hit"] for r in grp if isinstance(r.get(f"{test_type}_window_hit"), int)]
        mx = [r[f"{test_type}_max_log10p"] for r in grp if isinstance(r.get(f"{test_type}_max_log10p"), (int, float))]
        out.append({
            "n": n,
            "m": m,
            "levels": levels,
            "method": method,
            "test_type": test_type,
            "n_reps": len(grp),
            "top1_rate": float(np.mean(t1)) if t1 else "NA",
            "top10_rate": float(np.mean(t10)) if t10 else "NA",
            "window_hit_rate": float(np.mean(wh)) if wh else "NA",
            "mean_max_log10p": float(np.mean(mx)) if mx else "NA",
        })
    return out


def compute_runtime_summary(rows):
    groups = {}
    for r in rows:
        key = (r["n"], r["method"])
        groups.setdefault(key, []).append(r)
    out = []
    for (n, method), grp in groups.items():
        rt = [r["runtime_sec"] for r in grp if isinstance(r["runtime_sec"], (int, float))]
        out.append({
            "n": n,
            "method": method,
            "n_reps": len(grp),
            "mean_runtime": float(np.mean(rt)) if rt else "NA",
            "sd_runtime": float(np.std(rt, ddof=1)) if len(rt) > 1 else "NA",
            "min_runtime": float(np.min(rt)) if rt else "NA",
            "max_runtime": float(np.max(rt)) if rt else "NA",
        })
    return out


def compute_sigma2_summary(rows):
    sigma2_rows = [r for r in rows if r.get("sigma2_g") not in ("NA", None) and isinstance(r.get("sigma2_g"), (int, float))]
    if not sigma2_rows:
        return []
    groups = {}
    for r in sigma2_rows:
        key = (r["n"], r["method"])
        groups.setdefault(key, []).append(r)
    out = []
    for (n, method), grp in groups.items():
        sg = [r["sigma2_g"] for r in grp if isinstance(r["sigma2_g"], (int, float))]
        bd = [r["sigma2_g_boundary"] for r in grp if r["sigma2_g_boundary"] not in ("NA", None) and isinstance(r["sigma2_g_boundary"], (int, float))]
        out.append({
            "n": n,
            "method": method,
            "n_reps": len(grp),
            "mean_sigma2_g": float(np.mean(sg)) if sg else "NA",
            "sd_sigma2_g": float(np.std(sg, ddof=1)) if len(sg) > 1 else "NA",
            "min_sigma2_g": float(np.min(sg)) if sg else "NA",
            "max_sigma2_g": float(np.max(sg)) if sg else "NA",
            "boundary_rate": float(np.mean(bd)) if bd else "NA",
        })
    return out


def compute_method_status_summary(rows):
    groups = {}
    for r in rows:
        key = (r["n"], r["method"])
        groups.setdefault(key, []).append(r)
    out = []
    for (n, method), grp in groups.items():
        statuses = [r.get("method_status", "UNKNOWN") for r in grp]
        errors = [r["error_message"] for r in grp if r["error_message"]]
        out.append({
            "n": n,
            "method": method,
            "n_reps": len(grp),
            "status_OK": sum(1 for s in statuses if s == "OK"),
            "status_CONVERGED_WITH_WARNING": sum(1 for s in statuses if s == "CONVERGED_WITH_WARNING"),
            "status_NOT_CONVERGED": sum(1 for s in statuses if s == "NOT_CONVERGED"),
            "status_unsupported_dense_limit": sum(1 for s in statuses if s == "unsupported_dense_limit"),
            "error_count": len(errors),
        })
    return out


def generate_report(rows, lam_summary, power_summary, runtime_summary, sigma2_summary, method_status_summary):
    lines = [
        "# PSR-GWAS v0.5 Large-N Calibration Report",
        f"Date: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        "",
        "## Configuration",
        f"- Replicates: {CALIB_REPS}",
        f"- Sample sizes (N): {', '.join(str(x) for x in CALIB_N_LIST)}",
        f"- SNPs (M): {CALIB_M}",
        f"- Ordinal levels: {CALIB_LEVELS}",
        f"- Base seed: {CALIB_SEED}",
    ]
    methods = sorted(set(r["method"] for r in rows))
    lines.append(f"- Methods: {', '.join(methods)}")
    scenarios = sorted(set(r["scenario"] for r in rows))
    lines.append(f"- Scenarios: {', '.join(scenarios)}")
    if "power" in scenarios:
        lines.append(f"- Power effect size: {CALIB_EFFECT_SIZE}")
    lines.append("")

    verdict_key = {
        "persistent_conservative_bias": "lambda_GC < 0.5",
        "moderately_conservative": "0.5 <= lambda_GC < 0.8",
        "acceptable": "0.8 <= lambda_GC <= 1.2",
        "inflated": "lambda_GC > 1.3",
    }

    def verdict(mean_lam, n_reps):
        if n_reps < 2:
            return "needs_more_reps"
        if not isinstance(mean_lam, (int, float)):
            return "needs_more_reps"
        if mean_lam < 0.5:
            return "persistent_conservative_bias"
        if mean_lam < 0.8:
            return "moderately_conservative"
        if mean_lam <= 1.2:
            return "acceptable"
        if mean_lam > 1.3:
            return "inflated"
        return "acceptable"

    lines.append("## Lambda GC Summary (Null Scenarios)")
    lines.append("")
    null_lam = [r for r in lam_summary if r["scenario"] == "null"]
    if null_lam:
        hdr = ["N", "Method", "Reps", "Mean Lambda", "SD", "Min", "Max", "lt0.5 Rate", "Near1 Rate", "Convergence", "Verdict"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(null_lam, key=lambda x: (x["n"], x["method"])):
            mean_v = f"{sr['mean_lambda_gc']:.4f}" if isinstance(sr['mean_lambda_gc'], float) else "NA"
            sd_v = f"{sr['sd_lambda_gc']:.4f}" if isinstance(sr['sd_lambda_gc'], float) else "NA"
            lo_v = f"{sr['min_lambda_gc']:.4f}" if isinstance(sr['min_lambda_gc'], float) else "NA"
            hi_v = f"{sr['max_lambda_gc']:.4f}" if isinstance(sr['max_lambda_gc'], float) else "NA"
            lt05 = f"{sr['lambda_lt_0_5_rate'] * 100:.0f}%" if isinstance(sr['lambda_lt_0_5_rate'], float) else "NA"
            nr1 = f"{sr['lambda_near_1_rate'] * 100:.0f}%" if isinstance(sr['lambda_near_1_rate'], float) else "NA"
            conv = f"{sr['convergence_rate'] * 100:.0f}%" if isinstance(sr['convergence_rate'], float) else "NA"
            vd = verdict(sr.get("mean_lambda_gc"), sr["n_reps"])
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['n_reps']} | {mean_v} | {sd_v} | {lo_v} | {hi_v} | {lt05} | {nr1} | {conv} | {vd} |")
    else:
        lines.append("(No null-scenario lambda data)")
    lines.append("")

    power_rows = [r for r in rows if r["scenario"] == "power"]
    if power_summary:
        lines.append("## Power Summary (Power Scenarios)")
        lines.append("")
        hdr = ["N", "Method", "Reps", "Top1 Rate", "Top10 Rate", "Window Hit Rate", "Mean Max -log10(p)"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(power_summary, key=lambda x: (x["n"], x["method"])):
            t1 = f"{sr['top1_rate'] * 100:.0f}%" if isinstance(sr['top1_rate'], float) else "NA"
            t10 = f"{sr['top10_rate'] * 100:.0f}%" if isinstance(sr['top10_rate'], float) else "NA"
            wh = f"{sr['window_hit_rate'] * 100:.0f}%" if isinstance(sr['window_hit_rate'], float) else "NA"
            mx = f"{sr['mean_max_log10p']:.2f}" if isinstance(sr['mean_max_log10p'], float) else "NA"
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['n_reps']} | {t1} | {t10} | {wh} | {mx} |")
        lines.append("")

    if sigma2_summary:
        lines.append("## Sigma2_g Summary (Methods that estimate it)")
        lines.append("")
        hdr = ["N", "Method", "Reps", "Mean", "SD", "Min", "Max", "Boundary Rate"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(sigma2_summary, key=lambda x: (x["n"], x["method"])):
            mu = f"{sr['mean_sigma2_g']:.4f}" if isinstance(sr['mean_sigma2_g'], float) else "NA"
            sd = f"{sr['sd_sigma2_g']:.4f}" if isinstance(sr['sd_sigma2_g'], float) else "NA"
            lo = f"{sr['min_sigma2_g']:.4f}" if isinstance(sr['min_sigma2_g'], float) else "NA"
            hi = f"{sr['max_sigma2_g']:.4f}" if isinstance(sr['max_sigma2_g'], float) else "NA"
            br = f"{sr['boundary_rate'] * 100:.0f}%" if isinstance(sr['boundary_rate'], float) else "NA"
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['n_reps']} | {mu} | {sd} | {lo} | {hi} | {br} |")
        lines.append("")

    lines.append("## Method Status")
    lines.append("")
    if method_status_summary:
        hdr = ["N", "Method", "Reps", "OK", "Warnings", "Not Converged", "Dense-Skip", "Errors"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(method_status_summary, key=lambda x: (x["n"], x["method"])):
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['n_reps']} | {sr['status_OK']} | {sr['status_CONVERGED_WITH_WARNING']} | {sr['status_NOT_CONVERGED']} | {sr['status_unsupported_dense_limit']} | {sr['error_count']} |")
        lines.append("")

    lines.append("## Runtime Summary")
    lines.append("")
    if runtime_summary:
        hdr = ["N", "Method", "Reps", "Mean (s)", "SD", "Min", "Max"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(runtime_summary, key=lambda x: (x["n"], x["method"])):
            mn = f"{sr['mean_runtime']:.3f}" if isinstance(sr['mean_runtime'], float) else "NA"
            sd = f"{sr['sd_runtime']:.3f}" if isinstance(sr['sd_runtime'], float) else "NA"
            lo = f"{sr['min_runtime']:.3f}" if isinstance(sr['min_runtime'], float) else "NA"
            hi = f"{sr['max_runtime']:.3f}" if isinstance(sr['max_runtime'], float) else "NA"
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['n_reps']} | {mn} | {sd} | {lo} | {hi} |")
        lines.append("")

    lines.append("## Key Findings")
    lines.append("")
    has_finding = False
    for sr in null_lam or []:
        v = verdict(sr.get("mean_lambda_gc"), sr["n_reps"])
        if v == "persistent_conservative_bias":
            lines.append(f"- N={sr['n']}, {sr['method']}: persistent conservative bias (lambda_GC={sr['mean_lambda_gc']:.3f}). Conservativeness does NOT resolve at this sample size.")
            has_finding = True
        elif v == "moderately_conservative":
            lines.append(f"- N={sr['n']}, {sr['method']}: moderately conservative (lambda_GC={sr['mean_lambda_gc']:.3f}). Increasing N may resolve.")
            has_finding = True
        elif v == "acceptable":
            lines.append(f"- N={sr['n']}, {sr['method']}: acceptable calibration (lambda_GC={sr['mean_lambda_gc']:.3f}). Conservativeness resolved at this N.")
            has_finding = True
        elif v == "inflated":
            lines.append(f"- N={sr['n']}, {sr['method']}: inflated (lambda_GC={sr['mean_lambda_gc']:.3f}). Anti-conservative behavior detected.")
            has_finding = True
    if not has_finding:
        lines.append("- No definitive findings (insufficient data or all needs_more_reps).")
    lines.append("")
    lines.append("## Output Files")
    for f in [
        "large_n_calibration_runs.tsv",
        "large_n_lambda_summary.tsv",
        "large_n_power_summary.tsv",
        "large_n_runtime_summary.tsv",
        "large_n_sigma2_summary.tsv",
        "large_n_method_status.tsv",
        "large_n_calibration_report.md",
    ]:
        lines.append(f"- `{f}`")
    lines.append("")
    lines.append("*Report auto-generated by `run_psr_v05_large_n_calibration.py`.*")
    return "\n".join(lines)


def generate_score_diagnostic_report(rows, lam_summary, power_summary):
    lines = [
        "# PSR-GWAS v0.5 Large-N Score Test Diagnostic Report",
        f"Date: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        "",
        "## Configuration",
        f"- Test types requested: {CALIB_TESTS}",
        f"- Replicates: {CALIB_REPS}",
        f"- Sample sizes (N): {', '.join(str(x) for x in CALIB_N_LIST)}",
        f"- SNPs (M): {CALIB_M}",
        "",
        "## Lambda GC Comparison (Null Scenarios)",
        "",
    ]
    null_lam = [r for r in lam_summary if r["scenario"] == "null"]
    if null_lam:
        hdr = ["N", "Method", "Test Type", "Reps", "Mean Lambda", "SD", "Min", "Max", "lt0.5 Rate", "Near1 Rate", "Convergence"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(null_lam, key=lambda x: (x["n"], x["method"], x["test_type"])):
            mean_v = f"{sr['mean_lambda_gc']:.4f}" if isinstance(sr['mean_lambda_gc'], float) else "NA"
            sd_v = f"{sr['sd_lambda_gc']:.4f}" if isinstance(sr['sd_lambda_gc'], float) else "NA"
            lo_v = f"{sr['min_lambda_gc']:.4f}" if isinstance(sr['min_lambda_gc'], float) else "NA"
            hi_v = f"{sr['max_lambda_gc']:.4f}" if isinstance(sr['max_lambda_gc'], float) else "NA"
            lt05 = f"{sr['lambda_lt_0_5_rate'] * 100:.0f}%" if isinstance(sr['lambda_lt_0_5_rate'], float) else "NA"
            nr1 = f"{sr['lambda_near_1_rate'] * 100:.0f}%" if isinstance(sr['lambda_near_1_rate'], float) else "NA"
            conv = f"{sr['convergence_rate'] * 100:.0f}%" if isinstance(sr['convergence_rate'], float) else "NA"
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['test_type']} | {sr['n_reps']} | {mean_v} | {sd_v} | {lo_v} | {hi_v} | {lt05} | {nr1} | {conv} |")
    else:
        lines.append("(No null-scenario score lambda data)")
    lines.append("")

    if power_summary:
        lines.append("## Power Comparison (Power Scenarios)")
        lines.append("")
        hdr = ["N", "Method", "Test Type", "Reps", "Top1 Rate", "Top10 Rate", "Window Hit Rate", "Mean Max -log10(p)"]
        lines.append("| " + " | ".join(hdr) + " |")
        lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
        for sr in sorted(power_summary, key=lambda x: (x["n"], x["method"], x["test_type"])):
            t1 = f"{sr['top1_rate'] * 100:.0f}%" if isinstance(sr['top1_rate'], float) else "NA"
            t10 = f"{sr['top10_rate'] * 100:.0f}%" if isinstance(sr['top10_rate'], float) else "NA"
            wh = f"{sr['window_hit_rate'] * 100:.0f}%" if isinstance(sr['window_hit_rate'], float) else "NA"
            mx = f"{sr['mean_max_log10p']:.2f}" if isinstance(sr['mean_max_log10p'], float) else "NA"
            lines.append(f"| {sr['n']} | {sr['method']} | {sr['test_type']} | {sr['n_reps']} | {t1} | {t10} | {wh} | {mx} |")
        lines.append("")

    lines.append("## Key Findings")
    lines.append("")
    has_finding = False
    for sr in null_lam or []:
        lam = sr.get("mean_lambda_gc")
        if isinstance(lam, float):
            lines.append(f"- N={sr['n']}, {sr['method']}, {sr['test_type']}: lambda_GC = {lam:.4f}")
            has_finding = True
    if not has_finding:
        lines.append("- No score test findings (insufficient data).")
    lines.append("")

    lines.append("## Output Files")
    for fname in [
        "large_n_score_calibration_runs.tsv",
        "large_n_score_lambda_summary.tsv",
        "large_n_score_power_summary.tsv",
        "large_n_score_diagnostic_report.md",
    ]:
        lines.append(f"- `{fname}`")
    lines.append("")
    lines.append("*Report auto-generated by `run_psr_v05_large_n_calibration.py`.*")
    return "\n".join(lines)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    methods = resolve_methods()

    scenario_list = []
    if "null" in CALIB_SCENARIOS:
        scenario_list.append(("null", 0.0))
    if "power" in CALIB_SCENARIOS:
        scenario_list.append(("power", CALIB_EFFECT_SIZE))

    all_rows = []

    for n_val in CALIB_N_LIST:
        for rep in range(CALIB_REPS):
            for sidx, (scenario, effect) in enumerate(scenario_list):
                seed = CALIB_SEED + rep + 100 * sidx + 1000 * CALIB_N_LIST.index(n_val)
                rng = np.random.default_rng(seed)
                prefix = WORK_DIR / f"large_n.n{n_val}.{scenario}.rep{rep}"
                G, causal, payload = simulate_ordinal_data(rng, n_val, CALIB_M, effect, CALIB_LEVELS)
                smoke.write_plink(prefix, G)
                for method in methods:
                    config = ORDINAL_METHODS_MAP[method]
                    out_method = WORK_DIR / f"large_n.n{n_val}.{scenario}.rep{rep}.cosmic.{method}"
                    method_status = "OK"
                    method_error = ""

                    is_dense = config.get("stacked_engine") == "dense-zkz"
                    n_stack_check = n_val * (CALIB_LEVELS - 1)
                    if is_dense and n_stack_check > DENSE_LIMIT:
                        method_status = "unsupported_dense_limit"
                        pvals = np.full(CALIB_M, 1.0)
                        conv = ("NA", "NA", "NA")
                        extra = {
                            "ordinal_method": method,
                            "ordinal_stacked_engine": config.get("stacked_engine", "NA") or "NA",
                            "ordinal_vc_method": config.get("vc_method", "NA") or "NA",
                            "sigma2_g": "NA", "sigma2_g_boundary": "NA",
                            "has_blup": 0, "n_stack": "NA",
                            "reml_loglik": "NA", "n_grid_valid": "NA",
                            "vc_update": 1 if config.get("vc_method") == "grid-reml" else 0,
                        }
                        row = build_row(seed, scenario, n_val, CALIB_M, CALIB_LEVELS, rep, effect,
                                        causal, pvals, 0.0, conv, extra, method_status, method_error, score_data={})
                        row["ordinal_mode"] = config["mode"]
                        row["ordinal_implementation"] = f"stacked-psr-{config['stacked_engine']}"
                        if config.get("vc_method") == "grid-reml":
                            row["ordinal_implementation"] += "-vc-grid-reml"
                        all_rows.append(row)
                        continue

                    _start = time.time()
                    try:
                        pvals, conv, extra, score_data = run_ordinal_method(prefix, payload, out_method, method, config, CALIB_M)
                    except Exception as e:
                        pvals = np.full(CALIB_M, 1.0)
                        conv = (0, "ERROR", "ERROR")
                        extra = {
                            "ordinal_method": method,
                            "ordinal_stacked_engine": config.get("stacked_engine", "NA") or "NA",
                            "ordinal_vc_method": config.get("vc_method", "NA") or "NA",
                            "sigma2_g": "NA", "sigma2_g_boundary": "NA",
                            "has_blup": 0, "n_stack": "NA",
                            "reml_loglik": "NA", "n_grid_valid": "NA",
                            "vc_update": 1 if config.get("vc_method") == "grid-reml" else 0,
                        }
                        score_data = {}
                        method_error = str(e)
                        method_status = "ERROR"
                    runtime = time.time() - _start

                    if method_status != "ERROR" and isinstance(conv, tuple) and len(conv) == 3:
                        _, _, status_cls = conv
                        method_status = status_cls

                    row = build_row(seed, scenario, n_val, CALIB_M, CALIB_LEVELS, rep, effect,
                                    causal, pvals, runtime, conv, extra, method_status, method_error, score_data=score_data)
                    row["ordinal_mode"] = config["mode"]
                    row["ordinal_implementation"] = f"stacked-psr-{config['stacked_engine']}"
                    if config.get("vc_method") == "grid-reml":
                        row["ordinal_implementation"] += "-vc-grid-reml"
                    if config["mode"] == "experimental-cumulative-binary":
                        row["ordinal_implementation"] = "experimental-cumulative-binary"
                    all_rows.append(row)

    write_tsv(OUT_DIR / "large_n_calibration_runs.tsv", all_rows)

    lam_summary = compute_lambda_summary(all_rows)
    write_tsv(OUT_DIR / "large_n_lambda_summary.tsv", lam_summary)

    power_summary = compute_power_summary(all_rows)
    write_tsv(OUT_DIR / "large_n_power_summary.tsv", power_summary)

    runtime_summary = compute_runtime_summary(all_rows)
    write_tsv(OUT_DIR / "large_n_runtime_summary.tsv", runtime_summary)

    sigma2_summary = compute_sigma2_summary(all_rows)
    write_tsv(OUT_DIR / "large_n_sigma2_summary.tsv", sigma2_summary)

    method_status_summary = compute_method_status_summary(all_rows)
    write_tsv(OUT_DIR / "large_n_method_status.tsv", method_status_summary)

    report = generate_report(all_rows, lam_summary, power_summary, runtime_summary, sigma2_summary, method_status_summary)
    report_path = OUT_DIR / "large_n_calibration_report.md"
    report_path.write_text(report + "\n")
    print(f"[PSR-GWAS v0.5] Large-N calibration outputs written to {OUT_DIR}")
    print(f"[PSR-GWAS v0.5] Report: {report_path}")

    if CALIB_TESTS != "wald":
        score_rows = [r for r in all_rows if r.get("psr_score_lambda_gc") not in ("NA", None)]
        write_tsv(OUT_DIR / "large_n_score_calibration_runs.tsv", score_rows)

        score_lam_summary = compute_score_lambda_summary(score_rows)
        write_tsv(OUT_DIR / "large_n_score_lambda_summary.tsv", score_lam_summary)

        score_power = compute_score_power_summary(score_rows)
        write_tsv(OUT_DIR / "large_n_score_power_summary.tsv", score_power)

        diag_report = generate_score_diagnostic_report(score_rows, score_lam_summary, score_power)
        diag_path = OUT_DIR / "large_n_score_diagnostic_report.md"
        diag_path.write_text(diag_report + "\n")
        print(f"[PSR-GWAS v0.5] Score diagnostic report: {diag_path}")


if __name__ == "__main__":
    main()
