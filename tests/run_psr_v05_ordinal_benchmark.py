#!/usr/bin/env python3
"""PSR-GWAS v0.5 ordinal benchmark runner.

Benchmarks all ordinal methods on simulated ordinal data and outputs summary
tables and a markdown report.

Environment variables:
  PSR_V05_BENCH_REPS      Number of replicates (default: 3)
  PSR_V05_BENCH_N         Sample size (default: 80)
  PSR_V05_BENCH_M         Number of SNPs (default: 100)
  PSR_V05_BENCH_LEVELS    Number of ordinal levels (default: 3)
  PSR_V05_BENCH_METHODS   Comma-separated list or "all" (default: "all")
  PSR_V05_BENCH_OUTDIR    Output directory override (default: tests/results/psr_v05_ordinal_benchmark)
  PSR_V05_BENCH_SEED      Base seed (default: 20260601)
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
BENCH_REPS = int(os.environ.get("PSR_V05_BENCH_REPS", "3"))
BENCH_N = int(os.environ.get("PSR_V05_BENCH_N", "80"))
BENCH_M = int(os.environ.get("PSR_V05_BENCH_M", "100"))
BENCH_LEVELS = int(os.environ.get("PSR_V05_BENCH_LEVELS", "3"))
BENCH_SEED = int(os.environ.get("PSR_V05_BENCH_SEED", "20260601"))
RAW_METHODS = os.environ.get("PSR_V05_BENCH_METHODS", "all")
OUT_DIR = Path(
    os.environ.get(
        "PSR_V05_BENCH_OUTDIR",
        ROOT / "tests" / "results" / "psr_v05_ordinal_benchmark",
    )
)
WORK_DIR = OUT_DIR / "work"

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


def scenario_effect(scenario):
    return 0.0 if scenario.startswith("null") else 0.35 if scenario.endswith("low") else 0.7


def simulate_ordinal_data(rng, n, m, effect, levels=3):
    freqs = rng.uniform(0.08, 0.45, size=m)
    G = rng.binomial(2, freqs, size=(n, m)).astype(int)
    causal = m // 3
    pop = np.r_[np.zeros(n // 2), np.ones(n - n // 2)]
    eta = 0.25 * pop + 0.25 * rng.normal(size=n) + effect * (G[:, causal] - G[:, causal].mean())
    latent = eta + rng.normal(size=n)
    if levels == 4:
        thresh = [-0.7, 0.1, 0.9]
        y = np.ones(n, dtype=int)
        y[latent > thresh[0]] = 2
        y[latent > thresh[1]] = 3
        y[latent > thresh[2]] = 4
    elif levels == 2:
        y = (latent > 0).astype(int) + 1
    else:
        thresh = sorted(rng.uniform(-1.0, 1.0, size=levels - 1))
        y = np.ones(n, dtype=int)
        for k in range(levels - 1):
            y[latent > thresh[k]] = k + 2
    if len(set(y.tolist())) < min(3, levels):
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
    return pvals_arr, conv, ordinal_extra


def build_row(seed, scenario, n, m, rep, effect, causal, pvals, runtime, conv, extra, error_message=""):
    converged, warning_flag, status_class = conv
    if len(pvals) < m:
        padded = np.ones(m)
        padded[:len(pvals)] = pvals[:]
        pvals = padded
    top10 = np.argsort(pvals)[:10]
    causal_p = float(pvals[causal]) if causal is not None and causal < len(pvals) else "NA"
    causal_rank = int(np.where(np.argsort(pvals) == causal)[0][0] + 1) if causal is not None and causal < len(pvals) else "NA"
    top1_hit = int(int(np.argsort(pvals)[0]) == int(causal)) if causal is not None else "NA"
    top10_hit = 1 if causal is not None and causal in top10 else (0 if causal is not None else "NA")
    window_hit = int(any(abs(int(x) - int(causal)) <= 5 for x in top10)) if causal is not None else "NA"
    return {
        "engine": "cosmicscanner_psr",
        "seed": seed,
        "scenario": scenario,
        "n": n,
        "m": m,
        "rep": rep,
        "effect": effect,
        "runtime_sec": runtime,
        "converged": converged,
        "warning_flag": warning_flag,
        "status_class": status_class,
        "error_message": error_message,
        "causal_snp_index": causal,
        "causal_p": causal_p,
        "causal_rank": causal_rank,
        "top1_hit": top1_hit,
        "top10_hit": top10_hit,
        "window_hit": window_hit,
        "lambda_gc": lambda_gc(pvals),
        "valid_snp_count": len(pvals),
        "ordinal_method": extra.get("ordinal_method", "NA"),
        "ordinal_stacked_engine": extra.get("ordinal_stacked_engine", "NA"),
        "ordinal_vc_method": extra.get("ordinal_vc_method", "NA"),
        "sigma2_g": extra.get("sigma2_g", "NA"),
        "sigma2_g_boundary": extra.get("sigma2_g_boundary", "NA"),
        "has_blup": extra.get("has_blup", 0),
        "vc_update": extra.get("vc_update", 0),
        "n_stack": extra.get("n_stack", "NA"),
        "reml_loglik": extra.get("reml_loglik", "NA"),
    }


def write_tsv(path, rows):
    if not rows:
        path.write_text("")
        return
    fieldnames = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize_methods(rows):
    methods = sorted(set(r["ordinal_method"] for r in rows))
    summary_rows = []
    for meth in methods:
        mr = [r for r in rows if r["ordinal_method"] == meth]
        runtimes = [r["runtime_sec"] for r in mr if isinstance(r["runtime_sec"], (int, float))]
        conv_vals = [r["converged"] for r in mr if isinstance(r["converged"], (int, float))]
        valid_cnts = [r["valid_snp_count"] for r in mr if isinstance(r["valid_snp_count"], (int, float))]
        sigmas = [r["sigma2_g"] for r in mr if isinstance(r["sigma2_g"], (int, float))]
        boundaries = [r["sigma2_g_boundary"] for r in mr if r["sigma2_g_boundary"] not in ("NA", None) and isinstance(r["sigma2_g_boundary"], (int, float))]
        lams = [r["lambda_gc"] for r in mr if isinstance(r["lambda_gc"], (int, float))]
        top10_hits = [r["top10_hit"] for r in mr if isinstance(r["top10_hit"], int)]
        blups = [r["has_blup"] for r in mr if isinstance(r["has_blup"], int)]
        statuses = [r["status_class"] for r in mr]
        row = {
            "ordinal_method": meth,
            "n_runs": len(mr),
            "runtime_mean": float(np.mean(runtimes)) if runtimes else "NA",
            "runtime_sd": float(np.std(runtimes, ddof=1)) if len(runtimes) > 1 else "NA",
            "runtime_min": float(np.min(runtimes)) if runtimes else "NA",
            "runtime_max": float(np.max(runtimes)) if runtimes else "NA",
            "convergence_rate": float(np.mean(conv_vals)) if conv_vals else "NA",
            "valid_snp_rate": float(np.mean(valid_cnts) / BENCH_M) if valid_cnts else "NA",
            "sigma2_g_mean": float(np.mean(sigmas)) if sigmas else "NA",
            "sigma2_g_sd": float(np.std(sigmas, ddof=1)) if len(sigmas) > 1 else "NA",
            "sigma2_g_min": float(np.min(sigmas)) if sigmas else "NA",
            "sigma2_g_max": float(np.max(sigmas)) if sigmas else "NA",
            "sigma2_g_boundary_rate": float(np.mean(boundaries)) if boundaries else "NA",
            "lambda_gc_mean": float(np.mean(lams)) if lams else "NA",
            "top10_hit_rate": float(np.mean(top10_hits)) if top10_hits else "NA",
            "blup_finite_rate": float(np.mean(blups)) if blups else "NA",
        }
        for sc in sorted(set(statuses)):
            row[f"status_{sc}"] = sum(1 for s in statuses if s == sc)
        summary_rows.append(row)
    return summary_rows


def generate_report(rows, summary_rows):
    lines = [
        "# PSR-GWAS v0.5 Ordinal Benchmark Report",
        f"Date: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        "",
        "## Configuration",
        f"- Replicates: {BENCH_REPS}",
        f"- Sample size (N): {BENCH_N}",
        f"- SNPs (M): {BENCH_M}",
        f"- Ordinal levels: {BENCH_LEVELS}",
        f"- Base seed: {BENCH_SEED}",
    ]
    methods = sorted(set(r["ordinal_method"] for r in rows))
    lines.append(f"- Methods: {', '.join(methods)}")
    lines.append("- Scenarios: null, power_low (effect=0.35), power_high (effect=0.70)")
    lines.append("")
    lines.append("## Per-Method Summary")
    lines.append("")
    hdr = ["Method", "N", "Conv%", "Runtime(s)", "Sigma2_g", "Lambda_GC", "Top10%", "BLUP%"]
    lines.append("| " + " | ".join(hdr) + " |")
    lines.append("| " + " | ".join(["---"] * len(hdr)) + " |")
    for sr in summary_rows:
        meth = sr["ordinal_method"]
        n = sr["n_runs"]
        conv = f"{sr['convergence_rate'] * 100:.0f}%" if isinstance(sr['convergence_rate'], float) else "NA"
        rt = f"{sr['runtime_mean']:.3f}" if isinstance(sr['runtime_mean'], float) else "NA"
        sg = f"{sr['sigma2_g_mean']:.4f}" if isinstance(sr['sigma2_g_mean'], float) else "NA"
        lam = f"{sr['lambda_gc_mean']:.3f}" if isinstance(sr['lambda_gc_mean'], float) else "NA"
        top10 = f"{sr['top10_hit_rate'] * 100:.0f}%" if isinstance(sr['top10_hit_rate'], float) else "NA"
        blup = f"{sr['blup_finite_rate'] * 100:.0f}%" if isinstance(sr['blup_finite_rate'], float) else "NA"
        lines.append(f"| {meth} | {n} | {conv} | {rt} | {sg} | {lam} | {top10} | {blup} |")
    lines.append("")
    lines.append("## Runtime Comparison")
    lines.append("")
    lines.append("| Method | Mean (s) | SD (s) | Min (s) | Max (s) |")
    lines.append("| --- | --- | --- | --- | --- |")
    for sr in summary_rows:
        m = sr["ordinal_method"]
        mn = f"{sr['runtime_mean']:.3f}" if isinstance(sr['runtime_mean'], float) else "NA"
        sd = f"{sr['runtime_sd']:.3f}" if isinstance(sr['runtime_sd'], float) else "NA"
        lo = f"{sr['runtime_min']:.3f}" if isinstance(sr['runtime_min'], float) else "NA"
        hi = f"{sr['runtime_max']:.3f}" if isinstance(sr['runtime_max'], float) else "NA"
        lines.append(f"| {m} | {mn} | {sd} | {lo} | {hi} |")
    lines.append("")
    lines.append("## Convergence Comparison")
    lines.append("")
    lines.append("| Method | Converged | Total | Rate | Status Distribution |")
    lines.append("| --- | --- | --- | --- | --- |")
    for sr in summary_rows:
        meth = sr["ordinal_method"]
        n_conv = int(sr["convergence_rate"] * sr["n_runs"]) if isinstance(sr['convergence_rate'], float) else "?"
        n_tot = sr["n_runs"]
        rate = f"{sr['convergence_rate'] * 100:.1f}%" if isinstance(sr['convergence_rate'], float) else "NA"
        sc_dist = "; ".join(f"{k.replace('status_', '')}={v}" for k, v in sr.items() if k.startswith("status_"))
        lines.append(f"| {meth} | {n_conv} | {n_tot} | {rate} | {sc_dist} |")
    lines.append("")
    est_methods = [sr for sr in summary_rows if sr.get("sigma2_g_mean") not in ("NA", None)]
    if est_methods:
        lines.append("## Sigma2_g Behavior (Methods that estimate it)")
        lines.append("")
        lines.append("| Method | Mean | SD | Min | Max | Boundary Rate |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for sr in est_methods:
            m = sr["ordinal_method"]
            mu = f"{sr['sigma2_g_mean']:.4f}" if isinstance(sr['sigma2_g_mean'], float) else "NA"
            sd = f"{sr['sigma2_g_sd']:.4f}" if isinstance(sr['sigma2_g_sd'], float) else "NA"
            lo = f"{sr['sigma2_g_min']:.4f}" if isinstance(sr['sigma2_g_min'], float) else "NA"
            hi = f"{sr['sigma2_g_max']:.4f}" if isinstance(sr['sigma2_g_max'], float) else "NA"
            br = f"{sr['sigma2_g_boundary_rate'] * 100:.0f}%" if isinstance(sr['sigma2_g_boundary_rate'], float) else "NA"
            lines.append(f"| {m} | {mu} | {sd} | {lo} | {hi} | {br} |")
        lines.append("")
    lines.append("## Top SNP Recovery (Power Scenarios)")
    lines.append("")
    power_rows = [r for r in rows if r["scenario"].startswith("power")]
    if power_rows:
        lines.append("| Method | Top1 Hits | Top10 Hits | Window Hits | Total |")
        lines.append("| --- | --- | --- | --- | --- |")
        for meth in sorted(set(r["ordinal_method"] for r in power_rows)):
            mr = [r for r in power_rows if r["ordinal_method"] == meth]
            t1 = sum(r["top1_hit"] for r in mr if isinstance(r["top1_hit"], int))
            t10 = sum(r["top10_hit"] for r in mr if isinstance(r["top10_hit"], int))
            wh = sum(r["window_hit"] for r in mr if isinstance(r["window_hit"], int))
            lines.append(f"| {meth} | {t1} | {t10} | {wh} | {len(mr)} |")
        lines.append("")
    lines.append("## Key Observations")
    lines.append("")
    for sr in summary_rows:
        if sr.get("sigma2_g_mean") not in ("NA", None) and isinstance(sr.get("sigma2_g_boundary_rate"), float) and sr["sigma2_g_boundary_rate"] > 0.5:
            lines.append(f"- {sr['ordinal_method']}: high boundary rate ({sr['sigma2_g_boundary_rate'] * 100:.0f}%), sigma2_g may be hitting lower bound.")
    for sr in summary_rows:
        if isinstance(sr.get("convergence_rate"), float) and sr["convergence_rate"] < 1.0:
            lines.append(f"- {sr['ordinal_method']}: convergence rate is {sr['convergence_rate'] * 100:.0f}%, consider more iterations or alternative tuning.")
    for sr in summary_rows:
        if isinstance(sr.get("lambda_gc_mean"), float) and abs(sr["lambda_gc_mean"] - 1.0) > 0.1:
            lines.append(f"- {sr['ordinal_method']}: lambda_GC = {sr['lambda_gc_mean']:.3f} (deviates from expected 1.0 for null).")
    lines.append("")
    lines.append("## Output Files")
    for f in [
        "v05_ordinal_benchmark_runs.tsv",
        "v05_ordinal_method_summary.tsv",
        "v05_ordinal_runtime_summary.tsv",
        "v05_ordinal_convergence_summary.tsv",
        "v05_ordinal_sigma2_summary.tsv",
        "v05_ordinal_top_snp_summary.tsv",
        "v05_ordinal_benchmark_report.md",
    ]:
        lines.append(f"- `{f}`")
    lines.append("")
    lines.append("*Report auto-generated by `run_psr_v05_ordinal_benchmark.py`.*")
    return "\n".join(lines)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    methods = resolve_methods()
    scenarios = [("null", 0.0), ("power_low", 0.35), ("power_high", 0.7)]
    all_rows = []
    for rep in range(BENCH_REPS):
        for sidx, (scenario, effect) in enumerate(scenarios):
            seed = BENCH_SEED + rep + 100 * sidx
            rng = np.random.default_rng(seed)
            prefix = WORK_DIR / f"ordinal.{scenario}.rep{rep}"
            G, causal, payload = simulate_ordinal_data(rng, BENCH_N, BENCH_M, effect, BENCH_LEVELS)
            smoke.write_plink(prefix, G)
            for method in methods:
                config = ORDINAL_METHODS_MAP[method]
                out_method = WORK_DIR / f"ordinal.{scenario}.rep{rep}.cosmic.{method}"
                _start = time.time()
                method_error = ""
                try:
                    pvals, conv, extra = run_ordinal_method(prefix, payload, out_method, method, config, BENCH_M)
                except Exception as e:
                    pvals = np.full(BENCH_M, 1.0)
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
                    method_error = str(e)
                runtime = time.time() - _start
                row = build_row(seed, scenario, BENCH_N, BENCH_M, rep, effect, causal, pvals, runtime, conv, extra, method_error)
                all_rows.append(row)
    summary_rows = summarize_methods(all_rows)

    runtime_rows = [
        {k: r[k] for k in ("ordinal_method", "scenario", "rep", "runtime_sec", "converged", "error_message") if k in r}
        for r in all_rows
    ]
    conv_rows = [
        {k: r[k] for k in ("ordinal_method", "scenario", "rep", "converged", "warning_flag", "status_class", "error_message", "n_stack", "valid_snp_count") if k in r}
        for r in all_rows
    ]
    write_tsv(OUT_DIR / "v05_ordinal_benchmark_runs.tsv", all_rows)
    write_tsv(OUT_DIR / "v05_ordinal_method_summary.tsv", summary_rows)
    write_tsv(OUT_DIR / "v05_ordinal_runtime_summary.tsv", runtime_rows)
    write_tsv(OUT_DIR / "v05_ordinal_convergence_summary.tsv", conv_rows)

    sigma2_rows = [r for r in all_rows if r.get("sigma2_g") not in ("NA", None) and isinstance(r.get("sigma2_g"), (int, float))]
    if sigma2_rows:
        sigma2_export = [
            {k: r[k] for k in ("ordinal_method", "scenario", "rep", "sigma2_g", "sigma2_g_boundary", "reml_loglik", "vc_update") if k in r}
            for r in sigma2_rows
        ]
        write_tsv(OUT_DIR / "v05_ordinal_sigma2_summary.tsv", sigma2_export)

    power_rows = [r for r in all_rows if r["scenario"].startswith("power")]
    if power_rows:
        top_export = [
            {k: r[k] for k in ("ordinal_method", "scenario", "rep", "causal_snp_index", "causal_p", "causal_rank", "top1_hit", "top10_hit", "window_hit") if k in r}
            for r in power_rows
        ]
        write_tsv(OUT_DIR / "v05_ordinal_top_snp_summary.tsv", top_export)

    report = generate_report(all_rows, summary_rows)
    report_path = OUT_DIR / "v05_ordinal_benchmark_report.md"
    report_path.write_text(report + "\n")
    print(f"[PSR-GWAS v0.5] Ordinal benchmark outputs written to {OUT_DIR}")
    print(f"[PSR-GWAS v0.5] Report: {report_path}")


if __name__ == "__main__":
    main()
