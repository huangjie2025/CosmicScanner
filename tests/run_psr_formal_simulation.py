#!/usr/bin/env python3
"""Configurable PSR-GWAS formal simulation driver.

Defaults are intentionally small so the script is safe to smoke-test. Set
PSR_FORMAL_REPS, PSR_FORMAL_N, and PSR_FORMAL_M for paper-level runs.

Environment variables:
  PSR_FORMAL_REPS      Number of replicates (default: 2)
  PSR_FORMAL_N         Sample size (default: 80)
  PSR_FORMAL_M         Number of SNPs (default: 80)
  PSR_FORMAL_TRAITS    Comma-separated list (default: binary,binomial,poisson)
  PSR_FORMAL_SCENARIOS Comma-separated list to filter scenarios (default: all)
  PSR_FORMAL_SEED      Base seed (default: 20260624)
  PSR_FORMAL_OUTDIR    Output directory override (default: tests/results/psr_formal_simulation)
  PSR_FORMAL_ORDINAL_METHODS Ordinal methods: "all", comma-separated list, or single method (default: "fixed-gls")
                           Methods: cumulative-binary, fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml
                           When unset, falls back to PSR_FORMAL_ORDINAL_ENGINE / PSR_FORMAL_ORDINAL_VC_METHOD
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
OUT_DIR = Path(os.environ.get("PSR_FORMAL_OUTDIR", ROOT / "tests" / "results" / "psr_formal_simulation"))
WORK_DIR = OUT_DIR / "work"
ENGINE = "cosmicscanner_psr"
BASE_SEED = int(os.environ.get("PSR_FORMAL_SEED", "20260624"))
SCENARIO_FILTER = [x.strip() for x in os.environ.get("PSR_FORMAL_SCENARIOS", "").split(",") if x.strip()] or None

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
        "extra_args": ["--ordinal-vc-grid-min", "0", "--ordinal-vc-grid-max", "2", "--ordinal-vc-grid-points", "11"],
    },
}


def resolve_ordinal_methods():
    """Return list of ordinal method names to run.

    Respects PSR_FORMAL_ORDINAL_METHODS env var. When unset, falls back
    to the older PSR_FORMAL_ORDINAL_ENGINE / PSR_FORMAL_ORDINAL_VC_METHOD vars.
    """
    if "PSR_FORMAL_ORDINAL_METHODS" in os.environ:
        raw = os.environ["PSR_FORMAL_ORDINAL_METHODS"].strip()
        if raw == "all":
            return list(ORDINAL_METHODS_MAP.keys())
        return [x.strip() for x in raw.split(",") if x.strip()]
    engine = os.environ.get("PSR_FORMAL_ORDINAL_ENGINE")
    vc_method = os.environ.get("PSR_FORMAL_ORDINAL_VC_METHOD")
    if engine is not None or vc_method is not None:
        engine = engine or "fixed-gls"
        vc_method = vc_method or "fixed"
        if engine == "both":
            return ["fixed-gls", "dense-zkz-fixed"]
        if engine == "dense-zkz":
            return ["dense-zkz-grid-reml"] if vc_method == "grid-reml" else ["dense-zkz-fixed"]
        return ["fixed-gls"]
    return ["fixed-gls"]


def ordinal_display_mode(method_name, config):
    if config["mode"] == "experimental-cumulative-binary":
        return "experimental_cumulative_binary"
    engine_part = "fixed_gls" if config["stacked_engine"] == "fixed-gls" else "dense_zkz"
    vc_part = "_vc_grid_reml" if config.get("vc_method") == "grid-reml" else ""
    return f"stacked_psr_{engine_part}{vc_part}"


def norm_ppf(p):
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


def run_ordinal_method(prefix, payload, out_prefix, method_name, config, m):
    """Run CosmicScanner for one ordinal method and return parsed outputs.

    Returns (pvals, (converged, warning_flag, status_class), ordinal_extra_dict).
    Raises on subprocess failure.
    """
    pheno = smoke.add_ext(prefix, ".ordinal.pheno")
    smoke.write_one_col_pheno(pheno, payload["pheno"])
    cmd = [str(smoke.BIN), "--bfile", str(prefix), "--gwas", "--no-loco", "--maf", "0.0",
           "--psr-maxit", "20", "--ordinal-model", "--ordinal-pheno", str(pheno),
           "--link", "logit", "--out", str(out_prefix)]
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
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    pvals = []
    with smoke.add_ext(out_prefix, ".gwas.tsv").open("r", encoding="utf-8") as f:
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


def build_result_row(seed, trait_type, scenario, n, m, rep, effect, causal, pvals, runtime, conv, has_offset, link, payload, extra, error_message):
    converged, warning_flag, status_class = conv
    if len(pvals) < m:
        padded = np.ones(m)
        padded[:len(pvals)] = pvals[:]
        pvals = padded
    top = np.argsort(pvals)[:10]
    causal_snp_idx = causal + 1 if causal is not None else -1
    causal_snp_name = f"sim{causal_snp_idx:06d}" if causal is not None else "NA"
    causal_p = float(pvals[causal]) if causal is not None and causal < len(pvals) else "NA"
    causal_rank = int(np.where(np.argsort(pvals) == causal)[0][0] + 1) if causal is not None and causal < len(pvals) else "NA"
    top1_hit = int(int(np.argsort(pvals)[0]) == int(causal)) if causal is not None else "NA"
    top10_hit = 1 if causal is not None and causal in top else (0 if causal is not None else "NA")
    window_hit = int(any(abs(int(x) - int(causal)) <= 5 for x in top)) if causal is not None else "NA"
    power_val = int(causal in top) if scenario.startswith("power") and causal is not None else "NA"
    case_ratio = scenario_case_ratio_or_mean(trait_type, payload) if payload else "NA"
    trial_scheme = scenario_trial_scheme(trait_type, payload) if trait_type == "binomial" and payload else "NA"
    ord_levels = len(set(payload.get("pheno", []))) if trait_type == "ordinal" and payload else "NA"
    row = {
        "engine": ENGINE, "seed": seed, "trait_type": trait_type, "scenario": scenario,
        "n": n, "m": m, "rep": rep, "link": link,
        "runtime": runtime, "converged": converged,
        "warning_flag": warning_flag, "status_class": status_class,
        "error_message": error_message,
        "effect_size": effect,
        "case_ratio_or_mean": case_ratio,
        "trial_scheme": trial_scheme,
        "has_offset": int(has_offset),
        "ordinal_levels": ord_levels,
        "causal_snp": causal_snp_name,
        "causal_p": causal_p,
        "causal_rank": causal_rank,
        "top1_hit": top1_hit, "top10_hit": top10_hit, "window_hit": window_hit,
        "type1_1e_2": float(np.mean(pvals < 1e-2)),
        "type1_1e_3": float(np.mean(pvals < 1e-3)),
        "type1_rate": float(np.mean(pvals < 0.05)) if scenario.startswith("null") else "NA",
        "power": power_val,
        "lambda_gc": lambda_gc(pvals),
        "min_p": float(np.min(pvals)),
        "median_neglog10p": float(np.median(-np.log10(np.clip(pvals, 1e-300, 1.0)))),
        "pvals": pvals,
        "ordinal_mode": extra.get("ordinal_mode", "NA"),
        "ordinal_implementation": extra.get("ordinal_implementation", "NA"),
        "ordinal_method": extra.get("ordinal_method", "NA"),
        "ordinal_stacked_engine": extra.get("ordinal_stacked_engine", "NA"),
        "ordinal_vc_method": extra.get("ordinal_vc_method", "NA"),
        "sigma2_g": extra.get("sigma2_g", "NA"),
        "sigma2_g_boundary": extra.get("sigma2_g_boundary", "NA"),
        "has_blup": extra.get("has_blup", 0),
        "n_stack": extra.get("n_stack", "NA"),
        "reml_loglik": extra.get("reml_loglik", "NA"),
        "n_grid_valid": extra.get("n_grid_valid", "NA"),
        "vc_update": extra.get("vc_update", 0),
    }
    return row


def run_one(seed, trait_type, scenario, n, m, rep):
    rng = np.random.default_rng(seed)
    effect = scenario_effect(scenario)
    prefix = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}"
    out = WORK_DIR / f"{trait_type}.{scenario}.rep{rep}.cosmic"
    link = "probit"
    has_offset = False
    error_message = ""
    causal = None
    payload = {}
    G = None
    results = []

    try:
        if trait_type == "ordinal":
            G, causal, payload = simulate_ordinal_data(rng, n, m, effect)
        else:
            if scenario.startswith("link_logit"):
                link = "logit"
            G, causal, payload = smoke.simulate_data(rng, n, m, trait_type, effect)

        if trait_type == "binary":
            if scenario == "null_unbalanced":
                payload["pheno"] = rng.binomial(1, 0.1, size=n)
                if payload["pheno"].min() == payload["pheno"].max():
                    payload["pheno"][0] = 0
                    payload["pheno"][-1] = 1
            elif scenario == "null_rare":
                payload["pheno"] = rng.binomial(1, 0.05, size=n)
                if payload["pheno"].min() == payload["pheno"].max():
                    payload["pheno"][0] = 0
                    payload["pheno"][-1] = 1
            elif scenario == "null_balanced":
                payload["pheno"] = np.r_[np.zeros(n // 2, dtype=int), np.ones(n - n // 2, dtype=int)]
                rng.shuffle(payload["pheno"])

        if trait_type == "poisson":
            if scenario == "null_low_mean":
                payload["pheno"] = rng.poisson(0.2, size=n)
            elif scenario == "null_high_mean":
                payload["pheno"] = rng.poisson(10.0, size=n)
            elif scenario == "null_offset":
                payload["pheno"] = rng.poisson(2.0, size=n)
                has_offset = True

        if trait_type == "binomial":
            if scenario == "null_variable_trial":
                trials = rng.integers(3, 31, size=n)
                payload = {"events": rng.binomial(trials, 0.25), "trials": trials}
            elif scenario == "null_rare_prop":
                trials = rng.integers(10, 21, size=n)
                payload = {"events": rng.binomial(trials, 0.05), "trials": trials}
            elif scenario == "null_balanced_prop":
                trials = rng.integers(10, 21, size=n)
                payload = {"events": rng.binomial(trials, 0.5), "trials": trials}

        smoke.write_plink(prefix, G)

        if trait_type == "ordinal":
            methods = resolve_ordinal_methods()
            for method in methods:
                config = ORDINAL_METHODS_MAP[method]
                out_method = out.parent / f"{out.name}.{method}"
                _start = time.time()
                method_error = ""
                try:
                    pvals, conv, ord_extra = run_ordinal_method(prefix, payload, out_method, method, config, m)
                except Exception as e:
                    pvals = np.full(m, 1.0)
                    conv = (0, "ERROR", "ERROR")
                    ord_extra = {
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
                ord_mode = ordinal_display_mode(method, config)
                ord_extra["ordinal_mode"] = ord_mode
                ord_extra["ordinal_implementation"] = ord_mode
                row = build_result_row(seed, trait_type, scenario, n, m, rep, effect,
                                       causal, pvals, runtime, conv, has_offset, link,
                                       payload, ord_extra, method_error)
                results.append(row)
        else:
            pvals, runtime, conv = smoke.run_cosmic(prefix, trait_type, payload, out, link=link)
            runtime = max(runtime, time.time() - start_time)
            converged, warning_flag, status_class = conv
            row = build_result_row(seed, trait_type, scenario, n, m, rep, effect,
                                   causal, pvals, runtime, conv, has_offset, link,
                                   payload, {}, error_message)
            results.append(row)
    except Exception as e:
        runtime = time.time() - start_time
        error_message = str(e)
        pvals = np.full(m, 1.0)
        conv = (0, "ERROR", "ERROR")
        if trait_type == "ordinal":
            for method in resolve_ordinal_methods():
                config = ORDINAL_METHODS_MAP.get(method, {})
                ord_extra = {
                    "ordinal_method": method,
                    "ordinal_stacked_engine": config.get("stacked_engine", "NA") or "NA",
                    "ordinal_vc_method": config.get("vc_method", "NA") or "NA",
                    "sigma2_g": "NA", "sigma2_g_boundary": "NA",
                    "has_blup": 0, "n_stack": "NA",
                    "reml_loglik": "NA", "n_grid_valid": "NA",
                    "vc_update": 0,
                    "ordinal_mode": ordinal_display_mode(method, config) if config else "NA",
                    "ordinal_implementation": ordinal_display_mode(method, config) if config else "NA",
                }
                row = build_result_row(seed, trait_type, scenario, n, m, rep, effect,
                                       None, pvals, runtime, conv, False, "probit",
                                       {}, ord_extra, error_message)
                results.append(row)
        else:
            row = build_result_row(seed, trait_type, scenario, n, m, rep, effect,
                                   None, pvals, runtime, conv, False, "probit",
                                   {}, {}, error_message)
            results.append(row)

    for row in results:
        if trait_type != "ordinal":
            row["ordinal_method"] = "NA"
            row["ordinal_stacked_engine"] = "NA"
            row["ordinal_vc_method"] = "NA"
            row["ordinal_mode"] = "NA"
            row["ordinal_implementation"] = "NA"
            row["sigma2_g"] = "NA"
            row["sigma2_g_boundary"] = "NA"
            row["has_blup"] = 0
            row["n_stack"] = "NA"
            row["reml_loglik"] = "NA"
            row["n_grid_valid"] = "NA"
            row["vc_update"] = 0
    return results


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
    if not rows:
        path.write_text("")
        return
    clean = [{k: v for k, v in row.items() if k not in exclude} for row in rows]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=list(clean[0].keys()))
        writer.writeheader()
        writer.writerows(clean)


def generate_report(rows, trait_order):
    report_path = OUT_DIR / "psr_formal_simulation_report.md"
    null_rows = [r for r in rows if r["scenario"].startswith("null")]
    power_rows = [r for r in rows if r["scenario"].startswith("power")]
    lines = ["# PSR-GWAS Formal Simulation Report", f"Date: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}", ""]
    total_reps = len(set((r["trait_type"], r["scenario"], r["rep"]) for r in rows))
    config = f"Default smoke run (REPS={len(set(r['rep'] for r in rows))}, N={rows[0]['n']}, M={rows[0]['m']})"
    if int(os.environ.get("PSR_FORMAL_REPS", "2")) > 5:
        config = f"Formal run (REPS={os.environ.get('PSR_FORMAL_REPS')}, N={rows[0]['n']}, M={rows[0]['m']})"
    lines.append(f"Configuration: {config}")
    lines.append("")
    for tt in trait_order:
        tt_rows = [r for r in rows if r["trait_type"] == tt]
        if not tt_rows:
            continue
        scenarios = sorted(set(r["scenario"] for r in tt_rows))
        lines.append(f"## Trait: {tt}")
        lines.append(f"- Number of scenarios: {len(scenarios)}")
        lines.append(f"- Scenarios: {', '.join(scenarios)}")
        if tt == "ordinal":
            methods = sorted(set(r["ordinal_method"] for r in tt_rows if r["ordinal_method"] != "NA"))
            lines.append(f"- Ordinal methods ({len(methods)}): {', '.join(methods)}")
            for meth in methods:
                meth_rows = [r for r in tt_rows if r["ordinal_method"] == meth]
                n_conv = sum(r["converged"] for r in meth_rows)
                n_total = len(meth_rows)
                conv_rate = n_conv / n_total * 100 if n_total else 0
                meth_null = [r for r in null_rows if r["trait_type"] == tt and r["ordinal_method"] == meth]
                meth_power = [r for r in power_rows if r["trait_type"] == tt and r["ordinal_method"] == meth]
                lines.append(f"  ### Method: {meth}")
                lines.append(f"  - Convergence rate: {n_conv}/{n_total} ({conv_rate:.1f}%)")
                if meth_null:
                    type1_05 = np.mean([r["type1_rate"] for r in meth_null if isinstance(r["type1_rate"], (int, float))])
                    type1_1e2 = np.mean([r["type1_1e_2"] for r in meth_null if isinstance(r["type1_1e_2"], (int, float))])
                    lam_gc = np.mean([r["lambda_gc"] for r in meth_null if isinstance(r["lambda_gc"], (int, float))])
                    lines.append(f"  - Type I error @ 0.05: {type1_05:.4f}")
                    lines.append(f"  - Type I error @ 1e-2: {type1_1e2:.4f}")
                    lines.append(f"  - Mean lambda GC: {lam_gc:.4f}")
                if meth_power:
                    power_vals = [r["power"] for r in meth_power if isinstance(r["power"], (int, float))]
                    top1_hits = sum(r["top1_hit"] for r in meth_power if isinstance(r["top1_hit"], int))
                    top10_hits = sum(r["top10_hit"] for r in meth_power if isinstance(r["top10_hit"], int))
                    win_hits = sum(r["window_hit"] for r in meth_power if isinstance(r["window_hit"], int))
                    lines.append(f"  - Power (causal in top10): {sum(power_vals)}/{len(power_vals)} ({np.mean(power_vals)*100:.0f}%)")
                    lines.append(f"  - Top1 hits: {top1_hits}/{len(power_vals)}")
                    lines.append(f"  - Top10 hits: {top10_hits}/{len(power_vals)}")
                    lines.append(f"  - Window hits (causal +/- 5 SNPs): {win_hits}/{len(power_vals)}")
                meth_runtime = [r["runtime"] for r in meth_rows if isinstance(r["runtime"], (int, float))]
                if meth_runtime:
                    lines.append(f"  - Mean runtime: {np.mean(meth_runtime):.3f}s (range {np.min(meth_runtime):.3f}-{np.max(meth_runtime):.3f}s)")
                sigma2_vals = [r["sigma2_g"] for r in meth_rows if isinstance(r["sigma2_g"], (int, float))]
                if sigma2_vals:
                    lines.append(f"  - Mean sigma2_g: {np.mean(sigma2_vals):.4f} (range {np.min(sigma2_vals):.4f}-{np.max(sigma2_vals):.4f})")
                lines.append("")
        else:
            n_conv = sum(r["converged"] for r in tt_rows)
            n_total = len(tt_rows)
            conv_rate = n_conv / n_total * 100 if n_total else 0
            lines.append(f"- Convergence rate: {n_conv}/{n_total} ({conv_rate:.1f}%)")
            tt_null = [r for r in null_rows if r["trait_type"] == tt]
            if tt_null:
                type1_05 = np.mean([r["type1_rate"] for r in tt_null if isinstance(r["type1_rate"], (int, float))])
                type1_1e2 = np.mean([r["type1_1e_2"] for r in tt_null if isinstance(r["type1_1e_2"], (int, float))])
                lam_gc = np.mean([r["lambda_gc"] for r in tt_null if isinstance(r["lambda_gc"], (int, float))])
                lines.append(f"- Type I error @ 0.05: {type1_05:.4f}")
                lines.append(f"- Type I error @ 1e-2: {type1_1e2:.4f}")
                lines.append(f"- Mean lambda GC: {lam_gc:.4f}")
            tt_power = [r for r in power_rows if r["trait_type"] == tt]
            if tt_power:
                power_vals = [r["power"] for r in tt_power if isinstance(r["power"], (int, float))]
                top1_hits = sum(r["top1_hit"] for r in tt_power if isinstance(r["top1_hit"], int))
                top10_hits = sum(r["top10_hit"] for r in tt_power if isinstance(r["top10_hit"], int))
                win_hits = sum(r["window_hit"] for r in tt_power if isinstance(r["window_hit"], int))
                lines.append(f"- Power (causal in top10): {sum(power_vals)}/{len(power_vals)} ({np.mean(power_vals)*100:.0f}%)")
                lines.append(f"- Top1 hits: {top1_hits}/{len(power_vals)}")
                lines.append(f"- Top10 hits: {top10_hits}/{len(power_vals)}")
                lines.append(f"- Window hits (causal +/- 5 SNPs): {win_hits}/{len(power_vals)}")
            tt_runtime = [r["runtime"] for r in tt_rows if isinstance(r["runtime"], (int, float))]
            if tt_runtime:
                lines.append(f"- Mean runtime: {np.mean(tt_runtime):.3f}s (range {np.min(tt_runtime):.3f}-{np.max(tt_runtime):.3f}s)")
            lines.append("")

    has_issues = [r for r in rows if r["status_class"] != "OK"]
    lines.append("## Warning / Error Summary")
    if has_issues:
        lines.append(f"- Total runs with non-OK status: {len(has_issues)}/{len(rows)}")
        for sc in sorted(set(r["status_class"] for r in has_issues)):
            cnt = sum(1 for r in has_issues if r["status_class"] == sc)
            lines.append(f"  - {sc}: {cnt}")
        error_rows = [r for r in has_issues if r["error_message"]]
        if error_rows:
            lines.append("### Errors")
            for r in error_rows[:10]:
                lines.append(f"- {r['trait_type']}/{r['scenario']} rep{r['rep']} [{r.get('ordinal_method', 'NA')}]: {r['error_message']}")
    else:
        lines.append("- All runs OK")
    lines.append("")

    lines.append("## Output Files")
    for f in ["formal_type1_summary.tsv", "formal_power_summary.tsv", "formal_qq_points.tsv",
              "formal_runtime_summary.tsv", "formal_convergence_summary.tsv",
              "formal_ordinal_method_comparison.tsv", "formal_ordinal_sigma2_summary.tsv",
              "psr_formal_simulation_report.md"]:
        lines.append(f"- `{f}`")
    lines.append("")
    lines.append("## Notes")
    lines.append("- This report was auto-generated by `run_psr_formal_simulation.py`.")
    lines.append("- Results are on the latent/heritability scale used in simulation, not on the observed scale.")
    lines.append("- Formal simulation uses PSR/PQL approximation (pseudo-response + working weights + weighted LMM scan).")
    lines.append("- Binary/binomial/Poisson are stable traits; ordinal methods: cumulative-binary, fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml.")
    lines.append("- For paper-level results, set PSR_FORMAL_REPS=20+ and PSR_FORMAL_N=200+.")
    lines.append("- R baseline comparison is separate (see accuracy/speed benchmark).")
    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    n = int(os.environ.get("PSR_FORMAL_N", "80"))
    m = int(os.environ.get("PSR_FORMAL_M", "80"))
    reps = int(os.environ.get("PSR_FORMAL_REPS", "2"))
    traits = [x.strip() for x in os.environ.get("PSR_FORMAL_TRAITS", "binary,binomial,poisson").split(",") if x.strip()]
    scenarios = {
        "binary": ["null_balanced", "null_unbalanced", "null_rare", "power_low", "power_high", "link_logit_power_low", "link_logit_power_high"],
        "poisson": ["null_low_mean", "null_medium_mean", "null_high_mean", "null_offset", "power_low", "power_high"],
        "binomial": ["null_fixed_trial", "null_variable_trial", "null_rare_prop", "null_balanced_prop", "power_low", "power_high"],
        "ordinal": ["null_balanced", "power_low", "power_high"],
    }
    rows, qq_rows = [], []
    ordinal_method_comparison_rows = []
    trait_order = ["binary", "binomial", "poisson", "ordinal"]
    global start_time
    for trait_type in traits:
        if trait_type not in scenarios:
            raise SystemExit(f"Unknown PSR_FORMAL_TRAITS entry: {trait_type}")
        scenario_list = scenarios[trait_type]
        if SCENARIO_FILTER:
            scenario_list = [s for s in scenario_list if s in SCENARIO_FILTER]
        for sidx, scenario in enumerate(scenario_list):
            for rep in range(reps):
                seed = BASE_SEED + rep + 100 * sidx + 1000 * trait_order.index(trait_type)
                start_time = time.time()
                results_list = run_one(seed, trait_type, scenario, n, m, rep)
                for row in results_list:
                    rows.append(row)
                    if scenario.startswith("null"):
                        obs = np.sort(-np.log10(np.clip(row["pvals"], 1e-300, 1.0)))
                        exp = -np.log10((np.arange(1, len(obs) + 1) - 0.5) / len(obs))[::-1]
                        for k in np.linspace(0, len(obs) - 1, 50, dtype=int):
                            qq_rows.append({k2: row[k2] for k2 in ("engine", "seed", "trait_type", "scenario", "n", "m", "rep", "runtime", "converged", "status_class", "link", "ordinal_mode", "ordinal_implementation", "sigma2_g", "has_blup", "vc_update", "ordinal_method", "error_message")} | {"qq_index": int(k), "expected_neglog10p": float(exp[k]), "observed_neglog10p": float(obs[k])})
                    if trait_type == "ordinal":
                        ordinal_method_comparison_rows.append({k: row[k] for k in ("engine", "seed", "trait_type", "scenario", "n", "m", "rep", "ordinal_method", "ordinal_stacked_engine", "ordinal_vc_method", "sigma2_g", "sigma2_g_boundary", "has_blup", "n_stack", "reml_loglik", "n_grid_valid", "vc_update", "runtime", "converged", "lambda_gc", "type1_rate", "power", "top1_hit", "top10_hit", "window_hit", "error_message")})
    null_rows = [r for r in rows if r["scenario"].startswith("null")]
    power_rows = [r for r in rows if r["scenario"].startswith("power")]
    write(OUT_DIR / "formal_type1_summary.tsv", null_rows, exclude=("pvals", "causal_p", "causal_rank"))
    write(OUT_DIR / "formal_power_summary.tsv", power_rows, exclude=("pvals", "type1_1e_2", "type1_1e_3", "lambda_gc"))
    write(OUT_DIR / "formal_runtime_summary.tsv", rows, exclude=("pvals",))
    write(OUT_DIR / "formal_convergence_summary.tsv", rows, exclude=("pvals",))
    write(OUT_DIR / "formal_qq_points.tsv", qq_rows)
    if ordinal_method_comparison_rows:
        write(OUT_DIR / "formal_ordinal_method_comparison.tsv", ordinal_method_comparison_rows)
        ord_sigma2 = [r for r in ordinal_method_comparison_rows if r.get("sigma2_g") != "NA" and isinstance(r.get("sigma2_g"), (int, float))]
        if ord_sigma2:
            write(OUT_DIR / "formal_ordinal_sigma2_summary.tsv", ord_sigma2)
    report_path = generate_report(rows, trait_order)
    print(f"[PSR-GWAS] Formal simulation outputs written to {OUT_DIR}")
    print(f"[PSR-GWAS] Report: {report_path}")


start_time = time.time()

if __name__ == "__main__":
    main()
