#!/usr/bin/env python3
"""R31 BGEN AD final statistical/production acceptance gate.

Research/internal explicit harness only. The ordinary BGEN + additive_dominance
CLI guard is intentionally not modified.
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import stats

from run_bgen_ad_null_validation import phenotype_matrix, write_hardcalls
from run_bgen_dominance_production_gate import (
    generate_hardcalls, posterior_probabilities, require, write_bgen)
from run_bgen_ad_loco_calibration_root_cause import binmatrix


POLICY = "BGEN_DOM_V1_SHARED_AD_RESEARCH"
QC_POLICY_VERSION = "BGEN_DOM_V1"
AD_QC_POLICY_VERSION = "BGEN_AD_SHARED_ELIGIBILITY_RESEARCH_V1"
GATE_VERSION = "R31"
HEADER = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
          "BETA_A", "SE_A", "P_A", "LOG10P_A",
          "BETA_D", "SE_D", "P_D", "LOG10P_D",
          "P_JOINT", "LOG10P_JOINT"]


def read_meta(path):
    result = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[row["Item"]] = row["Value"]
    return result


def parse_assoc(path, chromosome=None):
    result = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == HEADER, f"unified header changed in {path}")
        for row in reader:
            if chromosome is not None and row["CHR"] != chromosome:
                continue
            result[row["SNP"]] = row
    return result


def read_null(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def replicate_metrics(pvalues, lam_key_df):
    p = np.asarray(pvalues, dtype=float)
    if len(p) == 0:
        return None
    statistic = stats.chi2.isf(np.clip(p, 1e-300, 1.0), lam_key_df)
    return {"type1": float(np.mean(p < .05)), "lambda_gc": float(np.median(statistic) / stats.chi2.ppf(.5, lam_key_df))}


def summarize_replicates(per_rep, stat_key, df):
    values_type = [row[stat_key] for row in per_rep]
    values_lambda = [row["lambda_gc"] for row in per_rep]
    rng = np.random.default_rng(4771)
    type_ci = []
    lambda_ci = []
    for _ in range(5000):
        index = rng.integers(0, len(per_rep), len(per_rep))
        type_ci.append(float(np.mean(np.asarray(values_type)[index])))
        lambda_ci.append(float(np.median(np.asarray(values_lambda)[index])))
    result = {
        "replicates": len(per_rep),
        "type1_mean": float(np.mean(values_type)),
        "type1_ci95": [float(v) for v in np.percentile(type_ci, [2.5, 97.5])],
        "lambda_gc_median_of_replicates": float(np.median(values_lambda)),
        "lambda_gc_ci95": [float(v) for v in np.percentile(lambda_ci, [2.5, 97.5])],
        "lambda_gc_mean": float(np.mean(values_lambda)),
    }
    return result


def lambda_gc(pvalues, df):
    p = np.asarray(pvalues, dtype=float)
    statistic = stats.chi2.isf(np.clip(p, 1e-300, 1.0), df)
    return float(np.median(statistic) / stats.chi2.ppf(.5, df))


def bound_ci(rate, reps, alpha=.05):
    lower = stats.beta.ppf(alpha / 2, rate * reps + .5, (1 - rate) * reps + .5)
    upper = stats.beta.ppf(1 - alpha / 2, rate * reps + .5, (1 - rate) * reps + .5)
    return float(lower), float(upper)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--harness", required=True)
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--calibration-driver", required=True)
    parser.add_argument("--null-driver", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--m", type=int, default=1200)
    parser.add_argument("--reps", type=int, default=100,
                        help="retained for compatibility; calibration-reps drives the strict gate")
    parser.add_argument("--calibration-reps", type=int, default=100)
    parser.add_argument("--diagnostic-reps", type=int, default=20)
    parser.add_argument("--partial-reps", type=int, default=100)
    parser.add_argument("--alternative-reps", type=int, default=20)
    parser.add_argument("--workers", type=int, default=48)
    parser.add_argument("--chunk", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--plink2", default=os.environ.get(
        "PLINK2_BIN", "/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2"))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if not args.smoke:
        require(args.n >= 1000, "production gate requires N >= 1000")
        require(args.m >= 1200, "production gate requires M >= 1200")
        require(args.calibration_reps >= 100, "production calibration requires >= 100 replicates")
        require(args.partial_reps >= 20, "partial calibration requires >= 20 replicates")
        require(args.alternative_reps >= 20, "mixed alternative requires >= 20 replicates")
    for executable in (args.harness, args.scanner, args.calibration_driver, args.null_driver, args.plink2):
        require(Path(executable).exists(), f"required executable missing: {executable}")
    if args.smoke:
        args.calibration_reps = min(args.calibration_reps, 4)
        args.diagnostic_reps = min(args.diagnostic_reps, 4)
        args.partial_reps = min(args.partial_reps, 4)
        args.alternative_reps = min(args.alternative_reps, 4)

    root = Path(args.workdir)
    root.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="r31-", dir=root))
    (root / "latest_run.txt").write_text(str(run) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def execute(command, quiet=False):
        started = time.time()
        proc = subprocess.run([str(v) for v in command], env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            if not quiet:
                raise AssertionError(proc.stdout[-8000:])
        return proc, time.time() - started

    start_time = time.time()
    samples = [f"R31_{i:05d}" for i in range(args.n)]
    calls, frequencies, causal = generate_hardcalls(args.n, args.m, args.seed)
    causal_id = f"rs{causal + 1:06d}"
    eligible_chr1 = np.flatnonzero((frequencies > 0.20) & (frequencies < 0.40) &
                                   (np.arange(args.m) < args.m // 2))
    if len(eligible_chr1):
        causal = int(eligible_chr1[len(eligible_chr1) // 2])
        causal_id = f"rs{causal + 1:06d}"
    rng = np.random.default_rng(args.seed + 1)
    noise = rng.normal(size=calls.shape)
    category = rng.choice(3, size=calls.shape, p=(0.50, 0.35, 0.15))
    quality_factor = np.choose(category, (0.55, 1.0, 1.8))
    missing = rng.random(calls.shape) < .01
    missing[0] = rng.random(args.n) < .08
    missing[causal] = False

    genotype_dir = run / "genotypes"
    genotype_dir.mkdir(parents=True)
    bed_prefix = genotype_dir / "panel"
    write_hardcalls(bed_prefix, np.where(missing, -1, calls), samples)
    execute([args.plink2, "--bfile", bed_prefix, "--make-pgen", "--out", genotype_dir / "panel_pgen"])
    execute([args.plink2, "--bfile", bed_prefix, "--export", "bgen-1.2", "--out", genotype_dir / "panel_bgen"])
    pgen_path = genotype_dir / "panel_pgen.pgen"
    onehot_bgen = genotype_dir / "panel_bgen.bgen"

    genotypes = {
        "bed_hardcall": ("bed", bed_prefix),
        "pgen_hardcall": ("pgen", pgen_path),
        "bgen_one_hot": ("bgen", onehot_bgen)}
    for quality, tau in (("bgen_high_confidence", .35),
                         ("bgen_medium_confidence", .65),
                         ("bgen_low_confidence", 1.0)):
        probs = posterior_probabilities(calls, frequencies, tau, noise, quality_factor)
        path = genotype_dir / f"{quality}.bgen"
        write_bgen(path, samples, probs, missing)
        genotypes[quality] = ("bgen", path)

    covariates = np.column_stack([np.ones(args.n), rng.normal(size=args.n)]).astype(np.float32).astype(float)
    covar_path = run / "covariates.txt"
    phenotype_matrix(covar_path, covariates)
    native_covar_path = run / "qcovar.txt"
    with native_covar_path.open("w") as handle:
        for sample, value in zip(samples, covariates[:, 1]):
            handle.write(f"{sample} {sample} {value:.17g}\n")

    # Fixed true-V encoded partitions used only to generate exact matched nulls.
    kernel_dir = run / "reference_kernels"
    kernel_dir.mkdir()
    dummy = run / "dummy.pheno"
    phenotype_matrix(dummy, np.ones((args.n, 1)))
    kernels = {}
    for quality, (fmt, source) in genotypes.items():
        parts = {}
        if quality in {"bed_hardcall", "pgen_hardcall"}:
            chromosome_by_label = {"global": "-", "loco1": "2", "loco2": "1"}
            for label, chromosome in chromosome_by_label.items():
                output = kernel_dir / f"{quality}_{label}"
                proc, _ = execute([args.null_driver, fmt, source, output, "-", chromosome, ".01", "-", "0"],
                                  quiet=True)
                require(proc.returncode == 0, f"null driver failed for {quality}/{label}: {proc.stdout[-3000:]}")
                parts[label] = (binmatrix(str(output) + ".ka.bin", args.n),
                                binmatrix(str(output) + ".kd.bin", args.n))
        else:
            for label, excluded, scan in (("global", "-", "1"), ("loco1", "1", "1"), ("loco2", "2", "2")):
                output = kernel_dir / f"{quality}_{label}"
                proc, _ = execute([args.calibration_driver, source, dummy, covar_path, output, scan, excluded,
                                   "fixed", ".25", ".20", ".55"], quiet=True)
                require(proc.returncode == 0,
                        f"calibration driver failed for {quality}/{label}: {proc.stdout[-3000:]}")
                parts[label] = (binmatrix(str(output) + ".ka.bin", args.n),
                                binmatrix(str(output) + ".kd.bin", args.n))
        kernels[quality] = parts

    # Direct cross-format kernel parity on chromosome 2 (outside partition for target 1).
    direct_dir = run / "direct_kernels"
    direct_dir.mkdir()
    for quality, (fmt, source) in genotypes.items():
        output = direct_dir / quality
        proc, _ = execute([args.null_driver, fmt, source, output, "-", "2", ".01", "-", "0"], quiet=True)
        require(proc.returncode == 0, f"direct kernel driver failed for {quality}")
    direct = {quality: (binmatrix(str(direct_dir / quality) + ".ka.bin", args.n),
                        binmatrix(str(direct_dir / quality) + ".kd.bin", args.n))
              for quality in genotypes}
    direct_errors = {}
    reference = direct["bed_hardcall"]
    for quality, pair in direct.items():
        direct_errors[quality] = {
            "ka_vs_bed": float(np.max(np.abs(pair[0] - reference[0]))),
            "kd_vs_bed": float(np.max(np.abs(pair[1] - reference[1])))}

    native_sources = {
        "bed_hardcall": ("--bfile", bed_prefix),
        "pgen_hardcall": ("--pfile", pgen_path.with_suffix(""))}

    truth = (.25, .20, .55)
    max_reps = max(args.calibration_reps, args.diagnostic_reps, args.partial_reps, args.alternative_reps)
    phenotype_noise = np.random.default_rng(args.seed + 2).standard_normal((args.n, max_reps))
    partial_noise = np.random.default_rng(args.seed + 3).standard_normal((args.n, max_reps))
    alternate_noise = np.random.default_rng(args.seed + 4).standard_normal((args.n, max_reps))

    def covariance(values, ka, kd, a, d, e):
        return a * ka + d * kd + e * np.eye(args.n)

    phenotype_cache = {}
    for quality, parts in kernels.items():
        global_cov = covariance(None, *parts["global"], *truth)
        noise_global = phenotype_noise
        phenotype_cache[(quality, "strict_global")] = (
            np.linalg.cholesky(global_cov) @ noise_global + .4 * covariates[:, 1, None], None, None)
        phenotype_cache[(quality, "strict_loco1")] = (
            np.linalg.cholesky(covariance(None, *parts["loco1"], *truth)) @ noise_global + .4 * covariates[:, 1, None],
            "1", True)
        phenotype_cache[(quality, "strict_loco2")] = (
            np.linalg.cholesky(covariance(None, *parts["loco2"], *truth)) @ noise_global + .4 * covariates[:, 1, None],
            "2", True)
        phenotype_cache[(quality, "partial_additive")] = (
            np.linalg.cholesky(covariance(None, *parts["loco1"], truth[0], 0.0, truth[2])) @ partial_noise + .4 * covariates[:, 1, None],
            "1", True)
        phenotype_cache[(quality, "partial_dominance")] = (
            np.linalg.cholesky(covariance(None, *parts["loco1"], 0.0, truth[1], truth[2])) @ partial_noise + .4 * covariates[:, 1, None],
            "1", True)
        base_alternate = (np.linalg.cholesky(covariance(None, *parts["loco1"], *truth)) @ alternate_noise
                          + .4 * covariates[:, 1, None])
        def marker_columns(j):
            observed = ~missing[j]
            dosage = np.where(observed, calls[j], 0.0).astype(float)
            heterozygote = np.where(observed & (calls[j] == 1), 1.0, 0.0)
            mean_dosage = dosage[observed].mean()
            mean_het = heterozygote[observed].mean()
            return (np.where(observed, dosage - mean_dosage, 0.0),
                    np.where(observed, heterozygote - mean_het, 0.0))

        causal_additive, causal_dominance = marker_columns(causal)
        for arm, beta in (("alternative_additive", (.6, 0.0)),
                          ("alternative_dominance", (0.0, .8)),
                          ("mixed_alternative", (.6, .8))):
            signal = beta[0] * causal_additive + beta[1] * causal_dominance
            phenotype_cache[(quality, arm)] = (base_alternate + signal[:, None], "1", True)

    bgen_qualities = [quality for quality in genotypes if quality.startswith("bgen_")]
    scenarios = {
        "strict_global": {"type": "calibration", "chromosome": None, "loco": False, "test_df": 2,
                          "reps": args.calibration_reps, "gate": True},
        "strict_loco1": {"type": "calibration", "chromosome": "1", "loco": True, "test_df": 2,
                         "reps": args.calibration_reps, "gate": True},
        "strict_loco2": {"type": "calibration", "chromosome": "2", "loco": True, "test_df": 2,
                         "reps": args.diagnostic_reps, "gate": False},
        "partial_additive": {"type": "conditional", "chromosome": "1", "loco": True, "test_df": 1, "field": "P_D",
                             "reps": args.partial_reps if bgen_qualities else args.diagnostic_reps,
                             "gate": True},
        "partial_dominance": {"type": "conditional", "chromosome": "1", "loco": True, "test_df": 1, "field": "P_A",
                              "reps": args.partial_reps if bgen_qualities else args.diagnostic_reps,
                              "gate": True},
        "mixed_alternative": {"type": "alternative", "chromosome": "1", "loco": True,
                              "reps": args.alternative_reps, "gate": False}}
    for name in ("partial_additive", "partial_dominance", "mixed_alternative"):
        if scenarios[name]["reps"] > args.diagnostic_reps:
            scenarios[name]["diagnostic_only_for_native"] = True
        else:
            scenarios[name]["diagnostic_only_for_native"] = False

    task_records = []
    for quality in genotypes:
        native = quality in {"bed_hardcall", "pgen_hardcall"}
        for scenario, spec in scenarios.items():
            if native and (spec["type"] == "alternative" or spec.get("diagnostic_only_for_native")):
                reps_for_quality = args.diagnostic_reps
            else:
                reps_for_quality = spec["reps"]
            values, _, _ = phenotype_cache[(quality, scenario)]
            for start in range(0, reps_for_quality, 1):
                end = min(reps_for_quality, start + 1)
                label = f"{quality}__{scenario}__rep{start:03d}_{end:03d}"
                directory = run / "runs" / label
                directory.mkdir(parents=True, exist_ok=True)
                if native:
                    command = ["__native__", quality, scenario, start, end, spec.get("chromosome") or "",
                               "true" if spec["loco"] else "false"]
                else:
                    phenopath = directory / "pheno.txt"
                    phenotype_matrix(phenopath, values[:, start:end])
                    config = {"maxit": "200", "tol": "1e-6"}
                    config["loco"] = "true" if spec["loco"] else "false"
                    if spec.get("chromosome"):
                        config["chromosome"] = spec["chromosome"]
                    cfg = directory / "config.tsv"
                    cfg.write_text("".join(f"{k}\t{v}\n" for k, v in config.items()))
                    command = [args.harness, genotypes[quality][1], phenopath, covar_path, directory, cfg]
                task_records.append({
                    "quality": quality, "scenario": scenario, "start": start, "end": end,
                    "directory": directory, "spec": spec, "command": command,
                    "runner": "native_cli" if native else "bgen_harness",
                    "label": label})

    collection = {}
    causal_rows = {quality: {"mixed_alternative": []} for quality in genotypes}
    null_entries = []
    runtimes = []
    completed = 0

    def process_task(record):
        local_rows = []
        local_null_rows = []
        started = time.time()
        if record["runner"] == "native_cli":
            quality = record["quality"]
            values, _, _ = phenotype_cache[(quality, record["scenario"])]
            flag, native_source = native_sources[quality]
            scan_chr = record["spec"].get("chromosome")
            loco = record["spec"]["loco"]
            for rep in range(record["start"], record["end"]):
                prefix = record["directory"] / f"rep{rep}"
                pheno = record["directory"] / f"rep{rep}.pheno"
                with pheno.open("w") as handle:
                    for sample, value in zip(samples, values[:, rep]):
                        handle.write(f"{sample} {sample} {value:.17g}\n")
                command = [args.scanner, flag, native_source, "--gwas-model", "additive_dominance",
                           "--pheno", pheno, "--qcovar", native_covar_path,
                           "--maf", ".01", "--reml-maxit", "200", "--reml-tol", "1e-6",
                           "--out", prefix]
                if not loco:
                    command.append("--no-loco")
                if scan_chr:
                    command.extend(["--chr", scan_chr])
                proc, _ = execute(command)
                require(proc.returncode == 0,
                        f"native CLI failed for {record['label']} rep {rep}: {proc.stdout[-3000:]}")
        else:
            proc, _ = execute(record["command"])
            require(proc.returncode == 0,
                    f"harness failed for {record['quality']}/{record['scenario']}: {proc.stdout[-3000:]}")
        for rep in range(record["start"], record["end"]):
            local = rep - record["start"]
            if record["runner"] == "native_cli":
                meta = read_meta(record["directory"] / f"rep{rep}.scan.meta.txt")
                assoc_path = record["directory"] / f"rep{rep}.assoc.ad.txt"
            else:
                meta = read_meta(record["directory"] / f"rep{local}.scan.meta.txt")
                assoc_path = record["directory"] / f"rep{local}.assoc.ad.txt"
            require(meta.get("RUN_STATUS", "SUCCESS") == "SUCCESS", f"failed scan metadata for {record['label']}")
            if "RUN_STATUS" not in meta:
                meta["RUN_STATUS"] = "SUCCESS" if meta.get("VCE_CONVERGED", "false") == "true" else "FAILED"
            if record["runner"] != "native_cli":
                require(meta["PER_SNP_NULL_REFITS"] == "0", "per-SNP null refit observed")
                require(meta["QC_POLICY"] == POLICY, "QC policy changed")
                require(meta["GENOTYPE_BACKEND"] == "BGEN", "BGEN backend metadata changed")
                require(meta["PROBABILITY_MODE"] == "true", "probability metadata changed")
            else:
                require(meta["MODEL"] == "additive_dominance", "native AD model metadata changed")
            rows = parse_assoc(assoc_path, chromosome=record["spec"].get("chromosome"))
            target_causal = parse_assoc(assoc_path, chromosome="1").get(causal_id)
            if record["runner"] != "native_cli":
                null_rows = read_null(record["directory"] / f"rep{local}.ad.null.tsv")
                for null in null_rows:
                    entry = dict(null)
                    entry.update({"quality": record["quality"], "scenario": record["scenario"], "replicate": rep})
                    local_null_rows.append(entry)
            local_rows.append((rep, meta, rows, target_causal))
        return record, time.time() - started, local_rows, local_null_rows

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_task, record) for record in task_records]
        for future in as_completed(futures):
            record, elapsed, local_rows, local_null_rows = future.result()
            completed += 1
            runtimes.append({"quality": record["quality"], "scenario": record["scenario"],
                             "start": record["start"], "end": record["end"], "seconds": elapsed})
            key = (record["quality"], record["scenario"])
            collection.setdefault(key, {"meta": [], "per_rep": [], "causal": []})
            for rep, meta, rows, causal_row in local_rows:
                fields = record["spec"].get("field")
                if record["spec"]["type"] in ("calibration", "conditional"):
                    pvalues = [float(row[fields]) for row in rows.values()] if fields else \
                              [float(row["P_JOINT"]) for row in rows.values()]
                    metric = replicate_metrics(pvalues, record["spec"].get("test_df", 2))
                    require(metric is not None and metric["type1"] is not None,
                            f"no tested markers for {key} replicate {rep}")
                    metric["replicate"] = rep
                    collection[key]["per_rep"].append(metric)
                meta_row = dict(meta)
                meta_row.update({"quality": record["quality"], "scenario": record["scenario"], "replicate": rep})
                collection[key]["meta"].append(meta_row)
                if record["spec"]["type"] == "alternative" and causal_row is not None:
                    rank = 1 + sum(float(r["P_JOINT"]) < float(causal_row["P_JOINT"]) for r in rows.values())
                    causal_rows[record["quality"]][record["scenario"]].append({
                        "replicate": rep,
                        "beta_a": float(causal_row["BETA_A"]), "beta_d": float(causal_row["BETA_D"]),
                        "se_a": float(causal_row["SE_A"]), "se_d": float(causal_row["SE_D"]),
                        "p_add": float(causal_row["P_A"]), "p_dom": float(causal_row["P_D"]),
                        "p_joint": float(causal_row["P_JOINT"]),
                        "joint_rank": rank})
            null_entries.extend(local_null_rows)
            # Keep metadata/assoc/null TSVs; remove huge replicated kernel dumps.
            for path in record["directory"].glob("*.bin"):
                path.unlink()
            if completed % 20 == 0 or completed == len(task_records):
                print(f"R31 progress: {completed}/{len(task_records)} harness chunks completed", flush=True)

    calibration_summary = {}
    conditional_summary = {}
    per_rep_rows = []
    for (quality, scenario), payload in sorted(collection.items()):
        slots = scenarios[scenario]
        per_meta = payload["meta"]
        convergence = all(item.get("NULL_CONVERGED", "false") == "true" for item in per_meta)
        max_gradient = max(float(item.get("NULL_GRADIENT_NORM", "nan")) for item in per_meta
                           if item.get("NULL_GRADIENT_NORM", "nan") != "nan") if any(
            item.get("NULL_GRADIENT_NORM", "nan") != "nan" for item in per_meta) else 0.0
        max_kkt = max(float(item.get("NULL_KKT_RESIDUAL", "nan")) for item in per_meta
                      if item.get("NULL_KKT_RESIDUAL", "nan") != "nan") if any(
            item.get("NULL_KKT_RESIDUAL", "nan") != "nan" for item in per_meta) else 0.0
        record_summary = {
            "quality": quality, "scenario": scenario,
            "runner": "native_cli" if quality in {"bed_hardcall", "pgen_hardcall"} else "bgen_harness",
            "replicates": len(per_meta),
            "null_convergence_rate": float(np.mean([
                item.get("NULL_CONVERGED", item.get("AD_NULL_CONVERGED_ALL", item.get("VCE_CONVERGED", "false"))) == "true"
                for item in per_meta])),
            "per_snp_null_refits": int(max(int(item.get("PER_SNP_NULL_REFITS", "0")) for item in per_meta)),
            "max_gradient_norm": max_gradient, "max_kkt_residual": max_kkt,
            "run_status_success_rate": float(np.mean([item.get("RUN_STATUS") == "SUCCESS" for item in per_meta])),
            "null_labels": sorted({item.get("AD_NULL_MODELS", "1") for item in per_meta}),
        }
        if "SIGMA_A2" in per_meta[0] or "SIGMA_A2_1" in per_meta[0]:
            record_summary["mean_sigma_a2"] = float(np.mean([
                float(item.get("SIGMA_A2_1", item.get("SIGMA_A2", "nan"))) for item in per_meta]))
            record_summary["mean_sigma_d2"] = float(np.mean([
                float(item.get("SIGMA_D2_1", item.get("SIGMA_D2", "nan"))) for item in per_meta]))
            record_summary["mean_sigma_e2"] = float(np.mean([
                float(item.get("SIGMA_E2_1", item.get("SIGMA_E2", "nan"))) for item in per_meta]))
        record_summary["boundary_status"] = sorted({item.get("BOUNDARY_STATUS", "UNKNOWN") for item in per_meta})
        record_summary["weak_information_null_fraction"] = float(np.mean([
            int(item.get("WEAK_INFORMATION_NULLS", "0")) > 0 for item in per_meta]))
        record_summary["d_local_se_weak_fraction"] = float(np.mean([
            "D_local_SE" in item for item in per_meta]))
        if slots["type"] == "calibration":
            summary = summarize_replicates(payload["per_rep"], "type1", 2)
            summary.update(record_summary)
            calibration_summary[(quality, scenario)] = summary
        elif slots["type"] == "conditional":
            summary = summarize_replicates(payload["per_rep"], "type1", 1)
            summary.update(record_summary)
            conditional_summary[(quality, scenario)] = summary
        for metric in payload["per_rep"]:
            per_rep_rows.append({"quality": quality, "scenario": scenario, "replicate": metric["replicate"],
                                 "type1": metric["type1"], "lambda_gc": metric["lambda_gc"]})

    alternative_summary = []
    for quality, arms in causal_rows.items():
        for arm, rows in arms.items():
            if not rows:
                continue
            p_joint = np.asarray([row["p_joint"] for row in rows])
            p_add = np.asarray([row["p_add"] for row in rows])
            p_dom = np.asarray([row["p_dom"] for row in rows])
            alternative_summary.append({
                "quality": quality, "scenario": arm, "replicates": len(rows),
                "mean_beta_a": float(np.mean([row["beta_a"] for row in rows])),
                "mean_beta_d": float(np.mean([row["beta_d"] for row in rows])),
                "mean_se_a": float(np.mean([row["se_a"] for row in rows])),
                "mean_se_d": float(np.mean([row["se_d"] for row in rows])),
                "power_a_at_0.05": float(np.mean(p_add < .05)),
                "power_d_at_0.05": float(np.mean(p_dom < .05)),
                "power_joint_at_0.05": float(np.mean(p_joint < .05)),
                "median_joint_rank": float(np.median([row["joint_rank"] for row in rows])),
                "lambda_gc_joint": lambda_gc(p_joint, 2)})

    strict_keys = [("global", "strict_global"), ("loco_target1", "strict_loco1")]
    strict_pass = True
    strict_details = []
    for quality in genotypes:
        for name, scenario in strict_keys:
            summary = calibration_summary.get((quality, scenario))
            require(summary is not None, f"missing calibration summary {quality}/{scenario}")
            require(summary["replicates"] >= (1 if args.smoke else 100),
                    f"strict calibration below required replicates for {quality}/{scenario}")
            pass_type = .025 <= summary["type1_mean"] <= .075
            pass_lambda = .85 <= summary["lambda_gc_median_of_replicates"] <= 1.15
            pass_converged = summary["null_convergence_rate"] == 1.0
            pass_kkt = summary["runner"] == "native_cli" or summary["max_kkt_residual"] <= 1e-7
            strict_pass = strict_pass and pass_type and pass_lambda and pass_converged and pass_kkt
            strict_details.append({"quality": quality, "scenario": scenario,
                                   "type1_pass": pass_type, "lambda_pass": pass_lambda,
                                   "convergence_pass": pass_converged, "kkt_pass": pass_kkt})
    diagnostic_details = []
    for quality in genotypes:
        summary = calibration_summary.get((quality, "strict_loco2"))
        require(summary is not None, f"missing diagnostic LOCO target 2 summary for {quality}")
        diagnostic_details.append({"quality": quality, "scenario": "strict_loco2",
                                   "replicates": summary["replicates"],
                                   "type1_mean": summary["type1_mean"],
                                   "lambda_gc_median_of_replicates": summary["lambda_gc_median_of_replicates"],
                                   "null_convergence_rate": summary["null_convergence_rate"]})
    partial_pass = True
    partial_details = []
    for quality in bgen_qualities:
        for scenario in ("partial_additive", "partial_dominance"):
            summary = conditional_summary.get((quality, scenario))
            require(summary is not None, f"missing conditional summary {quality}/{scenario}")
            require(summary["replicates"] >= (1 if args.smoke else 20),
                    f"partial calibration below required replicates for {quality}/{scenario}")
            pass_type = .025 <= summary["type1_mean"] <= .075
            pass_lambda = .85 <= summary["lambda_gc_median_of_replicates"] <= 1.15
            pass_converged = summary["null_convergence_rate"] == 1.0
            partial_pass = partial_pass and pass_type and pass_lambda and pass_converged
            partial_details.append({"quality": quality, "scenario": scenario,
                                    "replicates": summary["replicates"],
                                    "type1_pass": pass_type, "lambda_pass": pass_lambda,
                                    "convergence_pass": pass_converged})

    null_frame_columns = ["quality", "scenario", "replicate", "NULL", "SIGMA_A2", "SIGMA_D2", "SIGMA_E2",
                          "LOGL", "CONVERGED", "ITERATIONS", "LIKELIHOOD_DELTA", "GRADIENT_NORM",
                          "KKT_RESIDUAL", "BOUNDARY_STATUS", "INFORMATION_RANK", "D_LOCAL_SE"]
    with (run / "null_certificates.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=null_frame_columns, delimiter="\t")
        writer.writeheader()
        for row in null_entries:
            writer.writerow({key: row.get(key, "") for key in null_frame_columns})
    with (run / "per_replicate_calibration.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["quality", "scenario", "replicate", "type1", "lambda_gc"],
                                delimiter="\t")
        writer.writeheader(); writer.writerows(per_rep_rows)
    calibration_rows = []
    for key, summary in calibration_summary.items():
        row = dict(summary); row["key"] = f"{key[0]}|{key[1]}"; calibration_rows.append(row)
    with (run / "calibration.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(calibration_rows[0].keys()), delimiter="\t")
        writer.writeheader(); writer.writerows(calibration_rows)
    conditional_rows = []
    for key, summary in conditional_summary.items():
        row = dict(summary); row["key"] = f"{key[0]}|{key[1]}"; conditional_rows.append(row)
    with (run / "conditional_partial.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(conditional_rows[0].keys()), delimiter="\t")
        writer.writeheader(); writer.writerows(conditional_rows)
    if alternative_summary:
        with (run / "alternative_summary.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(alternative_summary[0].keys()), delimiter="\t")
            writer.writeheader(); writer.writerows(alternative_summary)

    gate = {
        "gate_version": GATE_VERSION,
        "run_directory": str(run),
        "created_at_epoch": time.time(),
        "runtime_seconds": time.time() - start_time,
        "n": args.n, "m": args.m, "seed": args.seed,
        "replicate_counts": {
            "strict_global_targets": args.calibration_reps,
            "strict_matched_loco_target1": args.calibration_reps,
            "strict_matched_loco_target2_diagnostic": args.diagnostic_reps,
            "partial_null_bgen": args.partial_reps,
            "partial_null_native_diagnostic": args.diagnostic_reps,
            "mixed_alternative": args.alternative_reps},
        "effective_causal_id": causal_id,
        "truth_sigma_a2": .25, "truth_sigma_d2": .20, "truth_sigma_e2": .55,
        "qc_policy": {"policy_name": POLICY, "qC_policy_version": QC_POLICY_VERSION,
                      "ad_shared_eligibility_policy": AD_QC_POLICY_VERSION,
                      "threshold_source": "E2.7_BGEN_DOM_V1_numerical_rules_shared_for_AD_research",
                      "new_ad_policy": AD_QC_POLICY_VERSION},
         "predeclared_thresholds": {
            "strict_type1": [.025, .075], "conditional_type1": [.025, .075],
            "lambda_gc": [.85, 1.15], "null_convergence_rate": 1.0,
            "max_kkt_residual": 1e-7,
            "cross_format_kernel_tolerance": 2e-6,
            "kernel_parity_scope": "BED/PGEN/BGEN one-hot only; uncertain BGEN kernels are expected to differ from BED hardcalls"},
        "direct_kernel_errors": direct_errors,
        "strict_calibration": {f"{k[0]}|{k[1]}": v for k, v in calibration_summary.items()},
        "conditional_partial": {f"{k[0]}|{k[1]}": v for k, v in conditional_summary.items()},
        "alternative_summary": alternative_summary,
        "strict_pass": strict_pass,
        "partial_pass": partial_pass,
        "strict_details": strict_details,
        "strict_loco2_diagnostic": diagnostic_details,
        "partial_details": partial_details,
        "task_runtimes": runtimes,
        "production_fitted_null_scale": {
            "n": args.n, "m": args.m,
            "strict_calibration_replicates": args.calibration_reps},
        "verdict": "PENDING"}
    if strict_pass and partial_pass and all(
            direct_errors[quality]["ka_vs_bed"] < 2e-6 and direct_errors[quality]["kd_vs_bed"] < 2e-6
            for quality in ("pgen_hardcall", "bgen_one_hot")):
        gate["verdict"] = "PRODUCTION READY WITH DOCUMENTED LIMITATIONS"
    else:
        gate["verdict"] = "BLOCKED"
    (run / "final_results.json").write_text(json.dumps(gate, indent=2) + "\n")
    print(json.dumps({"verdict": gate["verdict"], "run_directory": str(run),
                      "strict_pass": strict_pass, "partial_pass": partial_pass}, indent=2))
    require(gate["verdict"] != "BLOCKED", "R31 BGEN AD final gate blocked")


if __name__ == "__main__":
    main()
