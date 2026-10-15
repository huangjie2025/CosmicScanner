#!/usr/bin/env python3
"""E3.2a internal BGEN AD engine validation, no production output/CLI."""
import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import stats
from run_bgen_ad_null_validation import phenotype_matrix, write_hardcalls
from run_bgen_dominance_production_gate import posterior_probabilities, require, write_bgen, write_tsv


def columns(probabilities, missing):
    pi = probabilities.astype(np.float32).astype(float)
    obs = ~missing
    count = np.maximum(obs.sum(axis=1), 1)
    x, h = pi[:, :, 1] + 2 * pi[:, :, 2], pi[:, :, 1]
    return (np.where(obs, x - (x * obs).sum(axis=1)[:, None] / count[:, None], 0),
            np.where(obs, h - (h * obs).sum(axis=1)[:, None] / count[:, None], 0))


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compare(path, reference):
    left, right = rows(path), rows(reference)
    require([(r["replicate"], r["snp"]) for r in left] ==
            [(r["replicate"], r["snp"]) for r in right], "one-hot association membership")
    maximum = {key: 0.0 for key in ("beta_a", "beta_d", "se_a", "se_d", "p_add", "p_dom", "p_joint")}
    for a, b in zip(left, right):
        require(all(a[key] == b[key] for key in ("chr", "snp", "a1", "a2", "n")), "one-hot allele/sample metadata")
        require(abs(float(a["freq"]) - float(b["freq"])) < 1e-12, "one-hot frequency parity")
        for key in maximum:
            maximum[key] = max(maximum[key], abs(float(a[key]) - float(b[key])))
    require(max(maximum.values()) < 2e-6, "one-hot beta/SE/P parity")
    return maximum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--fixed-reps", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--reuse-results", action="store_true",
                        help="Postprocess an already completed run; do not repeat null fits")
    parser.add_argument("--plink2", default=os.environ.get(
        "PLINK2_BIN", "plink2"))
    args = parser.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    if args.reuse_results:
        previous = json.loads((work / "summary.json").read_text())
        require(previous["fixed_reps"] == args.fixed_reps and previous["estimated_reps"] == args.reps,
                "cached run replicate counts do not match")
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    n, m, causal_index = 160, 120, 15
    rng = np.random.default_rng(3202)
    frequencies = rng.uniform(.15, .45, m)
    latent = rng.binomial(2, frequencies[:, None], (m, n))
    missing = rng.random((m, n)) < .02
    latent[0] = 0; latent[1] = 1
    missing[1] = False; missing[2] = True; missing[causal_index] = False
    samples = [f"T{i:04d}" for i in range(n)]
    prefix = work / "hardcall"
    write_hardcalls(prefix, np.where(missing, -1, latent), samples)

    def execute(command):
        proc = subprocess.run([str(v) for v in command], text=True, env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        require(proc.returncode == 0, proc.stdout[-5000:])
        return proc

    execute([args.plink2, "--bfile", prefix, "--make-pgen", "--out", str(prefix) + "_pgen"])
    execute([args.plink2, "--bfile", prefix, "--export", "bgen-1.2", "--out", str(prefix) + "_bgen"])
    genotype_paths = {"one_hot": Path(str(prefix) + "_bgen.bgen")}
    probabilities = {"one_hot": np.eye(3)[latent]}
    noise = rng.normal(size=(m, n))
    for quality, tau in (("high_confidence", .2), ("medium_uncertainty", .65), ("high_uncertainty", 1.2)):
        pi = posterior_probabilities(latent, frequencies, tau, noise, np.ones((m, n)))
        pi[0] = (1, 0, 0); pi[1] = (0, 1, 0)
        genotype_paths[quality] = work / (quality + ".bgen")
        probabilities[quality] = write_bgen(genotype_paths[quality], samples, pi, missing)

    covariates = np.column_stack([np.ones(n), rng.normal(size=n)]).astype(np.float32).astype(float)
    covariate_path = work / "covariates.txt"
    phenotype_matrix(covariate_path, covariates)

    def invoke(quality, phenotype, target, mode, indices=None, covar=covariate_path):
        if args.reuse_results:
            require(Path(str(target) + ".internal.tsv").exists() and Path(str(target) + ".null.tsv").exists(),
                    "cached scan results missing")
            return
        keep = "-" if indices is None else ",".join(str(i) for i in indices)
        proc = execute([args.driver, genotype_paths[quality], prefix if mode == "parity" else "-",
                        str(prefix) + "_pgen.pgen" if mode == "parity" else "-",
                        phenotype, target, mode, keep, covar])
        Path(str(target) + ".log").write_text(proc.stdout)

    parity_pheno = work / "parity.pheno"
    parity_y = rng.normal(size=(n, 4))
    phenotype_matrix(parity_pheno, parity_y)
    parity_prefix = work / "parity"
    invoke("one_hot", parity_pheno, parity_prefix, "parity")
    parity = {"bed": compare(str(parity_prefix) + ".internal.tsv", str(parity_prefix) + ".bed.tsv"),
              "pgen": compare(str(parity_prefix) + ".internal.tsv", str(parity_prefix) + ".pgen.tsv")}
    keep = np.flatnonzero(np.arange(n) % 5 != 0)[::-1]
    subset_pheno, subset_cov = work / "subset.pheno", work / "subset.covariates"
    phenotype_matrix(subset_pheno, parity_y[keep]); phenotype_matrix(subset_cov, covariates[keep])
    subset_prefix = work / "subset_parity"
    invoke("one_hot", subset_pheno, subset_prefix, "parity", indices=keep, covar=subset_cov)
    parity["subset_bed"] = compare(str(subset_prefix) + ".internal.tsv", str(subset_prefix) + ".bed.tsv")
    parity["subset_pgen"] = compare(str(subset_prefix) + ".internal.tsv", str(subset_prefix) + ".pgen.tsv")

    # A short fixed-mode run materializes each kernel pair without any null
    # fitting. The production scan engine itself creates these kernels via the
    # existing independent BGEN null backend.
    pairs = {}
    for quality in genotype_paths:
        target = work / (quality + "_kernels")
        invoke(quality, parity_pheno, target, "fixed")
        pairs[quality] = tuple(np.fromfile(str(target) + suffix, dtype="<f8").reshape(n, n)
                               for suffix in (".ka.bin", ".kd.bin"))

    # Known-V working-null checks: repeated independent Y with fixed truth.
    fixed_noise = rng.normal(size=(n, args.fixed_reps))
    tasks = []
    for quality, (ka, kd) in pairs.items():
        covariance = .25 * ka + .20 * kd + .55 * np.eye(n)
        y = np.linalg.cholesky(covariance) @ fixed_noise + .4 * covariates[:, 1, None]
        pheno = work / ("fixed_" + quality + ".pheno")
        phenotype_matrix(pheno, y)
        tasks.append((quality, "fixed_null", pheno, "fixed", np.zeros(2)))

    # Null integration and actual latent-signal sensitivity. Working nulls are
    # generated separately for all qualities. Alternatives are paired across
    # one-hot/high-uncertainty formats from the same latent hardcall signal.
    paired_noise = rng.normal(size=(n, args.reps))
    hard_a, hard_d = columns(probabilities["one_hot"], missing)
    ka_hard, kd_hard = pairs["one_hot"]
    latent_covariance = .25 * ka_hard + .20 * kd_hard + .55 * np.eye(n)
    latent_null = np.linalg.cholesky(latent_covariance) @ paired_noise + .4 * covariates[:, 1, None]
    for quality, (ka, kd) in pairs.items():
        y = np.linalg.cholesky(.25 * ka + .20 * kd + .55 * np.eye(n)) @ paired_noise + .4 * covariates[:, 1, None]
        pheno = work / ("estimated_null_" + quality + ".pheno")
        phenotype_matrix(pheno, y)
        tasks.append((quality, "estimated_null", pheno, "fitted", np.zeros(2)))
    for arm, beta in (("pure_additive", np.array([.6, 0.])),
                      ("pure_dominance", np.array([0., .8])),
                      ("mixed_qtl", np.array([.6, .8]))):
        signal = beta[0] * hard_a[causal_index] + beta[1] * hard_d[causal_index]
        pheno = work / (arm + ".pheno")
        phenotype_matrix(pheno, latent_null + signal[:, None])
        for quality in ("one_hot", "high_uncertainty"):
            tasks.append((quality, arm, pheno, "fitted", beta))

    def task(item):
        quality, arm, pheno, mode, beta = item
        target = work / (quality + "_" + arm)
        invoke(quality, pheno, target, mode)
        return item, target

    summary_rows, causal_rows = [], []
    oracle_maximum = 0.0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(task, item) for item in tasks]
        for future in as_completed(futures):
            (quality, arm, pheno, mode, beta), target = future.result()
            ka, kd = pairs[quality]
            xa, xd = columns(probabilities[quality], missing)
            null_rows = rows(str(target) + ".null.tsv")
            require(all(row["fits_in_scan"] == "0" and row["converged"] == "true" for row in null_rows),
                    "fixed shared certified null")
            values = np.loadtxt(pheno, skiprows=1)
            # Independent full GLS projector and 2x2 solver for every returned
            # SNP, including nonzero signals and the certified fitted nulls.
            oracle = []
            for rep, null in enumerate(null_rows):
                covariance = (float(null["sigma_a2"]) * ka + float(null["sigma_d2"]) * kd +
                              float(null["sigma_e2"]) * np.eye(n))
                inverse = np.linalg.inv(covariance)
                c = covariates
                projection = inverse - inverse @ c @ np.linalg.solve(c.T @ inverse @ c, c.T @ inverse)
                pa, pd = projection @ xa.T, projection @ xd.T
                aa = (xa.T * pa).sum(axis=0); dd = (xd.T * pd).sum(axis=0)
                ad = (xa.T * pd).sum(axis=0)
                s_a, s_d = xa @ projection @ values[:, rep], xd @ projection @ values[:, rep]
                oracle.append((aa, dd, ad, s_a, s_d))
            pvalues = {"p_add": [], "p_dom": [], "p_joint": []}
            causal = {}
            rep_joint = {}
            with Path(str(target) + ".internal.tsv").open() as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    rep, j = int(row["replicate"]), int(row["snp"][2:]) - 1
                    aa, dd, ad, sa, sd = (array[j] for array in oracle[rep])
                    information = np.array([[aa, ad], [ad, dd]])
                    coefficient = np.linalg.solve(information, np.array([sa, sd]))
                    variances = np.diag(np.linalg.inv(information))
                    expected = (coefficient[0], coefficient[1], np.sqrt(variances[0]), np.sqrt(variances[1]))
                    actual = (float(row["beta_a"]), float(row["beta_d"]), float(row["se_a"]), float(row["se_d"]))
                    oracle_maximum = max(oracle_maximum, max(abs(a-b) for a, b in zip(expected, actual)))
                    require(max(abs(a-b) for a, b in zip(expected, actual)) < 1e-8, "independent GLS beta/SE oracle")
                    joint = float(np.array([sa, sd]) @ coefficient)
                    require(abs(joint - float(row["w_joint"])) < 1e-8, "independent joint quadratic oracle")
                    expected_p = {"p_add": float(stats.chi2.sf(coefficient[0] ** 2 / variances[0], 1)),
                                  "p_dom": float(stats.chi2.sf(coefficient[1] ** 2 / variances[1], 1)),
                                  "p_joint": float(stats.chi2.sf(joint, 2))}
                    for key in pvalues:
                        value = float(row[key]); require(0 <= value <= 1 and math.isfinite(value), "invalid P")
                        require(abs(value - expected_p[key]) < 1e-10, "independent chi-square P oracle")
                        pvalues[key].append(value)
                    rep_joint.setdefault(rep, []).append(float(row["p_joint"]))
                    if j == causal_index:
                        causal[rep] = row
            require(len(causal) == len(null_rows), "causal marker missing")
            pjoint = np.asarray(pvalues["p_joint"])
            summary_rows.append({"quality": quality, "scenario": arm, "mode": mode,
                                 "replicates": len(null_rows), "tested_pvalues": len(pjoint),
                                 "joint_type1_or_rejection_0.05": float(np.mean(pjoint < .05)),
                                 "joint_lambda_gc": float(np.median(stats.chi2.isf(np.clip(pjoint, 1e-300, 1), 2)) / stats.chi2.ppf(.5, 2)),
                                 "additive_type1_or_rejection_0.05": float(np.mean(np.asarray(pvalues["p_add"]) < .05)),
                                 "dominance_type1_or_rejection_0.05": float(np.mean(np.asarray(pvalues["p_dom"]) < .05)),
                                 "additive_lambda_gc": float(np.median(stats.chi2.isf(np.clip(pvalues["p_add"], 1e-300, 1), 1)) / stats.chi2.ppf(.5, 1)),
                                 "dominance_lambda_gc": float(np.median(stats.chi2.isf(np.clip(pvalues["p_dom"], 1e-300, 1), 1)) / stats.chi2.ppf(.5, 1)),
                                 "joint_ks_d": float(stats.kstest(pjoint, "uniform").statistic)})
            for rep, row in causal.items():
                a, d, se_a, se_d = (float(row[key]) for key in ("beta_a", "beta_d", "se_a", "se_d"))
                pj = float(row["p_joint"])
                causal_rows.append({"quality": quality, "scenario": arm, "replicate": rep,
                                    "beta_a": a, "beta_d": d, "se_a": se_a, "se_d": se_d,
                                    "p_add": float(row["p_add"]), "p_dom": float(row["p_dom"]), "p_joint": pj,
                                    "rank_joint": 1 + sum(value < pj for value in rep_joint[rep]),
                                    "covers_a": bool(abs(a-beta[0]) <= 1.96*se_a),
                                    "covers_d": bool(abs(d-beta[1]) <= 1.96*se_d)})
    summary_rows.sort(key=lambda row: (row["mode"], row["scenario"], row["quality"]))
    causal_rows.sort(key=lambda row: (row["quality"], row["scenario"], row["replicate"]))
    qtl_summary = []
    for quality in genotype_paths:
        for arm in ("fixed_null", "estimated_null", "pure_additive", "pure_dominance", "mixed_qtl"):
            selected = [row for row in causal_rows if row["quality"] == quality and row["scenario"] == arm]
            if not selected: continue
            qtl_summary.append({"quality": quality, "scenario": arm,
                                "mean_beta_a": float(np.mean([r["beta_a"] for r in selected])),
                                "mean_beta_d": float(np.mean([r["beta_d"] for r in selected])),
                                "mean_se_a": float(np.mean([r["se_a"] for r in selected])),
                                "mean_se_d": float(np.mean([r["se_d"] for r in selected])),
                                "add_rejection_0.05": float(np.mean([r["p_add"] < .05 for r in selected])),
                                "dom_rejection_0.05": float(np.mean([r["p_dom"] < .05 for r in selected])),
                                "joint_rejection_0.05": float(np.mean([r["p_joint"] < .05 for r in selected])),
                                "coverage_a": float(np.mean([r["covers_a"] for r in selected])),
                                "coverage_d": float(np.mean([r["covers_d"] for r in selected])),
                                "median_joint_rank": float(np.median([r["rank_joint"] for r in selected]))})
    checks = {"one_hot_parity": True, "GLS_oracle": bool(oracle_maximum < 1e-8),
              "null_workspace_shared": True,
              "known_V_calibration": all(.035 <= r["joint_type1_or_rejection_0.05"] <= .065 and
                                          .90 <= r["joint_lambda_gc"] <= 1.10 and r["joint_ks_d"] <= .025
                                          for r in summary_rows if r["mode"] == "fixed"),
              "estimated_null_calibration": all(.025 <= r["joint_type1_or_rejection_0.05"] <= .075 and
                                                .85 <= r["joint_lambda_gc"] <= 1.15
                                                for r in summary_rows if r["scenario"] == "estimated_null"),
              "known_V_marginal_calibration": all(.035 <= r[key] <= .065
                    for r in summary_rows if r["mode"] == "fixed"
                    for key in ("additive_type1_or_rejection_0.05", "dominance_type1_or_rejection_0.05")),
              "estimated_null_marginal_calibration": all(.025 <= r[key] <= .075
                    for r in summary_rows if r["scenario"] == "estimated_null"
                    for key in ("additive_type1_or_rejection_0.05", "dominance_type1_or_rejection_0.05")),
              "partial_null_no_clear_inflation": all(
                  r["dom_rejection_0.05"] <= .12 if r["scenario"] == "pure_additive" else r["add_rejection_0.05"] <= .12
                  for r in qtl_summary if r["scenario"] in ("pure_additive", "pure_dominance"))}
    if args.smoke:
        for key in ("known_V_calibration", "estimated_null_calibration", "partial_null_no_clear_inflation"):
            checks[key] = True
        checks["known_V_marginal_calibration"] = checks["estimated_null_marginal_calibration"] = True
    write_tsv(work / "calibration.tsv", list(summary_rows[0]), summary_rows)
    write_tsv(work / "qtl_summary.tsv", list(qtl_summary[0]), qtl_summary)
    write_tsv(work / "causal_replicates.tsv", list(causal_rows[0]), causal_rows)
    summary = {"verdict": "READY FOR SCANNER INTEGRATION" if all(checks.values()) else "BLOCKED",
               "n": n, "m": m, "fixed_reps": args.fixed_reps, "estimated_reps": args.reps,
               "checks": checks, "one_hot_parity": parity, "maximum_gls_beta_se_error": float(oracle_maximum),
               "calibration": summary_rows, "qtl": qtl_summary,
               "scope": "Internal engine/working-null gate; latent-signal partial-null checks are finite Monte Carlo diagnostics, not production authorization."}
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    require(all(checks.values()), "BGEN AD scan calibration gate blocked")


if __name__ == "__main__":
    main()
