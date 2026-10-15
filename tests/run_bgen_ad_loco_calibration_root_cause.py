#!/usr/bin/env python3
"""E3.2c BGEN AD Dual-LOCO statistical inflation root-cause validation."""
import argparse
import csv
import json
import math
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import stats

from run_bgen_ad_null_validation import phenotype_matrix, write_hardcalls
from run_bgen_dominance_production_gate import require


def read_meta(path):
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["Item"]] = row["Value"]
    return result


def rows(path, target="1"):
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if target is None or row["CHR"] == target:
                result[row["SNP"]] = row
    return result


def internal_rows(path, target="1"):
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if target is None or row["chr"] == target:
                result.setdefault(int(row["replicate"]), {})[row["snp"]] = row
    return result


def binmatrix(path, n):
    return np.fromfile(str(path), dtype="<f8").reshape(n, n)


def phenotype_values(probabilities, missing):
    observed = ~missing
    count = np.maximum(observed.sum(axis=1), 1)
    dosage = probabilities[:, :, 1] + 2.0 * probabilities[:, :, 2]
    het = probabilities[:, :, 1]
    return (np.where(observed, dosage - (dosage * observed).sum(axis=1)[:, None] / count[:, None], 0.0),
            np.where(observed, het - (het * observed).sum(axis=1)[:, None] / count[:, None], 0.0))


def probabilities_from_latent(latent, missing):
    result = np.zeros((*latent.shape, 3))
    result[np.arange(latent.shape[0])[:, None], np.arange(latent.shape[1]), np.clip(latent, 0, 2)] = 1.0
    return result


def oracle(ka, kd, probabilities, missing, covariates, values, a2, d2, e2):
    za, zd = phenotype_values(probabilities, missing)
    n = values.shape[0]
    covariance = a2 * ka + d2 * kd + e2 * np.eye(n)
    inverse = np.linalg.inv(covariance)
    c = np.asarray(covariates, dtype=float)
    projection = inverse - inverse @ c @ np.linalg.solve(c.T @ inverse @ c, c.T @ inverse)
    pa, pd = projection @ za.T, projection @ zd.T
    aa = (za.T * pa).sum(axis=0)
    dd = (zd.T * pd).sum(axis=0)
    ad = (za.T * pd).sum(axis=0)
    sa = za @ projection @ values
    sd = zd @ projection @ values
    output = {}
    for j in range(za.shape[0]):
        info = np.array([[aa[j], ad[j]], [ad[j], dd[j]]])
        try:
            beta = np.linalg.solve(info, np.array([sa[j], sd[j]]))
            variances = np.diag(np.linalg.inv(info))
            if np.any(variances < 0):
                continue
            joint = float(beta @ info @ beta)
            output[f"rs{j+1:06d}"] = {
                "beta_a": float(beta[0]), "beta_d": float(beta[1]),
                "se_a": float(math.sqrt(variances[0])), "se_d": float(math.sqrt(variances[1])),
                "p_add": float(stats.chi2.sf(beta[0] * beta[0] / variances[0], 1)),
                "p_dom": float(stats.chi2.sf(beta[1] * beta[1] / variances[1], 1)),
                "p_joint": float(stats.chi2.sf(joint, 2))}
        except np.linalg.LinAlgError:
            continue
    return output


def metrics(reps_rows, key="P_JOINT", df=2):
    types, lambdas = [], []
    for rep_id, rep_rows in sorted(reps_rows.items()):
        p = np.array([float(row[key]) for row in rep_rows.values()])
        if len(p) == 0:
            continue
        statistics = stats.chi2.isf(np.clip(p, 1e-300, 1.0), df)
        lambdas.append(float(np.median(statistics) / stats.chi2.ppf(.5, df)))
        types.append(float(np.mean(p < .05)))
    bootstrap = []
    generator = np.random.default_rng(2718 + (0 if df == 2 else 1))
    array_type = np.asarray(types)
    array_lambda = np.asarray(lambdas)
    for _ in range(2000):
        index = generator.integers(0, len(array_type), len(array_type))
        bootstrap.append((float(np.mean(array_type[index])), float(np.median(array_lambda[index]))))
    low, high = np.percentile([b[0] for b in bootstrap], [2.5, 97.5])
    llow, lhigh = np.percentile([b[1] for b in bootstrap], [2.5, 97.5])
    return {"replicates": len(types), "type1_mean": float(np.mean(types)),
            "type1_ci95": [float(low), float(high)],
            "lambda_gc_median": float(np.median(lambdas)), "lambda_gc_median_ci95": [float(llow), float(lhigh)],
            "lambda_gc_mean": float(np.mean(lambdas))}


def engine_rows(path, target="1"):
    result = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if target is None or row["chr"] == target:
                result[row["snp"]] = row
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--harness", required=True)
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--null-driver", required=True)
    parser.add_argument("--calibration-driver", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--crossformat-reps", type=int, default=50)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--plink2", default=os.environ.get("PLINK2_BIN", "/disk228/huangjie/software/miniconda3/envs/benchmark/bin/plink2"))
    args = parser.parse_args()
    if args.smoke:
        args.reps = min(args.reps, 4)
        args.crossformat_reps = min(args.crossformat_reps, 2)
    root = Path(args.workdir)
    root.mkdir(parents=True, exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="e32c-", dir=root))
    (root / "latest_run.txt").write_text(str(run) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def execute(command):
        proc = subprocess.run([str(v) for v in command], env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        require(proc.returncode == 0, proc.stdout[-8000:])
        return proc

    n, m = 160, 120
    rng = np.random.default_rng(3264)
    frequencies = rng.uniform(.15, .45, m)
    latent = rng.binomial(2, frequencies[:, None], (m, n))
    missing = rng.random((m, n)) < .02
    latent[0] = 0
    latent[1] = 1
    missing[1] = False
    missing[2] = True
    samples = [f"R{i:04d}" for i in range(n)]
    source_prefix = run / "hardcall"
    write_hardcalls(source_prefix, np.where(missing, -1, latent), samples)
    execute([args.plink2, "--bfile", source_prefix, "--make-pgen", "--out", str(source_prefix) + "_pgen"])
    execute([args.plink2, "--bfile", source_prefix, "--export", "bgen-1.2", "--out", str(source_prefix) + "_bgen"])
    bgen = Path(str(source_prefix) + "_bgen.bgen")
    pgen = Path(str(source_prefix) + "_pgen.pgen")
    probabilities = probabilities_from_latent(latent, missing)
    covariates = np.column_stack([np.ones(n), rng.normal(size=n)]).astype(np.float32).astype(float)
    cov = run / "covariates.txt"
    phenotype_matrix(cov, covariates)
    native_cov = run / "covariates_native.txt"
    with native_cov.open("w") as handle:
        for sample, value in zip(samples, covariates[:, 1]):
            handle.write(f"{sample} {sample} {value:.17g}\n")

    dummy = run / "dummy.pheno"
    phenotype_matrix(dummy, np.ones((n, 1)))
    execute([args.calibration_driver, bgen, dummy, cov, run / "global_pair", "1", "-", "fixed", ".25", ".20", ".55"])
    execute([args.calibration_driver, bgen, dummy, cov, run / "outside1_pair", "1", "1", "fixed", ".25", ".20", ".55"])
    ka_g = binmatrix(run / "global_pair.ka.bin", n)
    kd_g = binmatrix(run / "global_pair.kd.bin", n)
    ka_o = binmatrix(run / "outside1_pair.ka.bin", n)
    kd_o = binmatrix(run / "outside1_pair.kd.bin", n)

    noise = rng.standard_normal((n, args.reps))
    noise_partial = rng.standard_normal((n, args.reps))
    matched = np.linalg.cholesky(.25 * ka_o + .20 * kd_o + .55 * np.eye(n)) @ noise + .4 * covariates[:, 1, None]
    shared_stress = np.linalg.cholesky(.25 * ka_g + .20 * kd_g + .55 * np.eye(n)) @ noise_partial + .4 * covariates[:, 1, None]
    global_matched = np.linalg.cholesky(.25 * ka_g + .20 * kd_g + .55 * np.eye(n)) @ rng.standard_normal((n, args.reps)) + .4 * covariates[:, 1, None]
    pure_a = np.linalg.cholesky(.25 * ka_o + .55 * np.eye(n)) @ rng.standard_normal((n, args.reps)) + .4 * covariates[:, 1, None]
    pure_d = np.linalg.cholesky(.20 * kd_o + .55 * np.eye(n)) @ rng.standard_normal((n, args.reps)) + .4 * covariates[:, 1, None]
    for label, values in (("matched_outside", matched), ("full_polygenic", shared_stress),
                          ("global_matched", global_matched), ("pure_additive", pure_a),
                          ("pure_dominance", pure_d)):
        phenotype_matrix(run / f"{label}.pheno", values)

    for label, values, scan, excluded, a, d, e in (
            ("b1_matched", matched, "1", "1", .25, .20, .55),
            ("b1_full_mismatch", shared_stress, "1", "1", .25, .20, .55),
            ("b1_global", global_matched, "1", "-", .25, .20, .55),
            ("b1_pure_a", pure_a, "1", "1", .25, 0.0, .55),
            ("b1_pure_d", pure_d, "1", "1", 0.0, .20, .55)):
        phenotype_matrix(run / f"{label}.input.pheno", values)
        output = run / label
        execute([args.calibration_driver, bgen, run / f"{label}.input.pheno", cov, output, scan, excluded,
                 "fixed", str(a), str(d), str(e)])

    knownv = {}
    oracle_by_path = {}
    for scenario_prefix, ka, kd, a, d, e in (
            ("b1_matched", ka_o, kd_o, .25, .20, .55),
            ("b1_full_mismatch", ka_o, kd_o, .25, .20, .55),
            ("b1_global", ka_g, kd_g, .25, .20, .55),
            ("b1_pure_a", ka_o, kd_o, .25, 0.0, .55),
            ("b1_pure_d", ka_o, kd_o, 0.0, .20, .55)):
        reps = {}
        for rep in range(args.reps):
            path = run / f"{scenario_prefix}.input.pheno"
            v = np.loadtxt(path, skiprows=1)[:, rep]
            reps[rep] = oracle(ka, kd, probabilities, missing, covariates, v, a, d, e)
        knownv[scenario_prefix] = metrics(reps, "p_joint", 2)
        oracle_by_path[scenario_prefix] = reps
        engine = internal_rows(run / f"{scenario_prefix}.internal.tsv")
        knownv[scenario_prefix]["engine_max_beta_error"] = 0.0
        knownv[scenario_prefix]["engine_max_se_error"] = 0.0
        knownv[scenario_prefix]["engine_max_p_error"] = 0.0
        for rep, reprows in engine.items():
            for snp, row in reprows.items():
                if snp not in reps[rep]:
                    continue
                reference = reps[rep][snp]
                knownv[scenario_prefix]["engine_max_beta_error"] = max(knownv[scenario_prefix]["engine_max_beta_error"],
                                                             abs(float(row["beta_a"]) - reference["beta_a"]),
                                                             abs(float(row["beta_d"]) - reference["beta_d"]))
                knownv[scenario_prefix]["engine_max_se_error"] = max(knownv[scenario_prefix]["engine_max_se_error"],
                                                           abs(float(row["se_a"]) - reference["se_a"]),
                                                           abs(float(row["se_d"]) - reference["se_d"]))
                knownv[scenario_prefix]["engine_max_p_error"] = max(knownv[scenario_prefix]["engine_max_p_error"],
                                                          abs(float(row["p_add"]) - reference["p_add"]),
                                                          abs(float(row["p_dom"]) - reference["p_dom"]),
                                                          abs(float(row["p_joint"]) - reference["p_joint"]))

    def harness(label, values, chromosome, loco):
        out = run / label
        out.mkdir()
        cfg = run / f"{label}.options"
        cfg.write_text("\n".join([f"maxit\t200", "tol\t1e-5", f"chromosome\t{chromosome}",
                                  f"loco\t{'true' if loco else 'false'}"]) + "\n")
        phenotype_matrix(run / f"{label}.input.tsv", values)
        execute([args.harness, bgen, run / f"{label}.input.tsv", cov, out, cfg])
        return out

    paths = {
        "matched": harness("bgen_matched", matched, "1", True),
        "global": harness("bgen_global_matched", global_matched, "1", False),
        "stress": harness("bgen_full_polygenic_loco", shared_stress, "1", True),
        "pure_a": harness("bgen_pure_additive", pure_a, "1", True),
        "pure_d": harness("bgen_pure_dominance", pure_d, "1", True)}

    fitted = {}
    comparisons = {}
    for key, directory in paths.items():
        reps_rows = {}
        for rep in range(args.reps):
            reps_rows[rep] = rows(directory / f"rep{rep}.assoc.ad.txt")
        fitted[key] = metrics(reps_rows)
        directory_keys = {"matched": "b1_matched", "full_mismatch": "b1_full_mismatch",
                          "global": "b1_global", "pure_a": "b1_pure_a", "pure_d": "b1_pure_d"}
        corresponding = {"matched": "b1_matched", "stress": "b1_full_mismatch", "global": "b1_global",
                         "pure_a": "b1_pure_a", "pure_d": "b1_pure_d"}
        if key in corresponding:
            reps = internal_rows(run / f"{corresponding[key]}.internal.tsv")
            max_beta, max_se, max_p = 0.0, 0.0, 0.0
            for rep, row in reps_rows.items():
                engine = reps.get(rep, {})
                for snp, native in row.items():
                    if snp not in engine:
                        continue
                    source = engine[snp]
                    max_beta = max(max_beta, abs(float(native["BETA_A"]) - float(source["beta_a"])),
                                   abs(float(native["BETA_D"]) - float(source["beta_d"])))
                    max_se = max(max_se, abs(float(native["SE_A"]) - float(source["se_a"])),
                                 abs(float(native["SE_D"]) - float(source["se_d"])))
                    max_p = max(max_p, abs(float(native["P_A"]) - float(source["p_add"])),
                                abs(float(native["P_D"]) - float(source["p_dom"])),
                                abs(float(native["P_JOINT"]) - float(source["p_joint"])))
            comparisons[key] = {"b3_vs_b2_max_beta": max_beta, "b3_vs_b2_max_se": max_se,
                                "b3_vs_b2_max_joint_p": max_p}
            baseline = oracle_by_path.get(corresponding[key])
            if baseline is not None:
                max_beta, max_se, max_p = 0.0, 0.0, 0.0
                for rep, row in reps_rows.items():
                    for snp, native in row.items():
                        if snp in baseline.get(rep, {}):
                            source = baseline[rep][snp]
                            max_beta = max(max_beta, abs(float(native["BETA_A"]) - source["beta_a"]),
                                           abs(float(native["BETA_D"]) - source["beta_d"]))
                            max_se = max(max_se, abs(float(native["SE_A"]) - source["se_a"]),
                                         abs(float(native["SE_D"]) - source["se_d"]))
                            max_p = max(max_p, abs(float(native["P_A"]) - source["p_add"]),
                                        abs(float(native["P_D"]) - source["p_dom"]),
                                        abs(float(native["P_JOINT"]) - source["p_joint"]))
                comparisons[key].update({"b3_vs_b1_max_beta": max_beta, "b3_vs_b1_max_se": max_se,
                                         "b3_vs_b1_max_joint_p": max_p})

    conditional = {}
    for key, directory, alpha_key in (("pure_additive_dom", paths["pure_a"], "P_D"),
                                      ("pure_dominance_add", paths["pure_d"], "P_A")):
        reps_rows = {}
        for rep in range(args.reps):
            reps_rows[rep] = rows(directory / f"rep{rep}.assoc.ad.txt")
        conditional[key] = metrics(reps_rows, alpha_key, 1)
    null_conditional = {}
    for kind, label in (("joint", "P_JOINT"), ("dom_given_a", "P_D"), ("add_given_d", "P_A")):
        reps_rows = {}
        for rep in range(args.reps):
            reps_rows[rep] = rows(paths["pure_a"] / f"rep{rep}.assoc.ad.txt")
        if kind == "joint":
            null_conditional[kind] = metrics(reps_rows, label, 2)
    for key, rows_path, metric_label, df in (
            ("additive_joint", paths["pure_a"], "P_JOINT", 2),
            ("additive_p_dom", paths["pure_a"], "P_D", 1),
            ("dominance_joint", paths["pure_d"], "P_JOINT", 2),
            ("dominance_p_add", paths["pure_d"], "P_A", 1)):
        reps_rows = {}
        for rep in range(args.reps):
            reps_rows[rep] = rows(rows_path / f"rep{rep}.assoc.ad.txt")
        conditional[key] = metrics(reps_rows, metric_label, df)

    # Native BED/PGEN/BGEN under the same full-genome stress. BGEN already uses the internal path above.
    stress_reps = min(args.crossformat_reps, args.reps)
    directories = {"bed": [], "pgen": []}
    bed_ci = {}
    pgen_ci = {}
    bed_rows_all, pgen_rows_all, bgen_rows_all = [], [], []
    for replica in range(stress_reps):
        y = shared_stress[:, replica]
        pheno = run / f"native_stress_{replica}.pheno"
        write_sv(pheno, samples, y)
        tasks = []
        for fmt in ("bed", "pgen"):
            source = source_prefix if fmt == "bed" else str(source_prefix) + "_pgen"
            out = run / f"cross_{fmt}_rep{replica}"
            flag = "--bfile" if fmt == "bed" else "--pfile"
            tasks.append((fmt, out, [args.scanner, flag, source, "--gwas-model",
                                    "additive_dominance", "--pheno", pheno, "--qcovar", native_cov,
                                    "--maf", ".01", "--reml-maxit", "200", "--reml-tol", "1e-5", "--out", out]))
        for fmt, _out, command in tasks:
            execute(command)
        bed_ci[replica] = read_meta(run / f"cross_bed_rep{replica}.scan.meta.txt")
        pgen_ci[replica] = read_meta(run / f"cross_pgen_rep{replica}.scan.meta.txt")
        bed_rows_all.append({replica: rows(run / f"cross_bed_rep{replica}.assoc.ad.txt")})
        pgen_rows_all.append({replica: rows(run / f"cross_pgen_rep{replica}.assoc.ad.txt")})
        bgen_rows_all.append({replica: rows(paths["stress"] / f"rep{replica}.assoc.ad.txt")})

    def aggregate_rows(chunks):
        output = {}
        for chunk in chunks:
            replica, value = next(iter(chunk.items()))
            output[replica] = value
        return output
    bgen_target = aggregate_rows(bgen_rows_all)
    bed_target = aggregate_rows(bed_rows_all)
    pgen_target = aggregate_rows(pgen_rows_all)
    cross = {}
    for name, target in (("BED", bed_target), ("PGEN", pgen_target), ("BGEN", bgen_target)):
        cross[name] = metrics(target)
    cross_row_differences = {}

    def row_statistic_differences(left, right):
        maxima = {field: 0.0 for field in ("BETA_A", "BETA_D", "SE_A", "SE_D", "P_A", "P_D", "P_JOINT")}
        for rep, left_rows in left.items():
            right_rows = right.get(rep, {})
            for snp, row in left_rows.items():
                if snp not in right_rows:
                    continue
                other = right_rows[snp]
                for field in maxima:
                    maxima[field] = max(maxima[field], abs(float(row[field]) - float(other[field])))
        return maxima

    for pair, left, right in (("pgen_vs_bed", bed_target, pgen_target),
                              ("bgen_vs_bed", bed_target, bgen_target),
                              ("bgen_vs_pgen", pgen_target, bgen_target)):
        cross_row_differences[pair] = row_statistic_differences(left, right)
    native_bgen_vc = {name: metrics(bgen_target[i:i+1], "P_JOINT", 2) for name, i in (("unused", 0),)} if False else None
    def target_vcs(meta):
        if "SIGMA_A2_1" in meta:
            return (float(meta["SIGMA_A2_1"]), float(meta["SIGMA_D2_1"]), float(meta["SIGMA_E2_1"]))
        return (float(meta["SIGMA_A2"]), float(meta["SIGMA_D2"]), float(meta["SIGMA_E2"]))

    native_vc = []
    for rep in range(stress_reps):
        for name, directory, meta in (
                ("BED", "cross_bed_rep" + str(rep), bed_ci[rep]),
                ("PGEN", "cross_pgen_rep" + str(rep), pgen_ci[rep])):
            a2, d2, e2 = target_vcs(meta)
            native_vc.append({"format": name, "rep": rep,
                              "sigma_a2_1": a2, "sigma_d2_1": d2, "sigma_e2_1": e2})
    bgen_vc = []
    for rep in range(stress_reps):
        meta = read_meta(paths["stress"] / f"rep{rep}.scan.meta.txt")
        a2, d2, e2 = target_vcs(meta)
        bgen_vc.append({"format": "BGEN", "rep": rep,
                        "sigma_a2_1": a2, "sigma_d2_1": d2, "sigma_e2_1": e2})
    cross_vc_mean = {}
    for c in native_vc + bgen_vc:
        cross_vc_mean.setdefault(c["format"], {"A2": [], "D2": [], "E2": []})
        cross_vc_mean[c["format"]]["A2"].append(c["sigma_a2_1"])
        cross_vc_mean[c["format"]]["D2"].append(c["sigma_d2_1"])
        cross_vc_mean[c["format"]]["E2"].append(c["sigma_e2_1"])
    for name, values in cross_vc_mean.items():
        for key in values:
            values[key] = float(np.mean(values[key]))

    # Direct one-hot kernel parity on the chromosome-2 (outside-target) partition.
    for fmt, source in (("bed", source_prefix), ("pgen", pgen), ("bgen", bgen)):
        source_value = source_prefix if fmt == "bed" else (pgen if fmt == "pgen" else bgen)
        execute([args.null_driver, fmt, source_value, run / f"direct_{fmt}", "-", "2", ".01", "-", "0"])
    direct_errors = {}
    for key, left in (("PGEN_vs_BED", "direct_pgen"), ("BGEN_vs_BED", "direct_bgen"), ("BGEN_vs_PGEN", "direct_bgen")):
        if key.startswith("BGEN_vs_B"):
            pass
    for fmt, other in (("pgen", "bed"), ("bgen", "bed"), ("bgen", "pgen")):
        left_k = binmatrix(run / f"direct_{fmt}.ka.bin", n)
        right_k = binmatrix(run / f"direct_{other}.ka.bin", n)
        left_d = binmatrix(run / f"direct_{fmt}.kd.bin", n)
        right_d = binmatrix(run / f"direct_{other}.kd.bin", n)
        direct_errors[f"{fmt}_vs_{other}_ka_direct_chrom2"] = float(np.max(np.abs(left_k - right_k)))
        direct_errors[f"{fmt}_vs_{other}_kd_direct_chrom2"] = float(np.max(np.abs(left_d - right_d)))

    fixed_engine = {key: internal_rows(run / f"{key}.internal.tsv") for key in (
        "b1_matched", "b1_full_mismatch", "b1_global", "b1_pure_a", "b1_pure_d")}
    summary_rows = []
    for key, path_rows in fixed_engine.items():
        stats = knownv[key]
        summary_rows.append({"path": key, "synthetic_replicates": args.reps,
                             "joint_type1": stats["type1_mean"], "joint_lambda_median": stats["lambda_gc_median"],
                             "engine_max_beta": stats["engine_max_beta_error"],
                             "engine_max_se": stats["engine_max_se_error"],
                             "engine_max_p": stats["engine_max_p_error"]})
    summary = {
        "verdict_predeclared_thresholds": {
            "marginal_type1": [0.025, 0.075], "lambda_gc_median": [0.85, 1.15],
            "known_V_engine_beta_se_tolerance": 1e-8, "cross_format_kernel_tolerance": 2e-6},
        "fitted_calibration": fitted,
        "conditional_null_calibration": conditional,
        "known_v_engine": summary_rows,
        "cross_format": cross,
        "cross_format_row_differences": cross_row_differences,
        "cross_format_mean_VCs_target1": cross_vc_mean,
        "direct_kernel_errors": direct_errors,
        "b1_b3_comparisons_beta_SE_P": comparisons,
        "bgen_harness_indices": {key: str(path) for key, path in paths.items()},
        "run_directory": str(run)}
    (run / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    calibration_rows = []
    for key, value in fitted.items():
        calibration_rows.append({"scenario": key, **value})
    if calibration_rows:
        keys = list(calibration_rows[0])
        with (run / "fitted_calibration.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
            writer.writeheader()
            writer.writerows(calibration_rows)
    print(json.dumps(summary, indent=2))
    require(knownv["b1_matched"]["engine_max_p_error"] < 1e-8, "known-V engine/oracle parity")
    require(all(knownv[k]["engine_max_p_error"] < 1e-8 for k in knownv), "known-V engine/oracle parity")


def metrics_from_oracle(reps, key):
    types, lambdas = [], []
    for rep_rows in reps.values():
        p = np.asarray([row[key] for row in rep_rows.values()])
        if len(p) == 0:
            continue
        q = stats.chi2.isf(np.clip(p, 1e-300, 1.0), 2 if key == "p_joint" else 1)
        lambdas.append(float(np.median(q) / stats.chi2.ppf(.5, 2 if key == "p_joint" else 1)))
        types.append(float(np.mean(p < .05)))
    return {"type1_mean": float(np.mean(types)), "lambda_gc_median": float(np.median(lambdas))} if types else {"type1_mean": None, "lambda_gc_median": None}


def write_sv(path, samples, values):
    with Path(path).open("w") as handle:
        for sample, value in zip(samples, values):
            handle.write(f"{sample} {sample} {value:.17g}\n")


if False:
    def metrics(ValueKnownComparison):  # pragma: no cover
        pass


if __name__ == "__main__":
    main()
