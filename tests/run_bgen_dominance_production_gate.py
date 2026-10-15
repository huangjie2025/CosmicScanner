#!/usr/bin/env python3
"""Phase E2.7 production calibration and power gate for BGEN dominance."""

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import stats


POLICY_VERSION = "BGEN_DOM_V1"
CONDITIONS = ("bed", "pgen", "bgen_high", "bgen_medium", "bgen_low")
TAUS = {"bgen_high": 0.35, "bgen_medium": 0.65, "bgen_low": 1.00}
EXPECTED_HEADER = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
                   "BETA_A", "SE_A", "P_A", "LOG10P_A",
                   "BETA_D", "SE_D", "P_D", "LOG10P_D",
                   "P_JOINT", "LOG10P_JOINT"]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def pack_string2(value):
    raw = value.encode("utf-8")
    return struct.pack("<H", len(raw)) + raw


def pack_string4(value):
    raw = value.encode("utf-8")
    return struct.pack("<I", len(raw)) + raw


def generate_hardcalls(n, m, seed):
    rng = np.random.default_rng(seed)
    calls = np.empty((m, n), dtype=np.int8)
    frequencies = np.empty(m, dtype=float)
    previous = None
    for j in range(m):
        if j % 20 == 0:
            p = rng.uniform(0.12, 0.48)
            h1 = rng.binomial(1, p, n)
            h2 = rng.binomial(1, p, n)
        else:
            p = rng.uniform(0.12, 0.48)
            copy1 = rng.random(n) < 0.82
            copy2 = rng.random(n) < 0.82
            h1 = np.where(copy1, previous[0], rng.binomial(1, p, n))
            h2 = np.where(copy2, previous[1], rng.binomial(1, p, n))
        previous = (h1, h2)
        calls[j] = h1 + h2
        frequencies[j] = calls[j].mean() / 2.0
    eligible = np.flatnonzero((frequencies > 0.20) & (frequencies < 0.40))
    require(len(eligible) > 0, "no eligible causal marker")
    causal = int(eligible[len(eligible) // 2])
    return calls, frequencies, causal


def write_bed(prefix, calls, samples):
    m, n = calls.shape
    with Path(str(prefix) + ".fam").open("w") as out:
        for sample in samples:
            out.write(f"{sample} {sample} 0 0 1 -9\n")
    with Path(str(prefix) + ".bim").open("w") as out:
        for j in range(m):
            chrom = 1 if j < m // 2 else 2
            out.write(f"{chrom} rs{j + 1:06d} 0 {(j + 1) * 100} A C\n")
    bed = bytearray(b"\x6c\x1b\x01")
    for row in calls:
        packed = bytearray((n + 3) // 4)
        for i, value in enumerate(row):
            code = (0b11, 0b10, 0b00)[int(value)]
            packed[i // 4] |= code << (2 * (i % 4))
        bed.extend(packed)
    Path(str(prefix) + ".bed").write_bytes(bytes(bed))


def posterior_probabilities(calls, frequencies, tau, noise, quality_factor):
    states = np.arange(3, dtype=float)[:, None]
    m, n = calls.shape
    probs = np.empty((m, n, 3), dtype=np.float64)
    for j in range(m):
        p = frequencies[j]
        prior = np.array([(1.0 - p) ** 2, 2.0 * p * (1.0 - p), p ** 2])[:, None]
        sigma = tau * quality_factor[j][None, :]
        evidence = calls[j][None, :] + sigma * noise[j][None, :]
        weights = prior * np.exp(-0.5 * ((evidence - states) / sigma) ** 2)
        probs[j] = (weights / weights.sum(axis=0)).T
    return probs


def write_bgen(path, samples, probabilities, missing, bits=16):
    m, n, _ = probabilities.shape
    sample_block_len = 8 + sum(2 + len(s.encode("utf-8")) for s in samples)
    header_len = 20
    offset = header_len + sample_block_len
    denominator = (1 << bits) - 1
    flags = 1 | (2 << 2) | (1 << 31)
    decoded = np.zeros_like(probabilities)
    with Path(path).open("wb") as out:
        out.write(struct.pack("<IIII", offset, header_len, m, n))
        out.write(b"bgen")
        out.write(struct.pack("<I", flags))
        out.write(struct.pack("<II", sample_block_len, n))
        for sample in samples:
            out.write(pack_string2(sample))
        for j in range(m):
            variant = f"rs{j + 1:06d}"
            chrom = "1" if j < m // 2 else "2"
            out.write(pack_string2(variant))
            out.write(pack_string2(variant))
            out.write(pack_string2(chrom))
            out.write(struct.pack("<IH", (j + 1) * 100, 2))
            out.write(pack_string4("A"))
            out.write(pack_string4("C"))
            body = bytearray(struct.pack("<IHBB", n, 2, 2, 2))
            body.extend((130 if missing[j, i] else 2) for i in range(n))
            body.extend((0, bits))
            for i in range(n):
                if missing[j, i]:
                    q2 = q1 = 0
                else:
                    p0, p1, p2 = probabilities[j, i]
                    q2 = int(round(p2 * denominator))
                    q1 = int(round(p1 * denominator))
                    if q2 + q1 >= denominator:
                        if q2 >= q1:
                            q2 -= q2 + q1 - denominator + 1
                        else:
                            q1 -= q2 + q1 - denominator + 1
                    decoded[j, i, 2] = q2 / denominator
                    decoded[j, i, 1] = q1 / denominator
                    decoded[j, i, 0] = 1.0 - decoded[j, i, 1] - decoded[j, i, 2]
                body.extend(struct.pack("<HH", q2, q1))
            compressed = zlib.compress(bytes(body), 6)
            out.write(struct.pack("<II", len(compressed) + 4, len(body)))
            out.write(compressed)
    return decoded


def probability_metrics(probs, missing, calls, frequencies):
    m, n, _ = probs.shape
    info = np.zeros(m)
    entropy = np.zeros(m)
    max_prob = np.zeros(m)
    xd_var = np.zeros(m)
    z_columns = []
    denominator = 0.0
    non_one_hot = 0
    observed_total = 0
    for j in range(m):
        keep = ~missing[j]
        p0, p1, p2 = probs[j, keep].T
        p = (p1 + 2.0 * p2).mean() / 2.0
        q = 1.0 - p
        mu = -2.0 * p * p * p0 + 2.0 * p * q * p1 - 2.0 * q * q * p2
        denom = 4.0 * p * p * q * q
        info[j] = np.var(mu) / denom if denom > 0 else 0.0
        safe = np.clip(probs[j, keep], 1e-300, 1.0)
        entropy[j] = np.mean(-np.sum(safe * np.log(safe), axis=1) / math.log(3.0))
        max_prob[j] = np.mean(np.max(probs[j, keep], axis=1))
        xd_var[j] = np.var(p1)
        non_one_hot += int(np.sum(np.max(probs[j, keep], axis=1) < 1.0 - 1e-6))
        observed_total += int(keep.sum())
        z = np.zeros(n)
        z[keep] = mu
        z_columns.append(z)
        denominator += denom
    z = np.asarray(z_columns).T
    kd = z @ z.T / denominator
    dosage = probs[:, :, 1] + 2.0 * probs[:, :, 2]
    observed = ~missing
    dosage_corr = np.corrcoef(dosage[observed], calls[observed])[0, 1]
    dosage_rmse = math.sqrt(np.mean((dosage[observed] - calls[observed]) ** 2))
    one_hot_truth = np.eye(3)[calls]
    brier = np.mean(np.sum((probs[observed] - one_hot_truth[observed]) ** 2, axis=1))
    return {
        "non_one_hot_fraction": non_one_hot / observed_total,
        "mean_info_d": float(info.mean()),
        "p05_info_d": float(np.quantile(info, 0.05)),
        "mean_entropy": float(entropy.mean()),
        "mean_max_probability": float(max_prob.mean()),
        "dosage_correlation": float(dosage_corr),
        "dosage_rmse": dosage_rmse,
        "brier_score": float(brier),
        "mean_xd_variance": float(xd_var.mean()),
        "kd_mean_diagonal": float(np.diag(kd).mean()),
        "kd_max_abs": float(np.max(np.abs(kd))),
        "kd_symmetry_error": float(np.max(np.abs(kd - kd.T))),
        "finite": bool(np.isfinite(kd).all() and np.isfinite(info).all()),
    }, info, entropy


def read_meta(path):
    result = {}
    with Path(path).open() as handle:
        next(handle)
        for line in handle:
            key, value = line.rstrip("\n").split("\t", 1)
            result[key] = value
    return result


def read_assoc(path, causal_id):
    pvalues = []
    causal = None
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == EXPECTED_HEADER, "unified output header changed")
        for row in reader:
            pvalue = float(row["P_D"])
            pvalues.append(pvalue)
            if row["SNP"] == causal_id:
                causal = (float(row["BETA_D"]), float(row["SE_D"]), pvalue)
    require(causal is not None, f"causal marker {causal_id} was filtered")
    rank = 1 + sum(value < causal[2] for value in pvalues)
    return np.asarray(pvalues), causal, rank


def write_tsv(path, fieldnames, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--plink2", default=os.environ.get(
        "PLINK2_BIN", "plink2"))
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--m", type=int, default=1200)
    parser.add_argument("--reps", type=int, default=100)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260929)
    parser.add_argument("--effect", type=float, default=0.55)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if not args.smoke:
        require(args.n >= 1000, "production gate requires N >= 1000")
        require(args.m >= 1000, "production gate requires M >= 1000")
        require(args.reps >= 100, "production gate requires at least 100 replicates")
    require(Path(args.scanner).exists(), "scanner binary not found")
    args.plink2 = shutil.which(args.plink2) or args.plink2
    require(Path(args.plink2).exists(), "plink2 not found")

    work = Path(args.workdir)
    if work.exists():
        shutil.rmtree(work)
    (work / "fixtures").mkdir(parents=True)
    (work / "phenotypes").mkdir()
    (work / "runs").mkdir()
    samples = [f"E{i + 1:05d}" for i in range(args.n)]
    calls, frequencies, causal = generate_hardcalls(args.n, args.m, args.seed)
    causal_id = f"rs{causal + 1:06d}"
    bed_prefix = work / "fixtures" / "panel"
    write_bed(bed_prefix, calls, samples)
    subprocess.run([args.plink2, "--bfile", bed_prefix, "--make-pgen",
                    "--out", work / "fixtures" / "panel_pgen"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rng = np.random.default_rng(args.seed + 1)
    noise = rng.normal(size=calls.shape)
    category = rng.choice(3, size=calls.shape, p=(0.50, 0.35, 0.15))
    quality_factor = np.choose(category, (0.55, 1.0, 1.8))
    missing = rng.random(calls.shape) < 0.01
    missing[0] = rng.random(args.n) < 0.08
    missing[1] = rng.random(args.n) < 0.03

    bgen_paths = {}
    uncertainty_rows = []
    fixture_probabilities = {}
    fixture_info = {}
    fixture_entropy = {}
    for condition, tau in TAUS.items():
        raw_probs = posterior_probabilities(calls, frequencies, tau, noise, quality_factor)
        path = work / "fixtures" / f"{condition}.bgen"
        probs = write_bgen(path, samples, raw_probs, missing)
        metrics, info, entropy = probability_metrics(probs, missing, calls, frequencies)
        print(f"{condition} uncertainty metrics: {json.dumps(metrics, sort_keys=True)}")
        require(metrics["non_one_hot_fraction"] > 0.80,
                f"{condition} is not genuinely probabilistic")
        require(metrics["finite"] and metrics["kd_symmetry_error"] < 1e-12,
                f"{condition} Xd/KD stability failed")
        bgen_paths[condition] = path
        fixture_probabilities[condition] = probs
        fixture_info[condition] = info
        fixture_entropy[condition] = entropy
        uncertainty_rows.append({"condition": condition, "tau": tau, **metrics,
                                 "missing_fraction": float(missing.mean())})
    write_tsv(work / "uncertainty_qc.tsv", list(uncertainty_rows[0]), uncertainty_rows)

    examples = []
    representative = (("high_confidence", 0), ("medium_uncertainty", 1),
                      ("low_confidence", 2))
    medium_probs = fixture_probabilities["bgen_medium"]
    for label, cat in representative:
        locations = np.argwhere((category == cat) & (~missing))
        j, i = locations[len(locations) // 2]
        examples.append({"class": label, "snp": f"rs{j + 1:06d}", "sample": samples[i],
                         "latent_G": int(calls[j, i]), "P_AA": medium_probs[j, i, 2],
                         "P_AB": medium_probs[j, i, 1], "P_BB": medium_probs[j, i, 0],
                         "Xd": medium_probs[j, i, 1] - medium_probs[j, ~missing[j], 1].mean(),
                         "INFO_D": fixture_info["bgen_medium"][j],
                         "entropy": fixture_entropy["bgen_medium"][j]})
    j, i = np.argwhere(missing)[0]
    examples.append({"class": "hard_missing", "snp": f"rs{j + 1:06d}",
                     "sample": samples[i], "latent_G": int(calls[j, i]),
                     "P_AA": 0.0, "P_AB": 0.0, "P_BB": 0.0, "Xd": 0.0,
                     "INFO_D": fixture_info["bgen_medium"][j],
                     "entropy": fixture_entropy["bgen_medium"][j]})
    write_tsv(work / "probability_examples.tsv", list(examples[0]), examples)

    p = calls.mean(axis=1) / 2.0
    q = 1.0 - p
    z = np.where(calls == 0, -2.0 * p[:, None] ** 2,
                 np.where(calls == 1, 2.0 * p[:, None] * q[:, None],
                          -2.0 * q[:, None] ** 2)).T
    denom = np.sum(4.0 * p * p * q * q) - 4.0 * p[causal] ** 2 * q[causal] ** 2
    z[:, causal] = 0.0
    causal_x = (calls[causal] == 1).astype(float)
    causal_x -= causal_x.mean()
    phenotypes = []
    rng_pheno = np.random.default_rng(args.seed + 2)
    for rep in range(args.reps):
        effects = rng_pheno.normal(0.0, math.sqrt(0.25 / denom), args.m)
        effects[causal] = 0.0
        y_null = z @ effects + rng_pheno.normal(0.0, math.sqrt(0.75), args.n)
        for arm, values in (("null", y_null), ("alt", y_null + args.effect * causal_x)):
            path = work / "phenotypes" / f"rep{rep:03d}_{arm}.pheno"
            with path.open("w") as out:
                for sample, value in zip(samples, values):
                    out.write(f"{sample} {sample} {value:.12g}\n")
            phenotypes.append((rep, arm, path))

    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    (work / "config.json").write_text(json.dumps({
        **vars(args), "taus": TAUS, "qc_policy": POLICY_VERSION,
        "scanner_sha256": hashlib.sha256(Path(args.scanner).read_bytes()).hexdigest(),
        "probability_order": "G counts stored A; P_AA=P2, P_AB=P1, P_BB=P0",
    }, indent=2) + "\n")
    task_rows = []

    def run_task(condition, rep, arm, phenotype):
        prefix = work / "runs" / f"{condition}_{rep:03d}_{arm}"
        if condition == "bed":
            input_args = ["--bfile", bed_prefix]
        elif condition == "pgen":
            input_args = ["--pfile", work / "fixtures" / "panel_pgen"]
        else:
            input_args = ["--bgen", bgen_paths[condition]]
        command = ([args.scanner] + input_args + ["--pheno", phenotype, "--gwas",
                   "--gwas-model", "dominance", "--maf", "0.05", "--hwe", "0",
                   "--no-loco", "--out", prefix])
        started = time.monotonic()
        proc = subprocess.run([str(v) for v in command], text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        Path(str(prefix) + ".log").write_text(proc.stdout)
        require(proc.returncode == 0, f"task failed: {' '.join(map(str, command))}\n{proc.stdout[-4000:]}")
        assoc_path = Path(str(prefix) + ".assoc.dom.txt")
        meta_path = Path(str(prefix) + ".scan.meta.txt")
        pvalues, causal_stats, rank = read_assoc(assoc_path, causal_id)
        meta = read_meta(meta_path)
        require(meta.get("MODEL") == "dominance", "MODEL metadata mismatch")
        require(meta.get("D_NULL_CONVERGED_ALL") == "true", "null REML did not converge")
        require(np.isfinite(pvalues).all() and np.all((pvalues >= 0) & (pvalues <= 1)),
                "non-finite or invalid production P-values")
        require(all(math.isfinite(value) for value in causal_stats) and causal_stats[1] > 0,
                "invalid causal beta/SE/P")
        if condition.startswith("bgen"):
            require(meta.get("GENOTYPE_BACKEND") == "BGEN", "BGEN backend metadata missing")
            require(meta.get("PROBABILITY_MODE") == "true", "probability metadata missing")
            require(meta.get("DOMINANCE_ENCODING") == "posterior_mean_vitezica",
                    "dominance encoding metadata mismatch")
            require(meta.get("QC_POLICY_VERSION") == POLICY_VERSION, "QC policy metadata mismatch")
        if condition == "pgen":
            require(meta.get("PGEN_HARDCALL_ONLY") == "true", "PGEN hardcall contract missing")
        beta, se, pvalue = causal_stats
        return {"condition": condition, "replicate": rep, "arm": arm,
                "markers": len(pvalues), "elapsed_seconds": time.monotonic() - started,
                "causal_rank": rank, "causal_beta": beta, "causal_se": se,
                "causal_p": pvalue, "pvalues": pvalues,
                "null_converged": True,
                "sigma_d2": float(meta["SIGMA_D2"]), "sigma_e2": float(meta["SIGMA_E2"])}

    futures = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for condition in CONDITIONS:
            for rep, arm, phenotype in phenotypes:
                futures.append(executor.submit(run_task, condition, rep, arm, phenotype))
        for future in as_completed(futures):
            task_rows.append(future.result())
    task_rows.sort(key=lambda row: (row["condition"], row["replicate"], row["arm"]))

    min_markers = 1000 if not args.smoke else max(1, int(args.m * 0.75))
    require(all(row["markers"] >= min_markers for row in task_rows),
            "too few tested markers in at least one production scan")
    bed_rows = {(row["replicate"], row["arm"]): row for row in task_rows if row["condition"] == "bed"}
    pgen_rows = {(row["replicate"], row["arm"]): row for row in task_rows if row["condition"] == "pgen"}
    max_hardcall_diff = 0.0
    for key in bed_rows:
        require(bed_rows[key]["pvalues"].shape == pgen_rows[key]["pvalues"].shape,
                "BED/PGEN output marker counts differ")
        max_hardcall_diff = max(max_hardcall_diff, float(np.max(np.abs(
            bed_rows[key]["pvalues"] - pgen_rows[key]["pvalues"]))))
        for field in ("causal_beta", "causal_se", "causal_p"):
            max_hardcall_diff = max(max_hardcall_diff,
                                    abs(bed_rows[key][field] - pgen_rows[key][field]))
    require(max_hardcall_diff <= 2e-5, "BED/PGEN dominance regression mismatch")

    null_summary = []
    power_summary = []
    replicate_output = []
    for condition in CONDITIONS:
        null_rows = [row for row in task_rows if row["condition"] == condition and row["arm"] == "null"]
        alt_rows = [row for row in task_rows if row["condition"] == condition and row["arm"] == "alt"]
        pooled = np.concatenate([row["pvalues"] for row in null_rows])
        chisq = stats.chi2.isf(np.clip(pooled, 1e-300, 1.0), 1)
        null_summary.append({
            "condition": condition, "replicates": len(null_rows), "pvalues": len(pooled),
            "lambda_gc": float(np.median(chisq) / stats.chi2.ppf(0.5, 1)),
            "type1_0.05": float(np.mean(pooled < 0.05)),
            "type1_0.01": float(np.mean(pooled < 0.01)),
            "type1_0.001": float(np.mean(pooled < 0.001)),
            "type1_bonferroni": float(np.mean(pooled < 0.05 / null_rows[0]["markers"])),
            "ks_d": float(stats.kstest(pooled, "uniform").statistic),
            "p01": float(np.quantile(pooled, 0.01)), "p05": float(np.quantile(pooled, 0.05)),
            "p50": float(np.quantile(pooled, 0.50)), "p95": float(np.quantile(pooled, 0.95)),
            "p99": float(np.quantile(pooled, 0.99)),
        })
        bonf = [row["causal_p"] < 0.05 / row["markers"] for row in alt_rows]
        power_summary.append({
            "condition": condition, "replicates": len(alt_rows),
            "nominal_power": float(np.mean([row["causal_p"] < 0.05 for row in alt_rows])),
            "bonferroni_power": float(np.mean(bonf)),
            "median_rank": float(np.median([row["causal_rank"] for row in alt_rows])),
            "mean_log10p": float(np.mean([-math.log10(max(row["causal_p"], 1e-300)) for row in alt_rows])),
            "mean_beta": float(np.mean([row["causal_beta"] for row in alt_rows])),
            "mean_se": float(np.mean([row["causal_se"] for row in alt_rows])),
            "coverage_95": float(np.mean([abs(row["causal_beta"] - args.effect) <=
                                           1.96 * row["causal_se"] for row in alt_rows])),
        })
    for row in task_rows:
        replicate_output.append({key: value for key, value in row.items() if key != "pvalues"})
    write_tsv(work / "replicate_metrics.tsv", list(replicate_output[0]), replicate_output)
    write_tsv(work / "null_calibration_summary.tsv", list(null_summary[0]), null_summary)
    write_tsv(work / "power_summary.tsv", list(power_summary[0]), power_summary)

    verdict_checks = {
        "all_tasks_completed": len(task_rows) == len(CONDITIONS) * args.reps * 2,
        "minimum_marker_count": all(row["markers"] >= min_markers for row in task_rows),
        "bed_pgen_parity": max_hardcall_diff <= 2e-5,
        "fixture_non_one_hot": all(row["non_one_hot_fraction"] > 0.80 for row in uncertainty_rows),
        "null_lambda": all(0.90 <= row["lambda_gc"] <= 1.10 for row in null_summary),
        "null_type1": all(0.035 <= row["type1_0.05"] <= 0.065 for row in null_summary),
        "null_type1_0.01": all(0.006 <= row["type1_0.01"] <= 0.014 for row in null_summary),
        "null_type1_0.001": all(0.0004 <= row["type1_0.001"] <= 0.0016 for row in null_summary),
        "null_ks": all(row["ks_d"] <= 0.03 for row in null_summary),
    }
    power_by_condition = {row["condition"]: row for row in power_summary}
    verdict_checks["hardcall_power"] = power_by_condition["bed"]["bonferroni_power"] >= 0.75
    verdict_checks["uncertainty_power_trend"] = (
        power_by_condition["bgen_high"]["bonferroni_power"] + 0.05 >=
        power_by_condition["bgen_medium"]["bonferroni_power"] and
        power_by_condition["bgen_medium"]["bonferroni_power"] + 0.05 >=
        power_by_condition["bgen_low"]["bonferroni_power"])
    verdict_checks["uncertainty_se_trend"] = (
        power_by_condition["bgen_high"]["mean_se"] <=
        power_by_condition["bgen_medium"]["mean_se"] <=
        power_by_condition["bgen_low"]["mean_se"])
    if args.smoke:
        statistical_keys = ("null_lambda", "null_type1", "null_type1_0.01", "null_type1_0.001",
                            "null_ks", "hardcall_power",
                            "uncertainty_power_trend", "uncertainty_se_trend")
        for key in statistical_keys:
            verdict_checks[key] = True
    verdict = "PASS" if all(verdict_checks.values()) else "FAIL"
    summary = {"verdict": verdict, "n": args.n, "m": args.m, "replicates": args.reps,
               "effect": args.effect, "causal_snp": causal_id,
               "max_bed_pgen_causal_difference": max_hardcall_diff,
               "checks": verdict_checks, "uncertainty": uncertainty_rows,
               "null": null_summary, "power": power_summary}
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    require(verdict == "PASS", "BGEN dominance production gate failed")
    print("BGEN DOMINANCE PRODUCTION GATE PASSED")


if __name__ == "__main__":
    main()
