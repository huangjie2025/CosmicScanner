#!/usr/bin/env python3
"""E3.1: independent BGEN AD kernels/null VCE, with no GWAS scan."""

import argparse
import csv
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from run_bgen_dominance_production_gate import (
    posterior_probabilities, read_meta, require, write_bgen, write_tsv,
)


def write_hardcalls(prefix, calls, samples):
    m, n = calls.shape
    Path(str(prefix) + ".fam").write_text("".join(f"{s} {s} 0 0 1 -9\n" for s in samples))
    Path(str(prefix) + ".bim").write_text("".join(
        f"{1 if j < m // 2 else 2} rs{j + 1:06d} 0 {(j + 1) * 100} A C\n" for j in range(m)))
    body = bytearray(b"\x6c\x1b\x01")
    for row in calls:
        packed = bytearray((n + 3) // 4)
        for i, g in enumerate(row):
            code = 1 if g < 0 else (3, 2, 0)[int(g)]
            packed[i // 4] |= code << (2 * (i % 4))
        body.extend(packed)
    Path(str(prefix) + ".bed").write_bytes(body)


def phenotype_matrix(path, values):
    np.savetxt(path, values, header=f"{values.shape[0]} {values.shape[1]}", comments="", fmt="%.17g")


def read_estimates(prefix):
    with Path(str(prefix) + ".vce.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def probability_oracle(probabilities, missing):
    probabilities = probabilities.astype(np.float32).astype(float)
    obs = ~missing
    counts = obs.sum(axis=1)
    x = probabilities[:, :, 1] + 2 * probabilities[:, :, 2]
    p = (x * obs).sum(axis=1) / (2 * np.maximum(counts, 1))
    q = 1 - p
    selected = (counts > 0) & (np.minimum(p, q) >= .01)
    za = np.where(obs, x - 2 * p[:, None], 0.0)[selected].T
    mu = (-2 * p[:, None] ** 2 * probabilities[:, :, 0] +
          2 * p[:, None] * q[:, None] * probabilities[:, :, 1] -
          2 * q[:, None] ** 2 * probabilities[:, :, 2])
    zd = np.where(obs, mu, 0.0)[selected].T
    return (za @ za.T / np.sum(2 * p[selected] * q[selected]),
            zd @ zd.T / np.sum(4 * p[selected] ** 2 * q[selected] ** 2))


def diagnose(ka, kd):
    n = len(ka)
    outside = np.triu_indices(n, 1)
    corr = float(np.corrcoef(ka[outside], kd[outside])[0, 1])
    contrast = np.linalg.qr(np.ones((n, 1)), mode="complete")[0][:, 1:]
    a, d = contrast.T @ ka @ contrast, contrast.T @ kd @ contrast
    identity = np.eye(n - 1)
    v = .3 * a + .3 * d + .4 * identity
    li = np.linalg.inv(np.linalg.cholesky(v))
    kernels = [li @ k @ li.T for k in (a, d, identity)]
    fisher = np.array([[.5 * np.sum(x * y) for y in kernels] for x in kernels])
    scale = np.sqrt(np.diag(fisher))
    normalized = fisher / np.outer(scale, scale)
    return {"corr_ka_kd_offdiag": corr, "fisher_rank": int(np.linalg.matrix_rank(fisher)),
            "fisher_condition": float(np.linalg.cond(fisher)),
            "normalized_fisher_condition": float(np.linalg.cond(normalized)),
            "mean_ka_diagonal": float(np.diag(ka).mean()),
            "mean_kd_diagonal": float(np.diag(kd).mean())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--reps", type=int, default=30)
    parser.add_argument("--legacy-control", action="store_true",
                        help="Use unchanged EigenVceAD for the pre-robustification control")
    parser.add_argument("--plink2", default=os.environ.get(
        "PLINK2_BIN", "plink2"))
    args = parser.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    rng = np.random.default_rng(3101)
    n, m = 160, 120
    samples = [f"N{i:04d}" for i in range(n)]
    freqs = rng.uniform(.15, .45, m)
    latent = rng.binomial(2, freqs[:, None], (m, n))
    missing = rng.random((m, n)) < .02
    latent[0] = 0                  # monomorphic rejection
    missing[-1] = True            # all-missing rejection
    hardcalls = np.where(missing, -1, latent)

    def command(values):
        proc = subprocess.run([str(v) for v in values], env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        require(proc.returncode == 0, proc.stdout[-5000:])
        return proc

    def materials(prefix, calls, names):
        write_hardcalls(prefix, calls, names)
        command([args.plink2, "--bfile", prefix, "--make-pgen", "--out", str(prefix) + "_pgen"])
        command([args.plink2, "--bfile", prefix, "--export", "bgen-1.2", "--out", str(prefix) + "_bgen"])

    def run(format_name, source, target, keep=None, chromosome="-", info=0.0, phenotype="-", maxit=100):
        index_arg = "-" if keep is None else ",".join(str(i) for i in keep)
        driver_format = "bgen-legacy" if args.legacy_control and format_name == "bgen" else format_name
        proc = command([args.driver, driver_format, source, target, index_arg,
                        chromosome, info, phenotype, maxit])
        Path(str(target) + ".log").write_text(proc.stdout)
        meta = read_meta(str(target) + ".meta.tsv")
        count = int(meta["N"])
        ka = np.fromfile(str(target) + ".ka.bin", dtype="<f8").reshape(count, count, order="F")
        kd = np.fromfile(str(target) + ".kd.bin", dtype="<f8").reshape(count, count, order="F")
        require(np.isfinite(ka).all() and np.isfinite(kd).all(), "nonfinite dual kernels")
        require(float(meta["LOCO_MAX_ERROR"]) < 1e-10, "LOCO tolerance")
        return ka, kd, meta

    prefix = work / "hardcall"
    materials(prefix, hardcalls, samples)
    keep = np.array([i for i in range(n) if i % 5 != 0])[::-1]
    subset_prefix = work / "subset"
    materials(subset_prefix, hardcalls[:, keep], [samples[i] for i in keep])
    y = rng.normal(size=(n, 1))
    one_hot = np.eye(3)[latent]
    oracle_a, oracle_d = probability_oracle(one_hot, missing)
    mixed_rng = np.random.default_rng(3102)
    mixed_values = np.linalg.cholesky(.45 * oracle_a + .35 * oracle_d + .20 * np.eye(n)) @ mixed_rng.normal(size=(n, 3))
    y = np.column_stack([y, mixed_values])
    pheno = work / "parity.pheno"
    phenotype_matrix(pheno, y)
    sub_pheno = work / "parity_subset.pheno"
    phenotype_matrix(sub_pheno, y[keep])
    parity_rows = []
    global_pair = None
    for case, selection, chrom in (("all", None, "-"), ("sample_subset", keep, "-"),
                                    ("chromosome_subset", None, "1")):
        reference = subset_prefix if selection is not None else prefix
        outputs = {}
        for format_name, source in (("bed", reference),
                                   ("pgen", str(reference) + "_pgen.pgen"),
                                   ("bgen", str(prefix) + "_bgen.bgen")):
            target = work / f"{case}_{format_name}"
            outputs[format_name] = run(format_name, source, target,
                                      keep=selection if format_name == "bgen" else None,
                                      chromosome=chrom,
                                      phenotype=sub_pheno if selection is not None else pheno)
        for format_name in ("pgen", "bgen"):
            max_a = float(np.max(np.abs(outputs["bed"][0] - outputs[format_name][0])))
            max_d = float(np.max(np.abs(outputs["bed"][1] - outputs[format_name][1])))
            require(max(max_a, max_d) < 2e-6, f"one-hot kernel parity failed {case} {format_name}")
            require(outputs["bed"][2]["MARKERS"] == outputs[format_name][2]["MARKERS"], "marker membership count")
            bed_fits = read_estimates(work / f"{case}_bed")
            other_fits = read_estimates(work / f"{case}_{format_name}")
            vce_diff = max(abs(float(bed_fit[key]) - float(other_fit[key]))
                           for bed_fit, other_fit in zip(bed_fits, other_fits)
                           for key in ("sigma_a2", "sigma_d2", "sigma_e2"))
            require(vce_diff < 5e-4, f"one-hot VC parity failed {case} {format_name}: {vce_diff}")
            require(all(a["converged"] == b["converged"] for a, b in zip(bed_fits, other_fits)),
                    "VCE convergence parity")
            parity_rows.append({"case": case, "format": format_name, "ka_max_error": max_a,
                                "kd_max_error": max_d, "vce_max_error": vce_diff,
                                "loco_max_error": float(outputs[format_name][2]["LOCO_MAX_ERROR"])})
        if case == "all":
            global_pair = outputs["bgen"][:2]

    # Reordering/non-contiguous keep indices must equal physical subset BGEN.
    sub = run("bgen", str(subset_prefix) + "_bgen.bgen", work / "subset_materialized")
    indexed_a = np.fromfile(work / "sample_subset_bgen.ka.bin", dtype="<f8").reshape(len(keep), len(keep))
    indexed_d = np.fromfile(work / "sample_subset_bgen.kd.bin", dtype="<f8").reshape(len(keep), len(keep))
    require(np.max(np.abs(indexed_a - sub[0])) < 1e-12 and
            np.max(np.abs(indexed_d - sub[1])) < 1e-12, "materialized sample subset identity")

    noise = rng.normal(size=latent.shape)
    quality = np.ones(latent.shape)
    uncertainty_rows, uncertain_pairs = [], {}
    for level, tau in (("low_uncertainty", .20), ("medium_uncertainty", .65), ("high_uncertainty", 1.20)):
        probs = posterior_probabilities(latent, freqs, tau, noise, quality)
        probs[0, :, :] = np.array([1.0, 0.0, 0.0])
        for j in (7, 9):
            probs[j] = np.array([(1 - freqs[j]) ** 2, 2 * freqs[j] * (1 - freqs[j]), freqs[j] ** 2])
        source = work / f"{level}.bgen"
        decoded = write_bgen(source, samples, probs, missing)
        ka, kd, meta = run("bgen", source, work / level)
        expected_a, expected_d = probability_oracle(decoded, missing)
        oracle_error = max(float(np.max(np.abs(ka - expected_a))), float(np.max(np.abs(kd - expected_d))))
        require(oracle_error < 1e-10, f"uncertain probability KA/KD oracle failed: {oracle_error}")
        qc_a, qc_d, qc_meta = run("bgen", source, work / f"{level}_qc", info=.01)
        uncertainty_rows.append({"quality": level, "tau": tau, "markers": int(meta["MARKERS"]),
                                 "markers_info_d_0.01": int(qc_meta["MARKERS"]),
                                 "probability_oracle_error": oracle_error,
                                 "loco_max_error": float(meta["LOCO_MAX_ERROR"]), **diagnose(ka, kd)})
        uncertain_pairs[level] = (ka, kd, source)

    scenarios = (("sigma_d_zero", (.50, 0.0, .50), global_pair, str(prefix) + "_bgen.bgen"),
                 ("sigma_a_zero", (0.0, .50, .50), global_pair, str(prefix) + "_bgen.bgen"),
                 ("ad_mixed", (.30, .30, .40), global_pair, str(prefix) + "_bgen.bgen"),
                 ("high_uncertainty", (.30, .30, .40), uncertain_pairs["high_uncertainty"][:2],
                  uncertain_pairs["high_uncertainty"][2]))
    recovery_rows = []
    for label, truth, pair, source in scenarios:
        ka, kd = pair
        covariance = truth[0] * ka + truth[1] * kd + truth[2] * np.eye(n)
        values = np.linalg.cholesky(covariance) @ rng.normal(size=(n, args.reps))
        phenotype = work / f"{label}.pheno"
        phenotype_matrix(phenotype, values)
        target = work / f"recovery_{label}"
        run("bgen", source, target, phenotype=phenotype, maxit=200)
        estimates = read_estimates(target)
        fitted = np.array([[float(row[k]) for k in ("sigma_a2", "sigma_d2", "sigma_e2")]
                           for row in estimates])
        require(np.isfinite(fitted).all() and np.all(fitted >= 0), f"invalid recovery {label}")
        mean = fitted.mean(axis=0)
        recovery_rows.append({"scenario": label, "replicates": args.reps,
                              "true_a": truth[0], "true_d": truth[1], "true_e": truth[2],
                              "mean_a": mean[0], "mean_d": mean[1], "mean_e": mean[2],
                              "sd_a": fitted[:, 0].std(ddof=1), "sd_d": fitted[:, 1].std(ddof=1),
                              "sd_e": fitted[:, 2].std(ddof=1),
                              "converged": sum(row["converged"] == "true" for row in estimates),
                              "max_abs_mean_error": float(np.max(np.abs(mean - truth))),
                              "variance_recovery_within_0.20": bool(np.max(np.abs(mean - truth)) <= .20)})

    # A zero-iteration null must retain the actual non-convergence status.
    run("bgen", str(prefix) + "_bgen.bgen", work / "no_convergence", phenotype=pheno, maxit=0)
    require(read_estimates(work / "no_convergence")[0]["converged"] == "false", "false convergence")
    write_tsv(work / "one_hot_parity.tsv", list(parity_rows[0]), parity_rows)
    write_tsv(work / "uncertainty.tsv", list(uncertainty_rows[0]), uncertainty_rows)
    write_tsv(work / "variance_recovery.tsv", list(recovery_rows[0]), recovery_rows)
    recovery_gate = all(row["variance_recovery_within_0.20"] and
                        row["converged"] == row["replicates"] for row in recovery_rows)
    summary = {"solver_mode": "LEGACY_CONTROL" if args.legacy_control else "BGEN_ROBUST",
               "backend_verdict": "PASS",
               "variance_recovery_and_convergence_verdict": "PASS" if recovery_gate else "BLOCKED",
               "n": n, "m": m, "seed": 3101, "parity": parity_rows,
               "uncertainty": uncertainty_rows, "variance_recovery": recovery_rows,
               "notes": "Recovery under the Gaussian posterior-mean working model; no GWAS scan. "
                        "High-uncertainty weak identifiability and boundary optimizer results are reported, not hidden."}
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("BGEN AD NULL BACKEND/ONE-HOT/DUAL-LOCO VALIDATION PASSED; see recovery diagnostics")
    print("VARIANCE RECOVERY AND CONVERGENCE GATE: " + ("PASS" if recovery_gate else "BLOCKED"))


if __name__ == "__main__":
    main()
