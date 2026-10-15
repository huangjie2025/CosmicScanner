#!/usr/bin/env python3
"""BGEN-only boundary-aware REML validation; no association scan."""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from run_bgen_ad_null_validation import phenotype_matrix, read_estimates
from run_bgen_dominance_production_gate import require, write_tsv


def reference(kernels, y, start):
    def objective(theta):
        v = sum(t * k for t, k in zip(theta, kernels))
        try:
            factor = cho_factor(v, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return 1e100, np.zeros(3)
        vy = cho_solve(factor, y, check_finite=False)
        vi = cho_solve(factor, kernels[2], check_finite=False)
        value = np.log(np.diag(factor[0])).sum() + .5 * y @ vy
        gradient = np.array([.5 * (np.sum(vi * k.T) - vy @ k @ vy) for k in kernels])
        return float(value), gradient
    fits = [minimize(objective, seed, jac=True, method="L-BFGS-B",
                     bounds=((0.0, None), (0.0, None), (1e-8, None)),
                     options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
            for seed in (start, [.3, .3, .4], [0, 0, max(.1, np.var(y))])]
    return min(fits, key=lambda fit: fit.fun), objective


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    source = Path(args.source)
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if not (source / "summary.json").exists() or json.loads((source / "summary.json").read_text()).get("solver_mode") == "BGEN_ROBUST":
        source = work / "legacy_control"
        subprocess.run([sys.executable, str(Path(__file__).with_name("run_bgen_ad_null_validation.py")),
                        "--driver", args.driver, "--workdir", str(source), "--legacy-control"],
                       env=env, check=True, stdout=subprocess.DEVNULL)
    original = json.loads((source / "summary.json").read_text())
    n = original["n"]
    m = original["m"]
    # Fail before any driver call if the shared null fixture is incomplete.
    # The producer owns this directory; CTest also declares a dependency/resource lock.
    expected_bed_bytes = 3 + ((n + 3) // 4) * m
    require((source / "hardcall.bed").stat().st_size == expected_bed_bytes,
            "incomplete BGEN AD null fixture: hardcall.bed size mismatch")
    for required in ("hardcall_pgen.pgen", "hardcall_bgen.bgen", "high_uncertainty.bgen",
                     "parity.pheno"):
        require((source / required).is_file() and (source / required).stat().st_size > 0,
                f"incomplete BGEN AD null fixture: missing/empty {required}")
    contrast = np.linalg.qr(np.ones((n, 1)), mode="complete")[0][:, 1:]
    eye = np.eye(n - 1)
    scenarios = original["variance_recovery"]
    records, summaries, start_summaries = [], [], []

    def invoke(mode, genotype, phenotype, prefix, maxit=200):
        result = subprocess.run([args.driver, mode, str(genotype), str(prefix), "-", "-", "0",
                                 str(phenotype), str(maxit)], env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        Path(str(prefix) + ".log").write_text(result.stdout)
        require(result.returncode == 0, result.stdout[-5000:])

    # Real native-path controls: the test driver routes BED/PGEN directly to
    # the unmodified EigenVceAD, never to the BGEN robust adapter.
    for fmt, genotype in (("bed", source / "hardcall"), ("pgen", source / "hardcall_pgen.pgen")):
        target = work / f"unchanged_{fmt}"
        invoke(fmt, genotype, source / "parity.pheno", target, maxit=100)
        before = read_estimates(source / f"all_{fmt}")
        after = read_estimates(target)
        require(all(old[key] == new[key] for old, new in zip(before, after)
                    for key in ("sigma_a2", "sigma_d2", "sigma_e2", "logL", "converged", "iterations")),
                f"legacy {fmt} null results changed")

    def run_scenario(item):
        label = item["scenario"]
        genotype = source / ("high_uncertainty.bgen" if label == "high_uncertainty" else "hardcall_bgen.bgen")
        prefix = work / label
        invoke("bgen", genotype, source / f"{label}.pheno", prefix)
        # Preserve exact before snapshots independently of future test reruns.
        shutil.copyfile(source / f"recovery_{label}.vce.tsv", work / f"{label}.before.tsv")
        return item

    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(run_scenario, scenarios))

    for item in scenarios:
        label = item["scenario"]
        prefix = work / label
        before = read_estimates(source / f"recovery_{label}")
        after = read_estimates(prefix)
        ka, kd = [np.fromfile(str(prefix) + suffix, dtype="<f8").reshape(n, n)
                  for suffix in (".ka.bin", ".kd.bin")]
        kernels = (contrast.T @ ka @ contrast, contrast.T @ kd @ contrast, eye)
        values = contrast.T @ np.loadtxt(source / f"{label}.pheno", skiprows=1)
        with Path(str(prefix) + ".candidates.tsv").open() as handle:
            candidates = list(csv.DictReader(handle, delimiter="\t"))
        for rep, (old, new) in enumerate(zip(before, after)):
            theta = np.array([float(new[k]) for k in ("sigma_a2", "sigma_d2", "sigma_e2")])
            require(new["converged"] == "true", f"uncertified after null {label}/{rep}")
            require(float(new["gradient_norm"]) <= float(new["gradient_tolerance"]), "gradient gate")
            require(float(new["kkt_residual"]) <= float(new["gradient_tolerance"]), "KKT gate")
            require(float(new["likelihood_delta"]) <= float(new["likelihood_tolerance"]) * max(1, abs(float(new["logL"]))),
                    "likelihood delta gate")
            subset = [row for row in candidates if int(row["replicate"]) == rep]
            require({row["model"] for row in subset} == {"FULL", "A_ZERO", "D_ZERO", "AD_ZERO"}, "all faces evaluated")
            full = [row for row in subset if row["model"] == "FULL"]
            require({"balanced", "additive_dominated", "dominance_dominated", "near_boundary"}.issubset(
                {row["start"] for row in full}), "four starts missing")
            certified = [row for row in subset if row["face_converged"] == "true" and row["kkt_valid"] == "true"]
            best = max(float(row["logL"]) for row in certified)
            require(abs(best - float(new["logL"])) < 1e-10, "maximum certified candidate not selected")
            start_summaries.append({"scenario": label, "replicate": rep,
                                    "certified_candidate_count": len(certified),
                                    "certified_logL_spread": max(float(row["logL"]) for row in certified) -
                                        min(float(row["logL"]) for row in certified),
                                    "selected_start": new["selected_start"]})
            fit, objective = reference(kernels, values[:, rep], theta)
            value, raw_gradient = objective(theta)
            scale = float(values[:, rep] @ values[:, rep] / (n - 1))
            conversion = scale / np.maximum([k.trace() / (n - 1) for k in kernels], 1e-12)
            independent_projected = raw_gradient.copy()
            independent_projected[:2] = np.where(theta[:2] == 0, np.minimum(raw_gradient[:2], 0), raw_gradient[:2])
            independent_kkt = float(np.max(np.abs(independent_projected * conversion)) / (n - 1))
            require(independent_kkt <= 1.02 * float(new["gradient_tolerance"]), "independent KKT check")
            require(abs(-value - float(new["logL"])) < 1e-8, "REML value mismatch")
            gap = value - float(fit.fun)
            require(gap < 1e-6, f"reference likelihood gap {label}/{rep}: {gap}")
            gain = float(new["logL"]) - float(old["logL"])
            require(gain >= -1e-7, "robustification reduced likelihood")
            records.append({"scenario": label, "replicate": rep,
                            "before_converged": old["converged"], "after_converged": new["converged"],
                            "likelihood_gain": gain, "reference_gap": gap,
                            "independent_kkt_residual": independent_kkt,
                            "boundary_status": new["boundary_status"],
                            "before_a": float(old["sigma_a2"]), "before_d": float(old["sigma_d2"]),
                            "before_e": float(old["sigma_e2"]),
                            "after_a": theta[0], "after_d": theta[1], "after_e": theta[2]})
        selected = [row for row in records if row["scenario"] == label]
        old_mean = np.mean([[row["before_a"], row["before_d"], row["before_e"]] for row in selected], axis=0)
        new_mean = np.mean([[row["after_a"], row["after_d"], row["after_e"]] for row in selected], axis=0)
        truth = [item["true_a"], item["true_d"], item["true_e"]]
        summaries.append({"scenario": label, "replicates": len(selected), "truth": truth,
                          "before_converged": sum(row["before_converged"] == "true" for row in selected),
                          "after_converged": len(selected), "before_mean": old_mean.tolist(), "after_mean": new_mean.tolist(),
                          "before_bias": (old_mean - truth).tolist(), "after_bias": (new_mean - truth).tolist(),
                          "maximum_likelihood_gain": max(row["likelihood_gain"] for row in selected),
                          "maximum_reference_gap": max(row["reference_gap"] for row in selected),
                          "maximum_kkt_residual": max(row["independent_kkt_residual"] for row in selected),
                          "exact_a_zero": int(sum(row["after_a"] == 0 for row in selected)),
                          "exact_d_zero": int(sum(row["after_d"] == 0 for row in selected))})

    # Target a residual-only optimum to verify the exact A=D=0 face explicitly.
    ka, kd = [np.fromfile(source / ("all_bgen" + suffix), dtype="<f8").reshape(n, n)
              for suffix in (".ka.bin", ".kd.bin")]
    _, vectors = np.linalg.eigh(contrast.T @ (ka + kd) @ contrast)
    residual_y = np.sqrt(n - 1) * contrast @ vectors[:, 0]
    residual_pheno = work / "ad_zero.pheno"
    phenotype_matrix(residual_pheno, residual_y[:, None])
    invoke("bgen", source / "hardcall_bgen.bgen", residual_pheno, work / "ad_zero")
    corner = read_estimates(work / "ad_zero")[0]
    require(float(corner["sigma_a2"]) == 0 and float(corner["sigma_d2"]) == 0 and corner["converged"] == "true",
            "exact residual-only face failed")
    require(abs(float(corner["sigma_e2"]) - 1) < 1e-8, "analytic residual variance")
    invoke("bgen", source / "hardcall_bgen.bgen", source / "parity.pheno", work / "maxit_zero", maxit=0)
    require(all(row["converged"] == "false" for row in read_estimates(work / "maxit_zero")), "maxit=0 false success")

    write_tsv(work / "before_after.tsv", list(records[0]), records)
    write_tsv(work / "multistart.tsv", list(start_summaries[0]), start_summaries)
    summary = {"numerical_reml_verdict": "PASS", "n": n, "m": original["m"],
               "scenarios": summaries, "legacy_bed_pgen_unchanged": True,
               "exact_double_zero_verified": True, "maxit_zero_not_converged": True,
               "maximum_certified_start_spread": max(row["certified_logL_spread"] for row in start_summaries),
               "scope": "Numerical Gaussian working-null validation; high-uncertainty VC precision remains limited."}
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
