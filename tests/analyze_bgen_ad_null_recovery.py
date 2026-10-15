#!/usr/bin/env python3
"""Independent contrast-REML reference for saved E3.1 null estimates.

This is an offline diagnostic, not a replacement production optimizer or LRT.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from run_bgen_ad_null_validation import read_estimates
from run_bgen_dominance_production_gate import write_tsv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    work = Path(args.workdir)
    summary = json.loads((work / "summary.json").read_text())
    n = summary["n"]
    contrast = np.linalg.qr(np.ones((n, 1)), mode="complete")[0][:, 1:]
    rows = []
    summaries = []
    for scenario in summary["variance_recovery"]:
        label = scenario["scenario"]
        prefix = work / f"recovery_{label}"
        ka = np.fromfile(str(prefix) + ".ka.bin", dtype="<f8").reshape(n, n)
        kd = np.fromfile(str(prefix) + ".kd.bin", dtype="<f8").reshape(n, n)
        kernels = (contrast.T @ ka @ contrast, contrast.T @ kd @ contrast, np.eye(n - 1))
        values = contrast.T @ np.loadtxt(work / f"{label}.pheno", skiprows=1)
        estimates = read_estimates(prefix)
        for rep, estimate in enumerate(estimates):
            y = values[:, rep]

            def objective(theta):
                covariance = sum(coefficient * kernel for coefficient, kernel in zip(theta, kernels))
                try:
                    factor = cho_factor(covariance, lower=True, check_finite=False)
                except np.linalg.LinAlgError:
                    return 1e100, np.zeros(3)
                inverse_y = cho_solve(factor, y, check_finite=False)
                inverse = cho_solve(factor, kernels[2], check_finite=False)
                nll = np.log(np.diag(factor[0])).sum() + .5 * y @ inverse_y
                gradient = np.array([.5 * (np.sum(inverse * k.T) - inverse_y @ k @ inverse_y)
                                     for k in kernels])
                return float(nll), gradient

            original = np.array([float(estimate[k]) for k in ("sigma_a2", "sigma_d2", "sigma_e2")])
            starts = (original, np.array([.30, .30, .40]), np.array([0.0, 0.0, max(.1, np.var(y))]))
            solutions = [minimize(objective, start, jac=True, method="L-BFGS-B",
                                  bounds=((0.0, None), (0.0, None), (1e-8, None)),
                                  options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-7})
                         for start in starts]
            best = min(solutions, key=lambda result: result.fun)
            original_nll, _ = objective(original)
            rows.append({"scenario": label, "replicate": rep,
                         "solver_converged": estimate["converged"],
                         "reference_success": bool(best.success),
                         "likelihood_gap": float(original_nll - best.fun),
                         "sigma_a2_reference": float(best.x[0]),
                         "sigma_d2_reference": float(best.x[1]),
                         "sigma_e2_reference": float(best.x[2])})
        selected = [row for row in rows if row["scenario"] == label]
        summaries.append({"scenario": label,
                          "reference_successes": sum(row["reference_success"] for row in selected),
                          "maximum_likelihood_gap": max(row["likelihood_gap"] for row in selected),
                          "mean_a_reference": float(np.mean([row["sigma_a2_reference"] for row in selected])),
                          "mean_d_reference": float(np.mean([row["sigma_d2_reference"] for row in selected])),
                          "mean_e_reference": float(np.mean([row["sigma_e2_reference"] for row in selected]))})
    write_tsv(work / "reml_reference.tsv", list(rows[0]), rows)
    (work / "reml_reference_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
