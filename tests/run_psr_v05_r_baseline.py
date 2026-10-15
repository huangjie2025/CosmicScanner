#!/usr/bin/env python3
"""PSR-GWAS v0.5 R ordinal baseline runner.

Generates ordinal data, runs R fixed-effect ordinal GWAS methods (MASS::polr,
ordinal::clm) as external baselines, and outputs availability + summary TSVs.

R methods are optional: missing R or missing packages are recorded, not fatal.

Environment variables:
  PSR_V05_BENCH_REPS     Number of replicates (default: 1)
  PSR_V05_BENCH_N        Sample size (default: 80)
  PSR_V05_BENCH_M        Number of SNPs (default: 100)
  PSR_V05_BENCH_LEVELS   Ordinal levels (default: 3)
  PSR_V05_BENCH_OUTDIR   Output directory override
  PSR_V05_BENCH_SEED     Base seed (default: 20260601)
  PSR_V05_R_SCRIPT       Rscript binary path (default: auto-detect)
"""

import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_psr_simulation_summary as smoke


ROOT = Path(__file__).resolve().parents[1]
BENCH_REPS = int(os.environ.get("PSR_V05_BENCH_REPS", "1"))
BENCH_N = int(os.environ.get("PSR_V05_BENCH_N", "80"))
BENCH_M = int(os.environ.get("PSR_V05_BENCH_M", "100"))
BENCH_LEVELS = int(os.environ.get("PSR_V05_BENCH_LEVELS", "3"))
BENCH_SEED = int(os.environ.get("PSR_V05_BENCH_SEED", "20260601"))
OUT_DIR = Path(
    os.environ.get(
        "PSR_V05_BENCH_OUTDIR",
        ROOT / "tests" / "results" / "psr_v05_ordinal_benchmark",
    )
)
WORK_DIR = OUT_DIR / "work"

R_METHODS = ["MASS::polr", "ordinal::clm"]


def find_rscript():
    candidates = [
        os.environ.get("PSR_V05_R_SCRIPT", os.environ.get("RSCRIPT_BIN", str(Path(os.environ.get("CONDA_PREFIX", "/nonexistent")) / "bin/Rscript"))),
        "Rscript",
    ]
    for c in candidates:
        try:
            subprocess.run([c, "--version"], capture_output=True)
            return c
        except FileNotFoundError:
            pass
    return None


RSCRIPT = find_rscript()


def simulate_ordinal_data(rng, n, m, levels, effect=0.0):
    G = rng.binomial(2, 0.3, size=(n, m)).astype(np.float64)
    mu = 0.0
    if effect > 0.0 and m >= 3:
        causal_idx = max(1, m // 3) - 1
        G[:, causal_idx] = rng.binomial(2, 0.3, size=n).astype(np.float64)
        mu = G[:, causal_idx] * effect
    latent = mu + rng.normal(0, 1.0, size=n)
    spacing = 1.2
    thresholds = [(-(levels - 2) / 2.0 + k) * spacing for k in range(levels - 1)]
    y = np.zeros(n, dtype=int)
    for k in range(1, levels):
        y[latent > thresholds[k - 1]] = k
    return G, y


def scenario_effect(scenario):
    if scenario == "power_high":
        return 0.7
    if scenario == "power_low":
        return 0.35
    return 0.0


def write_plink(prefix, G, y, seed):
    smoke.write_plink(prefix, G)
    smoke.write_one_col_pheno(smoke.add_ext(prefix, ".pheno"), y)


def run_r_baseline(method, plink_prefix, pheno_file, out_prefix):
    if RSCRIPT is None:
        return {
            "method": method,
            "method_available": False,
            "runtime_sec": 0.0,
            "n_snps": 0,
            "valid_snp_count": 0,
            "top_snp": "NA",
            "max_log10p": "NA",
            "error_message": "Rscript binary not found",
        }
    try:
        subprocess.run(
            [
                RSCRIPT,
                str(ROOT / "tests" / "psr_external_benchmark.R"),
                "ordinal_gwas",
                method,
                str(plink_prefix),
                str(pheno_file),
                str(out_prefix),
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        return {
            "method": method,
            "method_available": None,
            "runtime_sec": 300.0,
            "n_snps": 0,
            "valid_snp_count": 0,
            "top_snp": "NA",
            "max_log10p": "NA",
            "error_message": "R process timed out after 300s",
        }
    except Exception as e:
        return {
            "method": method,
            "method_available": None,
            "runtime_sec": 0.0,
            "n_snps": 0,
            "valid_snp_count": 0,
            "top_snp": "NA",
            "max_log10p": "NA",
            "error_message": str(e),
        }
    summary_path = Path(str(out_prefix) + ".r_baseline_summary.tsv")
    if summary_path.exists():
        with summary_path.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row["method_available"] = row["method_available"] == "TRUE"
                row["runtime_sec"] = float(row["runtime_sec"])
                row["n_snps"] = int(row["n_snps"])
                row["valid_snp_count"] = int(row["valid_snp_count"])
                return row
    return {
        "method": method,
        "method_available": None,
        "runtime_sec": 0.0,
        "n_snps": 0,
        "valid_snp_count": 0,
        "top_snp": "NA",
        "max_log10p": "NA",
        "error_message": "R completed but no summary output found",
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


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(BENCH_SEED)

    print(f"[PSR-GWAS v0.5] R baseline runner")
    print(f"  reps={BENCH_REPS} n={BENCH_N} m={BENCH_M} levels={BENCH_LEVELS}")
    print(f"  Rscript={RSCRIPT or 'NOT FOUND'}")
    print(f"  R methods: {', '.join(R_METHODS)}")

    # Write availability
    if RSCRIPT is not None:
        avail_out = OUT_DIR / "v05_ordinal_r_baseline_availability.tsv"
        subprocess.run(
            [RSCRIPT, str(ROOT / "tests" / "psr_external_benchmark.R"),
             "ordinal_availability", str(avail_out)],
            capture_output=True,
        )
        if avail_out.exists():
            print(f"[PSR-GWAS] R baseline availability written to {avail_out}")

    scenarios = ["null", "power_low", "power_high"]
    r_rows = []
    for scenario in scenarios:
        effect = scenario_effect(scenario)
        causal_idx = max(1, BENCH_M // 3) - 1 if effect > 0.0 else -1
        for rep in range(BENCH_REPS):
            seed = BENCH_SEED + rep + 100 * scenarios.index(scenario)
            sub_rng = np.random.default_rng(seed)
            G, y = simulate_ordinal_data(sub_rng, BENCH_N, BENCH_M, BENCH_LEVELS, effect)
            work_prefix = WORK_DIR / f"ordinal_r_base_rep{rep}_{scenario}"
            pheno_file = smoke.add_ext(work_prefix, ".pheno")
            write_plink(work_prefix, G, y, seed)
            for method in R_METHODS:
                out_prefix = WORK_DIR / f"r_{method.replace('::', '_')}_rep{rep}_{scenario}"
                result = run_r_baseline(method, work_prefix, pheno_file, out_prefix)
                result["scenario"] = scenario
                result["rep"] = rep
                result["seed"] = seed
                result["n"] = BENCH_N
                result["m"] = BENCH_M
                r_rows.append(result)

    write_tsv(OUT_DIR / "v05_ordinal_r_baseline_summary.tsv", r_rows)
    print(f"[PSR-GWAS v0.5] R baseline summary written ({len(r_rows)} rows)")


if __name__ == "__main__":
    main()
