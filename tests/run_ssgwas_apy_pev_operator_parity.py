#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

import numpy as np


DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-7


def run(binary, out, core, blend, pev_mode):
    cmd = [
        binary,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "pcg", "--ssgwas-pcg-tol", "1e-10",
        "--ssgwas-pcg-maxit", "5000",
        "--ssgwas-pvalue-method", "apy-approx",
        "--ssgwas-apy-core-file", f"{CORE_DIR}/{core}.txt",
        "--ssgwas-grm-mode", "streamed",
        "--ssgwas-apy-inverse-mode", "operator",
        "--ssgwas-mme-apy-mode", "operator",
        "--ssgwas-pev-apy-mode", pev_mode,
        "--out", f"{out}/m2",
    ]
    if blend == "a22":
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    os.makedirs(out, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n{result.stdout[-2000:]}\n{result.stderr[-2000:]}")


def table(path):
    data = np.genfromtxt(path, names=True, delimiter="\t", dtype=None, encoding="utf-8")
    return np.atleast_1d(data)


def metadata(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            if "\t" in line:
                key, value = line.rstrip("\n").split("\t", 1)
                result[key] = value
    return result


def compare_column(reference, operator, column):
    a = np.asarray(reference[column], dtype=float)
    b = np.asarray(operator[column], dtype=float)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return float("inf")
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--outdir", default="tests/results/ssgwas_apy_pev_operator_parity")
    args = parser.parse_args()
    failures = []

    for core in ["core28", "core44"]:
        for blend in ["identity", "a22"]:
            stem = f"{args.outdir}/{core}_{blend}"
            run(args.bin, f"{stem}_dense", core, blend, "dense")
            run(args.bin, f"{stem}_operator", core, blend, "operator")
            dense = table(f"{stem}_dense/m2.ss_backsolve.pvalue.txt")
            operator = table(f"{stem}_operator/m2.ss_backsolve.pvalue.txt")
            meta = metadata(f"{stem}_operator/m2.ss_backsolve.meta.txt")
            if not np.array_equal(dense["SNP"], operator["SNP"]):
                failures.append(f"{core} {blend} SNP rows/order differ")
            metrics = {
                "SE": compare_column(dense, operator, "SE_ALLELE"),
                "VAR_EFFECT": compare_column(dense, operator, "VAR_EFFECT"),
                "Z": compare_column(dense, operator, "Z"),
                "P": compare_column(dense, operator, "P"),
                "LOG10P": compare_column(dense, operator, "LOG10P"),
            }
            print(f"{core} {blend}: " + " ".join(f"{k}={v:.3e}" for k, v in metrics.items()))
            failures.extend(
                f"{core} {blend} {name}={value:.3e}"
                for name, value in metrics.items() if value >= TOL)
            for key, expected in [
                ("PEV_APY_MODE", "operator"),
                ("PEV_APY_Q_REQUIRED", "false"),
                ("PEV_APY_CORE_MATRIX_REQUIRED", "false"),
                ("FULL_Q_REQUIRED", "false"),
                ("CORE_MATRIX_REQUIRED", "false"),
                ("FULL_GCN_STORED", "false"),
            ]:
                if meta.get(key) != expected:
                    failures.append(f"{core} {blend} {key}={meta.get(key, '<missing>')}")
            residual = float(meta.get("PEV_APY_MAX_RELATIVE_RESIDUAL", "inf"))
            if not np.isfinite(residual) or residual > 1e-10:
                failures.append(f"{core} {blend} residual={residual:.3e}")

    if failures:
        print("APY PEV OPERATOR DOWNSTREAM PARITY: FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("APY PEV OPERATOR DOWNSTREAM PARITY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
