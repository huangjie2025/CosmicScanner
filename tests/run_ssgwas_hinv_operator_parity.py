#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

import numpy as np


DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-7


def run(binary, out, core, hinv_mode):
    os.makedirs(out, exist_ok=True)
    cmd = [
        binary,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--h2-tau", "1.2", "--h2-omega", "0.7",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "pcg", "--ssgwas-pcg-tol", "1e-12",
        "--ssgwas-pcg-maxit", "5000",
        "--ssgwas-pvalue-method", "apy-approx",
        "--ssgwas-apy-core-file", f"{CORE_DIR}/{core}.txt",
        "--ssgwas-grm-mode", "streamed",
        "--ssgwas-apy-inverse-mode", "operator",
        "--ssgwas-mme-apy-mode", "operator",
        "--ssgwas-hinv-mode", hinv_mode,
        "--out", f"{out}/m2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n{result.stdout[-1600:]}\n{result.stderr[-1600:]}")


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
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--outdir", default="tests/results/ssgwas_hinv_operator_parity")
    args = parser.parse_args()
    failures = []

    for core in ["core28", "core44"]:
        stem = f"{args.outdir}/{core}"
        run(args.bin, f"{stem}_dense", core, "dense")
        run(args.bin, f"{stem}_operator", core, "operator")
        ebv_dense = table(f"{stem}_dense/m2.ss_backsolve.ebv.txt")
        ebv_operator = table(f"{stem}_operator/m2.ss_backsolve.ebv.txt")
        p_dense = table(f"{stem}_dense/m2.ss_backsolve.pvalue.txt")
        p_operator = table(f"{stem}_operator/m2.ss_backsolve.pvalue.txt")
        snp_dense = table(f"{stem}_dense/m2.ss_backsolve.snp.txt")
        snp_operator = table(f"{stem}_operator/m2.ss_backsolve.snp.txt")
        meta = metadata(f"{stem}_operator/m2.ss_backsolve.meta.txt")
        metrics = {
            "EBV": compare_column(ebv_dense, ebv_operator, "EBV"),
            "SNP_EFFECT": compare_column(snp_dense, snp_operator, "EFFECT_ALLELE"),
            "SE": compare_column(p_dense, p_operator, "SE_ALLELE"),
            "Z": compare_column(p_dense, p_operator, "Z"),
            "P": compare_column(p_dense, p_operator, "P"),
            "LOG10P": compare_column(p_dense, p_operator, "LOG10P"),
        }
        print(f"{core}: " + " ".join(f"{key}={value:.3e}" for key, value in metrics.items()))
        failures.extend(
            f"{core} {name}={value:.3e}" for name, value in metrics.items() if value >= TOL)
        for reference, operator, column in [
            (ebv_dense, ebv_operator, "EBV"),
            (snp_dense, snp_operator, "EFFECT_ALLELE"),
            (p_dense, p_operator, "LOG10P"),
        ]:
            corr = float(np.corrcoef(reference[column], operator[column])[0, 1])
            if not np.isfinite(corr) or corr < 1.0 - 1e-10:
                failures.append(f"{core} {column} correlation={corr:.12g}")
        for key, expected in [
            ("HINV_MODE", "operator"),
            ("HINV_MATERIALIZED", "false"),
            ("A22_MATERIALIZED", "false"),
            ("A22_INV_MATERIALIZED", "false"),
            ("UPSTREAM_DENSE_NG2_PRESENT", "false"),
        ]:
            if meta.get(key) != expected:
                failures.append(f"{core} {key}={meta.get(key, '<missing>')}")
        if int(meta.get("A22_OPERATOR_SOLVE_COUNT", "0")) <= 0:
            failures.append(f"{core} A22_OPERATOR_SOLVE_COUNT is not positive")
        if float(meta.get("A22_OPERATOR_MAX_RELATIVE_RESIDUAL", "inf")) >= 1e-10:
            failures.append(f"{core} A22 operator residual exceeds tolerance")

    if failures:
        print("A22/HINV OPERATOR DOWNSTREAM PARITY: FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("A22/HINV OPERATOR DOWNSTREAM PARITY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
