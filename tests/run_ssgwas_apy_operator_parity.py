#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

import numpy as np


DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-9


def run(binary, out, core, blend, inverse_mode):
    cmd = [
        binary,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "ldlt",
        "--ssgwas-pvalue-method", "apy-exact",
        "--ssgwas-apy-core-file", f"{CORE_DIR}/{core}.txt",
        "--ssgwas-grm-mode", "streamed",
        "--ssgwas-apy-inverse-mode", inverse_mode,
        "--out", f"{out}/m2",
    ]
    if blend == "a22":
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    os.makedirs(out, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n{result.stdout[-1200:]}\n{result.stderr[-1200:]}")


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
    parser.add_argument("--outdir", default="tests/results/ssgwas_apy_operator_parity")
    args = parser.parse_args()
    failures = []

    for core in ["core28", "core44"]:
        for blend in ["identity", "a22"]:
            stem = f"{args.outdir}/{core}_{blend}"
            run(args.bin, f"{stem}_dense", core, blend, "dense")
            run(args.bin, f"{stem}_operator", core, blend, "operator")
            ebv_dense = table(f"{stem}_dense/m2.ss_backsolve.ebv.txt")
            ebv_operator = table(f"{stem}_operator/m2.ss_backsolve.ebv.txt")
            p_dense = table(f"{stem}_dense/m2.ss_backsolve.pvalue.txt")
            p_operator = table(f"{stem}_operator/m2.ss_backsolve.pvalue.txt")
            snp_dense = table(f"{stem}_dense/m2.ss_backsolve.snp.txt")
            snp_operator = table(f"{stem}_operator/m2.ss_backsolve.snp.txt")
            operator_meta = metadata(f"{stem}_operator/m2.ss_backsolve.meta.txt")
            metrics = {
                "EBV": compare_column(ebv_dense, ebv_operator, "EBV"),
                "SNP_EFFECT": compare_column(snp_dense, snp_operator, "EFFECT_ALLELE"),
                "SE": compare_column(p_dense, p_operator, "SE_ALLELE"),
                "Z": compare_column(p_dense, p_operator, "Z"),
                "P": compare_column(p_dense, p_operator, "P"),
                "LOG10P": compare_column(p_dense, p_operator, "LOG10P"),
            }
            print(f"{core} {blend}: " + " ".join(f"{k}={v:.3e}" for k, v in metrics.items()))
            for name, value in metrics.items():
                if value >= TOL:
                    failures.append(f"{core} {blend} {name}={value:.3e}")
            if operator_meta.get("FULL_GCN_STORED") != "false":
                failures.append(
                    f"{core} {blend} FULL_GCN_STORED={operator_meta.get('FULL_GCN_STORED', '<missing>')}")

    if failures:
        print("APY OPERATOR DOWNSTREAM PARITY: FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("APY OPERATOR DOWNSTREAM PARITY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
