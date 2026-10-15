#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

import numpy as np


DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-7


def run(binary, out, core, blend, mme_mode):
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
        "--ssgwas-mme-apy-mode", mme_mode,
        "--out", f"{out}/m2",
    ]
    if blend == "a22":
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    os.makedirs(out, exist_ok=True)
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
    parser.add_argument("--outdir", default="tests/results/ssgwas_apy_mme_operator_parity")
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
            meta = metadata(f"{stem}_operator/m2.ss_backsolve.meta.txt")
            metrics = {
                "EBV": compare_column(ebv_dense, ebv_operator, "EBV"),
                "SNP_EFFECT": compare_column(snp_dense, snp_operator, "EFFECT_ALLELE"),
                "SE": compare_column(p_dense, p_operator, "SE_ALLELE"),
                "Z": compare_column(p_dense, p_operator, "Z"),
                "P": compare_column(p_dense, p_operator, "P"),
                "LOG10P": compare_column(p_dense, p_operator, "LOG10P"),
            }
            print(f"{core} {blend}: " + " ".join(f"{k}={v:.3e}" for k, v in metrics.items()))
            failures.extend(
                f"{core} {blend} {name}={value:.3e}"
                for name, value in metrics.items() if value >= TOL)
            for key, expected in [
                ("MME_APY_MODE", "operator"),
                ("MME_APY_Q_REQUIRED", "false"),
                ("MME_APY_CORE_BLOCK_REQUIRED", "false"),
                ("FULL_GCN_STORED", "false"),
            ]:
                if meta.get(key) != expected:
                    failures.append(f"{core} {blend} {key}={meta.get(key, '<missing>')}")

    if failures:
        print("APY MME OPERATOR DOWNSTREAM PARITY: FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("APY MME OPERATOR DOWNSTREAM PARITY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
