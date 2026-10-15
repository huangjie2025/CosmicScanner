#!/usr/bin/env python3
import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np


DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-11


def make_edge_fixture(outdir):
    os.makedirs(outdir, exist_ok=True)
    source = f"{DATASET}/sim"
    target = f"{outdir}/sim_edge"
    for suffix in [".bim", ".fam"]:
        shutil.copyfile(source + suffix, target + suffix)
    shutil.copyfile(f"{DATASET}/sim.pheno", target + ".pheno")
    shutil.copyfile(f"{DATASET}/sim.ped", target + ".ped")
    bed = bytearray(open(source + ".bed", "rb").read())
    n_samples = sum(1 for _ in open(source + ".fam"))
    bytes_per_snp = (n_samples + 3) // 4
    bed[3:3 + bytes_per_snp] = bytes([0xFF] * bytes_per_snp)
    second = 3 + bytes_per_snp
    for sample in range(0, n_samples, 5):
        byte_index = second + sample // 4
        shift = 2 * (sample % 4)
        bed[byte_index] = (bed[byte_index] & ~(3 << shift)) | (1 << shift)
    with open(target + ".bed", "wb") as handle:
        handle.write(bed)
    return target


def run(binary, out, bfile, pheno, ped, scan_mode, block_size, case):
    cmd = [
        binary,
        "--bfile", bfile,
        "--pheno", pheno,
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", ped,
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-pvalue-method", "apy-approx",
        "--ssgwas-grm-mode", "streamed",
        "--ssgwas-pev-scan-mode", scan_mode,
        "--ssgwas-pev-scan-block-size", str(block_size),
        "--out", f"{out}/m2",
    ]
    if case["random_core"]:
        cmd += ["--ssgwas-apy-core-size", "14", "--ssgwas-apy-seed", "2468"]
    else:
        cmd += ["--ssgwas-apy-core-file", f"{CORE_DIR}/{case['core']}.txt"]
    if case["operator"]:
        cmd += [
            "--ssgwas-solver", "pcg", "--ssgwas-pcg-tol", "1e-10",
            "--ssgwas-pcg-maxit", "5000",
            "--ssgwas-apy-inverse-mode", "operator",
            "--ssgwas-mme-apy-mode", "operator",
            "--ssgwas-pev-apy-mode", "operator",
            "--ssgwas-hinv-mode", "operator",
        ]
    else:
        cmd += ["--ssgwas-solver", "ldlt"]
    if case["a22"]:
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    os.makedirs(out, exist_ok=True)
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\n{result.stdout[-1800:]}\n{result.stderr[-1800:]}")


def table(path):
    return np.atleast_1d(np.genfromtxt(
        path, names=True, delimiter="\t", dtype=None, encoding="utf-8"))


def metadata(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            if "\t" in line:
                key, value = line.rstrip("\n").split("\t", 1)
                result[key] = value
    return result


def compare_numeric(reference, candidate, column):
    def numeric(values):
        return np.asarray([
            float(value) if str(value).upper() not in {"NA", "NAN"} else np.nan
            for value in np.atleast_1d(values)
        ], dtype=float)
    a = numeric(reference[column])
    b = numeric(candidate[column])
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        return math.inf
    finite = np.isfinite(a) & np.isfinite(b)
    return float(np.max(np.abs(a[finite] - b[finite]))) if np.any(finite) else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--outdir", default="tests/results/ssgwas_pev_scan_parity")
    args = parser.parse_args()
    edge = make_edge_fixture(f"{args.outdir}/edge_data")
    cases = [
        {"name": "dense_identity_fixed", "operator": False, "a22": False,
         "random_core": False, "core": "core_small", "block": 1, "edge": False},
        {"name": "dense_a22_random", "operator": False, "a22": True,
         "random_core": True, "core": "", "block": 7, "edge": False},
        {"name": "operator_identity_edge", "operator": True, "a22": False,
         "random_core": False, "core": "core_small", "block": 64, "edge": True},
        {"name": "operator_a22_fixed", "operator": True, "a22": True,
         "random_core": False, "core": "core_medium", "block": 256, "edge": False},
    ]
    failures = []
    for case in cases:
        if case["edge"]:
            bfile, pheno, ped = edge, edge + ".pheno", edge + ".ped"
        else:
            bfile = f"{DATASET}/sim"
            pheno = f"{DATASET}/sim.pheno"
            ped = f"{DATASET}/sim.ped"
        scalar_dir = f"{args.outdir}/{case['name']}_scalar"
        batched_dir = f"{args.outdir}/{case['name']}_batched"
        run(args.bin, scalar_dir, bfile, pheno, ped, "scalar", case["block"], case)
        run(args.bin, batched_dir, bfile, pheno, ped, "batched", case["block"], case)
        scalar = table(f"{scalar_dir}/m2.ss_backsolve.pvalue.txt")
        batched = table(f"{batched_dir}/m2.ss_backsolve.pvalue.txt")
        scalar_meta = metadata(f"{scalar_dir}/m2.ss_backsolve.meta.txt")
        batched_meta = metadata(f"{batched_dir}/m2.ss_backsolve.meta.txt")
        metrics = {column: compare_numeric(scalar, batched, column) for column in [
            "EFFECT_ALLELE", "SE_ALLELE", "VAR_EFFECT", "Z", "P", "LOG10P"]}
        print(case["name"] + ": " + " ".join(
            f"{column}={value:.3e}" for column, value in metrics.items()))
        failures.extend(
            f"{case['name']} {column}={value:.3e}"
            for column, value in metrics.items() if value > TOL)
        if not np.array_equal(scalar["STATUS"], batched["STATUS"]):
            failures.append(f"{case['name']} STATUS mismatch")
        expected_blocks = math.ceil(int(batched_meta["PEV_SCAN_VARIANTS_DECODED"]) /
                                    int(batched_meta["PEV_SCAN_BLOCK_SIZE"]))
        gates = [
            (batched_meta.get("PEV_SCAN_MODE"), "batched", "mode"),
            (batched_meta.get("PEV_SCAN_MULTIPLY_Z_V_CALLS"), "0", "multiply calls"),
            (batched_meta.get("PEV_SCAN_GENOTYPE_PASSES"), "1", "genotype passes"),
            (batched_meta.get("PEV_SCAN_BLOCKS"), str(expected_blocks), "block count"),
            (batched_meta.get("FULL_MXM_REQUIRED"), "false", "M-by-M gate"),
        ]
        for actual, expected, label in gates:
            if actual != expected:
                failures.append(f"{case['name']} {label}: {actual} != {expected}")
        if int(batched_meta["PEV_SCAN_BED_BYTES_READ_ESTIMATE"]) >= int(
                scalar_meta["PEV_SCAN_BED_BYTES_READ_ESTIMATE"]):
            failures.append(f"{case['name']} batched BED bytes were not reduced")
        if case["edge"]:
            statuses = set(np.asarray(batched["STATUS"], dtype=str))
            if "NA_MONOMORPHIC" not in statuses:
                failures.append(f"{case['name']} monomorphic status was not preserved")
            if int(batched_meta.get("PEV_SCAN_MISSING_TO_ZERO", "0")) <= 0:
                failures.append(f"{case['name']} missing-to-zero path was not exercised")

    if failures:
        print("PEV SCAN SCALAR/BATCHED PARITY: FAIL")
        for failure in failures:
            print(f"  {failure}")
        return 1
    print("PEV SCAN SCALAR/BATCHED PARITY: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
