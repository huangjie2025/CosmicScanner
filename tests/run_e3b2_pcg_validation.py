#!/usr/bin/env python3
"""Phase E3-B2 validation: Compare Jacobi vs LowRank preconditioners on PCG solver."""
import subprocess, sys, os
import numpy as np

BIN = "build/cosmicscanner"
DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"

def run_cs(out, core_file, preconditioner, rank=50):
    cmd = [
        BIN,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "ldlt",
        "--ssgwas-pvalue-method", "apy-exact",
        "--ssgwas-apy-core-file", core_file,
        "--ssgwas-apy-seed", "123",
        "--ssgwas-gcc-solver", "pcg",
        "--ssgwas-pcg-preconditioner", preconditioner,
        "--ssgwas-pcg-preconditioner-rank", str(rank),
        "--out", f"{out}/m2",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"run failed:\n{r.stderr[-1200:]}\n{r.stdout[-600:]}")
    return r

def load_vec(path, col=1):
    return np.array([float(l.split("\t")[col])
                     for l in open(path).read().splitlines()[1:]])

def maxdiff(a, b):
    return float(np.max(np.abs(a - b))) if len(a) == len(b) else float("inf")

def main():
     os.makedirs("tests/results/ssgwas_e3b2_precond_parity", exist_ok=True)

     core_file = f"{CORE_DIR}/core44.txt"

    print("Running Jacobi PCG...")
     run_cs("tests/results/ssgwas_e3b2_precond_parity/jacobi", core_file, "jacobi")

    print("Running LowRank PCG (rank=50)...")
     run_cs("tests/results/ssgwas_e3b2_precond_parity/lowrank50", core_file, "lowrank", 50)

    # Compare results
    ebv_j = load_vec("tests/results/ssgwas_e3b2_precond_parity/jacobi/m2.ss_backsolve.ebv.txt")
    ebv_lr = load_vec("tests/results/ssgwas_e3b2_precond_parity/lowrank50/m2.ss_backsolve.ebv.txt")
    d_ebv = maxdiff(ebv_j, ebv_lr)
     print(f"EBV max|diff| (Jacobi vs LowRank50) = {d_ebv:.3e}")

    for col, name in [(6, "EA"), (7, "SE"), (9, "Z"), (11, "LOG10P")]:
        p_j = load_vec("tests/results/ssgwas_e3b2_precond_parity/jacobi/m2.ss_backsolve.pvalue.txt", col)
        p_lr = load_vec("tests/results/ssgwas_e3b2_precond_parity/lowrank50/m2.ss_backsolve.pvalue.txt", col)
        dp = maxdiff(p_j, p_lr)
         print(f"  {name} max|diff| = {dp:.3e}")

    # Also compare with dense
    print("\nRunning Dense (reference)...")
     run_cs("tests/results/ssgwas_e3b2_precond_parity/dense", core_file, "dense")  # This will fail, dense doesn't use pcg

    # Instead, just run without --ssgwas-gcc-solver to get dense
    cmd_dense = [
        BIN,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "ldlt",
        "--ssgwas-pvalue-method", "apy-exact",
        "--ssgwas-apy-core-file", core_file,
        "--ssgwas-apy-seed", "123",
        "--out", "tests/results/ssgwas_e3b2_precond_parity/dense/m2",
    ]
    r = subprocess.run(cmd_dense, capture_output=True, text=True)
    if r.returncode != 0:
         raise RuntimeError(f"dense run failed:\n{r.stderr[-1200:]}")

    ebv_d = load_vec("tests/results/ssgwas_e3b2_precond_parity/dense/m2.ss_backsolve.ebv.txt")
    d_ebv_dense = maxdiff(ebv_d, ebv_j)
    print(f"EBV max|diff| (Dense vs Jacobi) = {d_ebv_dense:.3e}")
    d_ebv_dense_lr = maxdiff(ebv_d, ebv_lr)
    print(f"EBV max|diff| (Dense vs LowRank50) = {d_ebv_dense_lr:.3e}")

if __name__ == "__main__":
    main()
