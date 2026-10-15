#!/usr/bin/env python3
"""CLI regression: additive_dominance fail-fast guards (Phase 4.1).

Verifies that unsupported option combinations for --gwas-model additive_dominance
fail fast with a clear error instead of being silently ignored, and that the
supported AD invocation still runs.
"""
import math
import os
import subprocess
import sys

N = 40
M = 30


def next_double(s):
    s = (s * 1103515245 + 12345) & 0x7FFFFFFF
    return s / 0x7FFFFFFF, s


def generate_plink(workdir, seed=7):
    os.makedirs(workdir, exist_ok=True)
    s = seed
    n, m = N, M
    prefix = os.path.join(workdir, "sim")
    with open(prefix + ".bim", "w") as f:
        for j in range(m):
            a1 = "C" if (j % 2) else "T"
            a2 = "T" if a1 == "C" else "C"
            f.write(f"{j % 4 + 1}\tsnp{j}\t0\t{j * 500}\t{a1}\t{a2}\n")
    genotypes = [[0] * n for _ in range(m)]
    for j in range(m):
        p = 0.3 + 0.3 * ((j * 37) % 100) / 100.0
        for i in range(n):
            u, s = next_double(s)
            g = 0
            if u < p * p:
                g = 2
            elif u < p * p + 2 * p * (1 - p):
                g = 1
            genotypes[j][i] = g
    het = [1.0 if genotypes[5][i] == 1 else 0.0 for i in range(n)]
    y = [1.2 * het[i] + ((u := next_double(s)[0]) - 0.5) for i in range(n)]
    with open(os.path.join(workdir, "pheno.txt"), "w") as f:
        for i in range(n):
            f.write(f"F{i}\tI{i}\t{y[i]:.6f}\n")
    with open(prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"F{i}\tI{i}\t0\t0\t0\t1\n")
    with open(prefix + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01")
        byte_len = (n + 3) // 4
        for j in range(m):
            b = bytearray(byte_len)
            for i in range(n):
                code = {0: 0, 1: 2, 2: 3}[genotypes[j][i]]
                b[i >> 2] |= code << ((i & 3) * 2)
            f.write(bytes(b))
    return prefix


def run(scanner, prefix, out, extra):
    cmd = [scanner, "--bfile", prefix, "--out", out, "--gwas",
           "--gwas-model", "additive_dominance", "--no-loco", "--maf", "0.01"]
    cmd += list(extra)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--scanner", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    prefix = generate_plink(args.workdir)

    failures = []

    def expect_fail(extra, needle, label):
        r = run(args.scanner, prefix, os.path.join(args.workdir, label), extra)
        combined = (r.stdout or "") + (r.stderr or "")
        if r.returncode == 0:
            failures.append(f"{label}: expected non-zero exit, got {r.returncode}")
        elif needle not in combined:
            failures.append(f"{label}: error message missing '{needle}' "
                            f"(rc={r.returncode})\nstdout={r.stdout}\nstderr={r.stderr}")
        else:
            print(f"  PASS {label}: failed fast with '{needle}'")

    expect_fail(["--no-covar-proj"], "not compatible with --no-covar-proj", "no_covar_proj")
    expect_fail(["--no-covar-proj"], "requires covariate projection", "no_covar_proj_msg")
    expect_fail(["--out-gz"], "does not yet support --out-gz", "out_gz")
    expect_fail(["--out-ssf"], "does not yet support --out-ssf", "out_ssf")
    expect_fail(["--manhattan-tsv"], "--manhattan-tsv/--qq-tsv", "manhattan")
    expect_fail(["--qq-tsv"], "--manhattan-tsv/--qq-tsv", "qq")
    expect_fail(["--fastlmm"], "--fastlmm/--lmm-low-rank", "fastlmm")
    expect_fail(["--lmm-low-rank", "4"], "--fastlmm/--lmm-low-rank", "lmm_low_rank")

    r = run(args.scanner, prefix, os.path.join(args.workdir, "ok"), [])
    if r.returncode != 0:
        failures.append(f"baseline additive_dominance run failed: rc={r.returncode}\n"
                        f"stdout={r.stdout}\nstderr={r.stderr}")
    elif not os.path.exists(os.path.join(args.workdir, "ok.assoc.ad.txt")):
        failures.append("baseline additive_dominance run produced no .assoc.ad.txt")
    else:
        print("  PASS baseline: supported additive_dominance invocation runs")

    if failures:
        print("\n".join(failures))
        return 1
    print("All AD CLI guard tests PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())