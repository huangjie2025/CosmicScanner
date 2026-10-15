#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
import subprocess

import numpy as np


def run(command):
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError("Command failed:\n" + " ".join(map(str, command)) + "\n" + completed.stdout)


def read_selected(path):
    with path.open() as source:
        return {row["SNP"]: row for row in csv.DictReader(source, delimiter="\t")}


def close(observed, expected, tolerance, label):
    if not math.isfinite(observed) or abs(observed - expected) > tolerance:
        raise RuntimeError(f"{label}: observed={observed} expected={expected}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    binary = str(Path(args.bin).resolve())
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(20260914)
    n = 500
    causal = rng.binomial(2, 0.35, n)
    proxy = causal.copy()
    changed = rng.choice(n, size=70, replace=False)
    proxy[changed] = rng.binomial(2, 0.35, len(changed))
    removed = rng.binomial(2, 0.42, n)
    independent = rng.binomial(2, 0.28, n)
    genotype = np.column_stack((causal, removed, proxy, independent))
    ids = ["causal", "removed", "proxy", "independent"]
    positions = [1000, 2000, 3000, 4000]

    prefix = out / "reference"
    with Path(str(prefix) + ".fam").open("w") as target:
        for i in range(n):
            target.write(f"F{i + 1} I{i + 1} 0 0 0 -9\n")
    with Path(str(prefix) + ".bim").open("w") as target:
        for snp, position in zip(ids, positions):
            target.write(f"1 {snp} 0 {position} A C\n")
    bed = bytearray(b"\x6c\x1b\x01")
    encoding = {2: 0b00, 1: 0b10, 0: 0b11}
    for column in genotype.T:
        packed = bytearray((n + 3) // 4)
        for i, value in enumerate(column):
            packed[i // 4] |= encoding[int(value)] << (2 * (i % 4))
        bed.extend(packed)
    Path(str(prefix) + ".bed").write_bytes(bed)

    ld = np.corrcoef(genotype, rowvar=False)
    external = out / "reference.ld.tsv"
    with external.open("w") as target:
        target.write("SNP CHR BP A1 A2 " + " ".join(ids) + "\n")
        for i, snp in enumerate(ids):
            target.write(f"{snp} 1 {positions[i]} A C " +
                         " ".join(f"{value:.17g}" for value in ld[i]) + "\n")

    selected = [0, 2, 3]
    selected_ld = ld[np.ix_(selected, selected)]
    truth = np.array([5.0, 0.0, 4.0])
    marginal_z = selected_ld @ truth
    summary = out / "summary.tsv"
    with summary.open("w") as target:
        target.write("CHR SNP POS A1 A2 BETA SE P N AF\n")
        for local, source in enumerate(selected):
            z = marginal_z[local]
            p = math.erfc(abs(z) / math.sqrt(2.0))
            target.write(f"1 {ids[source]} {positions[source]} A C {0.1 * z:.17g} 0.1 {p:.17g} {n} 0.3\n")
    condition = out / "joint.snps"
    condition.write_text("causal\nproxy\nindependent\n", encoding="ascii")

    common = ["--cojo-joint", "--cojo-summary", str(summary), "--cojo-cond", str(condition),
              "--cojo-p", "0.5", "--maf", "0"]
    derived_prefix = out / "derived"
    supplied_prefix = out / "supplied"
    run([binary, "--bfile", str(prefix), *common, "--out", str(derived_prefix)])
    run([binary, "--cojo-ld-ref", str(external), *common, "--out", str(supplied_prefix)])

    derived = read_selected(Path(str(derived_prefix) + ".jma"))
    supplied = read_selected(Path(str(supplied_prefix) + ".jma"))
    solved = np.linalg.solve(selected_ld, marginal_z)
    covariance = np.linalg.inv(selected_ld)
    for i, snp in enumerate(("causal", "proxy", "independent")):
        expected_beta = 0.1 * solved[i]
        expected_se = 0.1 * math.sqrt(covariance[i, i])
        for label, rows in (("genotype", derived), ("external", supplied)):
            close(float(rows[snp]["BETA"]), expected_beta, 2e-8, f"{label} {snp} beta")
            close(float(rows[snp]["SE"]), expected_se, 2e-8, f"{label} {snp} SE")
        close(float(derived[snp]["BETA"]), float(supplied[snp]["BETA"]), 2e-8,
              f"derived/external {snp} beta")
        close(float(derived[snp]["SE"]), float(supplied[snp]["SE"]), 2e-8,
              f"derived/external {snp} SE")

    print("PASS COJO genotype-derived/external LD parity and NumPy multiple-regression oracle")


if __name__ == "__main__":
    main()
