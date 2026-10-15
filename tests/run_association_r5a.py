#!/usr/bin/env python3
"""Native SAIGE-style binary mixed score/SPA CLI smoke and benchmark."""

import argparse
import csv
import math
from pathlib import Path
import random
import sys
import tempfile

from run_association_r2 import Harness, fixture, require, table


CASES = ("saige_style", "benchmark")
HEADER = "MODEL CHR SNP BP A1 A2 FREQ N BETA SE P ITERATION PSEUDO_QTN"


def binary_phenotype(data, path):
    rng = random.Random(20260910)
    means = []
    for marker in (0, 3):
        observed = [value for value in data["calls"][marker] if value is not None]
        means.append(math.fsum(observed) / len(observed))
    outcomes = []
    for i in range(data["n"]):
        lead = means[0] if data["calls"][0][i] is None else data["calls"][0][i]
        secondary = means[1] if data["calls"][3][i] is None else data["calls"][3][i]
        probability = 1.0 / (1.0 + math.exp(-(-1.4 + 1.1 * lead + 0.35 * secondary)))
        outcomes.append(int(rng.random() < probability))
    require(0 in outcomes and 1 in outcomes, "SAIGE-style fixture needs cases and controls")
    path.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {outcomes[i]}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")


def validate_output(path, n):
    rows = table(path, HEADER)
    require(all(row["MODEL"] == "SAIGE_STYLE" for row in rows), "Wrong SAIGE-style MODEL")
    require(len(rows) == len({row["SNP"] for row in rows}), "Duplicate SAIGE-style SNP rows")
    for row in rows:
        require(int(row["N"]) == n, f"Wrong N for {row['SNP']}")
        require(math.isfinite(float(row["BETA"])), f"Nonfinite beta for {row['SNP']}")
        require(math.isfinite(float(row["SE"])) and float(row["SE"]) > 0,
                f"Invalid SE for {row['SNP']}")
        p = float(row["P"])
        require(math.isfinite(p) and 0 <= p <= 1, f"Invalid P for {row['SNP']}")
        require(int(row["ITERATION"]) >= 1 and row["PSEUDO_QTN"] == "0",
                f"Invalid null iteration metadata for {row['SNP']}")
    return rows


def saige_case(h, failures=True):
    phenotype = h.root / "binary.txt"
    binary_phenotype(h.data, phenotype)
    prefix = h.run("saige_style", ["--saige-style", "--biobank-mode",
                                    "--pheno", phenotype, "--sparse-cutoff", "0.05"])
    rows = validate_output(Path(str(prefix) + ".cosmic.assoc"), h.data["n"])
    require(len(rows) >= int(0.9 * h.data["m"]), "Too many SAIGE-style variants were filtered")
    metadata = Path(str(prefix) + ".saige_style.meta.tsv").read_text(encoding="ascii")
    require("model\tSAIGE_STYLE" in metadata and "status\texperimental" in metadata and
            "null_solver\tsparse_ldlt_pql" in metadata and
            "variance_component_method\tmoment_he" in metadata and
            "null_converged\t1" in metadata, "Missing SAIGE-style null metadata")
    fields = dict(line.split("\t", 1) for line in metadata.splitlines()[1:])
    require(float(fields["variance_ratio"]) > 0 and int(fields["variance_ratio_markers"]) > 0,
            "Invalid SAIGE-style variance-ratio calibration")
    require(int(fields["sparse_grm_nnz"]) >= h.data["n"], "Invalid sparse GRM metadata")
    if failures:
        h.run("saige_without_biobank", ["--saige-style", "--pheno", phenotype],
              error=r"requires --biobank-mode")
        invalid = h.root / "nonbinary.txt"
        invalid.write_text("".join(
            f"{fid} {iid} {2 if i == 0 else i % 2}\n"
            for i, (fid, iid) in enumerate(h.data["samples"])), encoding="ascii")
        h.run("saige_nonbinary", ["--saige-style", "--biobank-mode", "--pheno", invalid],
              error=r"phenotype.*0/1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", required=True, type=Path)
    parser.add_argument("--case", required=True, choices=CASES)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    binary = args.bin.expanduser().resolve()
    require(binary.is_file(), f"Missing executable: {binary}")
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"association_r5a_{args.case}_", dir=out) as directory:
        root = Path(directory)
        data = fixture(root, 120, 300) if args.case == "benchmark" else fixture(root, 80, 80)
        harness = Harness(binary, root, data)
        saige_case(harness, failures=args.case != "benchmark")
        if args.case == "benchmark":
            report = out / "benchmark.time.tsv"
            with report.open("w", encoding="ascii", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("command", "n", "m", "seconds", "returncode"))
                writer.writerows((name, n, m, f"{elapsed:.6f}", rc)
                                 for name, n, m, elapsed, rc in harness.timings)
            print(f"Timing report (no threshold): {report}")
        print(f"PASS {args.case}: {len(harness.timings)} commands, n={data['n']}, m={data['m']}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError) as error:
        print(f"FAIL association_r5a: {error}", file=sys.stderr)
        sys.exit(1)
