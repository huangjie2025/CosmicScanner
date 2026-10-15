#!/usr/bin/env python3
"""Native Firth logistic and conditional-joint CLI smoke/benchmark."""

import argparse
import csv
import math
from pathlib import Path
import random
import sys
import tempfile

from run_association_r2 import Harness, fixture, require, table


CASES = ("firth", "conditional_joint", "benchmark")
HEADER = "MODEL CHR SNP BP A1 A2 FREQ N BETA SE P ITERATION PSEUDO_QTN"


def binary_phenotype(data, path, separated):
    rng = random.Random(20260909)
    observed = [value for value in data["calls"][0] if value is not None]
    mean = math.fsum(observed) / len(observed)
    outcomes = []
    for value in data["calls"][0]:
        dosage = mean if value is None else value
        if separated:
            outcomes.append(int(dosage >= 1.0))
        else:
            probability = 1.0 / (1.0 + math.exp(-(-0.7 + 0.65 * dosage)))
            outcomes.append(int(rng.random() < probability))
    require(0 in outcomes and 1 in outcomes, "Binary fixture needs cases and controls")
    path.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {outcomes[i]}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")
    return outcomes


def joint_phenotype(data, path):
    columns = []
    for marker in (0, 3):
        observed = [value for value in data["calls"][marker] if value is not None]
        mean = math.fsum(observed) / len(observed)
        columns.append([mean if value is None else value for value in data["calls"][marker]])
    outcomes = [1.5 + 1.8 * columns[0][i] + 1.1 * columns[1][i] +
                0.08 * math.sin(i) for i in range(data["n"])]
    path.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {outcomes[i]:.12g}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")


def validate_output(path, model, n):
    rows = table(path, HEADER)
    require(all(row["MODEL"] == model for row in rows), f"Wrong MODEL in {path}")
    require(len({row["SNP"] for row in rows}) == len(rows), f"Duplicate SNP rows in {path}")
    for row in rows:
        require(int(row["N"]) == n, f"Wrong sample count for {row['SNP']}")
        require(math.isfinite(float(row["BETA"])), f"Nonfinite beta for {row['SNP']}")
        require(math.isfinite(float(row["SE"])) and float(row["SE"]) > 0,
                f"Invalid SE for {row['SNP']}")
        p = float(row["P"])
        require(math.isfinite(p) and 0 <= p <= 1, f"Invalid P for {row['SNP']}")
        require(int(row["ITERATION"]) >= 1 and row["PSEUDO_QTN"] in ("0", "1"),
                f"Invalid iteration/selection metadata for {row['SNP']}")
    return rows


def firth_case(h, failures=True):
    for separated in (False, True):
        kind = "separation" if separated else "normal"
        phenotype = h.root / f"{kind}_binary.txt"
        binary_phenotype(h.data, phenotype, separated)
        prefix = h.run(f"firth_{kind}", ["--firth-logistic", "--pheno", phenotype])
        rows = validate_output(Path(str(prefix) + ".cosmic.assoc"),
                               "FIRTH_LOGISTIC", h.data["n"])
        by_id = {row["SNP"]: row for row in rows}
        require("rs0001" in by_id, "Firth output omitted the simulated lead")
        if separated:
            require(float(by_id["rs0001"]["BETA"]) > 0,
                    "Firth separation estimate has the wrong direction")
        metadata = Path(str(prefix) + ".generalized.meta.tsv").read_text(encoding="ascii")
        require("model\tFIRTH_LOGISTIC" in metadata and "status\texperimental" in metadata,
                "Missing Firth metadata")
    if failures:
        invalid = h.root / "invalid_binary.txt"
        invalid.write_text("".join(
            f"{fid} {iid} {2 if i == 0 else i % 2}\n"
            for i, (fid, iid) in enumerate(h.data["samples"])), encoding="ascii")
        h.run("invalid_firth", ["--firth-logistic", "--pheno", invalid],
              error=r"Firth logistic phenotype.*0/1")


def conditional_joint_case(h, failures=True):
    condition_file = h.root / "joint_condition.txt"
    condition_file.write_text("rs0001\nrs0002\nrs0004\n", encoding="ascii")
    phenotype = h.root / "joint_phenotype.txt"
    joint_phenotype(h.data, phenotype)
    prefix = h.run("conditional_joint", ["--condition-snp", condition_file,
                                          "--conditional-joint", "--pheno", phenotype])
    rows = validate_output(Path(str(prefix) + ".conditional_joint.assoc"),
                           "conditional_joint", h.data["n"])
    selected = {row["SNP"] for row in rows if row["PSEUDO_QTN"] == "1"}
    require(selected == {"rs0001", "rs0004"},
            f"Unexpected LD/BIC conditional-joint selection: {sorted(selected)}")
    metadata = Path(str(prefix) + ".multilocus.meta.tsv").read_text(encoding="ascii")
    require("model\tconditional_joint" in metadata and "condition_count\t3" in metadata and
            "selected_pseudo_qtn\t2" in metadata, "Missing conditional-joint metadata")
    if failures:
        h.run("joint_without_leads", ["--conditional-joint", "--pheno", h.data["pheno"]],
              error=r"requires --condition-snp")


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
    with tempfile.TemporaryDirectory(prefix=f"association_r4a_{args.case}_", dir=out) as directory:
        root = Path(directory)
        data = fixture(root, 120, 300) if args.case == "benchmark" else fixture(root, 60, 48)
        harness = Harness(binary, root, data)
        if args.case == "firth":
            firth_case(harness)
        elif args.case == "conditional_joint":
            conditional_joint_case(harness)
        else:
            firth_case(harness, failures=False)
            conditional_joint_case(harness, failures=False)
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
        print(f"FAIL association_r4a: {error}", file=sys.stderr)
        sys.exit(1)
