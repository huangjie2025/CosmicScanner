#!/usr/bin/env python3
"""Native FarmCPU, BLINK, and conditional-association CLI smoke/benchmark."""

import argparse
import csv
import math
from pathlib import Path
import sys
import tempfile

from run_association_r2 import Harness, fixture, require, table


CASES = ("farmcpu", "blink", "conditional", "benchmark")
HEADER = "MODEL CHR SNP BP A1 A2 FREQ N BETA SE P ITERATION PSEUDO_QTN"


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


def multilocus_case(h, model, failures=True):
    flag = "--" + model
    extra = ["--farmcpu-max-iterations", "4", "--farmcpu-pseudo-qtn", "6"] if model == "farmcpu" else []
    prefix = h.run(model, [flag, "--pheno", h.data["pheno"], *extra])
    rows = validate_output(Path(str(prefix) + ".cosmic.assoc"), model, h.data["n"])
    ranked = sorted(rows, key=lambda row: float(row["P"]))
    require({"rs0001", "rs0002"}.intersection(row["SNP"] for row in ranked[:5]),
            f"{model} did not rank the simulated lead/proxy near the top")
    metadata = Path(str(prefix) + ".multilocus.meta.tsv").read_text(encoding="ascii")
    require(f"model\t{model}" in metadata and "status\texperimental" in metadata,
            f"Missing {model} metadata")
    if failures:
        h.run("conflicting_multilocus", ["--farmcpu", "--blink", "--pheno", h.data["pheno"]],
              error=r"choose one|multiple")
        if model == "farmcpu":
            h.run("bad_iterations", [flag, "--farmcpu-max-iterations", "0", "--pheno", h.data["pheno"]],
                  error=r"iterations")


def impute(values):
    mean = math.fsum(value for value in values if value is not None) / sum(value is not None for value in values)
    return [mean if value is None else value for value in values]


def conditional_beta(y, condition, marker):
    ym = math.fsum(y) / len(y)
    xm = math.fsum(condition) / len(condition)
    zm = math.fsum(marker) / len(marker)
    yc = [value - ym for value in y]
    xc = [value - xm for value in condition]
    zc = [value - zm for value in marker]
    xx = math.fsum(value * value for value in xc)
    y_on_x = math.fsum(a * b for a, b in zip(yc, xc)) / xx
    z_on_x = math.fsum(a * b for a, b in zip(zc, xc)) / xx
    yr = [a - y_on_x * b for a, b in zip(yc, xc)]
    zr = [a - z_on_x * b for a, b in zip(zc, xc)]
    return math.fsum(a * b for a, b in zip(zr, yr)) / math.fsum(value * value for value in zr)


def conditional_case(h, failures=True):
    condition_file = h.root / "condition.txt"
    condition_file.write_text("rs0001\n", encoding="ascii")
    prefix = h.run("conditional", ["--condition-snp", condition_file, "--pheno", h.data["pheno"]])
    rows = validate_output(Path(str(prefix) + ".conditional.assoc"), "conditional", h.data["n"])
    by_id = {row["SNP"]: row for row in rows}
    require("rs0002" not in by_id, "Perfect lead proxy was not removed by conditional collinearity")
    expected = conditional_beta(h.data["y"], impute(h.data["calls"][0]), impute(h.data["calls"][3]))
    observed = float(by_id["rs0004"]["BETA"])
    require(math.isclose(observed, expected, rel_tol=2e-6, abs_tol=2e-7),
            f"Conditional beta differs from OLS reference: {observed} vs {expected}")
    metadata = Path(str(prefix) + ".multilocus.meta.tsv").read_text(encoding="ascii")
    require("model\tconditional" in metadata and "condition_count\t1" in metadata,
            "Missing conditional metadata")
    if failures:
        missing = h.root / "missing_condition.txt"
        missing.write_text("not_present\n", encoding="ascii")
        h.run("missing_condition", ["--condition-snp", missing, "--pheno", h.data["pheno"]],
              error=r"condition SNP.*(absent|filtered|not found)")


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
    with tempfile.TemporaryDirectory(prefix=f"association_r3_{args.case}_", dir=out) as directory:
        root = Path(directory)
        data = fixture(root, 120, 300) if args.case == "benchmark" else fixture(root, 60, 48)
        harness = Harness(binary, root, data)
        if args.case == "farmcpu":
            multilocus_case(harness, "farmcpu")
        elif args.case == "blink":
            multilocus_case(harness, "blink")
        elif args.case == "conditional":
            conditional_case(harness)
        else:
            multilocus_case(harness, "farmcpu", failures=False)
            multilocus_case(harness, "blink", failures=False)
            conditional_case(harness, failures=False)
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
        print(f"FAIL association_r3: {error}", file=sys.stderr)
        sys.exit(1)
