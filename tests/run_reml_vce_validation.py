#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
import random
import subprocess
import tempfile

from run_association_r2 import fixture


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def standardized(values):
    mean = math.fsum(values) / len(values)
    centered = [x - mean for x in values]
    variance = math.fsum(x * x for x in centered) / (len(values) - 1)
    return [x / math.sqrt(max(1e-12, variance)) for x in centered]


def write_simulations(data, root):
    rng = random.Random(20260911)
    genetic = [0.0] * data["n"]
    for marker in range(5, 35):
        observed = [x for x in data["calls"][marker] if x is not None]
        mean = math.fsum(observed) / len(observed)
        effect = rng.gauss(0.0, 1.0)
        for sample in range(data["n"]):
            value = data["calls"][marker][sample]
            genetic[sample] += effect * ((mean if value is None else value) - mean)
    genetic = standardized(genetic)
    noise = standardized([rng.gauss(0.0, 1.0) for _ in range(data["n"])])
    paths = {}
    for name, h2 in (("low", 0.1), ("moderate", 0.5), ("high", 0.8)):
        phenotype = [math.sqrt(h2) * g + math.sqrt(1.0 - h2) * e
                     for g, e in zip(genetic, noise)]
        path = root / f"{name}.pheno"
        path.write_text("".join(
            f"{data['samples'][i][0]} {data['samples'][i][1]} {phenotype[i]:.12g}\n"
            for i in reversed(range(data["n"]))), encoding="ascii")
        paths[name] = path
    boundary = root / "boundary.pheno"
    boundary_values = standardized([rng.gauss(0.0, 1.0) for _ in range(data["n"])])
    boundary.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {boundary_values[i]:.12g}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")
    paths["boundary"] = boundary
    return paths


def run(binary, data, phenotype, prefix, method, production=True):
    route = ["--reml", "--reml-method", method] if production else ["--vc", "2"]
    command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(phenotype),
               *route, "--reml-maxit", "800", "--reml-tol", "1e-7",
               "--maf", "0", "--threads", "1", "--out", str(prefix)]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    require(result.returncode == 0, "Command failed:\n" + " ".join(command) + "\n" + result.stdout)


def read_reml(path):
    rows = list(csv.DictReader(path.open(encoding="ascii"), delimiter="\t"))
    required = {"COMPONENT", "VARIANCE", "SE", "PROPORTION", "HERITABILITY",
                "ITERATIONS", "CONVERGED", "MAX_SCORE", "LOG_LIKELIHOOD", "OPTIMIZER"}
    require(len(rows) == 2 and set(rows[0]) == required, f"Bad REML output: {path}")
    for row in rows:
        for key in ("VARIANCE", "SE", "PROPORTION", "HERITABILITY",
                    "MAX_SCORE", "LOG_LIKELIHOOD"):
            require(math.isfinite(float(row[key])), f"Non-finite {key}: {row}")
        require(float(row["VARIANCE"]) >= 0.0 and float(row["SE"]) >= 0.0,
                f"Negative variance/SE: {row}")
        require(row["CONVERGED"] == "1" and int(row["ITERATIONS"]) >= 1,
                f"REML did not converge: {row}")
    return float(rows[0]["HERITABILITY"]), rows


def read_old_summary(path):
    values = {}
    for line in path.read_text(encoding="ascii").splitlines():
        key, value = line.split("\t", 1)
        values[key] = value
    return float(values["sigma2_g"]), float(values["sigma2_e"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    binary = args.bin.resolve()
    args.out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="reml_vce_", dir=args.out) as directory:
        root = Path(directory)
        data = fixture(root, n=140, m=90)
        phenotypes = write_simulations(data, root)
        estimates = {"ai": {}, "em": {}}
        for method in estimates:
            for name in ("low", "moderate", "high", "boundary"):
                prefix = root / f"{method}_{name}"
                run(binary, data, phenotypes[name], prefix, method)
                h2, rows = read_reml(Path(str(prefix) + ".reml.tsv"))
                estimates[method][name] = h2
                require(all(row["OPTIMIZER"] == method for row in rows),
                        f"Wrong optimizer label for {method}")

        for method, values in estimates.items():
            require(values["high"] > values["low"],
                    f"{method} did not distinguish high and low h2: {values}")
            require(0.0 <= values["boundary"] <= 1.0,
                    f"{method} boundary estimate is uncontrolled: {values['boundary']}")
        for name in ("low", "moderate", "high"):
            require(abs(estimates["ai"][name] - estimates["em"][name]) < 2e-3,
                    f"AI/EM mismatch for {name}: {estimates}")

        old = root / "old"
        production = root / "production_brent"
        run(binary, data, phenotypes["moderate"], old, "brent", production=False)
        run(binary, data, phenotypes["moderate"], production, "brent")
        old_vg, old_ve = read_old_summary(Path(str(old) + ".vc.summary.txt"))
        _, rows = read_reml(Path(str(production) + ".reml.tsv"))
        new_vg = float(rows[0]["VARIANCE"])
        new_ve = float(rows[1]["VARIANCE"])
        require(abs(old_vg - new_vg) < 1e-9 and abs(old_ve - new_ve) < 1e-9,
                "Standalone REML does not match the internal Brent VCE path")

    print("PASS unified REML/VCE validation")


if __name__ == "__main__":
    main()
