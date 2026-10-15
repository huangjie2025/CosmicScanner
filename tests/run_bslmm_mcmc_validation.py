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


def write_phenotype(data, path, causal, seed, noise):
    rng = random.Random(seed)
    means = {}
    for marker, _ in causal:
        observed = [x for x in data["calls"][marker] if x is not None]
        means[marker] = math.fsum(observed) / len(observed)
    values = []
    for sample in range(data["n"]):
        value = rng.gauss(0.0, noise)
        for marker, effect in causal:
            genotype = data["calls"][marker][sample]
            value += effect * (means[marker] if genotype is None else genotype)
        values.append(value)
    path.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {values[i]:.12g}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")


def run(binary, data, phenotype, prefix, mode, burnin, samples, thin, seed):
    command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(phenotype),
               mode, "--mcmc-burnin", str(burnin), "--mcmc-samples", str(samples),
               "--mcmc-thin", str(thin), "--seed", str(seed), "--maf", "0",
               "--threads", "1", "--out", str(prefix)]
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    require(result.returncode == 0, "Command failed:\n" + " ".join(command) + "\n" + result.stdout)


def read_summary(path):
    rows = list(csv.reader(path.open(encoding="ascii"), delimiter="\t"))
    require(rows and rows[0] == ["METRIC", "VALUE"], f"Bad posterior header: {path}")
    values = dict(rows[1:])
    for key in ("PVE", "PGE", "PI", "RESIDUAL_VARIANCE", "POLYGENIC_VARIANCE",
                "INCLUSION_CHANGE_RATE", "PVE_SPLIT_HALF_DIFFERENCE"):
        require(math.isfinite(float(values[key])), f"Non-finite {key}: {path}")
    require(0.0 <= float(values["PVE"]) <= 1.0, "PVE outside [0,1]")
    require(0.0 <= float(values["PGE"]) <= 1.0, "PGE outside [0,1]")
    require(0.0 < float(values["PI"]) <= 1.0, "PI outside (0,1]")
    require(float(values["RESIDUAL_VARIANCE"]) > 0.0, "Residual variance is not positive")
    require(0.0 <= float(values["INCLUSION_CHANGE_RATE"]) <= 1.0,
            "Invalid inclusion change rate")
    return values


def read_effects(path):
    rows = list(csv.DictReader(path.open(encoding="ascii"), delimiter="\t"))
    require(rows and set(rows[0]) == {"SNP", "A1", "A2", "POSTERIOR_BETA", "PIP"},
            f"Bad effects schema: {path}")
    for row in rows:
        beta = float(row["POSTERIOR_BETA"])
        pip = float(row["PIP"])
        require(math.isfinite(beta) and math.isfinite(pip) and 0.0 <= pip <= 1.0,
                f"Invalid effect row: {row}")
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    binary = args.bin.resolve()
    args.out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="bslmm_mcmc_", dir=args.out) as directory:
        root = Path(directory)
        data = fixture(root, n=100, m=60)
        sparse = root / "sparse.pheno"
        null = root / "null.pheno"
        causal = [(7, 2.2), (19, -1.8)]
        write_phenotype(data, sparse, causal, 20260911, 0.55)
        write_phenotype(data, null, [], 20260912, 1.0)

        bslmm = root / "bslmm"
        bslmm_repeat = root / "bslmm_repeat"
        bvsr_short = root / "bvsr_short"
        bvsr_long = root / "bvsr_long"
        null_bvsr = root / "null_bvsr"
        run(binary, data, sparse, bslmm, "--bslmm-mcmc", 150, 450, 3, 91)
        run(binary, data, sparse, bslmm_repeat, "--bslmm-mcmc", 150, 450, 3, 91)
        run(binary, data, sparse, bvsr_short, "--bvsr", 150, 450, 3, 117)
        run(binary, data, sparse, bvsr_long, "--bvsr", 250, 900, 3, 117)
        run(binary, data, null, null_bvsr, "--bvsr", 150, 450, 3, 123)

        bslmm_summary = read_summary(Path(str(bslmm) + ".bslmm.posterior.tsv"))
        bvsr_summary = read_summary(Path(str(bvsr_short) + ".bslmm.posterior.tsv"))
        require(float(bslmm_summary["POLYGENIC_VARIANCE"]) > 0.0,
                "BSLMM did not sample a polygenic component")
        require(float(bvsr_summary["POLYGENIC_VARIANCE"]) == 0.0,
                "BVSR unexpectedly sampled a polygenic component")

        bslmm_effect_path = Path(str(bslmm) + ".bslmm.effects.tsv")
        repeat_effect_path = Path(str(bslmm_repeat) + ".bslmm.effects.tsv")
        require(bslmm_effect_path.read_bytes() == repeat_effect_path.read_bytes(),
                "Fixed-seed BSLMM effects are not reproducible")
        require(Path(str(bslmm) + ".bslmm.posterior.tsv").read_bytes() ==
                Path(str(bslmm_repeat) + ".bslmm.posterior.tsv").read_bytes(),
                "Fixed-seed BSLMM summaries are not reproducible")

        expected = {data["variants"][index][1] for index, _ in causal}
        short_rows = read_effects(Path(str(bvsr_short) + ".bvsr.effects.tsv"))
        long_rows = read_effects(Path(str(bvsr_long) + ".bvsr.effects.tsv"))
        top_short = {row["SNP"] for row in sorted(short_rows, key=lambda x: float(x["PIP"]), reverse=True)[:8]}
        top_long = {row["SNP"] for row in sorted(long_rows, key=lambda x: float(x["PIP"]), reverse=True)[:8]}
        require(expected.issubset(top_short) and expected.issubset(top_long),
                f"Causal SNPs not ranked near top: short={top_short}, long={top_long}")
        short_by_id = {row["SNP"]: float(row["POSTERIOR_BETA"]) for row in short_rows}
        long_by_id = {row["SNP"]: float(row["POSTERIOR_BETA"]) for row in long_rows}
        require(all(short_by_id[snp] * long_by_id[snp] > 0.0 for snp in expected),
                "Longer-chain causal effect directions are unstable")

        null_rows = read_effects(Path(str(null_bvsr) + ".bvsr.effects.tsv"))
        require(sum(float(row["PIP"]) >= 0.5 for row in null_rows) <= 3,
                "Null phenotype produced widespread high PIP")

        old = root / "old"
        command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(sparse),
                   "--bslmm", "2", "--bslmm-h2", "0.5", "--maf", "0", "--out", str(old)]
        result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, check=False)
        require(result.returncode == 0, "Old BSLMM path failed:\n" + result.stdout)
        metadata = Path(str(old) + ".bslmm.meta.txt").read_text(encoding="ascii")
        require("implementation\tdeterministic-ridge-gblup" in metadata,
                "Old deterministic BSLMM implementation contract changed")
        for suffix in (".param.txt", ".bv.txt", ".hyp.txt"):
            require(Path(str(old) + suffix).stat().st_size > 0, f"Missing old BSLMM output {suffix}")

    print("PASS BSLMM/BVSR MCMC validation")


if __name__ == "__main__":
    main()
