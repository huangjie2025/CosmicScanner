#!/usr/bin/env python3
"""Production validation for the dominance-only (d+D) GWAS path."""

import argparse
import csv
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys


def check(condition, message):
    if not condition:
        raise AssertionError(message)
    print("  PASS " + message)


def write_dataset(workdir, n, m, seed):
    rng = random.Random(seed)
    prefix = workdir / "dominance"
    genotypes = []
    for marker in range(m):
        p = 0.08 + 0.40 * ((marker * 37) % 101) / 100.0
        row = []
        for _ in range(n):
            u = rng.random()
            row.append(2 if u < p * p else (1 if u < p * p + 2 * p * (1 - p) else 0))
        genotypes.append(row)

    additive_index = 17
    dominance_index = 41
    monomorphic_index = 0
    all_heterozygous_index = 1
    genotypes[monomorphic_index] = [0] * n
    genotypes[all_heterozygous_index] = [1] * n
    genotypes[additive_index] = [i % 3 for i in range(n)]
    genotypes[dominance_index] = [i % 3 for i in range(n)]
    for i in range(0, n, 17):
        genotypes[2][i] = None

    with open(str(prefix) + ".bed", "wb") as output:
        output.write(b"\x6c\x1b\x01")
        for row in genotypes:
            block = bytearray((n + 3) // 4)
            for i, genotype in enumerate(row):
                code = 1 if genotype is None else {0: 3, 1: 2, 2: 0}[genotype]
                block[i // 4] |= code << (2 * (i % 4))
            output.write(block)
    with open(str(prefix) + ".bim", "w", encoding="utf-8") as output:
        for marker in range(m):
            output.write(f"{marker % 3 + 1}\tsnp{marker}\t0\t{1000 + marker * 10}\tA\tG\n")
    with open(str(prefix) + ".fam", "w", encoding="utf-8") as output:
        for i in range(n):
            output.write(f"F{i}\tI{i}\t0\t0\t0\t-9\n")

    covariate = [rng.gauss(0.0, 1.0) for _ in range(n)]
    with open(workdir / "covar.txt", "w", encoding="utf-8") as output:
        output.write("FID\tIID\tC1\n")
        for i, value in enumerate(covariate):
            output.write(f"F{i}\tI{i}\t{value:.12g}\n")

    additive = [genotypes[additive_index][i] - 1.0 for i in range(n)]
    dominance_raw = [1.0 if genotypes[dominance_index][i] == 1 else 0.0 for i in range(n)]
    dominance_mean = sum(dominance_raw) / n
    dominance = [value - dominance_mean for value in dominance_raw]
    phenotypes = {
        "null": [0.6 * covariate[i] + rng.gauss(0.0, 1.0) for i in range(n)],
        "pure": [0.6 * covariate[i] + 1.8 * dominance[i] + rng.gauss(0.0, 0.65)
                 for i in range(n)],
        "mixed": [0.6 * covariate[i] + 0.7 * additive[i] + 1.8 * dominance[i]
                  + rng.gauss(0.0, 0.65) for i in range(n)],
    }
    for label, values in phenotypes.items():
        with open(workdir / f"{label}.pheno", "w", encoding="utf-8") as output:
            output.write("FID\tIID\tTrait\n")
            for i, value in enumerate(values):
                output.write(f"F{i}\tI{i}\t{value:.12g}\n")
    return prefix, dominance_index, monomorphic_index, all_heterozygous_index


def run(scanner, prefix, phenotype, covariate, output, model, loco, extra=()):
    command = [scanner, "--bfile", str(prefix), "--pheno", str(phenotype),
               "--qcovar", str(covariate), "--out", str(output), "--gwas",
               "--gwas-model", model, "--maf", "0.01"]
    if not loco:
        command.append("--no-loco")
    command.extend(extra)
    environment = dict(os.environ)
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return subprocess.run(command, text=True, capture_output=True, env=environment)


def read_table(path):
    with open(path, encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        return reader.fieldnames, list(reader)


def by_snp(rows):
    return {row["SNP"]: row for row in rows}


def require_success(result, label):
    if result.returncode != 0:
        raise AssertionError(f"{label} failed\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def expect_failure(scanner, prefix, phenotype, covariate, workdir, extra, needle):
    result = run(scanner, prefix, phenotype, covariate, workdir / ("guard_" + str(len(extra))),
                 "dominance", False, extra)
    combined = (result.stdout or "") + (result.stderr or "")
    check(result.returncode != 0 and needle in combined,
          f"CLI guard rejects {' '.join(extra)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--n", type=int, default=500)
    parser.add_argument("--m", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--calibration-only", action="store_true")
    args = parser.parse_args()

    scanner = os.path.abspath(args.scanner)
    workdir = Path(args.workdir)
    if workdir.exists() and not args.keep:
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    prefix, causal, monomorphic, all_heterozygous = write_dataset(
        workdir, args.n, args.m, args.seed)
    covariate = workdir / "covar.txt"

    null_output = workdir / "null"
    result = run(scanner, prefix, workdir / "null.pheno", covariate,
                 null_output, "dominance", False)
    require_success(result, "dominance null calibration")
    header, null_rows = read_table(str(null_output) + ".assoc.dom.txt")
    expected_header = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
                       "BETA_A", "SE_A", "P_A", "LOG10P_A",
                       "BETA_D", "SE_D", "P_D", "LOG10P_D",
                       "P_JOINT", "LOG10P_JOINT"]
    check(header == expected_header, "dominance output header is complete")
    check(len(null_rows) >= args.m - 3, "eligible null markers are emitted")
    null_map = by_snp(null_rows)
    check(f"snp{monomorphic}" not in null_map, "monomorphic marker is filtered")
    check(f"snp{all_heterozygous}" not in null_map,
          "constant heterozygote indicator is rank deficient")
    numeric_names = ("FREQ", "BETA_D", "SE_D", "P_D", "LOG10P_D")
    check(all(all(math.isfinite(float(row[name])) for name in numeric_names) for row in null_rows),
          "all dominance output statistics are finite")
    check(all(row[name] == "NA" for row in null_rows
              for name in ("BETA_A", "SE_A", "P_A", "LOG10P_A",
                           "P_JOINT", "LOG10P_JOINT")),
          "unavailable additive and joint fields are NA")
    check(all(0.0 <= float(row["P_D"]) <= 1.0 for row in null_rows),
          "all P_D values are in [0,1]")
    p_values = sorted(float(row["P_D"]) for row in null_rows)
    median_p = p_values[len(p_values) // 2]
    fraction_005 = sum(value < 0.05 for value in p_values) / len(p_values)
    check(0.15 < median_p < 0.85, f"null median P is calibrated ({median_p:.4g})")
    check(fraction_005 < 0.15, f"null P<0.05 rate is bounded ({fraction_005:.4g})")

    if args.calibration_only:
        summary = workdir / "dominance_validation.summary.tsv"
        summary.write_text(
            "metric\tvalue\n"
            f"n\t{args.n}\n"
            f"m\t{args.m}\n"
            f"null_median_p\t{median_p:.12g}\n"
            f"null_fraction_p_lt_0.05\t{fraction_005:.12g}\n",
            encoding="utf-8")
        print("Dominance-only calibration validation passed")
        return 0

    pure_output = workdir / "pure_loco"
    result = run(scanner, prefix, workdir / "pure.pheno", covariate,
                 pure_output, "dominance", True)
    require_success(result, "dominance LOCO QTL scan")
    _, pure_rows = read_table(str(pure_output) + ".assoc.dom.txt")
    pure_map = by_snp(pure_rows)
    causal_row = pure_map[f"snp{causal}"]
    check(float(causal_row["P_D"]) < 1e-8, "pure dominance QTL is detected under LOCO")
    check(float(causal_row["BETA_D"]) > 1.0, "dominance QTL effect has the expected sign")
    meta = Path(str(pure_output) + ".scan.meta.txt").read_text(encoding="utf-8")
    for token in ("GWAS_MODEL\tdominance", "MODEL\tdominance", "N\t500",
                  "N_SEMANTICS\t", "FREQ_SEMANTICS\t", "KD_NORMALIZATION\tSUM_4P2Q2",
                  "LOCO\ttrue", "COVARIATE_PROJECTION\ttrue", "VCE_CONVERGED\t",
                  "SIGMA_D2\tPER_NULL_MODEL", "SIGMA_E2\tPER_NULL_MODEL",
                  "COVARIATE_COLUMNS\t2", "D_NULL_MODELS\t3"):
        check(token in meta, f"metadata contains {token}")

    mixed_output = workdir / "mixed_d"
    require_success(run(scanner, prefix, workdir / "mixed.pheno", covariate,
                        mixed_output, "dominance", False), "mixed-architecture dominance scan")
    _, mixed_rows = read_table(str(mixed_output) + ".assoc.dom.txt")
    mixed_causal = by_snp(mixed_rows)[f"snp{causal}"]
    check(float(mixed_causal["P_D"]) < 1e-8,
          "dominance QTL is detected with an additional additive QTL")
    check(float(mixed_causal["BETA_D"]) > 1.0,
          "mixed-architecture dominance effect has the expected sign")

    d_output = workdir / "pure_d"
    ad_output = workdir / "pure_ad"
    require_success(run(scanner, prefix, workdir / "pure.pheno", covariate,
                        d_output, "dominance", False), "D-only comparison scan")
    require_success(run(scanner, prefix, workdir / "pure.pheno", covariate,
                        ad_output, "additive_dominance", False), "AD comparison scan")
    _, d_rows = read_table(str(d_output) + ".assoc.dom.txt")
    _, ad_rows = read_table(str(ad_output) + ".assoc.ad.txt")
    d_causal = by_snp(d_rows)[f"snp{causal}"]
    ad_causal = by_snp(ad_rows)[f"snp{causal}"]
    d_logp = -math.log10(max(float(d_causal["P_D"]), 1e-300))
    ad_logp = -math.log10(max(float(ad_causal["P_D"]), 1e-300))
    check(abs(d_logp - ad_logp) < 2.0,
          f"D-only and AD conditional dominance evidence agree ({d_logp:.3g} vs {ad_logp:.3g})")
    check(abs(float(d_causal["BETA_D"]) - float(ad_causal["BETA_D"])) < 0.35,
          "D-only and AD dominance effect estimates agree")

    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--no-covar-proj"], "requires covariate projection")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--out-gz"], "does not yet support --out-gz")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--out-ssf"], "does not yet support --out-ssf")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--manhattan-tsv"], "does not yet support --manhattan-tsv/--qq-tsv")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--fastlmm"], "does not support --fastlmm/--lmm-low-rank")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--lmm-low-rank", "4"], "does not support --fastlmm/--lmm-low-rank")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--saige-style", "--biobank-mode"], "does not support --saige-style")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--lmm-gemma-mode", "4"], "does not support --lmm-gemma-mode/--lmm-gemma-refit")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--lmm-gemma-refit"], "does not support --lmm-gemma-mode/--lmm-gemma-refit")
    expect_failure(scanner, prefix, workdir / "pure.pheno", covariate, workdir,
                   ["--snp-model", "heterozygote"], "requires additive --snp-model")

    summary = workdir / "dominance_validation.summary.tsv"
    summary.write_text(
        "metric\tvalue\n"
        f"n\t{args.n}\n"
        f"m\t{args.m}\n"
        f"null_median_p\t{median_p:.12g}\n"
        f"null_fraction_p_lt_0.05\t{fraction_005:.12g}\n"
        f"loco_causal_p\t{causal_row['P_D']}\n"
        f"d_only_causal_log10p\t{d_logp:.12g}\n"
        f"ad_causal_log10p\t{ad_logp:.12g}\n",
        encoding="utf-8")
    print("Dominance-only production validation passed")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError) as error:
        print("FAIL:", error, file=sys.stderr)
        sys.exit(1)
