#!/usr/bin/env python3
"""Cross-format BED/PGEN equivalence for dominance and additive-dominance GWAS.

Phase E2.1. The production dominance model is defined on hardcalls
(AA=0, AB=1, BB=0, centered heterozygote indicator). BED and PGEN both supply
exact hardcall allele counts, so this test asserts that the two backends are
statistically indistinguishable:

  * K_D (Vitezica, algorithm 4) and K_A (VanRaden, algorithm 1) match
    element-wise between BED and PGEN.
  * Dominance and AD null models (variance components) match under both
    global and LOCO scans.
  * BETA_A / SE_A / P_A / BETA_D / SE_D / P_D / P_JOINT match per marker.
  * Per-marker filtering decisions agree (monomorphic, constant-heterozygote,
    rare, missing-genotype markers).
  * Statistical conclusions agree (QTL detection and null calibration).
  * BGEN dominance is accepted through its probability path, while BGEN
    additive-dominance remains fail-fast.

Usage:
  python3 run_pgen_dominance_ad_validation.py --scanner <path> --workdir <dir>
    [--plink2 <plink2 path>] [--n 500] [--m 600] [--keep]
"""

import argparse
import csv
import math
import os
import random
import shutil
import struct
import subprocess
import sys
from pathlib import Path

DEFAULT_PLINK2 = os.environ.get(
    "PLINK2_BIN", "plink2")

UNIFIED_HEADER = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
                  "BETA_A", "SE_A", "P_A", "LOG10P_A",
                  "BETA_D", "SE_D", "P_D", "LOG10P_D",
                  "P_JOINT", "LOG10P_JOINT"]
DOM_COLS = ["BETA_D", "SE_D", "P_D"]
AD_COLS = ["BETA_A", "SE_A", "P_A", "BETA_D", "SE_D", "P_D", "P_JOINT"]
# BED and PGEN accumulate the GRM through different floating-point block
# layouts, so agreement is float32-storage level rather than bit exact.
TOL = {"BETA_A": 1e-5, "SE_A": 1e-5, "P_A": 1e-5,
       "BETA_D": 1e-5, "SE_D": 1e-5, "P_D": 1e-5, "P_JOINT": 1e-5}

DOMINANCE_INDEX = 7
ADDITIVE_INDEX = 11
MIXED_ADDITIVE_INDEX = 23
MIXED_DOMINANCE_INDEX = 19
MONOMORPHIC_INDEX = 0
ULTRA_RARE_INDEX = 1
RARE_INDEX = 2
MISSING_INDEX = 5
ALL_HET_INDEX = 6
HWD_INDEX = 4


def check(condition, message):
    if not condition:
        raise AssertionError("FAIL: " + message)
    print("  ok: " + message)


def run(command, env=None):
    return subprocess.run([str(c) for c in command], text=True,
                          capture_output=True, env=env)


def single_thread_env():
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return env


def write_dataset(directory, n, m, seed, chromosomes):
    """Deterministic hardcall fixture with dominance/AD edge cases.

    Marker layout (index: purpose):
      0        monomorphic                (must be filtered by MAF)
      1        rare (p=0.02)              (exercises rare-variant handling)
      5        scattered missingness      (mean-imputed, must be tested)
      6        all heterozygous           (constant indicator, rank deficient)
      7        dominance QTL locus
      11       additive QTL locus
      19       mixed-architecture dominance QTL locus
      23       mixed-architecture additive QTL locus
    """
    rng = random.Random(seed)
    genotypes = []
    for j in range(m):
        if j == MONOMORPHIC_INDEX:
            freq = 0.0
        elif j == ULTRA_RARE_INDEX:
            freq = 0.02
        elif j == RARE_INDEX:
            freq = 0.12
        elif j % 7 == 3:
            freq = 0.5                      # maximum heterozygosity
        else:
            freq = 0.06 + 0.88 * ((j * 37) % 100) / 100.0
        column = []
        for _ in range(n):
            u = rng.random()
            column.append(2 if u < freq * freq
                          else (1 if u < freq * freq + 2 * freq * (1 - freq) else 0))
        genotypes.append(column)

    for i in range(n):
        genotypes[MONOMORPHIC_INDEX][i] = 0
    for i in range(n):
        genotypes[ALL_HET_INDEX][i] = 1
    # Strong heterozygote excess (HWD) so the exact-HWE filter actually rejects
    # this marker, which makes the BED/PGEN GRM SNP-membership comparison bite.
    for i in range(n):
        genotypes[HWD_INDEX][i] = 1 if (i % 3) != 0 else 0
    for i in range(0, n, 5):                 # ~20% missing on one marker
        genotypes[MISSING_INDEX][i] = -1
    for j in range(0, m, 29):                # light background missingness
        for i in range(0, n, 17):
            genotypes[j][i] = -1

    directory.mkdir(parents=True, exist_ok=True)
    prefix = directory / "pdom"
    with open(str(prefix) + ".bed", "wb") as handle:
        handle.write(bytes([0x6C, 0x1B, 0x01]))
        for column in genotypes:
            buffer = bytearray((n + 3) // 4)
            for i, value in enumerate(column):
                code = 1 if value < 0 else (3 if value == 0 else (2 if value == 1 else 0))
                buffer[i // 4] |= code << (2 * (i % 4))
            handle.write(bytes(buffer))
    with open(str(prefix) + ".bim", "w") as handle:
        for j in range(m):
            handle.write(f"{chromosome_of(j, m, chromosomes)}\trs{j}\t0\t"
                         f"{1000 + j * 100}\tA\tG\n")
    with open(str(prefix) + ".fam", "w") as handle:
        for i in range(n):
            handle.write(f"F{i}\tI{i}\t0\t0\t0\t-9\n")

    calls = {f"rs{j}": [("NA" if v < 0 else str(v)) for v in column]
             for j, column in enumerate(genotypes)}
    return prefix, calls


def chromosome_of(j, m, chromosomes):
    return 1 + min(chromosomes - 1, j * chromosomes // m)


def center(values):
    mean = sum(values) / len(values)
    return [v - mean for v in values]


def write_phenotypes(directory, calls, n, seed, chromosomes):
    rng = random.Random(seed + 991)
    cov1 = [rng.gauss(0.0, 1.0) for _ in range(n)]
    cov2 = [rng.gauss(0.0, 1.0) for _ in range(n)]
    cov3 = [1.0 if i % 2 == 0 else 0.0 for i in range(n)]

    directory.mkdir(parents=True, exist_ok=True)
    covariate = directory / "covar.txt"
    with open(covariate, "w") as handle:
        for i in range(n):
            handle.write(f"F{i}\tI{i}\t{cov1[i]:.10f}\t{cov2[i]:.10f}\t{cov3[i]:.1f}\n")

    def dosage(snp):
        return [0.0 if g == "NA" else float(g) for g in calls[snp]]

    def hetero(snp):
        return [1.0 if g == "1" else 0.0 for g in calls[snp]]

    additive_qtl = center(dosage(f"rs{ADDITIVE_INDEX}"))
    dominance_qtl = center(hetero(f"rs{DOMINANCE_INDEX}"))
    mixed_additive = center(dosage(f"rs{MIXED_ADDITIVE_INDEX}"))
    mixed_dominance = center(hetero(f"rs{MIXED_DOMINANCE_INDEX}"))

    scenarios = {
        "null": lambda i, s: 0.5 * s,
        "additive": lambda i, s: 0.5 * s + 0.95 * additive_qtl[i],
        "dominance": lambda i, s: 0.5 * s + 2.0 * dominance_qtl[i],
        "mixed": lambda i, s: 0.5 * s + 0.85 * mixed_additive[i] + 1.7 * mixed_dominance[i],
    }
    written = {}
    for name, fn in scenarios.items():
        path = directory / f"{name}.pheno"
        with open(path, "w") as handle:
            for i in range(n):
                s = 0.5 * cov1[i] + 0.3 * cov2[i] + 0.2 * cov3[i]
                handle.write(f"F{i}\tI{i}\t{fn(i, s) + rng.gauss(0.0, 0.5):.10f}\n")
        written[name] = path
    return written, covariate


def read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def read_meta(path):
    values = {}
    if not Path(path).exists():
        return values
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            values[key] = value
    return values


def read_grm(path):
    raw = Path(path).read_bytes()
    return struct.unpack(f"<{len(raw) // 4}f", raw)


def write_dosage_vcf(path, n, m, phased=False, seed=97, het_boost=0.15):
    """Write a VCF whose FORMAT carries fractional dosages.

    With `phased=False` the file is meant to be imported with `dosage=DS`, which
    makes plink2 emit a PGEN whose records carry a dosage payload. With
    `phased=True` a phased hardcall VCF is written, which makes plink2 emit the
    hardcall-phase flag instead. Both are non-hardcall-only record layouts that
    the hardcall-only contract must reject.
    """
    rng = random.Random(seed)
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Expected alternate allele dosage">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"F{i}_I{i}" for i in range(n)),
    ]
    for j in range(m):
        cols = []
        for _ in range(n):
            u = rng.random()
            if phased:
                gt = "0|1" if u < 0.4 else ("1|1" if u < 0.7 else "0|0")
                cols.append(f"{gt}:.")
            else:
                if u < 0.3:
                    gt, ds = "0/0", rng.random() * 0.4
                elif u < 0.7:
                    gt, ds = "0/1", 1.0 + (1.0 - (0.2 + 0.6 * rng.random()))
                else:
                    gt, ds = "1/1", 1.6 + rng.random() * 0.4
                cols.append(f"{gt}:{ds:.4f}")
        lines.append(f"1\t{1000 + j * 100}\trs{j}\tA\tG\t.\tPASS\t.\tGT:DS\t"
                     + "\t".join(cols))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def scan_command(scanner, in_flag, prefix, phenotype, covariate, out, model, loco):
    command = [scanner, in_flag, str(prefix), "--pheno", str(phenotype),
               "--qcovar", str(covariate),
               "--gwas", "--gwas-model", model, "--maf", "0.01", "--out", str(out)]
    if not loco:
        command.append("--no-loco")
    return command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--plink2", default=DEFAULT_PLINK2)
    parser.add_argument("--n", type=int, default=500)
    parser.add_argument("--m", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20260927)
    parser.add_argument("--chromosomes", type=int, default=4)
    parser.add_argument("--n-small", type=int, default=200,
                        help="sample count for the dosage/phased PGEN contract fixtures")
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args()

    scanner = os.path.abspath(args.scanner)
    if not os.path.exists(scanner):
        sys.exit("scanner not found: " + scanner)
    args.plink2 = shutil.which(args.plink2) or args.plink2
    if not os.path.exists(args.plink2):
        sys.exit("plink2 not found: " + args.plink2)

    workdir = Path(args.workdir)
    if workdir.exists() and not args.keep:
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    env = single_thread_env()

    print("=== fixture ===")
    bed_prefix, calls = write_dataset(workdir, args.n, args.m, args.seed,
                                      args.chromosomes)
    pgen_prefix = workdir / "pdom_pgen"
    result = run([args.plink2, "--bfile", str(bed_prefix), "--make-pgen",
                  "--out", str(pgen_prefix)], env)
    check(result.returncode == 0, "plink2 converts the BED fixture to PGEN hardcalls")
    bgen_prefix = workdir / "pdom_bgen"
    result = run([args.plink2, "--bfile", str(bed_prefix), "--export", "bgen-1.2",
                  "--out", str(bgen_prefix)], env)
    check(result.returncode == 0, "plink2 exports the BED fixture to BGEN")
    phenotypes, covariate = write_phenotypes(workdir, calls, args.n, args.seed,
                                             args.chromosomes)

    print("=== K_D / K_A element-wise GRM reference comparison ===")
    for algorithm, label in (("4", "K_D (Vitezica, algorithm 4)"),
                             ("1", "K_A (VanRaden, algorithm 1)")):
        matrices = {}
        for name, in_flag, prefix in (("bed", "--bfile", bed_prefix),
                                      ("pgen", "--pfile", pgen_prefix)):
            out = workdir / f"grm{algorithm}_{name}"
            result = run([scanner, in_flag, str(prefix), "--make-grm",
                          "--make-grm-alg", algorithm, "--maf", "0.01",
                          "--out", str(out)], env)
            check(result.returncode == 0, f"{label} builds from {name}")
            matrices[name] = read_grm(str(out) + ".grm.bin")
        check(len(matrices["bed"]) == len(matrices["pgen"]),
              f"{label} has identical element count")
        check(len(matrices["bed"]) > 0, f"{label} is non-empty")
        worst = max(abs(x - y) for x, y in zip(matrices["bed"], matrices["pgen"]))
        scale = max(max(abs(v) for v in matrices["bed"]), 1e-12)
        print(f"  {label}: max abs diff={worst:.3e} max|K|={scale:.6f}")
        check(worst <= 1e-5 * max(scale, 1.0),
              f"{label} matches element-wise between BED and PGEN")

    print("=== scan equivalence ===")
    summary_rows = []
    for trait, phenotype in phenotypes.items():
        for model, cols in (("dominance", DOM_COLS), ("additive_dominance", AD_COLS)):
            suffix = ".assoc.dom.txt" if model == "dominance" else ".assoc.ad.txt"
            for loco in (True, False):
                outputs = {}
                for name, in_flag, prefix in (("bed", "--bfile", bed_prefix),
                                              ("pgen", "--pfile", pgen_prefix)):
                    out = workdir / f"{trait}_{model}_{'loco' if loco else 'noloco'}_{name}"
                    result = run(scan_command(scanner, in_flag, prefix, phenotype,
                                              covariate, out, model, loco), env)
                    check(result.returncode == 0,
                          f"{trait}/{model}/{'LOCO' if loco else 'no-LOCO'} {name} completes")
                    outputs[name] = str(out) + suffix
                label = f"{trait}/{model}/{'LOCO' if loco else 'no-LOCO'}"

                header_bed, rows_bed = read_table(outputs["bed"])
                header_pgen, rows_pgen = read_table(outputs["pgen"])
                check(header_bed == UNIFIED_HEADER, f"{label} BED uses the unified schema")
                check(header_pgen == header_bed, f"{label} PGEN header matches BED")

                bed_map = {r["SNP"]: r for r in rows_bed}
                pgen_map = {r["SNP"]: r for r in rows_pgen}
                check(set(bed_map) == set(pgen_map),
                      f"{label} tested marker sets are identical "
                      f"(bed={len(bed_map)} pgen={len(pgen_map)})")

                worst = {c: 0.0 for c in cols}
                identity_mismatches = []
                for snp in sorted(set(bed_map) & set(pgen_map)):
                    a, b = bed_map[snp], pgen_map[snp]
                    for column in ("CHR", "BP", "A1", "A2", "N", "FREQ"):
                        if a[column] != b[column]:
                            identity_mismatches.append(
                                f"{snp}.{column}={a[column]}!={b[column]}")
                    for column in cols:
                        worst[column] = max(worst[column],
                                            abs(float(a[column]) - float(b[column])))
                check(not identity_mismatches,
                      f"{label} identity columns (CHR/BP/A1/A2/N/FREQ) match for all "
                      f"{len(bed_map)} markers")
                for column in cols:
                    check(worst[column] <= TOL[column],
                          f"{label} max |diff {column}| = {worst[column]:.3e} <= {TOL[column]:.0e}")

                meta_bed = read_meta(str(outputs["bed"]).replace(suffix, ".scan.meta.txt"))
                meta_pgen = read_meta(str(outputs["pgen"]).replace(suffix, ".scan.meta.txt"))
                check(meta_bed.get("GENOTYPE_FORMAT") == "PLINK_BED_HARDCALL",
                      f"{label} BED records PLINK_BED_HARDCALL")
                check(meta_pgen.get("GENOTYPE_FORMAT") == "PLINK_PGEN_HARDCALL",
                      f"{label} PGEN records PLINK_PGEN_HARDCALL")
                for key in ("MODEL", "N", "N_SEMANTICS", "FREQ_SEMANTICS", "LOCO",
                            "COVARIATE_PROJECTION", "VCE_CONVERGED"):
                    check(meta_bed.get(key) == meta_pgen.get(key),
                          f"{label} metadata {key} matches ({meta_pgen.get(key)})")
                sigma_keys = ["SIGMA_D2", "SIGMA_E2"]
                if model == "additive_dominance":
                    sigma_keys.insert(0, "SIGMA_A2")
                for key in sigma_keys:
                    check(key in meta_pgen, f"{label} PGEN reports {key}")
                    check(meta_bed.get(key) == meta_pgen.get(key),
                          f"{label} variance component {key} matches")
                for key in sorted(k for k in meta_pgen if k.startswith("SIGMA")):
                    if key in sigma_keys or meta_pgen[key] == "PER_NULL_MODEL":
                        continue
                    if not (key.startswith("SIGMA_A2") or key.startswith("SIGMA_D2")
                            or key.startswith("SIGMA_E2")):
                        continue
                    a, b = meta_bed.get(key), meta_pgen.get(key)
                    if a is None or b is None:
                        check(a == b, f"{label} per-chromosome {key} present in both")
                        continue
                    # Absolute floor so degenerate boundary solutions near zero
                    # (e.g. sigma_D^2 collapsed to ~1e-11) are not compared
                    # by a meaningless relative tolerance.
                    tol = 1e-5 * max(abs(float(a)), abs(float(b)), 1e-6)
                    check(abs(float(a) - float(b)) <= tol,
                          f"{label} per-chromosome {key} matches ({a} vs {b})")

                summary_rows.append((trait, model, loco, len(bed_map)))

    print("=== HWE-filtered GRM/SNP-membership parity ===")
    # The PGEN GRM paths gained the exact-HWE filter so that GRM SNP membership
    # matches BED. Exercise the filter explicitly, on a phenotype and a marker
    # set where HWE actually rejects variants (HWD at a random-effect locus).
    for model, suffix, cols in (("dominance", ".assoc.dom.txt", DOM_COLS),
                                ("additive_dominance", ".assoc.ad.txt", AD_COLS)):
        outputs = {}
        for name, in_flag, prefix in (("bed", "--bfile", bed_prefix),
                                      ("pgen", "--pfile", pgen_prefix)):
            out = workdir / f"hwe_{model}_{name}"
            result = run([scanner, in_flag, str(prefix), "--pheno", str(phenotypes["null"]),
                          "--qcovar", str(covariate), "--gwas", "--gwas-model", model,
                          "--maf", "0.01", "--hwe", "1e-4", "--out", str(out)], env)
            check(result.returncode == 0, f"HWE-filtered {model} {name} completes")
            outputs[name] = str(out) + suffix
        _, rows_bed = read_table(outputs["bed"])
        _, rows_pgen = read_table(outputs["pgen"])
        bed_map = {r["SNP"]: r for r in rows_bed}
        pgen_map = {r["SNP"]: r for r in rows_pgen}
        check(set(bed_map) == set(pgen_map),
              f"HWE-filtered {model} tested marker sets are identical "
              f"(bed={len(bed_map)} pgen={len(pgen_map)})")
        worst = {c: 0.0 for c in cols}
        for snp in set(bed_map) & set(pgen_map):
            for column in cols:
                worst[column] = max(worst[column],
                                    abs(float(bed_map[snp][column])
                                        - float(pgen_map[snp][column])))
        for column in cols:
            check(worst[column] <= TOL[column],
                  f"HWE-filtered {model} max |diff {column}| = "
                  f"{worst[column]:.3e} <= {TOL[column]:.0e}")
        unfiltered = {r["SNP"] for r in
                      read_table(str(workdir / f"null_{model}_loco_bed") + suffix)[1]}
        check(set(bed_map).issubset(unfiltered),
              f"HWE-filtered {model} marker set is a subset of the unfiltered BED run")
        check(f"rs{HWD_INDEX}" not in bed_map and f"rs{HWD_INDEX}" not in pgen_map,
              f"HWE-filtered {model} excludes the heterozygote-excess marker "
              f"rs{HWD_INDEX} in both BED and PGEN")
        if f"rs{HWD_INDEX}" in unfiltered:
            # dominance-only: the marker is otherwise testable, so the HWE filter
            # is what removes it.
            check(len(bed_map) == len(unfiltered) - 1,
                  f"HWE-filtered {model} removes exactly the HWE-failing marker "
                  f"({len(bed_map)} == {len(unfiltered)} - 1)")
        else:
            # additive_dominance: the marker is already rank deficient in the
            # joint two-column design, so HWE is not the binding constraint.
            print(f"  note: rs{HWD_INDEX} was already excluded for {model} "
                  f"rank deficiency; HWE parity still verified")

    print("=== PGEN hardcall-only contract ===")
    # plink2 can only strip dosage from a PGEN, never add it, so a
    # dosage-carrying PGEN is produced by importing a VCF with fractional
    # dosages (`--vcf <f> dosage=DS`). Its variant-record type byte sets bits
    # above the 3-bit hardcall storage type; a plain phased hardcall PGEN sets
    # bit 4. Both must be rejected rather than silently reduced to hardcalls.
    for label, make_target in (("dosage", "dosage"), ("phased", "phased")):
        vcf = workdir / f"{label}.vcf"
        write_dosage_vcf(vcf, args.n_small, 6, phased=(make_target == "phased"))
        out = workdir / f"{label}_pgen"
        result = run([args.plink2, "--vcf", str(vcf), "--double-id", "--make-pgen",
                      "--out", str(out)], env)
        check(result.returncode == 0, f"plink2 builds the {label} PGEN fixture")
        if make_target == "dosage":
            result = run([args.plink2, "--vcf", str(vcf), "dosage=DS", "--double-id",
                          "--make-pgen", "--out", str(out)], env)
            check(result.returncode == 0,
                  "plink2 builds the dosage-carrying PGEN fixture")
        # Sample IDs follow the VCF import convention (FID == IID).
        small_pheno = workdir / f"{label}.pheno"
        with open(small_pheno, "w") as handle:
            for i in range(args.n_small):
                sid = f"F{i}_I{i}"
                handle.write(f"{sid}\t{sid}\t{0.1 * ((i % 7) - 3):.6f}\n")
        for model, suffix in (("dominance", ".assoc.dom.txt"),
                              ("additive_dominance", ".assoc.ad.txt")):
            prefix = workdir / f"{label}_{model}"
            result = run([scanner, "--pfile", str(out), "--pheno", str(small_pheno),
                          "--gwas", "--gwas-model", model, "--out", str(prefix)], env)
            combined = (result.stdout or "") + (result.stderr or "")
            check(result.returncode != 0,
                  f"{label} PGEN + {model} is rejected")
            check("dosage records are not supported" in combined,
                  f"{label} PGEN + {model} states the dosage hardcall-only contract")
            check("hardcall PGEN or BED" in combined,
                  f"{label} PGEN + {model} names the accepted alternatives")
            check(not Path(str(prefix) + suffix).exists(),
                  f"{label} PGEN + {model} writes no association output")

    # A genuine hardcall PGEN must still be accepted, and must advertise the
    # contract in its metadata.
    for model, suffix in (("dominance", ".assoc.dom.txt"),
                          ("additive_dominance", ".assoc.ad.txt")):
        prefix = workdir / f"hardcall_contract_{model}"
        result = run([scanner, "--pfile", str(pgen_prefix), "--pheno",
                      str(phenotypes["null"]), "--qcovar", str(covariate), "--gwas",
                      "--gwas-model", model, "--maf", "0.01", "--out", str(prefix)], env)
        check(result.returncode == 0, f"hardcall PGEN + {model} is accepted")
        check(Path(str(prefix) + suffix).exists(),
              f"hardcall PGEN + {model} writes association output")
        meta = read_meta(str(prefix) + ".scan.meta.txt")
        check(meta.get("GENOTYPE_FORMAT") == "PLINK_PGEN_HARDCALL",
              f"hardcall PGEN + {model} reports GENOTYPE_FORMAT=PLINK_PGEN_HARDCALL")
        check(meta.get("PGEN_MODE") == "HARDCALL",
              f"hardcall PGEN + {model} reports PGEN_MODE=HARDCALL")
        check(meta.get("PGEN_HARDCALL_ONLY") == "true",
              f"hardcall PGEN + {model} reports PGEN_HARDCALL_ONLY=true")

    print("=== BED metadata is unchanged by the contract ===")
    for model in ("dominance", "additive_dominance"):
        meta = read_meta(str(workdir / f"null_{model}_loco_bed") + ".scan.meta.txt")
        check(meta.get("GENOTYPE_FORMAT") == "PLINK_BED_HARDCALL",
              f"BED {model} reports GENOTYPE_FORMAT=PLINK_BED_HARDCALL")
        check(meta.get("PGEN_MODE") == "NOT_APPLICABLE",
              f"BED {model} reports PGEN_MODE=NOT_APPLICABLE")
        check(meta.get("PGEN_HARDCALL_ONLY") == "NOT_APPLICABLE",
              f"BED {model} reports PGEN_HARDCALL_ONLY=NOT_APPLICABLE")

    print("=== extract / exclude GRM and LOCO contract ===")
    # Markers spread across chromosomes so that LOCO still has SNPs outside
    # each target chromosome after filtering.
    span = max(1, args.m // 4)
    keep_ids = [f"rs{1}", f"rs{span + 5}", f"rs{2 * span + 5}"]
    drop_ids = [f"rs{3 * span + 7}", f"rs{span + 40}"]
    extract_file = workdir / "extract.txt"
    exclude_file = workdir / "exclude.txt"
    extract_file.write_text("\n".join(keep_ids) + "\n", encoding="utf-8")
    exclude_file.write_text("\n".join(drop_ids) + "\n", encoding="utf-8")
    expected_markers = set(keep_ids) - set(drop_ids)
    check(len(expected_markers) > 0, "extract/exclude fixture leaves markers to test")

    for model, suffix, cols in (("dominance", ".assoc.dom.txt", DOM_COLS),
                                ("additive_dominance", ".assoc.ad.txt", AD_COLS)):
        for loco in (True, False):
            tag = "LOCO" if loco else "no-LOCO"
            outputs = {}
            for name, in_flag, prefix in (("bed", "--bfile", bed_prefix),
                                          ("pgen", "--pfile", pgen_prefix)):
                out = workdir / f"extract_{model}_{'loco' if loco else 'noloco'}_{name}"
                result = run([scanner, in_flag, str(prefix), "--pheno",
                              str(phenotypes["mixed"]), "--qcovar", str(covariate),
                              "--gwas", "--gwas-model", model, "--maf", "0.01",
                              "--extract", str(extract_file),
                              "--exclude", str(exclude_file),
                              "--out", str(out)] + ([] if loco else ["--no-loco"]), env)
                check(result.returncode == 0,
                      f"extract/exclude {model}/{tag} {name} completes")
                outputs[name] = str(out) + suffix
            _, rows_bed = read_table(outputs["bed"])
            _, rows_pgen = read_table(outputs["pgen"])
            bed_map = {r["SNP"]: r for r in rows_bed}
            pgen_map = {r["SNP"]: r for r in rows_pgen}
            check(set(bed_map) == set(pgen_map),
                  f"extract/exclude {model}/{tag} marker sets are identical")
            check(set(bed_map).issubset(expected_markers),
                  f"extract/exclude {model}/{tag} only tests extract-minus-exclude "
                  f"markers ({sorted(bed_map)} subset of {sorted(expected_markers)})")
            worst = {c: 0.0 for c in cols}
            for snp in set(bed_map) & set(pgen_map):
                for column in ("CHR", "BP", "A1", "A2", "N", "FREQ"):
                    check(bed_map[snp][column] == pgen_map[snp][column],
                          f"extract/exclude {model}/{tag} {snp}.{column} matches")
                for column in cols:
                    worst[column] = max(worst[column],
                                        abs(float(bed_map[snp][column])
                                            - float(pgen_map[snp][column])))
            for column in cols:
                check(worst[column] <= TOL[column],
                      f"extract/exclude {model}/{tag} max |diff {column}| = "
                      f"{worst[column]:.3e} <= {TOL[column]:.0e}")
            # LOCO denominators must also agree, which they only do if the
            # per-chromosome split GRM saw the same SNP set.
            meta_bed = read_meta(str(workdir / f"extract_{model}_"
                                          f"{'loco' if loco else 'noloco'}_bed")
                                 + ".scan.meta.txt")
            meta_pgen = read_meta(str(workdir / f"extract_{model}_"
                                           f"{'loco' if loco else 'noloco'}_pgen")
                                  + ".scan.meta.txt")
            if loco:
                check(meta_bed.get("D_NULL_MODELS") == meta_pgen.get("D_NULL_MODELS"),
                      f"extract/exclude {model} LOCO null count matches")
                shared = sorted(k for k in meta_bed
                                if k.startswith("SIGMA_") and k in meta_pgen)
                check(len(shared) > 0,
                      f"extract/exclude {model} LOCO reports per-chromosome sigma")
                for key in shared:
                    a, b = meta_bed[key], meta_pgen[key]
                    if a == "PER_NULL_MODEL":
                        continue
                    tol = 1e-5 * max(abs(float(a)), abs(float(b)), 1e-6)
                    check(abs(float(a) - float(b)) <= tol,
                          f"extract/exclude {model} LOCO {key} matches ({a} vs {b})")

    print("=== LOCO membership audit ===")
    # With the same chromosome assignment, the global and chromosome-excluded
    # marker sets and both GRMs must agree between BED and PGEN. The
    # extract/exclude runs above already compare the split-GRM consequences;
    # here the unfiltered runs confirm the same for the default configuration.
    for model, suffix in (("dominance", ".assoc.dom.txt"),
                          ("additive_dominance", ".assoc.ad.txt")):
        for name in ("bed", "pgen"):
            meta = read_meta(str(workdir / f"null_{model}_loco_{name}") + ".scan.meta.txt")
            check(meta.get("LOCO") == "true",
                  f"{model}/{name} LOCO runs report LOCO=true")
            null_key = "D_NULL_MODELS" if model == "dominance" else "AD_NULL_MODELS"
            check(null_key in meta,
                  f"{model}/{name} reports {null_key} covering every chromosome")
            per_chrom = [k for k in meta
                         if k.startswith("SIGMA_") and k.count("_") == 2]
            check(len(per_chrom) > 0,
                  f"{model}/{name} reports per-chromosome LOCO variance components "
                  f"({len(per_chrom)} entries)")
        bed_meta = read_meta(str(workdir / f"null_{model}_loco_bed") + ".scan.meta.txt")
        pgen_meta = read_meta(str(workdir / f"null_{model}_loco_pgen") + ".scan.meta.txt")
        for key in sorted(k for k in bed_meta
                          if k.startswith("SIGMA_") and k in pgen_meta):
            a, b = bed_meta[key], pgen_meta[key]
            if a == "PER_NULL_MODEL":
                check(a == b, f"{model} {key} marker matches")
                continue
            tol = 1e-5 * max(abs(float(a)), abs(float(b)), 1e-6)
            check(abs(float(a) - float(b)) <= tol,
                  f"LOCO membership {model} {key} matches ({a} vs {b})")

    print("=== edge-marker agreement ===")
    # The AD model additionally requires a jointly non-degenerate additive +
    # dominance design, so an ultra-rare marker is legitimately rank deficient
    # there while remaining testable in dominance-only. BED and PGEN must reach
    # identical decisions in every case.
    for model, suffix, expectations in (
            ("dominance", ".assoc.dom.txt",
             {f"rs{MONOMORPHIC_INDEX}": False,
              f"rs{ALL_HET_INDEX}": False,
              f"rs{ULTRA_RARE_INDEX}": True,
              f"rs{RARE_INDEX}": True,
              f"rs{MISSING_INDEX}": True}),
            ("additive_dominance", ".assoc.ad.txt",
             {f"rs{MONOMORPHIC_INDEX}": False,
              f"rs{ALL_HET_INDEX}": False,
              f"rs{ULTRA_RARE_INDEX}": False,
              f"rs{RARE_INDEX}": True,
              f"rs{MISSING_INDEX}": True})):
        rows = {}
        for name in ("bed", "pgen"):
            path = str(workdir / f"null_{model}_loco_{name}") + suffix
            _, table = read_table(path)
            rows[name] = {r["SNP"]: r for r in table}
        check(set(rows["bed"]) == set(rows["pgen"]),
              f"{model} tested marker set is identical between BED and PGEN")
        for marker, expected_present in expectations.items():
            for name in ("bed", "pgen"):
                present = marker in rows[name]
                check(present == expected_present,
                      f"{model}/{name} {marker} is "
                      f"{'tested' if expected_present else 'excluded'}")

    print("=== statistical conclusions agree ===")
    for name in ("bed", "pgen"):
        _, table = read_table(str(workdir / f"dominance_dominance_loco_{name}")
                              + ".assoc.dom.txt")
        row = {r["SNP"]: r for r in table}[f"rs{DOMINANCE_INDEX}"]
        check(float(row["P_D"]) < 1e-8, f"{name} detects the pure dominance QTL")
        check(float(row["BETA_D"]) > 0.5, f"{name} dominance effect has the expected sign")

        _, table = read_table(str(workdir / f"mixed_additive_dominance_loco_{name}")
                              + ".assoc.ad.txt")
        mixed = {r["SNP"]: r for r in table}
        dom_row = mixed[f"rs{MIXED_DOMINANCE_INDEX}"]
        add_row = mixed[f"rs{MIXED_ADDITIVE_INDEX}"]
        check(float(dom_row["P_D"]) < 1e-6, f"{name} detects the mixed dominance QTL")
        check(float(add_row["P_A"]) < 1e-6, f"{name} detects the mixed additive QTL")
        check(float(dom_row["P_JOINT"]) < 1e-6, f"{name} mixed dominance joint test is significant")
        check(float(add_row["P_JOINT"]) < 1e-6, f"{name} mixed additive joint test is significant")

    for model, suffix in (("dominance", ".assoc.dom.txt"),
                          ("additive_dominance", ".assoc.ad.txt")):
        p_values = {}
        for name in ("bed", "pgen"):
            _, table = read_table(str(workdir / f"null_{model}_loco_{name}") + suffix)
            values = sorted(float(r["P_D"] if model == "dominance" else r["P_JOINT"])
                            for r in table)
            median = values[len(values) // 2]
            rate = sum(1 for v in values if v < 0.05) / len(values)
            p_values[name] = (median, rate)
            check(0.15 < median < 0.85,
                  f"{model}/{name} null median P is calibrated ({median:.4g})")
            check(rate < 0.15, f"{model}/{name} null P<0.05 rate is bounded ({rate:.4g})")
        check(abs(p_values["bed"][0] - p_values["pgen"][0]) < 0.02,
              f"{model} null median P agrees between BED and PGEN")

    print("=== BGEN dominance-only boundary ===")
    out = workdir / "bgen_dominance"
    result = run([scanner, "--bgen", str(bgen_prefix) + ".bgen",
                  "--pheno", str(phenotypes["null"]), "--qcovar", str(covariate),
                  "--gwas", "--gwas-model", "dominance", "--maf", "0.01",
                  "--out", str(out)], env)
    check(result.returncode == 0, "BGEN + dominance is accepted")
    check(Path(str(out) + ".assoc.dom.txt").exists(),
          "BGEN + dominance writes dominance output")
    bgen_meta = read_meta(str(out) + ".scan.meta.txt")
    check(bgen_meta.get("GENOTYPE_FORMAT") == "BGEN_PROBABILITY",
          "BGEN dominance records probability input")
    check(bgen_meta.get("KD_ENCODING") == "POSTERIOR_MEAN_VITEZICA_OUTER_PRODUCT",
          "BGEN dominance records posterior-mean Vitezica K_D")

    out = workdir / "bgen_guard_additive_dominance"
    result = run([scanner, "--bgen", str(bgen_prefix) + ".bgen",
                  "--pheno", str(phenotypes["null"]), "--qcovar", str(covariate),
                  "--gwas", "--gwas-model", "additive_dominance", "--maf", "0.01",
                  "--out", str(out)], env)
    combined = (result.stdout or "") + (result.stderr or "")
    check(result.returncode != 0, "BGEN + additive_dominance is rejected")
    check("not supported for BGEN" in combined,
          "BGEN additive_dominance states the explicit model guard")
    check(not Path(str(out) + ".assoc.ad.txt").exists(),
          "BGEN additive_dominance produces no silent fallback output")

    print()
    print("marker counts by scenario (identical for BED and PGEN):")
    for trait, model, loco, count in summary_rows:
        print(f"  {trait:9s} {model:18s} {'LOCO' if loco else 'no-LOCO':8s} {count}")
    print()
    print("PGEN DOMINANCE / AD BED-EQUIVALENCE VALIDATION PASSED")
    print("workdir:", workdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
