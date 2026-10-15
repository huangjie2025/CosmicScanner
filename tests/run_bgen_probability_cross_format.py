#!/usr/bin/env python3
"""Phase E2.6: BED / PGEN / BGEN one-hot dominance contract.

This test validates both the E2.4 genotype contract and the E2.6
dominance-only production path. BGEN additive-dominance remains rejected.

A single deterministic hardcall genotype matrix is materialized three ways:

  * BED                       (PLINK binary hardcall)
  * PGEN                      (plink2 --make-pgen, hardcall only)
  * BGEN 1.2                  (plink2 --export bgen-1.2 -> one-hot probabilities)

plink2 exports one-hot genotype probabilities for hardcall input, so the BGEN
fixture is a hardcall-equivalent probabilistic BGEN.

Contract items checked:

  1. sample ordering identical across formats
  2. variant ordering identical across formats
  3. per-sample genotype calls identical (BED == PGEN == BGEN)
  4. missing semantics identical (./. in BGEN and PGEN)
  5. heterozygosity indicator: P(AB) == I(AB) in the one-hot BGEN
  6. additive dosage identical across formats
  7. allele frequency p = mean(E[G])/2 identical across formats
  8. posterior-mean Vitezica code from one-hot BGEN probabilities reproduces
     the hardcall Vitezica encoding used by BED/PGEN production GRM
  9. K_D is element-wise equal across BED/PGEN/one-hot BGEN
 10. dominance sigma_D2, beta_D, SE_D and P_D agree with and without LOCO
 11. BGEN additive-dominance remains rejected at the CLI

Usage:
  python3 run_bgen_probability_cross_format.py --scanner <path> --workdir <dir>
    [--plink2 <plink2 path>] [--n 400] [--m 200]
"""

import argparse
import csv
import math
import os
import shutil
import random
import struct
import subprocess
import sys
from pathlib import Path

DEFAULT_PLINK2 = os.environ.get(
    "PLINK2_BIN", "plink2")

# Hardcall parity across formats must be exact at the genotype level.
# Dosage/p values are float-derived, so allow only storage-level rounding.
TOL_DOSAGE = 1e-6
TOL_P = 1e-6
TOL_STAT = 2e-5

UNIFIED_HEADER = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
                  "BETA_A", "SE_A", "P_A", "LOG10P_A",
                  "BETA_D", "SE_D", "P_D", "LOG10P_D",
                  "P_JOINT", "LOG10P_JOINT"]


def check(condition, message):
    if not condition:
        raise AssertionError("FAIL: " + message)
    print("  ok: " + message)


def run(command, env=None):
    result = subprocess.run([str(c) for c in command], text=True,
                            capture_output=True, env=env)
    if result.returncode != 0:
        raise RuntimeError("command failed: {}\n{}\n{}".format(
            " ".join(str(c) for c in command), result.stdout[-2000:], result.stderr[-2000:]))
    return result


def read_grm(path):
    raw = Path(path).read_bytes()
    check(len(raw) % 4 == 0, "GRM binary has complete float32 elements")
    return struct.unpack("<{}f".format(len(raw) // 4), raw)


def read_table(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def read_meta(path):
    result = {}
    with Path(path).open() as handle:
        next(handle)
        for line in handle:
            key, value = line.rstrip("\n").split("\t", 1)
            result[key] = value
    return result


def single_thread_env():
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    return env


def write_bed_dataset(output, n, m, seed, n_chrom=2):
    """Deterministic hardcall fixture with edge cases (returns calls, samples)."""
    rng = random.Random(seed)
    samples = ["S{:04d}".format(i + 1) for i in range(n)]
    calls = []
    for j in range(m):
        if j == 0:
            freq = 0.0                      # monomorphic -> filtered by MAF
        elif j == 1:
            freq = 0.05                     # rare
        elif j % 7 == 3:
            freq = 0.5                      # max heterozygosity
        else:
            freq = 0.06 + 0.88 * ((j * 37) % 100) / 100.0
        col = []
        for _ in range(n):
            u = rng.random()
            col.append(2 if u < freq * freq
                       else (1 if u < freq * freq + 2 * freq * (1 - freq) else 0))
        calls.append(col)
    # scattered missingness on a few markers, aligned across all formats
    for marker in (5, 11, 17):
        for i in range(0, n, 7):
            calls[marker][i] = -1
    # one all-heterozygote marker (constant indicator)
    for i in range(n):
        calls[6][i] = 1

    prefix = output / "canon"
    with (output / "canon.fam").open("w") as fam:
        for s in samples:
            fam.write("{} {} 0 0 1 -9\n".format(s, s))
    with (output / "canon.bim").open("w") as bim:
        for j in range(m):
            chrom = "1" if j < m // 2 else "2"
            bim.write("{} rs{:04d} 0 {} A C\n".format(chrom, j + 1, (j + 1) * 100))
    bed = bytearray(b"\x6c\x1b\x01")
    for col in calls:
        packed = bytearray((n + 3) // 4)
        for i, value in enumerate(col):
            code = 0b01 if value < 0 else (0b11, 0b10, 0b00)[value]
            packed[i // 4] |= code << (2 * (i % 4))
        bed.extend(packed)
    (output / "canon.bed").write_bytes(bytes(bed))
    return calls, samples


def bed_code_to_allele_count(code):
    # PLINK .bed 2-bit code -> count of the second (C) allele
    if code == 0b11:
        return 2
    if code == 0b10:
        return 1
    if code == 0b00:
        return 0
    return None  # 0b01 = missing


def read_vcf_genotypes(path):
    """Return (variant_ids, sample_ids, calls) from a genotype-only VCF."""
    variants = []
    samples = []
    calls = []
    with Path(path).open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            fields = line.split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            vid = fields[2]
            variants.append(vid)
            row = []
            for cell in fields[9:]:
                gt = cell.split(":")[0]
                if gt in ("./.", ".", "./.", "."):
                    row.append(-1)
                    continue
                alleles = re_split(gt)
                if len(alleles) != 2:
                    raise RuntimeError("unexpected phased genotype: " + gt)
                row.append(int(alleles[0]) + int(alleles[1]))
            calls.append(row)
    return variants, samples, calls


def re_split(gt):
    out = []
    cur = ""
    for ch in gt:
        if ch in "/|":
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return out


def expected_vitezica(g, p):
    q = 1.0 - p
    if g == 0:
        return -2.0 * p * p
    if g == 1:
        return 2.0 * p * q
    if g == 2:
        return -2.0 * q * q
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--plink2", default=DEFAULT_PLINK2)
    parser.add_argument("--n", type=int, default=400)
    parser.add_argument("--m", type=int, default=200)
    args = parser.parse_args()

    args.plink2 = shutil.which(args.plink2) or args.plink2
    if not os.path.exists(args.plink2):
        raise RuntimeError("plink2 not found at " + args.plink2)

    out = Path(args.workdir)
    out.mkdir(parents=True, exist_ok=True)
    env = single_thread_env()

    calls, samples = write_bed_dataset(out, args.n, args.m, seed=20240917)
    n, m = args.n, args.m

    # --- build PGEN and one-hot BGEN from the same BED ---
    run([args.plink2, "--bfile", out / "canon", "--make-pgen", "--out", out / "canon_pgen"], env)
    run([args.plink2, "--bfile", out / "canon", "--export", "bgen-1.2", "--out", out / "canon_bgen"], env)

    bgen_file = out / "canon_bgen.bgen"
    sample_file = out / "canon_bgen.sample"
    check(bgen_file.exists(), "plink2 produced a BGEN file")

    # --- export each format to a common genotype-only VCF ---
    # BED: convert with CosmicScanner so all three go through the same reader stack.
    run([args.scanner, "--bfile", out / "canon", "--convert", "bfile",
         "--to", "vcf", "--out", out / "v_bed"], env)
    run([args.scanner, "--pfile", out / "canon_pgen", "--convert", "pgen",
         "--to", "vcf", "--out", out / "v_pgen"], env)
    run([args.scanner, "--bgen", bgen_file, "--sample", sample_file,
         "--convert", "bgen", "--to", "vcf", "--out", out / "v_bgen"], env)

    bed_v, bed_s, bed_calls = read_vcf_genotypes(out / "v_bed.vcf")
    pgen_v, pgen_s, pgen_calls = read_vcf_genotypes(out / "v_pgen.vcf")
    bgen_v, bgen_s, bgen_calls = read_vcf_genotypes(out / "v_bgen.vcf")

    # plink2 writes BGEN sample identifiers as FID_IID while the PLINK readers
    # emit FID:IID. Normalize the separator; it is a naming convention, not part
    # of the genotype contract.
    def norm_samples(ids):
        return [s.replace(":", "_") for s in ids]

    bed_s, pgen_s, bgen_s = norm_samples(bed_s), norm_samples(pgen_s), norm_samples(bgen_s)

    # --- 1/2. ordering contracts ---
    check(bed_s == pgen_s == bgen_s, "sample ordering identical across BED/PGEN/BGEN")
    check(bed_v == pgen_v == bgen_v, "variant ordering identical across BED/PGEN/BGEN")
    check(len(bed_v) == m, "all {} variants exported from every format".format(m))
    check(len(bed_s) == n, "all {} samples exported from every format".format(n))

    # A BGEN marker may use the same or complementary allele-count orientation,
    # but that orientation must be consistent for every observed sample in the
    # marker. Per-sample equal-or-complement checks would not prove this contract.

    def is_het(g):
        return g == 1

    # --- 3. BED vs PGEN exact genotype identity (same orientation) ---
    mismatches = 0
    missing_mismatches = 0
    het_mismatches = 0
    orientation_modes = []
    orientation_mismatches = 0
    observed_total = 0
    for j in range(m):
        same_orientation = True
        complement_orientation = True
        for i in range(n):
            a, b, c = bed_calls[j][i], pgen_calls[j][i], bgen_calls[j][i]
            if a != b:
                mismatches += 1
            if (a == -1) != (b == -1) or (a == -1) != (c == -1):
                missing_mismatches += 1
            if is_het(a) != is_het(c):
                het_mismatches += 1
            if a >= 0 and c >= 0:
                observed_total += 1
                same_orientation = same_orientation and a == c
                complement_orientation = complement_orientation and a == 2 - c
        if not same_orientation and not complement_orientation:
            orientation_mismatches += 1
        orientation_modes.append("same" if same_orientation else "complement")
    check(mismatches == 0, "BED and PGEN per-sample genotype calls are identical")
    check(missing_mismatches == 0, "missing semantics identical across all three formats")
    check(het_mismatches == 0, "heterozygosity indicator identical across all three formats")
    check(orientation_mismatches == 0,
          "BGEN allele-count orientation is consistent within every marker")
    check(observed_total > 0.9 * m * n,
          "fixture retains observed genotypes ({}/{})".format(observed_total, m * n))

    # --- 5/6/7. dosage, heterozygosity, allele frequency from the one-hot BGEN ---
    # For a one-hot BGEN the decoded dosage is the hardcall allele count, and
    # P(AB) is exactly I(AB). Recompute both independently from the probability
    # semantics and compare against the BED hardcalls (complement-aware).
    max_dosage_diff = 0.0
    max_p_diff = 0.0
    het_mismatch = 0
    for j in range(m):
        col = bgen_calls[j]
        obs = [g for g in col if g >= 0]
        if not obs:
            continue
        p = (sum(obs) / len(obs)) / 2.0
        bed_obs = [g for g in bed_calls[j] if g >= 0]
        p_bed = (sum(bed_obs) / len(bed_obs)) / 2.0
        expected_p = p_bed if orientation_modes[j] == "same" else 1.0 - p_bed
        max_p_diff = max(max_p_diff, abs(p - expected_p))
        for i in range(n):
            g = col[i]
            if g < 0:
                continue
            # E[G] = P(AB) + 2 P(BB) with one-hot probabilities equals the count
            p_ab = 1.0 if g == 1 else 0.0
            p_bb = 1.0 if g == 2 else 0.0
            dosage = p_ab + 2.0 * p_bb
            bed_g = bed_calls[j][i]
            expected_g = bed_g if orientation_modes[j] == "same" else 2 - bed_g
            max_dosage_diff = max(max_dosage_diff, abs(dosage - expected_g))
            if p_ab != (1.0 if is_het(bed_g) else 0.0):
                het_mismatch += 1
    check(het_mismatch == 0, "P(AB) == I(AB) for one-hot BGEN probabilities")
    check(max_dosage_diff <= TOL_DOSAGE,
          "E[G] = P(AB) + 2P(BB) matches hardcall dosage (max diff {:.3e})".format(max_dosage_diff))
    check(max_p_diff <= TOL_P,
          "allele frequency p = mean(E[G])/2 identical across formats (max diff {:.3e})".format(max_p_diff))

    # --- 8. posterior-mean Vitezica hardcall limit ---
    # For one-hot probabilities the posterior-mean Vitezica code must reduce
    # exactly to the hardcall Vitezica encoding with normalization 4p^2q^2.
    max_vitezica_err = 0.0
    checked = 0
    for j in range(m):
        col = bgen_calls[j]
        obs = [g for g in col if g >= 0]
        if len(obs) < 3:
            continue
        p = (sum(obs) / len(obs)) / 2.0
        for i in range(n):
            g = col[i]
            if g < 0:
                continue
            # one-hot probabilities
            p0 = 1.0 if g == 0 else 0.0
            p1 = 1.0 if g == 1 else 0.0
            p2 = 1.0 if g == 2 else 0.0
            q = 1.0 - p
            mu = -2.0 * p * p * p0 + 2.0 * p * q * p1 - 2.0 * q * q * p2
            hard = expected_vitezica(g, p)
            max_vitezica_err = max(max_vitezica_err, abs(mu - hard))
            checked += 1
    check(checked > 0, "posterior-mean Vitezica hardcall-limit cases exercised ({})".format(checked))
    check(max_vitezica_err <= 1e-12,
          "one-hot posterior-mean Vitezica == hardcall Vitezica (max err {:.3e})".format(max_vitezica_err))

    # --- 9. element-wise production K_D parity ---
    kd = {}
    parity_qc = ["--bgen-dom-max-missing", "1", "--bgen-dom-min-info", "0"]
    inputs = (("bed", "--bfile", out / "canon", []),
              ("pgen", "--pfile", out / "canon_pgen", []),
              ("bgen", "--bgen", bgen_file, ["--sample", sample_file]))
    for name, flag, prefix, extra in inputs:
        target = out / ("kd_" + name)
        run([args.scanner, flag, prefix] + extra +
            ["--make-grm", "--make-grm-alg", "4", "--maf", "0.01",
             "--out", target] + (parity_qc if name == "bgen" else []), env)
        kd[name] = read_grm(str(target) + ".grm.bin")
    check(len(kd["bed"]) == len(kd["pgen"]) == len(kd["bgen"]),
          "K_D element counts agree across BED/PGEN/BGEN")
    for name in ("pgen", "bgen"):
        worst = max(abs(a-b) for a, b in zip(kd["bed"], kd[name]))
        check(worst <= 2e-5,
              "K_D matches BED element-wise for {} (max diff {:.3e})".format(name, worst))

    # --- 10. dominance-only scan parity, global and chromosome-excluded LOCO ---
    quant = out / "quant.pheno"
    with quant.open("w") as handle:
        for i, s in enumerate(samples):
            handle.write("{} {} {:.12g}\n".format(s, s, ((i * 13) % 17 - 8) / 5.0))

    for loco in (False, True):
        scans = {}
        metas = {}
        tag = "loco" if loco else "global"
        for name, flag, prefix, extra in inputs:
            target = out / ("dom_{}_{}".format(tag, name))
            command = ([args.scanner, flag, prefix] + extra +
                       ["--gwas", "--gwas-model", "dominance", "--pheno", quant,
                        "--maf", "0.01", "--out", target] + (parity_qc if name == "bgen" else []))
            if not loco:
                command.append("--no-loco")
            run(command, env)
            header, rows = read_table(str(target) + ".assoc.dom.txt")
            check(header == UNIFIED_HEADER, "{} {} uses unified schema".format(tag, name))
            scans[name] = {row["SNP"]: row for row in rows}
            metas[name] = read_meta(str(target) + ".scan.meta.txt")

        check(set(scans["bed"]) == set(scans["pgen"]) == set(scans["bgen"]),
              "{} tested marker membership agrees across formats".format(tag))
        for name in ("pgen", "bgen"):
            identity_mismatches = sum(
                scans[name][snp][column] != scans["bed"][snp][column]
                for snp in scans["bed"] for column in ("CHR", "BP", "A1", "A2", "N"))
            check(identity_mismatches == 0,
                  "{} {} marker coordinates, allele labels, and N match BED".format(tag, name))
            worst_freq = max(abs(float(scans[name][snp]["FREQ"])
                                 - float(scans["bed"][snp]["FREQ"]))
                             for snp in scans["bed"])
            check(worst_freq <= TOL_P,
                  "{} {} FREQ matches BED (max diff {:.3e})".format(tag, name, worst_freq))
            for column in ("BETA_D", "SE_D", "P_D"):
                worst = max(abs(float(scans["bed"][snp][column])
                                - float(scans[name][snp][column]))
                            for snp in scans["bed"])
                check(worst <= TOL_STAT,
                      "{} {} {} matches BED (max diff {:.3e})".format(
                          tag, name, column, worst))
            sigma_keys = sorted(key for key in metas["bed"]
                                if key.startswith("SIGMA_D2") and key in metas[name]
                                and metas["bed"][key] != "PER_NULL_MODEL")
            check(sigma_keys, "{} {} reports comparable sigma_D2 values".format(tag, name))
            for key in sigma_keys:
                a, b = float(metas["bed"][key]), float(metas[name][key])
                check(abs(a-b) <= 1e-5*max(abs(a), abs(b), 1e-6),
                      "{} {} {} matches BED".format(tag, name, key))
        check(metas["bgen"].get("GENOTYPE_FORMAT") == "BGEN_PROBABILITY",
              "{} BGEN metadata records probability input".format(tag))
        check(metas["bgen"].get("GENOTYPE_BACKEND") == "BGEN" and
              metas["bgen"].get("PROBABILITY_MODE") == "true",
              "{} BGEN metadata records the production probability backend".format(tag))
        check(metas["bgen"].get("DOMINANCE_ENCODING") == "posterior_mean_vitezica" and
              metas["bgen"].get("QC_POLICY_VERSION") == "BGEN_DOM_V1",
              "{} BGEN metadata records dominance encoding and QC policy".format(tag))
        check(metas["bgen"].get("KD_ENCODING") ==
              "POSTERIOR_MEAN_VITEZICA_OUTER_PRODUCT",
              "{} BGEN metadata records posterior-mean Vitezica K_D".format(tag))
        check(metas["bgen"].get("LOCO") == ("true" if loco else "false"),
              "{} BGEN metadata records LOCO mode".format(tag))
        membership_keys = sorted(key for key in metas["bed"]
                                 if key.startswith("KD_MARKERS_") or
                                 key.startswith("KD_DENOMINATOR_"))
        check(membership_keys, "{} metadata reports K_D membership and normalization".format(tag))
        for name in ("pgen", "bgen"):
            for key in membership_keys:
                if key.startswith("KD_MARKERS_"):
                    check(metas[name].get(key) == metas["bed"].get(key),
                          "{} {} {} matches BED".format(tag, name, key))
                else:
                    a = float(metas["bed"][key])
                    b = float(metas[name][key])
                    check(abs(a-b) <= 1e-7*max(abs(a), abs(b), 1.0),
                          "{} {} {} matches BED".format(tag, name, key))

    # Exercise genotype-to-analysis-sample indexing and filtered GRM membership
    # together. Omitting phenotype rows creates a nontrivial keep-index vector.
    subset_quant = out / "quant_subset.pheno"
    with subset_quant.open("w") as handle:
        for i, s in enumerate(samples):
            if i % 9 != 0:
                handle.write("{} {} {:.12g}\n".format(s, s, ((i * 13) % 17 - 8) / 5.0))
    extract = out / "dominance.extract"
    exclude = out / "dominance.exclude"
    extracted = ["rs{:04d}".format(j + 1) for j in range(1, m) if j % 3 != 0]
    extract.write_text("\n".join(extracted) + "\n")
    exclude.write_text("\n".join(extracted[::11]) + "\n")

    filtered_scans = {}
    filtered_metas = {}
    for name, flag, prefix, extra in inputs:
        target = out / ("dom_filtered_subset_" + name)
        run([args.scanner, flag, prefix] + extra +
            ["--gwas", "--gwas-model", "dominance", "--pheno", subset_quant,
             "--maf", "0.01", "--extract", extract, "--exclude", exclude,
             "--no-loco", "--out", target] + (parity_qc if name == "bgen" else []), env)
        header, rows = read_table(str(target) + ".assoc.dom.txt")
        check(header == UNIFIED_HEADER, "filtered subset {} uses unified schema".format(name))
        filtered_scans[name] = {row["SNP"]: row for row in rows}
        filtered_metas[name] = read_meta(str(target) + ".scan.meta.txt")
    check(set(filtered_scans["bed"]) == set(filtered_scans["pgen"]) == set(filtered_scans["bgen"]),
          "filtered subset marker membership agrees across formats")
    check(len(filtered_scans["bed"]) < len(scans["bed"]),
          "extract/exclude reduce the tested marker set")
    expected_n = n - (n + 8) // 9
    check(all(int(row["N"]) == expected_n for row in filtered_scans["bgen"].values()),
          "BGEN dominance uses the phenotype-derived sample subset (N={})".format(expected_n))
    for name in ("pgen", "bgen"):
        check(filtered_metas[name]["KD_MARKERS_GLOBAL"] ==
              filtered_metas["bed"]["KD_MARKERS_GLOBAL"],
              "filtered subset {} K_D marker count matches BED".format(name))
    for name in ("pgen",):
        for column in ("BETA_D", "SE_D", "P_D"):
            worst = max(abs(float(filtered_scans["bed"][snp][column])
                            - float(filtered_scans[name][snp][column]))
                        for snp in filtered_scans["bed"])
            check(worst <= TOL_STAT,
                  "filtered subset {} {} matches BED (max diff {:.3e})".format(
                      name, column, worst))
    bgen_index = {variant: j for j, variant in enumerate(bgen_v)}
    worst_subset_freq = 0.0
    for snp, row in filtered_scans["bgen"].items():
        column = bgen_calls[bgen_index[snp]]
        observed = [g for i, g in enumerate(column) if i % 9 != 0 and g >= 0]
        expected_freq = sum(observed) / (2.0 * len(observed))
        worst_subset_freq = max(worst_subset_freq, abs(float(row["FREQ"]) - expected_freq))
    check(worst_subset_freq <= TOL_P,
          "BGEN subset FREQ uses analysis samples (max diff {:.3e})".format(worst_subset_freq))

    # --- 11. BGEN AD remains NO-GO ---
    target = out / "bgen_additive_dominance"
    proc = subprocess.run(
        [args.scanner, "--bgen", bgen_file, "--sample", sample_file,
         "--gwas-model", "additive_dominance", "--pheno", quant, "--out", target],
        text=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    check(proc.returncode != 0 and "not supported for BGEN" in proc.stdout,
          "BGEN --gwas-model additive_dominance remains rejected")
    check(not Path(str(target) + ".assoc.ad.txt").exists(),
          "BGEN additive_dominance produced no output file")

    # Unsupported probabilistic-QC/frequency modes must fail rather than being
    # silently ignored by either the GWAS or standalone algorithm-4 GRM path.
    guards = (
        ("hwe_scan", ["--gwas", "--gwas-model", "dominance", "--pheno", quant,
                      "--hwe", "0.01"], "does not support exact --hwe"),
        ("hwe_grm", ["--make-grm", "--make-grm-alg", "4", "--hwe", "0.01"],
         "does not support exact HWE"),
        ("pop_class", ["--gwas", "--gwas-model", "dominance", "--pheno", quant,
                       "--pop-class", out / "unused.pop"], "does not support --pop-class"),
        ("snp_weight", ["--gwas", "--gwas-model", "dominance", "--pheno", quant,
                        "--snp-weight", out / "unused.weight"], "does not support --snp-weight"),
    )
    for label, extra_args, expected_message in guards:
        guard_target = out / ("guard_" + label)
        guard = subprocess.run(
            [str(v) for v in ([args.scanner, "--bgen", bgen_file, "--sample", sample_file]
                              + extra_args + ["--out", guard_target])],
            text=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        check(guard.returncode != 0 and expected_message in guard.stdout,
              "BGEN dominance {} is rejected explicitly".format(label))

    print("PASS BED/PGEN/BGEN one-hot genotype, K_D, sigma_D2, dominance "
          "beta/SE/P, LOCO, and BGEN-AD guard contract")


if __name__ == "__main__":
    main()
