#!/usr/bin/env python3
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

from run_association_r2 import fixture


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_results(path, expected_method):
    with path.open(encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == ["CHR", "START", "END", "VARIANT_SET_ID", "METHOD",
                                      "STATISTIC", "P_VALUE"],
                f"Unexpected rare-variant schema: {reader.fieldnames}")
        rows = list(reader)
    require(rows, f"No rare-variant result rows: {path}")
    for row in rows:
        require(row["METHOD"] == expected_method, f"Unexpected method label: {row}")
        require(row["CHR"] and int(row["START"]) <= int(row["END"]),
                f"Invalid set coordinates: {row}")
        statistic = float(row["STATISTIC"])
        p_value = float(row["P_VALUE"])
        statistic_valid = not math.isnan(statistic) if expected_method == "ACAT_V" else math.isfinite(statistic)
        require(statistic_valid and math.isfinite(p_value) and 0.0 <= p_value <= 1.0,
                f"Invalid rare-variant inference: {row}")
    return {row["VARIANT_SET_ID"]: row for row in rows}


def run(binary, source, output, data, method):
    flag = "--" + method
    args = [str(binary), *map(str, source), "--out", str(output), "--pheno", str(data["pheno"]),
            "--set-file", str(data["set_file"]), "--maf", "0", "--hwe", "0",
            "--rare-max-maf", "0.5", "--rare-min-mac", "1", flag]
    if method in ("burden", "skat"):
        args.append("--gwas")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    completed = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               check=False, env=env)
    require(completed.returncode == 0, "Command failed:\n" + " ".join(args) + "\n" + completed.stdout)
    method_name = {"burden": "BURDEN", "skat": "SKAT", "skat-o": "SKAT_O_CANONICAL",
                   "acat-v": "ACAT_V"}[method]
    return read_results(Path(str(output) + ".rare.txt"), method_name)


def run_logistic(binary, data, root, flag, expected_model):
    phenotype = root / "binary.pheno"
    phenotype.write_text("".join(
        f"{data['samples'][i][0]} {data['samples'][i][1]} {i % 2}\n"
        for i in reversed(range(data["n"]))), encoding="ascii")
    output = root / expected_model.lower()
    command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(phenotype),
               "--out", str(output), "--maf", "0", *flag]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               check=False, env=env)
    require(completed.returncode == 0, "Command failed:\n" + " ".join(command) + "\n" + completed.stdout)
    with Path(str(output) + ".cosmic.assoc").open(encoding="ascii", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    require(rows and all(row["MODEL"] == expected_model for row in rows),
            f"{expected_model}: missing or mislabelled output")
    require(all(math.isfinite(float(row["BETA"])) and float(row["SE"]) > 0 and
                0 <= float(row["P"]) <= 1 for row in rows), f"{expected_model}: invalid inference")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="r8-b2-", dir=args.out) as temporary:
        root = Path(temporary)
        data = fixture(root)
        with data["set_file"].open("a", encoding="ascii") as stream:
            stream.write("empty_set not_present\n")
        run_logistic(args.bin, data, root, ("--logistic", "score"), "LOGISTIC_SCORE")
        run_logistic(args.bin, data, root, ("--logistic", "wald"), "LOGISTIC_WALD")
        run_logistic(args.bin, data, root, ("--logistic-glmm",), "LOGISTIC_GLMM")
        vcf = root / "fixture.vcf"
        lines = ["##fileformat=VCFv4.2",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                 "\t".join(iid for _, iid in data["samples"])]
        for j, (chrom, marker, position) in enumerate(data["variants"]):
            genotypes = []
            for value in data["calls"][j]:
                gt = "./." if value is None else ("0/0", "0/1", "1/1")[value]
                genotypes.append("20:" + gt)
            lines.append(f"{chrom}\t{position}\t{marker}\tG\tA\t.\tPASS\t.\tDP:GT\t" + "\t".join(genotypes))
        vcf.write_text("\n".join(lines) + "\n", encoding="ascii")
        vcf_pheno = root / "phenotype.vcf.txt"
        vcf_pheno.write_text("".join(
            f"0 {data['samples'][i][1]} {data['y'][i]:.8f}\n" for i in reversed(range(data["n"]))),
            encoding="ascii")
        vcf_data = dict(data, pheno=vcf_pheno)

        for method in ("burden", "skat", "skat-o", "acat-v"):
            bed = run(args.bin, ("--bfile", data["prefix"]), root / (method + "-bed"), data, method)
            direct_vcf = run(args.bin, ("--vcf", vcf), root / (method + "-vcf"), vcf_data, method)
            require("empty_set" not in bed and "empty_set" not in direct_vcf,
                    f"{method}: untested empty set was emitted")
            require(bed.keys() == direct_vcf.keys(), f"{method}: BED/VCF set mismatch")
            for set_id in bed:
                for field in ("CHR", "START", "END", "METHOD"):
                    require(bed[set_id][field] == direct_vcf[set_id][field],
                            f"{method} {set_id}: BED/VCF {field} mismatch")
                for field in ("STATISTIC", "P_VALUE"):
                    left, right = float(bed[set_id][field]), float(direct_vcf[set_id][field])
                    require(math.isclose(left, right, rel_tol=1e-7, abs_tol=1e-9),
                            f"{method} {set_id}: BED/VCF {field} mismatch: {left} vs {right}")
        print("R8-B2 BED/VCF rare-variant parity passed")


if __name__ == "__main__":
    main()
