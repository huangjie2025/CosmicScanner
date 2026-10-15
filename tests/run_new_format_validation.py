#!/usr/bin/env python3
import argparse
import math
import os
import pathlib
import struct
import subprocess


def run(command):
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError("Command failed:\n" + " ".join(map(str, command)) + "\n" + completed.stdout)
    return completed.stdout


def read_table(path):
    with path.open() as source:
        header = source.readline().split()
        return [dict(zip(header, line.split())) for line in source if line.strip()]


def compare_tables(left, right, keys, tolerance=1e-10):
    if len(left) != len(right):
        raise RuntimeError(f"Association row count mismatch: {len(left)} != {len(right)}")
    for row_left, row_right in zip(left, right):
        if row_left["SNP"] != row_right["SNP"]:
            raise RuntimeError("Association variant order mismatch")
        for key in keys:
            a = float(row_left[key])
            b = float(row_right[key])
            if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > tolerance * max(1.0, abs(a), abs(b)):
                raise RuntimeError(f"Association mismatch for {row_left['SNP']} {key}: {a} != {b}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--plink2", default=os.environ.get(
        "PLINK2_BIN", "plink2"))
    args = parser.parse_args()
    binary = pathlib.Path(args.bin).resolve()
    output = pathlib.Path(args.out).resolve()
    output.mkdir(parents=True, exist_ok=True)

    n_samples = 40
    n_variants = 30
    samples = [f"S{i + 1:03d}" for i in range(n_samples)]
    variants = []
    for j in range(n_variants):
        calls = [((i * 7 + j * 3 + j * j + i * j) % 3) for i in range(n_samples)]
        if j % 7 == 0:
            calls[(j * 5) % n_samples] = None
        variants.append(calls)

    hmp = output / "canonical.hmp.txt"
    with hmp.open("w") as target:
        target.write("rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode " + " ".join(samples) + "\n")
        for j, calls in enumerate(variants):
            chrom = "1" if j < n_variants // 2 else "2"
            encoded = ["NN" if value is None else ("AA", "AC", "CC")[value] for value in calls]
            target.write(f"rs{j + 1:03d} A/C {chrom} {(j + 1) * 100} + NA NA NA NA NA NA " + " ".join(encoded) + "\n")

    matrix = output / "canonical.geno"
    metadata = output / "canonical.var.info"
    with matrix.open("w") as target:
        target.write("VariantID " + " ".join(samples) + "\n")
        for j, calls in enumerate(variants):
            target.write(f"rs{j + 1:03d} " + " ".join("NA" if value is None else str(value) for value in calls) + "\n")
    with metadata.open("w") as target:
        target.write("ID CHROM POS REF ALT\n")
        for j in range(n_variants):
            target.write(f"rs{j + 1:03d} {'1' if j < n_variants // 2 else '2'} {(j + 1) * 100} A C\n")

    plink = output / "canonical"
    with pathlib.Path(str(plink) + ".fam").open("w") as target:
        for sample in samples:
            target.write(f"{sample} {sample} 0 0 0 -9\n")
    with pathlib.Path(str(plink) + ".bim").open("w") as target:
        for j in range(n_variants):
            target.write(f"{'1' if j < n_variants // 2 else '2'} rs{j + 1:03d} 0 {(j + 1) * 100} A C\n")
    bed = bytearray(b"\x6c\x1b\x01")
    for calls in variants:
        packed = bytearray((n_samples + 3) // 4)
        for i, value in enumerate(calls):
            code = 0b01 if value is None else (0b11, 0b10, 0b00)[value]
            packed[i // 4] |= code << (2 * (i % 4))
        bed.extend(packed)
    pathlib.Path(str(plink) + ".bed").write_bytes(bed)

    pgen = output / "canonical_pgen"
    bgen = output / "canonical_bgen"
    run([args.plink2, "--bfile", str(plink), "--make-pgen", "--out", str(pgen)])
    run([args.plink2, "--bfile", str(plink), "--export", "bgen-1.2", "--out", str(bgen)])

    vcf = output / "canonical.vcf"
    with vcf.open("w") as target:
        target.write("##fileformat=VCFv4.2\n")
        target.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for j, calls in enumerate(variants):
            encoded = ["./." if value is None else ("0/0", "0/1", "1/1")[value] for value in calls]
            target.write(f"{'1' if j < n_variants // 2 else '2'}\t{(j + 1) * 100}\trs{j + 1:03d}\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(encoded) + "\n")

    quantitative = output / "quantitative.pheno"
    binary_pheno = output / "binary.pheno"
    vcf_quantitative = output / "vcf.quantitative.pheno"
    vcf_binary = output / "vcf.binary.pheno"
    with quantitative.open("w") as qtarget, binary_pheno.open("w") as btarget, \
         vcf_quantitative.open("w") as vqtarget, vcf_binary.open("w") as vbtarget:
        for i, sample in enumerate(samples):
            first = variants[1][i]
            y = 0.8 * first + ((i * 13) % 17 - 8) / 5.0
            qtarget.write(f"{sample} {sample} {y:.12g}\n")
            btarget.write(f"{sample} {sample} {1 if i % 5 == 0 else 0}\n")
            vqtarget.write(f"0 {sample} {y:.12g}\n")
            vbtarget.write(f"0 {sample} {1 if i % 5 == 0 else 0}\n")

    formats = {
        "hmp": (["--hmp", str(hmp)], quantitative, binary_pheno),
        "text": (["--geno-matrix", str(matrix), "--geno-var", str(metadata)], quantitative, binary_pheno),
        "bed": (["--bfile", str(plink)], quantitative, binary_pheno),
        "pgen": (["--pfile", str(pgen)], quantitative, binary_pheno),
        "bgen": (["--bgen", str(bgen) + ".bgen", "--sample", str(bgen) + ".sample"], quantitative, binary_pheno),
        "vcf": (["--vcf", str(vcf)], vcf_quantitative, vcf_binary),
    }

    grm_values = {}
    for format_name, (format_args, _, _) in formats.items():
        prefix = output / f"{format_name}_grm"
        run([str(binary), *format_args, "--make-grm", "--maf", "0", "--out", str(prefix)])
        raw = pathlib.Path(str(prefix) + ".grm.bin").read_bytes()
        grm_values[format_name] = struct.unpack(f"<{len(raw) // 4}f", raw)
    reference_grm = grm_values["hmp"]
    for format_name, values in grm_values.items():
        if len(values) != len(reference_grm) or max(abs(a - b) for a, b in zip(reference_grm, values)) > 2e-5:
            raise RuntimeError(f"HapMap/{format_name} GRM mismatch")

    for name, extra, suffix, keys in [
        ("lm", ["--lm", "--lm-mode", "4", "--pheno", str(quantitative)], ".lm.assoc.txt", ["BETA", "SE", "P"]),
        ("lmm", ["--gwas", "--no-loco", "--pheno", str(quantitative)], ".assoc.txt", ["BETA_A", "SE_A", "P_A"]),
        ("loco", ["--gwas", "--pheno", str(quantitative)], ".assoc.txt", ["BETA_A", "SE_A", "P_A"]),
        ("saige", ["--saige-style", "--biobank-mode", "--pheno", str(binary_pheno)], ".cosmic.assoc", ["BETA", "SE", "P"]),
    ]:
        reference = None
        selected_formats = formats
        for format_name, (format_args, qpheno, bpheno) in selected_formats.items():
            prefix = output / f"{format_name}_{name}"
            format_extra = [str(qpheno) if value == str(quantitative)
                            else str(bpheno) if value == str(binary_pheno)
                            else value for value in extra]
            run([str(binary), *format_args, *format_extra, "--maf", "0", "--out", str(prefix)])
            table = read_table(pathlib.Path(str(prefix) + suffix))
            if reference is None:
                reference = table
            else:
                tolerance = 2e-5 if format_name in ("bed", "pgen", "bgen", "vcf") else 1e-10
                compare_tables(reference, table, keys, tolerance=tolerance)

    hardcall_formats = {key: value for key, value in formats.items() if key != "bgen"}
    for model in ("dominant", "recessive", "heterozygote"):
        for analysis, model_args, suffix in (
                ("lm", ["--lm", "--lm-mode", "1"], ".lm.assoc.txt"),
                ("lmm", ["--gwas", "--no-loco"], ".assoc.txt")):
            reference = None
            for format_name, (format_args, qpheno, _) in hardcall_formats.items():
                prefix = output / f"{format_name}_{analysis}_{model}"
                run([str(binary), *format_args, *model_args, "--snp-model", model,
                     "--pheno", str(qpheno), "--maf", "0", "--out", str(prefix)])
                table = read_table(pathlib.Path(str(prefix) + suffix))
                if reference is None:
                    reference = table
                else:
                    tolerance = 2e-5 if format_name in ("bed", "pgen", "vcf") else 1e-10
                    keys = ["BETA", "SE", "P"] if analysis == "lm" else ["BETA_A", "SE_A", "P_A"]
                    compare_tables(reference, table, keys, tolerance=tolerance)

    bgen_failure = subprocess.run(
        [str(binary), *formats["bgen"][0], "--lm", "--snp-model", "dominant",
         "--pheno", str(quantitative), "--out", str(output / "bgen_nonadditive")],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if bgen_failure.returncode == 0 or "hardcalls" not in bgen_failure.stdout:
        raise RuntimeError("BGEN non-additive coding was not rejected explicitly")

    for format_name in ("bed", "hmp", "text", "bgen"):
        format_args, qpheno, _ = formats[format_name]
        prefix = output / f"{format_name}_hwe_policy"
        run([str(binary), *format_args, "--gwas", "--no-loco", "--hwe", "1e-12",
             "--maf", "0", "--pheno", str(qpheno), "--out", str(prefix)])
        metadata = pathlib.Path(str(prefix) + ".scan.meta.txt").read_text()
        expected = "HWE_APPLIED\tfalse" if format_name == "bgen" else "HWE_APPLIED\ttrue"
        if expected not in metadata:
            raise RuntimeError(f"Incorrect explicit HWE policy for {format_name}")

    diagnostic = read_table(pathlib.Path(str(output / "hmp_saige") + ".saige_style.variants.tsv"))
    required = {"AF", "MAC", "MAC_BIN", "VARIANCE_RATIO", "LOG10P", "SPA_USED"}
    if not diagnostic or not required.issubset(diagnostic[0]):
        raise RuntimeError("SAIGE per-variant diagnostics are incomplete")
    metadata_text = pathlib.Path(str(output / "hmp_saige") + ".saige_style.meta.tsv").read_text()
    if "spa_marker_count\t" not in metadata_text or "mac_bin_0\t" not in metadata_text:
        raise RuntimeError("SAIGE run metadata is incomplete")

    print("PASS BED/PGEN/BGEN/VCF/HapMap/text GRM, LM, LMM, LOCO, and SAIGE parity")


if __name__ == "__main__":
    main()
