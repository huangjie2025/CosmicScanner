#!/usr/bin/env python3
import argparse
import csv
import math
import os
import pathlib
import subprocess


def run(command, expect_success=True):
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if (completed.returncode == 0) != expect_success:
        raise RuntimeError("Command status mismatch:\n" + " ".join(map(str, command)) + "\n" + completed.stdout)
    return completed.stdout


def read_rows(path):
    with path.open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--plink2", default=os.environ.get("PLINK2_BIN", "plink2"))
    args = parser.parse_args()
    binary = str(pathlib.Path(args.bin).resolve())
    out = pathlib.Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    n = 120
    variants = []
    for j in range(8):
        row = []
        for i in range(n):
            if j == 1:
                value = variants[0][i] if i % 10 else (variants[0][i] + 1) % 3
            elif j == 3:
                value = variants[2][i] if i % 8 else (variants[2][i] + 2) % 3
            else:
                value = (i * (j + 3) + i // (j + 2) + j * 7) % 3
            row.append(value)
        variants.append(row)

    prefix = out / "reference"
    with pathlib.Path(str(prefix) + ".fam").open("w") as target:
        for i in range(n):
            target.write(f"F{i + 1} I{i + 1} 0 0 0 -9\n")
    with pathlib.Path(str(prefix) + ".bim").open("w") as target:
        for j in range(len(variants)):
            target.write(f"1 rs{j + 1} 0 {1000 + j * 100} A C\n")
    bed = bytearray(b"\x6c\x1b\x01")
    for row in variants:
        packed = bytearray((n + 3) // 4)
        for i, value in enumerate(row):
            packed[i // 4] |= (0b11, 0b10, 0b00)[value] << (2 * (i % 4))
        bed.extend(packed)
    pathlib.Path(str(prefix) + ".bed").write_bytes(bed)

    pgen = out / "reference_pgen"
    bgen = out / "reference_bgen"
    run([args.plink2, "--bfile", str(prefix), "--make-pgen", "--out", str(pgen)])
    run([args.plink2, "--bfile", str(prefix), "--export", "bgen-1.2", "--out", str(bgen)])

    vcf = out / "reference.vcf"
    with vcf.open("w") as target:
        target.write("##fileformat=VCFv4.2\n")
        target.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                     "\t".join(f"I{i + 1}" for i in range(n)) + "\n")
        for j, row in enumerate(variants):
            calls = ("0/0", "0/1", "1/1")
            target.write(f"1\t{1000 + j * 100}\trs{j + 1}\tC\tA\t.\tPASS\t.\tGT\t" +
                         "\t".join(calls[value] for value in row) + "\n")

    betas = [0.8, 0.56, 0.7, 0.45, 0.03, -0.02, 0.01, 0.0]
    se = 0.1
    summary = out / "locus.cosmic.sumstats"
    swapped = out / "locus_swapped.cosmic.sumstats"
    strand = out / "locus_strand.cosmic.sumstats"
    strand_swapped = out / "locus_strand_swapped.cosmic.sumstats"
    ambiguous = out / "locus_ambiguous.cosmic.sumstats"
    header = "CHR\tSNP\tPOS\tA1\tA2\tBETA\tSE\tP\tN\tAF\tLOG10P\tMODEL\tZ\tWEIGHT\tINFO\tWARNING\n"
    orientations = ((summary, "base"), (swapped, "swap"), (strand, "strand"),
                    (strand_swapped, "strand_swap"), (ambiguous, "ambiguous"))
    for path, orientation in orientations:
        with path.open("w") as target:
            target.write("# CosmicScanner Summary Statistics v1.0\n#\n")
            target.write(header)
            for j, beta in enumerate(betas):
                value = -beta if orientation in ("swap", "strand_swap") else beta
                a1, a2 = {"base": ("A", "C"), "swap": ("C", "A"),
                          "strand": ("T", "G"), "strand_swap": ("G", "T")}.get(
                              orientation, ("A", "T") if j == 0 else ("A", "C"))
                z = value / se
                p = math.erfc(abs(z) / math.sqrt(2.0))
                target.write(f"1\trs{j + 1}\t{1000 + j * 100}\t{a1}\t{a2}\t{value}\t{se}\t{p}\t{n}\t0.3\t"
                             f"{-math.log10(p) if p else 300}\tLMM\t{z}\t100\t1\t.\n")

    common = ["--fine-map", "--finemap-region", "1:900-1800", "--finemap-prior-var", "0.04",
              "--credible-level", "0.95", "--finemap-max-signals", "2", "--finemap-signal-p", "1e-4"]
    runs = {
        "bed": ["--bfile", str(prefix), "--finemap-summary", str(summary)],
        "pgen": ["--pfile", str(pgen), "--finemap-summary", str(summary)],
        "bgen": ["--bgen", str(bgen) + ".bgen", "--finemap-summary", str(summary)],
        "vcf": ["--vcf", str(vcf), "--finemap-summary", str(summary)],
        "swap": ["--bfile", str(prefix), "--finemap-summary", str(swapped)],
        "strand": ["--bfile", str(prefix), "--finemap-summary", str(strand)],
        "strand_swap": ["--bfile", str(prefix), "--finemap-summary", str(strand_swapped)],
    }
    tables = {}
    credible_tables = {}
    for name, source in runs.items():
        output_prefix = out / name
        run([binary, *common, *source, "--out", str(output_prefix)])
        rows = read_rows(pathlib.Path(str(output_prefix) + ".finemap.tsv"))
        credible = read_rows(pathlib.Path(str(output_prefix) + ".credible_sets.tsv"))
        if not rows or not credible:
            raise RuntimeError(f"{name}: missing fine-mapping or credible-set rows")
        tables[name] = {(row["SIGNAL"], row["SNP"]): row for row in rows}
        credible_tables[name] = {(row["SIGNAL"], row["SNP"]): row for row in credible}

    parity_runs = ("pgen", "bgen", "vcf", "swap", "strand", "strand_swap")
    if any(tables["bed"].keys() != tables[name].keys() for name in parity_runs):
        raise RuntimeError("Fine-mapping format/orientation keys differ")
    if any(credible_tables["bed"].keys() != credible_tables[name].keys() for name in parity_runs):
        raise RuntimeError("Fine-mapping credible-set membership differs")
    for key, bed_row in tables["bed"].items():
        for name in parity_runs:
            row = tables[name][key]
            for column in ("PIP", "ABF", "CUM_POSTERIOR"):
                if abs(float(bed_row[column]) - float(row[column])) > 1e-9 * max(1.0, abs(float(bed_row[column]))):
                    raise RuntimeError(f"{name} parity failure for {key} {column}")
    for key, bed_row in credible_tables["bed"].items():
        for name in parity_runs:
            row = credible_tables[name][key]
            for column in ("PIP", "CUM_POSTERIOR", "LD_WITH_LEAD"):
                if abs(float(bed_row[column]) - float(row[column])) > 1e-9 * max(1.0, abs(float(bed_row[column]))):
                    raise RuntimeError(f"{name} credible-set parity failure for {key} {column}")

    first = tables["bed"][("1", "rs1")]
    variance = se * se
    shrinkage = 0.04 / (variance + 0.04)
    expected_abf = math.sqrt(variance / (variance + 0.04)) * math.exp(0.5 * (betas[0] / se) ** 2 * shrinkage)
    if abs(float(first["ABF"]) - expected_abf) > 1e-9 * expected_abf:
        raise RuntimeError("independent Wakefield ABF oracle mismatch")
    for signal in ("1", "2"):
        total = sum(float(row["PIP"]) for (row_signal, _), row in tables["bed"].items() if row_signal == signal)
        if abs(total - 1.0) > 1e-9:
            raise RuntimeError(f"signal {signal} PIP normalization failure")

    metadata = pathlib.Path(str(out / "bed") + ".finemap.meta.tsv").read_text()
    if "method\t" not in metadata or "prior_variance\t0.04" not in metadata or "ld_source\tBED" not in metadata:
        raise RuntimeError("fine-mapping metadata contract is incomplete")

    one_signal = out / "one_signal"
    run([binary, "--fine-map", "--bfile", str(prefix), "--finemap-summary", str(summary),
         "--finemap-region", "1:900-1800", "--finemap-prior-var", "0.04",
         "--finemap-max-signals", "1", "--finemap-signal-p", "1e-4", "--out", str(one_signal)])
    one_rows = read_rows(pathlib.Path(str(one_signal) + ".finemap.tsv"))
    if {row["SIGNAL"] for row in one_rows} != {"1"}:
        raise RuntimeError("One-signal fine-mapping produced an unexpected signal")

    ambiguous_prefix = out / "ambiguous"
    run([binary, *common, "--bfile", str(prefix), "--finemap-summary", str(ambiguous),
         "--out", str(ambiguous_prefix)])
    ambiguous_rows = read_rows(pathlib.Path(str(ambiguous_prefix) + ".finemap.tsv"))
    if any(row["SNP"] == "rs1" for row in ambiguous_rows):
        raise RuntimeError("Unsafe ambiguous SNP was not rejected")

    run([binary, "--fine-map", "--finemap-summary", str(summary), "--finemap-region", "1:900-1800",
         "--out", str(out / "missing_ld")], expect_success=False)
    print("PASS fine-mapping orientation, BED/VCF LD, one/two-signal, PIP, and credible-set parity")


if __name__ == "__main__":
    main()
