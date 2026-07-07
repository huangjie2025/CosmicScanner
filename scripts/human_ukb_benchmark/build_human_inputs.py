#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


TRAITS = [
    ("21001", "BMI"),
    ("21002", "Weight"),
    ("23105", "BasalMetabolicRate"),
    ("48", "WaistCircumference"),
    ("49", "HipCircumference"),
    ("50", "StandingHeight"),
    ("30000", "WBC"),
    ("30010", "RBC"),
    ("30080", "PlateletCount"),
    ("30150", "EosinophilCount"),
]


def read_psam(path):
    samples = []
    with path.open() as handle:
        header = handle.readline().strip().split()
        try:
            fid_idx = header.index("#FID")
        except ValueError:
            fid_idx = header.index("FID")
        iid_idx = header.index("IID")
        for line in handle:
            if not line.strip():
                continue
            fields = line.strip().split()
            samples.append((fields[fid_idx], fields[iid_idx]))
    return samples


def read_trait(path):
    values = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        value_col = [x for x in reader.fieldnames if x != "f.eid"][0]
        for row in reader:
            raw = row[value_col].strip()
            if raw in {"", "NA", "NaN", "nan"}:
                values[row["f.eid"]] = None
            else:
                values[row["f.eid"]] = float(raw)
    return values


def zscore(values):
    finite = [x for x in values if x is not None and math.isfinite(x)]
    mean = sum(finite) / len(finite)
    var = sum((x - mean) ** 2 for x in finite) / (len(finite) - 1)
    sd = math.sqrt(var)
    return [(None if x is None else (x - mean) / sd) for x in values], mean, sd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ucr-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    geno_dir = args.ucr_dir / "genotype_final"
    pheno_dir = args.ucr_dir / "phenotype"
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = read_psam(geno_dir / "CHR1_final.psam")
    trait_maps = {name: read_trait(pheno_dir / f"{field}.txt") for field, name in TRAITS}

    raw_rows = []
    for fid, iid in samples:
        vals = [trait_maps[name].get(iid) for _, name in TRAITS]
        raw_rows.append((fid, iid, vals))

    complete_rows = [row for row in raw_rows if all(v is not None and math.isfinite(v) for v in row[2])]
    columns = list(zip(*[row[2] for row in complete_rows]))
    z_columns = []
    stats = []
    for (_, name), values in zip(TRAITS, columns):
        z_values, mean, sd = zscore(values)
        z_columns.append(z_values)
        stats.append((name, len(values), mean, sd))

    pheno_path = out_dir / "human_ukb_10traits.standardized.pheno"
    with pheno_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["FID", "IID"] + [name for _, name in TRAITS])
        for row_idx, (fid, iid, _) in enumerate(complete_rows):
            writer.writerow([fid, iid] + [f"{z_columns[col_idx][row_idx]:.10g}" for col_idx in range(len(TRAITS))])

    raw_path = out_dir / "human_ukb_10traits.raw.pheno"
    with raw_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["FID", "IID"] + [name for _, name in TRAITS])
        for fid, iid, vals in complete_rows:
            writer.writerow([fid, iid] + [f"{v:.10g}" for v in vals])

    keep_path = out_dir / "human_ukb_complete_10traits.keep"
    with keep_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for fid, iid, _ in complete_rows:
            writer.writerow([fid, iid])

    pmerge_path = out_dir / "human_chr1_22.pmerge_list"
    with pmerge_path.open("w") as handle:
        for chrom in range(1, 23):
            handle.write(str(geno_dir / f"CHR{chrom}_final") + "\n")

    counts_path = out_dir / "human_input_audit.tsv"
    with counts_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["item", "value"])
        writer.writerow(["samples_in_psam", len(samples)])
        writer.writerow(["complete_10trait_samples", len(complete_rows)])
        writer.writerow(["dropped_for_any_missing_trait", len(samples) - len(complete_rows)])
        for name, n, mean, sd in stats:
            writer.writerow([f"{name}_n", n])
            writer.writerow([f"{name}_raw_mean", f"{mean:.10g}"])
            writer.writerow([f"{name}_raw_sd", f"{sd:.10g}"])

    print(f"wrote {pheno_path}")
    print(f"complete_10trait_samples={len(complete_rows)}")
    print(f"pmerge_list={pmerge_path}")


if __name__ == "__main__":
    main()
