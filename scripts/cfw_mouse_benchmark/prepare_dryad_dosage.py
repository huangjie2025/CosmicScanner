#!/usr/bin/env python3
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


TRAITS = [
    ("TA", "TA"),
    ("EDL", "EDL"),
    ("bw1", "bw1"),
    ("PPIweight", "PPIweight"),
    ("Tail length", "taillength"),
    ("tibia", "tibia"),
    ("Testis weight", "testisweight"),
    ("gastroc", "gastroc"),
    ("Plantaris", "plantaris"),
    ("Soleus", "soleus"),
]


def zscore(x):
    x = pd.to_numeric(x, errors="coerce")
    return (x - x.mean()) / x.std(ddof=1)


def read_geno_ids(path):
    ids = []
    discard = []
    with gzip.open(path, "rt") as f:
        header = f.readline().strip().split()
        snps = header[2:]
        for line in f:
            parts = line.split(maxsplit=2)
            ids.append(parts[0])
            discard.append(parts[1])
    return snps, pd.DataFrame({"id": ids, "geno_discard": discard})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--set-name", default="dryad_dosage")
    args = ap.parse_args()

    data = Path(args.data_dir)
    out = Path(args.out_dir) / args.set_name
    out.mkdir(parents=True, exist_ok=True)

    geno_path = data / "geno.txt.gz"
    map_path = data / "map.txt"
    pheno_path = data / "pheno.csv"
    if not geno_path.exists() or not map_path.exists():
        raise FileNotFoundError("Expected geno.txt.gz and map.txt in data dir")

    snps, geno_ids = read_geno_ids(geno_path)
    pheno = pd.read_csv(pheno_path, na_values=["NA", "", "nan"])
    pheno["id"] = pheno["id"].astype(str)
    trait_cols = [src for _, src in TRAITS]

    merged = geno_ids.merge(pheno[["id", "discard", "mixup"] + trait_cols], on="id", how="left")
    trait_complete = merged[trait_cols].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    # This reproduces the manuscript filtered sample count: pheno discard=no and all 10 traits observed.
    keep_mask = merged["discard"].eq("no") & trait_complete
    keep_ids = merged.loc[keep_mask, "id"].astype(str).tolist()
    keep_indices = np.where(keep_mask.to_numpy())[0]

    map_df = pd.read_csv(map_path, sep=r"\s+")
    if list(map_df["id"].astype(str)) != snps:
        raise RuntimeError("map.txt SNP order does not match geno.txt.gz header")

    n_samples = len(keep_ids)
    n_snps = len(snps)
    G = np.empty((n_samples, n_snps), dtype=np.float32)
    with gzip.open(geno_path, "rt") as f:
        f.readline()
        out_i = 0
        for row_i, line in enumerate(f):
            if row_i not in set(keep_indices):
                continue
            parts = line.strip().split()
            vals = np.array(parts[2:], dtype=np.float32)
            G[out_i, :] = vals
            out_i += 1
    if out_i != n_samples:
        raise RuntimeError(f"Loaded {out_i} samples, expected {n_samples}")

    ph = pheno.set_index("id").loc[keep_ids].reset_index()
    samples = pd.DataFrame({"FID": keep_ids, "IID": keep_ids})
    samples.to_csv(out / "samples.txt", sep="\t", header=False, index=False)
    with (out / "plink.psam").open("w", encoding="utf-8") as f:
        f.write("#FID\tIID\tSEX\n")
        for sid in keep_ids:
            f.write(f"{sid}\t{sid}\t1\n")

    pheno_out = samples.copy()
    raw_pheno_out = samples.copy()
    for target, src in TRAITS:
        raw_pheno_out[target] = pd.to_numeric(ph[src], errors="coerce").to_numpy()
        pheno_out[target] = zscore(ph[src]).to_numpy()
    pheno_out.to_csv(out / "pheno.tsv", sep="\t", index=False, na_rep="NA")
    raw_pheno_out.to_csv(out / "pheno.raw.tsv", sep="\t", index=False, na_rep="NA")
    for k in (2, 4, 6, 8, 10):
        cols = [t for t, _ in TRAITS[:k]]
        tdir = out / f"traits_{k}"
        tdir.mkdir(exist_ok=True)
        pheno_out[["FID", "IID"] + cols].to_csv(tdir / "pheno.tsv", sep="\t", index=False, na_rep="NA")
        pheno_out[cols].to_csv(tdir / "gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
        pheno_out[cols].to_csv(tdir / "emmv.P.csv", index=False, na_rep="NA")

    map_df.rename(columns={"id": "SNP", "chr": "CHR", "pos": "BP", "ref": "REF", "alt": "ALT"}).to_csv(out / "snp_info.tsv", sep="\t", index=False)
    with (out / "plink.map").open("w", encoding="utf-8") as f:
        for r in map_df.itertuples(index=False):
            f.write(f"{r.chr}\t{r.id}\t0\t{r.pos}\n")

    means = G.mean(axis=0)
    maf = np.minimum(means / 2.0, 1.0 - means / 2.0)
    pd.DataFrame({"SNP": snps, "CHR": map_df["chr"], "BP": map_df["pos"], "quality": map_df["quality"], "maf": maf, "mean_dosage": means}).to_csv(out / "snp_qc.tsv", sep="\t", index=False)

    # EMmvGWAS-compatible marker matrix: rows SNPs, columns samples, dosage centered at 1.
    with (out / "emmv.Z.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["SNP"] + keep_ids)
        for sid, row in zip(snps, (G - 1.0).T):
            writer.writerow([sid] + [f"{x:.8g}" for x in row])

    # PLINK2 dosage import format: variant rows, sample dosage columns.
    plink_sample_ids = [f"{sid}_{sid}" for sid in keep_ids]
    with (out / "plink.dosage.tsv").open("w", encoding="utf-8") as f:
        f.write("SNP\tA1\tA2\t" + "\t".join(plink_sample_ids) + "\n")
        for sid, ref, alt, row in zip(snps, map_df["ref"], map_df["alt"], G.T):
            # Dryad dosage is expected alt allele count, so A1 is ALT for PLINK2 import.
            f.write(f"{sid}\t{alt}\t{ref}\t" + "\t".join(f"{x:.8g}" for x in row) + "\n")

    audit = {
        "pheno_rows": int(len(pheno)),
        "geno_rows": int(len(geno_ids)),
        "geno_discard_no": int(geno_ids["geno_discard"].eq("no").sum()),
        "pheno_discard_no_in_geno": int(merged["discard"].eq("no").sum()),
        "complete_10_traits_all_geno": int(trait_complete.sum()),
        "filtered_samples_pheno_discard_no_complete10": int(n_samples),
        "snps": int(n_snps),
    }
    (out / "sample_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    manifest = {"set_name": args.set_name, "samples": n_samples, "snps": n_snps, "traits": [t for t, _ in TRAITS], "input": "Dryad geno.txt.gz dosage + map.txt", "sample_filter": "pheno discard=no and complete 10 target traits"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    report = ["# Dryad Dosage CFW Input Report", "", f"Samples: {n_samples}", f"SNPs: {n_snps}", "", "## Audit", "", "```json", json.dumps(audit, indent=2), "```"]
    (out / "input_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
