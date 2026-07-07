#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


TARGET_TRAITS = [
    "TA", "EDL", "bw1", "PPIweight", "Tail length", "tibia",
    "Testis weight", "gastroc", "Plantaris", "Soleus",
]


def norm_name(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def read_fam(path):
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["FID", "IID", "father", "mother", "sex", "pheno"], dtype=str)


def read_bim(path):
    return pd.read_csv(path, sep=r"\s+", header=None,
                       names=["CHR", "SNP", "CM", "BP", "A1", "A2"], dtype=str)


def plink_file(prefix, ext):
    prefix = Path(prefix)
    return prefix.parent / f"{prefix.name}.{ext}"


def find_id_column(pheno, fam):
    best = (None, None, -1)
    fid = set(fam["FID"].astype(str))
    iid = set(fam["IID"].astype(str))
    for col in pheno.columns:
        vals = set(pheno[col].dropna().astype(str))
        candidates = [("FID", len(vals & fid)), ("IID", len(vals & iid))]
        for kind, n in candidates:
            if n > best[2]:
                best = (col, kind, n)
    return best


def match_traits(columns):
    normalized = {norm_name(c): c for c in columns}
    aliases = {
        "Tail length": ["taillength", "tail_length", "tail length"],
        "Testis weight": ["testisweight", "testis_weight", "testis weight"],
        "Plantaris": ["plantaris"],
        "Soleus": ["soleus"],
    }
    rows = []
    for target in TARGET_TRAITS:
        keys = [norm_name(target)] + [norm_name(x) for x in aliases.get(target, [])]
        match = None
        score = 0
        for key in keys:
            if key in normalized:
                match = normalized[key]
                score = 100
                break
        if match is None:
            target_norm = norm_name(target)
            for key, col in normalized.items():
                common = len(set(target_norm) & set(key))
                s = common / max(1, max(len(set(target_norm)), len(set(key))))
                if s > score:
                    match, score = col, int(100 * s)
        rows.append({"target_trait": target, "matched_column": match or "", "score": score})
    return rows


def zscore(s):
    x = pd.to_numeric(s, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return x * np.nan
    return (x - mu) / sd


def write_markdown(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_tool_status(path):
    if path and Path(path).exists():
        return json.loads(Path(path).read_text(encoding="utf-8"))
    return {}


def decode_bed_to_emmv(bfile, fam_order, out_csv, max_snps=0):
    prefix = Path(bfile)
    fam = read_fam(plink_file(prefix, "fam"))
    bim = read_bim(plink_file(prefix, "bim"))
    keep_iids = list(fam_order["IID"].astype(str))
    index = {iid: i for i, iid in enumerate(fam["IID"].astype(str))}
    keep_idx = [index[iid] for iid in keep_iids]
    n = len(fam)
    bytes_per_snp = (n + 3) // 4
    limit = len(bim) if max_snps <= 0 else min(max_snps, len(bim))
    code_to_dosage = {0: 2.0, 2: 1.0, 3: 0.0, 1: np.nan}
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with plink_file(prefix, "bed").open("rb") as bed, out_csv.open("w", newline="", encoding="utf-8") as out:
        magic = bed.read(3)
        if magic != b"\x6c\x1b\x01":
            raise RuntimeError(f"Unsupported BED magic for {plink_file(prefix, 'bed')}")
        writer = csv.writer(out)
        writer.writerow(["SNP"] + keep_iids)
        for j in range(limit):
            raw = bed.read(bytes_per_snp)
            vals = np.empty(n, dtype=float)
            pos = 0
            for byte in raw:
                for shift in (0, 2, 4, 6):
                    if pos < n:
                        vals[pos] = code_to_dosage[(byte >> shift) & 3]
                    pos += 1
            row = vals[keep_idx]
            mean = np.nanmean(row)
            if not np.isfinite(mean):
                mean = 0.0
            row = np.where(np.isfinite(row), row, mean) - 1.0
            writer.writerow([bim.iloc[j]["SNP"]] + [f"{v:.8g}" for v in row])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--qc-prefix", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--tool-status-json")
    ap.add_argument("--audit-only", action="store_true")
    ap.add_argument("--export-emmv-z", action="store_true")
    ap.add_argument("--max-emmv-snps", type=int, default=0)
    args = ap.parse_args()

    data = Path(args.data_dir)
    results = Path(args.results_dir)
    processed = data / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    pheno = pd.read_csv(data / "pheno.csv", na_values=["NA", "", "nan"])
    fam_raw = read_fam(data / "cfw.fam")
    bim_raw = read_bim(data / "cfw.bim")
    id_col, id_kind, match_n = find_id_column(pheno, fam_raw)
    status = load_tool_status(args.tool_status_json)

    audit_lines = [
        "# CFW Mouse Data Audit", "",
        f"Repository: `{Path.cwd()}`",
        f"Data directory: `{data}`",
        f"Raw FAM samples: {len(fam_raw)}",
        f"Raw BIM SNPs: {len(bim_raw)}",
        f"pheno.csv rows: {len(pheno)}",
        f"pheno.csv columns: {', '.join(pheno.columns)}",
        f"Best phenotype ID match: `{id_col}` vs FAM `{id_kind}` ({match_n} samples)", "",
        "## Tool Status", "",
    ]
    for k, v in status.items():
        audit_lines.append(f"- {k}: {v}")
    write_markdown(results / "data_audit.md", audit_lines)

    if args.audit_only:
        return

    qc_prefix = Path(args.qc_prefix)
    fam_qc = read_fam(plink_file(qc_prefix, "fam"))
    bim_qc = read_bim(plink_file(qc_prefix, "bim"))
    key = fam_qc[id_kind].astype(str)
    pheno_idx = pheno.copy()
    pheno_idx[id_col] = pheno_idx[id_col].astype(str)
    pheno_idx = pheno_idx.drop_duplicates(id_col).set_index(id_col)
    matched = pheno_idx.reindex(key)
    trait_matches = match_traits(pheno.columns)
    (processed / "trait_matches.tsv").write_text(
        pd.DataFrame(trait_matches).to_csv(sep="\t", index=False), encoding="utf-8")

    trait_cols = [r["matched_column"] for r in trait_matches]
    if any(not c for c in trait_cols):
        missing = [r["target_trait"] for r in trait_matches if not r["matched_column"]]
        raise RuntimeError(f"Could not match target traits: {missing}")

    out = pd.DataFrame({"FID": fam_qc["FID"].astype(str), "IID": fam_qc["IID"].astype(str)})
    raw_out = out.copy()
    for row in trait_matches:
        target = row["target_trait"]
        col = row["matched_column"]
        out[target] = zscore(matched[col]).to_numpy()
        raw_out[target] = pd.to_numeric(matched[col], errors="coerce").to_numpy()
    out.to_csv(processed / "cfw_mouse_traits.tsv", sep="\t", index=False, na_rep="NA")
    raw_out.to_csv(processed / "cfw_mouse_traits.raw.tsv", sep="\t", index=False, na_rep="NA")

    trait_dir = processed / "trait_sets"
    single_dir = processed / "single_trait"
    trait_dir.mkdir(exist_ok=True)
    single_dir.mkdir(exist_ok=True)
    for t in TARGET_TRAITS:
        safe = norm_name(t)
        out[["FID", "IID", t]].to_csv(single_dir / f"{safe}.pheno.tsv", sep="\t", index=False, na_rep="NA")
    for k in (2, 4, 6, 8, 10):
        cols = TARGET_TRAITS[:k]
        with (trait_dir / f"traits_{k}.txt").open("w", encoding="utf-8") as f:
            f.write("\n".join(cols) + "\n")
        out[["FID", "IID"] + cols].to_csv(trait_dir / f"traits_{k}.pheno.tsv", sep="\t", index=False, na_rep="NA")
        out[cols].to_csv(trait_dir / f"traits_{k}.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
        out[cols].to_csv(trait_dir / f"traits_{k}.emmv.P.csv", index=False, na_rep="NA")

    with (processed / "first_1000_snps.txt").open("w", encoding="utf-8") as f:
        for snp in bim_qc["SNP"].head(1000):
            f.write(str(snp) + "\n")

    if args.export_emmv_z:
        decode_bed_to_emmv(qc_prefix, fam_qc, processed / "cfw.qc.emmv.Z.csv", args.max_emmv_snps)

    report = [
        "# CFW Mouse Preprocess Report", "",
        f"QC samples: {len(fam_qc)}",
        f"QC SNPs: {len(bim_qc)}", "",
        "## Trait Matching", "",
        "| target_trait | matched_column | score |",
        "|---|---|---|",
    ]
    for r in trait_matches:
        report.append(f"| {r['target_trait']} | {r['matched_column']} | {r['score']} |")
    report += [
        "", "## Outputs", "",
        f"- Standardized traits: `{processed / 'cfw_mouse_traits.tsv'}`",
        f"- Trait sets: `{trait_dir}`",
        f"- Single-trait phenotypes: `{single_dir}`",
    ]
    write_markdown(processed / "preprocess_report.md", report)


if __name__ == "__main__":
    main()
