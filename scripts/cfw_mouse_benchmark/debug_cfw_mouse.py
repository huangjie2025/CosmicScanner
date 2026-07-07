#!/usr/bin/env python3
import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

TARGET_TRAITS = ["TA", "EDL", "bw1", "PPIweight", "Tail length", "tibia", "Testis weight", "gastroc", "Plantaris", "Soleus"]


def norm_name(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def plink_file(prefix, ext):
    p = Path(prefix)
    return p.parent / f"{p.name}.{ext}"


def read_fam(path):
    return pd.read_csv(path, sep=r"\s+", header=None, names=["FID", "IID", "father", "mother", "sex", "pheno"], dtype=str)


def read_bim(path):
    return pd.read_csv(path, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"], dtype=str)


def match_traits(columns):
    aliases = {
        "Tail length": ["taillength", "tail_length"],
        "Testis weight": ["testisweight", "testis_weight"],
        "Plantaris": ["plantaris"],
        "Soleus": ["soleus"],
    }
    normalized = {norm_name(c): c for c in columns}
    rows = []
    for target in TARGET_TRAITS:
        keys = [norm_name(target)] + [norm_name(a) for a in aliases.get(target, [])]
        col = ""
        for key in keys:
            if key in normalized:
                col = normalized[key]
                break
        rows.append({"target_trait": target, "matched_column": col})
    return rows


def zscore(x):
    x = pd.to_numeric(x, errors="coerce")
    return (x - x.mean(skipna=True)) / x.std(skipna=True, ddof=1)


def write_md(path, lines):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(str(x) for x in lines) + "\n", encoding="utf-8")


def quantile_rows(name, s):
    x = pd.to_numeric(s, errors="coerce").dropna()
    qs = x.quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
    return {
        "metric": name, "n": len(x), "min": x.min(), "q01": qs.loc[0.01], "q05": qs.loc[0.05],
        "q25": qs.loc[0.25], "median": qs.loc[0.5], "mean": x.mean(), "q75": qs.loc[0.75],
        "q95": qs.loc[0.95], "q99": qs.loc[0.99], "max": x.max(),
    }


def load_assoc(path):
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    sep = "," if p.suffix == ".csv" else "\t"
    df = pd.read_csv(p, sep=sep)
    rename = {}
    p_selected = False
    preferred_snp = "SNP_ID" if "SNP_ID" in df.columns else ("SNP" if "SNP" in df.columns else ("rs" if "rs" in df.columns else None))
    for c in df.columns:
        lc = c.lower()
        if c == preferred_snp:
            rename[c] = "SNP"
        elif lc in ("p", "p_wald") and not p_selected:
            rename[c] = "P"
            p_selected = True
        elif lc in ("beta", "beta_1"):
            rename[c] = "BETA"
        elif lc in ("se", "se_1"):
            rename[c] = "SE"
        elif lc in ("wald", "wald_stat", "chisq", "chi2"):
            rename[c] = "STAT"
    df = df.rename(columns=rename)
    return df


def schema(path, out):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"), nrows=5)
    rows = []
    for i, c in enumerate(df.columns, 1):
        vals = df[c].dropna().astype(str).head(3).tolist()
        rows.append({"index": i, "column": c, "example_values": ";".join(vals)})
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)


def corr_rows(joined, pairs):
    out = []
    for label, a, b, transform in pairs:
        x = pd.to_numeric(joined[a], errors="coerce")
        y = pd.to_numeric(joined[b], errors="coerce")
        if transform == "logp":
            x = -np.log10(x.clip(lower=1e-300))
            y = -np.log10(y.clip(lower=1e-300))
        ok = np.isfinite(x) & np.isfinite(y)
        pearson = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
        if ok.sum() > 2:
            rx = pd.Series(x[ok]).rank(method="average").to_numpy()
            ry = pd.Series(y[ok]).rank(method="average").to_numpy()
            spearman = float(np.corrcoef(rx, ry)[0, 1])
        else:
            spearman = np.nan
        out.append({"metric": label, "n": int(ok.sum()), "pearson": pearson, "spearman": spearman})
    return out


def plot_scatter(joined, xcol, ycol, out, xlabel, ylabel, transform_logp=True):
    x = pd.to_numeric(joined[xcol], errors="coerce")
    y = pd.to_numeric(joined[ycol], errors="coerce")
    if transform_logp:
        x = -np.log10(x.clip(lower=1e-300))
        y = -np.log10(y.clip(lower=1e-300))
    ok = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(4.5, 4.2))
    plt.scatter(x[ok], y[ok], s=8, alpha=0.45)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"Pearson r={r:.4f}, n={ok.sum()}")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def decode_bed_to_emmv(bfile, keep_file, out_csv, max_snps=0):
    prefix = Path(bfile)
    fam = read_fam(plink_file(prefix, "fam"))
    bim = read_bim(plink_file(prefix, "bim"))
    keep = pd.read_csv(keep_file, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    index = {iid: i for i, iid in enumerate(fam["IID"].astype(str))}
    keep_iids = [iid for iid in keep["IID"].astype(str) if iid in index]
    keep_idx = [index[iid] for iid in keep_iids]
    n = len(fam)
    bytes_per_snp = (n + 3) // 4
    limit = len(bim) if max_snps <= 0 else min(max_snps, len(bim))
    code_to_dosage = {0: 2.0, 2: 1.0, 3: 0.0, 1: np.nan}
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with plink_file(prefix, "bed").open("rb") as bed, Path(out_csv).open("w", newline="", encoding="utf-8") as out:
        if bed.read(3) != b"\x6c\x1b\x01":
            raise RuntimeError("Unsupported BED magic")
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


def cmd_missing(args):
    res = Path(args.results_dir)
    raw = Path(args.raw_prefix)
    smiss = pd.read_csv(Path(args.missing_prefix).with_suffix(".smiss"), sep=r"\s+")
    vmiss = pd.read_csv(Path(args.missing_prefix).with_suffix(".vmiss"), sep=r"\s+")
    s_col = "F_MISS" if "F_MISS" in smiss.columns else smiss.columns[-1]
    v_col = "F_MISS" if "F_MISS" in vmiss.columns else vmiss.columns[-1]
    rows = [quantile_rows("per_sample_missing_rate", smiss[s_col]), quantile_rows("per_snp_missing_rate", vmiss[v_col])]
    out = res / "tables/missingness_summary.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    hi = smiss.sort_values(s_col, ascending=False).head(20)
    lo = smiss.sort_values(s_col, ascending=True).head(20)
    hi.to_csv(res / "tables/missingness_samples_highest20.tsv", sep="\t", index=False)
    lo.to_csv(res / "tables/missingness_samples_lowest20.tsv", sep="\t", index=False)
    readme = ""
    for p in [raw.parent / "README_for_pheno.txt", raw.parent / "README_for_geno.txt", raw.parent / "README_for_geno.txt.txt", raw.parent / "README_for_cfw.txt"]:
        if p.exists():
            readme += p.read_text(encoding="utf-8", errors="replace") + "\n"
    note = "README says genotypes are dosages 0-2 and can be real-valued when genotype uncertainty exists; PLINK BED stores hard calls/missing, so missingness filters on BED can be inappropriate for the paper-like imputed-dosage analysis."
    table = ["| metric | n | min | median | mean | max |", "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        table.append(f"| {r['metric']} | {r['n']} | {r['min']:.6g} | {r['median']:.6g} | {r['mean']:.6g} | {r['max']:.6g} |")
    write_md(res / "debug/missingness_debug.md", ["# Missingness Debug", "", f"Raw prefix: `{raw}`", "", "## Summary", "", *table, "", "## Interpretation", "", note])


def cmd_profiles(args):
    data = Path(args.data_dir)
    res = Path(args.results_dir)
    proc = data / "processed"
    paper = proc / "paper_like"
    strict = proc / "strict_qc"
    paper.mkdir(parents=True, exist_ok=True)
    strict.mkdir(parents=True, exist_ok=True)
    pheno = pd.read_csv(data / "pheno.csv", na_values=["NA", "", "nan"])
    fam = read_fam(data / "cfw.fam")
    bim = read_bim(data / "cfw.bim")
    matches = match_traits(pheno.columns)
    id_col = "id"
    pheno[id_col] = pheno[id_col].astype(str)
    ph = pheno.drop_duplicates(id_col).set_index(id_col)
    aligned = ph.reindex(fam["IID"].astype(str))
    out = pd.DataFrame({"FID": fam["FID"].astype(str), "IID": fam["IID"].astype(str)})
    raw = out.copy()
    trait_rows = []
    complete_mask = np.ones(len(out), dtype=bool)
    for m in matches:
        target, col = m["target_trait"], m["matched_column"]
        x_all = pd.to_numeric(pheno[col], errors="coerce")
        x = pd.to_numeric(aligned[col], errors="coerce")
        z = zscore(x)
        out[target] = z.to_numpy()
        raw[target] = x.to_numpy()
        complete_mask &= x.notna().to_numpy()
        trait_rows.append({
            "target_trait": target, "matched_column": col, "raw_nonmissing": int(x_all.notna().sum()),
            "fam_matched_nonmissing": int(x.notna().sum()), "raw_mean": x.mean(), "raw_sd": x.std(ddof=1),
            "raw_min": x.min(), "raw_max": x.max(), "z_mean": z.mean(), "z_sd": z.std(ddof=1),
        })
    out.to_csv(paper / "cfw_mouse_traits.tsv", sep="\t", index=False, na_rep="NA")
    raw.to_csv(paper / "cfw_mouse_traits.raw.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(trait_rows).to_csv(res / "tables/trait_matching_debug.tsv", sep="\t", index=False)
    complete = out.loc[complete_mask, ["FID", "IID"]]
    complete.to_csv(paper / "mouse_10traits_complete_samples.txt", sep="\t", index=False, header=False)
    out.loc[complete_mask].to_csv(paper / "traits_10.complete.pheno.tsv", sep="\t", index=False, na_rep="NA")
    out.loc[complete_mask, TARGET_TRAITS].to_csv(paper / "traits_10.complete.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
    out.loc[complete_mask, ["FID", "IID", "TA"]].to_csv(paper / "TA.complete.pheno.tsv", sep="\t", index=False, na_rep="NA")
    out.loc[complete_mask, ["TA"]].to_csv(paper / "TA.complete.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
    out.loc[complete_mask, ["FID", "IID", "TA", "EDL"]].to_csv(paper / "traits_2.complete.pheno.tsv", sep="\t", index=False, na_rep="NA")
    out.loc[complete_mask, ["TA", "EDL"]].to_csv(paper / "traits_2.complete.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
    out.loc[complete_mask, ["TA", "EDL"]].to_csv(paper / "traits_2.complete.emmv.P.csv", index=False, na_rep="NA")
    out.loc[complete_mask, TARGET_TRAITS].to_csv(paper / "traits_10.complete.emmv.P.csv", index=False, na_rep="NA")
    with (paper / "first_1000_raw_snps.txt").open("w", encoding="utf-8") as f:
        f.write("\n".join(bim["SNP"].head(1000).astype(str)) + "\n")
    lines = ["# Trait Matching Debug", "", "## pheno.csv Columns", "", ", ".join(pheno.columns), "", "## Matched Traits", "", "| target | column | raw_nonmissing | fam_matched_nonmissing | mean | sd | min | max | z_mean | z_sd |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in trait_rows:
        lines.append(f"| {r['target_trait']} | {r['matched_column']} | {r['raw_nonmissing']} | {r['fam_matched_nonmissing']} | {r['raw_mean']:.6g} | {r['raw_sd']:.6g} | {r['raw_min']:.6g} | {r['raw_max']:.6g} | {r['z_mean']:.3g} | {r['z_sd']:.6g} |")
    lines += ["", f"10-trait complete-case samples: {len(complete)}", "", "Article filtered sample size is ~1091. A small difference is expected if the paper used additional phenotype/sample flags or a slightly different trait missingness rule."]
    write_md(res / "trait_matching_debug.md", lines)
    prof = ["# CFW Mouse Preprocessing Profiles", "", "## paper_like", "", f"- Directory: `{paper}`", "- Uses raw PLINK BED/BIM/FAM without `--geno 0.05` or `--mind 0.05`.", "- SNP target: preserve the raw 92734 SNP IDs for paper-like comparisons; debug scans use a 1000-SNP extract only.", "- Samples are controlled by FAM/IID phenotype matching and 10-trait complete-case filtering.", f"- 10-trait complete-case samples: {len(complete)}", "", "## strict_qc", "", f"- Directory: `{strict}`", "- Represents robustness profile derived from the earlier fallback QC (`--maf 0.01 --geno 0.05`, no `--mind`).", "- The requested strict `--mind 0.05` is not paper-like for this BED because it removes all samples.", ""]
    write_md(res / "preprocessing_profiles.md", prof)


def cmd_join_univariate(args):
    gemma = load_assoc(args.gemma)
    cosmic = load_assoc(args.cosmic)
    if "SNP" not in gemma.columns and "rs" in gemma.columns:
        gemma = gemma.rename(columns={"rs": "SNP"})
    gemma = gemma.rename(columns={"P": "P_gemma", "BETA": "BETA_gemma", "STAT": "STAT_gemma", "beta": "BETA_gemma", "p_wald": "P_gemma"})
    cosmic = cosmic.rename(columns={"P": "P_cosmic", "BETA": "BETA_cosmic", "STAT": "STAT_cosmic"})
    keep_g = [c for c in ["SNP", "P_gemma", "BETA_gemma", "STAT_gemma"] if c in gemma.columns]
    keep_c = [c for c in ["SNP", "P_cosmic", "BETA_cosmic", "STAT_cosmic"] if c in cosmic.columns]
    joined = gemma[keep_g].merge(cosmic[keep_c], on="SNP", how="inner")
    if "STAT_gemma" not in joined and "BETA_gemma" in joined and "P_gemma" in joined:
        joined["STAT_gemma"] = np.nan
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.out, sep="\t", index=False)
    pairs = [("minus_log10_p", "P_gemma", "P_cosmic", "logp")]
    if "BETA_gemma" in joined and "BETA_cosmic" in joined:
        pairs.append(("beta", "BETA_gemma", "BETA_cosmic", "none"))
    rows = corr_rows(joined, pairs)
    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
    plot_scatter(joined, "P_gemma", "P_cosmic", args.plot, "GEMMA -log10(p)", "CosmicScanner -log10(p)")


def cmd_join_multi(args):
    gemma = load_assoc(args.gemma).rename(columns={"P": "P_gemma", "p_wald": "P_gemma"})
    emmv = load_assoc(args.emmv).rename(columns={"P": "P_emmv", "p": "P_emmv"})
    cosmic = load_assoc(args.cosmic).rename(columns={"P": "P_cosmic"})
    if "rs" in gemma.columns:
        gemma = gemma.rename(columns={"rs": "SNP"})
    if "SNP_ID" in emmv.columns:
        emmv = emmv.rename(columns={"SNP_ID": "SNP"})
    joined = gemma[["SNP", "P_gemma"]].merge(emmv[["SNP", "P_emmv"]], on="SNP", how="inner").merge(cosmic[["SNP", "P_cosmic"]], on="SNP", how="inner")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.out, sep="\t", index=False)
    rows = corr_rows(joined, [("EMmvGWAS_vs_GEMMA_logp", "P_gemma", "P_emmv", "logp"), ("CosmicScanner_vs_GEMMA_logp", "P_gemma", "P_cosmic", "logp"), ("CosmicScanner_vs_EMmvGWAS_logp", "P_emmv", "P_cosmic", "logp")])
    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
    plot_scatter(joined, "P_gemma", "P_emmv", args.plot_emmv, "GEMMA -log10(p)", "EMmvGWAS -log10(p)")
    plot_scatter(joined, "P_gemma", "P_cosmic", args.plot_cosmic, "GEMMA -log10(p)", "CosmicScanner -log10(p)")


def cmd_alignment(args):
    res = Path(args.results_dir)
    lines = ["# Alignment Debug", ""]
    for label, prefix in [("raw", args.raw_prefix), ("paper_like_debug", args.debug_prefix)]:
        fam = read_fam(plink_file(prefix, "fam"))
        bim = read_bim(plink_file(prefix, "bim"))
        lines += [f"## {label}", "", f"Samples: {len(fam)}", f"First 10 samples: {', '.join(fam['IID'].head(10))}", f"Last 10 samples: {', '.join(fam['IID'].tail(10))}", f"SNPs: {len(bim)}", f"First 10 SNPs: {', '.join(bim['SNP'].head(10))}", f"Last 10 SNPs: {', '.join(bim['SNP'].tail(10))}", ""]
    for label, path in [("GEMMA", args.gemma), ("EMmvGWAS", args.emmv), ("CosmicScanner", args.cosmic)]:
        df = load_assoc(path)
        snps = df["SNP"].astype(str) if "SNP" in df else pd.Series([], dtype=str)
        lines += [f"## {label} Output", "", f"Path: `{path}`", f"Rows: {len(df)}", f"Has SNP ID: {'SNP' in df.columns}", f"First 10 SNPs: {', '.join(snps.head(10))}", f"Last 10 SNPs: {', '.join(snps.tail(10))}", ""]
    write_md(res / "alignment_debug.md", lines)


def cmd_schema(args):
    schema(args.gemma, Path(args.results_dir) / "tables/result_schema_gemma.tsv")
    schema(args.emmv, Path(args.results_dir) / "tables/result_schema_emmvgwas.tsv")
    schema(args.cosmic, Path(args.results_dir) / "tables/result_schema_cosmicscanner.tsv")


def cmd_decode(args):
    decode_bed_to_emmv(args.bfile, args.keep, args.out, args.max_snps)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("missing"); p.add_argument("--raw-prefix", required=True); p.add_argument("--missing-prefix", required=True); p.add_argument("--results-dir", required=True); p.set_defaults(func=cmd_missing)
    p = sub.add_parser("profiles"); p.add_argument("--data-dir", required=True); p.add_argument("--results-dir", required=True); p.set_defaults(func=cmd_profiles)
    p = sub.add_parser("join-univariate"); p.add_argument("--gemma", required=True); p.add_argument("--cosmic", required=True); p.add_argument("--out", required=True); p.add_argument("--summary", required=True); p.add_argument("--plot", required=True); p.set_defaults(func=cmd_join_univariate)
    p = sub.add_parser("join-multi"); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic", required=True); p.add_argument("--out", required=True); p.add_argument("--summary", required=True); p.add_argument("--plot-emmv", required=True); p.add_argument("--plot-cosmic", required=True); p.set_defaults(func=cmd_join_multi)
    p = sub.add_parser("alignment"); p.add_argument("--results-dir", required=True); p.add_argument("--raw-prefix", required=True); p.add_argument("--debug-prefix", required=True); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic", required=True); p.set_defaults(func=cmd_alignment)
    p = sub.add_parser("schema"); p.add_argument("--results-dir", required=True); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic", required=True); p.set_defaults(func=cmd_schema)
    p = sub.add_parser("decode-emmv-z"); p.add_argument("--bfile", required=True); p.add_argument("--keep", required=True); p.add_argument("--out", required=True); p.add_argument("--max-snps", type=int, default=0); p.set_defaults(func=cmd_decode)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
