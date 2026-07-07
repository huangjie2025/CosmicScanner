#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PHENO_COLS = ["TA", "EDL", "bw1", "PPIweight", "Tail length", "tibia", "Testis weight", "gastroc", "Plantaris", "Soleus"]
TRAIT_SETS = {2: PHENO_COLS[:2], 4: PHENO_COLS[:4], 6: PHENO_COLS[:6], 8: PHENO_COLS[:8], 10: PHENO_COLS[:10]}


def plink_file(prefix, ext):
    p = Path(prefix)
    return p.parent / f"{p.name}.{ext}"


def read_fam(path):
    return pd.read_csv(path, sep=r"\s+", header=None, names=["FID", "IID", "father", "mother", "sex", "pheno"], dtype=str)


def read_bim(path):
    return pd.read_csv(path, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"], dtype=str)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def decode_bed(prefix, sample_iids, snp_ids):
    fam = read_fam(plink_file(prefix, "fam"))
    bim = read_bim(plink_file(prefix, "bim"))
    sample_index = {iid: i for i, iid in enumerate(fam["IID"].astype(str))}
    snp_index = {sid: i for i, sid in enumerate(bim["SNP"].astype(str))}
    keep_sample_idx = [sample_index[iid] for iid in sample_iids]
    keep_snp_idx = [snp_index[sid] for sid in snp_ids]
    n = len(fam)
    bytes_per_snp = (n + 3) // 4
    code_to_dosage = {0: 2.0, 2: 1.0, 3: 0.0, 1: np.nan}
    G = np.empty((len(sample_iids), len(snp_ids)), dtype=float)
    with plink_file(prefix, "bed").open("rb") as bed:
        if bed.read(3) != b"\x6c\x1b\x01":
            raise RuntimeError("Unsupported BED magic")
        for out_j, snp_j in enumerate(keep_snp_idx):
            bed.seek(3 + snp_j * bytes_per_snp)
            raw = bed.read(bytes_per_snp)
            vals = np.empty(n, dtype=float)
            pos = 0
            for byte in raw:
                for shift in (0, 2, 4, 6):
                    if pos < n:
                        vals[pos] = code_to_dosage[(byte >> shift) & 3]
                    pos += 1
            G[:, out_j] = vals[keep_sample_idx]
    return G, bim.set_index("SNP").loc[snp_ids].reset_index()


def write_bed(prefix, samples, snp_info, G):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with plink_file(prefix, "fam").open("w", encoding="utf-8") as f:
        for _, r in samples.iterrows():
            f.write(f"{r.FID}\t{r.IID}\t0\t0\t1\t-9\n")
    snp_info[["CHR", "SNP", "CM", "BP", "A1", "A2"]].to_csv(plink_file(prefix, "bim"), sep="\t", header=False, index=False)
    n = G.shape[0]
    pad_n = ((n + 3) // 4) * 4
    with plink_file(prefix, "bed").open("wb") as bed:
        bed.write(b"\x6c\x1b\x01")
        for j in range(G.shape[1]):
            g = np.rint(G[:, j]).astype(int)
            codes = np.zeros(pad_n, dtype=np.uint8)
            codes[:n] = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
            packed = (codes[0::4] | (codes[1::4] << 2) | (codes[2::4] << 4) | (codes[3::4] << 6)).astype(np.uint8)
            bed.write(packed.tobytes())


def build_set(args):
    out = Path(args.out_dir) / args.set_name
    out.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(args.samples, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    sample_iids = samples["IID"].astype(str).tolist()
    vmiss = pd.read_csv(args.vmiss, sep=r"\s+")
    snp_col = "ID" if "ID" in vmiss.columns else "SNP"
    miss_col = "F_MISS" if "F_MISS" in vmiss.columns else vmiss.columns[-1]
    vmiss = vmiss.rename(columns={snp_col: "SNP", miss_col: "F_MISS"})
    pool = vmiss[vmiss["F_MISS"] <= args.max_missing_rate].sort_values(["F_MISS", "SNP"])
    # Decode a larger candidate pool so MAF is computed on the fixed paper_like samples.
    max_candidates = len(pool) if args.n_snps <= 0 else min(len(pool), max(args.n_snps * 20, args.n_snps))
    candidate_snps = pool["SNP"].astype(str).head(max_candidates).tolist()
    G0, info0 = decode_bed(args.raw_prefix, sample_iids, candidate_snps)
    means0 = np.nanmean(G0, axis=0)
    maf0 = np.minimum(means0 / 2.0, 1.0 - means0 / 2.0)
    ok = np.where(maf0 >= args.min_maf)[0]
    if args.n_snps > 0:
        ok = ok[:args.n_snps]
    snps = [candidate_snps[i] for i in ok]
    G = G0[:, ok]
    snp_info = info0.iloc[ok].reset_index(drop=True)
    missing_before = np.isnan(G)
    imputed_counts = missing_before.sum(axis=0).astype(int)
    means = np.nanmean(G, axis=0)
    G_imp = np.where(np.isnan(G), means[None, :], G)
    maf = np.minimum(np.nanmean(G_imp, axis=0) / 2.0, 1.0 - np.nanmean(G_imp, axis=0) / 2.0)
    geno = pd.DataFrame(G_imp, columns=snps)
    geno.insert(0, "IID", sample_iids)
    geno_path = out / "genotype.tsv"
    geno.to_csv(geno_path, sep="\t", index=False, float_format="%.8g")
    samples.to_csv(out / "samples.txt", sep="\t", header=False, index=False)
    pd.DataFrame({"SNP": snps}).to_csv(out / "snps.txt", header=False, index=False)
    snp_info.to_csv(out / "snp_info.tsv", sep="\t", index=False)
    pd.DataFrame({"SNP": snps, "missing_rate_raw": pool.set_index("SNP").loc[snps]["F_MISS"].to_numpy(), "imputed_count": imputed_counts, "maf": maf, "mean_dosage": means}).to_csv(out / "snp_qc.tsv", sep="\t", index=False)
    z = G_imp - 1.0
    with (out / "emmv.Z.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["SNP"] + sample_iids)
        for sid, row in zip(snps, z.T):
            writer.writerow([sid] + [f"{x:.8g}" for x in row])
    Zc = G_imp - G_imp.mean(axis=0, keepdims=True)
    denom = np.sum(2 * (G_imp.mean(axis=0) / 2.0) * (1 - G_imp.mean(axis=0) / 2.0))
    K = (Zc @ Zc.T) / denom if denom > 0 else (Zc @ Zc.T) / G_imp.shape[1]
    pd.DataFrame(K, index=sample_iids, columns=sample_iids).to_csv(out / "kinship.tsv", sep="\t")
    write_bed(out / "plink", samples, snp_info, G_imp)
    pheno = pd.read_csv(args.pheno, sep="\t")
    pheno["IID"] = pheno["IID"].astype(str)
    ph = pheno.set_index("IID").loc[sample_iids].reset_index()
    for k, cols in TRAIT_SETS.items():
        tdir = out / f"traits_{k}"
        tdir.mkdir(exist_ok=True)
        keep = ["FID", "IID"] + cols
        ph[keep].to_csv(tdir / "pheno.tsv", sep="\t", index=False, na_rep="NA")
        ph[cols].to_csv(tdir / "gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
        ph[cols].to_csv(tdir / "emmv.P.csv", index=False, na_rep="NA")
    summaries = []
    for name, arr in [("missing_rate_raw", pool.set_index("SNP").loc[snps]["F_MISS"].to_numpy()), ("maf", maf), ("imputed_count", imputed_counts)]:
        s = pd.Series(arr)
        summaries.append({"metric": name, "n": len(s), "min": s.min(), "median": s.median(), "mean": s.mean(), "max": s.max()})
    pd.DataFrame(summaries).to_csv(out / "input_summary.tsv", sep="\t", index=False)
    hashes = {p.name: sha256(p) for p in [geno_path, out / "emmv.Z.csv", out / "kinship.tsv", plink_file(out / "plink", "bed"), plink_file(out / "plink", "bim"), plink_file(out / "plink", "fam")]}
    manifest = {"set_name": args.set_name, "samples": len(sample_iids), "snps": len(snps), "requested_snps": args.n_snps, "max_missing_rate": args.max_missing_rate, "min_maf": args.min_maf, "missing_after_imputation": int(np.isnan(G_imp).sum()), "hashes": hashes}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    lines = ["# Scaled Harmonized Input Report", "", f"Set: {args.set_name}", f"Samples: {len(sample_iids)}", f"Requested SNPs: {args.n_snps if args.n_snps > 0 else 'all'}", f"Actual SNPs: {len(snps)}", f"Max raw missing rate: {args.max_missing_rate}", f"Min harmonized MAF: {args.min_maf}", f"Missing after imputation: {int(np.isnan(G_imp).sum())}", "", "## Hashes", ""]
    lines += [f"- {k}: `{v}`" for k, v in hashes.items()]
    (out / "input_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_result(path, method):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"))
    if method == "gemma":
        return df.rename(columns={"rs": "SNP", "p_wald": "P", "beta_1": "BETA_1", "beta_2": "BETA_2"})
    if method == "emmv":
        return df.rename(columns={"SNP_ID": "SNP", "p": "P", "wald": "STAT"})
    return df.rename(columns={"WALD_STAT": "STAT"})


def corr(x, y):
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        return int(ok.sum()), np.nan, np.nan
    a = x[ok].to_numpy(); b = y[ok].to_numpy()
    return int(ok.sum()), float(np.corrcoef(a, b)[0, 1]), float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def summarize_combo(args):
    out = Path(args.out_dir)
    g = load_result(args.gemma, "gemma").add_suffix("_gemma").rename(columns={"SNP_gemma": "SNP"})
    e = load_result(args.emmv, "emmv").add_suffix("_emmv").rename(columns={"SNP_emmv": "SNP"})
    cg = load_result(args.cosmic_gemma_full, "cosmic").add_suffix("_cs_gemma_full").rename(columns={"SNP_cs_gemma_full": "SNP"})
    ce = load_result(args.cosmic_emmv, "cosmic").add_suffix("_cs_emmv").rename(columns={"SNP_cs_emmv": "SNP"})
    j = g.merge(e, on="SNP", how="inner").merge(cg, on="SNP", how="inner").merge(ce, on="SNP", how="inner")
    j.to_csv(out / "joined.tsv", sep="\t", index=False)
    rows = []
    for name, a, b in [("EMmvGWAS_vs_GEMMA_logp", "P_gemma", "P_emmv"), ("CosmicScanner_gemma_full_vs_GEMMA_logp", "P_gemma", "P_cs_gemma_full"), ("CosmicScanner_emmv_compatible_vs_EMmvGWAS_logp", "P_emmv", "P_cs_emmv"), ("CosmicScanner_gemma_full_vs_EMmvGWAS_logp", "P_emmv", "P_cs_gemma_full")]:
        x = -np.log10(pd.to_numeric(j[a], errors="coerce").clip(lower=1e-300))
        y = -np.log10(pd.to_numeric(j[b], errors="coerce").clip(lower=1e-300))
        n, p, s = corr(x, y)
        rows.append({"set": args.set_name, "traits": args.traits, "comparison": name, "n": n, "pearson": p, "spearman": s})
    pd.DataFrame(rows).to_csv(out / "correlation.tsv", sep="\t", index=False)
    counts = pd.DataFrame([{"set": args.set_name, "traits": args.traits, "input_snps": args.input_snps, "gemma_rows": len(g), "emmvgwas_rows": len(e), "cosmic_gemma_full_rows": len(cg), "cosmic_emmv_rows": len(ce), "inner_join_snps": len(j)}])
    counts.to_csv(out / "counts.tsv", sep="\t", index=False)


def collect(args):
    root = Path(args.out_dir)
    rows = []
    count_rows = []
    runtime_rows = []
    for corr_path in root.glob("*/traits_*/results/correlation.tsv"):
        rows.append(pd.read_csv(corr_path, sep="\t"))
    for count_path in root.glob("*/traits_*/results/counts.tsv"):
        count_rows.append(pd.read_csv(count_path, sep="\t"))
    corr_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    counts_df = pd.concat(count_rows, ignore_index=True) if count_rows else pd.DataFrame()
    if not corr_df.empty:
        corr_df.to_csv(root / "scaled_harmonized_summary.tsv", sep="\t", index=False)
    metrics = make_metrics(corr_df, counts_df)
    if not metrics.empty:
        metrics.to_csv(root / "scaled_harmonized_metrics.tsv", sep="\t", index=False)
    for tpath in root.glob("*/traits_*/logs/*.time"):
        parts = tpath.parts
        set_name = parts[-4]
        traits = int(parts[-3].split("_")[-1])
        method = tpath.stem
        vals = {"seconds": np.nan, "max_rss_kb": np.nan}
        for line in tpath.read_text(encoding="utf-8").splitlines():
            p = line.split()
            if len(p) == 2:
                vals[p[0]] = float(p[1])
        runtime_rows.append({"set": set_name, "traits": traits, "method": method, **vals})
    rt = pd.DataFrame(runtime_rows)
    if not rt.empty:
        rt.to_csv(root / "scaled_harmonized_runtime.tsv", sep="\t", index=False)
    plot_summary(root, corr_df, rt)
    plot_10traits(root)
    write_report(root, corr_df, rt, metrics)


def make_metrics(corr_df, counts_df):
    if corr_df.empty or counts_df.empty:
        return pd.DataFrame()
    rows = []
    key_comps = ["EMmvGWAS_vs_GEMMA_logp", "CosmicScanner_gemma_full_vs_GEMMA_logp", "CosmicScanner_emmv_compatible_vs_EMmvGWAS_logp"]
    for (set_name, traits), g in corr_df.groupby(["set", "traits"]):
        c = counts_df[(counts_df["set"] == set_name) & (counts_df["traits"] == traits)].iloc[0]
        vals = {r.comparison: r.pearson for r in g.itertuples()}
        join_ok = c.inner_join_snps >= 0.95 * c.input_snps
        corr_ok = all(vals.get(k, np.nan) >= 0.95 for k in key_comps)
        status = "PASS" if join_ok and corr_ok else ("WARN" if join_ok else "FAIL")
        rows.append({
            "set": set_name, "traits": traits, "input_snps": int(c.input_snps),
            "gemma_rows": int(c.gemma_rows), "emmvgwas_rows": int(c.emmvgwas_rows),
            "cosmic_gemma_full_rows": int(c.cosmic_gemma_full_rows), "cosmic_emmv_rows": int(c.cosmic_emmv_rows),
            "inner_join_snps": int(c.inner_join_snps),
            "EMmvGWAS_vs_GEMMA_logp_pearson": vals.get("EMmvGWAS_vs_GEMMA_logp", np.nan),
            "CS_gemma_full_vs_GEMMA_logp_pearson": vals.get("CosmicScanner_gemma_full_vs_GEMMA_logp", np.nan),
            "CS_emmv_vs_EMmvGWAS_logp_pearson": vals.get("CosmicScanner_emmv_compatible_vs_EMmvGWAS_logp", np.nan),
            "status": status,
        })
    return pd.DataFrame(rows)


def plot_summary(root, corr_df, rt):
    if not corr_df.empty:
        piv = corr_df.pivot_table(index="traits", columns=["set", "comparison"], values="pearson")
        plt.figure(figsize=(12, 5))
        x = np.arange(len(corr_df))
        labels = [f"{r.set}\nT{r.traits}\n{r.comparison.split('_vs_')[0]}" for r in corr_df.itertuples()]
        plt.bar(x, corr_df["pearson"])
        plt.axhline(0.95, color="red", lw=1)
        plt.xticks(x, labels, rotation=75, ha="right", fontsize=7)
        plt.ylabel("Pearson r (-log10 p)")
        plt.tight_layout()
        plt.savefig(root.parent / "plots/scaled_harmonized_consistency_grid.png", dpi=200)
        plt.close()
    if not rt.empty:
        plt.figure(figsize=(7, 4))
        for method, g in rt.groupby("method"):
            gg = g[g["set"] == sorted(g["set"].unique())[0]].sort_values("traits")
            plt.plot(gg["traits"], gg["seconds"] / 3600.0, marker="o", label=method)
        plt.xlabel("Number of traits")
        plt.ylabel("Runtime hours")
        plt.legend(frameon=False, fontsize=8)
        plt.tight_layout()
        plt.savefig(root.parent / "plots/scaled_harmonized_runtime.png", dpi=200)
        plt.close()


def load_p(path, method):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"))
    if method == "gemma":
        return df.rename(columns={"rs": "SNP", "p_wald": "P"})[["SNP", "P"]]
    if method == "emmv":
        return df.rename(columns={"SNP_ID": "SNP", "p": "P"})[["SNP", "P"]]
    return df[["SNP", "P"]]


def plot_10traits(root):
    candidates = sorted(root.glob("*/traits_10/results/joined.tsv"))
    if not candidates:
        return
    joined_path = candidates[-1]
    set_dir = joined_path.parents[2]
    res_dir = joined_path.parent
    joined = pd.read_csv(joined_path, sep="\t")
    specs = [("P_gemma", "P_emmv", "GEMMA vs EMmvGWAS"), ("P_gemma", "P_cs_gemma_full", "GEMMA vs CS gemma-full"), ("P_emmv", "P_cs_emmv", "EMmvGWAS vs CS emmv")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, (xcol, ycol, title) in zip(axes, specs):
        x = -np.log10(pd.to_numeric(joined[xcol], errors="coerce").clip(lower=1e-300))
        y = -np.log10(pd.to_numeric(joined[ycol], errors="coerce").clip(lower=1e-300))
        ok = np.isfinite(x) & np.isfinite(y)
        r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else np.nan
        ax.scatter(x[ok], y[ok], s=10, alpha=0.55)
        ax.set_title(f"{title}\nr={r:.4f}")
        ax.set_xlabel("x -log10(p)")
        ax.set_ylabel("y -log10(p)")
    plt.tight_layout()
    plt.savefig(root.parent / "plots/scaled_harmonized_scatter_10traits.png", dpi=200)
    plt.close()
    snp_info = pd.read_csv(set_dir / "snp_info.tsv", sep="\t")
    method_files = [
        ("GEMMA", res_dir / "gemma/output/gemma.assoc.txt", "gemma"),
        ("EMmvGWAS", res_dir / "emmvgwas/emmv_association_with_ids.csv", "emmv"),
        ("CS gemma-full", res_dir / "cosmicscanner_gemma_full/cosmic_gemma_full.assoc.mv.txt", "cosmic"),
        ("CS emmv", res_dir / "cosmicscanner_emmv/cosmic_emmv.assoc.mv.txt", "cosmic"),
    ]
    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    threshold = -np.log10(0.05 / max(1, len(snp_info)))
    for ax, (label, path, method) in zip(axes, method_files):
        df = load_p(path, method).merge(snp_info[["SNP", "CHR", "BP"]], on="SNP", how="left")
        y = -np.log10(pd.to_numeric(df["P"], errors="coerce").clip(lower=1e-300))
        x = np.arange(len(df))
        colors = ["#2458a6" if int(c) % 2 else "#222222" for c in df["CHR"].fillna(0)]
        ax.scatter(x, y, s=8, c=colors, alpha=0.75)
        ax.axhline(threshold, color="red", lw=1)
        ax.set_ylabel(label)
    axes[-1].set_xlabel("Harmonized SNP order")
    plt.tight_layout()
    plt.savefig(root.parent / "plots/scaled_harmonized_manhattan_10traits.png", dpi=200)
    plt.close()


def write_report(root, corr_df, rt, metrics):
    lines = ["# Scaled Harmonized Benchmark Report", "", "## Scope", "", "This benchmark uses harmonized genotype, phenotype, and derived PLINK/matrix inputs. Raw BED full benchmark is not used directly because high hard-call missingness causes method-specific SNP filtering and imputation differences.", "", "## Prior Parity", "", "The 200-SNP TA+EDL harmonized-input parity test passed: CosmicScanner gemma-full vs GEMMA and emmv-compatible vs EMmvGWAS both exceeded r=0.95, and p=1 multivariate parity passed.", "", "## Results", ""]
    if not corr_df.empty:
        lines += ["```tsv", corr_df.to_csv(sep="\t", index=False).strip(), "```", ""]
    if not metrics.empty:
        lines += ["## Metrics And Status", "", "```tsv", metrics.to_csv(sep="\t", index=False).strip(), "```", ""]
    if not rt.empty:
        lines += ["## Runtime", "", "```tsv", rt.to_csv(sep="\t", index=False).strip(), "```", ""]
    status = "PASS" if (not metrics.empty and (metrics["status"] == "PASS").all()) else "WARN"
    lines += ["## Decision", "", f"Status: **{status}** for completed scaled harmonized combinations.", "", "`lowmiss_5k` was not run by default. Run with `RUN_LOWMISS_5K=1` only after deciding whether WARN combinations are acceptable or after tightening the scaled SNP selection.", "", "Next: consider `geno.txt.gz`/dosage input before any paper-like full benchmark."]
    (root / "scaled_harmonized_benchmark_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build-set")
    p.add_argument("--raw-prefix", required=True); p.add_argument("--vmiss", required=True); p.add_argument("--samples", required=True); p.add_argument("--pheno", required=True); p.add_argument("--out-dir", required=True); p.add_argument("--set-name", required=True); p.add_argument("--n-snps", type=int, default=1000); p.add_argument("--max-missing-rate", type=float, default=0.0); p.add_argument("--min-maf", type=float, default=0.01); p.set_defaults(func=build_set)
    p = sub.add_parser("summarize-combo")
    p.add_argument("--out-dir", required=True); p.add_argument("--set-name", required=True); p.add_argument("--traits", type=int, required=True); p.add_argument("--input-snps", type=int, required=True); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic-gemma-full", required=True); p.add_argument("--cosmic-emmv", required=True); p.set_defaults(func=summarize_combo)
    p = sub.add_parser("collect"); p.add_argument("--out-dir", required=True); p.set_defaults(func=collect)
    args = ap.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
