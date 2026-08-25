#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


def write_bed_from_matrix(prefix, samples, snp_info, G):
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


def select_cmd(args):
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(args.samples, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    sample_iids = samples["IID"].astype(str).tolist()
    vmiss = pd.read_csv(args.vmiss, sep=r"\s+")
    snp_col = "ID" if "ID" in vmiss.columns else "SNP"
    miss_col = "F_MISS" if "F_MISS" in vmiss.columns else vmiss.columns[-1]
    vmiss = vmiss.rename(columns={snp_col: "SNP", miss_col: "F_MISS"})
    zero = vmiss[vmiss["F_MISS"] == 0].copy().sort_values(["F_MISS", "SNP"])
    lt01 = vmiss[vmiss["F_MISS"] < 0.01].copy().sort_values(["F_MISS", "SNP"])
    pools = [(zero, "missing rate = 0"), (lt01, "missing rate < 0.01"), (vmiss.sort_values(["F_MISS", "SNP"]), "lowest missing rate fallback")]
    selected = None
    G = None
    snp_info = None
    rule = ""
    for pool, base_rule in pools:
        candidate_snps = pool["SNP"].astype(str).head(max(args.n_snps * 50, args.n_snps)).tolist()
        G0, info0 = decode_bed(args.raw_prefix, sample_iids, candidate_snps)
        freq = np.nanmean(G0, axis=0) / 2.0
        maf = np.minimum(freq, 1.0 - freq)
        missing_counts0 = np.isnan(G0).sum(axis=0)
        ok = np.where((maf >= args.min_maf) & (missing_counts0 <= args.max_missing_count))[0]
        if len(ok) >= args.n_snps:
            take = ok[:args.n_snps]
            selected = pool.set_index("SNP").loc[[candidate_snps[i] for i in take]].reset_index()
            G = G0[:, take]
            snp_info = info0.iloc[take].reset_index(drop=True)
            rule = f"{base_rule}, harmonized MAF >= {args.min_maf}, missing count <= {args.max_missing_count}"
            break
    if selected is None:
        pool, base_rule = pools[-1]
        candidate_snps = pool["SNP"].astype(str).head(max(args.n_snps * 50, args.n_snps)).tolist()
        G0, info0 = decode_bed(args.raw_prefix, sample_iids, candidate_snps)
        freq = np.nanmean(G0, axis=0) / 2.0
        maf = np.minimum(freq, 1.0 - freq)
        missing_counts0 = np.isnan(G0).sum(axis=0)
        order = np.lexsort((-maf, missing_counts0))[:args.n_snps]
        selected = pool.set_index("SNP").loc[[candidate_snps[i] for i in order]].reset_index()
        G = G0[:, order]
        snp_info = info0.iloc[order].reset_index(drop=True)
        rule = f"lowest missing rate fallback with best MAF available; fewer than {args.n_snps} SNPs met MAF/missing gates"
    snps = selected["SNP"].astype(str).tolist()
    missing_before = np.isnan(G)
    imputed_counts = missing_before.sum(axis=0).astype(int)
    means = np.nanmean(G, axis=0)
    G_imp = np.where(np.isnan(G), means[None, :], G)
    geno = pd.DataFrame(G_imp, columns=snps)
    geno.insert(0, "IID", sample_iids)
    geno.to_csv(out / "harmonized_genotype.tsv", sep="\t", index=False, float_format="%.8g")
    pheno = pd.read_csv(args.pheno, sep="\t")
    pheno["IID"] = pheno["IID"].astype(str)
    ph = pheno.set_index("IID").loc[sample_iids].reset_index()
    ph = ph[["FID", "IID", "TA", "EDL"]] if "FID" in ph.columns else ph
    ph.to_csv(out / "harmonized_pheno_TA_EDL.tsv", sep="\t", index=False, na_rep="NA")
    ph[["TA", "EDL"]].to_csv(out / "harmonized_pheno_TA_EDL.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
    ph[["TA", "EDL"]].to_csv(out / "harmonized_pheno_TA_EDL.emmv.P.csv", index=False, na_rep="NA")
    ph[["FID", "IID", "TA"]].to_csv(out / "harmonized_pheno_TA.tsv", sep="\t", index=False, na_rep="NA")
    ph[["TA"]].to_csv(out / "harmonized_pheno_TA.gemma.no_header.tsv", sep="\t", index=False, header=False, na_rep="NA")
    pd.DataFrame({"SNP": snps}).to_csv(out / "harmonized_snps.txt", index=False, header=False)
    samples.to_csv(out / "harmonized_samples.txt", sep="\t", index=False, header=False)
    maf_final = np.minimum(means / 2.0, 1.0 - means / 2.0)
    pd.DataFrame({"SNP": snps, "missing_rate_raw": selected["F_MISS"].to_numpy(), "imputed_count_in_harmonized_samples": imputed_counts, "mean_dosage": means, "maf_harmonized": maf_final}).to_csv(out / "imputation_counts.tsv", sep="\t", index=False)
    snp_info.to_csv(out / "harmonized_snp_info.tsv", sep="\t", index=False)
    z = G_imp - 1.0
    with (out / "harmonized_emmv.Z.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["SNP"] + sample_iids)
        for sid, row in zip(snps, z.T):
            writer.writerow([sid] + [f"{x:.8g}" for x in row])
    Zc = G_imp - G_imp.mean(axis=0, keepdims=True)
    denom = np.sum(2 * (G_imp.mean(axis=0) / 2.0) * (1 - G_imp.mean(axis=0) / 2.0))
    K = (Zc @ Zc.T) / denom if denom > 0 else (Zc @ Zc.T) / G_imp.shape[1]
    pd.DataFrame(K, index=sample_iids, columns=sample_iids).to_csv(out / "harmonized_kinship.tsv", sep="\t")
    write_bed_from_matrix(out / "harmonized", samples, snp_info, G_imp)
    report = [
        "# Harmonized SNP Selection Report", "", f"Selection rule: {rule}", f"Target SNPs: {args.n_snps}",
        f"Selected SNPs: {len(snps)}", f"Raw zero-missing SNPs: {len(zero)}", f"Raw <0.01 missing SNPs: {len(lt01)}", "",
        f"Missing genotypes in harmonized matrix before imputation: {int(missing_before.sum())}", f"Missing genotypes after imputation: {int(np.isnan(G_imp).sum())}",
        f"Minimum harmonized MAF: {float(np.nanmin(maf_final)):.6g}",
    ]
    (out / "snp_selection_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    hashes = {p.name: sha256(p) for p in [out / "harmonized_genotype.tsv", out / "harmonized_pheno_TA_EDL.tsv", out / "harmonized_emmv.Z.csv", out / "harmonized_kinship.tsv", plink_file(out / "harmonized", "bed"), plink_file(out / "harmonized", "bim"), plink_file(out / "harmonized", "fam")]}
    input_report = ["# Harmonized Input Hash Report", "", f"Samples: {len(sample_iids)}", f"SNPs: {len(snps)}", "Traits: 2 (TA, EDL)", f"Genotype matrix has missing after imputation: {bool(np.isnan(G_imp).any())}", "Phenotypes are standardized z-scores from paper_like preprocessing.", "GEMMA, CosmicScanner, and exported matrices are derived from the same harmonized genotype/phenotype sample and SNP order. EMmvGWAS uses `harmonized_emmv.Z.csv` generated from the same matrix as dosage-1.", "", "## SHA256", ""]
    input_report += [f"- {k}: `{v}`" for k, v in hashes.items()]
    (out / "input_hash_report.md").write_text("\n".join(input_report) + "\n", encoding="utf-8")


def load_result(path, method):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"))
    if method == "gemma":
        out = df.rename(columns={"rs": "SNP", "chr": "CHR", "ps": "BP", "p_wald": "P", "beta_1": "BETA_1", "beta_2": "BETA_2"})
        out["STAT"] = np.nan
    elif method == "emmv":
        out = df.rename(columns={"SNP_ID": "SNP", "p": "P", "wald": "STAT", "TA": "BETA_1", "EDL": "BETA_2"})
    else:
        out = df.rename(columns={"WALD_STAT": "STAT"})
    cols = [c for c in ["SNP", "CHR", "BP", "P", "STAT", "BETA_1", "BETA_2"] if c in out.columns]
    return out[cols].copy()


def pearson_spearman(x, y):
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        return int(ok.sum()), np.nan, np.nan
    px = x[ok].to_numpy(); py = y[ok].to_numpy()
    return int(ok.sum()), float(np.corrcoef(px, py)[0, 1]), float(np.corrcoef(pd.Series(px).rank(), pd.Series(py).rank())[0, 1])


def compare_cmd(args):
    out_dir = Path(args.out_dir)
    gemma = load_result(args.gemma, "gemma").add_suffix("_gemma").rename(columns={"SNP_gemma": "SNP"})
    emmv = load_result(args.emmv, "emmv").add_suffix("_emmv").rename(columns={"SNP_emmv": "SNP"})
    cs_gf = load_result(args.cosmic_gemma_full, "cosmic").add_suffix("_cs_gemma_full").rename(columns={"SNP_cs_gemma_full": "SNP"})
    cs_em = load_result(args.cosmic_emmv, "cosmic").add_suffix("_cs_emmv").rename(columns={"SNP_cs_emmv": "SNP"})
    joined = gemma.merge(emmv, on="SNP", how="inner").merge(cs_gf, on="SNP", how="inner").merge(cs_em, on="SNP", how="inner")
    joined.to_csv(out_dir / "harmonized_joined_results.tsv", sep="\t", index=False)
    comps = [
        ("EMmvGWAS_vs_GEMMA_logp", "P_gemma", "P_emmv", "logp"),
        ("CosmicScanner_gemma_full_vs_GEMMA_logp", "P_gemma", "P_cs_gemma_full", "logp"),
        ("CosmicScanner_emmv_compatible_vs_EMmvGWAS_logp", "P_emmv", "P_cs_emmv", "logp"),
        ("CosmicScanner_gemma_full_vs_EMmvGWAS_logp", "P_emmv", "P_cs_gemma_full", "logp"),
        ("CosmicScanner_emmv_compatible_vs_GEMMA_logp", "P_gemma", "P_cs_emmv", "logp"),
    ]
    rows = []
    for name, a, b, mode in comps:
        x = joined[a]; y = joined[b]
        if mode == "logp":
            x = -np.log10(pd.to_numeric(x, errors="coerce").clip(lower=1e-300))
            y = -np.log10(pd.to_numeric(y, errors="coerce").clip(lower=1e-300))
        n, p, s = pearson_spearman(x, y)
        rows.append({"comparison": name, "n": n, "pearson": p, "spearman": s})
    for prefix in ["BETA_1", "BETA_2", "STAT"]:
        for label, a, b in [("emmv_vs_gemma", f"{prefix}_gemma", f"{prefix}_emmv"), ("cs_gemma_full_vs_gemma", f"{prefix}_gemma", f"{prefix}_cs_gemma_full"), ("cs_emmv_vs_emmv", f"{prefix}_emmv", f"{prefix}_cs_emmv")]:
            if a in joined.columns and b in joined.columns:
                n, p, s = pearson_spearman(joined[a], joined[b])
                rows.append({"comparison": f"{prefix}_{label}", "n": n, "pearson": p, "spearman": s})
    pd.DataFrame(rows).to_csv(out_dir / "harmonized_correlation_summary.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    plot_specs = [("P_gemma", "P_emmv", "EMmvGWAS vs GEMMA"), ("P_gemma", "P_cs_gemma_full", "CS gemma-full vs GEMMA"), ("P_emmv", "P_cs_emmv", "CS emmv vs EMmvGWAS")]
    for ax, (xcol, ycol, title) in zip(axes, plot_specs):
        x = -np.log10(pd.to_numeric(joined[xcol], errors="coerce").clip(lower=1e-300))
        y = -np.log10(pd.to_numeric(joined[ycol], errors="coerce").clip(lower=1e-300))
        ok = np.isfinite(x) & np.isfinite(y)
        r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else np.nan
        ax.scatter(x[ok], y[ok], s=12, alpha=0.6)
        ax.set_title(f"{title}\nr={r:.4f}, n={ok.sum()}")
        ax.set_xlabel("x -log10(p)")
        ax.set_ylabel("y -log10(p)")
    plt.tight_layout()
    plt.savefig(Path(args.plot), dpi=200)


def p1_cmd(args):
    gemma = load_result(args.gemma, "gemma").rename(columns={"P": "P_gemma", "BETA_1": "BETA_gemma"})
    uni = load_result(args.cosmic_uni, "cosmic").rename(columns={"P": "P_cosmic_uni", "BETA_1": "BETA_cosmic_uni"})
    mv = load_result(args.cosmic_mv_p1, "cosmic").rename(columns={"P": "P_cosmic_mv_p1", "BETA_1": "BETA_cosmic_mv_p1"})
    joined = gemma.merge(uni, on="SNP", how="inner").merge(mv, on="SNP", how="inner")
    rows = []
    for name, a, b in [("Cosmic_univariate_vs_GEMMA_logp", "P_gemma", "P_cosmic_uni"), ("Cosmic_mv_p1_vs_GEMMA_logp", "P_gemma", "P_cosmic_mv_p1"), ("Cosmic_mv_p1_vs_Cosmic_univariate_logp", "P_cosmic_uni", "P_cosmic_mv_p1")]:
        x = -np.log10(pd.to_numeric(joined[a], errors="coerce").clip(lower=1e-300))
        y = -np.log10(pd.to_numeric(joined[b], errors="coerce").clip(lower=1e-300))
        n, p, s = pearson_spearman(x, y)
        rows.append({"comparison": name, "n": n, "pearson": p, "spearman": s})
    pd.DataFrame(rows).to_csv(Path(args.out_dir) / "mv_p1_parity_summary.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    for ax, (name, a, b) in zip(axes, [("uni vs GEMMA", "P_gemma", "P_cosmic_uni"), ("mv p=1 vs GEMMA", "P_gemma", "P_cosmic_mv_p1"), ("mv p=1 vs uni", "P_cosmic_uni", "P_cosmic_mv_p1")]):
        x = -np.log10(pd.to_numeric(joined[a], errors="coerce").clip(lower=1e-300)); y = -np.log10(pd.to_numeric(joined[b], errors="coerce").clip(lower=1e-300))
        ok = np.isfinite(x) & np.isfinite(y)
        r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 else np.nan
        ax.scatter(x[ok], y[ok], s=12, alpha=0.6); ax.set_title(f"{name}\nr={r:.4f}"); ax.set_xlabel("x -log10(p)"); ax.set_ylabel("y -log10(p)")
    plt.tight_layout(); plt.savefig(args.plot, dpi=200)
    report = ["# Multivariate p=1 Parity Report", "", "```tsv", (Path(args.out_dir) / "mv_p1_parity_summary.tsv").read_text(encoding="utf-8").strip(), "```", ""]
    (Path(args.out_dir) / "mv_p1_parity_report.md").write_text("\n".join(report), encoding="utf-8")


def report_cmd(args):
    out = Path(args.out_dir)
    corr = pd.read_csv(out / "harmonized_correlation_summary.tsv", sep="\t")
    p1 = pd.read_csv(out / "mv_p1_parity_summary.tsv", sep="\t")
    joined_n = int(pd.read_csv(out / "harmonized_joined_results.tsv", sep="\t").shape[0])
    def get(df, key):
        r = df[df.iloc[:, 0] == key]
        return float(r["pearson"].iloc[0]) if len(r) else np.nan
    gates = {
        "inner_join_ge_190": joined_n >= 190,
        "mv_p1_vs_univariate_ge_0.99": get(p1, "Cosmic_mv_p1_vs_Cosmic_univariate_logp") >= 0.99,
        "cs_gemma_full_vs_gemma_ge_0.95": get(corr, "CosmicScanner_gemma_full_vs_GEMMA_logp") >= 0.95,
        "cs_emmv_vs_emmv_ge_0.95": get(corr, "CosmicScanner_emmv_compatible_vs_EMmvGWAS_logp") >= 0.95,
    }
    verdict = "PASS" if all(gates.values()) else "FAIL"
    if joined_n < 190:
        reason = "input / SNP join problem"
    elif not gates["mv_p1_vs_univariate_ge_0.99"]:
        reason = "p=1 multivariate engine problem"
    elif not gates["cs_gemma_full_vs_gemma_ge_0.95"] or not gates["cs_emmv_vs_emmv_ge_0.95"]:
        reason = "multivariate p-value/statistic or Wald implementation mismatch"
    else:
        reason = "all gates passed"
    lines = ["# Harmonized 2-Trait Parity Report", "", f"Verdict: **{verdict}**", f"Reason: {reason}", "", "## Why This Test", "", "The raw CFW BED has high hard-call missingness, so previous method comparisons mixed different SNP filtering and imputation behavior. This test fixes sample order, SNP order, phenotype matrix, and harmonized PLINK/matrix exports before running methods.", "", "## Inputs", "", (out / "input_hash_report.md").read_text(encoding="utf-8"), "", "## Correlations", "", "```tsv", (out / "harmonized_correlation_summary.tsv").read_text(encoding="utf-8").strip(), "```", "", "## p=1 Parity", "", "```tsv", (out / "mv_p1_parity_summary.tsv").read_text(encoding="utf-8").strip(), "```", "", "## Gates", ""]
    lines += [f"- {k}: {'PASS' if v else 'FAIL'}" for k, v in gates.items()]
    lines += ["", "## Next Fix", "", "If FAIL is driven by multivariate paths while p=1 passes, inspect CosmicScanner multivariate Wald statistic/GRM normalization against GEMMA/EMmvGWAS on the same harmonized matrix."]
    (out / "harmonized_2traits_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = Path(args.results_dir) / "debug_summary.md"
    prev = summary.read_text(encoding="utf-8") if summary.exists() else "# CFW Mouse Debug Summary\n"
    block = "\n## Harmonized 2-Trait Parity\n\n" + f"Verdict: **{verdict}** ({reason})\n\n" + "```tsv\n" + (out / "harmonized_correlation_summary.tsv").read_text(encoding="utf-8").strip() + "\n```\n\n" + "p=1:\n\n```tsv\n" + (out / "mv_p1_parity_summary.tsv").read_text(encoding="utf-8").strip() + "\n```\n"
    summary.write_text(prev.split("## Harmonized 2-Trait Parity")[0].rstrip() + block + "\n", encoding="utf-8")
    report = Path(args.results_dir) / "CFW_mouse_benchmark_report.md"
    rep = report.read_text(encoding="utf-8") if report.exists() else "# CFW Mouse Benchmark Report\n"
    block2 = "\n## Harmonized 2-Trait Parity Result\n\n" + f"- Verdict: **{verdict}**. Reason: {reason}.\n" + f"- Report: `{out / 'harmonized_2traits_report.md'}`\n" + "- Full benchmark remains blocked unless verdict is PASS.\n"
    report.write_text(rep.split("## Harmonized 2-Trait Parity Result")[0].rstrip() + block2, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("select"); p.add_argument("--raw-prefix", required=True); p.add_argument("--vmiss", required=True); p.add_argument("--samples", required=True); p.add_argument("--pheno", required=True); p.add_argument("--out-dir", required=True); p.add_argument("--n-snps", type=int, default=200); p.add_argument("--min-maf", type=float, default=0.01); p.add_argument("--max-missing-count", type=int, default=0); p.set_defaults(func=select_cmd)
    p = sub.add_parser("compare"); p.add_argument("--out-dir", required=True); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic-gemma-full", required=True); p.add_argument("--cosmic-emmv", required=True); p.add_argument("--plot", required=True); p.set_defaults(func=compare_cmd)
    p = sub.add_parser("p1"); p.add_argument("--out-dir", required=True); p.add_argument("--gemma", required=True); p.add_argument("--cosmic-uni", required=True); p.add_argument("--cosmic-mv-p1", required=True); p.add_argument("--plot", required=True); p.set_defaults(func=p1_cmd)
    p = sub.add_parser("report"); p.add_argument("--out-dir", required=True); p.add_argument("--results-dir", required=True); p.set_defaults(func=report_cmd)
    args = ap.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
