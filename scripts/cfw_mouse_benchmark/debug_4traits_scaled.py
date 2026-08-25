#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def corr(x, y):
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= 2:
        return int(ok.sum()), np.nan, np.nan, np.nan
    a = x[ok].to_numpy(); b = y[ok].to_numpy()
    pear = float(np.corrcoef(a, b)[0, 1])
    spear = float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])
    # Kendall tau-b without scipy; O(n^2), fine for <=439 SNPs.
    conc = disc = 0
    for i in range(len(a)):
        da = a[i+1:] - a[i]
        db = b[i+1:] - b[i]
        prod = da * db
        conc += int((prod > 0).sum())
        disc += int((prod < 0).sum())
    kend = (conc - disc) / (conc + disc) if (conc + disc) else np.nan
    return int(ok.sum()), pear, spear, float(kend)


def load_result(path, method):
    df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"))
    if method == "gemma":
        return df.rename(columns={"rs": "SNP", "chr": "CHR", "ps": "BP", "p_wald": "P_WALD", "p_lrt": "P_LRT", "p_score": "P_SCORE"})
    if method == "emmv":
        return df.rename(columns={"SNP_ID": "SNP", "p": "P", "wald": "STAT"})
    return df.rename(columns={"P": "P", "WALD_STAT": "STAT"})


def joined_4(args):
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    snp_qc = pd.read_csv(args.snp_qc, sep="\t")
    g = load_result(args.gemma, "gemma")
    e = load_result(args.emmv, "emmv")
    cg = load_result(args.cosmic_gemma_full, "cosmic")
    ce = load_result(args.cosmic_emmv, "cosmic")
    j = snp_qc[["SNP", "maf", "mean_dosage"]].copy()
    geno = pd.read_csv(args.genotype, sep="\t")
    var = geno.drop(columns=["IID"]).var(axis=0, ddof=1)
    j["geno_var"] = j["SNP"].map(var.to_dict())
    j = j.merge(g[[c for c in ["SNP", "CHR", "BP", "P_WALD", "P_LRT", "P_SCORE"] if c in g.columns]], on="SNP", how="inner")
    j = j.merge(e[["SNP", "P", "STAT"]].rename(columns={"P": "P_emmv", "STAT": "STAT_emmv"}), on="SNP", how="inner")
    j = j.merge(cg[["SNP", "P", "STAT"]].rename(columns={"P": "P_cs_gemma_full", "STAT": "STAT_cs_gemma_full"}), on="SNP", how="inner")
    j = j.merge(ce[["SNP", "P", "STAT"]].rename(columns={"P": "P_cs_emmv", "STAT": "STAT_cs_emmv"}), on="SNP", how="inner")
    j = j.rename(columns={"P_WALD": "P_gemma_wald", "P_LRT": "P_gemma_lrt", "P_SCORE": "P_gemma_score"})
    for c in ["P_gemma_wald", "P_gemma_lrt", "P_gemma_score", "P_emmv", "P_cs_gemma_full", "P_cs_emmv"]:
        if c in j.columns:
            j[c.replace("P_", "LOGP_")] = -np.log10(pd.to_numeric(j[c], errors="coerce").clip(lower=1e-300))
    j["absdiff_cs_gemma_full_vs_gemma"] = (j["LOGP_cs_gemma_full"] - j["LOGP_gemma_wald"]).abs()
    j["absdiff_emmv_vs_gemma"] = (j["LOGP_emmv"] - j["LOGP_gemma_wald"]).abs()
    j["absdiff_cs_emmv_vs_emmv"] = (j["LOGP_cs_emmv"] - j["LOGP_emmv"]).abs()
    j.sort_values("absdiff_cs_gemma_full_vs_gemma", ascending=False).to_csv(out / "joined_4traits_debug.tsv", sep="\t", index=False)
    j.sort_values("absdiff_cs_gemma_full_vs_gemma", ascending=False).head(30).to_csv(out / "top_discrepant_snps.tsv", sep="\t", index=False)
    robust_rows = []
    x0 = j["LOGP_gemma_wald"]; y0 = j["LOGP_cs_gemma_full"]
    d = j["absdiff_cs_gemma_full_vs_gemma"]
    for drop in [0, 1, 3, 5, 10]:
        keep = d.rank(method="first", ascending=False) > drop if drop else pd.Series(True, index=j.index)
        n, p, s, k = corr(x0[keep], y0[keep])
        diff = (y0[keep] - x0[keep]).abs()
        robust_rows.append({"comparison":"CS_gemma_full_vs_GEMMA_wald", "drop_top_absdiff":drop, "n":n, "pearson":p, "spearman":s, "kendall":k, "rmse_logp":float(np.sqrt(np.mean((y0[keep]-x0[keep])**2))), "median_absdiff_logp":float(diff.median()), "max_absdiff_logp":float(diff.max())})
    for name, x, y in [("EMmvGWAS_vs_GEMMA", j["LOGP_gemma_wald"], j["LOGP_emmv"]), ("CS_emmv_vs_EMmvGWAS", j["LOGP_emmv"], j["LOGP_cs_emmv"]), ("CS_gemma_full_vs_GEMMA_lrt", j.get("LOGP_gemma_lrt", j["LOGP_gemma_wald"]), j["LOGP_cs_gemma_full"]), ("CS_gemma_full_vs_GEMMA_score", j.get("LOGP_gemma_score", j["LOGP_gemma_wald"]), j["LOGP_cs_gemma_full"])]:
        n, p, s, k = corr(x, y)
        robust_rows.append({"comparison":name, "drop_top_absdiff":0, "n":n, "pearson":p, "spearman":s, "kendall":k, "rmse_logp":float(np.sqrt(np.mean((y-x)**2))), "median_absdiff_logp":float((y-x).abs().median()), "max_absdiff_logp":float((y-x).abs().max())})
    pd.DataFrame(robust_rows).to_csv(out / "robust_correlation_summary.tsv", sep="\t", index=False)
    lines = ["# Discrepant SNP Report", "", "Top discrepant SNPs are ranked by `abs(CS gemma-full - GEMMA Wald)` on -log10(p).", "", "## Summary", ""]
    top = j.sort_values("absdiff_cs_gemma_full_vs_gemma", ascending=False).head(30)
    lines += [f"- Top30 median MAF: {top['maf'].median():.6g}", f"- Top30 min MAF: {top['maf'].min():.6g}", f"- Top30 median genotype variance: {top['geno_var'].median():.6g}", f"- Max abs logp diff CS gemma-full vs GEMMA: {top['absdiff_cs_gemma_full_vs_gemma'].max():.6g}", f"- Median abs logp diff CS emmv-compatible vs EMmvGWAS in Top30: {top['absdiff_cs_emmv_vs_emmv'].median():.6g}"]
    (out / "discrepant_snp_report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    # schema
    schema_lines = ["# P-value Schema 4 Traits", "", "GEMMA columns used: `p_wald` as primary; additional columns `p_lrt` and `p_score` are also compared when present.", "CosmicScanner gemma-full columns: `WALD_STAT`, `P`; CLI requested `--mv-test wald,lrt,score`, but current assoc file exposes primary `P` and `WALD_STAT`.", "EMmvGWAS columns: `wald`, `p`.", ""]
    for label, path in [("GEMMA", args.gemma), ("EMmvGWAS", args.emmv), ("CS gemma-full", args.cosmic_gemma_full)]:
        df = pd.read_csv(path, sep=("," if str(path).endswith(".csv") else "\t"), nrows=3)
        schema_lines += [f"## {label}", "", ", ".join(df.columns), ""]
    (out / "pvalue_schema_4traits.md").write_text("\n".join(schema_lines)+"\n", encoding="utf-8")
    pd.DataFrame([r for r in robust_rows if "GEMMA_" in r["comparison"] or r["comparison"] == "CS_gemma_full_vs_GEMMA_wald"]).to_csv(out / "pvalue_column_correlation_4traits.tsv", sep="\t", index=False)


def covariance(args):
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    ph = pd.read_csv(args.full_pheno, sep="\t")
    sets = json.loads(Path(args.sets_json).read_text(encoding="utf-8"))
    rows = []
    report = ["# Trait Covariance Diagnostics", ""]
    for name, traits in sets.items():
        X = ph[traits].apply(pd.to_numeric, errors="coerce").dropna()
        cov = np.cov(X.to_numpy().T)
        ev = np.linalg.eigvalsh(cov)
        cond = float(ev.max() / ev.min()) if ev.min() > 0 else np.inf
        corrm = X.corr()
        maxcorr = corrm.where(~np.eye(len(traits), dtype=bool)).abs().max().max()
        rows.append({"trait_set": name, "traits": ",".join(traits), "n": len(X), "min_eigen": ev.min(), "max_eigen": ev.max(), "condition_number": cond, "max_abs_pairwise_corr": maxcorr})
        report += [f"## {name}", "", f"Traits: {', '.join(traits)}", f"N: {len(X)}", f"Eigenvalues: {', '.join(f'{v:.6g}' for v in ev)}", f"Condition number: {cond:.6g}", f"Max abs pairwise corr: {maxcorr:.6g}", ""]
    pd.DataFrame(rows).to_csv(out / "trait_covariance_diagnostics.tsv", sep="\t", index=False)
    (out / "trait_covariance_report.md").write_text("\n".join(report), encoding="utf-8")


def summarize_ablation(args):
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for combo_dir in sorted((out / "ablation_runs").glob("combo_*")):
        meta = json.loads((combo_dir / "meta.json").read_text(encoding="utf-8"))
        g = load_result(combo_dir / "gemma/output/gemma.assoc.txt", "gemma")
        e = load_result(combo_dir / "emmvgwas/emmv_association_with_ids.csv", "emmv")
        cg = load_result(combo_dir / "cosmic_gemma_full/cosmic_gemma_full.assoc.mv.txt", "cosmic")
        ce = load_result(combo_dir / "cosmic_emmv/cosmic_emmv.assoc.mv.txt", "cosmic")
        j = g.rename(columns={"P_WALD":"P_gemma"})[["SNP","P_gemma"]].merge(e.rename(columns={"P":"P_emmv"})[["SNP","P_emmv"]]).merge(cg.rename(columns={"P":"P_cg"})[["SNP","P_cg"]]).merge(ce.rename(columns={"P":"P_ce"})[["SNP","P_ce"]])
        for comp, a, b in [("EMmvGWAS_vs_GEMMA", "P_gemma", "P_emmv"), ("CS_gemma_full_vs_GEMMA", "P_gemma", "P_cg"), ("CS_emmv_vs_EMmvGWAS", "P_emmv", "P_ce")]:
            x = -np.log10(pd.to_numeric(j[a], errors="coerce").clip(lower=1e-300)); y = -np.log10(pd.to_numeric(j[b], errors="coerce").clip(lower=1e-300))
            n,p,s,k = corr(x,y)
            rows.append({"combo": meta["name"], "traits": ",".join(meta["traits"]), "comparison": comp, "n": n, "pearson": p, "spearman": s, "kendall": k})
    pd.DataFrame(rows).to_csv(out / "trait_ablation_metrics.tsv", sep="\t", index=False)


def intermediate(args):
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    top = pd.read_csv(out / "top_discrepant_snps.tsv", sep="\t").head(10)
    geno = pd.read_csv(args.genotype, sep="\t")
    ph = pd.read_csv(args.pheno, sep="\t")
    Y = ph[args.traits.split(",")].to_numpy(dtype=float)
    rows = []
    for _, r in top.iterrows():
        sid = r["SNP"]
        x = geno[sid].to_numpy(dtype=float)
        X = np.column_stack([np.ones(len(x)), x])
        B = np.linalg.lstsq(X, Y, rcond=None)[0]
        E = Y - X @ B
        S = (E.T @ E) / (len(x) - X.shape[1])
        xvar = float(np.var(x, ddof=1))
        beta = B[1, :]
        try:
            stat = float(beta @ np.linalg.inv(S / np.sum((x - x.mean())**2)) @ beta)
        except np.linalg.LinAlgError:
            stat = np.nan
        rows.append({"SNP": sid, "sample_size": len(x), "genotype_mean": x.mean(), "genotype_variance": xvar, "reference_beta_vector": ",".join(f"{v:.8g}" for v in beta), "reference_resid_cov": ";".join(",".join(f"{v:.8g}" for v in row) for row in S), "reference_wald_stat_ols": stat, "gemma_logp": r.get("LOGP_gemma_wald", np.nan), "cs_gemma_full_logp": r.get("LOGP_cs_gemma_full", np.nan), "emmv_logp": r.get("LOGP_emmv", np.nan), "cs_emmv_logp": r.get("LOGP_cs_emmv", np.nan)})
    pd.DataFrame(rows).to_csv(out / "intermediate_top10.tsv", sep="\t", index=False)
    (out / "intermediate_top10_report.md").write_text("# Intermediate Top10 Report\n\nA lightweight OLS multivariate Wald reference was computed for diagnostic direction only; it is not a GEMMA LMM replacement. GEMMA internal per-SNP covariance is not exported by GEMMA.\n", encoding="utf-8")


def final_report(args):
    out = Path(args.out_dir)
    robust = pd.read_csv(out / "robust_correlation_summary.tsv", sep="\t")
    ab = pd.read_csv(out / "trait_ablation_metrics.tsv", sep="\t") if (out / "trait_ablation_metrics.tsv").exists() else pd.DataFrame()
    cov = pd.read_csv(out / "trait_covariance_diagnostics.tsv", sep="\t") if (out / "trait_covariance_diagnostics.tsv").exists() else pd.DataFrame()
    main = robust[(robust.comparison=="CS_gemma_full_vs_GEMMA_wald") & (robust.drop_top_absdiff==0)].iloc[0]
    drop10 = robust[(robust.comparison=="CS_gemma_full_vs_GEMMA_wald") & (robust.drop_top_absdiff==10)].iloc[0]
    if drop10.pearson >= 0.95 and main.pearson < 0.95:
        verdict = "A"
        reason = "4-trait WARN is dominated by outlier/Pearson sensitivity."
    elif not ab.empty and (ab[(ab.comparison=="CS_gemma_full_vs_GEMMA") & (ab.pearson < 0.95)].shape[0] > 0):
        verdict = "B"
        reason = "4-trait WARN appears trait-set dependent across ablations."
    else:
        verdict = "C"
        reason = "CosmicScanner gemma-full multivariate path needs further calibration."
    lines = ["# Missing0 4-Trait Debug Report", "", f"Verdict option: **{verdict}**", f"Reason: {reason}", "", "## Robust Correlations", "", "```tsv", robust.to_csv(sep="\t", index=False).strip(), "```", ""]
    if not ab.empty:
        lines += ["## Trait Ablation", "", "```tsv", ab.to_csv(sep="\t", index=False).strip(), "```", ""]
    if not cov.empty:
        lines += ["## Covariance Diagnostics", "", "```tsv", cov.to_csv(sep="\t", index=False).strip(), "```", ""]
    (out / "debug_4traits_report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    for target in [Path(args.scaled_report), Path(args.debug_summary), Path(args.main_report)]:
        prev = target.read_text(encoding="utf-8") if target.exists() else ""
        block = f"\n## Missing0 4-Trait Debug\n\n- Verdict option: **{verdict}**. {reason}\n- Report: `{out / 'debug_4traits_report.md'}`\n- Robust CS gemma-full vs GEMMA Pearson: original={main.pearson:.4f}, drop_top10={drop10.pearson:.4f}.\n"
        target.write_text(prev.split("## Missing0 4-Trait Debug")[0].rstrip()+block, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p=sub.add_parser("joined4"); p.add_argument("--out-dir", required=True); p.add_argument("--snp-qc", required=True); p.add_argument("--genotype", required=True); p.add_argument("--gemma", required=True); p.add_argument("--emmv", required=True); p.add_argument("--cosmic-gemma-full", required=True); p.add_argument("--cosmic-emmv", required=True); p.set_defaults(func=joined_4)
    p=sub.add_parser("covariance"); p.add_argument("--out-dir", required=True); p.add_argument("--full-pheno", required=True); p.add_argument("--sets-json", required=True); p.set_defaults(func=covariance)
    p=sub.add_parser("summarize-ablation"); p.add_argument("--out-dir", required=True); p.set_defaults(func=summarize_ablation)
    p=sub.add_parser("intermediate"); p.add_argument("--out-dir", required=True); p.add_argument("--genotype", required=True); p.add_argument("--pheno", required=True); p.add_argument("--traits", required=True); p.set_defaults(func=intermediate)
    p=sub.add_parser("final-report"); p.add_argument("--out-dir", required=True); p.add_argument("--scaled-report", required=True); p.add_argument("--debug-summary", required=True); p.add_argument("--main-report", required=True); p.set_defaults(func=final_report)
    args=ap.parse_args(); args.func(args)


if __name__ == "__main__":
    main()
