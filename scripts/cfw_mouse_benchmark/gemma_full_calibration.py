#!/usr/bin/env python3
"""
gemma-full correlated-trait calibration.
Optimized: vectorized per-SNP GLS via precomputed per-sample sigma_inv.
"""

import numpy as np, pandas as pd, warnings, json, sys
from pathlib import Path
from scipy.linalg import cholesky, solve_triangular, cho_solve, eigh
from scipy.stats import chi2

BASE = Path("results/cfw_mouse_benchmark/scaled_harmonized/missing0_1k")
OUT = Path("results/cfw_mouse_benchmark/scaled_harmonized/debug_4traits/gemma_full_calibration")
OUT.mkdir(parents=True, exist_ok=True)

np.set_printoptions(precision=6, suppress=True, linewidth=120)

VG_GEMMA = np.array([[1.0095e-05,  2.93107e-06,  1.78963e-09, -1.29914e-09],
                     [2.93107e-06,  0.436493,     2.08598e-05, -3.80175e-05],
                     [1.78963e-09,  2.08598e-05,  9.97135e-06, -1.06616e-08],
                     [-1.29914e-09,-3.80175e-05, -1.06616e-08,  9.96961e-06]])
VE_GEMMA = np.array([[1.00967,  0.690724, 0.536205, 0.53871 ],
                     [0.690724, 0.991088, 0.581349, 0.565309],
                     [0.536205, 0.581349, 0.996293, 0.905716],
                     [0.53871,  0.565309, 0.905716, 0.995642]])

def load_data(trait_set="original"):
    G = pd.read_csv(BASE/"genotype.tsv", sep="\t", index_col=0).values.astype(np.float64)
    if trait_set == "original":
        Y = pd.read_csv(BASE/"traits_4/pheno.tsv", sep="\t").iloc[:,2:].values.astype(np.float64)
    elif trait_set == "alt_tail_tibia":
        Y = pd.read_csv(BASE/"traits_4/pheno.tsv", sep="\t").iloc[:,2:6].values.astype(np.float64)
        p6 = pd.read_csv(BASE/"traits_6/pheno.tsv", sep="\t")
        Y[:,2] = p6.iloc[:,4].values; Y[:,3] = p6.iloc[:,5].values
    elif trait_set == "alt_testis_gastroc":
        Y = pd.read_csv(BASE/"traits_4/pheno.tsv", sep="\t").iloc[:,2:6].values.astype(np.float64)
        p8 = pd.read_csv(BASE/"traits_8/pheno.tsv", sep="\t")
        Y[:,2] = p8.iloc[:,6].values; Y[:,3] = p8.iloc[:,7].values
    elif trait_set == "traits_2":
        Y = pd.read_csv(BASE/"traits_2/pheno.tsv", sep="\t").iloc[:,2:].values.astype(np.float64)
    elif trait_set == "traits_10":
        Y = pd.read_csv(BASE/"traits_10/pheno.tsv", sep="\t").iloc[:,2:].values.astype(np.float64)
    else:
        raise ValueError(trait_set)
    snp_ids = pd.read_csv(BASE/"genotype.tsv", sep="\t", index_col=0).columns.tolist()
    return G, Y, snp_ids

def precompute_sigma_inv(n, eigvals, Vg, Ve):
    """Precompute per-sample Sigma^{-1} as (n,d,d) tensor."""
    d = Vg.shape[0]
    D3 = eigvals[:, None, None]
    Sigma = D3 * Vg[None, :, :] + Ve[None, :, :]
    for i in range(n):
        np.fill_diagonal(Sigma[i], Sigma[i].diagonal() + 1e-8)
    SI = np.zeros((n, d, d))
    for i in range(n):
        L = cholesky(Sigma[i], lower=True)
        SI[i] = cho_solve((L, True), np.eye(d))
    return SI

def run_reference_scan(Y, G, eigvals, eigvecs, trait_name, Vg, Ve, label="ref"):
    n, d = Y.shape; m = G.shape[1]
    U = eigvecs
    UtY = U.T @ Y                    # n x d
    Ut_int = U.T @ np.ones(n)        # n
    Ut_all = (U.T @ G).T             # m x n — precompute all rotated SNPs

    # Precompute per-sample sigma_inv
    SI = precompute_sigma_inv(n, eigvals, Vg, Ve)

    # Intercept block (constant across SNPs)
    H_ii = np.tensordot(Ut_int**2, SI, axes=([0],[0]))
    np.fill_diagonal(H_ii, H_ii.diagonal() + 1e-8)
    rhs_i = np.einsum('i,ijk,ik->j', Ut_int, SI, UtY)

    ones_n = np.ones(n)
    rows = []
    for s in range(m):
        Utx = Ut_all[s]  # n
        H_is = np.einsum('i,i,ijk->jk', Ut_int, Utx, SI)
        H_ss = np.einsum('i,i,ijk->jk', ones_n, Utx**2, SI)
        rhs_s = np.einsum('i,ijk,ik->j', Utx, SI, UtY)

        H_full = np.block([[H_ii, H_is], [H_is.T, H_ss]])
        np.fill_diagonal(H_full, H_full.diagonal() + 1e-8)
        rhs_full = np.concatenate([rhs_i, rhs_s])

        try:
            Hinv = np.linalg.solve(H_full, np.eye(2*d))
        except np.linalg.LinAlgError:
            rows.append(dict(SNP=s, wald_stat=np.nan, p_wald=1.0)); continue
        beta_all = Hinv @ rhs_full
        if not np.all(np.isfinite(beta_all)):
            rows.append(dict(SNP=s, wald_stat=np.nan, p_wald=1.0)); continue

        beta = beta_all[d:2*d]; Vbeta = Hinv[d:2*d, d:2*d]
        if np.linalg.cond(Vbeta) > 1e12:
            rows.append(dict(SNP=s, wald_stat=np.nan, p_wald=1.0)); continue
        try:
            Vb_inv = np.linalg.solve(Vbeta, np.eye(d))
        except np.linalg.LinAlgError:
            rows.append(dict(SNP=s, wald_stat=np.nan, p_wald=1.0)); continue
        wald = (beta[:,None].T @ Vb_inv @ beta[:,None]).item()
        if not np.isfinite(wald) or wald < 0:
            rows.append(dict(SNP=s, wald_stat=np.nan, p_wald=1.0)); continue
        p = 1.0 - chi2.cdf(wald, d)
        if p == 0: p = np.finfo(np.float64).tiny

        r = dict(SNP=s, n=n, wald_stat=wald, p_wald=p,
                 genotype_mean=G[:,s].mean(), genotype_variance=G[:,s].var(ddof=0),
                 maf=min(G[:,s].mean()/2, 1-G[:,s].mean()/2))
        for ti in range(d):
            r[f"beta_{ti+1}"] = beta[ti]
            r[f"se_{ti+1}"] = np.sqrt(Vbeta[ti,ti]) if Vbeta[ti,ti] > 0 else np.nan
        rows.append(r)

    snp_ids = pd.read_csv(BASE/"genotype.tsv", sep="\t", index_col=0).columns.tolist()
    df = pd.DataFrame(rows)
    df["SNP"] = [snp_ids[s] if isinstance(s, int) else s for s in df["SNP"]]
    df.to_csv(OUT/f"reference_wald_{trait_name}_{label}.tsv", sep="\t", index=False)
    return df

# ── Loaders ──
def load_gemma_out(trait_dir="traits_4"):
    df = pd.read_csv(BASE/trait_dir/"results/gemma/output/gemma.assoc.txt", sep="\t")
    df.rename(columns={"rs":"SNP","p_wald":"P_GEMMA"}, inplace=True)
    return df[["SNP","P_GEMMA"]]

def load_cosmic_out(subdir):
    df = pd.read_csv(BASE/subdir/"cosmic_gemma_full.assoc.mv.txt", sep="\t")
    df.rename(columns={"P_WALD":"P_CS","WALD_STAT":"WALD_CS"}, inplace=True)
    return df[["SNP","P_CS","WALD_CS"]]

def load_emmv_out(subdir):
    df = pd.read_csv(BASE/subdir/"cosmic_emmv.assoc.mv.txt", sep="\t")
    df.rename(columns={"P":"P_CS_emmv","WALD_STAT":"WALD_CS_emmv"}, inplace=True)
    return df[["SNP","P_CS_emmv","WALD_CS_emmv"]]

# ── Diagnostics ──
def diagnose_covariance(Y, G, eigvals, eigvecs, Vg, Ve):
    n, d = Y.shape; m = G.shape[1]
    U = eigvecs; UtY = U.T @ Y; Ut_int = U.T @ np.ones(n)
    Ut_all = (U.T @ G).T
    SI = precompute_sigma_inv(n, eigvals, Vg, Ve)

    H_ii = np.tensordot(Ut_int**2, SI, axes=([0],[0]))
    np.fill_diagonal(H_ii, H_ii.diagonal()+1e-8)
    rhs_i = np.einsum('i,ijk,ik->j', Ut_int, SI, UtY)
    ones_n = np.ones(n)

    snp_ids = pd.read_csv(BASE/"genotype.tsv", sep="\t", index_col=0).columns.tolist()
    rows = []
    for s in range(m):
        Utx = Ut_all[s]
        H_is = np.einsum('i,i,ijk->jk', Ut_int, Utx, SI)
        H_ss = np.einsum('i,i,ijk->jk', ones_n, Utx**2, SI)
        rhs_s = np.einsum('i,ijk,ik->j', Utx, SI, UtY)
        H_full = np.block([[H_ii, H_is],[H_is.T, H_ss]])
        np.fill_diagonal(H_full, H_full.diagonal()+1e-8)
        try:
            Hinv = np.linalg.solve(H_full, np.eye(2*d))
        except: continue
        beta_all = Hinv @ np.concatenate([rhs_i, rhs_s])
        if not np.all(np.isfinite(beta_all)): continue
        Vbeta = Hinv[d:2*d, d:2*d]
        ev = np.linalg.eigvalsh(Vbeta)
        cn = ev[-1]/ev[0] if ev[0] > 1e-15 else np.inf
        rows.append(dict(SNP=snp_ids[s], cond_Vbeta=cn, det_Vbeta=np.linalg.det(Vbeta),
                         min_eig=ev[0], max_eig=ev[-1]))
    df = pd.DataFrame(rows)

    gemma = load_gemma_out().rename(columns={"P_GEMMA":"P_gemma"})
    cs = load_cosmic_out("traits_4/results/cosmicscanner_gemma_full")
    emmv = load_emmv_out("traits_4/results/cosmicscanner_emmv")
    ref = pd.read_csv(OUT/"reference_wald_original_ref.tsv", sep="\t").rename(columns={"p_wald":"P_ref"})

    all_df = gemma.merge(cs, on="SNP").merge(emmv, on="SNP").merge(ref, on="SNP").merge(df, on="SNP")
    all_df["logp_diff_cs_gemma"] = -np.log10(np.maximum(all_df["P_CS"],1e-300)) - (-np.log10(np.maximum(all_df["P_gemma"],1e-300)))
    all_df.to_csv(OUT/"covariance_stability_per_snp.tsv", sep="\t", index=False)

    ok = np.isfinite(all_df["cond_Vbeta"]) & np.isfinite(all_df["logp_diff_cs_gemma"])
    corr_cond = np.corrcoef(np.log10(all_df.loc[ok,"cond_Vbeta"]+1), np.abs(all_df.loc[ok,"logp_diff_cs_gemma"]))[0,1] if ok.sum()>3 else np.nan

    report = f"""# Covariance Stability Report
SNPs: {len(df)}, Traits: {d}
Ve condition number: {np.linalg.cond(Ve):.4f}
Vg condition number: {np.linalg.cond(Vg):.4f}
Vg+Ve condition number: {np.linalg.cond(Vg+Ve):.4f}
Vβ cond median={df['cond_Vbeta'].median():.2f}, 90%={df['cond_Vbeta'].quantile(0.9):.2f}, max={df['cond_Vbeta'].max():.2f}
Frac cond>100: {(df['cond_Vbeta']>100).mean():.4f}
Corr(log10(cond Vβ), |Δlog10p|): {corr_cond:.4f}
"""
    (OUT/"covariance_stability_report.md").write_text(report)
    return all_df

def run_inversion_sweep(Y, G, eigvals, eigvecs, Vg, Ve, trait_name):
    n, d = Y.shape; m = G.shape[1]
    U = eigvecs; UtY = U.T @ Y; Ut_int = U.T @ np.ones(n)
    Ut_all = (U.T @ G).T
    SI = precompute_sigma_inv(n, eigvals, Vg, Ve)

    H_ii = np.tensordot(Ut_int**2, SI, axes=([0],[0]))
    np.fill_diagonal(H_ii, H_ii.diagonal()+1e-8)
    rhs_i = np.einsum('i,ijk,ik->j', Ut_int, SI, UtY)
    ones_n = np.ones(n)

    snp_ids = pd.read_csv(BASE/"genotype.tsv", sep="\t", index_col=0).columns.tolist()
    gemma_dict = load_gemma_out().set_index("SNP")["P_GEMMA"].to_dict()

    strategies = {
        "solve": lambda V: np.linalg.solve(V, np.eye(V.shape[0])),
        "eigen_floor_1e_8": lambda V: _ef(V,1e-8),
        "eigen_floor_1e_6": lambda V: _ef(V,1e-6),
        "eigen_floor_1e_4": lambda V: _ef(V,1e-4),
        "ridge_1e_8": lambda V: np.linalg.solve(V+1e-8*np.trace(V)/V.shape[0]*np.eye(V.shape[0]), np.eye(V.shape[0])),
        "ridge_1e_6": lambda V: np.linalg.solve(V+1e-6*np.trace(V)/V.shape[0]*np.eye(V.shape[0]), np.eye(V.shape[0])),
        "ridge_1e_4": lambda V: np.linalg.solve(V+1e-4*np.trace(V)/V.shape[0]*np.eye(V.shape[0]), np.eye(V.shape[0])),
        "pinv": lambda V: np.linalg.pinv(V),
    }

    snp_vbeta = []
    for s in range(m):
        Utx = Ut_all[s]
        H_is = np.einsum('i,i,ijk->jk', Ut_int, Utx, SI)
        H_ss = np.einsum('i,i,ijk->jk', ones_n, Utx**2, SI)
        rhs_s = np.einsum('i,ijk,ik->j', Utx, SI, UtY)
        H_full = np.block([[H_ii, H_is],[H_is.T, H_ss]])
        np.fill_diagonal(H_full, H_full.diagonal()+1e-8)
        try:
            Hinv = np.linalg.solve(H_full, np.eye(2*d))
        except: continue
        beta_all = Hinv @ np.concatenate([rhs_i, rhs_s])
        if not np.all(np.isfinite(beta_all)): continue
        beta = beta_all[d:2*d]; Vbeta = Hinv[d:2*d, d:2*d]
        snp_vbeta.append((snp_ids[s], gemma_dict.get(snp_ids[s], np.nan), beta, Vbeta))

    rows = []
    for sid, p_ge, beta, Vbeta in snp_vbeta:
        for sname, sfunc in strategies.items():
            try:
                Vb_inv = sfunc(Vbeta)
                wald = (beta[:,None].T @ Vb_inv @ beta[:,None]).item()
                if not np.isfinite(wald) or wald < 0: raise ValueError
                p = 1.0 - chi2.cdf(wald, Vbeta.shape[0])
                if p == 0: p = np.finfo(np.float64).tiny
            except Exception:
                wald, p = np.nan, np.nan
            rows.append(dict(SNP=sid, strategy=sname, wald=wald, p=p, p_gemma=p_ge))

    df = pd.DataFrame(rows)
    df.to_csv(OUT/f"inversion_strategy_sweep_{trait_name}.tsv", sep="\t", index=False)

    # Per-strategy summary
    summ = df.groupby("strategy").agg(n=("p","count")).reset_index()
    ps = []
    for sname in df["strategy"].unique():
        sub = df[(df["strategy"]==sname)&np.isfinite(df["p"])&np.isfinite(df["p_gemma"])]
        if len(sub) < 3: continue
        lp_ref = -np.log10(np.maximum(sub["p"].values,1e-300))
        lp_ge = -np.log10(np.maximum(sub["p_gemma"].values,1e-300))
        ps.append(dict(strategy=sname, n=len(sub),
                       pearson=np.corrcoef(lp_ref,lp_ge)[0,1],
                       rmse=np.sqrt(np.mean((lp_ref-lp_ge)**2)),
                       max_abs_diff=np.max(np.abs(lp_ref-lp_ge)),
                       fail_rate=(~np.isfinite(df.loc[df["strategy"]==sname,"p"])).mean()))
    summ = pd.DataFrame(ps)
    summ.to_csv(OUT/f"inversion_strategy_summary_{trait_name}.tsv", sep="\t", index=False)

    report = f"# Inversion Strategy Report — {trait_name}\n\n```tsv\n{summ.to_string()}\n```\n"
    (OUT/f"inversion_strategy_report_{trait_name}.md").write_text(report)
    return df

def _ef(V, f):
    ev, evec = np.linalg.eigh(V); ev = np.maximum(ev, f)
    return evec @ np.diag(1.0/ev) @ evec.T

def run_whitening_test():
    G, Y, _ = load_data("original")
    K = np.loadtxt(BASE/"traits_4/results/gemma/kinship/output/kinship.cXX.txt")
    ev, evec = eigh(K); ev = np.maximum(ev, 0.0)

    V_total = VG_GEMMA + VE_GEMMA
    L = cholesky(V_total, lower=True)
    L_inv = solve_triangular(L, np.eye(4), lower=True)

    Y_w = Y @ L_inv.T
    Vg_w = L_inv @ VG_GEMMA @ L_inv.T
    Ve_w = L_inv @ VE_GEMMA @ L_inv.T

    df_w = run_reference_scan(Y_w, G, ev, evec, "whitened_original", Vg_w, Ve_w, label="whitened")
    df_orig = pd.read_csv(OUT/"reference_wald_original_ref.tsv", sep="\t")
    merged = df_orig.rename(columns={"p_wald":"p_orig","wald_stat":"w_orig"}).merge(
        df_w.rename(columns={"p_wald":"p_w","wald_stat":"w_w"}), on="SNP")
    ok = np.isfinite(merged["w_orig"]) & np.isfinite(merged["w_w"])

    gemma = load_gemma_out().set_index("SNP")["P_GEMMA"]
    def _pearson_vs_ge(pcol):
        m = merged.merge(gemma.to_frame(), left_on="SNP", right_index=True)
        lp1 = -np.log10(np.maximum(m[pcol].values,1e-300))
        lp2 = -np.log10(np.maximum(m["P_GEMMA"].values,1e-300))
        ok2 = np.isfinite(lp1)&np.isfinite(lp2)
        return np.corrcoef(lp1[ok2],lp2[ok2])[0,1] if ok2.sum()>3 else np.nan

    metrics = dict(
        pearson_wald_vs_ge_orig=_pearson_vs_ge("p_orig"),
        pearson_wald_vs_ge_whitened=_pearson_vs_ge("p_w"),
        max_abs_diff_wald=float(np.max(np.abs(merged.loc[ok,"w_orig"]-merged.loc[ok,"w_w"]))),
        rmse_wald=float(np.sqrt(np.mean((merged.loc[ok,"w_orig"]-merged.loc[ok,"w_w"])**2))),
    )
    pd.DataFrame([metrics]).to_csv(OUT/"whitening_test_metrics.tsv", sep="\t", index=False)
    report = f"# Whitening Test Report\n\n{json.dumps(metrics,indent=2)}\n"
    (OUT/"whitening_test_report.md").write_text(report)
    return metrics

# ── Main ──
if __name__ == "__main__":
    print("=== gemma-full correlated-trait calibration ===\n")
    G, Y, snp_ids = load_data("original")
    K = np.loadtxt(BASE/"traits_4/results/gemma/kinship/output/kinship.cXX.txt")
    eigvals, eigvecs = eigh(K); eigvals = np.maximum(eigvals, 0.0)
    print(f"Data: {Y.shape[0]} samples, {G.shape[1]} SNPs, {Y.shape[1]} traits")
    print(f"Ve cond={np.linalg.cond(VE_GEMMA):.2f}, Vg+Ve cond={np.linalg.cond(VG_GEMMA+VE_GEMMA):.2f}")

    # Stage 2: Reference Wald
    print("\n--- Stage 2: Reference Wald ---")
    df_ref = run_reference_scan(Y, G, eigvals, eigvecs, "original", VG_GEMMA, VE_GEMMA)

    # Stage 3: Comparison
    print("--- Stage 3: Comparisons ---")
    gemma = load_gemma_out()
    cs = load_cosmic_out("traits_4/results/cosmicscanner_gemma_full")
    emmv = load_emmv_out("traits_4/results/cosmicscanner_emmv")
    ref = pd.read_csv(OUT/"reference_wald_original_ref.tsv", sep="\t").rename(columns={"p_wald":"P_ref"})

    # Build unified table
    unified = gemma.merge(cs, on="SNP").merge(emmv, on="SNP").merge(ref, on="SNP")
    unified.to_csv(OUT/"intermediate_all439.tsv", sep="\t", index=False)

    comparisons = [
        ("GEMMA",  gemma.set_index("SNP")["P_GEMMA"]),
        ("CS_gemma_full", cs.set_index("SNP")["P_CS"]),
        ("CS_emmv", emmv.set_index("SNP")["P_CS_emmv"]),
        ("Ref_GEMMA_GRM", ref.set_index("SNP")["P_ref"]),
    ]
    comp_rows = []
    for li, si in comparisons:
        for lj, sj in comparisons:
            if comparisons.index((li,si)) >= comparisons.index((lj,sj)): continue
            common = si.dropna().to_frame("a").join(sj.dropna().to_frame("b"), how="inner")
            if len(common) < 3: continue
            lp1 = -np.log10(np.maximum(common["a"].values,1e-300))
            lp2 = -np.log10(np.maximum(common["b"].values,1e-300))
            ok = np.isfinite(lp1)&np.isfinite(lp2)
            comp_rows.append(dict(comparison=f"{li}_vs_{lj}", n=int(ok.sum()),
                                  pearson=np.corrcoef(lp1[ok],lp2[ok])[0,1],
                                  spearman=pd.Series(lp1[ok]).corr(pd.Series(lp2[ok]),method="spearman"),
                                  rmse_logp=np.sqrt(np.mean((lp1[ok]-lp2[ok])**2)),
                                  max_abs_diff_logp=np.max(np.abs(lp1[ok]-lp2[ok]))))
    comp_df = pd.DataFrame(comp_rows)
    comp_df.to_csv(OUT/"intermediate_comparison_summary.tsv", sep="\t", index=False)
    print(comp_df.to_string(index=False))

    # Additional reference: VanRaden GRM + scaled Vg/Ve
    print("\n  VanRaden GRM reference...")
    K_van = pd.read_csv(BASE/"kinship.tsv", sep="\t", index_col=0).values.astype(np.float64)
    ratio = np.trace(K)/np.trace(K_van)  # GEMMA_VanRaden scaling ratio
    Vg_van = VG_GEMMA / ratio
    Ve_van = VE_GEMMA  # Ve is independent of GRM scaling
    ev_van, evec_van = eigh(K_van); ev_van = np.maximum(ev_van, 0.0)
    ref_van = run_reference_scan(Y, G, ev_van, evec_van, "original", Vg_van, Ve_van, label="vanraden")
    ref_van2 = pd.read_csv(OUT/"reference_wald_original_vanraden.tsv", sep="\t").rename(columns={"p_wald":"P_ref_van"})
    # Compare VanRaden reference vs CS_gemma_full and vs GEMMA
    merged_van = gemma.merge(cs, on="SNP").merge(ref_van2, on="SNP")
    lp_ge = -np.log10(np.maximum(merged_van["P_GEMMA"].values,1e-300))
    lp_cs = -np.log10(np.maximum(merged_van["P_CS"].values,1e-300))
    lp_rv = -np.log10(np.maximum(merged_van["P_ref_van"].values,1e-300))
    ok = np.isfinite(lp_ge)&np.isfinite(lp_cs)&np.isfinite(lp_rv)
    print(f"  VanRaden ref vs GEMMA: r={np.corrcoef(lp_rv[ok],lp_ge[ok])[0,1]:.4f}")
    print(f"  VanRaden ref vs CS:    r={np.corrcoef(lp_rv[ok],lp_cs[ok])[0,1]:.4f}")
    # Add to unified table
    unified = unified.merge(ref_van2, on="SNP", how="outer")
    unified.to_csv(OUT/"intermediate_all439.tsv", sep="\t", index=False)

    # Stage 4: Covariance diagnostics
    print("\n--- Stage 4: Covariance Diagnostics ---")
    stab_df = diagnose_covariance(Y, G, eigvals, eigvecs, VG_GEMMA, VE_GEMMA)
    print(f"Median Vβ cond: {stab_df['cond_Vbeta'].median():.2f}")

    # Stage 5: Inversion sweep
    print("\n--- Stage 5: Inversion Sweep ---")
    for tset in ["original", "alt_tail_tibia", "alt_testis_gastroc", "traits_2", "traits_10"]:
        print(f"  {tset}...", end=" ")
        try:
            g2, y2, _ = load_data(tset)
            nt = y2.shape[1]
            vg2 = VG_GEMMA[:nt,:nt]; ve2 = VE_GEMMA[:nt,:nt]
            run_inversion_sweep(y2, g2, eigvals, eigvecs, vg2, ve2, tset)
            summ = OUT/f"inversion_strategy_summary_{tset}.tsv"
            if summ.exists():
                print(f"best Pearson: {pd.read_csv(summ,sep='\t')['pearson'].max():.4f}")
            else:
                print("done")
        except Exception as e:
            print(f"ERROR: {e}")

    # Stage 6: Whitening test
    print("\n--- Stage 6: Whitening Test ---")
    wt = run_whitening_test()
    print(json.dumps(wt, indent=2))

    print("\n=== Calibration complete ===")
