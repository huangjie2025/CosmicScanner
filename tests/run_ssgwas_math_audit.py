#!/usr/bin/env python3
"""Phase 6.5 independent mathematical red-team audit of the ssGBLUP back-solving
workflow (--single-step-method backsolve).

This script deliberately does NOT reuse the Phase 6 analytic reference as the
ground truth; it re-derives every quantity independently:

GATE 1  MME direction from an explicit-variance-component reference:
        R = Ve*I, CoVar(a) = H*Va, so the MME is
        [X'R^-1 X  X'R^-1 Z;  Z'R^-1 X  Z'R^-1 Z + H^-1/Va][b;a] = [X'R^-1y; Z'R^-1y]
        solved canonically (numpy, full system, intercept kept). No lambda
        shortcut and no h2 anywhere. Two direction-distinguishing VCE pairs:
        (Va=0.4, Ve=0.6) -> lambda = Ve/Va = 1.5  (h2 = 0.4)
        (Va=0.8, Ve=0.2) -> lambda = Ve/Va = 0.25 (h2 = 0.8)
        If the C++ lambda direction were wrong the EBV would disagree here;
        h2 = 0.5 would NOT distinguish the direction and is deliberately avoided.

GATE 2  tau/omega are tuning weights, not heritability fractions: Hinv =
        Ainv + tau*Ginv - omega*A22^-1. Runs at tau=1.5 and omega=1.2 must
        match explicit references built with the SAME tau/omega, and must
        differ from the tau=omega=1 run (sensitivity), proving they are
        independent of lambda.

GATE 3  Canonical (intercept kept) vs mu-absorbed solvers across several
        datasets: the C++ LDLT path eliminates the intercept column leftmost,
        which reproduces the mu-absorbed system exactly. This gate checks that
        (a) canonical and mu-absorbed solutions agree, and (b) diagnoses the
        rank/conditioning structure (Hinv*1, min singular value of the
        mu-absorbed matrix) to decide whether a free gauge actually exists.

GATE 4  Back-solving reconstruction: a_g = G q, effects alpha = (M' q)/S.
        (i)    M @ EFFECT_ALLELE must reconstruct EBV(typed) (dosage scale)
        (ii)   Z @ SNP_EFFECT_STD must reconstruct EBV(typed) (standardized
               scale) if the current SNP_EFFECT_STD column is reconstruction-
               consistent; a systematic scale mismatch is reported.
        Reported for both: Pearson, Spearman, RMSE, max_abs, slope, intercept.
        (iii)  OLS sanity: regress EBV(typed) on M; EFFECT_ALLELE should align
               with per-allele dosage regression coefficients.

GATE 5  SNP variance scale: Var(M_j * EA_j) over animals vs 2p(1-p) EA_j^2
        must agree up to the n/(n-1) sample-variance factor (slope check).

GATE 6  Window-variance definitions under LD. Definition A (meta formula):
        sum_j 2p_j(1-p_j) alpha_j^2 over the window. Definition B:
        Var(M_win alpha_win) (genomic variance realized by the window). On a
        dedicated high-LD dataset (r ~ 0.9) the two are very different; the
        gate verifies the C++ WINDOW_VARIANCE column equals Definition A and
        quantifies B/A.

GATE 7  PCT denominator and metadata semantics: PCT_GENETIC_VARIANCE_j =
        SNP_VARIANCE_j / TOTAL_SNP_VARIANCE is a per-marker marginal share,
        not a fraction of additive genetic variance; asserted numerically.

No implementation file (src/) is modified by this audit.
"""
import argparse
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_ssgwas_backsolve_analytic import (
    build_A, build_Ainv, make_pedigree, ridge_inv, write_plink,
)

REL = 1e-4
ABSF = 1e-8


def relok(d, ref):
    return d <= REL * max(1.0, float(np.abs(ref).max())) + ABSF


def run_binary(bin_path, pref, tag_dir, vc, tau=1.0, omega=1.0, solver="ldlt",
               window=5):
    os.makedirs(tag_dir, exist_ok=True)
    va, ve = vc
    cmd = [bin_path, "--bfile", pref, "--pheno", pref + ".pheno", "--gwas",
           "--single-step", "--single-step-method", "backsolve",
           "--ped", pref + ".ped",
           "--reml-vg", str(va), "--reml-ve", str(ve),
           "--h2-tau", str(tau), "--h2-omega", str(omega),
           "--ssgwas-solver", solver,
           "--no-loco", "--ssgwas-window-snps", str(window),
           "--out", os.path.join(tag_dir, "m2")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise SystemExit("FAIL: cosmicscanner returned %d" % r.returncode)
    return tag_dir


def load_ebv(rd):
    ebv = {}
    with open(os.path.join(rd, "m2.ss_backsolve.ebv.txt")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            ebv[p[0]] = float(p[1])
    return ebv


def load_snp(rd):
    cols = {}
    with open(os.path.join(rd, "m2.ss_backsolve.snp.txt")) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            cols[d["SNP"]] = d
    return cols


def load_win(rd):
    rows = []
    with open(os.path.join(rd, "m2.ss_backsolve.window.txt")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append(p)
    return rows


def load_meta(rd):
    meta = {}
    with open(os.path.join(rd, "m2.ss_backsolve.meta.txt")) as f:
        for l in f.read().splitlines():
            if "\t" in l:
                k, v = l.split("\t", 1)
                meta[k] = v
    return meta


def explicit_sigma_mme(ped, idx, typed, Ainv, A22_inv, Ginv, pheno_ids,
                       y_vals, tau, omega, va, ve):
    """Full canonical MME with intercept, built ONLY from variance components:
    R = Ve I, CoVar(a) = H Va. No lambda, no h2 shortcut.
    Returns (mu, a_full, a_mu, Hinv)."""
    n_ped = len(ped)
    Hinv = Ainv.copy()
    for gi, i in enumerate(typed):
        for gj, j in enumerate(typed):
            Hinv[i, j] += tau * Ginv[gi, gj] - omega * A22_inv[gi, gj]
    Hinv = 0.5 * (Hinv + Hinv.T)

    y_full = np.full(n_ped, np.nan)
    for a in pheno_ids:
        y_full[idx[a]] = y_vals[a]
    rec = [i for i in range(n_ped) if math.isfinite(y_full[i])]
    n_r = len(rec)
    y = y_full[rec]
    Zd = np.zeros((n_r, n_ped))
    for k, i in enumerate(rec):
        Zd[k, i] = 1.0

    rinv = 1.0 / ve
    Z1 = Zd.sum(axis=0)
    C = np.zeros((1 + n_ped, 1 + n_ped))
    C[0:1, 0:1] = rinv * n_r
    C[0:1, 1:] = rinv * Z1
    C[1:, 0:1] = rinv * Z1.reshape(-1, 1)
    C[1:, 1:] = rinv * (Zd.T @ Zd) + (1.0 / va) * Hinv
    rhs = np.concatenate([rinv * np.array([y.sum()]), rinv * (Zd.T @ y)])
    sol = np.linalg.solve(C, rhs)
    mu, a_full = sol[0], sol[1:]

    Z1 = Zd.sum(axis=0)
    lam = ve / va
    C_mu = (Zd.T @ Zd + lam * Hinv) - np.outer(Z1, Z1) / n_r
    rhs_mu = Zd.T @ y - Z1 * y.mean()
    a_mu = np.linalg.solve(C_mu, rhs_mu)
    return mu, a_full, a_mu, Hinv


def stats_report(tag, x, y, name):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    pear = float(np.corrcoef(x, y)[0, 1])
    sp = 1.0 - 6.0 * float(np.sum((rx - ry) ** 2)) / (len(x) * (len(x) ** 2 - 1))
    rmse = float(np.sqrt(np.mean((x - y) ** 2)))
    mabs = float(np.abs(x - y).max())
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    print(f"  [{tag}] {name}: pearson={pear:.6f} spearman={sp:.6f} "
          f"rmse={rmse:.3e} max_abs={mabs:.3e} slope={slope:.6f} "
          f"intercept={intercept:.6f}")
    return dict(pearson=pear, spearman=sp, rmse=rmse, max_abs=mabs,
                slope=slope, intercept=intercept)


def gate3_check(ped, idx, typed, pheno_ids, y_vals, Ainv, A22_inv, Ginv,
                tau, omega, va, ve, print_res=True):
    mu, a_full, a_mu, Hinv = explicit_sigma_mme(
        ped, idx, typed, Ainv, A22_inv, Ginv, pheno_ids, y_vals, tau, omega,
        va, ve)
    d = float(np.abs(a_full - a_mu).max())
    h1 = float(np.abs(Hinv @ np.ones(len(ped))).max())
    lam = ve / va
    n_ped = len(ped)
    rec = [i for i in range(n_ped) if math.isfinite(y_vals.get(ped[i][0],
                                                               np.nan))]
    n_r = len(rec)
    Zd = np.zeros((n_r, n_ped))
    for k, i in enumerate(rec):
        Zd[k, i] = 1.0
    W = Zd.T @ Zd + lam * Hinv
    Z1 = Zd.sum(axis=0)
    C_mu = W - np.outer(Z1, Z1) / n_r
    s = np.linalg.svd(C_mu, compute_uv=False)
    if print_res:
        print(f"  [gate3] canonical vs mu-absorbed max|diff|={d:.3e}  "
              f"max|Hinv*1|={h1:.3e}  min_sv(C_mu)={s[-1]:.3e}  "
              f"cond(C_mu)={s[0]/s[-1]:.3e}")
    return dict(diff=d, h1=h1, min_sv=s[-1], ratio=s[0] / s[-1])


def zquant(p):
    """Standard normal quantile by bisection (no scipy dependency)."""
    lo, hi = -8.0, 8.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def correlated_snps(rng, n, m, p, corr):
    """m dosages in {0,1,2} with MAF p and pairwise latent correlation ~corr."""
    latent = rng.normal(0.0, 1.0, n)
    z0 = zquant((1.0 - p) ** 2)
    z1 = zquant(1.0 - p * p)
    X = np.zeros((n, m), dtype=int)
    for j in range(m):
        u = corr * latent + math.sqrt(1.0 - corr * corr) * rng.normal(0.0, 1.0, n)
        x = np.zeros(n)
        x[u > z0] = 1.0
        x[u > z1] = 2.0
        X[:, j] = x
    return X


def defA_windows(chrom, v, W):
    """Definition A per sliding window (same algorithm as the C++ output)."""
    out = []
    s = 0
    m = len(v)
    while s < m:
        ch = chrom[s]
        e = s
        while e + 1 < m and chrom[e + 1] == ch:
            e += 1
        L = e - s + 1
        if L <= W:
            out.append((ch, float(v[s:e + 1].sum())))
        else:
            for k in range(s, e - W + 2):
                out.append((ch, float(v[k:k + W].sum())))
        s = e + 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--gates", default="1,2,3,4,5,6,7",
                    help="comma list of gates to run")
    args = ap.parse_args()
    gates = {int(g) for g in args.gates.split(",")}
    rng = np.random.default_rng(args.seed)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    results = {}

    ped = make_pedigree()
    idx = {p[0]: i for i, p in enumerate(ped)}
    n_ped = len(ped)
    A = build_A(ped, idx)
    Ainv = build_Ainv(ped, idx, A)
    typed = [idx[f"g2_{i:02d}"] for i in range(24)] + \
            [idx[f"f{i:02d}"] for i in range(12)] + \
            [idx[f"g3_{i:02d}"] for i in range(8)]
    typed_names = [ped[i][0] for i in typed]
    pheno_ids = [f"f{i:02d}" for i in range(12)] + \
                [f"g2_{i:02d}" for i in range(28)] + \
                [f"g3_{i:02d}" for i in range(14)]
    m1, m2 = 50, 40
    m = m1 + m2
    chrom = ["1"] * m1 + ["2"] * m2
    snp_ids = [f"rs{j:04d}" for j in range(m)]
    bp = [1000 * (j + 1) for j in range(m)]
    X = np.zeros((len(typed), m), dtype=int)
    for j in range(m):
        p = rng.uniform(0.15, 0.5)
        X[:, j] = rng.binomial(2, p, size=len(typed))
    n_g = len(typed)

    pref = os.path.join(outdir, "sim")
    write_plink(pref, typed_names, X, chrom, snp_ids, bp)
    a_true = rng.normal(0.0, 1.0, n_ped)
    y_vals = {}
    for a in pheno_ids:
        y = 5.0 + a_true[idx[a]] + rng.normal(0.0, math.sqrt(1.5))
        y_vals[a] = round(float(y), 6)
    with open(pref + ".pheno", "w") as f:
        for a in pheno_ids:
            f.write(f"0 {a} {y_vals[a]}\n")
    with open(pref + ".ped", "w") as f:
        for a, s, d in ped:
            f.write(f"{a} {s} {d}\n")

    means = X.mean(axis=0)
    s2pq = 2.0 * (means / 2.0) * (1.0 - means / 2.0)
    S = s2pq.sum()
    M = X.astype(float) - means
    G = (M @ M.T) / S
    BLEND_ALPHA = 0.05
    Gb = (1.0 - BLEND_ALPHA) * G + BLEND_ALPHA * np.eye(n_g)
    GinvC = ridge_inv(Gb, "G")
    A22 = A[np.ix_(typed, typed)]
    A22_inv = ridge_inv(A22, "A22")
    lam_default = 2.0 / 1.0
    print(f"data: n_ped={n_ped} n_typed={n_g} n_pheno={len(pheno_ids)} "
          f"S={S:.6f}")

    # ---------------- GATE 1 + GATE 4 + GATE 5 + GATE 7 ----------------
    vcs = [(0.4, 0.6), (0.8, 0.2)]
    for (va, ve) in vcs:
        lam = ve / va
        tag = f"va{va}_ve{ve}"
        if 1 in gates or 4 in gates or 5 in gates or 7 in gates:
            rd = run_binary(args.bin, pref, os.path.join(outdir, tag),
                            (va, ve))
            ebv = load_ebv(rd)
            snp = load_snp(rd)
            meta = load_meta(rd)
            _, a_r, a_mu_r, _ = explicit_sigma_mme(
                ped, idx, typed, Ainv, A22_inv, GinvC, pheno_ids, y_vals,
                tau=1.0, omega=1.0, va=va, ve=ve)
            a_c = np.array([ebv[a] for a, _, _ in ped], dtype=float)
            if 1 in gates:
                d_can = float(np.abs(a_c - a_r).max())
                d_mu = float(np.abs(a_c - a_mu_r).max())
                print(f"[GATE1] ({tag}, lambda=Ve/Va={lam:.6f}, h2={va/(va+ve):.4f}) "
                      f"explicit-sigma canonical max|diff|={d_can:.3e} "
                      f"mu-absorbed max|diff|={d_mu:.3e}")
                ok = relok(d_can, a_r) and relok(d_mu, a_mu_r)
                print(f"[GATE1] ({tag}) {'PASS' if ok else 'FAIL'}")
                results[f"gate1_{tag}"] = ok
                if not ok:
                    raise SystemExit("FAIL: GATE1 direction mismatch")
            if 4 in gates:
                a_g = a_c[typed]
                q = GinvC @ a_g
                EA = (M.T @ q) / S
                sd = np.sqrt(np.where(s2pq > 1e-9, s2pq, 1.0))
                Zs = M / sd
                rec_d = M @ EA
                ea_col = np.array([float(snp[s]["EFFECT_ALLELE"])
                                   for s in snp_ids], dtype=float)
                std_col = np.array([float(snp[s]["SNP_EFFECT_STD"])
                                    for s in snp_ids], dtype=float)
                stats_report(f"gate4.dosage.{tag}", rec_d, a_g,
                             "M@EFFECT_ALLELE vs GEBV(typed)")
                stats_report(f"gate4.std.{tag}", Zs @ std_col, a_g,
                             "Z@SNP_EFFECT_STD vs GEBV(typed)")
                ols_b = np.linalg.lstsq(np.hstack([M, np.ones((n_g, 1))]),
                                        a_g, rcond=None)[0][:m]
                rec_ols = M @ ols_b
                stats_report(f"gate4.ols.{tag}", rec_ols, a_g,
                             "M@ols_beta(EBV~M) vs GEBV(typed)")
                stats_report(f"gate4.olsalpha.{tag}", ols_b, ea_col,
                             "EA vs OLS beta")
                d_ea = float(np.abs(ea_col - EA).max())
                ok_rec = relok(d_ea, EA)
                std_ref = EA * sd
                d_std = float(np.abs(std_col - std_ref).max())
                print(f"[GATE4] ({tag}) EFFECT_ALLELE vs (M'q)/S "
                      f"max|diff|={d_ea:.3e} "
                      f"{'PASS' if ok_rec else 'FAIL'}")
                print(f"[GATE4] ({tag}) SNP_EFFECT_STD vs EA*sigma_j "
                      f"(reconstruction-consistent scale) max|diff|="
                      f"{d_std:.3e}")
                results[f"gate4_ea_{tag}"] = ok_rec
            if 5 in gates:
                sv_ref = s2pq * (np.array([float(snp[s]["EFFECT_ALLELE"])
                                           for s in snp_ids])) ** 2
                sv = np.array([float(snp[s]["SNP_VARIANCE"])
                               for s in snp_ids])
                d_sv = float(np.abs(sv - sv_ref).max())
                ok_sv = relok(d_sv, sv_ref)
                ea_col = np.array([float(snp[s]["EFFECT_ALLELE"])
                                   for s in snp_ids], dtype=float)
                obs_var = np.array([np.sum((M[:, j] * ea_col[j]) ** 2) /
                                    (n_g - 1) for j in range(m)])
                sel = sv_ref > 1e-12
                wmean = float(obs_var[sel].sum() / sv_ref[sel].sum())
                ok_n = abs(wmean - 1.0) <= 0.05
                print(f"[GATE5] ({tag}) SNP_VARIANCE vs 2pq*EA^2 "
                      f"max|diff|={d_sv:.3e} "
                      f"{'PASS' if ok_sv else 'FAIL'}")
                print(f"[GATE5] ({tag}) Var_s(M_j*EA_j) vs 2pq_j*EA_j^2 "
                      f"weighted-mean ratio={wmean:.6f} (expect 1.0, "
                      f"n={n_g}) {'PASS' if ok_n else 'FAIL'}")
                results[f"gate5_{tag}"] = ok_sv and ok_n
            if 7 in gates:
                sv = np.array([float(snp[s]["SNP_VARIANCE"])
                               for s in snp_ids])
                pct = np.array([float(snp[s]["PCT_GENETIC_VARIANCE"])
                                for s in snp_ids])
                tot = float(meta["TOTAL_SNP_VARIANCE"])
                ok_den = max(np.abs(pct - sv / tot).max(),
                             abs(tot - sv.sum())) <= \
                    REL * max(1.0, tot) + ABSF
                print(f"[GATE7] ({tag}) PCT=SNP_VAR/TOTAL "
                      f"{'PASS' if ok_den else 'FAIL'}  "
                      f"TOTAL_SNP_VARIANCE={tot:.6e}")
                results[f"gate7_{tag}"] = ok_den
            if tag == "va0.4_ve0.6":
                print("[GATE1] meta keys: "
                      f"BACKSOLVE_FORMULA={meta.get('BACKSOLVE_FORMULA')!r} "
                      f"S_SNP={meta.get('S_SNP')!r} TAU={meta.get('TAU')!r} "
                      f"OMEGA={meta.get('OMEGA')!r} "
                      f"VCE_MODE={meta.get('VCE_MODE')!r}")

    # ---------------- GATE 2 -----------------
    if 2 in gates:
        tau_rd = run_binary(args.bin, pref, os.path.join(outdir, "tau1p5"),
                            (0.4, 0.6), tau=1.5)
        om_rd = run_binary(args.bin, pref, os.path.join(outdir, "om1p2"),
                           (0.4, 0.6), omega=1.2)
        base_rd = os.path.join(outdir, "va0.4_ve0.6")
        _, a_ref_tau, _, _ = explicit_sigma_mme(
            ped, idx, typed, Ainv, A22_inv, GinvC, pheno_ids, y_vals,
            tau=1.5, omega=1.0, va=0.4, ve=0.6)
        _, a_ref_om, _, _ = explicit_sigma_mme(
            ped, idx, typed, Ainv, A22_inv, GinvC, pheno_ids, y_vals,
            tau=1.0, omega=1.2, va=0.4, ve=0.6)
        ebv_tau = load_ebv(tau_rd)
        ebv_om = load_ebv(om_rd)
        ebv_base = load_ebv(base_rd)
        a_tau = np.array([ebv_tau[a] for a, _, _ in ped])
        a_om = np.array([ebv_om[a] for a, _, _ in ped])
        a_base = np.array([ebv_base[a] for a, _, _ in ped])
        d_tau = float(np.abs(a_tau - a_ref_tau).max())
        d_om = float(np.abs(a_om - a_ref_om).max())
        d_sens_t = float(np.abs(a_tau - a_base).max())
        d_sens_o = float(np.abs(a_om - a_base).max())
        ok = relok(d_tau, a_ref_tau) and relok(d_om, a_ref_om) and \
            d_sens_t > 1e-3 and d_sens_o > 1e-3
        print(f"[GATE2] tau=1.5 ref max|diff|={d_tau:.3e}  "
              f"omega=1.2 ref max|diff|={d_om:.3e}  "
              f"sensitivity |tau-run - base|={d_sens_t:.3e}  "
              f"|omega-run - base|={d_sens_o:.3e}")
        print(f"[GATE2] {'PASS' if ok else 'FAIL'} "
              f"(tau/omega are independent tuning weights, not h2)")
        results["gate2"] = ok
        if not ok:
            raise SystemExit("FAIL: GATE2")

    # ---------------- GATE 3 -----------------
    if 3 in gates:
        print("[GATE3] canonical vs mu-absorbed across random datasets "
              "(independent numpy solve)")
        worst = 0.0
        for s in range(5):
            r3 = np.random.default_rng(1000 + s)
            p3 = make_pedigree()
            i3 = {p[0]: i for i, p in enumerate(p3)}
            A3 = build_A(p3, i3)
            Ainv3 = build_Ainv(p3, i3, A3)
            t3 = [i3[f"g2_{i:02d}"] for i in range(24)] + \
                 [i3[f"f{i:02d}"] for i in range(12)] + \
                 [i3[f"g3_{i:02d}"] for i in range(8)]
            ph3 = [f"f{i:02d}" for i in range(12)] + \
                  [f"g2_{i:02d}" for i in range(28)] + \
                  [f"g3_{i:02d}" for i in range(14)]
            X3 = r3.binomial(2, r3.uniform(0.15, 0.5, 30),
                             size=(len(t3), 30))
            y3 = {a: float(np.round(r3.normal(0.0, 1.0), 6)) for a in ph3}
            n3 = len(t3)
            mean3 = X3.mean(axis=0)
            s3 = 2.0 * (mean3 / 2.0) * (1.0 - mean3 / 2.0)
            M3 = X3.astype(float) - mean3
            G3 = (M3 @ M3.T) / s3.sum()
            G3b = 0.95 * G3 + 0.05 * np.eye(n3)
            Ginv3 = ridge_inv(G3b, "G")
            A22i3 = ridge_inv(A3[np.ix_(t3, t3)], "A22")
            res = gate3_check(p3, i3, t3, ph3, y3, Ainv3, A22i3, Ginv3,
                              tau=1.0, omega=1.0, va=0.4, ve=0.6)
            worst = max(worst, res["diff"])
        ok = worst <= 1e-7
        print(f"[GATE3] worst canonical-vs-mu-absorbed diff={worst:.3e}  "
              f"{'PASS' if ok else 'FAIL (check gauge story)'}")
        results["gate3"] = ok
        if not ok:
            raise SystemExit("FAIL: GATE3")

    # ---------------- GATE 6 -----------------
    if 6 in gates:
        r6 = np.random.default_rng(5000)
        n6, m6 = 200, 10
        p6 = 0.3
        X6 = correlated_snps(r6, n6, m6, p6, 0.95)
        mean6 = X6.mean(axis=0)
        s6 = 2.0 * (mean6 / 2.0) * (1.0 - mean6 / 2.0)
        M6 = X6.astype(float) - mean6
        alpha6 = np.ones(m6)
        a6 = M6 @ alpha6 + r6.normal(0.0, 0.1, n6)
        EA6 = np.linalg.lstsq(M6, a6, rcond=None)[0]
        A_def = float((s6 * EA6 ** 2).sum())
        B_def = float(np.var(M6 @ EA6))
        print(f"[GATE6] LD example (latent corr=0.95, m=10, MAF=0.3): "
              f"Definition A (sum 2pq*alpha^2)={A_def:.4f}  "
              f"Definition B (Var(M_win*alpha_win))={B_def:.4f}  "
              f"B/A={B_def / A_def:.3f}  (Definitions A and B are NOT "
              f"interchangeable under LD)")
        rd = run_binary(args.bin, pref, os.path.join(outdir, "g6win"),
                        (0.4, 0.6))
        win = load_win(rd)
        snp6 = load_snp(rd)
        ea6 = np.array([float(snp6[s]["EFFECT_ALLELE"]) for s in snp_ids])
        v6 = s2pq * ea6 ** 2
        win_ref = defA_windows(chrom, v6, 5)
        win_cpp = np.array([float(r[6]) for r in win])
        if len(win_cpp) == len(win_ref):
            d_w = float(np.abs(win_cpp - np.array([w[1] for w in win_ref]))
                        .max())
            ok = d_w <= REL * max(1.0, float(abs(win_ref[0][1]))) + ABSF
        else:
            d_w = float("inf")
            ok = False
        tot_win = float(win_cpp.sum())
        tot_snp = float(v6.sum())
        print(f"[GATE6] C++ WINDOW_VARIANCE per window vs Definition A "
              f"(n_windows={len(win_cpp)} vs {len(win_ref)}) "
              f"max|diff|={d_w:.3e}  sums {tot_win:.6e} vs {tot_snp:.6e} "
              f"(sliding overlap: sum != TOTAL by construction)  "
              f"{'PASS' if ok else 'FAIL'}")
        meta6 = load_meta(rd)
        print(f"[GATE6] meta WINDOW_SIZE={meta6.get('WINDOW_SIZE')!r}")
        results["gate6"] = ok

    print("\n--- GATE SUMMARY ---")
    for k in sorted(results):
        print(f"  {k}: {'PASS' if results[k] else 'FAIL'}")
    if not all(results.values()):
        raise SystemExit("FAIL: one or more audit gates failed")
    print("PASS: Phase 6.5 mathematical red-team audit gates "
          f"({','.join(sorted(str(g) for g in gates))})")
    return 0


if __name__ == "__main__":
    sys.exit(main())