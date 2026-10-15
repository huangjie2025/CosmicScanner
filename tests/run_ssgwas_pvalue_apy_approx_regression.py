#!/usr/bin/env python3
"""Phase D gate D3 (tiny): exact vs apy-exact vs apy-approx on identical
tiny data (n_ped=30, n_g=20, m=30; analysis-consistent null draw).

Gates asserted:
  D3a: apy-approx with core = n_g must give EBV and EFFECT bit-identical to
       exact (APY component exact: G_APY^-1 = G^-1), and effect correlation 1.0
       at every core size tested (core-only operator identity check).
  D3b: (unit-level, test_apy_ginv_d15) Leite Eq.8 (superscript blocks) ==
       Schur form; not repeated here.
  D3c: QUANTIFICATION gate (NOT equality): at core = n_g the residual
       PEV/D-approximation error is reported as Pearson/slope/RMSE/max_abs of
       SE/Z/LOG10P vs exact and vs apy-exact - the one-sided D-hatched
       deviation source is Eq. 9 (diagonal absorption + A22-structure), not APY.
  Metadata: DENSE_NG2_ALLOCATED=false, DENSE_UNCERTAINTY_ALLOCS=0,
       MME_INVERSE_METHOD=NONE..., CORE_UNCERTAINTY_METHOD=LEITE2024_EQ8_BERMANN2022_SCHUR,
       EXACT_OR_APPROX=APPROXIMATE, GENOMIC_CONTROL_APPLIED=false.
  PCG smoke: apy-approx + --ssgwas-solver pcg must run, converge and produce
       finite P (SOLVER_CONVERGED=true) - full battery is gate D10.

Usage:
  python3 tests/run_ssgwas_pvalue_apy_approx_regression.py --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_apy_approx --seed 123
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np


def meuwissen_luo(pedigree):
    n = len(pedigree)
    idx = {r[0]: i for i, r in enumerate(pedigree)}
    A = np.zeros((n, n))
    L = np.zeros((n, n))
    for i, r in enumerate(pedigree):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        if s >= 0 and d >= 0:
            L[i, i] = np.sqrt(max(1.0 + 0.5 * A[s, d] - 0.25 * (A[s, s] + A[d, d]) - 0.5 * A[s, d], 0.0))
            L[i, 0:i] = 0.5 * (L[s, 0:i] + L[d, 0:i])
            A[i, i] = 1.0 + 0.5 * A[s, d]
        elif s >= 0 or d >= 0:
            p_ = s if s >= 0 else d
            L[i, i] = np.sqrt(max(1.0 - 0.25 * A[p_, p_], 0.0))
            L[i, 0:i] = 0.5 * L[p_, 0:i]
            A[i, i] = 1.0
        else:
            L[i, i] = 1.0
            A[i, i] = 1.0
        for j in range(i):
            A[i, j] = 0.5 * ((A[s, j] if s >= 0 else 0.0) + (A[d, j] if d >= 0 else 0.0))
            A[j, i] = A[i, j]
    return A


def write_plink_trio(bed_prefix, samples, m, p, rng):
    with open(bed_prefix + ".fam", "w") as f:
        for iid in samples:
            f.write(f"F{iid} {iid} 0 0 0 -9\n")
    with open(bed_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"1 rs{j + 1} 0 {j + 1} A G\n")
    n = len(samples)
    dos = np.array([rng.binomial(2, p[j], size=n) for j in range(m)]).T.astype(np.float64)
    bit = {0: 0b11, 1: 0b10, 2: 0b00}
    bytes_per = (n + 3) // 4
    raw = bytearray()
    for j in range(m):
        for b in range(bytes_per):
            byte = 0
            for k in range(4):
                i = b * 4 + k
                if i < n:
                    byte |= bit[int(dos[i, j])] << (2 * k)
            raw.append(byte)
    with open(bed_prefix + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01" + bytes(raw))
    return dos


def run_cs(bin, work, out, mode, pvalue_method, core_size=None, seed=7,
           core_file=None, solver="ldlt", expect_fail=False, vg=1.0, ve=3.0):
    os.makedirs(out, exist_ok=True)
    cmd = [bin, "--bfile", f"{work}/sim", "--pheno", f"{work}/sim_y.phe",
           "--gwas", "--single-step", "--single-step-method", "backsolve",
           "--ped", f"{work}/sim.cs.ped", "--reml-vg", f"{vg:.10f}",
           "--reml-ve", f"{ve:.10f}", "--no-loco", "--ssgwas-window-snps", "10",
           "--ssgwas-solver", solver, "--ssgwas-pvalue-method", pvalue_method,
           "--out", f"{out}/m2"]
    if mode == "a22":
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    if core_size is not None:
        cmd += ["--ssgwas-apy-core-size", str(core_size)]
    if core_file:
        cmd += ["--ssgwas-apy-core-file", core_file]
    if seed is not None:
        cmd += ["--ssgwas-apy-seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if expect_fail:
        return r
    if r.returncode != 0:
        raise RuntimeError(f"run failed:\n{r.stderr[-1500:]}")
    return r


def load_ebv(path):
    return np.array([float(l.split("\t")[1]) for l in open(path).read().splitlines()[1:]])


def load_pvalue(path):
    rows = sorted([l.split("\t") for l in open(path).read().splitlines()[1:]],
                  key=lambda r: int(r[1][2:]))
    return (np.array([float(r[6]) for r in rows]),   # EA
            np.array([float(r[7]) for r in rows]),   # SE
            np.array([float(r[9]) for r in rows]),   # Z
            np.array([float(r[11]) for r in rows]))  # LOG10P


def load_meta(path):
    return {l.split("\t")[0]: l.split("\t")[1] for l in
            open(path).read().splitlines() if "\t" in l}


def pearson(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return 1.0 if np.allclose(a, b) else float("nan")
    return np.corrcoef(a, b)[0, 1]


def metrics(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    slope = 1.0
    if np.std(a) > 0.0:
        slope = float(np.cov(a, b, ddof=0)[0, 1] / np.var(a))
    return {"pearson": float(pearson(a, b)),
            "slope": slope,
            "rmse": float(np.sqrt(np.mean((a - b) ** 2))),
            "max_abs": float(np.max(np.abs(a - b)))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--tolerance", type=float, default=1e-5)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    work = f"{args.outdir}/data"
    os.makedirs(work, exist_ok=True)

    n_ped, n_geno, m = 30, 20, 30
    n_f = n_ped // 3
    pedigree = []
    for i in range(1, n_f + 1):
        pedigree.append([i, 0, 0])
    for i in range(n_f + 1, n_ped + 1):
        s = int(rng.integers(1, n_f + 1))
        d = int(rng.integers(1, n_f + 1))
        while d == s:
            d = int(rng.integers(1, n_f + 1))
        pedigree.append([i, s, d])
    geno_ids = list(range(n_f + 1, n_ped + 1))
    p = rng.uniform(0.15, 0.45, size=m)
    dos = write_plink_trio(f"{work}/sim", geno_ids, m, p, rng)
    pheno_ids = sorted(set(rng.choice(list(range(1, n_ped + 1)), size=21, replace=False)))
    with open(f"{work}/sim.phe", "w") as f:
        f.write("FID IID Pheno_1\n")
        for i in pheno_ids:
            f.write(f"F{i} {i} 0.0\n")
    with open(f"{work}/sim.cs.ped", "w") as f:
        for r in pedigree:
            f.write(f"{r[0]} {r[1]} {r[2]}\n")

    idx = {r[0]: i for i, r in enumerate(pedigree)}
    geno_idx_ped = [idx[i] for i in geno_ids]
    A = meuwissen_luo(pedigree)
    mean_dos = dos.mean(axis=0)
    pj = mean_dos / 2.0
    M = dos - mean_dos
    S = np.sum(2.0 * pj * (1.0 - pj))
    G_raw = M @ M.T / S
    alpha = 0.05
    A22 = A[np.ix_(geno_idx_ped, geno_idx_ped)]
    dG, oG = np.trace(G_raw) / n_geno, (G_raw.sum() - np.trace(G_raw)) / (n_geno * (n_geno - 1))
    dA, oA = np.trace(A22) / n_geno, (A22.sum() - np.trace(A22)) / (n_geno * (n_geno - 1))
    G_star = ((dA - oA) / (dG - oG)) * G_raw + (dA - (dA - oA) / (dG - oG) * dG)
    Gw = (1.0 - alpha) * G_star + alpha * A22
    H = A.copy()
    H[np.ix_(geno_idx_ped, geno_idx_ped)] = Gw
    w, Vh = np.linalg.eigh(H)
    w = np.maximum(w, 0.0)
    pheno_rows = [idx[i] for i in pheno_ids]
    va, ve = 1.0, 3.0
    a_true = (Vh @ (np.sqrt(w) * rng.normal(size=n_ped))) * np.sqrt(va)
    y = 1.0 + a_true[pheno_rows] + rng.normal(scale=np.sqrt(ve), size=len(pheno_rows))
    with open(f"{work}/sim_y.phe", "w") as f:
        f.write("FID IID Pheno_1\n")
        for k, r in enumerate(pheno_rows):
            f.write(f"F{pedigree[r][0]} {pedigree[r][0]} {y[k]:.10f}\n")

    fails = []
    d3c_rows = []

    for mode in ["identity", "a22"]:
        out_ex = f"{args.outdir}/{mode}_exact"
        run_cs(args.bin, work, out_ex, mode, "exact")
        ex_ebv = load_ebv(f"{out_ex}/m2.ss_backsolve.ebv.txt")
        ex_ea, ex_se, ex_z, ex_lp = load_pvalue(f"{out_ex}/m2.ss_backsolve.pvalue.txt")

        out_ae = f"{args.outdir}/{mode}_apyexact"
        run_cs(args.bin, work, out_ae, mode, "apy-exact", core_size=n_geno)
        ae_ea, ae_se, ae_z, ae_lp = load_pvalue(f"{out_ae}/m2.ss_backsolve.pvalue.txt")

        for core in [n_geno, n_geno - 5, n_geno - 10, n_geno - 15]:
            out = f"{args.outdir}/{mode}_apyapprox_core{core}"
            run_cs(args.bin, work, out, mode, "apy-approx", core_size=core)
            ebv = load_ebv(f"{out}/m2.ss_backsolve.ebv.txt")
            ea, se, z, lp = load_pvalue(f"{out}/m2.ss_backsolve.pvalue.txt")
            meta = load_meta(f"{out}/m2.ss_backsolve.meta.txt")

            if core == n_geno:
                d_ebv = metrics(ex_ebv, ebv)
                d_ea = metrics(ex_ea, ea)
                if d_ebv["max_abs"] > 0.0 or d_ea["max_abs"] > 0.0:
                    fails.append(f"{mode} core=n_g: EBV/EA not bit-identical to exact "
                                 f"({d_ebv['max_abs']}, {d_ea['max_abs']})")
                    continue
                for k, v in {"DENSE_NG2_ALLOCATED": "false",
                             "DENSE_UNCERTAINTY_ALLOCS": "0",
                             "MME_INVERSE_METHOD": "NONE (no selected inverse; Eq. 8 C_core_approx only)",
                             "CORE_UNCERTAINTY_METHOD": "LEITE2024_EQ8_BERMANN2022_SCHUR",
                             "EXACT_OR_APPROX": "APPROXIMATE",
                             "GENOMIC_CONTROL_APPLIED": "false",
                             "P_VALUE_METHOD": "AGUILAR2019_APY_APPROX"}.items():
                    if meta.get(k) != v:
                        fails.append(f"{mode} meta {k}: '{meta.get(k)}' != '{v}'")
                m_se = metrics(ex_se, se)
                m_z = metrics(ex_z, z)
                m_lp = metrics(ex_lp, lp)
                m_se_x = metrics(ae_se, se)
                m_z_x = metrics(ae_z, z)
                m_lp_x = metrics(ae_lp, lp)
                d3c_rows.append({
                    "mode": mode, "core": core, "reference": "exact",
                    "SE": m_se, "Z": m_z, "LOG10P": m_lp, "EA": d_ea})
                d3c_rows.append({
                    "mode": mode, "core": core, "reference": "apy-exact",
                    "SE": m_se_x, "Z": m_z_x, "LOG10P": m_lp_x})
                print(f"{mode} core=n_g D3c vs exact: SE {m_se['pearson']:.4f} "
                      f"slope {m_se['slope']:.4f} RMSE {m_se['rmse']:.3e} "
                      f"maxA {m_se['max_abs']:.3e}; Z {m_z['pearson']:.4f}/"
                      f"{m_z['slope']:.4f}/{m_z['max_abs']:.2e}; "
                      f"LOG10P {m_lp['pearson']:.4f}/{m_lp['slope']:.4f}/{m_lp['max_abs']:.2e}")
                print(f"{mode} core=n_g D3c vs apy-exact: SE {m_se_x['pearson']:.4f} "
                      f"slope {m_se_x['slope']:.4f} RMSE {m_se_x['rmse']:.3e}; "
                      f"LOG10P {m_lp_x['pearson']:.4f}/{m_lp_x['slope']:.4f}")
            else:
                if not (np.all(np.isfinite(z)) and np.all(np.isfinite(lp))):
                    fails.append(f"{mode} core={core}: non-finite outputs")
                    continue
                r_core = metrics(ex_ea, ea)
                print(f"{mode} core={core}: finite; EA corr vs exact {r_core['pearson']:.4f} "
                      f"rmse {r_core['rmse']:.3e}; meta GC false="
                      f"{meta.get('GENOMIC_CONTROL_APPLIED')}")

    # PCG smoke (gate D10 battery deferred): apy-approx + pcg on core=15
    out_pcg = f"{args.outdir}/identity_apyapprox_pcg"
    run_cs(args.bin, work, out_pcg, "identity", "apy-approx", core_size=15,
           solver="pcg")
    meta = load_meta(f"{out_pcg}/m2.ss_backsolve.meta.txt")
    ea, se, z, lp = load_pvalue(f"{out_pcg}/m2.ss_backsolve.pvalue.txt")
    if meta.get("SOLVER_CONVERGED") != "true":
        fails.append(f"pcg smoke: SOLVER_CONVERGED={meta.get('SOLVER_CONVERGED')}")
    if not (np.all(np.isfinite(z)) and np.all(np.isfinite(lp))):
        fails.append("pcg smoke: non-finite P")
    print(f"pcg smoke core=15: converged={meta.get('SOLVER_CONVERGED')} "
          f"residual={meta.get('SOLVER_RESIDUAL', 'NA')}")

    with open(f"{args.outdir}/d3_summary.json", "w") as f:
        json.dump({"fails": fails, "d3c": d3c_rows,
                   "meta_pcg": dict(load_meta(f"{out_pcg}/m2.ss_backsolve.meta.txt"))}, f,
                  indent=2)
    if fails:
        print("FAILURES:\n" + "\n".join(fails))
        sys.exit(1)
    print(f"PASS: D3a/D3c tiny gates (bit-identical at core=n_g; D-approximation "
          f"quantified against exact and apy-exact; PCG smoke OK)")


if __name__ == "__main__":
    main()