#!/usr/bin/env python3
# Phase 4 / Method A record-level RRM-GWAS simulation (tests).
#
# Generates an unbalanced longitudinal PLINK dataset whose causal SNP effect is
# a continuous trajectory beta(t) = phi(t)'gamma (record-level Method A target),
# runs Cosmicscanner --rrm-unbalance --rrm-gwas-method record, and validates:
#   - calibration of the rank-aware joint P under the null (scenario A)
#   - truth recovery / joint significance (scenarios B..G)
#   - trajectory recovery: reconstructed beta(t) vs the injected trajectory
#   - no silent drops: every scanned SNP carries RANK/COND/STATUS diagnostics
#
# Data model (per record r of subject s(r) at time t_r):
#   y_r = b0 + phi(t_r)' u_{s(r)}   <- random genetic trajectory
#         + phi(t_r)' p_{s(r)}      <- permanent-environment trajectory
#         + x_{s(r),j} * phi(t_r)' gamma_j  <- SNP fixed trajectory (Method A)
#         + e_r
# with u ~ N(0, Ku (x) I), p ~ N(0, Kp (x) I), e ~ N(0, sigma2e).
import argparse
import math
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def rng_seed(scenario, base=20260815):
    return base + sum(ord(c) for c in scenario)


def draw_genotypes(n, m):
    rng = np.random.default_rng(0)
    g = np.empty((n, m), dtype=np.int8)
    for j in range(m):
        maf = 0.15 + 0.35 * ((j * 37) % 101) / 100.0
        p0 = (1 - maf) ** 2
        p1 = 2 * maf * (1 - maf)
        u = rng.random(n)
        g[:, j] = np.where(u < p0, 0, np.where(u < p1, 1, 2)).astype(np.int8)
    return g


def write_plink(prefix, g, n, m):
    with (prefix.with_suffix(".fam")).open("w", encoding="utf-8") as f:
        for i in range(n):
            f.write(f"F{i+1:04d}\tI{i+1:04d}\t0\t0\t0\t-9\n")
    with (prefix.with_suffix(".bim")).open("w", encoding="utf-8") as f:
        for j in range(m):
            chrom = min(3, 1 + j // (m // 3))
            snp = f"rs_causal" if j == 0 else f"rs_bg_{j+1:04d}"
            f.write(f"{chrom}\t{snp}\t0\t{1000 + j * 10}\tA\tG\n")
    bed = bytearray(b"\x6c\x1b\x01")
    code = {0: 0b00, 1: 0b10, 2: 0b11}
    bpb = (n + 3) // 4
    for j in range(m):
        block = bytearray(bpb)
        for i in range(n):
            block[i // 4] |= code[int(g[i, j])] << ((i % 4) * 2)
        bed.extend(block)
    prefix.with_suffix(".bed").write_bytes(bed)


def scenario_spec(scenario, seed):
    """Return (n, times_pool, missing_rate, n_records_range, gamma_true, extra_notes)."""
    r = random.Random(seed)
    if scenario == "A":    # null
        n, times, miss, rec_range, gamma = 250, [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0], 0.2, (5, 8), None
    elif scenario == "B":  # constant
        n, times, miss, rec_range, gamma = 250, [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0], 0.2, (5, 8), [0.8, 0.0, 0.0, 0.0]
    elif scenario == "C":  # increasing
        n, times, miss, rec_range, gamma = 250, [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0], 0.2, (5, 8), [0.25, 0.55, 0.0, 0.0]
    elif scenario == "D":  # crossing (sign-switching linear)
        n, times, miss, rec_range, gamma = 250, [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0], 0.2, (5, 8), [0.30, -0.60, 0.0, 0.0]
    elif scenario == "E":  # curved
        n, times, miss, rec_range, gamma = 250, [0.0, 0.15, 0.3, 0.45, 0.6, 0.8, 1.0], 0.2, (5, 8), [0.2, 0.1, 0.5, 0.0]
    elif scenario == "F":  # severe unbalance
        n, times, miss, rec_range, gamma = 300, [0.0, 0.12, 0.25, 0.4, 0.55, 0.7, 0.9], 0.5, (2, 6), [0.7, 0.0, 0.0, 0.0]
    else:                  # G: rank-deficient (few subjects, few distinct times)
        n, times, miss, rec_range, gamma = 60, [0.0, 0.5, 1.0], 0.3, (2, 5), [0.9, 0.0, 0.0, 0.0]
    if gamma is not None and len(gamma) < 4:
        gamma = gamma + [0.0] * (4 - len(gamma))
    # Scenarios where causal-signal/trajectory recovery is NOT an objective:
    # A=null calibration, G=diagnostics-only rank-deficiency exercise.
    diagnostics_only = scenario in {"A", "G"}
    return n, times, miss, rec_range, gamma, diagnostics_only


def legendre_phi(t, order=3):
    # Same as buildLegendreMatrix: shifted Legendre on [-1,1], normalized sqrt((2k+1)/2).
    tn = 2.0 * t - 1.0
    phi = np.zeros(order + 1)
    phi[0] = 1.0
    phi[1] = tn
    for k in range(2, order + 1):
        phi[k] = ((2 * k - 1) * tn * phi[k - 1] - (k - 1) * phi[k - 2]) / k
    for k in range(order + 1):
        phi[k] *= math.sqrt((2 * k + 1) / 2.0)
    return phi


def generate_dataset(data_dir, scenario):
    seed = rng_seed(scenario)
    n, times_pool, miss, rec_range, gamma, diagnostics_only = scenario_spec(scenario, seed)
    m = 300
    order = 3
    prefix = data_dir / scenario
    g = draw_genotypes(n, m)
    write_plink(prefix, g, n, m)

    r = random.Random(seed)
    g_causal = g[:, 0].astype(float)
    g_causal = g_causal - g_causal.mean()

    Ku = np.array([[1.5, 0.10, -0.02, 0.001], [0.10, 0.35, 0.02, -0.001],
                   [-0.02, 0.02, 0.08, 0.0], [0.001, -0.001, 0.0, 0.02]])
    Kp = 0.5 * np.array([[1.0, 0.05, 0.0, 0.0], [0.05, 0.4, 0.01, 0.0],
                         [0.0, 0.01, 0.1, 0.0], [0.0, 0.0, 0.0, 0.05]])
    sigma2e = 0.35

    rows = []
    times = np.array(times_pool, dtype=float)
    for i in range(n):
        n_rec = r.randint(rec_range[0], rec_range[1])
        t_sub = [r.choice(times) for _ in range(n_rec)]
        rngi = np.random.default_rng(seed + i)
        Lk = np.linalg.cholesky(Ku + 1e-6 * np.eye(4))
        Lp = np.linalg.cholesky(Kp + 1e-8 * np.eye(4))
        u = Lk @ rngi.normal(size=4)
        p = Lp @ rngi.normal(size=4)
        for t in t_sub:
            if r.random() < miss:
                continue
            phi = legendre_phi(t, order)
            mu_t = 1.0 + 0.3 * (2.0 * t - 1.0)
            snp_part = 0.0 if gamma is None else g_causal[i] * float(np.dot(gamma, phi))
            y = mu_t + float(np.dot(phi, u)) + float(np.dot(phi, p)) + snp_part + r.gauss(0, math.sqrt(sigma2e))
            rows.append((f"F{i+1:04d}", f"I{i+1:04d}", t, y))
    pheno = pd.DataFrame(rows, columns=["FID", "IID", "TIME", "Y"])
    out = data_dir / f"{scenario}.long.pheno"
    pheno.to_csv(out, sep="\t", index=False, header=["FID", "IID", "TIME", "Y"])
    np.save(data_dir / f"{scenario}.gamma.npy", np.asarray(gamma, dtype=float) if gamma is not None else None)
    with (data_dir / f"{scenario}.meta.txt").open("w", encoding="utf-8") as f:
        f.write(f"scenario={scenario}\nN={n}\nM={m}\nTIME_POOL={','.join(map(str, times_pool))}\n")
        f.write(f"MISS={miss}\nREC_RANGE={rec_range[0]}-{rec_range[1]}\nRECORDS={len(rows)}\n")
        f.write(f"GAMMA={','.join(map(str, gamma)) if gamma is not None else 'None'}\n")
    return pheno


def run_cosmic(bin_path, data_dir, scenario, result_dir, extra=()):
    cmd = [
        str(bin_path),
        "--bfile", str(data_dir / scenario),
        "--pheno", str(data_dir / f"{scenario}.long.pheno"),
        "--tpheno", "1",
        "--mpheno", "2",
        "--rrm-unbalance",
        "--rrm-gwas-method", "record",
        "--rrm-order", "3",
        "--rrm-ai-iter", "8",
        "--rrm-curve-grid", "7",
        "--rrm-top", "10",
        "--rrm-sig-threshold", "0.05",
        "--maf", "0",
        "--gwas",
        "--no-loco",
        "--out", str(result_dir / scenario),
    ] + list(extra)
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return result_dir / f"{scenario}.rrm_unbal.rec.assoc.txt"


def validate(bin_path, data_dir, scenario, result_dir, diagnostics_only):
    out = []

    # Calibration under the null (scenario A): joint P of background SNPs.
    if scenario == "A":
        assoc = pd.read_csv(result_dir / "A.rrm_unbal.rec.assoc.txt", sep="\t")
        p = assoc["P_JOINT"].astype(float).dropna()
        fp05 = float((p < 0.05).mean())
        lam = float(np.median(assoc["STAT_JOINT"].astype(float).dropna()) / np.median(np.random.chisquare(4, 200000)))
        out.append(("N_SNP", len(assoc)))
        out.append(("FPR_0.05", fp05))
        out.append(("INFLATION", lam))
        ok = (len(assoc) > 0 and 0.005 <= fp05 <= 0.15 and 0.5 <= lam <= 1.8
              and (assoc["STATUS"] == "OK").mean() > 0.99)
        out.append(("Pass", bool(ok)))
        return dict(out)

    assoc = pd.read_csv(result_dir / f"{scenario}.rrm_unbal.rec.assoc.txt", sep="\t")
    gamma = np.load(data_dir / f"{scenario}.gamma.npy")
    meta = {}
    with (data_dir / f"{scenario}.meta.txt").open("r", encoding="utf-8") as f:
        for line in f:
            k, _, v = line.strip().partition("=")
            meta[k] = v
    times_pool = [float(x) for x in meta["TIME_POOL"].split(",")]

    # 1. Causal SNP recovered and joint-significant.
    causal = assoc[assoc["SNP"] == "rs_causal"]
    rank = None
    if len(causal) != 1:
        out.append(("causal_found", False))
    else:
        out.append(("causal_found", True))
        rank = int((assoc["P_JOINT"].astype(float) < causal.iloc[0]["P_JOINT"]).sum()) + 1
        out.append(("causal_rank", rank))
        out.append(("causal_joint_p", float(causal.iloc[0]["P_JOINT"])))
        out.append(("causal_status", causal.iloc[0]["STATUS"]))

    # 2. Trajectory recovery for the causal SNP: reconstructed beta(t) vs the
    #    injected trajectory. A global sign flip is allowed (allele-coding
    #    convention is arbitrary), so use |Pearson| for variable trajectories.
    #    For a constant trajectory (scenario B) Pearson is undefined; check that
    #    the reconstruction is flat at the correct magnitude instead.
    gamma_flat = np.asarray(gamma, dtype=float)
    cb = float("nan")
    flat = float("nan")
    amp = float("nan")
    true_scale = float("nan")
    if len(causal) == 1 and gamma_flat is not None:
        df = pd.read_csv(result_dir / f"{scenario}.rrm_unbal.rec.dynamic_curve.txt", sep="\t")
        cur = df[df["SNP"] == "rs_causal"]
        if len(cur) > 1:
            times = list(dict.fromkeys(cur["Time"].astype(float)))
            beta_hat = np.asarray([cur[cur["Time"] == t]["Beta"].iloc[0] for t in times])
            beta_true = np.asarray([float(np.dot(gamma_flat, legendre_phi(t, 3))) for t in times])
            if np.std(beta_true) > 1e-6:
                cb = float(abs(np.corrcoef(beta_hat, beta_true)[0, 1]))
            else:
                amp = float(np.max(np.abs(beta_hat)))
                flat = (float(np.max(beta_hat)) - float(np.min(beta_hat))) / max(amp, 1e-9)
                true_scale = float(abs(beta_true[0]))
    out.append(("trajectory_corr", cb))
    out.append(("trajectory_flat", flat))

    # 3. Diagnostics: no silent drops, all P finite.
    n_snp = len(assoc)
    n_finite = int(assoc["P_JOINT"].notna().sum())
    statuses = assoc["STATUS"].value_counts().to_dict()
    out.append(("N_SNP", n_snp))
    out.append(("N_FINITE_P", n_finite))
    out.append(("STATUS_counts", statuses))
    no_numeric_failure = statuses.get("NUMERICAL_FAILURE", 0) == 0
    ok_frac = statuses.get("OK", 0) / max(1, n_snp)
    if diagnostics_only:
        # The rank-deficient scenario intentionally produces RANK_DEFICIENT /
        # ILL_CONDITIONED statuses; the requirement is that every SNP is still
        # reported with a finite P and no silent NUMERICAL_FAILURE.
        rank_diag_ok = bool(n_finite == n_snp and no_numeric_failure and n_snp > 0)
        out.append(("Pass", rank_diag_ok))
        return dict(out)
    rank_ok = bool(ok_frac > 0.9 and n_finite == n_snp and n_snp > 0)

    pass_causal = bool(
        len(causal) == 1
        and causal.iloc[0]["P_JOINT"] < 1e-3
        and rank is not None
        and rank <= 5
    )
    if cb == cb:
        pass_traj = bool(cb > 0.5)
    elif flat == flat:
        pass_traj = bool(flat <= 0.6 and 0.1 * true_scale <= amp <= 3.0 * true_scale)
    else:
        pass_traj = False
    pass_diag = rank_ok
    ok = bool(pass_causal and pass_traj and pass_diag)
    out.append(("Pass", ok))
    return dict(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="/disk228/huangjie/PanOmics/CosmicScanner/build/cosmicscanner")
    ap.add_argument("--data-dir", default="/tmp/opencode/rrm_rec_method_a_data")
    ap.add_argument("--result-dir", default="/tmp/opencode/rrm_rec_method_a_results")
    ap.add_argument("--scenarios", default="A,G", help="comma-separated A-G or 'all'")
    args = ap.parse_args()

    scenarios = list("ABCDEFG") if args.scenarios == "all" else [s.strip() for s in args.scenarios.split(",")]
    for s in scenarios:
        if s not in "ABCDEFG":
            raise SystemExit(f"Unknown scenario {s}")

    data_dir = Path(args.data_dir)
    result_dir = Path(args.result_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    bin_path = Path(args.bin)
    if not bin_path.exists():
        raise SystemExit(f"Missing binary {bin_path}")

    summary_rows = []
    for scenario in scenarios:
        print(f"=== Scenario {scenario} ===")
        seed = rng_seed(scenario)
        _, _, _, _, _, diagnostics_only = scenario_spec(scenario, seed)
        generate_dataset(data_dir, scenario)
        assoc = run_cosmic(bin_path, data_dir, scenario, result_dir)
        print(f"  scanned {len(pd.read_csv(assoc, sep='\t'))} SNPs")
        res = validate(bin_path, data_dir, scenario, result_dir, diagnostics_only)
        summary_rows.append({"Scenario": scenario, **{k: str(v) for k, v in res.items()}})
        print(f"  {res}")

    pd.DataFrame(summary_rows).to_csv(result_dir / "method_a_simulation_summary.tsv",
                                      sep="\t", index=False)
    print(f"\nSummary: {result_dir / 'method_a_simulation_summary.tsv'}")
    if any(r["Pass"] != "True" for r in summary_rows):
        raise SystemExit("Method A simulation: one or more scenarios failed")


if __name__ == "__main__":
    main()
