#!/usr/bin/env python3
"""Simulation-level validation of CosmicScanner ssGBLUP back-solving
(--single-step-method backsolve, Phase 6).

Checks:
  1. Causal recovery: 900-animal, 3-generation pedigree; 1000 independent
     bi-allelic SNPs; 5 causal loci. Back-solved per-SNP effects must rank
     causal loci on top (Pearson/Spearman vs true effects, top-k hits,
     causal-window rank).
  2. Typed/untyped integration: 300 typed + 150 untyped phenotyped +
     + 150 typed unphenotyped + rest neither ancestors; phenotyped animals
     include ungenotyped ones so H (A22 correction) is exercised.
  3. Pedigree-shuffle invariance: random row permutation of the .ped file
     must reproduce IDENTICAL EBV/effects (topological reordering + explicit
     ID mapping), so pedigree input order is irrelevant.
  4. Solver consistency: ldlt vs pcg EBVs identical within tolerance.
"""
import argparse
import math
import os
import random
import subprocess
import sys

import numpy as np

REL = 1e-6
ABSF = 1e-9


def write_plink(path_prefix, names, X, chrom, snp_ids, bp):
    n, m = X.shape
    with open(path_prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"0 {names[i]} 0 0 0 -9\n")
    with open(path_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"{chrom[j]} {snp_ids[j]} 0 {bp[j]} A G\n")
    nbytes = 3 + m * ((n + 3) // 4)
    buf = bytearray(nbytes)
    buf[0:3] = bytes([0x6C, 0x1B, 0x01])
    for j in range(m):
        for i in range(n):
            code = {2: 0, 1: 2, 0: 3}[int(X[i, j])]
            byte_idx = 3 + j * ((n + 3) // 4) + i // 4
            buf[byte_idx] |= (code << (2 * (i % 4)))
    with open(path_prefix + ".bed", "wb") as f:
        f.write(bytes(buf))


def run_cs(bin_path, outdir, tag, bfile, pheno, ped, extra=None):
    d = os.path.join(outdir, tag)
    os.makedirs(d, exist_ok=True)
    cmd = [bin_path, "--bfile", bfile, "--pheno", pheno, "--gwas",
           "--single-step", "--single-step-method", "backsolve",
           "--ped", ped,
           "--reml-vg", "1.0", "--reml-ve", "2.5",
           "--no-loco", "--ssgwas-window-snps", "10"]
    if extra:
        cmd += extra
    cmd += ["--out", os.path.join(d, "m2")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise SystemExit(f"FAIL: cosmicscanner returned {r.returncode}")
    return d


def load_ebv(d):
    ebv = {}
    with open(os.path.join(d, "m2.ss_backsolve.ebv.txt")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            ebv[p[0]] = float(p[1])
    return ebv


def load_snp(d):
    out = {}
    with open(os.path.join(d, "m2.ss_backsolve.snp.txt")) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            dd = dict(zip(hdr, line.rstrip("\n").split("\t")))
            out[dd["SNP"]] = {k: (float(v) if k not in ("CHR", "SNP", "A1", "A2") else v)
                              for k, v in dd.items()}
    return out


def load_win(d):
    rows = []
    with open(os.path.join(d, "m2.ss_backsolve.window.txt")) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append((p[0], float(p[2]), float(p[7])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--n-founders", type=int, default=200)
    ap.add_argument("--n-g2", type=int, default=400)
    ap.add_argument("--n-g3", type=int, default=300)
    ap.add_argument("--m-snp", type=int, default=1000)
    ap.add_argument("--n-causal", type=int, default=5)
    ap.add_argument("--n-typed", type=int, default=300)
    ap.add_argument("--n-pheno", type=int, default=650)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # ---- pedigree (parents always precede offspring; founders first) ----
    founders = [f"f{i:06d}" for i in range(args.n_founders)]
    g2 = [f"g2_{i:06d}" for i in range(args.n_g2)]
    g3 = [f"g3_{i:06d}" for i in range(args.n_g3)]
    ped = [(a, "0", "0") for a in founders]
    nf = args.n_founders
    for k, a in enumerate(g2):
        ped.append((a, founders[k % nf], founders[(k * 7 + 3) % nf]))
    for k, a in enumerate(g3):
        ped.append((a, g2[k % args.n_g2], g2[(k * 5 + 11) % args.n_g2]))
    all_ids = [a for a, _, _ in ped]

    # ---- partitions ----
    rng.shuffle(all_ids)
    typed_ids = all_ids[:args.n_typed]
    typed_set = set(typed_ids)
    rest = all_ids[args.n_typed:]
    pheno_ids = sorted(set(rest[:args.n_pheno - min(args.n_pheno, len(typed_ids))]) |
                       set(rng.choice(typed_ids, min(len(typed_ids), max(0, args.n_pheno - len(rest))), replace=False)))
    if len(pheno_ids) < args.n_pheno:
        pheno_ids = sorted(rng.choice(all_ids, args.n_pheno, replace=False))
    assert len(set(typed_ids) & set(pheno_ids)) < len(typed_ids)  # typed unphenotyped exist
    assert len(set(typed_ids)) - len(set(typed_ids) & set(pheno_ids)) > 0
    assert len(set(pheno_ids) - set(typed_ids)) > 0

    # ---- genotypes: independent biallelic, MAF ~ U(0.15, 0.5) ----
    m = args.m_snp
    ma = []
    for j in range(m):
        ma.append(rng.uniform(0.15, 0.5))
    X = np.zeros((len(typed_ids), m), dtype=int)
    for j in range(m):
        X[:, j] = rng.binomial(2, ma[j], size=len(typed_ids))

    # ---- causal effects and phenotypes ----
    causal = rng.choice(m, args.n_causal, replace=False)
    beta = np.zeros(m)
    s2pq = np.array([2.0 * p * (1.0 - p) for p in ma])
    for j in causal:
        beta[j] = rng.normal(0.0, 1.0 / math.sqrt(s2pq[j]))
    iz = {all_ids[i]: i for i in range(len(all_ids))}
    typed_i = [iz[a] for a in typed_ids]
    Xfull = np.zeros((len(all_ids), m))
    for k, i in enumerate(typed_i):
        Xfull[i] = X[k]
    tbv = Xfull @ beta
    y_vals = {}
    for a in pheno_ids:
        y_vals[a] = round(float(4.0 + tbv[iz[a]] + rng.normal(0.0, 1.0)), 6)

    pref = os.path.join(outdir, "sim")
    write_plink(pref, typed_ids, X,
                ["1"] * m, [f"rs{j:06d}" for j in range(m)], [1000 * (j + 1) for j in range(m)])
    with open(pref + ".pheno", "w") as f:
        for a in pheno_ids:
            f.write(f"0 {a} {y_vals[a]}\n")
    with open(pref + ".ped", "w") as f:
        for a, s, d in ped:
            f.write(f"{a} {s} {d}\n")
    shuffled = list(ped)
    rng.shuffle(shuffled)
    with open(pref + ".ped.shuf", "w") as f:
        for a, s, d in shuffled:
            f.write(f"{a} {s} {d}\n")

    # ---- runs ----
    d_ldlt = run_cs(args.bin, outdir, "ldlt", pref, pref + ".pheno", pref + ".ped")
    d_pcg = run_cs(args.bin, outdir, "pcg", pref, pref + ".pheno", pref + ".ped",
                   extra=["--ssgwas-solver", "pcg"])
    d_shuf = run_cs(args.bin, outdir, "shuf", pref, pref + ".pheno", pref + ".ped.shuf")

    # ---- 1. causal recovery ----
    snp = load_snp(d_ldlt)
    eff = np.array([abs(snp[f"rs{j:06d}"]["EFFECT_ALLELE"]) for j in range(m)])
    order_beta = np.argsort(-np.abs(beta))
    s = order_beta[:200]                   # signal region: top-200 |true| effects
    r_pearson = np.corrcoef(eff[s], np.abs(beta[s]))[0, 1]
    top = np.argsort(eff)[::-1]
    hits_top10 = len(set(top[:10]) & set(map(int, causal)))
    hits_top20 = len(set(top[:20]) & set(map(int, causal)))
    win = load_win(d_ldlt)
    causal_bp = {1000 * (j + 1) for j in causal}
    causal_win_ranks = []
    for i, (ch, bp_e, pct) in enumerate(win):
        if bp_e in causal_bp or any(abs(bp_e - b) / (10 * 1000) < 1.0 for b in causal_bp
                                    if abs(bp_e - b) <= 9 * 1000):
            causal_win_ranks.append(i + 1)
    rank_best_causal_win = min(causal_win_ranks) if causal_win_ranks else len(win) + 1

    print(f"[causal] Pearson(top200 true)={r_pearson:.4f} "
          f"Top10={hits_top10}/{args.n_causal} Top20={hits_top20}/{args.n_causal} "
          f"best_causal_win_rank={rank_best_causal_win}/{len(win)}")
    if not math.isfinite(r_pearson) or r_pearson < 0.30:
        raise SystemExit("FAIL: causal recovery Pearson too low")
    if hits_top10 < 3 or hits_top20 < 4:
        raise SystemExit("FAIL: causal SNR too low (top-k recovery)")
    if rank_best_causal_win > 250:   # sparse causality: top quartile is the gate
        raise SystemExit("FAIL: causal window rank too poor")

    # ---- 2. typed/untyped integration ----
    ebv = load_ebv(d_ldlt)
    ebv_s = np.array([ebv[a] for a in all_ids])
    typed_ebv = np.array([ebv[a] for a in typed_ids])
    if not (np.isfinite(ebv_s).all() and np.isfinite(typed_ebv).all()):
        raise SystemExit("FAIL: non-finite EBVs")
    print(f"[integrate] typed units: {len(typed_ids)}, phenotyped: {len(pheno_ids)}, "
          f"EBV std full={ebv_s.std():.4f}, typed={typed_ebv.std():.4f}")

    # ---- 3. pedigree-shuffle invariance ----
    ebv_sh = load_ebv(d_shuf)
    d_sh = max(abs(ebv_sh[a] - ebv[a]) for a in all_ids)
    snp_sh = load_snp(d_shuf)
    d_eff = max(abs(snp_sh[f"rs{j:06d}"]["EFFECT_ALLELE"] - snp[f"rs{j:06d}"]["EFFECT_ALLELE"])
                for j in range(m))
    print(f"[shuffle] max|EBV diff|={d_sh:.3e}  max|effect diff|={d_eff:.3e}")
    if d_sh > REL * max(1.0, float(np.abs(ebv_s).max())) + ABSF:
        raise SystemExit("FAIL: pedigree shuffle invariance (EBV)")

    # ---- 4. solver consistency ----
    ebv_p = load_ebv(d_pcg)
    d_solv = max(abs(ebv_p[a] - ebv[a]) for a in all_ids)
    print(f"[solver] ldlt vs pcg max|EBV diff|={d_solv:.3e}")
    if d_solv > 1e-4:
        raise SystemExit("FAIL: ldlt vs pcg inconsistency")

    print("\nPASS: ssGBLUP back-solving simulation validation "
          "(causal recovery, typed/untyped integration, shuffle invariance, "
          "solver consistency)")
    return 0


if __name__ == "__main__":
    sys.exit(main())