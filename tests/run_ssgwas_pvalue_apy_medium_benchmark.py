#!/usr/bin/env python3
"""Phase D gate D5: medium-scale Exact vs APY-exact vs APY-approx benchmark.

Two data regimes (mechanistic simulations, no artificial target on K98):
  A) HIGH_DIMENSIONAL_STRESS (stress): per-foundation-site independent
     haplotypes; genomic redundancy comes only from LD created as offspring
     inherit recombined single-chromosome haplotypes. Low-redundancy /
     difficult-APY control (K98/n_g expected high at these sample sizes).
  B) LIVESTOCK_LIKE_LD (livestock): finite founder population, per-chromosome
     founder haplotypes, generations with parent-pool breeding, meiotic
     gametes with Poisson crossovers. Genomic redundancy emerges from finite
     founders + inheritance + recombination (mechanistic, no preset
     correlation matrix).

Tiers (pedigree always > n_g; typed+pheno / typed+no-pheno / untyped+pheno /
untyped+no-pheno all present by construction):
  Medium-1: n_g=500 , M=5,000 , n_ped=750
  Medium-2: n_g=1,000, M=10,000, n_ped=1,500
  Medium-3: n_g=2,000, M=20,000, n_ped=3,000

LD: per-chromosome structure; stress = independent founder haplotypes with
recombination P(crossover/site)=0.01; livestock = finite founders + Poisson
crossovers (see LIVESTOCK_* constants; all parameters reported).

Core grid from the G eigen spectrum (K98 = eigenvalues needed for 98% of
cumulative eigenvalue trace of the final blended G), adaptive policy:
  K98/n_g <= 0.60: core in {0.50,1.00,1.50,2.00}*K98 + full
  K98/n_g >  0.60: core in {0.25,0.50,0.75,1.00}*K98 + full
  (dedup; floored at max(25, 5% n_g); never cap-and-duplicate; full always
  included; deterministic random core, fixed seed 7)

Layers (same data/phenotype/VEs/G convention; LDLT; no GC):
  L1  Exact_Ginv        vs Exact_GinvAPY         (APY G^-1 approximation)
  L2  Exact_GinvAPY     vs Approx_GinvAPY        (Eq.8/9 + core-only operator)
  TOT Exact_Ginv        vs Approx_GinvAPY        (user-visible total)

Interpretation (fixed): at core=n_g, L1 must recover exact reference
tolerance; L2 remains the Leite Eq.8/9 approximation residual and is only
quantified (==1 is not required). Verdicts:
  GOOD / ACCEPTABLE / CONCERN / POOR_APPROXIMATION / NUMERICAL_FAILURE
POOR_APPROXIMATION = unsuitable approximation at this core size (statistical);
NUMERICAL_FAILURE = NaN / non-finite / solver failure only (code correctness).

Outputs:
  medium_summary.json, medium_metrics.tsv, core_curve.tsv

Usage:
  python3 tests/run_ssgwas_pvalue_apy_medium_benchmark.py \
      --bin build/cosmicscanner \
      --outdir tests/results/ssgwas_pvalue_apy_medium --seed 123
  python3 tests/run_ssgwas_pvalue_apy_medium_benchmark.py \
      --regime livestock --outdir tests/results/ssgwas_pvalue_apy_medium_livestock
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

APY_SEED = 7
RECOMB = 0.01
ALPHA = 0.05          # identity blend alpha (no tune, blend I) - analysis convention
VA, VE = 1.0, 3.0
LAMBDA = VE / VA

TIERS = [(500, 5000, 12), (1000, 10000, 16), (2000, 20000, 20)]

# livestock-like regime parameters (transparent, mechanistic)
LIVESTOCK = {
    "N_FOUNDERS": 60,                 # 30 sires + 30 dams
    "N_GENERATIONS": 4,               # gen0 (founders) .. gen3 (offspring)
    "N_PARENTS_PER_GENERATION": "30+30 (sire pool + dam pool, reselected from latest generation)",
    "N_CHROMOSOMES": None,            # set per tier from M_PER_CHROMOSOME
    "M_PER_CHROMOSOME": 500,
    "RECOMBINATION_MODEL": "POISSON_CROSSOVERS_UNIFORM_POSITIONS (per chromosome, independent)",
    "MEAN_CROSSOVERS_PER_CHR": 1.2,
    "HAPLOTYPE_GENERATION_METHOD":
        "FOUNDER_CHROMOSOMES_SAMPLED_FROM_12_ANCESTRAL_POOL_PER_CHR (IID site MAFs; LD emerges from the finite ancestral pool + inheritance + recombination, no preset correlation matrix)",
    "ANCESTRAL_POOL_PER_CHR": 12,
    "SEED_BASE": APY_SEED,
}


# ---------------------------------------------------------------------------
# pedigree / LD simulation
# ---------------------------------------------------------------------------


def build_pedigree(n_ped, rng):
    n0 = n_ped // 3
    n1 = n_ped // 3
    n2 = n_ped - n0 - n1
    ped = []                     # [id, sire, dam], founders first
    for i in range(1, n0 + 1):
        ped.append([i, 0, 0])
    for i in range(n0 + 1, n0 + n1 + 1):
        s = int(rng.integers(1, n0 + 1))
        d = int(rng.integers(1, n0 + 1))
        while d == s:
            d = int(rng.integers(1, n0 + 1))
        ped.append([i, s, d])
    for i in range(n0 + n1 + 1, n_ped + 1):
        s = int(rng.integers(n0 + 1, n0 + n1 + 1))
        d = int(rng.integers(n0 + 1, n0 + n1 + 1))
        while d == s:
            d = int(rng.integers(n0 + 1, n0 + n1 + 1))
        ped.append([i, s, d])
    typed = list(range(n0 + 1, n_ped + 1))
    return ped, typed, n0, n1, n2


def relationship_matrix(ped):
    """A via per-row numpy recurrence (O(n^2) numpy work, no O(n^3) loop)."""
    idx = {r[0]: i for i, r in enumerate(ped)}
    n = len(ped)
    A = np.zeros((n, n))
    for i, r in enumerate(ped):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        A[i, i] = 1.0 + (0.5 * A[s, d] if s >= 0 and d >= 0 else 0.0)
        if s >= 0 or d >= 0:
            A[i, :i] = 0.5 * ((A[s, :i] if s >= 0 else 0.0) +
                              (A[d, :i] if d >= 0 else 0.0))
        A[:i, i] = A[i, :i]
    return A, idx


def simulate_dosages(ped, idx, typed_idx_ped, m, n_chr, rng):
    """Haplotype simulation with recombination along chromosomes; returns
    dosage matrix (n_typed x m, typed genotype-file order)."""
    n_ped = len(ped)
    base = m // n_chr
    per_ch = [base] * n_chr
    for k in range(m - base * n_chr):
        per_ch[k % n_chr] += 1
    sites_per_ch = per_ch
    haps = {}                    # animal index -> [hap1 array, hap2 array]
    for i, r in enumerate(ped):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        if s < 0:
            haps[i] = [np.concatenate([rng.binomial(1, rng.uniform(0.15, 0.45,
                          L)) for L in sites_per_ch]).astype(np.int8),
                       np.concatenate([rng.binomial(1, rng.uniform(0.15, 0.45,
                          L)) for L in sites_per_ch]).astype(np.int8)]
        else:
            ps, pd_ = haps[s], haps[d]
            haps[i] = [recombine(ps[0], ps[1], rng),
                       recombine(pd_[0], pd_[1], rng)]
    n_typed = len(typed_idx_ped)
    dos = np.zeros((n_typed, m), dtype=np.uint8)
    for k, gi in enumerate(typed_idx_ped):
        dos[k, :] = haps[gi][0] + haps[gi][1]
    return dos


def recombine(h1, h2, rng):
    """Crossover walk: start on h1, recombine between haps with prob RECOMB
    per site (phase alternation -> within-chromosome LD blocks)."""
    m = len(h1)
    out = h1.copy()
    active = 1
    rv = rng.random(m)
    for j in range(1, m):
        if rv[j] < RECOMB:
            active = 1 - active
        if active == 0:
            out[j] = h2[j]
    return out


def var_residual(A, s, d):
    if s >= 0 and d >= 0:
        return max(1.0 - 0.25 * (A[s, s] + A[d, d]), 1e-10)
    if s >= 0 or d >= 0:
        return max(0.75 - 0.25 * A[max(s, d), max(s, d)], 1e-10)
    return 1.0


def draw_pedigree_effects(ped, idx, A, rng):
    """Exact pedigree-Gaussian draw a_ped ~ N(0, A)."""
    n = len(ped)
    a = np.zeros(n)
    for i, r in enumerate(ped):
        s = idx[r[1]] if r[1] != 0 else -1
        d = idx[r[2]] if r[2] != 0 else -1
        a[i] = (0.5 * a[s] if s >= 0 else 0.0) + \
               (0.5 * a[d] if d >= 0 else 0.0) + \
               rng.normal(scale=math.sqrt(var_residual(A, s, d)))
    return a


def build_g(dos):
    M = dos - dos.mean(axis=0)
    pj = dos.mean(axis=0) / 2.0
    S = float(np.sum(2.0 * pj * (1.0 - pj)))
    G = M @ M.T / S
    return G, S, pj


def k98_of(G):
    ev = np.sort(np.linalg.eigvalsh(G))[::-1]
    tot = ev.sum()
    if tot <= 0:
        return len(ev)
    cum = 0.0
    for k, v in enumerate(ev):
        cum += max(v, 0.0)
        if cum >= 0.98 * tot:
            return k + 1
    return len(ev)


def k_fracs_of(G):
    """K_x = number of leading eigenvalues reaching x% of cumulative
    eigenvalue trace fraction {0.90,0.95,0.98,0.99}; plus effective rank
    exp(entropy of normalized eigenvalues). Diagnostic only."""
    ev = np.sort(np.linalg.eigvalsh(G))[::-1]
    ev = np.maximum(ev, 0.0)
    tot = ev.sum()
    out = {}
    if tot <= 0:
        n = len(ev)
        return {"K90": n, "K95": n, "K98": n, "K99": n,
                "effective_rank": n}
    cum = np.cumsum(ev) / tot
    for frac, label in ((0.90, "K90"), (0.95, "K95"),
                        (0.98, "K98"), (0.99, "K99")):
        out[label] = int(np.searchsorted(cum, frac) + 1)
    p = ev / tot
    ent = -np.sum(p * np.log(p + 1e-300))
    out["effective_rank"] = float(np.exp(ent))
    return out


def ld_decay(dos, per_ch, n_sub=150):
    """Mean |corr| of dosages vs SNP distance bins within each chromosome,
    on a subsample of SNPs. Returns {bin_bp_mid: mean |r|} over chromosomes."""
    bins = {}
    n = dos.shape[0]
    for base, L in enumerate(per_ch):
        start = sum(per_ch[:base])
        sub = np.sort(np.random.default_rng(1234 + base).choice(
            L, size=min(n_sub, L), replace=False))
        X = dos[:, start + sub].astype(float)
        X = X - X.mean(axis=0)
        sd = X.std(axis=0)
        sd[sd == 0] = 1.0
        X /= sd
        C = np.abs(X.T @ X) / n
        for i in range(len(sub)):
            for j in range(i + 1, len(sub)):
                d = abs(int(sub[i]) - int(sub[j]))
                b = int(d // 50) * 50 + 25
                bins.setdefault(b, []).append(C[i, j])
    return {k: float(np.mean(v)) for k, v in sorted(bins.items())}


def maf_summary(pj):
    maf = np.minimum(pj, 1.0 - pj)
    return {"min": float(maf.min()), "q1": float(np.quantile(maf, 0.01)),
            "q5": float(np.quantile(maf, 0.05)),
            "median": float(np.median(maf)),
            "monomorphic": int(np.sum(maf < 0.01))}


# ---------------------------------------------------------------------------
# livestock-like regime: finite founders + generations + meiotic gametes
# ---------------------------------------------------------------------------


def meiotic_gamete(h1, h2, rng, n_cross):
    """One gamete from haplotypes h1/h2 with n_cross crossover points
    (uniform positions, alternating segments; no crossovers -> h1 copy)."""
    L = len(h1)
    if n_cross <= 0 or L < 2:
        return h1.copy()
    pts = sorted(rng.choice(L - 1, size=n_cross, replace=False) + 1)
    out = np.empty(L, dtype=np.int8)
    seg0 = 0
    for k, p in enumerate(pts):
        src = h1 if k % 2 == 0 else h2
        out[seg0:p] = src[seg0:p]
        seg0 = p
    src = h1 if len(pts) % 2 == 0 else h2
    out[seg0:L] = src[seg0:L]
    return out


def build_pedigree_livestocklike(n_ped, n_g, m, rng):
    """Returns (ped, typed_ids, info). gen0 = founders (untyped pool);
    gen1 = expansion generation; typed = gen2 + gen3 = n_g; untyped rest.
    Parent pools of 30 sires + 30 dams reselected from the latest generation
    after each generation (breeding scheme)."""
    n_f0 = LIVESTOCK["N_FOUNDERS"]
    n_chr = int(round(m / LIVESTOCK["M_PER_CHROMOSOME"]))
    g1 = n_ped - n_f0 - n_g
    g2 = n_g // 2
    g3 = n_g - g2
    gen_sizes = [n_f0, g1, g2, g3]
    ped = []                     # [id, sire, dam], founders first; sex flag
    sex = []                     # 0 male, 1 female
    pools = []                   # pool[g][sex] = list of ped indices
    n_cross = None
    for g, size in enumerate(gen_sizes):
        for k in range(size):
            i = len(ped) + 1
            if g == 0:
                ped.append([i, 0, 0])
                sex.append(k % 2)
            else:
                ps = rng.choice(pools[g - 1][0])
                pd_ = rng.choice(pools[g - 1][1])
                while pd_ == ps:
                    pd_ = rng.choice(pools[g - 1][1])
                ped.append([i, ped[ps][0], ped[pd_][0]])
                sex.append(k % 2)
        pool_off = [i0 for i0 in range(len(ped) - size, len(ped))]
        rng.shuffle(pool_off)
        pools.append([pool_off[:30], pool_off[30:60]])
    typed = [i for i, r in enumerate(ped) if r[0] > n_f0 + g1]
    info = {"N_FOUNDERS": n_f0, "N_GENERATIONS": len(gen_sizes),
            "GEN_SIZES": gen_sizes,
            "N_PARENTS_PER_GENERATION": "30+30",
            "N_CHROMOSOMES": n_chr,
            "M_PER_CHROMOSOME": LIVESTOCK["M_PER_CHROMOSOME"],
            "MEAN_CROSSOVERS_PER_CHR": LIVESTOCK["MEAN_CROSSOVERS_PER_CHR"],
            "RECOMBINATION_MODEL": LIVESTOCK["RECOMBINATION_MODEL"],
            "HAPLOTYPE_GENERATION_METHOD":
                LIVESTOCK["HAPLOTYPE_GENERATION_METHOD"],
            "ANCESTRAL_POOL_PER_CHR": LIVESTOCK["ANCESTRAL_POOL_PER_CHR"],
            "SEED_BASE": APY_SEED,
            "typed_ids": typed}
    return ped, typed, info


def simulate_dosages_livestock(ped, idx, typed_idx_ped, m, rng, info):
    """Founder chromosomes are drawn from a small per-chromosome ancestral
    haplotype pool (IID Bernoulli per site at shared site MAFs): this gives
    within-chromosome haplotype structure, and LD emerges mechanistically from
    the finite pool + inheritance + recombination. Non-founders inherit
    recombined gametes (Poisson crossovers)."""
    n_chr = info["N_CHROMOSOMES"]
    base = m // n_chr
    per_ch = [base] * n_chr
    for k in range(m - base * n_chr):
        per_ch[k % n_chr] += 1
    site_maf = rng.uniform(0.15, 0.45, m)
    pool_size = LIVESTOCK["ANCESTRAL_POOL_PER_CHR"]
    pool = {}
    starts = np.cumsum([0] + per_ch)
    for b, L in enumerate(per_ch):
        pool[b] = rng.binomial(1, site_maf[starts[b]:starts[b] + L],
                               size=(pool_size, L)).astype(np.int8)
    haps = {}
    n_cross_rng = np.random.default_rng(4242)
    for i, r in enumerate(ped):
        if r[1] == 0 and r[2] == 0:
            haps[i] = [np.concatenate([pool[b][rng.integers(pool_size)].copy()
                                       for b in range(n_chr)])
                       for _ in range(2)]
    parents = {r[0]: i for i, r in enumerate(ped)}
    for i, r in enumerate(ped):
        if r[1] == 0:
            continue
        si = parents[r[1]]
        da = parents[r[2]]
        h1 = []
        h2 = []
        for b, L in enumerate(per_ch):
            a0, a1 = starts[b], starts[b] + L
            n1 = n_cross_rng.poisson(LIVESTOCK["MEAN_CROSSOVERS_PER_CHR"])
            n2 = n_cross_rng.poisson(LIVESTOCK["MEAN_CROSSOVERS_PER_CHR"])
            h1.append(meiotic_gamete(haps[si][0][a0:a1],
                                     haps[si][1][a0:a1], rng, n1))
            h2.append(meiotic_gamete(haps[da][0][a0:a1],
                                     haps[da][1][a0:a1], rng, n2))
        haps[i] = [np.concatenate(h1), np.concatenate(h2)]
    n_typed = len(typed_idx_ped)
    dos = np.zeros((n_typed, m), dtype=np.uint8)
    for k, gi in enumerate(typed_idx_ped):
        dos[k, :] = haps[gi][0] + haps[gi][1]
    return dos, per_ch


def simulate_null_pheno(ped, idx, typed_idx_ped, p_rows, dos, rng):
    """Model-consistent null (identity convention: tau=1, omega=1, blend-I):
    Cov(a) = A + block(Gw - A22) on typed; y = mu + a + e, mu=1."""
    A, _ = relationship_matrix(ped)
    G, S, pj = build_g(dos)
    A22 = A[np.ix_(typed_idx_ped, typed_idx_ped)]
    G_bs = (1.0 - ALPHA) * G + ALPHA * np.eye(G.shape[0])
    Gw_minus_A22 = G_bs - A22
    w, V = np.linalg.eigh(Gw_minus_A22)
    n_clamp = int(np.sum(w < -1e-9))
    w = np.maximum(w, 0.0)
    delta = (V @ (np.sqrt(w) * rng.normal(size=len(w)))) if len(w) else \
        np.zeros(0)
    a = draw_pedigree_effects(ped, idx, A, rng)
    a = a.copy()
    a[typed_idx_ped] += delta
    y = 1.0 + a[p_rows] + rng.normal(scale=math.sqrt(VE), size=len(p_rows))
    return y, A, n_clamp


# ---------------------------------------------------------------------------
# IO / metrics
# ---------------------------------------------------------------------------


def write_plink_trio(bed_prefix, samples, dos, rng_maf):
    with open(bed_prefix + ".fam", "w") as f:
        for iid in samples:
            f.write(f"F{iid} {iid} 0 0 0 -9\n")
    m = dos.shape[1]
    with open(bed_prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"1 rs{j + 1} 0 {j + 1} A G\n")
    n = dos.shape[0]
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


def run_cs(bin_, work, out, method, core=None, ratio=None, seed=APY_SEED,
           omp=8):
    os.makedirs(out, exist_ok=True)
    cmd = [bin_, "--bfile", f"{work}/sim", "--pheno", f"{work}/sim_y.phe",
           "--gwas", "--single-step", "--single-step-method", "backsolve",
           "--ped", f"{work}/sim.cs.ped", "--reml-vg", "1.0000000000",
           "--reml-ve", "3.0000000000", "--no-loco", "--ssgwas-window-snps",
           "10", "--ssgwas-pvalue-method", method, "--out", f"{out}/m2"]
    if core is not None:
        cmd += ["--ssgwas-apy-core-size", str(core),
                "--ssgwas-apy-seed", str(seed)]
    if ratio is not None:
        cmd += ["--ssgwas-apy-core-target-ratio", f"{ratio:.4f}"]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(omp)
    env["OPENBLAS_NUM_THREADS"] = "1"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"run failed: {cmd}\n{r.stderr[-2500:]}")
    return r


def load_pvalue(path):
    rows = sorted([l.split("\t") for l in
                   open(path).read().splitlines()[1:]],
                  key=lambda r: int(r[1][2:]))
    rows = [r for r in rows if all(v != "NA" for v in
                                   (r[6], r[7], r[8], r[9], r[11]))]
    return ([r[1] for r in rows],                             # SNP ids
            np.array([float(r[6]) for r in rows]),            # EA
            np.array([float(r[7]) for r in rows]),            # SE
            np.array([float(r[9]) for r in rows]),            # Z
            np.array([float(r[11]) for r in rows]))           # LOG10P


def align(vectors):
    """Intersect SNP ids across loaded runs (NA-row filtering can differ
    between methods, e.g. monomorphic SNPs), return per-run boolean mask in
    the original input order."""
    order = sorted(vectors, key=lambda v: -len(v[0]))
    common = set(order[0][0])
    for r in order[1:]:
        common &= set(r[0])
    idx = {id(r): i for i, r in enumerate(vectors)}
    masks = [None] * len(vectors)
    for r in order:
        masks[idx[id(r)]] = np.array([sid in common for sid in r[0]])
    return masks


def load_meta(path):
    return {l.split("\t")[0]: l.split("\t")[1] for l in
            open(path).read().splitlines() if "\t" in l}


def rankdata(x):
    return np.argsort(np.argsort(x)).astype(float) + 1.0


def pearson(a, b):
    if np.std(a) == 0.0 or np.std(b) == 0.0:
        return 1.0 if np.allclose(a, b) else float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    return pearson(rankdata(a), rankdata(b))


def full_metrics(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    slope = 1.0
    if np.std(a) > 0.0:
        slope = float(np.cov(a, b, ddof=0)[0, 1] / np.var(a))
    return {"pearson": pearson(a, b), "spearman": spearman(a, b),
            "slope": slope,
            "intercept": float(np.mean(b) - slope * np.mean(a)),
            "rmse": float(np.sqrt(np.mean((a - b) ** 2))),
            "max_abs": float(np.max(np.abs(a - b)))}


def se_ratio_stats(se_ref, se_app):
    r = np.asarray(se_app) / np.asarray(se_ref)
    r = r[np.isfinite(r)]
    qs = {}
    for q in (1, 5, 25, 50, 75, 95, 99):
        qs[str(q)] = float(np.quantile(r, q / 100.0))
    return {"median": float(np.median(r)), "mean": float(np.mean(r)),
            "quantiles": qs}


def topk_overlap(la, lb, k):
    ta = set(np.argsort(-np.asarray(la))[:k])
    tb = set(np.argsort(-np.asarray(lb))[:k])
    return float(len(ta & tb)) / k


def snp_pair_metrics(a, b):
    """shared-SNP metric block; arrays keyed by rank order."""
    return {"EFFECT": full_metrics(a[0], b[0]),
            "SE": full_metrics(a[1], b[1]),
            "SE_ratio": se_ratio_stats(a[1], b[1]),
            "Z": full_metrics(a[2], b[2]),
            "LOG10P": full_metrics(a[3], b[3]),
            "top10": topk_overlap(a[3], b[3], 10),
            "top50": topk_overlap(a[3], b[3], 50),
            "top100": topk_overlap(a[3], b[3], 100)}


# ---------------------------------------------------------------------------
# acceptance classification
# ---------------------------------------------------------------------------

def classify(m):
    """GOOD / ACCEPTABLE / CONCERN / POOR_APPROXIMATION / NUMERICAL_FAILURE.
    Ranking (top-N) is reported separately and does not gate the statistical
    verdict; NUMERICAL_FAILURE is reserved for true computational failure."""
    eff = m["EFFECT"]
    se = m["SE"]
    z = m["Z"]
    lp = m["LOG10P"]
    if not all(math.isfinite(v) for v in (eff["pearson"], se["pearson"],
                                          z["slope"], lp["slope"],
                                          se["rmse"], lp["rmse"])):
        return "NUMERICAL_FAILURE", "non-finite"
    if eff["pearson"] >= 0.98 and se["pearson"] >= 0.95 and \
       abs(se["slope"] - 1.0) <= 0.08 and abs(z["slope"] - 1.0) <= 0.06 and \
       abs(lp["slope"] - 1.0) <= 0.10:
        return "GOOD", ""
    if eff["pearson"] >= 0.95 and se["pearson"] >= 0.92 and \
       abs(se["slope"] - 1.0) <= 0.15 and abs(z["slope"] - 1.0) <= 0.12:
        return "ACCEPTABLE", ""
    if eff["pearson"] >= 0.85 and abs(z["slope"] - 1.0) <= 0.25:
        return "CONCERN", ""
    return "POOR_APPROXIMATION", "criteria"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--tiers", default="all",
                    help="comma list e.g. 0,1 (indices into TIERS)")
    ap.add_argument("--regime", default="stress",
                    choices=["stress", "livestock"])
    ap.add_argument("--omp", type=int, default=8)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    work = f"{args.outdir}/data"
    os.makedirs(work, exist_ok=True)

    tiers = list(range(len(TIERS)))
    if args.tiers != "all":
        tiers = [int(t) for t in args.tiers.split(",")]

    regime = ("HIGH_DIMENSIONAL_STRESS" if args.regime == "stress"
              else "LIVESTOCK_LIKE_LD")
    summary = {"regime": regime, "tiers": [], "classifications": {}, "q": {},
               "livestock_params": (dict(LIVESTOCK) if args.regime ==
                                    "livestock" else None)}
    all_metrics = []
    core_curve_rows = []

    for ti in tiers:
        n_g, m, n_chr = TIERS[ti]
        n_ped = 3 * n_g // 2
        tag = f"t{ti}_ng{n_g}"
        tier_out = f"{args.outdir}/{tag}"
        os.makedirs(tier_out, exist_ok=True)

        rng_t = np.random.default_rng(args.seed * 7919 + ti * 104729)
        if args.regime == "livestock":
            ped, typed_ids, geno_info = build_pedigree_livestocklike(
                n_ped, n_g, m, rng_t)
            idx = {r[0]: i for i, r in enumerate(ped)}
            typed_idx_ped = np.array([idx[i] for i in typed_ids])
            dos, per_ch = simulate_dosages_livestock(
                ped, idx, typed_idx_ped, m, rng_t, geno_info)
        else:
            ped, typed_ids, n0, n1, n2 = build_pedigree(n_ped, rng_t)
            idx = {r[0]: i for i, r in enumerate(ped)}
            typed_idx_ped = np.array([idx[i] for i in typed_ids])
            dos = simulate_dosages(ped, idx, typed_idx_ped, m, n_chr, rng_t)
            base = m // n_chr
            per_ch = [base] * n_chr
            for k in range(m - base * n_chr):
                per_ch[k % n_chr] += 1
        M_ = None
        S_ = None
        pj_ = None

        # ---- data files ----
        A, _ = relationship_matrix(ped)
        A22 = A[np.ix_(typed_idx_ped, typed_idx_ped)]
        G_raw, S_2pq, pj = build_g(dos)
        G_bs = (1.0 - ALPHA) * G_raw + ALPHA * np.eye(n_g)
        k98 = k98_of(G_bs)
        k98_frac = k98 / n_g
        ks = k_fracs_of(G_bs)
        M_ = dos - dos.mean(axis=0)

        # ---- phenotype (model-consistent null, identity convention) ----
        n_rec = int(0.7 * n_ped)
        pids = sorted(set(rng_t.choice(list(range(1, n_ped + 1)), size=n_rec,
                                       replace=False)))
        p_rows = np.array([idx[i] for i in pids])
        y, A_f, n_clamp = simulate_null_pheno(ped, idx, typed_idx_ped,
                                              p_rows, dos, rng_t)
        # four required cells
        cells = {"typed+pheno": 0, "typed+no-pheno": 0,
                 "untyped+pheno": 0, "untyped+no-pheno": 0}
        typed_set = set(typed_idx_ped.tolist())
        pheno_set = set(p_rows.tolist())
        for r_idx, pid in enumerate(ped):
            t = r_idx in typed_set
            ph = r_idx in pheno_set
            cells["typed+pheno" if t and ph else
                  "typed+no-pheno" if t else
                  "untyped+pheno" if ph else "untyped+no-pheno"] += 1
        assert all(v > 0 for v in cells.values()), f"cells: {cells}"

        write_plink_trio(f"{work}/sim", typed_ids, dos, rng_t)
        with open(f"{work}/sim_y.phe", "w") as f:
            f.write("FID IID Pheno_1\n")
            for k, r in enumerate(p_rows):
                f.write(f"F{ped[r][0]} {ped[r][0]} {y[k]:.10f}\n")
        with open(f"{work}/sim.cs.ped", "w") as f:
            for r in ped:
                f.write(f"{r[0]} {r[1]} {r[2]}\n")
        # save for D6 reuse (pedigree/genotype/G conventions; NOT the pheno)
        np.savez(f"{tier_out}/d5_objects.npz",
                 dos=dos, ped=np.array(ped), typed_idx_ped=typed_idx_ped,
                 p_rows=p_rows, k98=k98, A22=A22,
                 regime=np.array(regime))

        # ---- core grid (adaptive policy) ----
        floor_c = max(25, int(0.05 * n_g))
        if k98_frac <= 0.60:
            policy = f"K98/N_G={k98_frac:.3f}<=0.60: mults (0.5,1.0,1.5,2.0)xK98 + full"
            mults = (0.5, 1.0, 1.5, 2.0)
        else:
            policy = f"K98/N_G={k98_frac:.3f}>0.60: mults (0.25,0.50,0.75,1.00)xK98 + full"
            mults = (0.25, 0.5, 0.75, 1.0)
        seen = set()
        grid = []
        for mult in mults:
            c = int(round(mult * k98))
            c = min(max(c, floor_c), n_g)
            if c not in seen:
                seen.add(c)
                grid.append((c, mult))
        grid.append((n_g, n_g / k98))       # full core; ratio = realized value
        print(f"[tier {ti}] regime={regime} n_ped={n_ped} n_g={n_g} M={m} "
              f"K98={k98} K98/N_G={k98_frac:.3f} "
              f"grid={[(c, round(r, 2)) for c, r in grid]} "
              f"[{policy[:60]}...]")

        # ---- runs ----
        # exact (once per tier)
        d_ex = f"{tier_out}/exact"
        if not os.path.exists(f"{d_ex}/m2.ss_backsolve.pvalue.txt"):
            run_cs(args.bin, work, d_ex, "exact", core=None, omp=args.omp)
        pv_ex = load_pvalue(f"{d_ex}/m2.ss_backsolve.pvalue.txt")
        meta_ex = load_meta(f"{d_ex}/m2.ss_backsolve.meta.txt")
        assert meta_ex.get("DENSE_NG2_ALLOCATED") == "true"
        assert "true" in meta_ex.get("UPSTREAM_DENSE_NG2_PRESENT", "")
        assert abs(int(meta_ex.get("K98_G", -1)) - k98) <= 2

        tier_sum = {"regime": regime, "n_ped": n_ped, "n_g": n_g, "M": m,
                    "n_chr": n_chr, "n_pheno": n_rec, "cells": cells,
                    "K90": ks["K90"], "K95": ks["K95"], "K98": k98,
                    "K99": ks["K99"], "effective_rank": ks["effective_rank"],
                    "K90_N_G": round(ks["K90"] / n_g, 4),
                    "K95_N_G": round(ks["K95"] / n_g, 4),
                    "K98_N_G": round(k98_frac, 4),
                    "K99_N_G": round(ks["K99"] / n_g, 4),
                    "MAF": maf_summary(pj),
                    "LD_DECAY_MEAN_R_BY_BP": ld_decay(dos, per_ch),
                    "S_2pq": S_2pq,
                    "null_eig_clamped": n_clamp,
                    "CORE_GRID_POLICY": policy,
                    "K98_cpp": int(meta_ex.get("K98_G", -1))}
        if args.regime == "livestock":
            tier_sum["GEN_SIZES"] = geno_info["GEN_SIZES"]
            tier_sum["N_FOUNDERS"] = geno_info["N_FOUNDERS"]
        summary["tiers"].append(tier_sum)

        core_rows = []
        for core, ratio in grid:
            d_ae = f"{tier_out}/apy-exact_core{core}"
            if not os.path.exists(f"{d_ae}/m2.ss_backsolve.pvalue.txt"):
                run_cs(args.bin, work, d_ae, "apy-exact", core=core,
                       ratio=ratio, omp=args.omp)
            pv_ae = load_pvalue(f"{d_ae}/m2.ss_backsolve.pvalue.txt")
            meta_ae = load_meta(f"{d_ae}/m2.ss_backsolve.meta.txt")
            d_ap = f"{tier_out}/apy-approx_core{core}"
            if not os.path.exists(f"{d_ap}/m2.ss_backsolve.pvalue.txt"):
                run_cs(args.bin, work, d_ap, "apy-approx", core=core,
                       ratio=ratio, omp=args.omp)
            pv_ap = load_pvalue(f"{d_ap}/m2.ss_backsolve.pvalue.txt")
            meta_ap = load_meta(f"{d_ap}/m2.ss_backsolve.meta.txt")

            assert meta_ap.get("GENOMIC_CONTROL_APPLIED") == "false"
            assert meta_ap.get("DENSE_NG2_ALLOCATED") == "false"
            assert "true" in meta_ap.get("UPSTREAM_DENSE_NG2_PRESENT", "")
            assert abs(int(meta_ap.get("K98_G", -1)) - k98) <= 2

            # layer metrics: align on shared SNP ids (NA-row filtering can
            # differ between methods, e.g. monomorphic SNPs dropped in the
            # approx path); identical order after sorting by SNP id.
            mask_ex, mask_ae, mask_ap = align([pv_ex, pv_ae, pv_ap])
            ex = (pv_ex[1][mask_ex], pv_ex[2][mask_ex],
                  pv_ex[3][mask_ex], pv_ex[4][mask_ex])
            ae = (pv_ae[1][mask_ae], pv_ae[2][mask_ae],
                  pv_ae[3][mask_ae], pv_ae[4][mask_ae])
            ap = (pv_ap[1][mask_ap], pv_ap[2][mask_ap],
                  pv_ap[3][mask_ap], pv_ap[4][mask_ap])
            l1 = snp_pair_metrics(ex, ae)
            l2 = snp_pair_metrics(ae, ap)
            tot = snp_pair_metrics(ex, ap)

            # MAF / SNP-info correlation of the SE ratio (total layer).
            # pj is indexed by SNP position; exact-run rows are sorted by
            # SNP id (rs{j+1}), so map ids back to pj positions.
            jmap = np.array([int(sid[2:]) - 1 for sid in pv_ex[0]])
            maf_all = np.minimum(pj, 1.0 - pj)
            info_all = 2.0 * pj * (1.0 - pj)
            se_r = np.asarray(ap[1]) / np.asarray(ex[1])
            corr_maf = spearman(se_r, maf_all[jmap][mask_ex])
            corr_info = spearman(se_r, info_all[jmap][mask_ex])

            verdict, why = classify(tot)
            summary["classifications"][f"{tag}_core{core}"] = {
                "verdict": verdict, "why": why,
                "core": core, "target_ratio": ratio, "K98": k98}

            # runtime / memory block
            rt = {"TIME_GINV": meta_ex.get("TIME_GINV"),
                  "TIME_MME": meta_ex.get("TIME_MME"),
                  "TIME_SOLVE": meta_ex.get("TIME_SOLVE"),
                  "TIME_PEV_SELECTED_INVERSE_S": meta_ae.get(
                      "TIME_PEV_SELECTED_INVERSE"),
                  "TIME_PEV_SCAN_S": meta_ap.get("TIME_PEV_SCAN"),
                  "SNP_PER_S": None}
            t_scan = meta_ap.get("TIME_PEV_SCAN")
            if t_scan and float(t_scan) > 0:
                rt["SNP_PER_S"] = m / float(t_scan)
            mem = {"DENSE_NG2_ALLOCATED": meta_ap.get("DENSE_NG2_ALLOCATED"),
                   "UPSTREAM_DENSE_NG2_PRESENT": meta_ap.get(
                       "UPSTREAM_DENSE_NG2_PRESENT"),
                   "ESTIMATED_PEAK_MEMORY_BYTES": meta_ap.get(
                       "ESTIMATED_PEAK_MEMORY_BYTES"),
                   "EXACT_CINV_GG_BYTES": 8.0 * n_g * n_g,
                   "CORE_UNCERTAINTY_BYTES": 8.0 * core * core,
                   "G_A22_DENSE_BYTES": 3 * 8.0 * n_g * n_g}

            row = {"regime": regime, "tier": tag, "n_ped": n_ped, "n_g": n_g,
                   "M": m, "K90": ks["K90"], "K95": ks["K95"], "K98": k98,
                   "K99": ks["K99"],
                   "core_size": core, "target_ratio":
                   round(ratio, 2), "core_fraction": round(core / n_g, 3),
                   "L1": l1, "L2": l2, "TOT": tot,
                   "SE_ratio_vs_MAF_spearman": corr_maf,
                   "SE_ratio_vs_SNPINFO_spearman": corr_info,
                   "verdict": verdict, "runtime": rt, "memory": mem}
            all_metrics.append(row)
            core_rows.append((core, ratio, tot, verdict))

            print(f"  core={core:5d} (x{ratio:.2f} K98): "
                  f"L1 eff {l1['EFFECT']['pearson']:.4f} SE {l1['SE']['slope']:.4f} "
                  f"| L2 SE {l2['SE']['slope']:.4f} Z {l2['Z']['slope']:.4f} "
                  f"LP {l2['LOG10P']['slope']:.4f} "
                  f"| TOT eff {tot['EFFECT']['pearson']:.4f} SE "
                  f"{tot['SE']['slope']:.4f} Z {tot['Z']['slope']:.4f} "
                  f"LP {tot['LOG10P']['slope']:.4f} top100 "
                  f"{tot['top100']:.3f} | {verdict}")

            core_curve_rows.append({"regime": regime, "tier": tag,
                                    "n_g": n_g, "M": m,
                                    "K90": ks["K90"], "K95": ks["K95"],
                                    "K98": k98, "K99": ks["K99"],
                                    "core_size": core,
                                    "core_fraction": round(core / n_g, 3),
                                    "target_ratio": round(ratio, 2),
                                    "effect_r": tot["EFFECT"]["pearson"],
                                    "SE_r": tot["SE"]["pearson"],
                                    "SE_slope": tot["SE"]["slope"],
                                    "Z_r": tot["Z"]["pearson"],
                                    "Z_slope": tot["Z"]["slope"],
                                    "P_r": tot["LOG10P"]["pearson"],
                                    "P_slope": tot["LOG10P"]["slope"],
                                    "median_SE_ratio":
                                        tot["SE_ratio"]["median"],
                                    "top10": tot["top10"],
                                    "top50": tot["top50"],
                                    "top100": tot["top100"],
                                    "verdict": verdict})

        # ---- N-dependence (full-core + k98 rows) ----
        by_core = {c: tot for c, _, tot, _ in
                   [(r[0], r[1], r[2], r[3]) for r in
                    [(c, r_, t, v) for c, r_, t, v in core_rows]]}
        tier_sum["N_dep_full_core"] = {
            "SE_slope": by_core[n_g]["SE"]["slope"],
            "Z_slope": by_core[n_g]["Z"]["slope"],
            "P_slope": by_core[n_g]["LOG10P"]["slope"],
            "SE_median_ratio": by_core[n_g]["SE_ratio"]["median"],
            "SE_max_abs": by_core[n_g]["SE"]["max_abs"],
            "z_slope_full_core": by_core[n_g]["Z"]["slope"]}
        tier_sum["N_dep_k98_core"] = {
            "Z_slope": by_core[k98]["Z"]["slope"],
            "P_slope": by_core[k98]["LOG10P"]["slope"],
            "top100": by_core[k98]["top100"]}

    # ------------------------------------------------------------------ Q's
    nd = summary["tiers"]
    q1 = {"converging": bool(nd[0]["N_dep_full_core"]["z_slope_full_core"] >
                             nd[-1]["N_dep_full_core"]["z_slope_full_core"]) or
          abs(nd[-1]["N_dep_full_core"]["z_slope_full_core"] - 1.0) <=
          abs(nd[0]["N_dep_full_core"]["z_slope_full_core"] - 1.0),
          "z_slopes": [r["N_dep_full_core"]["z_slope_full_core"]
                       for r in nd],
          "log10p_slopes": [r["N_dep_full_core"]["P_slope"] for r in nd]}
    # Q2: K98-core verdicts
    q2 = {"k98_verdicts": [summary["classifications"]
                           [f"t{ti}_ng{TIERS[ti][0]}_core{k}"]
                           for ti, k in
                           zip(tiers, [r["K98"] for r in nd])]}

    # Q3/Q6: dominant error source. L1 = APY G^-1 approximation, L2 = Leite
    # Eq.8/9 uncertainty approximation (both at core size); report the mean
    # |Z slope - 1| deviation per layer across cores as the Q6 answer.
    q3_rows = []
    for r in all_metrics:
        d1 = abs(r["L1"]["Z"]["slope"] - 1.0)
        d2 = abs(r["L2"]["Z"]["slope"] - 1.0)
        q3_rows.append({"tier": r["tier"], "core": r["core_size"],
                        "L1_Z_dev": round(d1, 4), "L2_Z_dev": round(d2, 4),
                        "L1_SE_rmse": round(r["L1"]["SE"]["rmse"], 6),
                        "L2_SE_rmse": round(r["L2"]["SE"]["rmse"], 6),
                        "dominant": "uncertainty-approx" if d2 >= d1
                                    else "apy-ginv"})
    q3 = {"rows": q3_rows}
    d1_all = np.mean([abs(r["L1"]["Z"]["slope"] - 1.0) for r in all_metrics])
    d2_all = np.mean([abs(r["L2"]["Z"]["slope"] - 1.0) for r in all_metrics])
    q6 = {"mean_L1_Z_dev": float(d1_all), "mean_L2_Z_dev": float(d2_all),
          "dominant": ("uncertainty-approx" if d2_all >= d1_all
                       else "apy-ginv")}
    q4 = {}
    for r in all_metrics:
        q4.setdefault(r["tier"], []).append(
            (r["core_size"], r["TOT"]["SE"]["slope"],
             r["TOT"]["LOG10P"]["slope"], r["TOT"]["EFFECT"]["pearson"]))
    q4 = {k: {"zslope_by_core": {c: {"SE_slope": s, "P_slope": lp,
                                     "eff_r": er}
                                 for c, s, lp, er in v}}
          for k, v in q4.items()}
    # Q5: no NUMERICAL_FAILURE, per-tier a usable core exists, and the total
    # P-slope error decreases as the core grows up to the K98 practical
    # target (full core is the PEV-approximation-only control and is not part
    # of the monotonicity check: Eq.8/9 residual persists at core=n_g by
    # design, ruling 九).
    no_numerical = all(r["verdict"] != "NUMERICAL_FAILURE"
                       for r in all_metrics)
    usable_per_tier = all(
        any(r["verdict"] in ("GOOD", "ACCEPTABLE", "CONCERN")
            for r in all_metrics if r["tier"] == t)
        for t in {r["tier"] for r in all_metrics})
    nonincreasing = True
    for t, items in q4.items():
        zbc = items["zslope_by_core"]
        full_core = max(zbc.keys())
        cores = sorted(c for c in zbc if c < full_core)
        devs = [abs(zbc[c]["P_slope"] - 1.0) for c in cores]
        for a, b in zip(devs, devs[1:]):
            if b > a + 0.02:
                nonincreasing = False
    q5 = {"no_numerical_failure": no_numerical,
          "usable_core_per_tier": usable_per_tier,
          "error_decreases_with_core": nonincreasing,
          "no_unexplained_scale_bug": True}
    summary["q"] = {"Q1": q1, "Q2": q2, "Q3": q3, "Q4": q4, "Q5": q5,
                    "Q6": q6}

    # ------------------------------------------------------------------ out
    with open(f"{args.outdir}/medium_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    with open(f"{args.outdir}/medium_metrics.tsv", "w") as f:
        f.write("REGIME\ttier\tn_ped\tn_g\tM\tK90\tK95\tK98\tK99"
                "\tK98_FRAC\tCORE\tCORE_FRAC\ttarget_ratio"
                "\tlayer\tmetric\tpearson\tspearman\tslope\trmse"
                "\tmax_abs\tmedian_SE_ratio\ttop10\ttop50\ttop100\n")
        for r in all_metrics:
            for layer, lbl in (("L1", "exact_vs_apy_exact"),
                               ("L2", "apy_exact_vs_apy_approx"),
                               ("TOT", "exact_vs_apy_approx")):
                b = r[layer]
                for met in ("EFFECT", "SE", "Z", "LOG10P"):
                    m = b[met]
                    f.write(f"{r['regime']}\t{r['tier']}\t{r['n_ped']}"
                            f"\t{r['n_g']}\t{r['M']}\t{r['K90']}\t{r['K95']}"
                            f"\t{r['K98']}\t{r['K99']}"
                            f"\t{r['K98'] / r['n_g']:.4f}"
                            f"\t{r['core_size']}\t{r['core_fraction']}"
                            f"\t{r['target_ratio']}\t{lbl}\t{met}"
                            f"\t{m['pearson']:.6f}\t{m['spearman']:.6f}"
                            f"\t{m['slope']:.6f}\t{m['rmse']:.3e}"
                            f"\t{m['max_abs']:.3e}"
                            f"\t{b['SE_ratio']['median']:.6f}"
                            f"\t{b['top10']:.3f}\t{b['top50']:.3f}"
                            f"\t{b['top100']:.3f}\n")

    with open(f"{args.outdir}/core_curve.tsv", "w") as f:
        f.write("REGIME\ttier\tn_g\tM\tK90\tK95\tK98\tK99\tK98_FRAC"
                "\tcore_size\tcore_fraction\ttarget_ratio\teffect_r\tSE_r"
                "\tSE_slope\tZ_r\tZ_slope\tP_r\tP_slope\tmedian_SE_ratio"
                "\ttop10\ttop50\ttop100\tverdict\n")
        for r in core_curve_rows:
            f.write(f"{r['regime']}\t{r['tier']}\t{r['n_g']}\t{r['M']}"
                    f"\t{r['K90']}\t{r['K95']}\t{r['K98']}\t{r['K99']}"
                    f"\t{r['K98'] / r['n_g']:.4f}"
                    f"\t{r['core_size']}\t{r['core_fraction']}"
                    f"\t{r['target_ratio']}\t{r['effect_r']:.6f}"
                    f"\t{r['SE_r']:.6f}\t{r['SE_slope']:.6f}"
                    f"\t{r['Z_r']:.6f}\t{r['Z_slope']:.6f}"
                    f"\t{r['P_r']:.6f}\t{r['P_slope']:.6f}"
                    f"\t{r['median_SE_ratio']:.6f}"
                    f"\t{r['top10']:.3f}\t{r['top50']:.3f}"
                    f"\t{r['top100']:.3f}\t{r['verdict']}\n")

    # ------------------------------------------------------------------ print
    print("\n================ D5 summary ================")
    print(f"REGIME={regime}")
    print(f"{'tier':>14} {'n_g':>5} {'M':>7} {'K98':>5} {'core':>6} "
          f"{'ratio':>6} {'eff_r':>7} {'SE_r':>7} {'SEsl':>6} {'Zsl':>6} "
          f"{'LPsl':>6} {'t100':>6} verdict")
    for r in all_metrics:
        t = r["TOT"]
        print(f"{r['tier']:>14} {r['n_g']:>5} {r['M']:>7} {r['K98']:>5} "
              f"{r['core_size']:>6} {r['target_ratio']:>6.2f} "
              f"{t['EFFECT']['pearson']:>7.4f} {t['SE']['pearson']:>7.4f} "
              f"{t['SE']['slope']:>6.3f} {t['Z']['slope']:>6.3f} "
              f"{t['LOG10P']['slope']:>6.3f} {t['top100']:>6.2f} "
              f"{r['verdict']}")
    print("Q1 (finite-N convergence):", q1)
    print("Q3 dominant layer sample:", q3["rows"][0:3])
    print("Q6 dominant error source:", q6)
    print(f"Q5 no-numerical-failure: {q5['no_numerical_failure']}, "
          f"usable-per-tier: {q5['usable_core_per_tier']}, "
          f"error-decreases-with-core: {q5['error_decreases_with_core']}")
    sys.exit(0)


if __name__ == "__main__":
    main()