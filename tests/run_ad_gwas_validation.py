#!/usr/bin/env python3
"""Phase 2B Additive-Dominance GWAS validation driver.

Generates PLINK BED/bim/fam data with an additive quantitative-trait locus and a
separate dominance-segregating locus, runs CosmicScanner under the default
--gwas-model additive and the new --gwas-model additive_dominance, then checks:

  (a) Null-dominance similarity: when sigma_D ~ 0 (no dominance noise OR a pure
      additive causal locus), -log10(P_A) from the AD scan tracks the standard
      additive P across shared SNPs.
  (b) Dominance QTL detectability: a pure-dominance locus is detected by the AD
      scan (small P_D / P_JOINT) while the additive-alone scan misses it.
  (c) Output format: header, column count, finite numerics, no NaN rows kept.

Usage:
  python3 run_ad_gwas_validation.py --scanner <path> [--workdir <dir>]
"""
import argparse
import math
import os
import shutil
import struct
import subprocess
import sys

WORKDIR_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "ad_gwas_p2b")
M = 400
N = 400
AD_ADD_SNP = 17     # additive causal SNP index (0-based)
AD_DOM_SNP = 41     # dominance causal SNP index (0-based)
SEED = 1234


def rng(seed):
    return seed * 69621 + 7


def next_double(s):
    s = (s * 1103515245 + 12345) & 0x7FFFFFFF
    return s / 0x7FFFFFFF, s


def generate_plink(workdir, seed=SEED):
    os.makedirs(workdir, exist_ok=True)
    s = seed
    n, m = N, M
    prefix = os.path.join(workdir, "sim")
    with open(prefix + ".bim", "w") as f:
        for j in range(m):
            a1 = "G" if s % 2 else "A"
            s = (s * 1103515245 + 12345) & 0x7FFFFFFF
            a2 = "A" if a1 == "G" else "G"
            f.write(f"{j % 3 + 1}\tsnp{j}\t0\t{j * 1000}\t{a1}\t{a2}\n")
    genotypes = [[0] * n for _ in range(m)]
    for j in range(m):
        p = 0.3 + 0.3 * ((j * 37) % 100) / 100.0
        for i in range(n):
            u, s = next_double(s)
            g = 0
            if u < p * p:
                g = 2
            elif u < p * p + 2 * p * (1 - p):
                g = 1
            genotypes[j][i] = g
    # Dominance causal SNP: heterozygotes over-represented -> clear HWD
    for i in range(n):
        genotypes[AD_DOM_SNP][i] = 0 if (i % 4 == 0) else (1 if (i % 4 in (1, 2)) else 2)
    # Additive causal SNP: keep all three genotypic classes so the joint
    # additive+dominance design is full rank (AA/Aa/aa present).
    for i in range(n):
        genotypes[AD_ADD_SNP][i] = (i % 3)
    dosages = [[float(genotypes[j][i]) for i in range(n)] for j in range(m)]

    def centered(v):
        return [v[i] - sum(v) / n for i in range(n)]

    g_add = centered([dosages[AD_ADD_SNP][i] for i in range(n)])
    het = [1.0 if genotypes[AD_DOM_SNP][i] == 1 else 0.0 for i in range(n)]
    g_dom = centered(het)
    var_add = sum(x * x for x in g_add) / (n - 1)
    var_dom = sum(x * x for x in g_dom) / (n - 1)

    # Dataset A: pure additive trait (sigma_D = 0)
    y_av = []
    for i in range(n):
        u, s = next_double(s)
        e = (u - 0.5) * 2.0 * math.sqrt(0.6)
        y_av.append(g_add[i] + e)
    # Dataset B: additive + strong dominance QTL
    y_dv = []
    for i in range(n):
        u, s = next_double(s)
        e = (u - 0.5) * 2.0 * math.sqrt(0.6)
        y_dv.append(g_add[i] + 2.0 * g_dom[i] + e)

    with open(os.path.join(workdir, "add.pheno"), "w") as f:
        for i in range(n):
            f.write(f"F{i}\tI{i}\t{y_av[i]:.6f}\n")
    with open(os.path.join(workdir, "ad.pheno"), "w") as f:
        for i in range(n):
            f.write(f"F{i}\tI{i}\t{y_dv[i]:.6f}\n")
    with open(prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"F{i}\tI{i}\t0\t0\t0\t1\n")
    with open(prefix + ".bed", "wb") as f:
        f.write(b"\x6c\x1b\x01")
        remain = n
        byte_len = (n + 3) // 4
        for j in range(m):
            b = bytearray(byte_len)
            for i in range(n):
                # PLINK BED 2-bit codes: 00=hom ref, 01=missing, 10=het, 11=hom alt
                code = {0: 0, 1: 2, 2: 3}[genotypes[j][i]]
                b[i >> 2] |= code << ((i & 3) * 2)
            f.write(bytes(b))
    return prefix


def read_lines(path):
    with open(path) as f:
        return [l.rstrip("\n") for l in f if l.rstrip("\n")]


def parse_assoc(prefix):
    lines = read_lines(prefix)
    header = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:]]
    return header, rows


def col(header, name):
    return header.index(name)


def run_scanner(scanner, prefix, out_prefix, pheno=None, extra=(), loco=False):
    cmd = [scanner, "--bfile", prefix, "--out", out_prefix, "--gwas",
           "--maf", "0.01"]
    if not loco:
        cmd += ["--no-loco"]
    if pheno:
        cmd += ["--pheno", pheno]
    cmd += list(extra)
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    return r


def check(cond, msg):
    if not cond:
        raise AssertionError("FAIL: " + msg)
    print("  ok: " + msg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scanner", required=True)
    ap.add_argument("--workdir", default=WORKDIR_DEFAULT)
    ap.add_argument("--keep", action="store_true", help="keep generated data")
    args = ap.parse_args()
    if not os.path.exists(args.scanner):
        sys.exit("scanner not found: " + args.scanner)
    if not args.keep and os.path.exists(args.workdir):
        shutil.rmtree(args.workdir)
    prefix = generate_plink(args.workdir, SEED)

    print("=== Run default additive (pure-additive trait) ===")
    add_pheno = os.path.join(args.workdir, "add.pheno")
    ra = run_scanner(args.scanner, prefix, os.path.join(args.workdir, "add"), pheno=add_pheno)
    print(ra.stdout[-2500:])
    if ra.returncode != 0:
        sys.exit("additive run failed:\n" + ra.stderr)
    add_file = os.path.join(args.workdir, "add.assoc.txt")
    header_a, rows_a = parse_assoc(add_file)
    expected = ["CHR", "SNP", "BP", "A1", "A2", "FREQ", "N",
                "BETA_A", "SE_A", "P_A", "LOG10P_A",
                "BETA_D", "SE_D", "P_D", "LOG10P_D",
                "P_JOINT", "LOG10P_JOINT"]
    check(header_a == expected, "additive header matches unified schema")
    check(all(r_[col(header_a, name)] == "NA" for r_ in rows_a
              for name in ("BETA_D", "SE_D", "P_D", "LOG10P_D",
                           "P_JOINT", "LOG10P_JOINT")),
          "additive unavailable fields are NA")

    print("=== Run additive_dominance (pure-additive trait) ===")
    rb = run_scanner(args.scanner, prefix, os.path.join(args.workdir, "ad"),
                     pheno=add_pheno,
                     extra=["--gwas-model", "additive_dominance"])
    print(rb.stdout[-2500:])
    if rb.returncode != 0:
        sys.exit("additive_dominance run failed:\n" + rb.stderr)
    ad_file = os.path.join(args.workdir, "ad.assoc.ad.txt")
    if not os.path.exists(ad_file):
        sys.exit("AD output missing: " + ad_file)
    header_b, rows_b = parse_assoc(ad_file)

    print("=== (c) output format checks ===")
    check(header_b == expected, "AD header matches spec")
    check(len(rows_b) > 0, "AD output non-empty")
    for r_ in rows_b:
        check(len(r_) == len(header_b), "row length == header length")
    for name in ("BETA_A", "SE_A", "P_A", "LOG10P_A", "BETA_D", "SE_D", "P_D", "LOG10P_D",
                 "P_JOINT", "LOG10P_JOINT"):
        vals = [float(r_[col(header_b, name)]) for r_ in rows_b]
        check(all(math.isfinite(v) for v in vals), name + " all finite")
    pj = [float(r_[col(header_b, "P_JOINT")]) for r_ in rows_b]
    check(all(0.0 <= v <= 1.0 for v in pj), "P_JOINT in [0,1]")
    logp = [float(r_[col(header_b, "LOG10P_JOINT")]) for r_ in rows_b]
    check(all(v >= 0.0 for v in logp), "LOG10P_JOINT >= 0")

    print("=== (a) additive null-dominance similarity ===")
    import math as _m
    snp_col_a = col(header_a, "SNP")
    p_col_a = col(header_a, "P_A")
    p_add_a = {}
    for r_ in rows_a:
        p_add_a[r_[snp_col_a]] = float(r_[p_col_a])
    p_add_b = {}
    for r_ in rows_b:
        p_add_b[r_[col(header_b, "SNP")]] = float(r_[col(header_b, "P_A")])
    common = [s for s in p_add_a if s in p_add_b]
    check(len(common) > M * 0.9, f"shared SNPs across run types ({len(common)}/{M})")
    la = [_m.log10(p_add_a[s]) for s in common]
    lb = [_m.log10(p_add_b[s]) for s in common]
    corr = corr_of(la, lb)
    check(corr > 0.95, f"-log10 P_A corr vs standard LMM = {corr:.4f} > 0.95")

    print("=== (b) dominance QTL detectability (LOCO: avoid GRM absorbing own causal signal) ===")
    dom_pheno = os.path.join(args.workdir, "ad.pheno")
    rd = run_scanner(args.scanner, prefix, os.path.join(args.workdir, "add_only_on_dom_trait"),
                     pheno=dom_pheno, loco=True)
    print(rd.stdout[-800:])
    if rd.returncode != 0:
        sys.exit("additive run (dominance trait) failed:\n" + rd.stderr)
    _, rows_a_dom = parse_assoc(os.path.join(args.workdir, "add_only_on_dom_trait.assoc.txt"))
    p_add_a_dom = {r_[col(header_a, "SNP")]: float(r_[col(header_a, "P_A")]) for r_ in rows_a_dom}

    rc = run_scanner(args.scanner, prefix, os.path.join(args.workdir, "ad_dom"),
                     pheno=dom_pheno, loco=True,
                     extra=["--gwas-model", "additive_dominance"])
    print(rc.stdout[-800:])
    if rc.returncode != 0:
        sys.exit("additive_dominance run (dominance trait) failed:\n" + rc.stderr)
    header_c, rows_c = parse_assoc(os.path.join(args.workdir, "ad_dom.assoc.ad.txt"))

    dom_rows = {r_[col(header_c, "SNP")]: r_ for r_ in rows_c}
    ad_dom = dom_rows.get("snp" + str(AD_DOM_SNP))
    check(ad_dom is not None, "dominance SNP row present")
    p_dom_val = float(ad_dom[col(header_c, "P_D")])
    p_joint_val = float(ad_dom[col(header_c, "P_JOINT")])
    check(p_dom_val < 1e-4, f"dominance SNP P_D significant ({p_dom_val:.3e})")
    check(p_joint_val < 1e-4, f"dominance SNP P_JOINT significant ({p_joint_val:.3e})")
    add_only_p_dom = p_add_a_dom.get("snp" + str(AD_DOM_SNP), 1.0)
    # HWD locus: dosage and het-indicator correlate, so additive-only has some signal;
    # the AD dominance-conditional test must be decisively (>= 1e4x) more significant.
    check(p_dom_val < add_only_p_dom / 1e4,
          f"AD dominance test beats additive-only by >=1e4x "
          f"(P_D={p_dom_val:.3e} vs additive-only P={add_only_p_dom:.3e})")

    pj_all = {r_[col(header_c, "SNP")]: float(r_[col(header_c, "P_JOINT")]) for r_ in rows_c}
    top = sorted(pj_all, key=pj_all.get)[:10]
    check("snp" + str(AD_DOM_SNP) in top, "dominance SNP in top-10 joint hits")
    check("snp" + str(AD_ADD_SNP) in top, "additive SNP in top-10 joint hits")

    for prefix, model in ((os.path.join(args.workdir, "add"), "additive"),
                          (os.path.join(args.workdir, "ad"), "additive_dominance")):
        meta = open(prefix + ".scan.meta.txt", encoding="utf-8").read()
        for token in (f"MODEL\t{model}", f"N\t{N}", "N_SEMANTICS\t",
                      "FREQ_SEMANTICS\t", "LOCO\tfalse",
                      "COVARIATE_PROJECTION\ttrue", "VCE_CONVERGED\t"):
            check(token in meta, f"{model} metadata contains {token}")
    ad_meta = open(os.path.join(args.workdir, "ad.scan.meta.txt"), encoding="utf-8").read()
    for token in ("SIGMA_A2\t", "SIGMA_D2\t", "SIGMA_E2\t"):
        check(token in ad_meta, f"AD metadata contains {token}")

    print("= PHASE 2B AD GWAS VALIDATION PASSED =")
    print("workdir:", args.workdir)


def corr_of(a, b):
    n = len(a)
    ma = sum(a) / n
    mb = sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    if den <= 0:
        return 0.0
    return num / den


if __name__ == "__main__":
    main()
