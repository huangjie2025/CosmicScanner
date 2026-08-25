#!/usr/bin/env python3
"""Generate pedigree-based animal-model GWAS simulation data.

Simulates a 3-generation pedigree:
  gen0: 2*N_F founders (untyped, no recorded parents)
  gen1: 2*N_F children from N_F founder couples (untyped)
  gen2: N_TYPE typed analysis individuals, parents drawn from gen1

Additive breeding values are simulated recursively over the pedigree
(founders ~ N(0, sigma2_a); offspring = 0.5*(a_sire+a_dam) + m, with
Mendelian sampling variance adjusted for unknown parents). The causal SNP is
independently drawn (HWE); phenotype y = mu + beta*(g-2p) + a + e.

Outputs (in OUT_DIR):
  <prefix>.bed/.bim/.fam/.pheno   typed individuals (BED written with the
                                  CosmicScanner plink conventions used by
                                  tests/data/gemma_accuracy)
  pedigree.txt                    ID SIRE DAM (0 = unknown)
  truth.txt                       simulation truths for validation

Options control the tested scenarios: sigma2_a = 0 (no polygenic signal),
high heritability, incomplete pedigrees (missing one/both parents), shuffled
pedigree row order, and no-ancestor pedigrees.
"""
import argparse
import random
import struct
import sys
from pathlib import Path


def plink_code(dosage):
    if dosage == 2:
        return 0b00
    if dosage == 1:
        return 0b10
    if dosage == 0:
        return 0b11
    return 0b01


def write_plink(prefix, samples, variants, genotypes, phenotypes):
    n = len(samples)
    with prefix.with_suffix(".fam").open("w", encoding="utf-8") as f:
        for (fid, iid), y in zip(samples, phenotypes):
            f.write(f"{fid}\t{iid}\t0\t0\t1\t{y:.8f}\n")
    with prefix.with_suffix(".pheno").open("w", encoding="utf-8") as f:
        for (fid, iid), y in zip(samples, phenotypes):
            f.write(f"{fid}\t{iid}\t{y:.8f}\n")
    with prefix.with_suffix(".bim").open("w", encoding="utf-8") as f:
        for chrom, snp, bp, a1, a2 in variants:
            f.write(f"{chrom}\t{snp}\t0\t{bp}\t{a1}\t{a2}\n")
    bytes_per_snp = (n + 3) // 4
    bed = bytearray(b"\x6c\x1b\x01")
    for row in genotypes:
        block = bytearray(bytes_per_snp)
        for i, g in enumerate(row):
            block[i // 4] |= plink_code(g) << ((i % 4) * 2)
        bed.extend(block)
    prefix.with_suffix(".bed").write_bytes(bed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="animal_model")
    ap.add_argument("--n-type", type=int, default=150, help="typed analysis individuals (gen2)")
    ap.add_argument("--n-founders", type=int, default=15, help="founder couples (gen0)")
    ap.add_argument("--m", type=int, default=600, help="number of SNPs")
    ap.add_argument("--h2", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=20260815)
    ap.add_argument("--causal-effect", type=float, default=0.6)
    ap.add_argument("--causal-maf", type=float, default=0.3)
    ap.add_argument("--missing-one-parent", type=float, default=0.0,
                    help="prob. a typed individual lacks one recorded parent")
    ap.add_argument("--missing-both-parents", type=float, default=0.0,
                    help="prob. a typed individual lacks both recorded parents")
    ap.add_argument("--shuffle-pedigree", action="store_true",
                    help="write pedigree rows in shuffled order (parents may follow offspring)")
    ap.add_argument("--no-ancestors", action="store_true",
                    help="pedigree contains only the typed individuals")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    prefix = out / args.prefix

    n_type = args.n_type
    n_f = args.n_founders
    n_gen0 = 2 * n_f
    n_gen1 = 2 * n_f
    sigma2_a = args.h2
    sigma2_e = 1.0 - args.h2

    # --- pedigree: (id, sire_id or None, dam_id or None) ---
    pedigree = []  # gen0 founders first
    for i in range(n_gen0):
        pedigree.append((f"F{i:04d}", None, None))
    for i in range(n_f):  # gen1 children of founder couples
        pedigree.append((f"M{i*2:04d}", f"F{i:04d}", f"F{i+n_f:04d}"))
        pedigree.append((f"M{i*2+1:04d}", f"F{i:04d}", f"F{i+n_f:04d}"))
    gen1_ids = [f"M{i:04d}" for i in range(n_gen1)]
    for t in range(n_type):
        sire = dam = None
        if rng.random() >= args.missing_both_parents:
            sire = gen1_ids[rng.randrange(n_gen1)]
            if rng.random() >= args.missing_one_parent:
                dam = gen1_ids[rng.randrange(n_gen1)]
        if args.no_ancestors:
            sire = dam = None
        pedigree.append((f"T{t:04d}", sire, dam))
    n_ped = len(pedigree)

    # --- recursive breeding values (parents listed before offspring by construction) ---
    a = {}
    for pid, sire, dam in pedigree:
        has_s = sire is not None
        has_d = dam is not None
        term = 0.0
        if has_s:
            term += 0.5 * a[sire]
        if has_d:
            term += 0.5 * a[dam]
        m_var = sigma2_a * (1.0 - 0.25 * (1 if has_s else 0) - 0.25 * (1 if has_d else 0))
        a[pid] = term + rng.gauss(0.0, m_var ** 0.5)

    # --- genotypes for typed individuals (HWE), causal SNP first ---
    typed = [pid for pid, _, _ in pedigree if pid.startswith("T")]
    assert len(typed) == n_type
    variants = []
    genotypes = []
    for j in range(args.m):
        chrom = str((j % 4) + 1)
        snp = f"rs{j+1:06d}"
        bp = 1000 + j * 10
        p = args.causal_maf if j == 0 else 0.10 + 0.30 * ((j * 37) % 100) / 99.0
        a1 = "C" if p <= 0.5 else "A"
        a2 = "A" if p <= 0.5 else "C"
        # keep the effect allele as the minor allele (a1) for the causal SNP
        row = []
        for _ in typed:
            g = (1 if rng.random() < p else 0) + (1 if rng.random() < p else 0)
            row.append(g)
        variants.append((chrom, snp, bp, a1, a2))
        genotypes.append(row)

    # --- phenotypes ---
    mu = 0.0
    beta = args.causal_effect
    p_causal = args.causal_maf
    phenotypes = []
    for t, pid in enumerate(typed):
        y = mu + beta * (genotypes[0][t] - 2.0 * p_causal) + a[pid]
        y += rng.gauss(0.0, sigma2_e ** 0.5)
        phenotypes.append(y)

    # --- sample records (fam) and shuffled pedigree file ---
    samples = [("P", pid) for pid in typed]
    write_plink(prefix, samples, variants, genotypes, phenotypes)

    ped_rows = []
    for pid, sire, dam in pedigree:
        ped_rows.append(f"{pid}\t{sire if sire is not None else '0'}\t{dam if dam is not None else '0'}")
    if args.shuffle_pedigree:
        rng.shuffle(ped_rows)
    (out / "pedigree.txt").write_text("\n".join(ped_rows) + "\n", encoding="utf-8")

    with (out / "truth.txt").open("w", encoding="utf-8") as f:
        f.write(f"n_pedigree\t{n_ped}\n")
        f.write(f"n_type\t{n_type}\n")
        f.write(f"sigma2_a\t{sigma2_a}\n")
        f.write(f"sigma2_e\t{sigma2_e}\n")
        f.write(f"h2\t{args.h2}\n")
        f.write(f"causal_snp\trs1\n")
        f.write(f"causal_beta\t{beta}\n")
        f.write(f"causal_maf\t{p_causal}\n")

    print(f"Generated animal-model simulation data in {out}")
    print(f"  pedigree n={n_ped} (typed={n_type}, untyped ancestors={n_ped - n_type})")
    print(f"  SNPs={args.m} causal=rs1 beta={beta} h2={args.h2}")


if __name__ == "__main__":
    main()