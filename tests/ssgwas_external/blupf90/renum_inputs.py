#!/usr/bin/env python3
"""Convert CosmicScanner ssGBLUP inputs (PLINK bed/bim/fam + pedigree) to
BLUPF90 renumbered input + pedigree files (readiness kit; BLUPF90 not present
here, so this is only exercised by hand when the reference is available).

Usage:
  python3 renum_inputs.py --bed PREFIX --ped PED [--pheno PHE] --out DIR
"""
import argparse
import os
import sys
import numpy as np


def read_bed(bed_path, n_ind, n_snp):
    with open(bed_path, "rb") as f:
        magic = f.read(3)
        if magic != b"\x6c\x1b\x01":
            raise SystemExit(f"{bed_path}: not a valid PLINK1 bed (magic {magic.hex()})")
        bytes_per_ind = (n_snp + 3) // 4
        nbytes = n_ind * bytes_per_ind
        data = np.frombuffer(f.read(nbytes), dtype=np.uint8).reshape(n_ind, bytes_per_ind)
    codes = (data[:, :, None] >> np.array([0, 2, 4, 6])) & 3
    dos = 2 - codes.reshape(n_ind, n_snp)
    return np.where(dos > 2, np.nan, dos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed")
    ap.add_argument("--ped")
    ap.add_argument("--pheno")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    fam = np.loadtxt(args.bed + ".fam", dtype=str, usecols=(0, 1))
    bim = np.loadtxt(args.bed + ".bim", dtype=str, usecols=(0, 1, 3, 4, 5))
    n_ind, n_snp = len(fam), len(bim)
    dos = read_bed(args.bed + ".bed", n_ind, n_snp)

    ped = [l.split() for l in open(args.ped) if l.strip()]
    if ped[0][0].lower() == "fid":
        ped = ped[1:]
    id2ped = {int(r[0]): r for r in ped}
    id2row = {i: k for k, i in enumerate(fam[:, 1].astype(int))}

    with open(f"{args.out}/input.dat", "w") as f:
        for i, (fid, iid) in enumerate(fam):
            f.write(f"{i + 1} {i + 1} 1 {iid}\n")
    with open(f"{args.out}/input.ped", "w") as f:
        f.write("id sire dam\n")
        for r in ped:
            f.write(" ".join(r[0:3]) + "\n")
    if args.pheno:
        rows = [l.split() for l in open(args.pheno) if l.strip()]
        if rows[0][0].lower() in ("fid", "id"):
            rows = rows[1:]
        with open(f"{args.out}/input.pheno", "w") as f:
            f.write("id trait\n")
            for r in rows:
                f.write(f"{int(r[0]) + 1} {r[1]}\n")
    print(f"wrote {args.out}/input.dat (n={n_ind}), input.ped (n={len(ped)}), input.pheno (n={n_snp} snp via bim)")


if __name__ == "__main__":
    main()
