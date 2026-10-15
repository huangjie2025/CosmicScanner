#!/usr/bin/env python3
"""Phase E2 downstream regression: dense vs block-streamed GRM mode on the
frozen benchmark dataset.

Runs the same cosmicscanner command twice (dense default vs --ssgwas-grm-mode
streamed) for identity-blend and a22-blend modes, then compares:
  - EBV (max |diff|)
  - pvalue (EA, SE, Z, LOG10P: max |diff|)
  - snp (max |diff|)
  - meta fields: FULL_G_MATERIALIZED, FULL_GNN_MATERIALIZED, G_MATERIALIZATION

Tolerance 1e-9. Bit-exact expected for the 44-animal frozen dataset because
the block-streamed builder reproduces the dense Gcc/Gcn/diag blocks exactly.
"""
import argparse, json, os, subprocess, sys
import numpy as np

BIN = None
DATASET = "benchmarks/ssgwas_blupf90_parity/dataset"
CORE_DIR = "benchmarks/ssgwas_blupf90_parity/core_ids"
TOL = 1e-9


def run_cs(out, mode, core_file=None, grm_mode=None):
    cmd = [
        BIN,
        "--bfile", f"{DATASET}/sim",
        "--pheno", f"{DATASET}/sim.pheno",
        "--gwas", "--single-step", "--single-step-method", "backsolve",
        "--ped", f"{DATASET}/sim.ped",
        "--reml-vg", "1.0", "--reml-ve", "2.0", "--no-loco",
        "--ssgwas-window-snps", "10",
        "--ssgwas-solver", "ldlt",
        "--ssgwas-pvalue-method", "apy-exact",
        "--ssgwas-apy-core-file", core_file,
        "--ssgwas-apy-seed", "123",
        "--out", f"{out}/m2",
    ]
    if mode == "a22":
        cmd += ["--ssgwas-tune-g", "--ssgwas-blend-a22"]
    if grm_mode is not None:
        cmd += ["--ssgwas-grm-mode", grm_mode]
    os.makedirs(out, exist_ok=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"CMD: {' '.join(cmd)}", file=sys.stderr)
        raise RuntimeError(f"run failed:\nSTDERR:\n{r.stderr[-1200:]}\nSTDOUT:\n{r.stdout[-600:]}")
    return r


def load_vec(path, col=1):
    return np.array([float(l.split("\t")[col])
                     for l in open(path).read().splitlines()[1:]])


def load_meta(path):
    out = {}
    for l in open(path).read().splitlines():
        if "\t" in l:
            k, v = l.split("\t", 1)
            out[k] = v
    return out


def maxdiff(a, b):
    return float(np.max(np.abs(a - b))) if len(a) == len(b) else float("inf")


def main():
    global BIN
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", required=True)
    ap.add_argument("--outdir", default="tests/results/ssgwas_e2_stream_dense_parity")
    ap.add_argument("--core", default="core44", help="core28 or core44")
    args = ap.parse_args()
    BIN = args.bin
    os.makedirs(args.outdir, exist_ok=True)

    core_file = f"{CORE_DIR}/{args.core}.txt"
    fails = []

    for blend in ["identity", "a22"]:
        out_d = f"{args.outdir}/{blend}_dense"
        out_s = f"{args.outdir}/{blend}_streamed"
        print(f"--- {blend} blend ---")
        print(f"  running dense ...")
        r = run_cs(out_d, blend, core_file=core_file, grm_mode=None)
        print(f"  dense done, out_d={out_d}, files={os.listdir(out_d)}")
        print(f"  running streamed ...")
        r = run_cs(out_s, blend, core_file=core_file, grm_mode="streamed")
        print(f"  streamed done, out_s={out_s}, files={os.listdir(out_s)}")

        # EBV
        ebv_d = load_vec(f"{out_d}/m2.ss_backsolve.ebv.txt")
        ebv_s = load_vec(f"{out_s}/m2.ss_backsolve.ebv.txt")
        d_ebv = maxdiff(ebv_d, ebv_s)
        print(f"  EBV max|diff| = {d_ebv:.3e}", "PASS" if d_ebv < TOL else "FAIL")
        if d_ebv >= TOL:
            fails.append(f"{blend} EBV diff={d_ebv}")

        # pvalue
        for col, name in [(6, "EA"), (7, "SE"), (9, "Z"), (11, "LOG10P")]:
            p_d = load_vec(f"{out_d}/m2.ss_backsolve.pvalue.txt", col)
            p_s = load_vec(f"{out_s}/m2.ss_backsolve.pvalue.txt", col)
            dp = maxdiff(p_d, p_s)
            tag = "PASS" if dp < TOL else "FAIL"
            print(f"  {name} max|diff| = {dp:.3e} {tag}")
            if dp >= TOL:
                fails.append(f"{blend} {name} diff={dp}")

        # snp (compare EFFECT_ALLELE column = col 6)
        snp_d = load_vec(f"{out_d}/m2.ss_backsolve.snp.txt", 6)
        snp_s = load_vec(f"{out_s}/m2.ss_backsolve.snp.txt", 6)
        dsnp = maxdiff(snp_d, snp_s)
        print(f"  SNP effect max|diff| = {dsnp:.3e}", "PASS" if dsnp < TOL else "FAIL")
        if dsnp >= TOL:
            fails.append(f"{blend} SNP diff={dsnp}")

        # meta
        meta_d = load_meta(f"{out_d}/m2.ss_backsolve.meta.txt")
        meta_s = load_meta(f"{out_s}/m2.ss_backsolve.meta.txt")
        for k in ["FULL_G_MATERIALIZED", "FULL_GNN_MATERIALIZED", "G_MATERIALIZATION"]:
            vd = meta_d.get(k, "<missing>")
            vs = meta_s.get(k, "<missing>")
            ok = vd == vs
            print(f"  meta {k}: dense={vd} streamed={vs} {'PASS' if ok else 'MISMATCH'}")
        # Streamed-specific checks
        vs_fullg = meta_s.get("FULL_G_MATERIALIZED", "?")
        vs_gnn = meta_s.get("FULL_GNN_MATERIALIZED", "?")
        vs_mat = meta_s.get("G_MATERIALIZATION", "?")
        if vs_fullg != "false" or vs_gnn != "false" or "BLOCK_STREAMED" not in vs_mat:
            fails.append(f"{blend} meta G_MATERIALIZATION fields wrong: {vs_fullg}/{vs_gnn}/{vs_mat}")

    print()
    if fails:
        print(f"FAILURES ({len(fails)}):")
        for f in fails:
            print(f"  {f}")
        sys.exit(1)
    print("DOWNSTREAM_SSGWAS_PARITY: PASS (dense vs streamed, identity+a22 blend, frozen 44-animal dataset)")


if __name__ == "__main__":
    main()
