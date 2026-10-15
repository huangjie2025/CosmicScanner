#!/usr/bin/env python3
"""Targeted CLI checks for tail serialization, failed writers and rare boundaries."""
import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

from run_association_r2 import fixture, require


def rows(path):
    with path.open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    binary = args.bin.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    report = {}
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE", MKL_DYNAMIC="FALSE")
    with tempfile.TemporaryDirectory(prefix="r27_", dir=args.out) as directory:
        root = Path(directory)
        data = fixture(root, n=240, m=200)

        def run(name, flags, failure=None, obstruction=None):
            prefix = root / name
            if obstruction:
                suffix, target = obstruction
                path = Path(str(prefix) + suffix)
                if target is None:
                    path.mkdir()
                else:
                    path.symlink_to(target)
            command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(data["pheno"]),
                       "--maf", "0", "--hwe", "0", "--threads", "1", "--out", str(prefix),
                       *map(str, flags)]
            completed = subprocess.run(command, env=env, capture_output=True, text=True)
            log = completed.stdout + completed.stderr
            if failure:
                require(completed.returncode > 0 and failure.lower() in log.lower(),
                        f"{name}: expected clean failure {failure!r}, rc={completed.returncode}\n{log}")
            else:
                require(completed.returncode == 0, f"{name}: {log}")
            return prefix

        phenotype = root / "strong_binary.pheno"
        values = [int(value is not None and value >= 1) for value in data["calls"][0]]
        for i in range(0, data["n"], 17):
            values[i] = 1 - values[i]
        phenotype.write_text("".join(f"{fid} {iid} {values[i]}\n"
                                     for i, (fid, iid) in enumerate(data["samples"])))
        saige_flags = ["--saige-style", "--biobank-mode", "--no-loco", "--sparse-cutoff", "1",
                       "--pheno", phenotype]
        prefix = run("saige_tiny", saige_flags)
        main_rows = rows(Path(str(prefix) + ".cosmic.assoc"))
        metadata = {r["SNP"]: r for r in rows(Path(str(prefix) + ".saige_style.variants.tsv"))}
        require(main_rows, "SAIGE precision fixture produced no rows")
        pvalues = []
        for row in main_rows:
            p = float(row["P"])
            require(math.isfinite(p) and 0 < p <= 1 and float(row["SE"]) > 0,
                    f"Lost positive P/SE: {row}")
            expected_logp = float(metadata[row["SNP"]]["LOG10P"])
            require(math.isclose(-math.log10(p), expected_logp, rel_tol=2e-11, abs_tol=2e-11),
                    f"Main/metadata tail mismatch: {row}")
            pvalues.append(p)
        require(min(pvalues) < 1e-8, "Fixture did not exercise six-decimal tail loss")
        report["saige_min_positive_p"] = min(pvalues)

        # Writers must fail on open and on late flush/close errors. /dev/full
        # exercises the latter even for short gzip output buffered until close.
        run("saige_open_failure", saige_flags, "Cannot write SAIGE-style", (".cosmic.assoc", None))
        lmm_flags = ["--gwas", "--no-loco"]
        run("lmm_open_failure", lmm_flags, "Cannot open association output", (".assoc.txt", None))
        for compressed in (False, True):
            flags = lmm_flags + (["--out-gz"] if compressed else [])
            if Path("/dev/full").exists():
                run(f"lmm_full_{compressed}", flags, "Failed", (".assoc.txt", "/dev/full"))
        if Path("/dev/full").exists():
            run("saige_full", saige_flags, "Failed while writing SAIGE-style", (".cosmic.assoc", "/dev/full"))
        plain = run("lmm_plain", lmm_flags)
        zipped = run("lmm_gzip", lmm_flags + ["--out-gz"])
        require(Path(str(plain) + ".assoc.txt").read_bytes() ==
                gzip.decompress(Path(str(zipped) + ".assoc.txt").read_bytes()),
                "gzip hardening changed association content")
        report["writer_open_flush_close"] = "PASS" if Path("/dev/full").exists() else "OPEN_ONLY"
        report["gzip_plain_parity"] = "PASS"

        covariates = root / "duplicate.qcovar"
        covariates.write_text("".join(f"{fid} {iid} {i % 7} {i % 7}\n"
                                      for i, (fid, iid) in enumerate(data["samples"])))
        # Zero usable SNPs includes absent, monomorphic and entirely missing variants.
        bed_path = Path(str(data["prefix"]) + ".bed")
        bed = bytearray(bed_path.read_bytes())
        stride = (data["n"] + 3) // 4
        bed[3 + 4*stride:3 + 5*stride] = bytes([0xff]) * stride
        bed[3 + 5*stride:3 + 6*stride] = bytes([0x55]) * stride
        bed_path.write_bytes(bed)
        empty_sets = root / "empty.sets"
        empty_sets.write_text("absent unknown\nmono rs0005\nmissing rs0006\n")
        single_set = root / "single.set"
        single_set.write_text("single rs0001\n")
        for method in ("burden", "skat", "skat-o", "acat-v"):
            common = ["--" + method, "--gwas", "--rare-max-maf", ".5", "--rare-min-mac", "1"]
            run(method + "_rank", [*common, "--set-file", single_set, "--qcovar", covariates],
                "rank deficient")
            empty = run(method + "_empty", [*common, "--set-file", empty_sets])
            require(rows(Path(str(empty) + ".rare.txt")) == [], f"{method}: fabricated empty result")
            one = run(method + "_single", [*common, "--set-file", single_set])
            result = rows(Path(str(one) + ".rare.txt"))
            require(len(result) == 1 and result[0]["VARIANT_SET_ID"] == "single" and
                    math.isfinite(float(result[0]["P_VALUE"])) and 0 <= float(result[0]["P_VALUE"]) <= 1,
                    f"{method}: invalid single-variant boundary")
            run(method + "_output_failure", [*common, "--set-file", single_set],
                "Cannot write rare-variant results", (".rare.txt", None))
        report["rare_rank_empty_single_boundaries"] = "PASS_ALL_FOUR_METHODS"
    (args.out / "r27_closure_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
