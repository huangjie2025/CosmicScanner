#!/usr/bin/env python3
"""Production QC edge cases, numerical KD oracle, and sample-subset identity."""

import argparse
import csv
import os
import struct
import subprocess
import zlib
from pathlib import Path

import numpy as np

from run_bgen_dominance_production_gate import read_meta, require, write_bgen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2707)
    n, m = 200, 80
    samples = [f"Q{i:04d}" for i in range(n)]
    calls = rng.binomial(2, 0.30, (m, n))
    prior = np.array([0.49, 0.42, 0.09])
    probs = 0.80 * np.eye(3)[calls] + 0.20 * prior
    missing = np.zeros((m, n), dtype=bool)
    missing[0, :11] = True       # > 5%: filter
    missing[1, :10] = True       # exactly 5%: retain
    probs[2] = prior             # uninformative: filter by INFO_D
    probs[3] = 0.05 * np.eye(3)[calls[3]] + 0.95 * prior  # INFO_D < .01: filter
    probs[4] = 0.25 * np.eye(3)[calls[4]] + 0.75 * prior  # warning-only retained
    fixture = work / "qc.bgen"
    decoded = write_bgen(fixture, samples, probs, missing).astype(np.float32)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def run(command):
        result = subprocess.run([str(v) for v in [args.scanner] + command],
                                env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        require(result.returncode == 0, result.stdout[-4000:])
        return result

    grm = work / "oracle"
    run(["--bgen", fixture, "--make-grm", "--make-grm-alg", "4", "--out", grm])
    z_columns, denominator, accepted = [], 0.0, set()
    for j in range(m):
        obs = ~missing[j]
        p0, p1, p2 = decoded[j, obs].astype(float).T
        freq = (p1 + 2 * p2).mean() / 2
        q = 1 - freq
        denom = 4 * freq ** 2 * q ** 2
        mu = -2 * freq ** 2 * p0 + 2 * freq * q * p1 - 2 * q ** 2 * p2
        info = np.var(mu) / denom
        if missing[j].sum() / n > .05 or info < .01 or min(freq, q) < .01:
            continue
        accepted.add(f"rs{j + 1:06d}")
        z = np.zeros(n)
        z[obs] = mu.astype(np.float32)
        z_columns.append(z)
        denominator += denom
    z = np.asarray(z_columns).T
    expected = z @ z.T / denominator
    actual = np.fromfile(str(grm) + ".grm.bin", dtype="<f4")
    error = float(np.max(np.abs(actual - expected[np.tril_indices(n)])))
    require(error < 2e-6, f"production KD differs from posterior probability oracle: {error}")

    phenotype = work / "qc.pheno"
    values = rng.normal(size=n)
    phenotype.write_text("".join(f"{s} {s} {v:.12g}\n" for s, v in zip(samples, values)))
    output = work / "scan"
    run(["--bgen", fixture, "--pheno", phenotype, "--gwas-model", "dominance",
         "--no-loco", "--out", output])
    with Path(str(output) + ".assoc.dom.txt").open() as handle:
        scanned = {row["SNP"] for row in csv.DictReader(handle, delimiter="\t")}
    require(scanned == accepted, "QC marker sets differ between scan and KD")
    require("rs000002" in scanned, "exact 5% missingness boundary incorrectly filtered")
    require("rs000005" in scanned, "entropy/INFO warning incorrectly used as hard filter")
    meta = read_meta(str(output) + ".scan.meta.txt")
    require(int(meta["BGEN_DOM_QC_REJECTED_MISSING"]) == 1, "missingness rejection count")
    require(int(meta["BGEN_DOM_QC_REJECTED_INFO_D"]) == 2, "INFO_D rejection count")
    require(int(meta["BGEN_DOM_QC_INFO_D_WARNINGS"]) > 0, "INFO_D warning not reported")
    require(int(meta["BGEN_DOM_QC_ENTROPY_WARNINGS"]) > 0, "entropy warning not reported")

    # The subset scan must equal an independently materialized subset BGEN,
    # including KD normalization, global/LOCO variance components, and Wald tests.
    keep = np.flatnonzero(np.arange(n) % 4 != 0)
    sub_samples = [samples[i] for i in keep]
    sub_fixture = work / "subset.bgen"
    write_bgen(sub_fixture, sub_samples, probs[:, keep], missing[:, keep])
    sub_pheno = work / "subset.pheno"
    sub_pheno.write_text("".join(f"{samples[i]} {samples[i]} {values[i]:.12g}\n" for i in keep))
    for loco in (False, True):
        outputs = []
        for label, source in (("indexed", fixture), ("materialized", sub_fixture)):
            target = work / f"{label}_{loco}"
            run(["--bgen", source, "--pheno", sub_pheno, "--gwas-model", "dominance",
                 "--out", target] + ([] if loco else ["--no-loco"]))
            with Path(str(target) + ".assoc.dom.txt").open() as handle:
                outputs.append(list(csv.DictReader(handle, delimiter="\t")))
            sub_meta = read_meta(str(target) + ".scan.meta.txt")
            require(sub_meta["FREQUENCY_SAMPLE_SET"] == "ANALYSIS_SAMPLES", "frequency metadata")
        require(outputs[0] == outputs[1], f"sample subset scan differs from materialized BGEN, LOCO={loco}")

    # A one-byte decompressed body must be rejected before header dereferences.
    raw = fixture.read_bytes()
    offset = 4 + struct.unpack_from("<I", raw, 0)[0]
    for _ in range(3):
        length = struct.unpack_from("<H", raw, offset)[0]
        offset += 2 + length
    offset += 6
    for _ in range(2):
        length = struct.unpack_from("<I", raw, offset)[0]
        offset += 4 + length
    compressed = zlib.compress(b"\x00")
    malformed = work / "truncated_header.bgen"
    malformed.write_bytes(raw[:offset] + struct.pack("<II", len(compressed) + 4, 1) + compressed)
    bad = subprocess.run([args.scanner, "--bgen", str(malformed), "--make-grm",
                          "--make-grm-alg", "4", "--out", str(work / "bad")],
                         env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    require(bad.returncode != 0 and "Truncated BGEN genotype header" in bad.stdout,
            "truncated probability header not rejected safely")
    print(f"BGEN DOMINANCE QC PASSED: KD max error={error:.3e}; {len(accepted)} accepted markers; "
          "filter/warning/boundary decisions and global/LOCO sample-subset identity verified")


if __name__ == "__main__":
    main()
