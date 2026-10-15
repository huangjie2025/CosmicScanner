#!/usr/bin/env python3
"""Production-path calibration check for probabilistic BGEN dominance."""

import argparse
import csv
import math
import os
import random
import shutil
import struct
import subprocess
import zlib
from pathlib import Path


def pack_string2(value):
    data = value.encode("utf-8")
    return struct.pack("<H", len(data)) + data


def pack_string4(value):
    data = value.encode("utf-8")
    return struct.pack("<I", len(data)) + data


def normal_pdf_delta(delta, tau):
    return math.exp(-0.5*(delta/tau)**2)


def posterior_triplet(g, noise, p, tau):
    evidence = g + tau*noise
    prior = ((1-p)**2, 2*p*(1-p), p*p)
    weights = [prior[k]*normal_pdf_delta(evidence-k, tau) for k in range(3)]
    total = sum(weights)
    return tuple(value/total for value in weights)


def write_bgen(path, sample_ids, freqs, latent, noise, tau, bits=16):
    n = len(sample_ids)
    m = len(freqs)
    sample_block_len = 8 + sum(2+len(s.encode("utf-8")) for s in sample_ids)
    header_len = 20
    flags = 1 | (2 << 2) | (1 << 31)  # zlib, Layout 2, embedded sample IDs
    offset = header_len + sample_block_len
    denominator = (1 << bits)-1

    with Path(path).open("wb") as out:
        out.write(struct.pack("<IIII", offset, header_len, m, n))
        out.write(b"bgen")
        out.write(struct.pack("<I", flags))
        out.write(struct.pack("<II", sample_block_len, n))
        for sid in sample_ids:
            out.write(pack_string2(sid))

        for j, p in enumerate(freqs):
            variant = "u{:04d}".format(j+1)
            chrom = "1" if j < m//2 else "2"
            out.write(pack_string2(variant))
            out.write(pack_string2(variant))
            out.write(pack_string2(chrom))
            out.write(struct.pack("<IH", (j+1)*100, 2))
            out.write(pack_string4("A"))
            out.write(pack_string4("C"))

            body = bytearray(struct.pack("<IHBB", n, 2, 2, 2))
            body.extend([2]*n)
            body.extend((0, bits))  # unphased, bit depth
            for i in range(n):
                p0, p1, p2 = posterior_triplet(latent[j][i], noise[j][i], p, tau)
                # Stored A0A0 is G=2 because the reader's G counts stored A0.
                q2 = int(round(p2*denominator))
                q1 = int(round(p1*denominator))
                # Leave at least one quantization unit for the derived remainder;
                # the strict reader rejects a tiny negative p0 from float addition.
                if q2+q1 >= denominator:
                    if q2 >= q1:
                        q2 -= q2+q1-denominator+1
                    else:
                        q1 -= q2+q1-denominator+1
                body.extend(struct.pack("<HH", q2, q1))
            compressed = zlib.compress(bytes(body), 6)
            out.write(struct.pack("<II", len(compressed)+4, len(body)))
            out.write(compressed)


def write_phenotype(path, sample_ids, values):
    with Path(path).open("w", encoding="utf-8") as out:
        for sid, value in zip(sample_ids, values):
            out.write("{} {} {:.12g}\n".format(sid, sid, value))


def run_scan(scanner, bgen, phenotype, output, env):
    command = [scanner, "--bgen", bgen, "--pheno", phenotype,
               "--gwas", "--gwas-model", "dominance", "--maf", "0.05",
               "--no-loco", "--out", output]
    result = subprocess.run([str(v) for v in command], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env=env)
    if result.returncode != 0:
        raise RuntimeError("scan failed: {}\n{}".format(" ".join(map(str, command)),
                                                        result.stdout[-4000:]))


def read_causal(path, causal_id):
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["SNP"] == causal_id:
                return float(row["BETA_D"]), float(row["SE_D"]), float(row["P_D"])
    raise RuntimeError("causal marker {} missing from {}".format(causal_id, path))


def check(condition, message):
    if not condition:
        raise AssertionError("FAIL: " + message)
    print("  ok: " + message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scanner", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--m", type=int, default=100)
    parser.add_argument("--reps", type=int, default=60)
    args = parser.parse_args()

    workdir = Path(args.workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

    rng = random.Random(2606)
    samples = ["U{:04d}".format(i+1) for i in range(args.n)]
    freqs = [0.10 + 0.80*((j*37) % 97)/96.0 for j in range(args.m)]
    causal = min(12, args.m//4)
    causal_id = "u{:04d}".format(causal+1)
    latent = []
    evidence_noise = []
    for p in freqs:
        column = []
        noise_column = []
        for _ in range(args.n):
            u = rng.random()
            column.append(0 if u < (1-p)**2 else (1 if u < 1-p*p else 2))
            noise_column.append(rng.gauss(0.0, 1.0))
        latent.append(column)
        evidence_noise.append(noise_column)
    causal_h = [1.0 if value == 1 else 0.0 for value in latent[causal]]

    levels = (("low", 0.15), ("medium", 0.50), ("high", 1.20))
    bgen_paths = {}
    for label, tau in levels:
        path = workdir / (label + ".bgen")
        write_bgen(path, samples, freqs, latent, evidence_noise, tau)
        bgen_paths[label] = path

    effect = 1.0
    summaries = {}
    # Pair every phenotype realization across uncertainty levels.
    residuals = [[rng.gauss(0.0, 1.0) for _ in range(args.n)]
                 for _ in range(args.reps)]
    for label, _ in levels:
        null_p = []
        alt_p = []
        alt_beta = []
        alt_se = []
        for rep in range(args.reps):
            for arm in ("null", "alt"):
                values = residuals[rep] if arm == "null" else [
                    residuals[rep][i] + effect*causal_h[i] for i in range(args.n)]
                phenotype = workdir / ("{}_{}_{}.pheno".format(label, arm, rep))
                output = workdir / ("{}_{}_{}".format(label, arm, rep))
                write_phenotype(phenotype, samples, values)
                run_scan(args.scanner, bgen_paths[label], phenotype, output, env)
                beta, se, pvalue = read_causal(str(output)+".assoc.dom.txt", causal_id)
                if arm == "null":
                    null_p.append(pvalue)
                else:
                    alt_p.append(pvalue)
                    alt_beta.append(beta)
                    alt_se.append(se)
        summaries[label] = {
            "type1": sum(value < 0.05 for value in null_p)/args.reps,
            "power": sum(value < 0.05 for value in alt_p)/args.reps,
            "beta": sum(alt_beta)/args.reps,
            "se": sum(alt_se)/args.reps,
        }

    print("quality  type-I   power  mean-beta  mean-SE")
    for label, _ in levels:
        row = summaries[label]
        print("{:<7} {:7.4f} {:7.4f} {:10.4f} {:8.4f}".format(
            label, row["type1"], row["power"], row["beta"], row["se"]))

    for label, _ in levels:
        check(summaries[label]["type1"] <= 0.15,
              "{} uncertainty type-I is bounded ({:.3f})".format(
                  label, summaries[label]["type1"]))
        check(abs(summaries[label]["beta"]-effect) <= 0.45,
              "{} uncertainty beta is calibrated ({:.3f})".format(
                  label, summaries[label]["beta"]))
    check(summaries["low"]["power"] >= 0.80,
          "low-uncertainty power is high ({:.3f})".format(summaries["low"]["power"]))
    check(summaries["low"]["power"] >= summaries["high"]["power"],
          "power does not improve from low to high uncertainty")
    check(summaries["low"]["se"] < summaries["high"]["se"],
          "SE increases from low to high uncertainty")
    print("BGEN DOMINANCE UNCERTAINTY VALIDATION PASSED")


if __name__ == "__main__":
    main()
