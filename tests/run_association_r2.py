#!/usr/bin/env python3
"""Stdlib-only CLI regression harness for association expansion R2.

Example: python3 tests/run_association_r2.py --bin build/cosmicscanner \
    --case genetic_model --out build/test-results/association-r2

Fixtures and command outputs live in a fresh TemporaryDirectory under --out.
Command failures include the log before that directory is cleaned.
The benchmark leaves benchmark.time.tsv; timings are informational, not gates.
SKAT-O checks cover the fixed-grid dependence-adjusted implementation.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import random
import re
import shlex
import subprocess
import sys
import tempfile
import time


CASES = ("pca", "qc_filter", "ld_prune", "ld_clump", "skat_o", "acat_v",
         "genetic_model", "benchmark")
MODELS = ("additive", "dominant", "recessive", "heterozygote")
LD_ARGS = ["--ld-window", "250000", "--ld-step", "1", "--ld-r2", "0.8"]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def table(path, header=None):
    with path.open(encoding="ascii", newline="") as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    require(bool(rows), f"Empty table: {path}")
    names = rows.pop(0)
    require(len(names) == len(set(names)), f"Duplicate columns: {path}")
    if header is not None:
        require(names == header.split(), f"Unexpected header in {path}: {names}")
    require(bool(rows), f"No data rows: {path}")
    require(all(len(row) == len(names) for row in rows), f"Ragged table: {path}")
    return [dict(zip(names, row)) for row in rows]


def probability(value, context):
    value = float(value)
    require(math.isfinite(value) and 0 <= value <= 1, f"Invalid P in {context}: {value}")


def fixture(root, n=40, m=24):
    rng = random.Random(20260908)
    calls = []
    for j in range(m):
        column = [rng.choices((0, 1, 2), (0.35, 0.45, 0.20))[0] for _ in range(n)]
        # Every SNP supports all four models, even after the QC sample subset.
        column[n - 5:n - 2] = [0, 1, 2]
        if j >= 2:
            column[0] = None
            column[1 + (7 * j + 11) % (n - 8)] = None
        calls.append(column)
    calls[0][5] = None
    calls[1] = calls[0][:]
    calls[2][:n // 4] = [None] * (n // 4)
    # Identical across chromosomes must not be pruned/clumped together.
    calls[m // 2] = calls[0][:]
    samples = [(f"F{i:03d}", f"I{i:03d}") for i in range(n)]
    variants = [(str(1 + j // (m // 2)), f"rs{j + 1:04d}",
                 1000 + (j % (m // 2)) * 10000) for j in range(m)]
    fam = [f"{fid} {iid} 0 0 {1 + i % 2} -9" for i, (fid, iid) in enumerate(samples)]
    bim = [f"{chrom} {snp} 0 {bp} A G" for chrom, snp, bp in variants]
    prefix = root / "fixture"
    Path(str(prefix) + ".fam").write_text("\n".join(fam) + "\n", encoding="ascii")
    Path(str(prefix) + ".bim").write_text("\n".join(bim) + "\n", encoding="ascii")
    stride = (n + 3) // 4
    packed = bytearray(b"\x6c\x1b\x01") + bytearray(stride * m)
    # CosmicScanner counts BIM A1: BED 00 -> 2, 10 -> 1, 11 -> 0, 01 -> missing.
    codes = {2: 0, 1: 2, 0: 3, None: 1}
    for j, column in enumerate(calls):
        for i, value in enumerate(column):
            packed[3 + j * stride + i // 4] |= codes[value] << (2 * (i % 4))
    Path(str(prefix) + ".bed").write_bytes(packed)
    y = [round(3 + 0.65 * (calls[0][i] if calls[0][i] is not None else 0.8)
               + 0.35 * (calls[3][i] == 1) + rng.gauss(0, 1), 8) for i in range(n)]
    pheno = root / "phenotype.txt"
    # Reverse order catches accidental positional phenotype alignment.
    pheno.write_text("".join(f"{samples[i][0]} {samples[i][1]} {y[i]:.8f}\n"
                             for i in reversed(range(n))), encoding="ascii")
    pvalues = [0.02 + j / (2 * m) for j in range(m)]
    pvalues[0], pvalues[1], pvalues[m // 2] = 1e-9, 1e-12, 1e-8
    summary = root / "summary.tsv"
    summary.write_text("SNP\tP\tCHR\tBP\n" + "".join(
        f"{snp}\t{pvalues[j]:.17g}\t{chrom}\t{bp}\n"
        for j, (chrom, snp, bp) in reversed(list(enumerate(variants)))), encoding="ascii")
    sets = {"duplicate_set": [0, 1, 3], "missing_set": [2, 4, 5],
            "chr2_set": [m // 2, m // 2 + 1, m // 2 + 2]}
    set_file = root / "sets.txt"
    set_file.write_text("".join(f"{name} {variants[j][1]}\n"
                                for name, indices in sets.items() for j in indices), encoding="ascii")
    return dict(prefix=prefix, calls=calls, samples=samples, variants=variants,
                fam=fam, bim=bim, y=y, pheno=pheno, summary=summary,
                pvalues=pvalues, sets=sets, set_file=set_file, n=n, m=m)


class Harness:
    def __init__(self, binary, root, data):
        self.binary, self.root, self.data = binary, root, data
        self.timings = []

    def run(self, name, args, error=None, occupied=None, pgen=False):
        prefix = self.root / f"{len(self.timings):03d}_{name}"
        if occupied:
            Path(str(prefix) + occupied).write_bytes(b"exclusive-create sentinel\n")
        source = ["--pfile", self.root / "unsupported_pgen"] if pgen else [
            "--bfile", self.data["prefix"]]
        command = [str(x) for x in [self.binary, *source, "--out", prefix,
                                   "--threads", "1", "--maf", "0", "--hwe", "0", *args]]
        log = Path(str(prefix) + ".command.log")
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                   MKL_NUM_THREADS="1", LC_ALL="C")
        start = time.perf_counter()
        # Direct file redirection, never a shell or PIPE; inspect every return code.
        with log.open("wb") as stream:
            result = subprocess.run(command, cwd=self.root, env=env,
                                    stdout=stream, stderr=subprocess.STDOUT, check=False)
        elapsed = time.perf_counter() - start
        self.timings.append((name, self.data["n"], self.data["m"], elapsed, result.returncode))
        text = log.read_text(encoding="utf-8", errors="replace")
        diagnostic = f"{shlex.join(command)}\nreturncode={result.returncode}\n{text}"
        if error is None:
            require(result.returncode == 0, diagnostic)
        else:
            require(result.returncode > 0, f"Expected a clean nonzero error exit:\n{diagnostic}")
            require("unknown option" not in text.lower(), f"Wrong failure (unknown flag):\n{diagnostic}")
            require(re.search(error, text, re.IGNORECASE), f"Wrong error; expected {error!r}:\n{diagnostic}")
        if occupied:
            require(Path(str(prefix) + occupied).read_bytes() == b"exclusive-create sentinel\n",
                    f"Existing output was overwritten: {prefix}{occupied}")
        return prefix


def check_qc(prefix, data, samples, variants):
    qc = str(prefix) + ".qc"
    require(Path(qc + ".fam").read_text().splitlines() == [data["fam"][i] for i in samples],
            "QC FAM subset/order/metadata mismatch")
    require(Path(qc + ".bim").read_text().splitlines() == [data["bim"][j] for j in variants],
            "QC BIM subset/order/metadata mismatch")
    raw = Path(qc + ".bed").read_bytes()
    stride = (len(samples) + 3) // 4
    require(raw[:3] == b"\x6c\x1b\x01" and len(raw) == 3 + stride * len(variants),
            "QC BED magic or exact size mismatch")
    decode = (2, None, 1, 0)
    for k, j in enumerate(variants):
        observed = [decode[(raw[3 + k * stride + i // 4] >> (2 * (i % 4))) & 3]
                    for i in range(len(samples))]
        require(observed == [data["calls"][j][i] for i in samples],
                f"QC altered hardcalls/missingness for variant {j}")


def pca_case(h, failures=True):
    prefix = h.run("pca", ["--pca", "--pca-k", "3"])
    rows = table(Path(str(prefix) + ".eigenvec"), "FID IID PC1 PC2 PC3")
    require([(r["FID"], r["IID"]) for r in rows] == h.data["samples"], "PCA sample IDs/order")
    columns = [[float(r[f"PC{k}"]) for r in rows] for k in range(1, 4)]
    require(all(math.isfinite(x) for col in columns for x in col), "Nonfinite PCA scores")
    for i in range(3):
        for j in range(3):
            dot = math.fsum(a * b for a, b in zip(columns[i], columns[j]))
            require(abs(dot - int(i == j)) < 1e-5, "PCA eigenvectors not orthonormal")
    values = [float(x) for x in Path(str(prefix) + ".eigenval").read_text().split()]
    require(len(values) == 3 and all(math.isfinite(x) and x > 0 for x in values),
            "PCA expected three finite positive eigenvalues")
    require(values == sorted(values, reverse=True), "PCA eigenvalues not descending")
    if failures:
        for k in ("0", "-1", str(h.data["n"] + 1)):
            h.run("bad_pca_k", ["--pca", "--pca-k", k], error=r"PCA.*(k|sample)")
        h.run("conflicting_standalone", ["--pca", "--ld-prune"], error="conflict")
        h.run("conflicting_lm", ["--pca", "--lm", "--lm-mode", "1", "--pheno", h.data["pheno"]],
              error="conflict")
        h.run("pgen_pca", ["--pca"], error=r"BED|PGEN.*unsupported", pgen=True)


def qc_case(h, failures=True):
    d = h.data
    keep = h.root / "keep.txt"
    remove = h.root / "remove.txt"
    keep.write_text("".join(f"{fid} {iid}\n" for fid, iid in d["samples"][:-2])
                    + f"wrong_fid {d['samples'][-1][1]}\n", encoding="ascii")
    remove.write_text(" ".join(d["samples"][4]) + "\n"
                      + f"wrong_fid {d['samples'][6][1]}\n", encoding="ascii")
    for name, geno, mind, select in (("qc_geno", 0.1, 1.0, False),
                                     ("qc_mind", 1.0, 0.2, False),
                                     ("qc_combined", 0.1, 0.2, True)):
        args = ["--geno", str(geno), "--mind", str(mind)]
        candidates = list(range(d["n"]))
        if select:
            args += ["--keep", keep, "--remove", remove]
            candidates = [i for i in range(d["n"] - 2) if i != 4]
        samples = [i for i in candidates
                   if sum(col[i] is None for col in d["calls"]) / d["m"] <= mind]
        variants = [j for j, col in enumerate(d["calls"])
                    if sum(col[i] is None for i in samples) / len(samples) <= geno]
        require(bool(samples) and bool(variants), "Invalid QC fixture")
        prefix = h.run(name, args)
        check_qc(prefix, d, samples, variants)
    if failures:
        for flag in ("--geno", "--mind"):
            for value in ("-0.1", "1.1", "nan", "inf"):
                h.run("bad_qc_threshold", [flag, value], error=r"geno|mind|finite|range")
        remove_all = h.root / "remove_all.txt"
        remove_all.write_text("".join(f"{fid} {iid}\n" for fid, iid in d["samples"]), encoding="ascii")
        h.run("zero_samples", ["--remove", remove_all], error=r"zero samples|all samples|no samples")
        for suffix in (".qc.bed", ".qc.bim", ".qc.fam"):
            h.run("exclusive_qc", ["--geno", "1"], error=r"exist|exclusive", occupied=suffix)
        h.run("pgen_qc", ["--geno", "0.1"], error=r"BED|PGEN.*unsupported", pgen=True)


def ld_case(h, clump=False, failures=True):
    d = h.data
    normalized = []
    for col in d["calls"]:
        called = [x for x in col if x is not None]
        mean = math.fsum(called) / len(called)
        centered = [0.0 if x is None else x - mean for x in col]
        norm = math.sqrt(math.fsum(x * x for x in centered))
        normalized.append([x / norm for x in centered])

    def linked(i, j):
        ci, _, pi = d["variants"][i]
        cj, _, pj = d["variants"][j]
        return ci == cj and abs(pi - pj) <= 250000 and math.fsum(
            a * b for a, b in zip(normalized[i], normalized[j])) ** 2 >= 0.8

    ids = [v[1] for v in d["variants"]]
    args = ["--ld-clump", d["summary"], "--clump-p", "1"] if clump else ["--ld-prune"]
    prefix = h.run("ld_clump" if clump else "ld_prune", [*args, *LD_ARGS])
    if clump:
        expected = {}
        for i in sorted(range(d["m"]), key=lambda j: (d["pvalues"][j], j)):
            if ids[i] in expected:
                continue
            expected[ids[i]] = ids[i]
            for j in range(d["m"]):
                if ids[j] not in expected and linked(i, j):
                    expected[ids[j]] = ids[i]
        rows = table(Path(str(prefix) + ".clump.tsv"), "SNP P CHR BP LEAD_ID")
        leads = table(Path(str(prefix) + ".lead.tsv"), "SNP P CHR BP LEAD_ID")
        require([r["SNP"] for r in rows] == ids, "Clump membership/order mismatch")
        require({r["SNP"]: r["LEAD_ID"] for r in rows} == expected, "Clump differs from LD reference")
        lead_ids = [snp for snp in ids if expected[snp] == snp]
        require([r["SNP"] for r in leads] == lead_ids, "Clump lead selection/order mismatch")
        require(expected[ids[0]] == ids[1] and ids[1] in lead_ids
                and ids[0] not in lead_ids and ids[d["m"] // 2] in lead_ids,
                "Duplicate LD pair / cross-chromosome lead expectation")
        metadata = {snp: (chrom, bp, d["pvalues"][j])
                    for j, (chrom, snp, bp) in enumerate(d["variants"])}
        for row in rows + leads:
            chrom, bp, p = metadata[row["SNP"]]
            require(row["CHR"] == chrom and int(row["BP"]) == bp
                    and math.isclose(float(row["P"]), p, rel_tol=1e-10), "Clump summary fields changed")
            require(row["LEAD_ID"] == expected[row["SNP"]], "Incorrect lead metadata")
    else:
        kept = []
        for j in range(d["m"]):
            if not any(linked(i, j) for i in kept):
                kept.append(j)
        observed_in = Path(str(prefix) + ".prune.in").read_text().split()
        observed_out = Path(str(prefix) + ".prune.out").read_text().split()
        require(observed_in == [ids[j] for j in kept], "LD prune keep list differs from reference")
        require(observed_out == [ids[j] for j in range(d["m"]) if j not in kept],
                "LD prune removal list differs from reference")
        require(ids[0] in observed_in and ids[1] in observed_out
                and ids[d["m"] // 2] in observed_in, "Duplicate LD pair / chromosome boundary")
    if failures:
        for flag, value in (("--ld-r2", "-0.1"), ("--ld-r2", "1.1"),
                            ("--ld-r2", "nan"), ("--ld-r2", "inf"),
                            ("--ld-step", "0"), ("--ld-window", "-1")):
            h.run("bad_ld_threshold", [*args, flag, value], error=r"ld|window|step")
        if clump:
            h.run("bad_clump_p", ["--ld-clump", d["summary"], "--clump-p", "1.1"], error="clump")
        h.run("pgen_ld", args, error=r"BED|PGEN.*unsupported", pgen=True)


def rare_case(h, method):
    d = h.data
    flag = "--" + method.replace("_", "-")
    prefix = h.run(method, [flag, "--set-file", d["set_file"], "--pheno", d["pheno"]])
    rows = table(Path(str(prefix) + ".rare.txt"),
                 "CHR START END VARIANT_SET_ID METHOD STATISTIC P_VALUE")
    require(len(rows) == len(d["sets"]) and
            {r["VARIANT_SET_ID"] for r in rows} == set(d["sets"]),
            "Missing or duplicate rare-variant sets")
    method_id = "SKAT_O_CANONICAL" if method == "skat_o" else "ACAT_V"
    for row in rows:
        set_id = row["VARIANT_SET_ID"]
        members = [d["variants"][index] for index in d["sets"][set_id]]
        require(row["CHR"] == members[0][0] and int(row["START"]) == min(x[2] for x in members)
                and int(row["END"]) == max(x[2] for x in members),
                f"Wrong coordinates for {set_id}")
        require(row["METHOD"] == method_id, f"Wrong method label for {set_id}")
        statistic = float(row["STATISTIC"])
        require(not math.isnan(statistic) if method == "acat_v" else math.isfinite(statistic),
                f"Invalid statistic for {set_id}")
        probability(row["P_VALUE"], set_id)
    metadata = Path(str(prefix) + ".rare.meta.tsv").read_text(encoding="ascii")
    fields = [line.split("\t") for line in metadata.splitlines() if line]
    require(len(fields) >= 2 and all(len(row) >= 2 for row in fields), "Missing TSV method metadata")
    require(any(method_id in row for row in fields[1:]), f"Missing method label {method_id}")
    require("method" in metadata.lower(), "Rare metadata does not identify its method")
    if method == "skat_o":
        require(all(name in fields[0] for name in ("RHO_GRID", "BEST_RHO", "P_METHOD", "FALLBACK_USED")),
                "Canonical SKAT-O audit metadata is incomplete")
        require("DEPENDENCE_ADJUSTED_GAUSSIAN_QMC" in metadata, "SKAT-O dependence adjustment not reported")


def genetic_model_case(h, failures=True):
    d = h.data
    args = ["--lm", "--lm-mode", "1", "--pheno", d["pheno"]]
    default = h.run("lm_default", args)
    transforms = {"additive": lambda x: x, "dominant": lambda x: int(x >= 1),
                  "recessive": lambda x: int(x == 2), "heterozygote": lambda x: int(x == 1)}
    ymean = math.fsum(d["y"]) / d["n"]
    for model in MODELS:
        prefix = h.run("lm_" + model, [*args, "--snp-model", model])
        path = Path(str(prefix) + ".lm.assoc.txt")
        if model == "additive":
            require(path.read_bytes() == Path(str(default) + ".lm.assoc.txt").read_bytes(),
                    "Default and explicit additive results differ byte-for-byte")
        rows = table(path, "CHR SNP BP A1 A2 FREQ BETA SE P")
        require([r["SNP"] for r in rows] == [v[1] for v in d["variants"]],
                f"Missing/duplicate/reordered LM variants for {model}")
        for j, row in enumerate(rows):
            transformed = [None if x is None else transforms[model](x) for x in d["calls"][j]]
            called = [x for x in transformed if x is not None]
            mean = math.fsum(called) / len(called)
            # Transform observed hardcalls FIRST, then impute their transformed mean.
            x = [mean if value is None else value for value in transformed]
            denominator = math.fsum((value - mean) ** 2 for value in x)
            beta = math.fsum((value - mean) * (y - ymean) for value, y in zip(x, d["y"])) / denominator
            observed = float(row["BETA"])
            require(math.isfinite(observed) and math.isclose(observed, beta, rel_tol=2e-6, abs_tol=2e-7),
                    f"{model} {row['SNP']}: BETA={observed}, intercept OLS reference={beta}")
            require(math.isfinite(float(row["SE"])) and float(row["SE"]) > 0, "Invalid LM SE")
            probability(row["P"], f"{model} {row['SNP']}")
    if failures:
        h.run("invalid_model", [*args, "--snp-model", "unsupported"], error=r"snp.model|model")
        # --snp-model controls marker coding; standard LM/LMM accept non-additive codings.
        for model in MODELS[1:]:
            lmm_prefix = h.run("lmm_" + model, ["--lmm-gemma-mode", "1", "--pheno", d["pheno"],
                                               "--snp-model", model])
            lmm_rows = table(Path(str(lmm_prefix) + ".assoc.txt"))
            require(bool(lmm_rows), f"No LMM rows for {model}")
            for row in lmm_rows:
                require(math.isfinite(float(row["BETA"])), f"invalid LMM BETA for {model}")
                require(math.isfinite(float(row["SE"])), f"invalid LMM SE for {model}")
                require(math.isfinite(float(row["P"])), f"invalid LMM P for {model}")
            for mode in (1, 2, 3, 4):
                prefix = h.run("lm_mode_" + str(mode) + "_" + model,
                               ["--lm", "--lm-mode", str(mode), "--pheno", d["pheno"],
                                "--snp-model", model])
                rows = table(Path(str(prefix) + ".lm.assoc.txt"))
                require(bool(rows), f"No LM rows for {model} mode {mode}")
                for row in rows:
                    require(math.isfinite(float(row["BETA"])), f"invalid BETA for {model} mode {mode}")
                    require(math.isfinite(float(row["SE"])), f"invalid SE for {model} mode {mode}")
                    require(math.isfinite(float(row["P"])), f"invalid P for {model} mode {mode}")
                    if mode == 4:
                        for key in ("P_WALD", "P_LRT", "P_SCORE"):
                            require(math.isfinite(float(row[key])), f"invalid {key} for {model}")
        # --gwas-model dominance/additive_dominance requires additive marker coding.
        h.run("gwas_model_dominance_nonadditive", ["--gwas-model", "dominance", "--pheno", d["pheno"],
                                                   "--snp-model", "dominant"],
              error=r"additive --snp-model")
        h.run("gwas_model_ad_nonadditive", ["--gwas-model", "additive_dominance", "--pheno", d["pheno"],
                                            "--snp-model", "dominant"],
              error=r"additive --snp-model")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", required=True, type=Path, help="CosmicScanner executable")
    parser.add_argument("--case", required=True, choices=CASES)
    parser.add_argument("--out", required=True, type=Path, help="Parent directory for temporary runs")
    args = parser.parse_args()
    binary = args.bin.expanduser().resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        parser.error(f"Not an executable file: {binary}")
    out = args.out.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"association_r2_{args.case}_", dir=out) as directory:
        root = Path(directory)
        data = fixture(root, 150, 500) if args.case == "benchmark" else fixture(root)
        h = Harness(binary, root, data)
        if args.case == "benchmark":
            pca_case(h, failures=False)
            qc_case(h, failures=False)
            ld_case(h, failures=False)
            ld_case(h, clump=True, failures=False)
            genetic_model_case(h, failures=False)
            timing = out / "benchmark.time.tsv"
            with timing.open("w", encoding="ascii", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("command", "n", "m", "seconds", "returncode"))
                writer.writerows((name, n, m, f"{elapsed:.6f}", rc)
                                 for name, n, m, elapsed, rc in h.timings)
            print(f"Timing report (no time threshold): {timing}")
        elif args.case == "pca":
            pca_case(h)
        elif args.case == "qc_filter":
            qc_case(h)
        elif args.case in ("ld_prune", "ld_clump"):
            ld_case(h, clump=args.case == "ld_clump")
        elif args.case in ("skat_o", "acat_v"):
            rare_case(h, args.case)
        else:
            genetic_model_case(h)
        print(f"PASS {args.case}: {len(h.timings)} commands, n={data['n']}, m={data['m']}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError) as exc:
        print(f"FAIL association_r2: {exc}", file=sys.stderr)
        sys.exit(1)
