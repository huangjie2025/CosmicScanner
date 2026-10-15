#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
import subprocess
import tempfile
import time

from run_association_r2 import fixture


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(binary, data, phenotype, output, rank=0, tolerance=1.0, no_loco=True,
        auto_threshold=None, auto_min=1, auto_max=256):
    command = [str(binary), "--bfile", str(data["prefix"]), "--pheno", str(phenotype),
               "--gwas", "--maf", "0", "--threads", "1",
               "--out", str(output)]
    if no_loco:
        command.append("--no-loco")
    if rank:
        command.extend(["--lmm-low-rank", str(rank),
                        "--lmm-low-rank-oversample", "10",
                        "--lmm-low-rank-seed", "20260911",
                         "--lmm-low-rank-error-tol", str(tolerance)])
    if auto_threshold is not None:
        command.extend(["--fastlmm", "--fastlmm-trace-threshold", str(auto_threshold),
                        "--fastlmm-min-rank", str(auto_min),
                        "--fastlmm-max-rank", str(auto_max),
                        "--lmm-low-rank-seed", "20260911"])
    start = time.perf_counter()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    elapsed = time.perf_counter() - start
    require(result.returncode == 0, "Command failed:\n" + " ".join(command) + "\n" + result.stdout)
    return elapsed


def read_assoc(path):
    with path.open(encoding="ascii") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    result = {}
    for row in rows:
        if row.get("SNP") and row.get("P_A"):
            values = [float(row[key]) for key in ("BETA_A", "SE_A", "P_A")]
            require(all(math.isfinite(value) for value in values), f"Non-finite result: {row}")
            result[row["SNP"]] = row
    return result


def correlation(left, right):
    mean_left = math.fsum(left) / len(left)
    mean_right = math.fsum(right) / len(right)
    numerator = math.fsum((x - mean_left) * (y - mean_right) for x, y in zip(left, right))
    denominator = math.sqrt(math.fsum((x - mean_left) ** 2 for x in left) *
                            math.fsum((y - mean_right) ** 2 for y in right))
    return numerator / denominator


def read_meta(path):
    with path.open(encoding="ascii") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    binary = args.bin.resolve()
    args.out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="low_rank_lmm_", dir=args.out) as directory:
        root = Path(directory)
        data = fixture(root, n=96, m=240)
        phenotype = data["pheno"]
        exact = root / "exact"
        low_rank = root / "low_rank"
        repeat = root / "repeat"
        fallback = root / "fallback"
        loco = root / "loco"
        auto = root / "auto"
        auto_repeat = root / "auto_repeat"
        auto_fallback = root / "auto_fallback"
        exact_seconds = run(binary, data, phenotype, exact)
        low_rank_seconds = run(binary, data, phenotype, low_rank, rank=80)
        run(binary, data, phenotype, repeat, rank=80)
        run(binary, data, phenotype, fallback, rank=40, tolerance=1e-12)
        run(binary, data, phenotype, loco, rank=80, no_loco=False)
        run(binary, data, phenotype, auto, auto_threshold=0.8, auto_min=10, auto_max=95)
        run(binary, data, phenotype, auto_repeat, auto_threshold=0.8, auto_min=10, auto_max=95)
        run(binary, data, phenotype, auto_fallback, auto_threshold=0.99,
            auto_min=5, auto_max=10)

        exact_rows = read_assoc(Path(str(exact) + ".assoc.txt"))
        low_rank_rows = read_assoc(Path(str(low_rank) + ".assoc.txt"))
        require(exact_rows.keys() == low_rank_rows.keys() and len(exact_rows) == 240,
                "Exact and low-rank tested-marker sets differ")
        exact_logp = [-math.log10(max(float(exact_rows[key]["P_A"]), 1e-300)) for key in exact_rows]
        low_rank_logp = [-math.log10(max(float(low_rank_rows[key]["P_A"]), 1e-300)) for key in exact_rows]
        require(correlation(exact_logp, low_rank_logp) > 0.995,
                "Low-rank and exact -log10(P) correlation is below 0.995")
        max_logp_error = max(abs(x - y) for x, y in zip(exact_logp, low_rank_logp))
        require(max_logp_error < 0.1,
                f"Low-rank maximum -log10(P) error exceeds 0.1: {max_logp_error}")
        exact_top = set(sorted(exact_rows, key=lambda key: float(exact_rows[key]["P_A"]))[:10])
        low_rank_top = set(sorted(low_rank_rows, key=lambda key: float(low_rank_rows[key]["P_A"]))[:10])
        require(len(exact_top & low_rank_top) >= 9, "Low-rank Top10 overlap is below 9")
        require(Path(str(low_rank) + ".assoc.txt").read_bytes() ==
                Path(str(repeat) + ".assoc.txt").read_bytes(),
                "Fixed low-rank seed is not byte reproducible")

        metadata = read_meta(Path(str(low_rank) + ".low_rank.meta.tsv"))
        require(len(metadata) == 1 and metadata[0]["FALLBACK"] == "none",
                f"Unexpected low-rank metadata: {metadata}")
        require(metadata[0]["REQUESTED_RANK"] == "80" and
                metadata[0]["EFFECTIVE_RANK"] == "80",
                "Requested/effective rank metadata is incorrect")
        for key in ("RESIDUAL_EIGENVALUE", "RELATIVE_FROBENIUS_ERROR",
                    "DECOMPOSITION_SECONDS"):
            require(math.isfinite(float(metadata[0][key])), f"Non-finite metadata {key}")
        fallback_metadata = read_meta(Path(str(fallback) + ".low_rank.meta.tsv"))
        require(fallback_metadata[0]["FALLBACK"] == "error_tolerance",
                "Strict approximation tolerance did not trigger exact fallback")
        require(Path(str(exact) + ".assoc.txt").read_bytes() ==
                Path(str(fallback) + ".assoc.txt").read_bytes(),
                "Exact fallback changed stable exact LMM output")
        loco_metadata = read_meta(Path(str(loco) + ".low_rank.meta.tsv"))
        require(len(loco_metadata) == 2 and all(row["FALLBACK"] == "none" for row in loco_metadata),
                f"LOCO low-rank metadata is incomplete: {loco_metadata}")
        require(len(read_assoc(Path(str(loco) + ".assoc.txt"))) == 240,
                "LOCO low-rank scan did not return all markers")
        auto_metadata = read_meta(Path(str(auto) + ".low_rank.meta.tsv"))[0]
        selected_rank = int(auto_metadata["SELECTED_RANK"])
        require(auto_metadata["RANK_SELECTION"] == "AUTO_TRACE" and
                auto_metadata["REQUESTED_RANK"] == "95" and
                10 <= selected_rank <= 95 and
                float(auto_metadata["TRACE_FRACTION"]) >= 0.8 and
                auto_metadata["FALLBACK"] == "none",
                f"Automatic trace rank metadata is invalid: {auto_metadata}")
        require(Path(str(auto) + ".assoc.txt").read_bytes() ==
                Path(str(auto_repeat) + ".assoc.txt").read_bytes(),
                "Automatic trace-ranked FastLMM is not deterministic")
        auto_fallback_metadata = read_meta(Path(str(auto_fallback) + ".low_rank.meta.tsv"))[0]
        require(auto_fallback_metadata["FALLBACK"] == "trace_threshold",
                "Unreachable automatic trace target did not trigger exact fallback")
        require(Path(str(exact) + ".assoc.txt").read_bytes() ==
                Path(str(auto_fallback) + ".assoc.txt").read_bytes(),
                "Automatic-rank exact fallback changed exact LMM output")

        scale_root = root / "scale"
        scale_root.mkdir()
        scale_data = fixture(scale_root, n=500, m=700)
        scale_exact = root / "scale_exact"
        scale_low_rank = root / "scale_low_rank"
        scale_exact_seconds = run(binary, scale_data, scale_data["pheno"], scale_exact)
        scale_low_rank_seconds = run(binary, scale_data, scale_data["pheno"],
                                     scale_low_rank, rank=80)
        scale_exact_rows = read_assoc(Path(str(scale_exact) + ".assoc.txt"))
        scale_rows = read_assoc(Path(str(scale_low_rank) + ".assoc.txt"))
        require(len(scale_rows) == 700 and scale_exact_rows.keys() == scale_rows.keys(),
                "Scaled exact and low-rank marker sets differ")
        scale_exact_logp = [-math.log10(max(float(scale_exact_rows[key]["P_A"]), 1e-300))
                            for key in scale_exact_rows]
        scale_low_rank_logp = [-math.log10(max(float(scale_rows[key]["P_A"]), 1e-300))
                               for key in scale_exact_rows]
        scale_correlation = correlation(scale_exact_logp, scale_low_rank_logp)
        require(scale_correlation > 0.95,
                f"Scaled low-rank and exact correlation is below 0.95: {scale_correlation}")
        scale_metadata = read_meta(Path(str(scale_low_rank) + ".low_rank.meta.tsv"))[0]
        low_memory = int(scale_metadata["ESTIMATED_PEAK_BYTES"])
        exact_memory = int(scale_metadata["FULL_EIGEN_ESTIMATED_PEAK_BYTES"])
        require(low_memory < 0.7 * exact_memory,
                f"Low-rank decomposition memory estimate is not materially lower: {low_memory}/{exact_memory}")

        print(f"PASS low-rank LMM validation corr={correlation(exact_logp, low_rank_logp):.6f} "
              f"top10={len(exact_top & low_rank_top)}/10 memory_ratio={low_memory/exact_memory:.3f} "
              f"small_wall_exact={exact_seconds:.3f} small_wall_low_rank={low_rank_seconds:.3f} "
              f"scale_corr={scale_correlation:.6f} scale_wall_exact={scale_exact_seconds:.3f} "
              f"scale_wall_low_rank={scale_low_rank_seconds:.3f}")


if __name__ == "__main__":
    main()
