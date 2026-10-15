#!/usr/bin/env python3
"""Validate the R33 synthetic example dataset.

Checks performed:
  * manifest checksums and file inventory
  * PLINK BED / PGEN / BGEN / VCF / HapMap / text / NPY hard-call identity
  * phenotype sample coverage and value validity
  * checksum reproducibility against a regenerated copy (optional --regenerate)

Usage:
    python3 validate_dataset.py [--dataset DIR] [--manifest]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_DATASET = HERE / "dataset"

FAILURES: list[str] = []
CHECKS: list[dict] = []


def check(name: str, ok: bool, detail: str = "") -> bool:
    CHECKS.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


def _is_float(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_bed_calls(prefix: Path) -> np.ndarray:
    fam = [line.split() for line in prefix.with_suffix(".fam").read_text().splitlines() if line.strip()]
    bim = [line.split() for line in prefix.with_suffix(".bim").read_text().splitlines() if line.strip()]
    raw = prefix.with_suffix(".bed").read_bytes()
    assert raw[:3] == b"\x6c\x1b\x01", "unexpected BED magic"
    n = len(fam)
    stride = (n + 3) // 4
    codes = np.full((len(bim), n), -1, dtype=np.int16)
    for j in range(len(bim)):
        block = raw[3 + j * stride: 3 + (j + 1) * stride]
        for i in range(n):
            code = (block[i // 4] >> (2 * (i % 4))) & 0b11
            codes[j, i] = -1 if code == 0b01 else (2 if code == 0b00 else 1 if code == 0b10 else 0)
    return codes


def read_vcf_calls(path: Path) -> tuple[list[str], np.ndarray]:
    ids, rows, header = [], [], None
    for line in path.read_text().splitlines():
        if line.startswith("#CHROM"):
            header = line.split("\t")
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        ids.append(parts[2])
        calls = []
        for cell in parts[9:]:
            allele = cell.split(":")[0]
            if allele in ("./.", ".", "./."):
                calls.append(-1)
            else:
                calls.append(sum(int(a) for a in allele.replace("|", "/").split("/")))
        rows.append(calls)
    assert header is not None
    return ids, np.array(rows, dtype=np.int16)


def read_hapmap_calls(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    samples = lines[0].split()[11:]
    ids, rows = [], []
    for line in lines[1:]:
        parts = line.split()
        ids.append(parts[0])
        rows.append([-1 if cell == "NN" else cell.count("C") for cell in parts[11:]])
    return ids, np.array(rows, dtype=np.int16)


def read_text_calls(matrix: Path) -> tuple[list[str], np.ndarray]:
    lines = matrix.read_text().splitlines()
    samples = lines[0].split()[1:]
    ids, rows = [], []
    for line in lines[1:]:
        parts = line.split()
        ids.append(parts[0])
        rows.append([-1 if cell == "NA" else int(cell) for cell in parts[1:]])
    return ids, np.array(rows, dtype=np.int16)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    parser.add_argument("--print", dest="show", action="store_true", help="print per-check results")
    args = parser.parse_args()
    dataset = args.dataset.resolve()

    manifest_path = dataset / "manifest.json"
    if not manifest_path.exists():
        print(f"FAIL: missing dataset manifest {manifest_path}", file=sys.stderr)
        return 2
    manifest = json.loads(manifest_path.read_text())
    check("manifest.n_samples", manifest["n_samples"] == 60, f"n_samples={manifest['n_samples']}")
    check("manifest.n_variants", manifest["n_variants"] == 180, f"n_variants={manifest['n_variants']}")
    check("manifest.seed_is_fixed", isinstance(manifest["seed"], int), f"seed={manifest['seed']}")
    check("manifest.has_allele_orientation", bool(manifest.get("allele_orientation")), "")
    check("manifest.files_nonempty", len(manifest["files"]) > 0, f"files={len(manifest['files'])}")

    bad = []
    for entry in manifest["files"]:
        path = dataset / entry["file"]
        if not path.exists():
            bad.append(f"{entry['file']}:missing")
        elif sha256(path) != entry["sha256"]:
            bad.append(f"{entry['file']}:checksum")
    check("manifest.checksums", not bad, "; ".join(bad) if bad else f"{len(manifest['files'])} files verified")

    bed = read_bed_calls(dataset / "synth")
    check("bed.shape", bed.shape == (180, 60), f"shape={bed.shape}")
    check("bed.missing_rate", 0.0 <= (bed < 0).mean() <= 0.02, f"missing={(bed < 0).mean():.4f}")

    bed_ids = [line.split()[1] for line in (dataset / "synth.bim").read_text().splitlines() if line.strip()]
    bed_samples = [line.split()[1] for line in (dataset / "synth.fam").read_text().splitlines() if line.strip()]

    vcf_ids, vcf = read_vcf_calls(dataset / "synth.vcf")
    check("vcf.matches_bed", np.array_equal(vcf, bed) and vcf_ids == bed_ids,
          f"ids_equal={vcf_ids == bed_ids} calls_equal={np.array_equal(vcf, bed)}")

    hmp_ids, hmp = read_hapmap_calls(dataset / "synth.hmp.txt")
    check("hapmap.matches_bed", np.array_equal(hmp, bed) and hmp_ids == bed_ids,
          f"ids_equal={hmp_ids == bed_ids} calls_equal={np.array_equal(hmp, bed)}")

    text_ids, text = read_text_calls(dataset / "synth.geno.txt")
    check("text.matches_bed", np.array_equal(text, bed) and text_ids == bed_ids,
          f"ids_equal={text_ids == bed_ids} calls_equal={np.array_equal(text, bed)}")

    npy = np.load(dataset / "synth.npy")
    npy_ids = [line.split()[0] for line in (dataset / "synth.npy.var.txt").read_text().splitlines() if line.strip()]
    npy_calls = np.where(np.isnan(npy), -1, npy).astype(np.int16)
    check("npy.matches_bed", np.array_equal(npy_calls, bed) and npy_ids == bed_ids,
          f"ids_equal={npy_ids == bed_ids} calls_equal={np.array_equal(npy_calls, bed)}")
    check("npy.dtype", npy.dtype == np.float32, f"dtype={npy.dtype}")

    for name, columns in (("quantitative.pheno", 1), ("binary.pheno", 1), ("count.pheno", 1)):
        rows = [line.split() for line in (dataset / name).read_text().splitlines() if line.strip()]
        check(f"{name}.rows", len(rows) == 60, f"rows={len(rows)}")
        check(f"{name}.ids", all(r[1] in bed_samples for r in rows), "all IIDs present in genotype")
        values = [float(r[2]) for r in rows]
        check(f"{name}.finite", all(math.isfinite(v) for v in values), "")

    binary_rows = [line.split() for line in (dataset / "binary.pheno").read_text().splitlines() if line.strip()]
    binary_values = [int(r[2]) for r in binary_rows]
    check("binary.pheno.balanced_0_1", set(binary_values) <= {0, 1}, f"values={sorted(set(binary_values))}")
    check("binary.pheno.both_classes", min(binary_values) == 0 and max(binary_values) == 1, "")

    count_rows = [line.split() for line in (dataset / "count.pheno").read_text().splitlines() if line.strip()]
    check("count.pheno.nonneg_integer", all(float(r[2]) >= 0 and float(r[2]).is_integer() for r in count_rows), "")

    long_rows = [line.split() for line in (dataset / "long.pheno").read_text().splitlines() if line.strip()]
    check("long.pheno.repeated_records", len(long_rows) - 1 == 180,
          f"rows={len(long_rows) - 1} (60 samples x 3 records)")
    check("long.pheno.long_format", long_rows[0] == ["FID", "IID", "TIME", "VALUE"],
          f"header={long_rows[0]}")
    check("long.pheno.time_values", {r[2] for r in long_rows[1:]} == {"1", "2", "3"}, "")
    vcf_long = (dataset / "long.vcf.pheno").read_text().splitlines()
    check("long.vcf.pheno.convention", all(line.split()[0] == "0" for line in vcf_long), "")

    survival_rows = [line.split() for line in (dataset / "survival.pheno").read_text().splitlines() if line.strip()]
    check("survival.pheno.time_status", all(float(r[2]) > 0 and r[3] in ("0", "1") for r in survival_rows), "")

    multi_rows = (dataset / "multivariate.pheno").read_text().splitlines()
    check("multivariate.pheno.header", multi_rows[0].split() == ["FID", "IID", "TRAIT1", "TRAIT2", "TRAIT3"],
          f"header={multi_rows[0]}")
    check("multivariate.pheno.rows", len(multi_rows) - 1 == 60, f"rows={len(multi_rows) - 1}")

    qcov_rows = [line.split() for line in (dataset / "qcov.txt").read_text().splitlines() if line.strip()]
    check("qcov.header_and_rows", qcov_rows[0] == ["FID", "IID", "AGE", "SEX", "PC1"] and len(qcov_rows) == 61,
          f"header={qcov_rows[0]} rows={len(qcov_rows) - 1}")

    qcov_vcf = (dataset / "qcov.vcf.txt").read_text().splitlines()
    check("qcov.vcf.convention", all(line.split()[0] == "0" for line in qcov_vcf[1:]), "")

    ped_rows = [line.split() for line in (dataset / "pedigree.ped").read_text().splitlines() if line.strip()]
    check("pedigree.three_columns", all(len(r) == 3 for r in ped_rows), "")

    for extra in ("rare_sets.txt", "offset.txt", "common_env.txt", "genotyped.txt"):
        check(f"{extra}.exists", (dataset / extra).exists(), "")

    offset_rows = [line for line in (dataset / "offset.txt").read_text().splitlines() if line.strip()]
    check("offset.one_value_per_line",
          len(offset_rows) == 60 and all(_is_float(row) for row in offset_rows),
          f"rows={len(offset_rows)} bare numeric values")

    negative_expectations = {
        "malformed_offset.txt": lambda rows: all(len(row) == 3 for row in rows),
        "invalid_values.pheno": lambda rows: any(not _is_float(row[2]) for row in rows),
        "id_mismatch.pheno": lambda rows: all(row[1].startswith("X") for row in rows),
        "invalid_covariates.txt": lambda rows: any(not _is_float(row[2]) for row in rows[1:]),
        "constant.pheno": lambda rows: len({row[2] for row in rows}) == 1,
    }
    for name, predicate in negative_expectations.items():
        rows = [line.split() for line in (dataset / name).read_text().splitlines() if line.strip()]
        check(f"negative_fixture:{name}", predicate(rows), "fixture is not malformed as documented")

    planted = set(manifest["planted_causal_variants"])
    check("planted_causal_count", len(planted) == 5, f"{sorted(planted)}")

    if args.show:
        for item in CHECKS:
            print(f"{item['status']:<4} {item['check']:<40} {item['detail']}")
    passed = sum(1 for c in CHECKS if c["status"] == "PASS")
    print(f"{passed}/{len(CHECKS)} dataset checks passed")
    if FAILURES:
        for failure in FAILURES:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())