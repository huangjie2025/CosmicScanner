#!/usr/bin/env python3
"""Generate the fully synthetic R33 example dataset.

Everything produced here is machine-generated from a fixed seed. No third-party
human data, no consent/privacy question, and no absolute paths are involved.

Usage:
    python3 generate_dataset.py [--out DIR] [--seed 20260301] [--force]

The PLINK binary files are the canonical source of truth. The PGEN and BGEN
copies are produced with plink2 (auto-detected) so that the example dataset
matches the on-disk formats used by CosmicScanner users.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

SEED = 20260301
N_SAMPLES = 60
N_VARIANTS = 180
SPLIT = 90  # chr1: rs00001..rs00090, chr2: rs00091..rs00180

# Planted additive QTL (variant index, effect on the primary quantitative trait).
CAUSAL = ((4, 1.50), (23, 1.20), (44, -1.00), (95, 1.30), (130, 0.90))

MISSING_RATE = 0.005
MISSING_CODE = 0b01


def find_plink2() -> str | None:
    candidate = shutil.which("plink2")
    if candidate:
        return candidate
    env = Path.home() / "software" / "miniconda3" / "envs" / "benchmark" / "bin" / "plink2"
    return str(env) if env.exists() else None


def sample_ids() -> list[str]:
    return [f"S{i + 1:03d}" for i in range(N_SAMPLES)]


def variant_table() -> list[tuple[str, str, int, str, str]]:
    variants = []
    for j in range(N_VARIANTS):
        chrom = "1" if j < SPLIT else "2"
        variants.append((chrom, f"rs{j + 1:05d}", 10_000 + (j % SPLIT) * 5_000, "A", "C"))
    return variants


def draw_genotypes(rng: np.random.Generator) -> np.ndarray:
    """Return an (M, N) matrix of hard calls in {0, 1, 2}, -1 marks missing."""
    freq = rng.uniform(0.12, 0.48, size=(1, N_VARIANTS))
    dosage = rng.binomial(2, freq, size=(N_SAMPLES, N_VARIANTS)).astype(np.int16)
    missing = rng.random((N_SAMPLES, N_VARIANTS)) < MISSING_RATE
    dosage[missing] = -1
    return dosage


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    spread = values.std()
    return (values - values.mean()) / (spread if spread > 0 else 1.0)


def build_traits(genotypes: np.ndarray, rng: np.random.Generator) -> dict[str, np.ndarray]:
    age = rng.normal(35, 8, size=N_SAMPLES)
    sex = rng.integers(0, 2, size=N_SAMPLES)
    age = standardize(age)
    poly = np.column_stack([genotypes[:, index].astype(float) for index, _ in CAUSAL])
    poly[genotypes[:, [index for index, _ in CAUSAL]] == -1] = np.nan
    poly = np.nan_to_num(poly)

    genetic = np.zeros(N_SAMPLES)
    for (index, effect), column in zip(CAUSAL, poly.T):
        genetic += effect * column
    residual = rng.normal(0.0, 0.5, size=N_SAMPLES)
    quantitative = genetic + 0.25 * age + 0.30 * sex + residual

    pc_source = poly[:, 0] + poly[:, 1] + poly[:, 2]
    covariates = np.column_stack([age, sex.astype(float), standardize(pc_source)])

    liability = standardize(quantitative)
    binary = (liability > np.median(liability)).astype(int)

    count_mean = np.exp(0.45 + 0.30 * poly[:, 0] + 0.05 * age)
    count = rng.poisson(count_mean)

    trait2 = 0.60 * standardize(quantitative) + 0.45 * poly[:, 3] + rng.normal(0, 0.6, N_SAMPLES)
    trait3 = 0.50 * standardize(quantitative) - 0.80 * poly[:, 2] + rng.normal(0, 0.6, N_SAMPLES)

    permanent = rng.normal(0.0, 0.8, size=N_SAMPLES)
    record_residual = rng.normal(0.0, 0.4, size=(N_SAMPLES, 3))
    long_values = (np.repeat(genetic[:, None], 3, axis=1) + permanent[:, None] + record_residual).reshape(-1)

    ordinal = np.digitize(standardize(quantitative), [-0.35, 0.55])
    survival_time = np.round(np.exp(3.0 + 0.25 * standardize(quantitative) + rng.normal(0, 0.2, N_SAMPLES)), 3)
    survival_status = (rng.random(N_SAMPLES) < 0.65).astype(int)

    return {
        "age": age,
        "sex": sex,
        "covariates": covariates,
        "quantitative": standardize(quantitative),
        "binary": binary,
        "count": count,
        "trait2": standardize(trait2),
        "trait3": standardize(trait3),
        "long": long_values,
        "ordinal": ordinal,
        "survival_time": survival_time,
        "survival_status": survival_status,
        "genetic_score": genetic,
    }


def write_bed(prefix: Path, samples: list[str], variants, genotypes: np.ndarray) -> None:
    with prefix.with_suffix(".fam").open("w", encoding="ascii") as fam:
        for sample in samples:
            fam.write(f"{sample} {sample} 0 0 {1 if int(sample[1:]) % 2 else 2} -9\n")
    with prefix.with_suffix(".bim").open("w", encoding="ascii") as bim:
        for chrom, name, position, allele1, allele2 in variants:
            bim.write(f"{chrom}\t{name}\t0\t{position}\t{allele1}\t{allele2}\n")
    bed = bytearray(b"\x6c\x1b\x01")
    bytes_per_variant = (N_SAMPLES + 3) // 4
    for row in genotypes.T:
        packed = bytearray(bytes_per_variant)
        for i, call in enumerate(row):
            code = MISSING_CODE if call < 0 else (0b11, 0b10, 0b00)[int(call)]
            packed[i // 4] |= code << (2 * (i % 4))
        bed.extend(packed)
    prefix.with_suffix(".bed").write_bytes(bytes(bed))


def write_vcf(path: Path, samples: list[str], variants, genotypes: np.ndarray) -> None:
    with path.open("w", encoding="ascii") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=CosmicScanner-R33-synthetic\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for j, (chrom, name, position, ref, alt) in enumerate(variants):
            calls = []
            for call in genotypes[:, j]:
                if call < 0:
                    calls.append("./.")
                else:
                    calls.append(("0/0", "0/1", "1/1")[int(call)])
            vcf.write(f"{chrom}\t{position}\t{name}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")


def write_hapmap(path: Path, samples: list[str], variants, genotypes: np.ndarray) -> None:
    with path.open("w", encoding="ascii") as target:
        target.write("rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode "
                     + " ".join(samples) + "\n")
        for j, (chrom, name, position, ref, alt) in enumerate(variants):
            encoded = []
            for call in genotypes[:, j]:
                if call < 0:
                    encoded.append("NN")
                else:
                    encoded.append(("AA", "AC", "CC")[int(call)])
            target.write(f"{name} {ref}/{alt} {chrom} {position} + NA NA NA NA NA NA " + " ".join(encoded) + "\n")


def write_text_genotype(path: Path, var_path: Path, samples, variants, genotypes) -> None:
    with path.open("w", encoding="ascii") as target:
        target.write("VariantID " + " ".join(samples) + "\n")
        for j, (_, name, _, _, _) in enumerate(variants):
            calls = ["NA" if call < 0 else str(int(call)) for call in genotypes[:, j]]
            target.write(f"{name} " + " ".join(calls) + "\n")
    with var_path.open("w", encoding="ascii") as target:
        target.write("ID CHROM POS REF ALT\n")
        for chrom, name, position, ref, alt in variants:
            target.write(f"{name} {chrom} {position} {ref} {alt}\n")


def write_npy(path: Path, var_path: Path, sample_path: Path, samples, variants, genotypes) -> None:
    matrix = np.where(genotypes < 0, np.nan, genotypes).astype(np.float32).T.copy()
    np.save(path, matrix)
    with var_path.open("w", encoding="ascii") as target:
        for chrom, name, position, ref, alt in variants:
            target.write(f"{name}\t{chrom}\t{position}\t{ref}\t{alt}\n")
    with sample_path.open("w", encoding="ascii") as target:
        for sample in samples:
            target.write(f"{sample}\t{sample}\n")


def write_phenotypes(directory: Path, samples: list[str], traits: dict) -> None:
    def rows(values, decimals=6):
        for sample, value in zip(samples, values):
            yield f"{sample} {sample} {value:.{decimals}f}\n"

    (directory / "quantitative.pheno").write_text("".join(rows(traits["quantitative"])), encoding="ascii")
    (directory / "binary.pheno").write_text(
        "".join(f"{s} {s} {int(v)}\n" for s, v in zip(samples, traits["binary"])), encoding="ascii")
    (directory / "count.pheno").write_text(
        "".join(f"{s} {s} {int(v)}\n" for s, v in zip(samples, traits["count"])), encoding="ascii")
    (directory / "ordinal.pheno").write_text(
        "".join(f"{s} {s} {int(v)}\n" for s, v in zip(samples, traits["ordinal"])), encoding="ascii")
    (directory / "survival.pheno").write_text(
        "".join(f"{s} {s} {t:.3f} {int(st)}\n"
                for s, t, st in zip(samples, traits["survival_time"], traits["survival_status"])), encoding="ascii")
    (directory / "multivariate.pheno").write_text(
        "FID\tIID\tTRAIT1\tTRAIT2\tTRAIT3\n" + "".join(
            f"{s}\t{s}\t{traits['quantitative'][i]:.6f}\t{traits['trait2'][i]:.6f}\t{traits['trait3'][i]:.6f}\n"
            for i, s in enumerate(samples)), encoding="ascii")
    (directory / "long.pheno").write_text(
        "FID\tIID\tTIME\tVALUE\n" + "".join(
            f"{s}\t{s}\t{record + 1}\t{traits['long'][3 * i + record]:.6f}\n"
            for i, s in enumerate(samples) for record in range(3)), encoding="ascii")
    (directory / "long.vcf.pheno").write_text(
        "".join(f"0 {s} {record + 1} {traits['long'][3 * i + record]:.6f}\n"
                for i, s in enumerate(samples) for record in range(3)), encoding="ascii")
    (directory / "multivariate.vcf.pheno").write_text(
        "FID\tIID\tTRAIT1\tTRAIT2\tTRAIT3\n" + "".join(
            f"0\t{s}\t{traits['quantitative'][i]:.6f}\t{traits['trait2'][i]:.6f}\t{traits['trait3'][i]:.6f}\n"
            for i, s in enumerate(samples)), encoding="ascii")
    # VCF phenotypes use the 0/FID convention used by the reader's VCF sample identity.
    (directory / "quantitative.vcf.pheno").write_text(
        "".join(f"0 {s} {v:.6f}\n" for s, v in zip(samples, traits["quantitative"])), encoding="ascii")
    (directory / "binary.vcf.pheno").write_text(
        "".join(f"0 {s} {int(v)}\n" for s, v in zip(samples, traits["binary"])), encoding="ascii")
    with (directory / "qcov.txt").open("w", encoding="ascii") as target:
        target.write("FID\tIID\tAGE\tSEX\tPC1\n")
        for i, sample in enumerate(samples):
            age, sex, pc1 = traits["covariates"][i]
            target.write(f"{sample}\t{sample}\t{age:.6f}\t{sex:.0f}\t{pc1:.6f}\n")
    with (directory / "qcov_nosex.txt").open("w", encoding="ascii") as target:
        target.write("FID\tIID\tAGE\tPC1\n")
        for i, sample in enumerate(samples):
            target.write(f"{sample}\t{sample}\t{traits['covariates'][i][0]:.6f}\t{traits['covariates'][i][2]:.6f}\n")
    (directory / "qcov.vcf.txt").write_text(
        "FID\tIID\tAGE\tSEX\tPC1\n" + "".join(
            f"0\t{s}\t{traits['covariates'][i][0]:.6f}\t{traits['covariates'][i][1]:.0f}\t{traits['covariates'][i][2]:.6f}\n"
            for i, s in enumerate(samples)), encoding="ascii")
    (directory / "qcov_nosex.vcf.txt").write_text(
        "FID\tIID\tAGE\tPC1\n" + "".join(
            f"0\t{s}\t{traits['covariates'][i][0]:.6f}\t{traits['covariates'][i][2]:.6f}\n"
            for i, s in enumerate(samples)), encoding="ascii")
    # CountModel::loadOffset reads exactly one numeric value per line, in genotype order.
    (directory / "offset.txt").write_text(
        "".join(f"{4.0 + 0.001 * i:.3f}\n" for i in range(len(samples))), encoding="ascii")
    (directory / "malformed_offset.txt").write_text(
        "".join(f"{s}\t{s}\t{4.0:.3f}\n" for s in samples), encoding="ascii")
    with (directory / "common_env.txt").open("w", encoding="ascii") as target:
        for i, sample in enumerate(samples):
            target.write(f"{sample} {sample} PEN{1 + i % 3}\n")
    with (directory / "pedigree.ped").open("w", encoding="ascii") as target:
        for i, sample in enumerate(samples):
            sire = samples[(i + 7) % N_SAMPLES]
            dam = samples[(i + 19) % N_SAMPLES]
            target.write(f"{sample} {sire} {dam}\n")
    (directory / "genotyped.txt").write_text("".join(f"{s}\n" for s in samples), encoding="ascii")


def write_negative_fixtures(directory: Path, samples: list[str]) -> None:
    """Intentionally invalid inputs used by the negative-test section of the suite."""
    (directory / "invalid_values.pheno").write_text(
        "".join(f"{s} {s} {'not_a_number' if i == 3 else f'{1.0 + 0.05 * i:.4f}'}\n"
                for i, s in enumerate(samples)),
        encoding="ascii")
    (directory / "id_mismatch.pheno").write_text(
        "".join(f"X{i + 1:03d} X{i + 1:03d} 1.0\n" for i in range(len(samples))), encoding="ascii")
    (directory / "invalid_covariates.txt").write_text(
        "FID\tIID\tAGE\tSEX\tPC1\n" + "".join(
            f"{s}\t{s}\t{'not_a_number' if i == 2 else '1.0'}\t1\t0.2\n" for i, s in enumerate(samples)),
        encoding="ascii")
    (directory / "constant.pheno").write_text(
        "".join(f"{s} {s} 2.0\n" for s in samples), encoding="ascii")


def write_sets(path: Path, variants) -> None:
    """Set files are parsed as `SET_ID<whitespace>SNP_ID` one SNP per line.

    Verified against SetFileReader::read in src/rare_variant_test.cpp: each line
    contributes exactly one SNP to the set named in the first column.
    """
    with path.open("w", encoding="ascii") as target:
        for set_id, span in (("SET_A", range(0, 30)), ("SET_B", range(40, 70)), ("SET_C", range(100, 130))):
            for i in span:
                target.write(f"{set_id}\t{variants[i][1]}\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "dataset")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--force", action="store_true", help="regenerate even when data already exists")
    args = parser.parse_args()

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    plink_prefix = out / "synth"

    if (out / "manifest.json").exists() and not args.force:
        print(f"dataset already present: {out} (use --force to regenerate)")
        return 0

    rng = np.random.default_rng(args.seed)
    samples = sample_ids()
    variants = variant_table()
    genotypes = draw_genotypes(rng)
    traits = build_traits(genotypes, rng)

    write_bed(plink_prefix, samples, variants, genotypes)
    write_vcf(out / "synth.vcf", samples, variants, genotypes)
    write_hapmap(out / "synth.hmp.txt", samples, variants, genotypes)
    write_text_genotype(out / "synth.geno.txt", out / "synth.var.txt", samples, variants, genotypes)
    write_npy(out / "synth.npy", out / "synth.npy.var.txt", out / "synth.npy.sample.txt", samples, variants, genotypes)
    write_phenotypes(out, samples, traits)
    write_sets(out / "rare_sets.txt", variants)
    write_negative_fixtures(out, samples)

    plink2 = find_plink2()
    if plink2 is None:
        print("WARNING: plink2 not found; PGEN/BGEN mirrors were not generated", file=sys.stderr)
    else:
        subprocess.run([plink2, "--bfile", str(plink_prefix), "--make-pgen", "--out", str(out / "synth")],
                       check=True, stdout=subprocess.DEVNULL)
        subprocess.run([plink2, "--bfile", str(plink_prefix), "--export", "bgen-1.2", "--out", str(out / "synth")],
                       check=True, stdout=subprocess.DEVNULL)

    entries = []
    purpose = {
        "synth.bed": "PLINK 1.x binary hard calls (canonical example genotype)",
        "synth.bim": "PLINK 1.x variant metadata",
        "synth.fam": "PLINK 1.x sample metadata",
        "synth.pgen": "PLINK 2.0 hard-call mirror of the canonical genotype",
        "synth.pvar": "PLINK 2.0 variant metadata",
        "synth.psam": "PLINK 2.0 sample metadata",
        "synth.bgen": "BGEN Layout 2 (zlib) probability mirror of the canonical genotype",
        "synth.sample": "Oxford sample file for the BGEN mirror",
        "synth.vcf": "VCFv4.2 hard-call mirror of the canonical genotype",
        "synth.hmp.txt": "HapMap text mirror of the canonical genotype",
        "synth.geno.txt": "Text genotype matrix (variants x samples)",
        "synth.var.txt": "Text genotype variant metadata",
        "synth.npy": "NPY float32 genotype matrix (variants x samples)",
        "synth.npy.var.txt": "NPY variant metadata",
        "synth.npy.sample.txt": "NPY sample metadata",
        "quantitative.pheno": "Continuous trait with five planted additive QTL",
        "binary.pheno": "Balanced binary trait derived from the quantitative liability",
        "count.pheno": "Integer count trait with offset support",
        "ordinal.pheno": "Three-level ordinal trait (experimental route)",
        "survival.pheno": "FID IID TIME STATUS survival trait (experimental route)",
        "multivariate.pheno": "FID IID TRAIT1 TRAIT2 TRAIT3 for mvLMM",
        "long.pheno": "FID IID TIME VALUE with three repeated records per individual (repeatability / RRM)",
        "long.vcf.pheno": "Long phenotype in the 0/FID convention used for VCF input",
        "multivariate.vcf.pheno": "Multivariate phenotype in the 0/FID convention used for VCF input",
        "quantitative.vcf.pheno": "Quantitative trait in the 0/FID convention used for VCF input",
        "binary.vcf.pheno": "Binary trait in the 0/FID convention used for VCF input",
        "qcov.txt": "AGE/SEX/PC1 quantitative covariates",
        "qcov_nosex.txt": "AGE/PC1 covariates without sex",
        "qcov.vcf.txt": "Covariates in the 0/FID convention used for VCF input",
        "qcov_nosex.vcf.txt": "AGE/PC1 covariates without sex, 0/FID convention for VCF input",
        "offset.txt": "Constant log-exposure offsets, one numeric value per line (count model)",
        "malformed_offset.txt": "NEGATIVE fixture: offset file with FID/IID columns instead of bare values",
        "common_env.txt": "FID IID GROUP common-environment assignment",
        "pedigree.ped": "ID SIRE DAM pedigree for pedigree-based examples",
        "genotyped.txt": "Genotyped-ID list used by pedigree-based examples",
        "rare_sets.txt": "Variant sets (one SNP per line) for SKAT/Burden/SKAT-O/ACAT-V examples",
        "invalid_values.pheno": "NEGATIVE fixture: one non-numeric phenotype value",
        "id_mismatch.pheno": "NEGATIVE fixture: sample IDs absent from the genotype",
        "invalid_covariates.txt": "NEGATIVE fixture: one non-numeric covariate value",
        "constant.pheno": "NEGATIVE fixture: zero-variance phenotype (degenerate information)",
    }
    for name in sorted(purpose):
        path = out / name
        if path.exists():
            entries.append({
                "file": name,
                "format": path.suffix.lstrip(".") or "text",
                "n_samples": N_SAMPLES if name.split(".")[0] in {"synth"} else None,
                "n_variants": N_VARIANTS if name.split(".")[0] in {"synth"} else None,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "purpose": purpose[name],
            })

    manifest = {
        "dataset": "CosmicScanner R33 synthetic example dataset",
        "source": "100% synthetic, generated by examples/data/r33_synthetic/generate_dataset.py",
        "generator": "generate_dataset.py",
        "seed": args.seed,
        "n_samples": N_SAMPLES,
        "n_variants": N_VARIANTS,
        "chromosomes": ["1", "2"],
        "allele_orientation": "BIM/VCF REF=A ALT=C; dosage counts the ALT(C) allele; PLINK BED code maps 0->AAB,1->AB,2->BB",
        "missing_encoding": "PLINK 0b01, VCF ./., HapMap NN, text NA, NPY NaN",
        "planted_causal_variants": [variants[index][1] for index, _ in CAUSAL],
        "planted_effects": {variants[index][1]: effect for index, effect in CAUSAL},
        "plink2_used_for": ["synth.pgen", "synth.bgen"],
        "files": entries,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    print(f"generated {len(entries)} files under {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())