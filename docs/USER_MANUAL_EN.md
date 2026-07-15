# CosmicScanner v1.1.1 User Manual

## Overview

CosmicScanner v1.1.1 is a focused source release supporting four core workflows
for the first manuscript:

- **Single-trait GWAS**: Linear mixed model association analysis for a single
  continuous trait
- **Multivariate GWAS**: Joint association testing of multiple correlated
  continuous traits
- **Repeatability GWAS**: Mixed-model association analysis for repeated-measure
  phenotypes, with exact record-level analysis
- **A/G/H Matrix Construction**: Pedigree-based A⁻¹, genomic G, and H⁻¹ matrix
  building

All workflows share a unified infrastructure for genotype reading, sample
matching, covariate design, variant filtering, and GRM construction.

## Installation

### System Requirements

- Linux x86_64
- CMake 3.15 or newer
- C++17 compiler (GCC 8+, Clang 10+)
- Eigen 3
- OpenMP
- zlib and zstd

### Building from Source

```bash
git clone --branch v1.1.1 --depth 1 https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner

# Default build (Eigen native solvers)
bash scripts/build.sh

# Optional BLAS/LAPACK acceleration
COSMICSCANNER_ENABLE_BLAS=ON bash scripts/build.sh

# Custom build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cosmicscanner -j$(nproc)
```

If CMake cannot find Eigen3, set CMAKE_PREFIX_PATH to the directory
containing share/eigen3/cmake/Eigen3Config.cmake.

Prebuilt Linux binaries are not currently published. Use the source build above;
the GitHub release provides the archived source corresponding to each tag.

### Build Verification

```bash
bash scripts/validate_release.sh
```

## Input File Formats

### Genotype Data

**PLINK 1 binary (recommended):** `--bfile prefix` for .bed/.bim/.fam.
**PLINK 2 PGEN:** `--pfile prefix` for .pgen/.pvar/.psam.
**BGEN:** `--bgen file.bgen --sample file.sample`.

### Phenotype Files

Single-trait/multivariate: `FID IID TRAIT1 [TRAIT2 ...]`.
Repeatability: `FID IID TIME PHENOTYPE`.
Missing values: NA or -9.

### Covariates

`--qcovar`: quantitative. `--covar`: categorical (dummy-coded).

## CLI Parameters

### Input/Output

| Parameter | Description | Required |
|---|---|---|
| --bfile <prefix> | PLINK 1 binary prefix | One of three |
| --pfile <prefix> | PLINK 2 PGEN prefix | One of three |
| --bgen <file> --sample <f> | BGEN input | Paired |
| --pheno <file> | Phenotype file | Yes |
| --qcovar <file> | Quantitative covariates | No |
| --covar <file> | Categorical covariates | No |
| --out <prefix> | Output prefix | Yes |

### Single-trait GWAS

| Parameter | Description | Default |
|---|---|---|
| --gwas | Enable GWAS scan | Required |
| --scan-tier approx\|semi\|exact | Scan tier | Legacy default |
| --no-loco | Disable LOCO | LOCO enabled |
| --lm <1\|2\|3\|4> | Linear model mode | Off |
| --lmm-gemma-mode <1\|2\|3\|4> | LMM Wald/LRT/Score | Wald |

### Multivariate GWAS

| Parameter | Description | Required |
|---|---|---|
| --multivariate | Enable multivariate | Yes |
| --mpheno-list c1,c2,... | Trait columns | Yes |
| --scan-tier approx\|semi\|exact | Scan tier | Recommended |
| --mv-method <method> | Method (emmv, etc.) | emmv |
| --mv-output-vbeta | Beta covariance | No |

### Repeatability GWAS

| Parameter | Description | Required |
|---|---|---|
| --repeatability-model | Enable model | Yes |
| --tpheno <n> | Time column | Yes |
| --mpheno <n> | Trait column | Yes |
| --repeatability-mode auto\|exact\|fast | Mode | auto |
| --scan-tier approx\|semi\|exact | Scan tier | Recommended |

### A/G/H Matrices

| Parameter | Description | Required |
|---|---|---|
| --make-agh | Build matrices | One of --gwas |
| --ped <file> | Pedigree file | Yes |
| --genotyped-list <file> | Genotyped IDs | No |
| --blend-g <float> | G blend weight | 0.05 |

### Filters

`--maf <x>`, `--hwe <x>`, `--extract <file>`, `--exclude <file>`,
`--chr <chr>`, `--from-bp <pos>`, `--to-bp <pos>`.

### Runtime

`--threads <n>`, `--lmm-null-solver <solver>`, `--reml-maxit <n>`,
`--reml-tol <float>`, `--help`.

## Examples

### Single-trait

```bash
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 \
  --gwas --scan-tier semi --out results/single
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 \
  --gwas --scan-tier semi --no-loco --out results/global
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier approx --out app
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi --out semi
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --out exact
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi \
  --maf 0.01 --hwe 1e-6 --chr 1 --extract keep.txt --threads 8 --out filt
```

### Multivariate

```bash
cosmicscanner --bfile data --pheno traits.txt \
  --mpheno-list 1,2 --multivariate --gwas --scan-tier semi --no-loco \
  --out results/mv
cosmicscanner --bfile data --pheno traits.txt \
  --mpheno-list 1,2,3 --multivariate --gwas --scan-tier semi --no-loco \
  --out results/mv3
cosmicscanner --bfile data --pheno traits.txt \
  --mpheno-list 1,2 --multivariate --gwas --scan-tier exact --no-loco \
  --out results/mv_exact
```

### Repeatability

```bash
cosmicscanner --bfile data --pheno repeated.txt \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact \
  --gwas --scan-tier semi --out results/rep
cosmicscanner --bfile data --pheno repeated.txt \
  --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode fast \
  --gwas --out results/rep_fast
```

### A/G/H

```bash
cosmicscanner --bfile data --make-agh \
  --ped pedigree.txt --out results/agh
cosmicscanner --bfile data --make-agh \
  --ped pedigree.txt \
  --genotyped-list genotyped.txt \
  --blend-g 0.05 --out results/agh
```

## Scan Tiers

--scan-tier approx|semi|exact:

| Tier | Handling | Use case |
|---|---|---|
| approx | Fixed null parameters | High-throughput screening |
| semi | Fixed ratios, profile scale | Recommended main scan |
| exact | Per-SNP full refit | Small data / fine-mapping |

Limitations: multivariate/repeatability exact require --no-loco.
Repeatability exact requires PLINK input, no record weights.

## Output Formats

### assoc.txt (single-trait, repeatability)

Columns: CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P

### assoc.mv.txt (multivariate)

Columns: CHR SNP BP A1 A2 FREQ WALD_STAT P BETA_1 SE_1 BETA_2 SE_2 ...

### repeatability.null.txt

Key-value: n_records, n_subjects, Vg, Vpe, Ve, Vp, h2, r2, timing, scan mode.

### scan.meta.txt

analysis_model, requested_scan_tier, resolved_scan_tier, snp_vc_strategy.

### Matrix .bin/.id

Lower-triangle row-by-row, 32-bit float, native endian. .id = one ID per line.

### agh.meta.txt

n_pedigree, n_genotyped, grm_algorithm, h_tau, h_omega, blend_g,
id_mapping_verified, hinv_symmetry_max_abs.

## GEMMA Compatibility

--lm <1|2|3|4>, --lmm-gemma-mode <1|2|3|4>, --lmm-gemma-refit,
--lmm-null-solver brent (same as GEMMA default).

## Version

v1.1.1. See CITATION.cff for citation metadata.
