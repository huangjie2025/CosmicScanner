# CosmicScanner Official User Guide (Tutorial)

**Author**: Jie Huang
**Development-tree documentation version**: 1.1.0
**GitHub**: [https://github.com/huangjie2025/CosmicScanner](https://github.com/huangjie2025/CosmicScanner)
**Updated**: 2026-08-24

> This guide describes extended capabilities in the development tree. Paper
> inclusion and software maturity are separate. The release-candidate maturity
> labels are `Stable`, `Experimental validated`, and
> `Experimental approximate`; see the
> [paper model support matrix](PAPER_MODEL_SUPPORT_MATRIX.md).

---

## Table of Contents

- [1. Introduction](#1-introduction)
- [2. Core Features](#2-core-features)
- [3. Installation](#3-installation)
- [4. Data Preparation](#4-data-preparation)
- [5. Quick Start](#5-quick-start)
- [6. Detailed Usage & Examples](#6-detailed-usage--examples)
  - [6.1 Single-Trait LMM GWAS](#61-single-trait-lmm-gwas)
  - [6.2 Multivariate LMM GWAS](#62-multivariate-lmm-gwas)
  - [6.3 Repeatability Model GWAS](#63-repeatability-model-gwas)
  - [6.4 Random Regression Model (RRM) GWAS](#64-random-regression-model-rrm-gwas)
  - [6.5 Maternal Effects Model GWAS](#65-maternal-effects-model-gwas)
  - [6.6 Single-Step GWAS (ssGWAS)](#66-single-step-gwas-ssgwas)
  - [6.7 Threshold/Probit Model GWAS](#67-thresholdprobit-model-gwas)
  - [6.8 Count/Binomial PSR-GWAS](#68-countbinomial-psr-gwas)
  - [6.9 Ordinal PSR-GWAS](#69-ordinal-psr-gwas)
  - [6.10 Survival Cox Model GWAS](#610-survival-cox-model-gwas)
  - [6.11 GxE Reaction Norm Model GWAS](#611-gxe-reaction-norm-model-gwas)
  - [6.12 Rare Variant Analysis (SKAT / Burden)](#612-rare-variant-analysis-skat--burden)
  - [6.13 Biobank-Scale Analysis](#613-biobank-scale-analysis)
  - [6.14 Empirical Bayes EB-GWAS](#614-empirical-bayes-eb-gwas)
  - [6.15 Common-Environment Model GWAS](#615-common-environment-model-gwas-experimental)
  - [6.16 Pedigree Animal Model GWAS](#616-pedigree-animal-model-gwas)
- [7. Output File Formats](#7-output-file-formats)
- [8. Complete Parameter Reference](#8-complete-parameter-reference)
- [9. Performance Tuning](#9-performance-tuning)
- [10. FAQ](#10-faq)
- [11. Related Documentation](#11-related-documentation)

---

## 1. Introduction

**CosmicScanner** is a high-performance computational tool built for Genome-Wide Association Studies (GWAS) and is a core component of the **Cosmic** pan-omics analysis suite. It shatters traditional bottlenecks in computing speed, memory consumption, and model versatility, easily handling ultra-large cohort genotype data and complex phenotypes.

**Supported Platforms**:
- **Linux**: Native support for major Linux distributions (x86_64 architecture).
- **Windows**: Recommended to run via WSL2 (Windows Subsystem for Linux).

---

## 2. Core Features

### GWAS Model Suite (Core Models and Compatibility Routes)

| # | Model | Use Case | Key Parameters |
|---|-------|----------|---------------|
| 1 | Single-Trait LMM GWAS | Standard single-phenotype association; GEMMA `-lmm` mode/refit output | `--gwas` / `--lmm-gemma-mode` / `--lmm-gemma-refit` |
| 2 | GEMMA-style LM GWAS | Ordinary linear model without random effects; GEMMA `-lm` parity | `--lm --lm-mode 1\|2\|3\|4` |
| 3 | Multivariate LMM GWAS | Multi-phenotype joint testing | `--multivariate --mv-method` |
| 4 | Repeatability Model GWAS | Multiple measurements per individual | `--repeatability-model` |
| 5 | Random Regression (RRM) GWAS | Longitudinal/dynamic traits | `--long-pheno --rrm-order` |
| 6 | Maternal Effects Model GWAS | Birth/weaning weight traits | `--maternal-model` |
| 7 | Single-Step GWAS (ssGWAS) | Pedigree + genomic joint analysis | `--single-step --ped` |
| 8 | Threshold/Probit Model GWAS | Binary traits | `--threshold-model --link` |
| 9 | Count/Binomial PSR-GWAS | Poisson/Binomial traits | `--count-model --family` |
| 10 | Ordinal PSR-GWAS | 3+ ordered category traits | `--ordinal-model` |
| 11 | Survival Cox Model GWAS | Time-to-event traits | `--survival-model` |
| 12 | GxE Reaction Norm GWAS | Environment gradient interaction | `--gxe-model --env-covar` |
| 13 | GEMMA-style utilities | VC, summary-stat VC/CI, GRM eigen, SNP correlation, effect prediction, multi-locus outputs | `--vc` / `--ci` / `--gs` / `--gq` / `--beta` / `--cor` / `--cat` / `--mcat` / `--wcat` / `--wsnp` / `--study` / `--ref` / `--mstudy` / `--mref` / `--gemma-eigen` / `--calccor` / `--predict-effect` / `--bslmm` |
| 14 | Common-Environment GWAS (Experimental) | Repeated records with group-level random intercept | `--common-env-model --common-env-file` |

### Release-Candidate Status for Paper Models

| Model/path | Software status |
|---|---|
| General LMM; multivariate LMM | Stable |
| Repeatability; RRM Method A; animal; maternal; common environment | Experimental validated |
| H-matrix single-step GWAS (Method 1) | Experimental validated |
| ssGBLUP back-solving GWAS (Method 2) | Experimental validated |
| Exact back-solved SNP P-value | Stable |
| APY-exact back-solved SNP P-value | Experimental validated |
| APY-approx back-solved SNP P-value | Experimental approximate |

Method 1's fixed-SNP coefficient and Method 2's back-solved marker effect are
different estimands and must not be reported under an unqualified `ssGWAS`
label. APY-approx raw P-values do not guarantee strict type-I error control;
CosmicScanner applies no genomic control or P-value rescaling to this route.

### Engine & Infrastructure

| # | Feature | Key Parameters |
|---|---------|---------------|
| 14 | Rare Variant Analysis (SKAT/Burden/SPA) | `--skat` / `--burden` / `--spa` |
| 14 | Biobank 2-step (REGENIE-style) | `--biobank-mode` |
| 15 | IDUL/AIDUL Accelerated Solvers | `--lmm-null-solver` |
| 16 | Input Adaptation Layer MVP | `--input-file --input-model` |
| 17 | Ultra-Fast GRM Construction | `--make-grm --make-grm-alg` |
| 18 | Empirical Bayes EB-GWAS (validated auxiliary route) | `--model lmm-eb` |

---

## 3. Installation

### 3.1 Pre-compiled Binaries (Recommended)

After publication, the binary will be available from the corresponding Release page. The current local release-candidate package is:

```bash
tar -zxvf CosmicGWAS-v1.1.0-linux-x86_64.tar.gz
cd CosmicGWAS-v1.1.0-linux-x86_64
bash install.sh
bin/cosmicscanner --help
```

### 3.2 Compiling from Source

Ensure `CMake (>= 3.15)` and a C++17 compatible compiler (e.g., GCC 9+) are installed:

```bash
git clone --recurse-submodules https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
```

The executable will be located at `build/cosmicscanner`.

### 3.3 Verify Installation

```bash
./cosmicscanner --help
```

If the full parameter list is displayed, installation is successful.

---

## 4. Data Preparation

### 4.1 Genotype Files

| Format | Parameter | Files |
|--------|-----------|-------|
| PLINK 1.9 | `--bfile <prefix>` | `.bed` / `.bim` / `.fam` |
| PLINK 2.0 | `--pfile <prefix>` or `--pgen`/`--pvar`/`--psam` | `.pgen` / `.pvar` / `.psam` |
| BGEN | `--bgen <file> --sample <file>` | `.bgen` / `.sample` (Oxford) |

**Variant filtering**:
- `--maf <float>`: Minimum MAF filter (default: 0.01)
- `--hwe <float>`: HWE exact test p-value filter (default: 0.0, disabled)
- `--extract <file>`: Include only listed SNP IDs
- `--exclude <file>`: Exclude listed SNP IDs
- `--chr <chr>`: Restrict scan to one chromosome
- `--from-bp <pos>` / `--to-bp <pos>`: Restrict scan range

### 4.2 Phenotype Files

**Standard format** (single/multi-trait):

```
FID    IID    Pheno1    Pheno2    Pheno3
IND001 IND001 1.23      2.45      NA
IND002 IND002 0.88      NA        1.67
```

- Missing values: `NA` or `-9`.
- `--mpheno <n>`: Phenotype column number (1-indexed), default 1.
- `--mpheno-list <list>`: Multi-trait column numbers (e.g., `1,2,3`).

**Long format** (Repeatability / RRM):

```
FID    IID    TIME    VALUE
IND001 IND001 1       1.23
IND001 IND001 2       1.35
```

- `--tpheno <n>`: TIME column number.

**Survival format** (Cox model):

```
FID    IID    TIME    STATUS
IND001 IND001 365     1
IND002 IND002 720     0
```

- `STATUS`: 1 = event occurred, 0 = censored.

**Maternal effects format**:

```
FID    IID    DAM_ID    PHENO
IND001 IND001 DAM001    1.35
```

- Column 3 is maternal ID (`--maternal-id <col>` to specify).

**Binomial event/trial format**:

```
FID    IID    EVENT    TRIAL
IND001 IND001 45       100
```

- Use `--binomial-pheno <file>`. Do not input proportions only.

### 4.3 Covariate Files

| Parameter | Type | Format |
|-----------|------|--------|
| `--covar <file>` | Categorical (sex, batch) | `FID IID CAT1 CAT2 ...` |
| `--qcovar <file>` | Quantitative (age, weight) | `FID IID VAL1 VAL2 ...` |

### 4.4 Input Adaptation Layer MVP

Converts `csv/tsv/txt` business tables to auditable prepared input for long-format models.

**Basic usage**:

```bash
./cosmicscanner \
  --bfile data/genotype \
  --input-file data/semen.csv \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration,motility \
  --record-qcovars age,temp \
  --subject-ccovars farm,line \
  --time-parser date \
  --time-index-multiplier auto \
  --prepare-input-only \
  --prepared-input-dir prep/semen
```

**Reuse prepared input**:

```bash
./cosmicscanner \
  --bfile data/genotype \
  --prepared-input-dir prep/semen \
  --input-trait volume \
  --repeatability-model \
  --gwas \
  --out results/volume_repeatability
```

**Using config file**:

```bash
# Generate template
./cosmicscanner --input-file data/semen.csv --write-template-config input_prepare.json

# Use config
./cosmicscanner --bfile data/genotype --input-config input_prepare.json --prepare-input-only
```

**Model integration**:

| `--input-model` | Files Used | Description |
|----------------|-----------|-------------|
| `repeatability` | `traits/<trait>.pheno.tsv` + record covariate | Record-level repeatability |
| `rrm` | `rrm/<trait>.pheno.tsv` | Trajectory time = `TIME_BASE` |
| `lmm` | `subject/<trait>.pheno.tsv` | Mean aggregation per subject |
| `multivariate` | `subject/multivariate.pheno.tsv` | Multi-trait mean aggregation |

> **Note**: Excel (`xlsx/xls`) is not supported in MVP. Save as `csv/tsv` first.

See [docs/INPUT_ADAPTATION_SPEC_CN.md](INPUT_ADAPTATION_SPEC_CN.md) for full spec.

---

## 5. Quick Start

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --covar data/my_covar.txt \
  --gwas \
  --out results/standard_gwas
```

Output: `results/standard_gwas.assoc.txt`

---

## 6. Detailed Usage & Examples

### 6.1 Single-Trait LMM GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --covar data/my_covar.txt \
  --qcovar data/my_qcovar.txt \
  --gwas \
  --scan-tier semi \
  --out results/standard_gwas
```

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--no-loco` | Disable LOCO, use global null | LOCO enabled |
| `--scan-tier <m>` | `approx`/`semi`/`exact`; `semi` is recommended genome-wide | Explicit |
| `--fast-lmm` | Legacy single-trait fixed-null alias; use `--scan-tier approx` in new commands | Off |
| `--lmm-null-solver <s>` | `brent`/`idul`/`idul_plus`/`aidul` | `brent` |
| `--maf <f>` | Minimum MAF | 0.01 |
| `--hwe <f>` | HWE p-value filter | 0.0 (off) |
| `--threads <n>` | OpenMP threads | Max available |

**Output**: `*.assoc.txt`

### 6.2 Multivariate LMM GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2,3 \
  --multivariate \
  --gwas \
  --scan-tier semi \
  --out results/multi_gwas
```

The unified multivariate tiers are fixed-null Wald (`approx`), fixed covariance ratio with per-SNP residual covariance (`semi`), and per-SNP full covariance refit (`exact`). Multivariate `exact` currently requires `--no-loco`. Do not combine `--scan-tier` and `--mv-method` in one command.

**Advanced method selection**:

| Method | Description | Use Case |
|--------|-------------|----------|
| `emmv` (default) | EMmvGWAS-style semi-exact Wald | Fast genome-wide scan |
| `wald` | Fixed-null basic Wald | Debugging |
| `gemma` | GEMMA-style per-SNP refit Wald | Small validation |
| `gemma-full` | GEMMA-style Wald/LRT/Score with `--mv-gemma-mode 1/2/3/4` | GEMMA mvLMM parity |
| `gemma-exact` | Full alternative refit with GEMMA EM/NR controls | Exact-style validation |
| `gemma-gxe` | Multivariate SNP-by-environment Wald | GEMMA-style multivariate GxE scan |
| `gemma-gxe-full` | Multivariate GxE Wald/LRT/Score | GEMMA-style multivariate GxE multi-test |
| `hybrid` (recommended) | EMmv fast scan + GEMMA-full refine | Large genome-wide |
| `eb` | Empirical Bayes random SNP | Validated auxiliary signal enhancement |

`gemma-exact` follows the GEMMA exact-style route with EM/NR control options, mode selection, and per-SNP full alternative refit. Add `--mv-crt` to enable GEMMA-style Edgeworth/PCRT small-sample correction; the output includes `P_WALD_CRT`, `P_LRT_CRT`, and `P_SCORE_CRT`, and the main `P` column follows the corrected p-value selected by `--mv-gemma-mode`. This is a native CosmicScanner/Eigen implementation rather than a direct link to GEMMA's GSL code.

**Hybrid method (recommended)**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method hybrid \
  --mv-test wald,lrt,score \
  --mv-refine-top 100 \
  --mv-refine-p 1e-3 \
  --mv-refine-score-prepass 500 \
  --mv-refine-discrepancy-top 50 \
  --mv-refine-discrepancy 1.0 \
  --gwas \
  --out results/mv_hybrid_loco
```

**Hybrid outputs**:
- `*.assoc.mv.txt`: Stage 1 genome-wide scan
- `*.refine.snplist`: Refined SNP list
- `*.refined.assoc.mv.txt`: Refined results
- `*.hybrid.summary.tsv`: Hybrid summary
- `*.score_prepass.refined.assoc.mv.txt`: Score prepass results (if enabled)

**GEMMA-full refit**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2,3 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --mv-refit-p 0.05 \
  --no-loco \
  --gwas \
  --out results/mv_gemma_full_refit
```

**EMmvGWAS R compatibility**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2 \
  --emmvgwas \
  --emmvgwas-output \
  --out results/emmvgwas_compatible
```

`--emmvgwas` is a direct alias for the EMmvGWAS-compatible multivariate route (`--multivariate --mv-method emmv --mv-compatible emmvgwas --no-loco`). Use `--emmvgwas-output` to write EMmvGWAS-style `G.csv`, `R.csv`, `beta.csv`, and association CSV files.

For the native EMmvGWAS two-step workflow, use `--emmvgwas-vc-only` for `vc_estimation`, then pass the generated matrices to `--emmvgwas --emmvgwas-G <G.csv> --emmvgwas-R <R.csv>` for `semi_gwas`.

**VCE solver control**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--mv-vce-solver <s>` | `em` or `aireml` | `em` |
| `--mv-em-min-iter <n>` | Min EM iterations | 10 |
| `--mv-ai-max-iter <n>` | Max AI-REML iterations | 50 |
| `--mv-ai-tol <f>` | AI-REML convergence | 1e-6 |

### 6.3 Repeatability Model GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat_pheno.txt \
  --tpheno 1 --mpheno 2 \
  --repeatability-model \
  --repeatability-mode exact \
  --gwas \
  --scan-tier semi \
  --repeatability-grm-cache grm_cache \
  --out results/repeat_gwas
```

> **Command convention**: multivariate and repeatability model flags currently start scans implicitly for compatibility, but canonical commands explicitly include `--gwas`, phenotype columns, and `--scan-tier`.

**Strategy** (`--repeatability-mode <m>`):

| Strategy | Description | Use Case |
|----------|-------------|----------|
| `auto` | Automatic (switches to STCG when records > 80000) | Default |
| `exact` | Record-level exact GLS + Woodbury | **Production recommended** |
| `fast` | Mean aggregation approximation | Quick screening |
| `stcg` | STCG matrix-free solver | Biobank-scale longitudinal |

> **Note**: `refit-validation` is not a valid value for `--repeatability-mode`. Use the separate flag `--repeatability-refit-validation` for per-SNP VCE refit.

**VCE mode** (`--vce-mode <m>`):

| Mode | Description |
|------|-------------|
| `auto` (default) | Scalable subject-level Woodbury VCE |
| `fast` | Fast VCE, skip exact REML refinement |
| `subject-level-exact` | Force subject-level exact REML refinement |
| `full-record-exact` | Full record-level exact VCE (small data validation only) |

**Key parameters**:

| Parameter | Description |
|-----------|-------------|
| `--repeatability-grm-cache <prefix>` | Disk cache for subject GRM via `mmap` |
| `--repeatability-mixed-precision` | Float32 to reduce memory |
| `--repeatability-weight <file>` | Record weights |
| `--repeatability-refit-validation` | Per-SNP VCE refit |
| `--repeatability-refit-limit <n>` | Refit limit (0 = all) |

**Outputs**: `*.assoc.txt`, `*.refit.assoc.txt`, `*.stcg.assoc.txt`

### 6.4 Random Regression Model (RRM) GWAS

**Balanced longitudinal**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --long-pheno data/long_pheno.txt \
  --mpheno 2 \
  --tpheno 1 \
  --rrm-order 2 \
  --gwas \
  --out results/rrm_gwas
```

**Unbalanced longitudinal**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/unbal_pheno.txt \
  --rrm-unbalance \
  --rrm-gwas-method record \   # or coefficient (default)
  --tpoint time_point \
  --tcovar age,bw \
  --tclass farm,parity \
  --rrm-order 2 \
  --rrm-basis legendre \
  --out results/rrm_unbal_gwas
```

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rrm-order <n>` | Basis order (1-4) | — |
| `--rrm-basis <m>` | `legendre` or `bspline` | `legendre` |
| `--rrm-gwas-method <m>` | `record` (record-level `x⊗Φ` joint test, Type B, `Experimental validated`) or `coefficient` (BLUP coefficient-space multivariate Wald, GMAT-aligned, default) | `coefficient` |
| `--rrm-curve-grid <n>` | Curve grid points | 21 |
| `--rrm-top <n>` | Top SNP curves | 100 |
| `--rrm-ai-iter <n>` | AI refinement iterations (-1=auto, 0=fast) | -1 |
| `--rrm-sig-threshold <p>` | Dynamic time window p | 0.05 |
| `--rrm-loco-global-null` | Fast approximate LOCO | Off |
| `--tpoint <col>` | Unbalanced time column | — |
| `--tcovar <cols>` | Time-dependent q-covariates | — |
| `--tclass <cols>` | Time-dependent c-covariates | — |
| `--rrm-compatible <m>` | `gma-fixreg` compatibility | — |
| `--rrm-compatible-varcom <file>` | GMA `.var` file | — |

**Outputs**: `*.rrm_joint.assoc.txt`, `*.rrm.top_curve.txt`, `*.rrm_unbal.assoc.txt`

The two unbalanced RRM scan paths differ statistically:

- `--rrm-gwas-method record` runs a record-level joint test with SNP design `x_{subject(r)} × φ_k(t_r)` (genotype × Legendre/B-spline basis at each record), a true time-varying SNP effect test (Type B), with explicit `STATUS/DF` diagnostics (`OK/RANK_DEFICIENT/ILL_CONDITIONED/NUMERICAL_FAILURE`). It is `Experimental validated`, supported by balanced-path equivalence (Pearson ≥ 0.99999) and scenarios A–G in `run_test_rrm_method_a_simulation.sh`.
- `--rrm-gwas-method coefficient` (the default) first BLUPs individual Legendre coefficients, then runs a Kronecker multivariate Wald scan over the coefficient traits; it aligns with GMAT (`run_gmat_unbalanced_alignment.sh`, autosome Pearson = 0.952, top SNP concordant) and is the historical-compatible path.

See `docs/RRM_GWAS_METHOD_AUDIT.md` for the full five-layer method audit, accuracy table, and the X-chromosome methodological note.

### 6.5 Maternal Effects Model GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/birth_weight.txt \
  --maternal-model \
  --maternal-id 3 \
  --covar data/covar.txt \
  --gwas \
  --out results/maternal_gwas
```

**Phenotype format**: `FID IID DAM_ID PHENO`

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--maternal-model` | Enable maternal effects model | — |
| `--maternal-id <col>` | Maternal ID column (1-indexed) | 3 |
| `--maternal-snp-test <mode>` | Model-specific SNP test: `direct` (offspring genotype, default), `maternal` (dam genotype single-column), `joint` (direct+maternal 2-df) | `direct` |

**Variance components**: Direct additive (Va), Maternal additive (Vm), Direct-maternal covariance (Cov_am), Maternal permanent environment (Vpe), Residual (Ve).

**Model-specific SNP tests (`--maternal-snp-test`, experimental)**:

```bash
# maternal: test the DAM genotype for a maternal SNP effect (H0: beta_m = 0)
./cosmicscanner --bfile data/my_geno --pheno data/birth_weight.txt \
  --maternal-model --maternal-id 3 --covar data/covar.txt --gwas \
  --maternal-snp-test maternal --out results/maternal_snp
# joint: 2-df Wald test H0: beta_d = beta_m = 0 at each SNP
./cosmicscanner --bfile data/my_geno --pheno data/birth_weight.txt \
  --maternal-model --maternal-id 3 --covar data/covar.txt --gwas \
  --maternal-snp-test joint --out results/joint_snp
```

- Eligible samples are animals with a phenotype **and** a genotyped dam (the dam must be a subject in the PLINK family file); non-genotyped dams are explicitly filtered and counted in `N_TOTAL / N_WITH_DAM_GENOTYPE / N_DROPPED_MISSING_DAM` in the output and in `*.maternal.snp.meta.txt`.
- The maternal scan uses the dam's ALT allele dosage (`x_dam`) centered per-SNP; the null covariance mirrors the maternal VCE: `σ²ad·K(offspring) + σ²am·K(dams) + σ²e·I`, factored once per chromosome and reused across SNPs. LOCO follows the same policy as the direct scan.
- Default `direct` reproduces the legacy `--maternal-model` behaviour byte-for-byte (verified regression).
- Outputs: `direct` → `*.assoc.txt`; `maternal` → `*.maternal.assoc.txt` (`CHR SNP BP A1 A2 FREQ N BETA_MATERNAL SE_MATERNAL P LOG10P` + N columns); `joint` → `*.joint.assoc.txt` (also `BETA_DIRECT SE_DIRECT BETA_MATERNAL SE_MATERNAL COV_BETA_DM CORR_BETA_DM STAT_JOINT DF P_JOINT COND RANK STATUS`; `status=COLLINEAR` when the joint information matrix is rank-deficient/ill-conditioned, e.g. offspring==dam).
- BGEN/PGEN inputs are not yet wired for the `maternal|joint` modes (PLINK only; an error is raised).

**Method note (statistical definition)**: the default (`direct`) model is a **maternal-adjusted GWAS**. The random-effect background is `Xb + Z_d a_d + Z_m a_m + e` with direct and maternal additive components (currently `G_m = G_d`; `Cov(a_d,a_m)` is not yet estimated). The SNP scan uses the **individual's own genotype** as a single direct fixed SNP effect (`x_j β_j`) tested with `H0: β_j = 0`, and outputs `BETA/SE/P` exactly like standard LMM. With `--maternal-snp-test maternal|joint`, the dam genotype enters the SNP design (maternal SNP effect `x_dam β_m` and 2-df joint test). Software status: `Experimental validated`, supported by a 10-scenario simulation suite including null calibration and direct/maternal/joint truth recovery; no external maternal-SNP reference is claimed. See [GWAS_MODEL_DEFINITIONS.md](GWAS_MODEL_DEFINITIONS.md#5-maternal-model--maternal-model).

**Output**: `*.assoc.txt`, and with `maternal`/`joint` also `*.maternal.assoc.txt` / `*.joint.assoc.txt` + `*.maternal.snp.meta.txt`.

### 6.6 Single-Step GWAS (ssGWAS)

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --single-step \
  --ped data/pedigree.txt \
  --genotyped-list data/genotyped_ids.txt \
  --gwas \
  --out results/ssgwas_out
```

This route uses the standard A22-corrected H⁻¹ block for pedigree-genomic joint association scanning.

**Method note (naming precision)**: two routes are provided.

- **Method 1 — H-matrix single-step GWAS** (default; `Experimental validated`). The scan tests the fixed SNP effect `H0: β_j = 0` with genotype `x_j` as the design and the H-blocks (built from `H⁻¹ = A⁻¹ + [0; τG⁻¹ − ωA22⁻¹]`) as the background covariance on the genotyped analysis sample, producing `BETA/SE/P`.
- **Method 2 — ssGBLUP back-solving GWAS** (`--single-step-method backsolve`; `Experimental validated`). It solves the single-step MME over the full pedigree (including non-genotyped animals), gives EBV/GEBV, then back-solves SNP effects `ĝ = Z'G⁻¹ û_g` (VanRaden 2008; Strandén-Garrick 2009; Wang 2012) and sliding-window genetic variance (`%Var_window`, 10-SNP windows, never crossing chromosomes). Validated analytically (max EBV diff 5e-7 vs brute-force reference) and in simulation (5/5 causal loci in top-10 SNPs; pedigree-shuffle invariance exact). No iterative weighted ssGBLUP is implemented.
- **Optional SNP P-values for Method 2**: `exact` is `Stable`; `apy-exact` is an `Experimental validated` scalable APY uncertainty path requiring a core; `apy-approx` is `Experimental approximate` (`STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`, raw P-values only, no genomic control/rescale, and scalability limited to the uncertainty layer). Do not report `apy-approx` as calibrated inference.

See [GWAS_MODEL_DEFINITIONS.md](GWAS_MODEL_DEFINITIONS.md#8-ssgwas--single-step).

Back-solving command:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --single-step --single-step-method backsolve \
  --ped data/pedigree.txt \
  --genotyped-list data/genotyped_ids.txt \
  --reml-vg 0.5 --reml-ve 1.0 \
  --out results/ssgwas_bs
```

Outputs: `*.ss_backsolve.ebv.txt` (all-animal EBV/GEBV), `*.ss_backsolve.snp.txt` (per-SNP effects: `EFFECT_ALLELE = (z'q/σ_z)/S` allele-substitution scale, `EFFECT_STD = EFFECT_ALLELE·σ_z` standardized scale with Z·STD = M·ALLELE = GEBV, `SNP_VAR = 2p(1−p)·EFFECT_ALLELE²`), `*.ss_backsolve.window.txt` (WINDOW_VARIANCE = Σ_window 2p(1−p)·EFFECT_ALLELE²), `*.ss_backsolve.meta.txt`.

**Pedigree format** (`--ped`, no header): `ID SIRE DAM`

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--single-step` | Enable ssGWAS | — |
| `--single-step-method lmm\|backsolve` | Method 1 (H-kernel LMM association) or Method 2 (breeding-value back-solving) | `lmm` |
| `--make-agh` | Build and write A-inverse, G, and H-inverse without phenotypes | — |
| `--write-agh` | Write A-inverse, G, and H-inverse during ssGWAS | Off |
| `--ped <file>` | Pedigree file | — |
| `--genotyped-list <file>` | Genotyped individual ID list | — |
| `--h2-tau <f>` | Tau scaling for A22 blend (`--single-step-method backsolve`); with `--h2-omega` also fixes the heritability ratio τ = h², ω = (1−h²)τ | 1.0 |
| `--h2-omega <f>` | Omega scaling | 1.0 |
| `--blend-g <f>` | G blend weight with I | 0.05 |
| `--ssgwas-window-snps <n>` | Sliding-window SNP count for `%Var_window` (backsolve) | 10 |
| `--ssgwas-solver ldlt\|pcg` | MME solver (backsolve) | auto |
| `--ssgwas-pvalue-method off\|exact\|apy-exact\|apy-approx` | Optional backsolve SNP P-value route; `apy-approx` is experimental approximate and not strict type-I controlled | `off` |
| `--ssgwas-apy-core-size <n>` | Required APY core size for `apy-exact` / `apy-approx` | — |

Matrix-only command:

```bash
./cosmicscanner --bfile data/my_geno --make-agh \
  --ped data/pedigree.txt --genotyped-list data/genotyped_ids.txt \
  --out results/agh
```

Matrix outputs are `*.Ainv.bin/.id`, `*.G.bin/.id`, `*.Hinv.bin/.id`, and `*.agh.meta.txt`. ssGWAS outputs remain `*.assoc.txt` and `*_ssgwas_null.txt`. Phenotyped animals in VCE/scanning must currently be genotyped; additional pedigree animals may contribute ancestry.

### 6.7 Threshold/Probit Model GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/binary.pheno \
  --mpheno 1 \
  --gwas \
  --threshold-model \
  --link probit \
  --out results/psr_binary
```

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--threshold-model` | Enable binary PSR/probit | — |
| `--link probit\|logit` | Link function | `probit` |
| `--spa` | Saddle Point Approximation | Off |

> **Tip**: Add `--spa` for improved calibration on binary traits.

### 6.8 Count/Binomial PSR-GWAS

**Poisson count**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/count.pheno \
  --mpheno 1 \
  --gwas \
  --count-model \
  --family poisson \
  --offset data/count.log_offset \
  --out results/psr_poisson
```

**Binomial proportion**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --count-model \
  --family binomial \
  --binomial-pheno data/binomial_event_trial.pheno \
  --out results/psr_binomial
```

> **Note**: Use event/trial columns, not proportions only.

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--count-model` | Enable count/binomial PSR | — |
| `--family poisson\|binomial\|negbin` | Distribution family | `poisson` |
| `--offset <file>` | Poisson log-offset | — |
| `--binomial-pheno <file>` | Binomial event/trial file | — |

> **Note**: `--family negbin` is experimental; overdispersion not fully validated.

### 6.9 Ordinal PSR-GWAS

**Recommended: stacked proportional-odds PSR**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --ordinal-model \
  --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --out results/psr_ordinal
```

**Baseline: cumulative-binary**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --ordinal-model \
  --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode experimental-cumulative-binary \
  --out results/psr_ordinal_baseline
```

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--ordinal-model` | Enable ordinal PSR | — |
| `--ordinal-pheno <file>` | Ordinal phenotype `FID IID ORDINAL` | — |
| `--ordinal-mode <name>` | `experimental-cumulative-binary` / `stacked-psr` | `experimental-cumulative-binary` |
| `--ordinal-stacked-engine <name>` | `fixed-gls` / `dense-zkz` | `fixed-gls` |
| `--ordinal-vc-method <name>` | `fixed` / `grid-reml` | `fixed` |
| `--ordinal-test <name>` | `wald`/`score`/`robust-score`/`all` | `wald` |
| `--ordinal-allow-fallback` | Fall back to cumulative-binary | Off |
| `--ordinal-sigma2-g <f>` | Fixed genetic variance | — |
| `--ordinal-vc-grid-min <f>` | Grid search lower bound | — |
| `--ordinal-vc-grid-max <f>` | Grid search upper bound | — |
| `--ordinal-vc-grid-points <n>` | Grid search points | — |

**Outputs**: `*.gwas.tsv`, `*.ordinal.score.tsv` (with `--ordinal-test all`)

### 6.10 Survival Cox Model GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/survival.pheno \
  --survival-model \
  --covar data/covar.txt \
  --gwas \
  --out results/survival_gwas
```

**Phenotype format**: `FID IID TIME STATUS` (STATUS: 1=event, 0=censored)

> **Implementation note**: The current survival model is a simplified implementation — it uses Nelson-Aalen cumulative hazard to estimate martingale residuals, then runs a fast-LMM scan on the residuals, rather than a full Cox PH likelihood solve. Suitable for screening; significant loci should be validated with R `coxme`/`GMMAT`.

**Output**: `*.assoc.txt`

### 6.11 GxE Reaction Norm Model GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/gxe_pheno.txt \
  --gxe-model \
  --env-covar 3 \
  --covar data/covar.txt \
  --gwas \
  --out results/gxe_gwas
```

**Phenotype format**: `FID IID ENV PHENO`

- `--env-covar <col>`: Environment variable column (1-indexed).

**Output**: `*.gxe.assoc.txt`

**Columns**: `CHR SNP BP A1 A2 FREQ N BETA_G SE_G BETA_GxE SE_GxE JOINT_CHISQ JOINT_P`

### 6.12 Rare Variant Analysis (SKAT / Burden)

**SKAT**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --set-file data/gene_sets.txt \
  --skat \
  --gwas \
  --out results/skat_gwas
```

**Burden**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --set-file data/gene_sets.txt \
  --burden \
  --gwas \
  --out results/burden_gwas
```

> **Limitation**: SKAT/Burden are only available in the default single-trait LMM path and cannot be combined with other model flags (`--multivariate`, `--repeatability-model`, etc.).

**Set file format**: `GENE1 rs123 rs456 rs789`

**SPA for binary traits**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/binary.pheno \
  --gwas \
  --spa \
  --out results/spa_gwas
```

### 6.13 Biobank-Scale Analysis

**Full 2-step**:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --block-size 1000 \
  --gwas \
  --out results/biobank_gwas
```

**Step-by-step**:

```bash
# Step 1: Generate LOCO predictions
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --regenie-step 1 \
  --write-pred results/loco_pred.txt \
  --gwas \
  --out results/biobank_step1

# Step 2: Association using predictions
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --regenie-step 2 \
  --use-pred results/loco_pred.txt \
  --gwas \
  --out results/biobank_step2
```

> **Limitation**: Biobank mode is only available in the default single-trait LMM path and cannot be combined with other model flags.

**Key parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--biobank-mode` | Enable REGENIE-style 2-step | — |
| `--block-size <n>` | Ridge block size | 1000 |
| `--regenie-model-snps <file>` | Step 1 model SNP list | — |
| `--write-pred <file>` | Output Step 1 LOCO predictions | — |
| `--use-pred <file>` | Use existing predictions | — |
| `--regenie-step <n>` | 0=both, 1=predictions, 2=association | 0 |

**Output**: `*.biobank.assoc.txt`

### 6.14 Empirical Bayes EB-GWAS

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --model lmm-eb \
  --test eb-wald \
  --adjust eb-mbonf \
  --gwas \
  --out results/eb_gwas
```

**Key parameters**:

| Parameter | Description |
|-----------|-------------|
| `--model lmm-eb` | Enable EB-LMM (alias: `--gwas-method eb-lmm`) |
| `--test eb-wald` | Enable EB Wald test |
| `--adjust eb-mbonf` | Experimental threshold in `.eb.summary.txt` |

> **Note**: EB-GWAS is a validated auxiliary route. Always report standard LMM results and ordinary Bonferroni/FDR alongside.

### 6.15 Common-Environment Model GWAS (Experimental)

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/long_pheno.txt \
  --tpheno 2 \
  --common-env-model \
  --common-env-file data/group.txt \
  --gwas \
  --out results/ce_gwas
```

**Model**: `y = Xb + x_j * beta_j + Z_a * a + Z_c * c + e`, with
`a ~ N(0, K * sigma2_a)` (genomic GRM across subjects) and
`c ~ N(0, I * sigma2_c)` (i.i.d. common-environment effects across groups), so
records sharing a group share a `sigma2_c` correlation. Suitable when records
of one subject come from different environments and a random intercept per
environment (not per subject) is desired.

**Phenotype format**: long format with `FID IID [TIME] VALUE` columns
(`--tpheno <col>` selects the value column, 1-indexed).

**Group file** (`--common-env-file <f>`, whitespace-separated):

- Subject-level (default): one line per subject, `FID IID GROUP`; every record
  of that subject is assigned to the group.
- Record-level (`--common-env-record-level`): one line per record,
  `FID IID GROUP` (optionally `FID IID TIME GROUP`); the first
  columns must match the phenotype file record order.

| Parameter | Description |
|-----------|-------------|
| `--common-env-model` | Enable the common-environment model |
| `--common-env-file <f>` | Group assignment file (subject- or record-level) |
| `--common-env-record-level` | Interpret the group file as record-level |
| `--reml-vg/-vc/-ve <v>` | Fix variance components, skipping VCE |
| `--no-loco` | Disable leave-one-chromosome-out scans |

**Output**:
- `*.assoc.txt`: association results (columns `CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`).
- `*.common_env.null.txt`: null-model summary (`sigma2_a / sigma2_common / sigma2_e / h2 / common_environment_fraction`).

> **Note**: experimental model in the development tree; not part of the
> focused CosmicGWAS v1.1.0 release scope. PLINK/BGEN/PGEN inputs are
> supported; MAF/HWE filters and LOCO follow the same options as the
> single-trait LMM path.

---

### 6.16 Pedigree Animal Model GWAS

The simple animal model replaces the genomic GRM with a pedigree-based
additive relationship matrix `A`:

```
y = Xb + x_j * beta_j + Z*a + e,   a ~ N(0, A * sigma2_a)
```

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --animal-model \
  --ped data/pedigree.txt \
  --gwas \
  --out results/animal_assoc
```

**Pedigree file** (`--ped <f>`, whitespace-separated, no header):

| Column | Content |
|--------|---------|
| 1 | Individual ID (must include every analysis sample) |
| 2 | Sire ID (`NA`, `.`, or `0` if unknown) |
| 3 | Dam ID (`NA`, `.`, or `0` if unknown) |

IDs are matched against samples as `FID:IID` first, then bare `IID`
(falls back when the `.fam` FID is a placeholder). The pedigree may
include extra untyped ancestors; only the analysis subset is used to
build the `A22` block.

**Options**:

| Parameter | Description |
|-----------|-------------|
| `--animal-model` | Enable the pedigree animal model |
| `--ped <f>` | Pedigree file `ID SIRE DAM` (required) |
| `--write-ped-kernel` | Dump the exact analysis-subset `A22` matrix to `*.animal_model.k.txt` + `*.animal_model.k.id` (pass directly to GEMMA `-k`) |
| `--no-loco` | Disable LOCO (needed to silence the LOCO-not-applicable warning on multi-chromosome data) |

**Output**:
- `*.assoc.txt`: association results (columns `CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`).
- `*.animal_model.null.txt`: null-model summary (`sigma2_a / sigma2_e / h2`, plus pedigree parse/`Ainv`/`A22`/VCE/scan timings).
- `*.animal_model.k.txt` / `*.animal_model.k.id`: exact dense `A22` kernel and sample order (only with `--write-ped-kernel`).

**Notes**:
- The animal model shares the same `PedigreeRelationship` parser and
  `Ainv`/`A22` construction with ssGWAS (`buildA22` solves the Schur
  complement of `Ainv`; no dense `A` is materialized). With
  `--h2-tau 0 --h2-omega 0 --no-loco`, ssGWAS reduces to this model and
  produces per-SNP identical BETA/SE/-log10P.
- LOCO is not applicable (only one relationship matrix and no chromosome
  separation); the tool prints a warning on multi-chromosome data unless
  `--no-loco` is given.
- For a GEMMA-style external check, GEMMA 0.98.5 reads the dense
  `*.animal_model.k.txt` via `-k` and `-lmm 4`; no GRM normalization is
  required (scale is absorbed into `sigma2_g`).
- Validation: `tests/run_test_animal_model.sh` (simulation suite) and
  `tests/run_animal_model_gemma_validation.sh` (external GEMMA `-k`
  benchmark; SKIP if no GEMMA binary).

---

## 7. Output File Formats

### 7.1 Standard LMM GWAS — `*.assoc.txt`

| Column | Description |
|--------|-------------|
| CHR | Chromosome |
| SNP | SNP ID |
| BP | Physical position (bp) |
| A1 | Effect allele |
| A2 | Non-effect allele |
| FREQ | A1 allele frequency |
| N | Sample size |
| BETA | Effect estimate |
| SE | Standard error |
| P | p-value |
| LOG10P | -log10(P) |

**EB-LMM extra columns**: `EB_BETA`, `EB_SE`, `EB_WALD`

### 7.2 Multivariate GWAS — `*.assoc.mv.txt`

| Column | Description |
|--------|-------------|
| CHR / SNP / BP / A1 / A2 | Basic info |
| FREQ | A1 frequency |
| WALD_STAT | Wald statistic |
| P | p-value |

GEMMA-full mode adds LRT, Score columns.

### 7.3 RRM GWAS — `*.rrm_joint.assoc.txt`

| Column | Description |
|--------|-------------|
| CHR / SNP / BP / A1 / A2 | Basic info |
| FREQ | A1 frequency |
| N | Sample size |
| BETA_L0 | Baseline effect |
| SE_L0 | Baseline SE |
| JOINT_CHISQ | Joint test chi-square |
| JOINT_P | Joint test p-value |

### 7.4 GxE GWAS — `*.gxe.assoc.txt`

| Column | Description |
|--------|-------------|
| CHR / SNP / BP / A1 / A2 / FREQ / N | Basic info |
| BETA_G | Main effect (genotype) |
| SE_G | Main effect SE |
| BETA_GxE | Interaction effect |
| SE_GxE | Interaction SE |
| JOINT_CHISQ | 2-df joint Wald chi-square |
| JOINT_P | Joint p-value |

### 7.5 Ordinal PSR-GWAS — `*.gwas.tsv`

Default: corrected Wald / pseudo-score. `--ordinal-test all` adds `*.ordinal.score.tsv`.

### 7.6 Biobank GWAS — `*.biobank.assoc.txt`

Same format as standard `*.assoc.txt`.

### 7.7 EBI GWAS-SSF Format

Add `--out-ssf` to output EBI GWAS-SSF standard format.

---

## 8. Complete Parameter Reference

### 8.1 Input/Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--bfile <prefix>` | PLINK 1.9 binary prefix | — |
| `--bgen <file>` | BGEN genotype file | — |
| `--sample <file>` | BGEN sample file (Oxford) | — |
| `--pfile <prefix>` | PLINK 2.0 binary prefix | — |
| `--pgen <file>` / `--pvar <file>` / `--psam <file>` | PLINK 2.0 individual files | — |
| `--out <prefix>` | Output prefix | `cosmic_gwas` |
| `--out-ssf` | Output EBI GWAS-SSF format | Off |
| `--extract <file>` | Include only listed SNPs | — |
| `--exclude <file>` | Exclude listed SNPs | — |

### 8.2 Phenotype/Covariates

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pheno <file>` | Phenotype file | — |
| `--mpheno <n>` | Phenotype column (1-indexed) | 1 |
| `--mpheno-list <list>` | Multi-trait columns (e.g., `1,2,3`) | — |
| `--tpheno <n>` | TIME column (Repeatability/RRM) | — |
| `--covar <file>` | Categorical covariates | — |
| `--qcovar <file>` | Quantitative covariates | — |
| `--long-pheno <file>` | Long-format phenotype (enables RRM) | — |
| `--binomial-pheno <file>` | Binomial event/trial file | — |
| `--offset <file>` | Poisson log-offset | — |

### 8.3 Variant Filtering

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--maf <f>` | Minimum MAF | 0.01 |
| `--hwe <f>` | HWE p-value filter | 0.0 (off) |
| `--chr <chr>` | Restrict to chromosome | — |
| `--from-bp <pos>` | Scan start position | — |
| `--to-bp <pos>` | Scan end position | — |

### 8.4 LMM Engine

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--gwas` | Enable GWAS scan (**required for all models**) | Off |
| `--no-loco` | Disable LOCO | LOCO enabled |
| `--fast-lmm` | Fixed-delta approximate acceleration | Off |
| `--lmm-null-solver <s>` | `brent`/`idul`/`idul_plus`/`aidul` | `brent` |
| `--no-covar-proj` | Disable covariate projection | Enabled |
| `--force-slow` | Force slow path | Off |
| `--spa` | Saddle Point Approximation (binary) | Off |

> **Important**: `--gwas` is the mandatory switch to trigger scanning. Specifying only model flags (`--repeatability-model`, `--multivariate`, `--maternal-model`, etc.) without `--gwas` will exit normally but produce no association output. Exception: `--rrm-unbalance` implicitly enables `--gwas`.

### 8.5 REML / VCE Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reml-maxit <n>` | Max VCE iterations | 50 |
| `--reml-tol <f>` | VCE convergence | 1e-6 |
| `--reml-vg <f>` | Fixed genetic variance | — |
| `--reml-vpe <f>` | Fixed permanent environment variance | — |
| `--reml-ve <f>` | Fixed residual variance | — |

### 8.6 Multivariate GWAS

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--multivariate` | Enable multivariate GWAS | — |
| `--mv-method <m>` | `emmv`/`wald`/`eb`/`gemma`/`gemma-full`/`hybrid` | `emmv` |
| `--mv-compatible <name>` | Compatibility mode (`emmvgwas`) | — |
| `--mv-test <list>` | GEMMA-full tests: `wald,lrt,score` | All |
| `--mv-refine-top <n>` | Refine SNP limit | 100 |
| `--mv-refine-p <p>` | Refine p threshold | 1e-3 |
| `--mv-refine-score-prepass <n>` | Score prepass top n | — |
| `--mv-refine-discrepancy-top <n>` | Discrepancy SNP refine count | — |
| `--mv-refine-discrepancy <x>` | Min \|log10P diff\| | 1.0 |
| `--mv-refit-p <p>` | GEMMA-full refit threshold | 1.0 |
| `--mv-vce-solver <s>` | `em` or `aireml` | `em` |
| `--mv-em-min-iter <n>` | Min EM iterations | 10 |
| `--mv-ai-max-iter <n>` | Max AI-REML iterations | 50 |
| `--mv-ai-tol <f>` | AI-REML convergence | 1e-6 |

### 8.7 Repeatability Model

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--repeatability-model` | Enable repeatability model | — |
| `--repeatability-mode <m>` | `auto`/`exact`/`fast`/`stcg` | `auto` |
| `--repeatability-grm-cache <prefix>` | GRM disk cache | — |
| `--repeatability-mixed-precision` | Float32 memory reduction | Off |
| `--repeatability-weight <file>` | Record weights | — |
| `--repeatability-refit-validation` | Per-SNP VCE refit | Off |
| `--repeatability-refit-limit <n>` | Refit limit (0=all) | 0 |
| `--vce-mode <m>` | `auto`/`fast`/`subject-level-exact`/`full-record-exact` | `auto` |

### 8.8 Random Regression Model (RRM)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rrm-order <n>` | Basis order (1-4) | — |
| `--rrm-basis <m>` | `legendre` or `bspline` | `legendre` |
| `--rrm-curve-grid <n>` | Curve grid points | 21 |
| `--rrm-top <n>` | Top SNP curves | 100 |
| `--rrm-ai-iter <n>` | AI iterations (-1=auto, 0=fast) | -1 |
| `--rrm-sig-threshold <p>` | Dynamic time window p | 0.05 |
| `--rrm-unbalance` | Unbalanced longitudinal mode | Off |
| `--rrm-loco-global-null` | Global null fast LOCO | Off |
| `--tpoint <col>` | Unbalanced time column | — |
| `--tcovar <cols>` | Time-dependent q-covariates | — |
| `--tclass <cols>` | Time-dependent c-covariates | — |
| `--rrm-compatible <m>` | `gma-fixreg` compatibility | — |
| `--rrm-compatible-varcom <file>` | GMA `.var` file | — |

### 8.9 Maternal Effects Model

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--maternal-model` | Enable maternal effects model | — |
| `--maternal-id <col>` | Maternal ID column (1-indexed) | 3 |

### 8.10 Single-Step GWAS (ssGWAS)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--single-step` | Enable ssGWAS | — |
| `--ped <file>` | Pedigree file `ID SIRE DAM` | — |
| `--genotyped-list <file>` | Genotyped ID list | — |
| `--write-ebv` | Reserved switch (legacy); `--single-step-method backsolve` always writes `*.ss_backsolve.ebv.txt` | Off |
| `--h2-tau <f>` | Tau scaling | 1.0 |
| `--h2-omega <f>` | Omega scaling | 1.0 |
| `--blend-g <f>` | G blend weight with I | 0.05 |

### 8.11 Threshold/Probit Model

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--threshold-model` | Enable binary PSR/probit | — |
| `--link probit\|logit` | Link function | `probit` |

### 8.12 Count/Binomial PSR-GWAS

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--count-model` | Enable count/binomial PSR | — |
| `--family poisson\|binomial\|negbin` | Distribution family | `poisson` |
| `--offset <file>` | Poisson log-offset | — |
| `--binomial-pheno <file>` | Binomial event/trial | — |

### 8.13 Ordinal PSR-GWAS

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--ordinal-model` | Enable ordinal PSR | — |
| `--ordinal-pheno <file>` | Ordinal phenotype file | — |
| `--ordinal-mode <name>` | `experimental-cumulative-binary` / `stacked-psr` | `experimental-cumulative-binary` |
| `--ordinal-stacked-engine <name>` | `fixed-gls` / `dense-zkz` | `fixed-gls` |
| `--ordinal-vc-method <name>` | `fixed` / `grid-reml` | `fixed` |
| `--ordinal-test <name>` | `wald`/`score`/`robust-score`/`all` | `wald` |
| `--ordinal-allow-fallback` | Fall back to cumulative-binary | Off |
| `--ordinal-sigma2-g <f>` | Fixed genetic variance | — |
| `--ordinal-vc-grid-min <f>` | Grid search lower bound | — |
| `--ordinal-vc-grid-max <f>` | Grid search upper bound | — |
| `--ordinal-vc-grid-points <n>` | Grid search points | — |

### 8.14 PSR Iteration Control

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--psr-maxit <n>` | PSR/PQL max iterations | 50 |
| `--psr-tol <f>` | Convergence tolerance | 1e-6 |
| `--psr-mu-eps <f>` | Mu clipping epsilon | 1e-6 |
| `--psr-min-weight <f>` | Min working weight | 1e-8 |
| `--psr-max-weight <f>` | Max working weight | 1e8 |

### 8.15 Survival Cox Model

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--survival-model` | Enable Cox model | — |

### 8.16 GxE Reaction Norm Model

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--gxe-model` | Enable GxE model | — |
| `--env-covar <col>` | Environment variable column (1-indexed) | — |

### 8.17 Rare Variant

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skat` | SKAT test | — |
| `--burden` | Burden test | — |
| `--set-file <file>` | Variant set definition | — |

### 8.18 Biobank Mode

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--biobank-mode` | REGENIE-style 2-step | — |
| `--block-size <n>` | Ridge block size | 1000 |
| `--regenie-model-snps <file>` | Step 1 model SNPs | — |
| `--write-pred <file>` | Output LOCO predictions | — |
| `--use-pred <file>` | Use existing predictions | — |
| `--regenie-step <n>` | 0=both, 1=predictions, 2=association | 0 |

### 8.19 GRM Construction

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--make-grm` | Calculate GRM and exit | — |
| `--make-grm-alg <n>` | 1=VanRaden, 2=Zeng, 3=Yang, 4=Vitezica, 5=TrueVar, 6=EMmvGWAS-compatible | 1 |
| `--make-grm-sparse` | Sparse GRM | Off |
| `--sparse-cutoff <f>` | Sparse cutoff | 0.05 |
| `--snp-weight <file>` | SNP weights | — |
| `--pop-class <file>` | Population stratification | — |
| `--force-slow` | Force float slow path (verification) | Off |

> **Note**: `--make-grm-alg 6` is only enabled with `--mv-compatible emmvgwas` (auto-set by the multivariate model); using `--make-grm --make-grm-alg 6` alone does not trigger EMmvGWAS encoding.

### 8.20 EB-GWAS

| Parameter | Description |
|-----------|-------------|
| `--model lmm-eb` | Enable EB-LMM (alias: `--gwas-method eb-lmm`) |
| `--test eb-wald` | EB Wald test |
| `--adjust eb-mbonf` | Experimental threshold |

### 8.21 Input Adaptation Layer

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input-file <file>` | Raw csv/tsv/txt table | — |
| `--input-config <json>` | Load from JSON config | — |
| `--write-template-config <json>` | Generate template and exit | — |
| `--input-format <fmt>` | `auto`/`csv`/`tsv`/`txt` | `auto` |
| `--input-model <model>` | `repeatability`/`rrm`/`lmm`/`multivariate` | `repeatability` |
| `--prepared-input-dir <dir>` | Prepared input directory | — |
| `--prepare-input-only` | Prepare input and exit | — |
| `--id-col <name>` | ID column name | — |
| `--iid-col <name>` | IID column (when id-col is FID) | — |
| `--time-col <name>` | Time column name | — |
| `--trait-cols <list>` | Trait columns (comma-sep) | — |
| `--input-trait <name>` | Trait to run after preparation | First |
| `--record-qcovars <list>` | Record-level q-covariates | — |
| `--record-ccovars <list>` | Record-level c-covariates | — |
| `--subject-qcovars <list>` | Subject-level q-covariates | — |
| `--subject-ccovars <list>` | Subject-level c-covariates | — |
| `--time-parser <mode>` | `auto`/`date`/`integer`/`batch` | `auto` |
| `--time-unique-policy <p>` | `strict`/`within-time-index` | `within-time-index` |
| `--time-index-multiplier <n\|auto>` | TIME encoding multiplier | `auto` |
| `--preview-schema` | Output schema preview and exit | — |
| `--infer-columns` | Column inference suggestions | — |
| `--wide-to-long` | Wide-to-long conversion | — |

### 8.22 Engineering

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--threads <n>` | OpenMP threads | Max available |
| `--omp <n>` | Alias for `--threads` | — |

---

## 9. Performance Tuning

### 9.1 Accelerate Null Model Solving

```bash
# IDUL (medium-scale data)
--lmm-null-solver idul

# IDUL+ or AIDUL (large-scale data)
--lmm-null-solver aidul
```

See [docs/IDUL_AIDUL_ROADMAP.md](IDUL_AIDUL_ROADMAP.md).

### 9.2 Repeatability Model Optimization

```bash
# GRM disk cache reuse
--repeatability-grm-cache grm_cache

# Float32 memory reduction
--repeatability-mixed-precision
```

### 9.3 Multi-threading

```bash
--threads 16
```

### 9.4 Biobank Large-Scale

```bash
--biobank-mode --block-size 1000
```

### 9.5 Fast Screening

```bash
# Fast LMM approximation
--fast-lmm

# Repeatability fast mode
--repeatability-mode fast
```

---

## 10. FAQ

### Q1: What do PSR-GWAS convergence warnings mean?

Convergence warnings are diagnostics, **not failures**:

| `status_class` | Meaning |
|----------------|---------|
| `CONVERGED` | Normal convergence |
| `CONVERGED_WITH_WARNING` | Converged but with boundary variance, non-monotone quasi-likelihood, or clipping |
| `NOT_CONVERGED` | Not converged; use cautiously in formal analysis |

**Warning flags**:

| Flag | Meaning |
|------|---------|
| `WARNING_PSR_MAXIT_REACHED` | Iteration limit reached |
| `WARNING_PSR_CLIPPED` | mu or working weight hit clipping bounds |
| `WARNING_PSR_SIGMA_BOUNDARY` | Genetic variance near zero |
| `WARNING_PSR_QUASILIK_DECREASED` | Quasi-likelihood proxy decreased |

### Q2: How to handle difficult-to-converge or extremely sparse traits?

Adjust PSR/PQL parameters:

```bash
--psr-maxit 100 \
--psr-tol 1e-8 \
--psr-mu-eps 1e-8 \
--psr-min-weight 1e-10 \
--psr-max-weight 1e10
```

### Q3: Should PSR-GWAS significant loci be validated?

**Yes**. PSR-GWAS is not exact GLMM. Validate significant loci with candidate refinement or external exact GLMM (GMMAT/POLMM).

### Q4: How to choose multivariate method?

- **Large genome-wide**: `--mv-method hybrid` (recommended)
- **Small validation**: `--mv-method gemma-full`
- **Fast scan**: `--mv-method emmv` (default)

### Q5: How to benchmark against EMmvGWAS R?

Use `--mv-compatible emmvgwas --no-loco` to enable EMmvGWAS-compatible GRM (`--make-grm-alg 6`).

### Q6: Why can't I input binomial proportions only?

Proportion-only input loses trial size information, causing incorrect variance estimation. Use event/trial columns.

### Q7: How to output EBI GWAS-SSF format?

Add `--out-ssf`.

### Q8: Does the input adaptation layer support Excel?

Not in MVP. Save as `csv/tsv` first.

### Q9: How to calculate GRM only without GWAS?

```bash
./cosmicscanner --bfile data/my_geno --make-grm --out results/grm
```

### Q10: Why did I specify a model flag but got no association output?

**You must explicitly add `--gwas` to trigger scanning.** For example, the following command produces no output:

```bash
# Wrong: missing --gwas
./cosmicscanner --bfile geno --pheno pheno.txt --repeatability-model --out out
```

Correct:

```bash
./cosmicscanner --bfile geno --pheno pheno.txt --repeatability-model --gwas --out out
```

Exception: `--rrm-unbalance` implicitly enables `--gwas`.

### Q11: How to restrict scan to one chromosome?

```bash
--chr 1
```

Or restrict range:

```bash
--chr 1 --from-bp 1000000 --to-bp 2000000
```

---

## 11. Related Documentation

| Document | Description |
|----------|-------------|
| [Tutorial_CN.md](Tutorial_CN.md) | Chinese tutorial |
| [PSR_GWAS_ROADMAP.md](PSR_GWAS_ROADMAP.md) | PSR-GWAS roadmap |
| [PSR_GWAS_VALIDATION_REPORT_V0.5.md](PSR_GWAS_VALIDATION_REPORT_V0.5.md) | PSR-GWAS v0.5 validation report |
| [PSR_GWAS_FULL_ORDINAL_DESIGN.md](PSR_GWAS_FULL_ORDINAL_DESIGN.md) | Full ordinal PSR design |
| [PSR_GWAS_BENCHMARK_GUIDE.md](PSR_GWAS_BENCHMARK_GUIDE.md) | PSR-GWAS benchmark guide |
| [tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md](tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md) | Ordinal PSR tutorial |
| [REPEATABILITY_GWAS_RELEASE_CN.md](REPEATABILITY_GWAS_RELEASE_CN.md) | Repeatability GWAS release notes |
| [REPEATABILITY_GWAS_DEVELOPMENT_SUMMARY_20260529_CN.md](REPEATABILITY_GWAS_DEVELOPMENT_SUMMARY_20260529_CN.md) | Repeatability GWAS development summary |
| [REPEATABILITY_GWAS_MANUSCRIPT_CN.md](REPEATABILITY_GWAS_MANUSCRIPT_CN.md) | Repeatability GWAS manuscript |
| [RRM_GWAS_VALIDATION_CN.md](RRM_GWAS_VALIDATION_CN.md) | RRM methods & validation |
| [RRM_GWAS_INPUT_DESIGN_CN.md](RRM_GWAS_INPUT_DESIGN_CN.md) | RRM input reuse design |
| [IDUL_AIDUL_ROADMAP.md](IDUL_AIDUL_ROADMAP.md) | IDUL/AIDUL acceleration roadmap |
| [EB_GWAS_MODULE_CN.md](EB_GWAS_MODULE_CN.md) | EB-GWAS module |
| [INPUT_ADAPTATION_SPEC_CN.md](INPUT_ADAPTATION_SPEC_CN.md) | Input adaptation spec |
| [models/NEW_MODELS_DESIGN.md](models/NEW_MODELS_DESIGN.md) | New models design |
| [GWAS_RELEASE_CLOSEOUT_CN.md](GWAS_RELEASE_CLOSEOUT_CN.md) | GWAS release closeout checklist |
