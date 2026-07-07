# GWAS IDUL/AIDUL Acceleration Roadmap

## Scope

CosmicScanner keeps `brent` as the default exact null-model solver and exposes faster experimental solvers via:

```bash
--lmm-null-solver brent|idul|idul_plus|aidul
```

These solvers estimate the single LMM variance ratio `eta = sigma_g2 / sigma_e2` after GRM eigendecomposition and data rotation. They accelerate the null-model variance-ratio update, not GRM construction, eigendecomposition, genotype decoding, or output I/O.

## Applicability Matrix

| Model/path | Status | Rationale |
|---|---|---|
| Standard single-trait LMM-GWAS | Implemented | Model is one GRM plus residual; IDUL/AIDUL directly applies. |
| Standard LOCO LMM-GWAS | Implemented | Each chromosome uses a one-GRM LOCO null model; solver option is reused per chromosome. |
| Repeatability fast/aggregated scan | Implemented | Aggregated subject-level fast path reduces to one-GRM Eigen LMM; exact record-level repeatability remains separate. |
| Multi-phenotype independent scans | Planned | Reuse one GRM eigendecomposition and rotated genotype across many phenotypes; this is the main IDUL/OmiGA speed target. |
| SNP-specific variance refit | Planned/experimental only | IDUL/AIDUL can replace per-SNP Newton/REML if a refit mode is explicitly requested; default GWAS should keep null-model reuse. |
| molQTL / high-throughput omics | Planned | Needs a dedicated batch phenotype engine: shared rotated genotype/covariates, phenotype blocks, null-cache per phenotype, and ACAT/Clipper-style multiple testing. |
| Multivariate LMM | Not directly applicable | Current model estimates trait covariance matrices, not a single scalar eta. Needs a separate matrix-variate update. |
| Exact record-level repeatability/RRM | Not directly applicable | Models include genetic, permanent environmental/time-basis, and residual components; require record-level mixed-model updates. |

## Four Implementation Tracks

1. LOCO reuse: already connected. Use `RUN_LOCO=1 LMM_NULL_SOLVER=aidul bash CosmicScanner/tests/run_gemma_speed_benchmark.sh` to benchmark.
2. Multi-phenotype independent GWAS: implement a batch univariate engine that prepares GRM/eigenvectors once, rotates genotype blocks once, and loops over phenotype blocks with IDUL/AIDUL null caches.
3. SNP-specific refit: add an explicit opt-in validation mode, never default, where selected SNPs update eta with IDUL_PLUS/AIDUL and compare against Brent/GEMMA.
4. molQTL: build on the batch univariate engine, add cis-window scanning and trait-level multiple-testing support.

## Biobank-Scale GWAS Route

The UKB-scale path must not use dense exact LMM as its main engine. The primary route is REGENIE-style two-step regression:

1. Step 1: block-wise ridge predictors from common model SNPs, with LOCO predictions.
2. Step 2: streaming single-SNP association tests after residualizing phenotype and genotype against covariates plus the chromosome-specific LOCO prediction.

Current implementation status:

- `--biobank-mode` supports PLINK `--bfile` quantitative-trait smoke runs.
- Step 2 no longer calls dense dummy LMM; it performs REGENIE-style residualized linear regression by chromosome/LOCO prediction.
- Step 1 can be restricted to model SNPs with `--regenie-model-snps <snplist>`.
- LOCO predictions can be cached with `--write-pred <file>` and reused with `--use-pred <file>`, allowing Step 2 reruns without retraining ridge.
- PGEN/BGEN biobank Step 2 is intentionally blocked until the PLINK path is validated.
- Accuracy against exact LMM is not expected on tiny toy data if Step 1 trains on all SNPs or if chromosomes are not grouped naturally; production use needs an LD-pruned model-SNP set and sorted chromosome-aware LOCO blocks.

Near-term upgrades required before treating `--biobank-mode` as production:

- Add chromosome/chunk filters for Step 2 distributed runs.
- Add block-size guidance and auto-tuning; small blocks can be much faster for small/medium tests.
- Add calibration benchmarks against REGENIE/BOLT/fastGWA, measuring GC lambda, top-hit overlap, and p-value rank correlation, not exact LMM equality.

Example smoke workflow:

```bash
cosmicscanner --bfile data --pheno pheno.txt --gwas --biobank-mode \
  --regenie-model-snps model.snplist --write-pred data.pred --out step1_and_step2

cosmicscanner --bfile data --pheno pheno.txt --gwas --biobank-mode \
  --use-pred data.pred --out step2_only
```
