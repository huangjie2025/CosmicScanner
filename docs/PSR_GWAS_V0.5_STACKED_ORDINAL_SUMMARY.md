# PSR-GWAS v0.5 Stacked Ordinal Summary

## Status: Stage 1 — Fixed-VC Weighted GLS First Pass

The `stacked-psr` mode has advanced from **design-only** (v0.4.1) to a **first-pass runnable C++ implementation** (v0.5). This is explicitly Stage 1 / first pass — not the full random-effect ordinal PSR-GLMM.

## What Was Implemented (Stage 1)

| Component | Implementation |
| --- | --- |
| Stacked data builder | Builds `z_stack` (length `n * (C-1)`) with cumulative indicators `I(y <= k)`, design matrix `X_stack` with threshold-dummy encoding |
| Fixed-VC PSR/PQL loop | Reweighted IRLS using `V = I` (identity working covariance), `sigma2_g = 0.5` fixed |
| Per-SNP Wald scan | `P = I - X (X'X)^{-1} X'` on the sqrt-weighted scale; each SNP's genotype is stacked across C-1 thresholds per individual |
| Covariate support | Reads `--qcovar` / `--covar`; places covariates in X_stack alongside threshold dummies |
| Fallback logic | `--ordinal-mode stacked-psr` fails with clear error if internal failure; `--ordinal-allow-fallback` re-routes to cumulative-binary |
| Output files | `.gwas.tsv` (header + per-SNP row with `ordinal_implementation=fixed_vc_weighted_gls_first_pass`), `.psr.convergence.tsv`, `.ordinal.thresholds.tsv` |

## Current Limitations (Documented)

- `sigma2_g = 0.5` fixed — no variance component estimation (even in Stage 2A; Stage 2B adds grid REML)
- `K_stack = Z K Z'` is now implemented in the `dense-zkz` engine, but `sigma2_g` is still fixed
- BLUP is now output for the `dense-zkz` engine: `u_hat = sigma2_g * G * Z' * V^{-1} * (y_tilde - X * beta_hat)`
- IRLS converges immediately when no covariates are present (initial GLS equals starting values for intercept-only threshold model)
- The output `ordinal_fallback_used` field is always `FALSE` when stacked succeeds

These limitations are the design boundary of **Stage 1 (fixed-VC weighted GLS first pass) and Stage 2A (dense ZKZ' BLUP first pass)**. Stage 2B added grid/profile REML variance component estimation for sigma2_g. Stage 2C will add AI-REML and large-sample support.

## Validation

| Test | Status |
| --- | --- |
| Basic smoke (N=60, M=50, levels=4) | 50/50 SNPs finite effect/se/p-value |
| Covariate smoke (N=60, M=30, levels=3, 1 qcovar) | 30/30 SNPs finite |
| Fallback error path | Correctly errors on bad data |
| No-NaN regression | All effect/se/wald/p_value/log10p are finite across all scenarios |
| Cumulative-binary unchanged | `--ordinal-mode experimental-cumulative-binary` continues to work |

## Stage 2A (Completed — 2026-06-23)

- [x] `--ordinal-stacked-engine dense-zkz` CLI option
- [x] `K_stack = Z K Z'` — full `n*(C-1) x n*(C-1)` stacked kinship, directly constructed as `K_stack[r,s] = K[row2indiv[r], row2indiv[s]]`
- [x] `V = sigma2_g * K_stack + diag(1/w)` in PSR iteration for GLS beta update
- [x] Full `P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}` Wald scan using LDLT Cholesky
- [x] Individual-level BLUP: `u_hat = sigma2_g * G * Z' * V^{-1} * (y_tilde - X * beta_hat)`
- [x] Output `dense_zkz_blup_first_pass` implementation label
- [x] `N_stack <= 2000` guard
- [x] Validation wrapper `tests/run_psr_ordinal_dense_zkz_validation.sh`
- [x] Parity matrix update (`run_psr_full_ordinal_design_check.sh` now reports `stage2a_dense_zkz_blup`)

## Stage 2B (Completed — 2026-06-23)

- [x] `--ordinal-vc-method fixed|grid-reml` CLI option
- [x] `--ordinal-sigma2-g <double>` initial/override sigma2_g value
- [x] `--ordinal-vc-grid-min <double>`, `--ordinal-vc-grid-max <double>`, `--ordinal-vc-grid-points <int>` grid search params
- [x] Grid/profile REML sigma2_g estimation within each PQL iteration
- [x] `.ordinal.vc.tsv` output (VC log per iteration)
- [x] `vc_method`, `vc_update`, `sigma2_g` fields in `.gwas.tsv`
- [x] `ordinal_stacked_engine`, `vc_method`, `vc_update` fields in `.psr.convergence.tsv`
- [x] sigma2_g = 0 boundary allowed
- [x] BLUP uses updated sigma2_g
- [x] `--ordinal-stacked-engine fixed-gls` with `--ordinal-vc-method grid-reml` error
- [x] Not AI-REML (Stage 2C)
- [x] Only sigma2_g estimated (sigma2_e via R_stack, fixed)
- [x] Validation: 15/15 tests pass

## Stage 2C (Planned)

- AI-REML or faster optimization for sigma2_g estimation
- Sparse/low-rank K_stack approximation for large N
- sigma2_e estimation
- Better conditioning / jitter diagnostics for dense V
- Accuracy comparison against cumulative-binary and R `MASS::polr`
- Speed benchmark against Stage 1/2A/2B and cumulative-binary

## Files Changed in v0.5

```
src/models/ordinal_model.h     — declare runStackedPsr, buildStackedData, buildGrmMatrix, etc.
src/models/ordinal_model.cpp   — full stacked PSR implementation (Stage 1)
src/gwas_app.cpp               — include ordinal_model.h
src/config.h                   — ordinal_mode, ordinal_allow_fallback fields
src/main.cpp                   — --ordinal-mode, --ordinal-allow-fallback CLI parsing
tests/run_psr_ordinal_stacked_smoke.sh — new: basic, covariate, and fallback tests
.gitignore                     — ignore generated PSR-GWAS validation artifacts
```

## Usage

```bash
# Stage 1 stacked-psr (fixed-VC weighted GLS)
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr

# Stage 2A stacked-psr (dense ZKZ' BLUP) — experimental
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz

# With covariates
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr --qcovar data/qcovar.txt

# Fallback to cumulative-binary if stacked fails
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr --ordinal-allow-fallback
```
