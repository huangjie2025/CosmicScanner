# PSR-GWAS v0.5 Stage 2A — Dense ZKZ' BLUP First Pass

## Status: Stage 2A — Experimental

The `stacked-psr` mode now supports a **dense ZKZ' BLUP first pass** engine
(`--ordinal-stacked-engine dense-zkz`), building on Stage 1's fixed-VC weighted
GLS foundation. This is Stage 2A — `sigma2_g` remains fixed at 0.5, and no
variance component estimation is implemented yet.

## What Was Added in Stage 2A

| Component | Implementation |
| --- | --- |
| `--ordinal-stacked-engine` CLI | `fixed-gls` (default, Stage 1) or `dense-zkz` (Stage 2A) |
| `K_stack` construction | `K_stack[r,s] = G[row2indiv[r], row2indiv[s]]` — full `n*(C-1) x n*(C-1)` kinship, NOT block-diagonal |
| PSR GLS update | `V = sigma2_g * K_stack + diag(1/w)`; `beta = (X'V^{-1}X)^{-1} X'V^{-1} y_tilde` via LDLT |
| Wald scan | `P = V^{-1} - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}`; per-SNP `gtPg` and `gtPy` using precomputed LDLT factor |
| BLUP output | `u_hat = sigma2_g * G * Z' * V^{-1} * (y_tilde - X * beta_hat)` written to `.ordinal.blup.tsv` (FID, IID, blup, sigma2_g) |
| N_stack guard | Error if `N = n * (C-1) > 2000` (dense LDLT memory/time limit) |
| sigma2_g | Fixed at 0.5 (same as Stage 1) |
| VC update | Not implemented (Stage 2B) |

## Output Files

| File | Content |
| --- | --- |
| `.gwas.tsv` | `ordinal_implementation = dense_zkz_blup_first_pass` |
| `.psr.convergence.tsv` | Per-iteration sigma_g2, quasi-likelihood, weights, clipping |
| `.ordinal.thresholds.tsv` | Threshold estimates theta_k |
| `.ordinal.blup.tsv` | Individual-level BLUP (FID, IID, blup, sigma2_g) |

## Comparison: Stage 1 vs Stage 2A

| Aspect | Stage 1 (fixed-gls) | Stage 2A (dense-zkz) |
|--------|---------------------|----------------------|
| Variance structure | `V = diag(1/w)` | `V = sigma2_g * K_stack + diag(1/w)` |
| Beta update | WLS: beta = (X'WX)^{-1} X'Wy_tilde | GLS: beta = (X'V^{-1}X)^{-1} X'V^{-1}y_tilde |
| Wald projection | `P_w = I - X_w(X_w'X_w)^{-1}X_w'` (weighted) | `P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}` (full) |
| BLUP | None | `u_hat = sigma2_g * G * Z' * V^{-1} * (y - X*beta)` |
| Implementation label | `fixed_vc_weighted_gls_first_pass` | `dense_zkz_blup_first_pass` |
| K_stack | Not built | Built as `n*(C-1) x n*(C-1)` dense |
| Runtime | Fast (no K_stack) | Slower (LDLT factorization + per-SNP solve) |
| Sample limit | None practical | N_stack <= 2000 |

## Validation

| Test | Status |
| --- | --- |
| dense-zkz basic (N=48, M=36, levels=3) | Finite effect/se/p, BLUP finite |
| dense-zkz + covariate (N=40, M=20, levels=3) | Finite effect/se/p, BLUP finite |
| fixed-gls on same data | Unchanged, all finite |
| fixed-gls vs dense-zkz comparison | log10p correlation, top10 overlap |
| Schema check | `ordinal_implementation=dense_zkz_blup_first_pass` |
| N_stack > 2000 guard | Correctly errors |
| stacked smoke (21 tests) | All pass (fixed-gls + dense-zkz) |
| ordinal validation | cumulative-binary + stacked fixed-gls + stacked dense-zkz |

## Limitations

1. **sigma2_g fixed at 0.5** — not estimated from data. This is the key
   limitation of Stage 2A. The BLUP and Wald statistics use `sigma2_g = 0.5`
   regardless of the true genetic variance. Stage 2B adds grid/profile REML
   estimation to address this.
2. **N_stack <= 2000** — the dense `(n*(C-1)) x (n*(C-1))` LDLT factorization
   limits applicability to small/medium datasets. For n=500 with C=5,
   N_stack=2000 is at the boundary. For UK Biobank-scale data, use fixed-gls.
3. **No VC update** — `sigma2_g` is not iteratively re-estimated within the
   PSR loop. This may affect convergence and BLUP accuracy when the true
   heritability differs substantially from 0.5.
4. **Experimental** — Stage 2A is marked `experimental` and may change in
   future releases. Always validate results against Stage 1 or R baselines.

## Stage 2B (Completed — 2026-06-23)

Stage 2B added grid/profile REML sigma2_g estimation (see `docs/PSR_GWAS_V0.5_STAGE2B_VC_UPDATE_SUMMARY.md`). The `--ordinal-vc-method grid-reml` option now estimates sigma2_g within each PSR iteration via grid search on the REML objective. `--ordinal-vc-method fixed` preserves Stage 2A behavior exactly.

## Stage 2C Roadmap (Planned)

- AI-REML or faster optimization for sigma2_g estimation
- Sparse/low-rank K_stack approximation for large N
- Better conditioning diagnostics (jitter for near-singular V)
- Wald vs Score comparison across stages

## Usage

```bash
# Stage 2A dense-zkz
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz

# With covariates
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --qcovar data/qcovar.txt
```

## Files Changed (Stage 2A)

```
src/config.h                           — ordinal_stacked_engine field
src/main.cpp                           — --ordinal-stacked-engine CLI + help
src/models/ordinal_model.h             — Options + writeStackedBlup declaration
src/models/ordinal_model.cpp           — dense-zkz PSR/GLS/Wald/BLUP implementation
tests/run_psr_ordinal_stacked_smoke.sh — dense-zkz basic + covariate tests
tests/run_psr_ordinal_dense_zkz_validation.sh — standalone dense-zkz validation (new)
tests/run_psr_ordinal_validation.sh    — added stacked fixed-gls + dense-zkz
tests/run_psr_full_parity_validation.sh — added dense-zkz validation step
tests/run_gwas_core_regression.sh      — added dense-zkz under --with-psr
tests/run_psr_full_ordinal_design_check.sh — status=stage2a_dense_zkz_blup
tests/run_psr_formal_simulation.py     — PSR_FORMAL_ORDINAL_ENGINE support
tests/run_psr_accuracy_speed_benchmark.py — stacked engine rows in report
docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md   — Stage 2A update
docs/PSR_GWAS_ROADMAP.md               — Stage 2A/2B roadmap update
docs/PSR_GWAS_V0.5_STACKED_ORDINAL_SUMMARY.md — Stage 2A section
docs/PSR_GWAS_BENCHMARK_GUIDE.md       — updated for stacked engines
docs/PSR_GWAS_V0.5_STAGE2A_DENSE_ZKZ_SUMMARY.md — this document
```
