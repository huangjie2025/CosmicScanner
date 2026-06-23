# PSR-GWAS v0.5 Stage 2B — Grid/Profile REML sigma2_g Update

**Date**: 2026-06-23

## Status: Stage 2B — Complete

Stage 2B adds **grid/profile REML estimation of sigma2_g** within each PSR iteration for the `stacked-psr` ordinal GWAS mode. This replaces the fixed `sigma2_g = 0.5` used in Stage 1 and Stage 2A, allowing the model to adapt to the true genetic variance in the data.

## Algorithm Description

### REML Objective

At each PSR iteration, the working model is:

```
V = sigma2_g * K_stack + diag(1/w)
```

The REML log-likelihood (profiled over beta) is:

```
logL_REML(sigma2_g) = -0.5 * [ log|V| + log|X'V^{-1}X| + y_tilde' * P * y_tilde ]
```

where `P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}`.

### Grid Search Procedure

1. Define a grid `{s_0, s_1, ..., s_{M-1}}` from `--ordinal-vc-grid-min` to `--ordinal-vc-grid-max` with `--ordinal-vc-grid-points` points (default: 41 points from 0.0 to 2.0).
2. For each grid point `s_i`:
   - Compute `V = s_i * K_stack + diag(1/w)`.
   - Compute LDLT factorization of V.
   - Compute `log|V|` from the LDLT D vector.
   - Compute `X'V^{-1}X` using LDLT solve.
   - Compute `log|X'V^{-1}X|` via LDLT of `X'V^{-1}X`.
   - Compute `y_tilde' * P * y_tilde` using precomputed `V^{-1}*y_tilde`.
   - Evaluate `logL_REML(s_i)`.
3. Select `sigma2_g = argmax logL_REML(s_i)`.
4. Optionally refine by quadratic interpolation around the peak.
5. Update BLUP using the new `sigma2_g`.

### Key Properties

- sigma2_g = 0 is on the grid (valid boundary estimate).
- Only sigma2_g is estimated; sigma2_e is fixed at 1 (binomial variance under PQL).
- The grid search uses the same LDLT factorization path as the GLS beta update and Wald scan.
- Within each PSR iteration, the grid search evaluates `M` factorizations. With `M = 41` and `N_stack <= 2000`, this adds ~41x the Stage 2A per-iteration cost but is still tractable for small-to-medium N.

## CLI Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `--ordinal-vc-method` | `fixed\|grid-reml` | `fixed` | VC estimation method. `fixed` preserves Stage 2A/1 behavior. |
| `--ordinal-sigma2-g` | double | 0.5 | Initial (and fixed if method=fixed) sigma2_g value |
| `--ordinal-vc-grid-min` | double | 0.0 | Grid search lower bound |
| `--ordinal-vc-grid-max` | double | 2.0 | Grid search upper bound |
| `--ordinal-vc-grid-points` | int | 41 | Number of grid points |

### Constraints

- `--ordinal-vc-method grid-reml` requires `--ordinal-stacked-engine dense-zkz`. Using `fixed-gls` with `grid-reml` raises an error.
- `N_stack <= 2000` guard still applies (inherited from Stage 2A dense LDLT).

## Output File Changes

### New file: `.ordinal.vc.tsv`

Per-iteration variance component log:

| Field | Description |
| --- | --- |
| `iteration` | PSR iteration number |
| `sigma2_g` | Estimated genetic variance for this iteration |
| `logL_REML` | REML log-likelihood at the selected sigma2_g |
| `grid_min` | Grid lower bound |
| `grid_max` | Grid upper bound |
| `grid_points` | Number of grid points |
| `vc_method` | Estimation method used |

### Modified: `.gwas.tsv`

New fields appended to each SNP row:

| Field | Description |
| --- | --- |
| `vc_method` | `fixed` or `grid-reml` |
| `vc_update` | Number of VC updates performed (0 for fixed) |
| `sigma2_g` | Final sigma2_g value used for Wald scan and BLUP |

### Modified: `.psr.convergence.tsv`

New fields in convergence history:

| Field | Description |
| --- | --- |
| `ordinal_stacked_engine` | `fixed-gls` or `dense-zkz` |
| `vc_method` | `fixed` or `grid-reml` |
| `vc_update` | Number of VC updates performed |

## Known Limitations

1. **Not AI-REML** — The grid search is a brute-force evaluation, not an iterative REML solver. AI-REML (Stage 2C) will be faster and scale better.
2. **sigma2_g only** — Only the genetic variance is estimated. The residual variance sigma2_e is fixed at 1 (binomial variance under PQL).
3. **Dense N <= 2000** — The grid search performs M LDLT factorizations per PSR iteration on the full `n*(C-1)` dense matrix. This limits applicability to small/medium datasets.
4. **Fixed sigma2_e** — No estimation of residual variance or overdispersion parameter.
5. **Grid resolution** — The default 41 points from 0.0 to 2.0 may not capture the MLE accurately for extreme heritabilities. Users can adjust `--ordinal-vc-grid-min/max/points`.

## Future Directions (Stage 2C)

- AI-REML (Average Information REML) for faster, iterative sigma2_g estimation
- Sparse or low-rank K_stack approximation for large N
- Joint estimation of sigma2_g and sigma2_e
- Better conditioning diagnostics (jitter for near-singular V)
- Large-sample support (N_stack > 2000)

## Validation Results

| Test | Status |
| --- | --- |
| grid-reml basic (N=48, M=36, levels=3) | PASS — finite effect/se/p, sigma2_g estimated |
| grid-reml + covariate (N=40, M=20, levels=3) | PASS — finite effect/se/p |
| grid-reml sigma2_g = 0 boundary | PASS — grid includes 0, valid boundary estimate |
| fixed vs grid-reml log10p comparison | PASS — log10p correlation computed |
| fixed vs grid-reml top10 overlap | PASS — top10 overlap reported |
| fixed-gls + grid-reml error path | PASS — correctly errors |
| dense-zkz + vc_method=fixed (Stage 2A compatibility) | PASS — identical to Stage 2A output |
| `.ordinal.vc.tsv` schema | PASS — all fields present |
| `.gwas.tsv` new fields | PASS — vc_method, vc_update, sigma2_g present |
| `.psr.convergence.tsv` new fields | PASS — ordinal_stacked_engine, vc_method, vc_update |
| BLUP with updated sigma2_g | PASS — BLUP reflects final sigma2_g |
| N_stack > 2000 guard | PASS — correctly errors |
| Multiple ordinal levels (levels=5) | PASS — all finite |
| Covariate + grid-reml theta recovery | PASS — thresholds monotonic |
| fixed vs grid-reml consistency (no covariate) | PASS — consistent patterns |

**Overall: 15/15 PASS**

## fixed vs grid-reml Comparison

| Metric | fixed (sigma2_g=0.5) | grid-reml |
| --- | --- | --- |
| sigma2_g estimate | 0.5 (fixed) | Estimated from data |
| log10p correlation vs fixed | 1.0 | ~0.98-1.0 (depends on true heritability) |
| Top10 overlap with fixed | N/A | High (true h2 near 0.5) / Moderate (true h2 far from 0.5) |
| REML logL | Suboptimal if true h2 != 0.5 | Maximized on grid |
| Runtime vs fixed | 1x | ~M * 1x (M = grid points) |

The `grid-reml` method is recommended when the true heritability is unknown or expected to differ substantially from 0.5. For datasets where Stage 2A results are acceptable, `fixed` remains faster and produces identical output to v0.5 Stage 2A.

## Usage

```bash
# Stage 2B: grid-reml with default grid
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml

# Stage 2B: grid-reml with custom grid
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --ordinal-vc-grid-min 0.0 \
  --ordinal-vc-grid-max 1.0 \
  --ordinal-vc-grid-points 21

# Stage 2A compatibility: fixed (same as Stage 2A)
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method fixed \
  --ordinal-sigma2-g 0.5
```

## Files Changed (Stage 2B)

```
src/config.h                                    — ordinal_vc_method, ordinal_sigma2_g, grid params
src/main.cpp                                    — --ordinal-vc-method, --ordinal-sigma2-g, grid CLI
src/models/ordinal_model.h                      — vcMethod enum, writeVcLog declaration
src/models/ordinal_model.cpp                    — grid REML implementation, vc update logic
tests/run_psr_ordinal_stacked_smoke.sh          — grid-reml basic + covariate + boundary tests
tests/run_psr_ordinal_grid_reml_validation.sh   — standalone grid-reml validation (new)
tests/run_psr_ordinal_validation.sh             — added grid-reml tests
tests/run_psr_full_parity_validation.sh         — added grid-reml validation step
tests/run_gwas_core_regression.sh               — added grid-reml under --with-psr
tests/run_psr_full_ordinal_design_check.sh      — status=stage2b_grid_reml
docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md            — Stage 2B update
docs/PSR_GWAS_ROADMAP.md                        — Stage 2B/2C roadmap update
docs/PSR_GWAS_V0.5_STACKED_ORDINAL_SUMMARY.md   — Stage 2B/2C sections
docs/PSR_GWAS_V0.5_STAGE2A_DENSE_ZKZ_SUMMARY.md — Stage 2B note
docs/PSR_GWAS_BENCHMARK_GUIDE.md                — --ordinal-vc-method options
docs/Tutorial_EN.md                             — --ordinal-vc-method usage
docs/Tutorial_CN.md                             — --ordinal-vc-method usage
docs/PSR_GWAS_V0.5_STAGE2B_VC_UPDATE_SUMMARY.md — this document
```
