# Full Ordinal PSR-GWAS Design

## Status

**Current implementation**: `stacked_psr` Stage 2B (v0.5 grid/profile REML sigma2_g update, landed) + `experimental_cumulative_binary` (v0.4 MVP, default)
**Stage 1** (v0.5, landed): Fixed-VC weighted GLS first pass — `sigma2_g = 0.5` fixed, `V = I`, no random effects
**Stage 2A** (v0.5, landed): Dense `ZKZ'` engine — `V = sigma2_g * K_stack + diag(1/w)`, GLS update, full Wald scan, individual BLUP output
**Stage 2B** (v0.5, landed): Grid/profile REML sigma2_g estimation — `--ordinal-vc-method grid-reml`, VC update within PSR iteration, `.ordinal.vc.tsv`. Wald p-value bug fixed; focused large-N calibration passed (lambda GC ~1.0 at N=200/500).
**Stage 2C** (planned): AI-REML or faster optimization, sparse/low-rank K_stack for large N, sigma2_e estimation

See `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md` for the v0.5 landing calibration results.

## Limitations of Current Experimental Cumulative-Binary

The v0.4 experimental ordinal runs C-1 independent binary PSR scans (one per
cumulative threshold `y <= k`), then merges SNP results by `min(P)`. This has
known limitations:

1. No proportional-odds constraint — each threshold has its own SNP effect
   estimate, so the model does not share information across thresholds.
2. No stacked pseudo-response — working weights and pseudo-responses from
   different thresholds are modelled separately, reducing effective sample size.
3. No joint likelihood — model comparison across thresholds is ad-hoc (minP),
   not a single ordinal likelihood.
4. No threshold separation constraint — in theory thresholds could cross
   (though this is unlikely with the cumulative-binary formulation).

## Full Ordinal PSR Target Model

### Cumulative Logit/Probit

For a C-level ordinal outcome `y_i ∈ {1, ..., C}`:

```
P(y_i <= k) = g^{-1}(theta_k - x_i' * beta - g_i)
```

where:
- `g` is the link function (logit or probit)
- `theta_k` are C-1 threshold/intercept parameters (theta_1 < theta_2 < ... < theta_{C-1})
- `beta` is the shared SNP effect vector (proportional-odds assumption)
- `g_i` is the random polygenic effect: `g ~ N(0, K * sigma_g^2)`

### Stacked Pseudo-Response

The PSR/PQL linearisation creates (C-1) pseudo-observations per individual:

```
y_tilde_ik = eta_ik + (y_ik - p_ik) / (p_ik * (1 - p_ik)) * d_eta_ik / d_mu_ik
```

where `y_tilde_ik` is the working pseudo-response for individual `i` at
threshold `k`, `eta_ik = theta_k - x_i' * beta - g_i` is the linear predictor,
`p_ik = P(y_i <= k)` is the cumulative probability, and `y_ik = I(y_i <= k)` is
the cumulative indicator.

The stacked vector is:

```
y_tilde = [y_tilde_11, ..., y_tilde_1{C-1}, y_tilde_21, ..., y_tilde_n{C-1}]'
```

### Working Weights

Diagonal working weight matrix `W` of size `n*(C-1) x n*(C-1)`, where each
block corresponds to individual `i` at threshold `k`:

```
w_ik = p_ik^2 * (1 - p_ik)^2 / var(y_ik)
```

For binomial cumulative: `var(y_ik) = p_ik * (1 - p_ik)`, so
`w_ik = p_ik * (1 - p_ik)`.

### Kinship for Stacked Observations

Define incidence matrix `Z` of size `n*(C-1) x n` where `Z[j,i] = 1` if
stacked row `j` corresponds to individual `i` (each individual maps to C-1
rows). Then:

```
K_stack = Z * K * Z'
```

Equivalently:

```
K_stack[r, s] = K[row2indiv[r], row2indiv[s]]
```

**K_stack is NOT block-diagonal.** Each of the C-1 thresholds maps to the
same individual effect, so `K_stack` contains (C-1) x (C-1) identical
blocks of `K`. This means the off-diagonal *across-threshold* entries are
the same individual's kinship with themselves, and cross-individual
entries are the full kinship. The matrix has an `n*(C-1)` square shape and
rank `min(n, C-1)`.

The variance structure is:

```
V_stack = K_stack * sigma_g^2 + W^{-1} * sigma_e^2
```

with `sigma_e^2 = 1` (binomial variance fixed to 1 under PQL).

### PSR Iteration

The PSR loop updates `eta`, `p`, `y_tilde`, `W`, and `K_stack` at each
iteration, then solves the mixed model equations:

```
(X' * R^{-1} * X + K_stack^{-1} / sigma_g^2) * beta_hat = X' * R^{-1} * y_tilde
```

where `R = W^{-1}`.

After convergence, the final weighted LMM scan structure is:

```
beta_SNP_hat = (S' * V_stack^{-1} * S)^{-1} * S' * V_stack^{-1} * y_tilde_res
```

where `S` is the stacked SNP dosage vector (replicated C-1 times per
individual), and `y_tilde_res` is the pseudo-response residualised for fixed
effects.

### Output Interpretation

- SNP effects are on the latent (linear predictor) scale.
- A positive `beta` means higher SNP dosage increases the probability of
  being in a higher category (under proportional odds with logit link,
  `exp(beta)` is the cumulative odds ratio).
- Thresholds `theta_k` are nuisance parameters and are reported in the
  thresholds output file.

## Difference from Experimental Cumulative-Binary

| Aspect | Experimental Cumulative-Binary | Full Stacked PSR |
|--------|-------------------------------|-------------------|
| SNP effect | one per threshold (C-1 values) | single shared beta |
| Merging strategy | minP across thresholds | joint likelihood |
| PSR iterations | independent per threshold | single joint PSR loop |
| Effective N | n per threshold | n * (C-1) stacked |
| Kinship | n x n per threshold | n*(C-1) x n*(C-1) stacked |
| Threshold constraint | not enforced | theta_1 < ... < theta_{C-1} |
| link | logit only | logit or probit |
| Output | `ordinal_mode=experimental_cumulative_binary` | `ordinal_mode=stacked_psr` |

## Validation Plan

### R Baseline Methods

- **MASS::polr**: Cumulative logit/probit with proportional odds, no random
  effects. Can serve as fixed-effect baseline for beta comparison.
- **ordinal::clm**: More flexible cumulative link models, supports various
  link functions. No random effects in the base version.
- **ordinal::clmm2**: Cumulative link mixed model with one random intercept.
  Closest R analogue to stacked PSR (but with approximate likelihood).
- **POLMM**: POLMM R package, designed for ordinal GWAS with saddlepoint
  approximation. If available, best external reference.
- **PSR-GLMM/R**: Official PSR-GLMM/R tool (http://106.13.144.206/psrglmm/).
  Not currently reachable; tracked separately.

### Validation Steps

1. **Fixed-effect parity**: Compare stacked PSR beta vs `MASS::polr` and
   `ordinal::clm` (no random effect, small N).
2. **Mixed-model parity**: Compare stacked PSR vs `ordinal::clmm2` and
   optional POLMM for small N with GRM.
3. **Type I / Power**: Same formal simulation as binary/binomial/Poisson but
   with ordinal outcomes. Compare against R baselines.
4. **Threshold recovery**: Check that `theta_k` are ordered correctly and
   increase monotonically.
5. **Proportional-odds diagnostic**: If proportional odds is violated in
   simulated data, check whether stacked PSR still gives reasonable results
   (robustness check).

## Implementation Plan

### Phase 1: Design (this document)
- [x] Document model, stacked structure, variance, output.

### Phase 2: Scaffolding — C++ interface
- [x] Add `--ordinal-mode` CLI argument with values:
  - `experimental-cumulative-binary` (default, existing)
  - `stacked-psr` (new, experimental)
- [x] Add `--ordinal-allow-fallback` flag to allow fallback to experimental
      mode if stacked-psr fails.
- [x] Add `OrdinalModel::Options` with `mode` field.

### Phase 3: Stacked PSR core loop (Stage 1 — fixed-VC, first pass)
- [x] `buildStackedDesign()` — build y_tilde, W, Z, stacked dosage S.
- [x] `buildStackedKinship()` — build K_stack = Z * K * Z' (implemented for Stage 2A).
- [x] `runStackedPSR()` — PSR iteration with fixed sigma_g^2 = 0.5, V = I.
- [x] `runStackedGWAS()` — weighted GLS scan on stacked pseudo-response (P = I - X (X'X)^{-1}X' on sqrtW scale).
- [x] Output `<out>.gwas.tsv` with `ordinal_mode=stacked_psr`, `ordinal_implementation=fixed_vc_weighted_gls_first_pass`.
- [x] Output `<out>.ordinal.thresholds.tsv` with theta_k estimates.
- [x] Output `<out>.psr.convergence.tsv` with iteration history.
- [x] **Stage 2A**: dense ZKZ' engine + BLUP (sigma2_g fixed, 2026-06-23).
- [x] Stage 2B: grid/profile REML sigma2_g estimation, `--ordinal-vc-method grid-reml`, VC update within PSR iteration.
- [ ] Stage 2C: AI-REML or faster optimization, large-sample support.

### Phase 4: Validation (Stage 1 done; Stage 2 pending)
- [ ] Add R baseline parity script using MASS::polr and ordinal::clm (deferred — R packages unavailable).
- [x] Add stacked-psr smoke test with 3+ levels (basic, covariate, fallback).
- [ ] Add stacked-psr vs experimental cumulative-binary comparison (Stage 2).
- [ ] Add formal simulation with ordinal outcomes (Stage 2).

### Phase 5: Documentation (Stage 1 done)
- [x] Update `docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md` with Stage 1 implementation details.
- [x] Update Tutorial docs with stacked-psr usage.
- [x] Stage 1/2A/2B details recorded in `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md`.
- [x] Update parity matrix: ordinal → stacked_psr (Stage 2B).

## Risks

1. **Computational cost**: Stacked V is `n*(C-1) x n*(C-1)`, which is
   (C-1)^2 times larger per iteration than the univariate case. With C=5
   and n=10000, the stacked matrix is 40000 x 40000.
2. **Memory**: K_stack is block-diagonal with C-1 identical blocks, but the
   stacked structure requires storing `n*(C-1)` pseudo-observations and
   weights.
3. **Convergence**: Ordinal PSR can have non-monotone quasi-likelihood,
   especially when thresholds are close together or categories have few
   observations.
4. **Proportional-odds violation**: If SNP effects differ substantially
   across thresholds, a single beta is a poor summary. A diagnostic test
   (comparing threshold-specific betas from experimental cumulative-binary)
   should be provided.
5. **Threshold separation**: If `theta_k ≈ theta_{k+1}`, the cumulative
   binary for threshold k+1 has almost no cases, leading to boundary
   variance estimates.
6. **Weights**: Extreme weights can arise when `p_ik` is near 0 or 1. The
   existing `--psr-min-weight` / `--psr-max-weight` clamps apply.

## Current Status (v0.5)

```
full_ordinal_psr_status = stage2b_grid_reml
```

Stage 1 of stacked-psr is implemented and validated (see `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md`). The implementation uses fixed `sigma2_g = 0.5`, no stacked kinship `K_stack`, and a weighted GLS scan (`V = I`).

Stage 2A (dense `ZKZ'` BLUP first pass) is implemented and validated (see `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md`). The dense-zkz engine (`--ordinal-stacked-engine dense-zkz`) builds the full `K_stack = Z * K * Z'`, uses `V = sigma2_g * K_stack + diag(1/w)` for GLS beta update and full Wald scan, and outputs individual-level BLUP to `.ordinal.blup.tsv`. `sigma2_g` remains fixed at 0.5; VC update (Stage 2B) is not yet implemented.

Stage 2B (grid/profile REML sigma2_g estimation) is implemented and validated (see `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md`). `--ordinal-vc-method grid-reml` performs a grid search over sigma2_g at each PSR iteration to maximize the REML log-likelihood, supporting sigma2_g = 0 boundary. `--ordinal-vc-method fixed` preserves Stage 2A/1 behavior. The Wald p-value bug has been fixed and focused large-N calibration passed (lambda GC ~1.0 at N=200/500).

Stage 2C (AI-REML) is planned but not yet started.

## References

- Breslow & Clayton (1993). Approximate inference in generalized linear mixed
  models. *JASA* 88(421):9-25.
- Wolfinger & O'Connell (1993). Generalized linear mixed models: a
  pseudo-likelihood approach. *J. Statist. Comput. Simul.* 48:233-243.
- Agresti (2010). *Analysis of Ordinal Categorical Data*. Wiley.
- Chen et al. (2023). POLMM: an efficient analysis tool for ordinal GWAS.
  *Bioinformatics* 39(11):btad667.
