# PSR-GWAS Benchmark Guide

The v0.4 benchmark framework compares CosmicScanner PSR-GWAS against deterministic R baselines and optional GLMM tools when available.

## Commands

```bash
bash tests/run_psr_gwas_validation.sh
bash tests/run_psr_normal_special_case.sh
bash tests/run_psr_ordinal_validation.sh
python3 tests/run_psr_accuracy_speed_benchmark.py
python3 tests/compare_psr_benchmarks.py
python3 tests/run_psr_formal_simulation.py          # v0.4.1+
bash tests/run_psr_full_ordinal_design_check.sh      # v0.4.1+
bash tests/run_psr_full_parity_validation.sh         # all steps combined
```

## Baselines

| Trait | Baseline |
| --- | --- |
| normal | Existing CosmicScanner LMM equivalence path |
| binary | R `glm`; optional `logistf`, `GMMAT` availability |
| binomial | R `glm(cbind(event, trial-event))`; optional `lme4`, `glmmTMB` availability |
| poisson | R `glm(..., family=poisson)`; optional `lme4`, `glmmTMB` availability |
| ordinal | Optional `MASS::polr`, `ordinal::clm`, `POLMM` availability |

Missing R packages are written as `method_available=FALSE` with `error_message`; they do not fail CI.

## Metrics

The main summary includes `runtime_sec`, `beta_corr`, `se_corr`, `wald_corr`, `log10p_corr`, `top10_overlap`, `top20_overlap`, `top_peak_match`, `lambda_gc`, `converged_rate`, `status_class`, and method availability fields.

Ordinal rows support three modes:
- `experimental_cumulative_binary` (v0.4 MVP, default): independent binary PSR per threshold
- `fixed_vc_weighted_gls_first_pass` (Stage 1): stacked-psr with V=I
- `dense_zkz_blup_first_pass` (Stage 2A, experimental): stacked-psr with V = sigma2_g * K_stack + diag(1/w), produces individual BLUP

Use `--ordinal-stacked-engine fixed-gls|dense-zkz` to select between stacked engines.

Stage 2B adds `--ordinal-vc-method fixed|grid-reml` to control variance component estimation:
- `fixed` (default): sigma2_g stays at the initial value (`--ordinal-sigma2-g`, default 0.5).
- `grid-reml`: estimates sigma2_g via grid/profile REML within each PSR iteration. Only works with `--ordinal-stacked-engine dense-zkz`.

Additional grid params: `--ordinal-vc-grid-min <double>` (default 0.0), `--ordinal-vc-grid-max <double>` (default 2.0), `--ordinal-vc-grid-points <int>` (default 41).

## v0.5 Additions

### R Ordinal Baseline
```bash
bash tests/run_psr_v05_r_baseline_validation.sh
```
Runs `MASS::polr` and `ordinal::clm` as fixed-effect ordinal GWAS baselines. Missing R or missing packages are recorded, not fatal.

### Ordinal Benchmark Runner
```bash
bash tests/run_psr_v05_validation_benchmark.sh
```
Benchmarks all 4 ordinal methods (cumulative-binary, fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml) with configurable N, M, reps. Outputs TSV summary + markdown report.

### Large-N Calibration (Stacked PSR Conservativeness)
```bash
bash tests/run_psr_v05_large_n_calibration_smoke.sh
```
Tests whether stacked PSR lambda GC ≈ 0 resolves at larger N. Sweeps across N values in N_LIST. See `docs/PSR_GWAS_V0.5_LARGE_N_CALIBRATION_PLAN.md` for details.

### Cross-commit Comparison
```bash
python3 tests/compare_psr_benchmarks.py --old <old_dir> --new <new_dir>
```

### Current Known Issue
Stacked PSR methods show lambda GC ≈ 0 at N=200. The large-N calibration branch (`psr-gwas-v05-large-n-calibration`) investigates whether this resolves at larger sample sizes. Stage 2C AI-REML/Brent is not expected to fix this — the issue is in PQL approximation, not null model solver.
