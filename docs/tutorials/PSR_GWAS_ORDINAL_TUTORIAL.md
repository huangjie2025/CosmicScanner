# PSR-GWAS Ordinal Traits Tutorial

## Overview

PSR-GWAS supports ordinal traits (3+ ordered categories) through two modes: the **experimental cumulative-binary** approximation (v0.4) and the **stacked proportional-odds PSR** pipeline (v0.5). For v0.5 production use, prefer the dense stacked PSR path with grid/profile REML.

Recommended v0.5 command:

```bash
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --out results/ordinal_reml
```

The default `.gwas.tsv` statistic is the corrected Wald / pseudo-score test. Use `--ordinal-test all` only when you need `<out>.ordinal.score.tsv` diagnostics; use `--ordinal-test robust-score` only for cluster-robust sensitivity analysis.

## Ordinal Methods Reference

### 1. cumulative-binary (Baseline)
- **CLI**: `--ordinal-mode experimental-cumulative-binary` (default)
- **Engine**: Runs C-1 separate binary PSR scans (one per cumulative threshold y <= k), merges by min(P)
- **No stacked PSR** — each threshold has its own SNP effect
- **Pros**: Fast, simple, no assumptions about variance components
- **Cons**: No proportional-odds constraint, reduced power, no BLUP

### 2. fixed-gls (Fastest Stacked PSR)
- **CLI**: `--ordinal-mode stacked-psr --ordinal-stacked-engine fixed-gls`
- **Engine**: Stage 1 stacked PSR with V=I (identity working covariance), sigma2_g=0.5 fixed
- **Weighted GLS** scan on stacked pseudo-response
- **Pros**: Fastest ordinal GWAS, no dense matrix factorization
- **Cons**: No BLUP, no sigma2_g estimation

### 3. dense-zkz-fixed (BLUP Available)
- **CLI**: `--ordinal-mode stacked-psr --ordinal-stacked-engine dense-zkz --ordinal-vc-method fixed`
- **Engine**: Stage 2A — full K_stack = Z K Z' construction, V = sigma2_g * K_stack + diag(1/w)
- **Full Wald scan** using LDLT Cholesky of V
- **BLUP output**: `.ordinal.blup.tsv` with individual-level BLUP coefficients
- **Pros**: BLUP available, sigma2_g fixed at 0.5
- **Cons**: N_stack <= 2000 guard, slower than fixed-gls

### 4. dense-zkz-grid-reml (Full REML)
- **CLI**: `--ordinal-mode stacked-psr --ordinal-stacked-engine dense-zkz --ordinal-vc-method grid-reml`
- **Engine**: Stage 2B — grid/profile REML sigma2_g estimation within each PSR iteration
- **VC output**: `.ordinal.vc.tsv` with per-iteration variance component log
- **Full Wald scan + BLUP with updated sigma2_g**
- **Pros**: REML-estimated sigma2_g, BLUP available
- **Cons**: Moderate speed cost (grid search over M=41+ points per iteration)

## How to Choose Between Methods

| Criterion | Recommended Method |
|-----------|-------------------|
| Baseline comparison | cumulative-binary |
| Quick genome-wide scan | fixed-gls |
| BLUP coefficients needed | dense-zkz-fixed or dense-zkz-grid-reml |
| Unknown heritability | dense-zkz-grid-reml |
| Production v0.5 ordinal scan | dense-zkz-grid-reml |
| Score diagnostics | dense-zkz-grid-reml + `--ordinal-test all` |
| N*(C-1) > 2000 | fixed-gls (dense-zkz engine unavailable) |

## Examples

```bash
# Baseline: cumulative-binary (default)
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --out results/ordinal_baseline

# Fast stacked PSR scan
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --out results/ordinal_fast

# With BLUP and fixed variance
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method fixed \
  --out results/ordinal_blup

# Full REML variance estimation
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --out results/ordinal_reml

# Score diagnostics (keeps .gwas.tsv as Wald; writes .ordinal.score.tsv)
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --ordinal-test all \
  --out results/ordinal_score_diagnostics

# With covariates
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --qcovar data/qcovar.txt \
  --out results/ordinal_covar

# Fallback to cumulative-binary if stacked fails
./cosmicscanner --bfile data/geno --gwas \
  --ordinal-model --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-allow-fallback \
  --out results/ordinal_fallback
```

## Output Files

| File | Description |
|------|-------------|
| `.gwas.tsv` | Per-SNP association results. `ordinal_implementation` field identifies the algorithm used: `fixed_vc_weighted_gls_first_pass` (fixed-gls), `dense_zkz_blup_first_pass` (dense-zkz), or `experimental_cumulative_binary` |
| `.psr.convergence.tsv` | PSR/PQL iteration diagnostics per threshold (cumulative-binary) or per stacked iteration |
| `.ordinal.thresholds.tsv` | Threshold/intercept estimates theta_k |
| `.ordinal.blup.tsv` | Individual-level BLUP coefficients (dense-zkz engine only) |
| `.ordinal.vc.tsv` | Per-iteration variance component log (grid-reml only) |
| `.ordinal.score.tsv` | Optional score diagnostics when `--ordinal-test score`, `robust-score`, or `all` is used |

## CLI Parameter Reference

| Parameter | Values | Default | Description |
|-----------|--------|---------|-------------|
| `--ordinal-mode` | `experimental-cumulative-binary`, `stacked-psr` | `experimental-cumulative-binary` | Ordinal algorithm |
| `--ordinal-stacked-engine` | `fixed-gls`, `dense-zkz` | `fixed-gls` | Engine for stacked PSR |
| `--ordinal-vc-method` | `fixed`, `grid-reml` | `fixed` | VC estimation method |
| `--ordinal-sigma2-g` | float | 0.5 | Initial/fixed sigma2_g value |
| `--ordinal-vc-grid-min` | float | 0.0 | Grid search lower bound |
| `--ordinal-vc-grid-max` | float | 2.0 | Grid search upper bound |
| `--ordinal-vc-grid-points` | int | 41 | Grid search points |
| `--ordinal-test` | `wald`, `score`, `robust-score`, `all` | `wald` | Primary test / diagnostics for stacked PSR. `all` writes all diagnostics while keeping Wald as `.gwas.tsv` primary output |
| `--ordinal-allow-fallback` | flag | — | Fall back to cumulative-binary if stacked fails |

## Validation Scripts

### Formal Simulation
```bash
# Run all 4 ordinal methods in formal simulation
PSR_FORMAL_TRAITS=ordinal PSR_FORMAL_ORDINAL_METHODS=all \
  PSR_FORMAL_REPS=5 PSR_FORMAL_N=200 \
  python3 tests/run_psr_formal_simulation.py
```

### Standalone Benchmark
```bash
# Run ordinal benchmark with integrity checks
PSR_V05_BENCH_REPS=5 PSR_V05_BENCH_N=200 \
  bash tests/run_psr_v05_validation_benchmark.sh
```

### Score-Test Regression
```bash
bash tests/run_psr_ordinal_score_test_validation.sh
bash tests/run_psr_v05_score_calibration_smoke.sh
```

### Compare Across Commits
```bash
python3 tests/compare_psr_benchmarks.py \
  --old tests/results/psr_v05_ordinal_benchmark_OLD \
  --new tests/results/psr_v05_ordinal_benchmark
```

## References

- `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md` — v0.5 validation report (includes Stage 1/2A/2B details and landing calibration)
- `docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md` — Full ordinal PSR design document
- `docs/PSR_GWAS_ROADMAP.md` — PSR-GWAS roadmap and acceptance criteria
- `docs/Tutorial_EN.md` — General CosmicScanner tutorial
