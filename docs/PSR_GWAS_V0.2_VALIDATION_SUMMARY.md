# PSR-GWAS v0.2 Validation Summary

This document records the current PSR-GWAS v0.2 validation state for CosmicScanner.

## Supported Trait Types

- Binary threshold GWAS: `--threshold-model --link probit|logit`
- Poisson count GWAS: `--count-model --family poisson`, optional `--offset`
- Binomial event/trial GWAS: `--count-model --family binomial --binomial-pheno FID IID EVENT TRIAL`

## PSR/PQL Parameters

- `--psr-maxit <n>`: maximum PSR/PQL iterations, default `50`
- `--psr-tol <f>`: convergence tolerance, default `1e-6`
- `--psr-mu-eps <f>`: `mu` clipping epsilon, default `1e-6`
- `--psr-min-weight <f>`: minimum working weight, default `1e-8`
- `--psr-max-weight <f>`: maximum working weight, default `1e8`

## Validation Commands

```bash
make -C /disk228/huangjie/PanOmics/CosmicScanner/build cosmicscanner -j4
bash /disk228/huangjie/PanOmics/CosmicScanner/tests/run_psr_gwas_validation.sh
bash /disk228/huangjie/PanOmics/CosmicScanner/tests/run_gwas_core_regression.sh --fast --with-psr
```

All three commands passed on 2026-06-23.

## Weighted LMM Equivalence

`tests/test_weighted_lmm_equivalence.py` verifies that working weights are handled through the weighted LMM transform rather than an unweighted pseudo phenotype scan.

Latest smoke output:

| Test | Max Absolute Difference | Status |
|---|---:|---|
| `weighted_gls` | `8.88e-16` | pass |
| `weighted_lmm_fixed_v` | `4.44e-16` | pass |
| `weighted_lmm_scan_fixed_vc` | `4.44e-16` | pass |

## Smoke Benchmark Summary

The validation script compares CosmicScanner PSR-GWAS against a local R GLM baseline for smoke-level ranking and effect consistency.

Latest smoke output:

| Trait | Overlap | beta corr | SE corr | -log10(P) corr | Top20 overlap |
|---|---:|---:|---:|---:|---:|
| binary | 240 | 0.99958 | 0.99176 | 0.99955 | 20 |
| poisson | 240 | 0.98707 | 0.98158 | 0.97117 | 17 |
| binomial | 240 | 0.99958 | 0.99176 | 0.99955 | 20 |

## Convergence Diagnostics

Each PSR run writes `<out>.psr.convergence.tsv`. The benchmark summarizer also writes `psr_convergence_summary.tsv` with a `status_class` field:

- `OK`: converged and no warning flags.
- `CONVERGED_WITH_WARNING`: converged, but warning flags are present.
- `NOT_CONVERGED`: did not meet the convergence criterion.

Warning flags currently include:

- `WARNING_PSR_MAXIT_REACHED`: iteration limit reached.
- `WARNING_PSR_CLIPPED`: `mu` or working weight clipping occurred.
- `WARNING_PSR_SIGMA_BOUNDARY`: final genetic variance is close to the lower boundary.
- `WARNING_PSR_QUASILIK_DECREASED`: quasi-likelihood proxy decreased in at least one iteration.

Converged runs with warnings are not automatic failures, but formal analyses should inspect `.psr.convergence.tsv` and `psr_convergence_summary.tsv`.

## Simulation Smoke

`tests/run_psr_simulation_summary.py` now performs actual CosmicScanner PSR-GWAS runs and writes:

- `psr_type1_summary.tsv`
- `psr_power_summary.tsv`
- `psr_qq_points.tsv`

The output includes `engine=cosmicscanner_psr` to distinguish it from earlier simplified score-smoke summaries.

## Known Limits

- `--family negbin` remains experimental and is not recommended for production claims.
- Official PSR-GLMM/R is not integrated because the reference website timed out from this development host and license status is not confirmed.
- Candidate refinement is still script-based. The default method is R `glm`; binary `method=firth` is available only when the optional `logistf` package is installed.
- Current simulation smoke is designed for regression monitoring, not publication-grade Type I error or power claims.

## Next Steps

- v0.3: upgrade candidate refinement outputs and optional methods while keeping the workflow script-based.
- Later: add larger paper-level simulations for case-control imbalance, low-count Poisson traits, and varying-trial binomial traits.
- Defer ordinal, repeated-measure, and random-regression PSR-GWAS until binary/binomial/Poisson validation is stable.
