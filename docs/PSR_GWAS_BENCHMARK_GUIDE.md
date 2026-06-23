# PSR-GWAS Benchmark Guide

The v0.4 benchmark framework compares CosmicScanner PSR-GWAS against deterministic R baselines and optional GLMM tools when available.

## Commands

```bash
bash tests/run_psr_gwas_validation.sh
bash tests/run_psr_normal_special_case.sh
bash tests/run_psr_ordinal_validation.sh
python3 tests/run_psr_accuracy_speed_benchmark.py
python3 tests/compare_psr_benchmarks.py
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

Ordinal rows are marked `experimental_cumulative_binary` until a full stacked proportional-odds PSR/PQL model is implemented and validated.
