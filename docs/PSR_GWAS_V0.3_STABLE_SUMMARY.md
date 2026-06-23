# PSR-GWAS v0.3 Stable Summary

This stable sprint closes the current binary/binomial/Poisson PSR-GWAS module as a verifiable, regression-tested CosmicScanner component. It does not add ordinal, repeated-measure, random-regression, or exact GLMM models.

## Scope

Supported trait types:

- Binary: `--threshold-model --link probit|logit`
- Poisson: `--count-model --family poisson`, optional `--offset`
- Binomial: `--count-model --family binomial --binomial-pheno FID IID EVENT TRIAL`

Model contract:

- PSR-GWAS is a PSR/PQL approximation, not an exact GLMM.
- The backend is pseudo response + working weights + weighted LMM scan.
- Candidate refinement remains script-based and is not merged into the C++ main program.

## Diagnostics

Each PSR run writes `<out>.psr.convergence.tsv`. The validation summary writes `psr_convergence_summary.tsv` and classifies each trait as:

- `OK`: converged without warning flags.
- `CONVERGED_WITH_WARNING`: converged, but diagnostics such as boundary variance, clipping, or non-monotone quasi-likelihood were observed.
- `NOT_CONVERGED`: convergence criterion was not met; formal analyses should treat results cautiously.

Warning flags include:

- `WARNING_PSR_MAXIT_REACHED`
- `WARNING_PSR_CLIPPED`
- `WARNING_PSR_SIGMA_BOUNDARY`
- `WARNING_PSR_QUASILIK_DECREASED`

## Candidate Refinement

`tests/psr_candidate_refine.R` supports:

- `glm`: default method for binary, binomial, and Poisson.
- `firth`: binary-only, optional `logistf` dependency.
- `glmer`: optional `lme4` dependency.
- `glmmTMB`: optional `glmmTMB` dependency.

Missing optional R packages do not crash the suite. The script writes per-SNP rows with `converged=FALSE` and a clear `error_message`.

Unified output columns:

- `trait_type`
- `method`
- `snp`
- `beta`
- `se`
- `wald_or_stat`
- `p_value`
- `runtime`
- `converged`
- `error_message`

## Simulation Smoke

`tests/run_psr_simulation_summary.py` performs actual CosmicScanner PSR-GWAS calls for binary, Poisson, and binomial traits under null and power scenarios.

Outputs:

- `psr_type1_summary.tsv`
- `psr_power_summary.tsv`
- `psr_qq_points.tsv`
- `method_metrics.tsv`

Each summary includes `engine=cosmicscanner_psr`, seed, trait type, scenario, sample size, SNP count, replicate, runtime, convergence flag, and `status_class`.

## Formal Simulation Skeleton

`tests/run_psr_formal_simulation.py` is a configurable real-CosmicScanner simulation driver for larger validation runs. Defaults are intentionally small; paper-level runs can set:

```bash
PSR_FORMAL_REPS=100 PSR_FORMAL_N=500 PSR_FORMAL_M=2000 python3 tests/run_psr_formal_simulation.py
```

Outputs:

- `formal_type1_summary.tsv`
- `formal_power_summary.tsv`
- `formal_qq_points.tsv`
- `formal_runtime_summary.tsv`
- `formal_convergence_summary.tsv`

The formal driver covers binary balanced/unbalanced scenarios, Poisson low/medium/high mean scenarios, and binomial fixed/variable trial scenarios. It remains separate from CI because paper-level repetitions are intentionally expensive.

## Stable CI Wrapper

`tests/run_psr_stable_ci.sh` wraps the fast core regression path:

```bash
bash tests/run_psr_stable_ci.sh --build
```

This runs `tests/run_gwas_core_regression.sh --fast --with-psr` after an optional build.

## Unified Report

`tests/summarize_psr_gwas_benchmark.py` now writes:

- `psr_gwas_benchmark.summary.tsv`
- `psr_convergence_summary.tsv`
- `psr_validation_report.md`

The unified TSV and Markdown report include convergence diagnostics, R baseline comparison, weighted equivalence results, candidate refinement status, and simulation smoke metrics.

## Input Validation Regression

`tests/test_psr_input_validation.sh` covers:

- Binary traits with one class after filtering.
- Binomial `EVENT > TRIAL`.
- Binomial `TRIAL <= 0`.
- Negative Poisson counts.
- Offset files with too few rows.
- Offset files containing `NaN`/`Inf`.
- `method=firth` without `logistf` package, which must not crash.

## Validation Commands

```bash
make -C /disk228/huangjie/PanOmics/CosmicScanner/build cosmicscanner -j4
bash /disk228/huangjie/PanOmics/CosmicScanner/tests/run_psr_gwas_validation.sh
bash /disk228/huangjie/PanOmics/CosmicScanner/tests/run_gwas_core_regression.sh --fast --with-psr
```

## Remaining Limits

- Negative binomial remains experimental and is not recommended for formal reporting.
- Official PSR-GLMM/R is not integrated because source/license availability has not been confirmed from the development host.
- Current simulations are smoke/regression scale, not publication-grade Type I error and power studies.
- Ordinal, repeated-measure, and random-regression PSR-GWAS remain roadmap items.
# v0.4 Follow-Up Note

The v0.3 stable scope remains binary/binomial/Poisson PSR-GWAS. v0.4 builds on this by adding normal as an LMM special case, an experimental ordinal cumulative-binary MVP, official PSR-GLMM/R availability tracking, and an accuracy/speed benchmark framework. The ordinal MVP is not full PSR-GLMM/R-equivalent ordinal inference.
