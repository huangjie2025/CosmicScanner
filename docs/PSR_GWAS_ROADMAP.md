# PSR-GWAS Roadmap And Acceptance Criteria

Goal: make binary/binomial/Poisson/ordinal PSR-GWAS statistically auditable and reproducible while keeping repeated-measure PSR-GLMM and random-regression PSR-GWAS as separate model tracks.

## Current Scope

Implemented trait types:

- Binary: `--threshold-model --link probit|logit`
- Poisson: `--count-model --family poisson`, optional `--offset`
- Binomial: `--count-model --family binomial --binomial-pheno FID IID EVENT TRIAL`
- Ordinal: `--ordinal-model`; v0.5 recommended path is `--ordinal-mode stacked-psr --ordinal-stacked-engine dense-zkz --ordinal-vc-method grid-reml`

Current model contract:

- PSR/PQL linearization produces a pseudo response.
- Working residual variance is passed to the scanner as `obs_weights = sqrt(W)`.
- The SNP scan uses weighted LMM, not ordinary LMM on an unweighted pseudo phenotype.
- Candidate refinement is external/scripted for now.

## v0.2a Diagnostics And Input Validation

Implemented acceptance criteria:

- Each PSR run writes `<out>.psr.convergence.tsv`.
- Convergence file includes iteration, genetic variance, quasi-likelihood proxy, mean/min/max `mu`, min/max weight, clipping count, and convergence flag.
- CLI tuning is available through `--psr-maxit`, `--psr-tol`, `--psr-mu-eps`, `--psr-min-weight`, and `--psr-max-weight`.
- `tests/summarize_psr_gwas_benchmark.py` writes `psr_convergence_summary.tsv` and merges convergence status into `psr_gwas_benchmark.summary.tsv`.
- Warning flags currently include `WARNING_PSR_MAXIT_REACHED`, `WARNING_PSR_CLIPPED`, `WARNING_PSR_SIGMA_BOUNDARY`, and `WARNING_PSR_QUASILIK_DECREASED`.
- Convergence summaries classify each trait as `OK`, `CONVERGED_WITH_WARNING`, or `NOT_CONVERGED`.
- `docs/PSR_GWAS_V0.2_VALIDATION_SUMMARY.md` records the current validation commands, outputs, and known limitations.
- Binary traits require both cases and controls after filtering.
- Binomial records require finite `0 <= EVENT <= TRIAL` and `TRIAL > 0`.
- Poisson/negative-binomial traits require finite non-negative phenotype values.
- Offset files must be readable, finite, and have at least one value per sample.

## v0.2b Weighted LMM Equivalence

Implemented smoke criteria:

- `tests/test_weighted_lmm_equivalence.py` verifies weighted GLS equivalence.
- It also verifies fixed-variance weighted LMM equivalence for `V = K sigma_g2 + R` versus transformed `Vw = sqrt(W) K sqrt(W) sigma_g2 + I`.
- It includes a fixed-variance multi-SNP scan smoke to guard the per-SNP weighted LMM path.
- Current tolerance target: absolute difference below `1e-8` for beta, SE, Wald, and p-value.

## v0.2c Simulation Summary MVP

Implemented smoke outputs:

- `psr_type1_summary.tsv`
- `psr_power_summary.tsv`
- `psr_qq_points.tsv`

These are lightweight simulation summaries from actual CosmicScanner PSR-GWAS runs for CI-style trend monitoring. They are not a final publication-grade PSR benchmark. A full method validation should scale the same null/alternative simulations and compare against official PSR-GLMM/R or exact GLMM baselines when available.

## v0.3 Candidate Refinement Upgrade

Current MVP:

- `tests/psr_candidate_refine.R` supports GLM refit for binary, binomial, and Poisson top SNPs.
- The script now emits `method`, `converged`, and `error_message` columns and accepts optional `method=glm|firth`; `firth` is binary-only and requires the optional R package `logistf`.
- Stable sprint update: `tests/psr_candidate_refine.R` now accepts optional `method=glm|firth|glmer|glmmTMB` while keeping `glm` as the default. Missing optional R packages are reported per SNP without crashing.
- Unified candidate refinement output columns are `trait_type`, `method`, `snp`, `beta`, `se`, `wald_or_stat`, `p_value`, `runtime`, `converged`, and `error_message`.
- `tests/run_psr_candidate_refinement_suite.sh` tests the default GLM path and optional dependency degradation paths.
- `tests/test_psr_input_validation.sh` covers invalid PSR inputs and the missing-package refinement degradation path.
- `tests/summarize_psr_gwas_benchmark.py` now emits `psr_validation_report.md` and folds convergence, weighted equivalence, R baseline, refinement, and simulation summaries into `psr_gwas_benchmark.summary.tsv`.
- `docs/PSR_GWAS_V0.3_STABLE_SUMMARY.md` records the stable sprint state and limitations.
- `tests/run_psr_stable_ci.sh` provides a one-command fast CI wrapper around `run_gwas_core_regression.sh --fast --with-psr`.
- `tests/run_psr_formal_simulation.py` provides a configurable real-CosmicScanner formal simulation skeleton for paper-level validation, kept separate from CI-scale smoke tests.

Planned acceptance criteria:

- Add at least two of: Firth logistic, `glmer`, `glmmTMB`, GMMAT adapter. Firth, `glmer`, and `glmmTMB` hooks are now present as optional script-level methods; GMMAT remains deferred.
- Keep refinement script-based until dependencies and data contracts are stable.
- Output beta, SE, Wald, p-value, runtime, method, and convergence status.

## Later Roadmap

Completed:

- v0.4: experimental ordinal cumulative-binary PSR-GWAS MVP.
- v0.5: stacked proportional-odds ordinal PSR-GWAS with dense ZKZ, grid/profile REML, corrected Wald / pseudo-score calibration, score diagnostics, and large-N validation.

Next:

- Repeated-measure PSR-GWAS with additive genetic + permanent environment + working residual variance.
- Random-regression PSR-GWAS experimental entry point.
- Larger ordinal calibration grid and external real-data validation for manuscript reporting.
