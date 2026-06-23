# PSR-GWAS Roadmap And Acceptance Criteria

Goal: make binary/binomial/Poisson PSR-GWAS a statistically auditable, reproducible CosmicScanner module before expanding to ordinal or repeated-measure PSR-GLMM.

## Current Scope

Implemented trait types:

- Binary: `--threshold-model --link probit|logit`
- Poisson: `--count-model --family poisson`, optional `--offset`
- Binomial: `--count-model --family binomial --binomial-pheno FID IID EVENT TRIAL`

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

- v0.4: ordinal PSR-GWAS, validated against POLMM or another ordinal baseline.
- v0.5: stacked ordinal PSR-GWAS Stage 1 (fixed-VC weighted GLS first pass). Completed: stacked data design, threshold encoding, weighted IRLS, per-SNP Wald scan, covariate support, no-NaN validation, fallback logic, documentation.
- v0.5+: stacked ordinal PSR-GWAS Stage 2A (dense ZKZ' BLUP first pass). **Completed**: `--ordinal-stacked-engine dense-zkz`, full `K_stack = Z * K * Z'`, `V = sigma2_g * K_stack + diag(1/w)` GLS update, full `P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}` Wald scan, `.ordinal.blup.tsv` individual BLUP output, `dense_zkz_blup_first_pass` implementation label, `N_stack <= 2000` safeguard, standalone validation wrapper, parity matrix update.
- v0.5++: stacked ordinal PSR-GWAS Stage 2B (grid/profile REML for sigma2_g). **Completed**: `--ordinal-vc-method fixed|grid-reml`, grid REML sigma2_g estimation within each PQL iteration, `.ordinal.vc.tsv`, new output fields, sigma2_g=0 boundary, BLUP with updated sigma2_g, validation (15/15 pass).
- v0.5+++: stacked ordinal PSR-GWAS Stage 2C (AI-REML or faster optimization for sigma2_g, large-sample support). **Planned**.
- v0.6: repeated-measure PSR-GWAS with additive genetic + permanent environment + working residual variance.
- v0.7: random-regression PSR-GWAS experimental entry point.
# v0.4 Parity-Ready Sprint

v0.4 adds a parity-ready validation surface for normal, binary, binomial, Poisson, and ordinal traits.

- Normal traits are handled as the Gaussian identity-link special case using existing LMM GWAS and validated by `tests/run_psr_normal_special_case.sh`.
- Binary, binomial, and Poisson remain v0.3 stable PSR/PQL implementations.
- Ordinal traits are implemented as `--ordinal-model --ordinal-pheno` using an experimental cumulative-binary approximation. This is explicitly not full stacked proportional-odds PSR-GLMM parity.
- Official PSR-GLMM/R availability is tracked by `tests/check_psr_glmm_r_availability.sh`; unavailability is not a CI blocker.
- Alternative benchmarks are generated by `tests/run_psr_accuracy_speed_benchmark.py` and optional R package availability is recorded without failing CI.

# v0.4.1 Validation Infrastructure

v0.4.1 extends v0.4 parity-ready into a reportable validation infrastructure:

- Formal simulation (`tests/run_psr_formal_simulation.py`) expanded with richer scenarios: binary null_rare (5% prevalence), link_logit variants, poisson null_offset, binomial null_rare_prop/null_balanced_prop, plus `link`, `ordinal_mode`, `has_offset`, `error_message` output fields and `psr_formal_simulation_report.md` generation.
- Accuracy/speed benchmark (`tests/run_psr_accuracy_speed_benchmark.py`) now generates `psr_accuracy_speed_report.md` with trait coverage, accuracy, runtime, and method availability summaries.
- `--ordinal-mode` and `--ordinal-allow-fallback` CLI flags added to `config.h` and `main.cpp`. Current values: `experimental-cumulative-binary` (default), `stacked-psr` (design-only, raises clear error). See `docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md`.
- Full ordinal PSR design document (`docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md`) documents stacked proportional-odds PSR, C-1 thresholds, shared beta, stacked pseudo-response, W, K_stack, V_stack, and 5-stage engineering plan. Status: `design_only`.
- Full ordinal design check (`tests/run_psr_full_ordinal_design_check.sh`) confirms design document presence and reports `full_ordinal_psr_status=design_only`.
- Full parity validation wrapper (`tests/run_psr_full_parity_validation.sh`) includes all new steps: formal simulation, accuracy/speed report, ordinal design check.
- Commit tag `psr-gwas-v0.4-parity-ready` frozen at v0.4 stable point.
- `.gitignore` added to ignore `build/`, `tests/__pycache__/`, `tests/results/`.

Repeated-measure PSR-GLMM, random-regression PSR-GWAS, exact GLMM inference, and stable negative-binomial reporting remain future roadmap items.

## Recent Progress (2026-06-23)

### Item 1: PSR-GWAS v0.5 Ordinal Methods Full Validation (2026-06-23)
- Extended formal simulation (`run_psr_formal_simulation.py`) to support `PSR_FORMAL_ORDINAL_METHODS=all`, running all 4 ordinal methods per rep
- Created `run_psr_v05_ordinal_benchmark.py` standalone benchmark runner with per-method summary outputs and markdown report
- Created `run_psr_v05_validation_benchmark.sh` wrapper with integrity checks
- Created `compare_psr_benchmarks.py` for comparing benchmark results across commits
- Created `docs/PSR_GWAS_V0.5_VALIDATION_REPORT.md` validation report
- Verified: all 4 methods produce output (cumulative-binary, fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml)
- Verified: formal simulation with PSR_FORMAL_ORDINAL_METHODS=all generates per-method report sections and comparison TSVs
- Verified: backward compatibility with PSR_FORMAL_ORDINAL_ENGINE / PSR_FORMAL_ORDINAL_VC_METHOD env vars
- Verified: R baseline (MASS::polr, ordinal::clm) available; both produce valid 30/30 SNP results
- Tagged: `psr-gwas-v0.5-r-baseline-medium-validation`
- Remaining: real phenotype validation, paper-level benchmark runs (REPS=20+, N=200+)

### Item 2: Stacked PSR Large-N Calibration (2026-06-23)
- Known issue: all stacked PSR methods show lambda GC ≈ 0 at N=200, M=500 — conservative test statistics, weak power
- Created `tests/run_psr_v05_large_n_calibration.py` for systematic N-sweep calibration (N=100→1000)
- Created `tests/run_psr_v05_large_n_calibration_smoke.sh` validation wrapper
- Created `docs/PSR_GWAS_V0.5_LARGE_N_CALIBRATION_PLAN.md` documenting possible causes and diagnosis
- Diagnostic verdict rules: `persistent_conservative_bias`, `sample_size_dependent_pql_conservativeness`, `inflated`, `acceptable_or_needs_more_reps`
- N_stack guard for dense-zkz (n * (levels-1) > 2000 → `unsupported_dense_limit`, not a crash)
- Unknown: whether conservativeness resolves at N ≥ 500
- Next: run calibration (smoke: PSR_CALIB_REPS=1 N_LIST=80,120; medium: REPS=3 N_LIST=200,500,1000)
- Note: Stage 2C AI-REML/Brent is an optimization improvement, not guaranteed to fix lambda GC issue
