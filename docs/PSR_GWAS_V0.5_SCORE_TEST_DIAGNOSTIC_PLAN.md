# PSR-GWAS v0.5 Score Test Diagnostic Plan

## 0. Landing Decision

For v0.5 delivery, the recommended ordinal stacked-PSR path is:

```
--ordinal-mode stacked-psr \
--ordinal-stacked-engine dense-zkz \
--ordinal-vc-method grid-reml
```

The primary association statistic is the standard `.gwas.tsv` Wald statistic. The previous lambda-GC collapse was caused by an implementation bug in the two-sided Wald p-value calculation (`2*Phi(|z|)` instead of `erfc(|z|/sqrt(2))`), not by a fundamental failure of the stacked-PSR test.

`--ordinal-test all` remains useful for diagnostics. `--ordinal-test robust-score` makes the cluster-robust score the primary `.gwas.tsv` statistic, but it is no longer the recommended default because the corrected Wald / pseudo-score statistic is better calibrated and more powerful in the focused v0.5 calibration.

Focused calibration used for the v0.5 landing decision:

| N | M | Reps | Method | Test | Mean Lambda GC | Power Top1 | Power Top10 | Convergence |
|---|---|------|--------|------|----------------|------------|-------------|-------------|
| 200 | 500 | 5 | dense-zkz-grid-reml | Wald / psr_score | 1.0141 | 60% | 100% | 100% |
| 500 | 500 | 5 | dense-zkz-grid-reml | Wald / psr_score | 0.9352 | 100% | 100% | 100% |
| 200 | 500 | 5 | dense-zkz-grid-reml | cluster_score | 0.9213 | 20% | 100% | 100% |
| 500 | 500 | 5 | dense-zkz-grid-reml | cluster_score | 0.7516 | 100% | 100% | 100% |

Interpretation: the p-value fix resolves the apparent ordinal calibration failure. Wald and pseudo-score are algebraically identical and are now acceptable in focused large-N calibration. Cluster-robust score is retained as an experimental diagnostic / sensitivity statistic.

## 1. Motivation

The earlier large-N calibration showed that all stacked PSR ordinal methods (fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml) produced **lambda GC = 0.000** with the standard Wald test. Follow-up diagnostics identified the root cause as an incorrect two-sided p-value formula in the stacked ordinal Wald output, not as a failure of the PSR approximation.

The Wald test uses: `Wald = (g'Py)² / (g'Pg)`. The pseudo-LMM score test is algebraically equivalent: `Score = (g'Py)² / (g'Pg)`. Therefore the pseudo-LMM score cannot be expected to perform differently.

Alternative test statistics that may differ:

| Test | Formula | Expected Difference |
|------|---------|-------------------|
| Wald (current) | (g'Py)² / (g'Pg) | Current baseline |
| Pseudo-LMM score | (g'Py)² / (g'Pg) | Algebraically identical to Wald |
| Raw working score | (g'W r_raw)² / (g'W g) | Uses raw residual r_raw = z - mu instead of working vector y_tilde |
| Cluster-robust score | (ΣS_i)² / (ΣS_i²) | Aggregates by individual, uses empirical variance |

## 2. CLI

```
--ordinal-test wald|score|robust-score|all   (default: wald)
```

- `wald`: corrected Wald statistic in `.gwas.tsv` (default)
- `score`: pseudo-LMM score in `.gwas.tsv` plus diagnostics (algebraically identical to Wald)
- `robust-score`: cluster-robust score in `.gwas.tsv` plus diagnostics
- `all`: diagnostic mode; `.gwas.tsv` remains Wald and `<out>.ordinal.score.tsv` contains all score variants

Only effective with `--ordinal-mode stacked-psr`.

## 3. Test Statistics

### A. Pseudo-LMM Score

```
U = g' P y_tilde
I = g' P g
score_stat = U² / I = (g'Py)² / (g'Pg) = Wald
```

Expected: **identical to Wald** within numerical precision. If this is not the case, there is a bug.

### B. Raw Working Score

```
r_raw = z - mu              # raw cumulative residual
g_resid = P_wls g           # WLS-projected genotype
U_raw = g_resid' W r_raw
I_raw = g_resid' W g_resid
score_stat_raw = U_raw² / I_raw
```

Uses the raw residual (observed minus fitted cumulative probability) instead of the PSR-adjusted working vector. For fixed-gls, the projection is WLS (`raw_score_projection=wls`); for dense-zkz, it uses the dense V projection (`raw_score_projection=dense_v`).

### C. Cluster-Robust Score

```
For each individual i:
  S_i = Σ_k g_resid_ik * r_raw_ik

U_cluster = Σ_i S_i
Var_cluster = Σ_i S_i²
cluster_score_stat = U_cluster² / Var_cluster
```

Aggregates the score contribution by original individual (across ordinal thresholds). The empirical variance `Σ S_i²` is robust to within-individual correlation among the stacked rows.

## 4. Output

When `--ordinal-test` is not `wald`, the file `<out>.ordinal.score.tsv` is written with columns:

| Column | Description |
|--------|-------------|
| wald_stat, wald_p, wald_log10p | Standard Wald test |
| psr_score_stat, psr_score_p, psr_score_log10p | Pseudo-LMM score (identical to Wald) |
| raw_score_stat, raw_score_p, raw_score_log10p | Raw working residual score |
| cluster_score_stat, cluster_score_p, cluster_score_log10p | Cluster-robust score |
| wald_psr_score_equivalent | TRUE/FALSE (TRUE confirms equivalence) |
| raw_score_projection | "wls" or "dense_v" |
| n_clusters | Number of individuals aggregated |

## 5. Scoring Criteria

Calibration is interpreted as:

| Lambda GC Range | Interpretation |
|----------------|---------------|
| 0.8 - 1.2 | Well-calibrated |
| 0.5 - 0.8 | Moderately conservative |
| < 0.5 | Persistent conservative bias |
| > 1.3 | Inflated / anti-conservative |

The production default remains Wald / pseudo-score when its lambda GC is in the 0.8-1.2 range. Cluster-score is a sensitivity statistic and can be used with `--ordinal-test robust-score` if future datasets show Wald inflation.

If Wald / pseudo-score again show lambda GC near zero, first verify the p-value range (`0 <= p <= 1`) and the `wald_psr_score_equivalent` diagnostic before revisiting the model approximation.

## 6. Files

- Source: `src/config.h`, `src/main.cpp`, `src/models/ordinal_model.h`, `src/models/ordinal_model.cpp`
- Validation: `tests/run_psr_ordinal_score_test_validation.sh`
- Calibration: `tests/run_psr_v05_large_n_calibration.py` (via `PSR_CALIB_TESTS=all`)
- Smoke: `tests/run_psr_v05_score_calibration_smoke.sh`
- Output: `tests/results/psr_v05_large_n_calibration/large_n_score_diagnostic_report.md`

## References

- Breslow & Clayton (1993). Approximate inference in generalized linear mixed models. JASA.
- Liang & Zeger (1986). Longitudinal data analysis using generalized linear models. Biometrika.
- PSR-GWAS large-N calibration plan: `docs/PSR_GWAS_V0.5_LARGE_N_CALIBRATION_PLAN.md`
