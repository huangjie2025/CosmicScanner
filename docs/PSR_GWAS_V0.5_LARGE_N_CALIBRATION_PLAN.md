# PSR-GWAS v0.5 Large-N Calibration Plan

## 1. Why Large-N Calibration?

The N=200 M=500 benchmark showed that all stacked PSR ordinal methods (fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml) produce **lambda GC ≈ 0**, indicating conservative test statistics (p-values systematically too large). Correspondingly, power is weak — top10 hit rate for the causal SNP is near 0%.

This calibration study systematically tests whether the conservativeness resolves at larger sample sizes, or whether it is a structural issue with the PSR/PQL approximation for ordinal traits.

## 2. Observations from N=200 Benchmark

| Method | Lambda GC | Top10% | Window% |
|--------|-----------|--------|---------|
| cumulative-binary (baseline) | ~2.0 (inflated) | 50% | 75% |
| fixed-gls (stacked PSR) | ~0.0 | 3% | 20% |
| dense-zkz-fixed (stacked PSR) | ~0.0 | 0% | 15% |
| dense-zkz-grid-reml (stacked PSR) | ~0.0 | 0% | 10% |

**Key**: cumulative-binary (per-SNP cumulative logit) has inflated lambda GC at N=200, while stacked PSR methods have lambda GC ≈ 0. The two approaches are on opposite sides of calibration.

## 3. Possible Causes

1. **PQL approximation conservativeness**: PQL is known to be conservative for binary/ordinal data, especially with small clusters and few levels. The working weights and pseudo-response may not adequately capture the ordinal information.

2. **Wald SE overestimation**: The Wald test on the stacked PSR null model may produce standard errors that are too large, shrinking chi-square statistics toward zero.

3. **Latent scale vs. observed scale**: PSR operates on the latent scale through the ordinal link function; p-values from the latent-scale scan may not be directly comparable to the cumulative-binary approach.

4. **V matrix / working weights absorbing signal**: The pseudo-response variance-covariance structure may absorb part of the genetic signal into residual variation.

5. **Small-sample category imbalance**: With 3 levels and N=200, some categories may have few observations, making the PQL approximation less stable.

## 4. Calibration Design

We test multiple N values with the same M (marker count) and ordinal levels:

| N | M | Levels | Reps | Methods | Scenarios |
|---|----|--------|------|---------|-----------|
| 100 | 200 | 3 | 5 | all 4 | null, power |
| 200 | 200 | 3 | 5 | all 4 | null, power |
| 500 | 200 | 3 | 5 | all 4 | null, power |
| 1000 | 200 | 3 | 5 | all 4 | null, power |

If lambda GC increases from ~0 toward 1 as N grows from 100 to 1000, the issue is **sample-size-dependent PQL conservativeness** (can be managed by requiring larger samples).

If lambda GC remains < 0.5 even at N=1000, the issue is **persistent conservative bias** (needs algorithmic correction).

## 5. Diagnostic Verdict Rules

| Condition | Verdict |
|-----------|---------|
| lambda_gc < 0.5 at max(N) | `persistent_conservative_bias` |
| lambda_gc increases toward 1 as N grows | `sample_size_dependent_pql_conservativeness` |
| lambda_gc > 1.3 at any N | `inflated` |
| lambda_gc in [0.8, 1.2] at max(N) | `acceptable_or_needs_more_reps` |

Note: This verdict is a diagnostic engineering assessment, **not** a paper-level conclusion.

## 6. Method Comparison

- **fixed-gls**: Fastest, weakest BLUP support. If lambda GC improves with N, it's usable for screening.
- **dense-zkz-fixed**: BLUP available, sigma2_g fixed at 0.5. If conservativeness persists, the fixed VC may be the wrong value.
- **dense-zkz-grid-reml**: Full REML VC estimation. Most informative for understanding whether VC estimation affects calibration.
- **cumulative-binary**: External baseline. Expected to have proper or slightly inflated lambda GC at moderate N.

Note: dense-zkz methods have `N_stack = n * (levels - 1) ≤ 2000` limitation. For levels=3 and N=1000, N_stack=2000, which is the limit. For levels=4, N_stack=3000 > 2000, so dense-zkz will be recorded as `unsupported_dense_limit`.

## 7. Next Steps if Conservativeness Persists

If lambda GC remains < 0.5 at N ≥ 500:

1. **PQL correction**: Investigate Breslow & Clayton-type corrections for the PQL approximation in ordinal models.
2. **Robust/sandwich SE**: Replace model-based SE with robust (sandwich/Huber-White) SE in the Wald test.
3. **Score test**: Replace Wald test with score (Lagrange multiplier) test, which may be less sensitive to variance misspecification.
4. **Firth-type correction**: Apply Firth's penalized likelihood to reduce small-sample bias.
5. **Stage 2C AI-REML / Brent**: These are optimization improvements for the null model, not guaranteed to fix the lambda GC issue.

If the issue resolves at N ≥ 500:
- Document recommended minimum sample size for stacked PSR ordinal GWAS.
- Consider adding an automatic small-sample warning when N < 500.

## 8. Branch and Scripts

- Branch: `psr-gwas-v05-large-n-calibration`
- Runner: `tests/run_psr_v05_large_n_calibration.py`
- Smoke test: `tests/run_psr_v05_large_n_calibration_smoke.sh`
- Output dir: `tests/results/psr_v05_large_n_calibration/`

## References

- Breslow & Clayton (1993). Approximate inference in generalized linear mixed models. JASA.
- Lin & Breslow (1996). Bias correction in generalized linear mixed models with multiple components of dispersion. JASA.
- PSR-GWAS v0.5 validation report: `docs/PSR_GWAS_V0.5_VALIDATION_REPORT.md`
