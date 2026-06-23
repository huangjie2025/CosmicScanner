# PSR-GWAS v0.5 Landing Summary

Date: 2026-06-23
Branch: `psr-gwas-v05-score-test-calibration`

## Final Recommendation

Ordinal PSR-GWAS v0.5 is ready to land with the stacked proportional-odds PSR path:

```bash
--ordinal-mode stacked-psr \
--ordinal-stacked-engine dense-zkz \
--ordinal-vc-method grid-reml
```

The primary result is the corrected Wald / pseudo-score statistic in `.gwas.tsv`.

Use `--ordinal-test all` only for diagnostic output in `.ordinal.score.tsv`. Use `--ordinal-test robust-score` only for cluster-robust sensitivity analysis.

## Root Cause Fixed

The apparent stacked ordinal lambda-GC collapse was caused by an implementation bug in the two-sided Wald p-value conversion. The code used a lower-tail expression equivalent to `2*Phi(|z|)`, which can produce invalid p-values greater than 1 and made lambda GC collapse toward zero.

The fix uses the correct upper-tail two-sided normal p-value:

```text
p = erfc(|z| / sqrt(2))
```

After this fix, Wald and pseudo-score are algebraically equivalent and calibrated in focused large-N validation.

## Final Calibration

Focused calibration command:

```bash
PSR_CALIB_REPS=5 \
PSR_CALIB_N_LIST=200,500 \
PSR_CALIB_M=500 \
PSR_CALIB_METHODS=dense-zkz-grid-reml \
PSR_CALIB_TESTS=all \
python3 tests/run_psr_v05_large_n_calibration.py
```

Results:

| N | M | Reps | Method | Lambda GC | Top1 | Top10 | Convergence | Verdict |
|---|---|------|--------|-----------|------|-------|-------------|---------|
| 200 | 500 | 5 | dense-zkz-grid-reml | 1.0141 | 60% | 100% | 100% | acceptable |
| 500 | 500 | 5 | dense-zkz-grid-reml | 0.9352 | 100% | 100% | 100% | acceptable |

## Validation Passed

```bash
make -C build cosmicscanner -j4
git diff --check
bash tests/run_psr_ordinal_score_test_validation.sh
bash tests/run_psr_v05_score_calibration_smoke.sh
bash tests/run_psr_v05_large_n_calibration_smoke.sh
bash tests/run_psr_v05_validation_benchmark.sh
bash tests/run_psr_full_parity_validation.sh
```

Observed pass criteria:

- Score validation: 30/30 passed
- Score calibration smoke: 8/8 passed
- Large-N calibration smoke: 11/11 passed
- v0.5 validation benchmark: 8/8 passed
- Full parity validation: passed

## Commits

```text
394e89d Fix ordinal PSR score calibration
857605f Document ordinal PSR v0.5 landing path
```

This landing summary is the final handoff note for the v0.5 ordinal PSR-GWAS implementation.

## Remaining Work

No blocker remains for v0.5 delivery.

Recommended post-landing work:

- Larger calibration grid for manuscript reporting (`reps >= 20`, N up to 1000+).
- External real-data validation for ordinal phenotypes.
- Optional exact ordinal mixed-model comparison when POLMM / ordinal mixed-model baselines are available.
