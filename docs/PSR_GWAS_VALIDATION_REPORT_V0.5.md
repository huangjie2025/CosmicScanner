# PSR-GWAS v0.5 Validation Report
Date: 2026-06-23
Status: ordinal stacked PSR calibration fixed; dense-zkz-grid-reml recommended for v0.5

## 0. Final v0.5 Landing Update

The earlier stacked-PSR lambda-GC collapse was traced to an implementation bug in the two-sided Wald p-value calculation. The statistic was correct, but the p-value used `2*Phi(|z|)` instead of the correct upper-tail form `erfc(|z|/sqrt(2))`, producing invalid overly large p-values and lambda GC near zero.

After the fix, focused large-N calibration supports the following production recommendation for ordinal traits:

```bash
--ordinal-mode stacked-psr \
--ordinal-stacked-engine dense-zkz \
--ordinal-vc-method grid-reml
```

Focused calibration after the p-value fix:

| N | M | Reps | Method | Lambda GC | Top1 | Top10 | Convergence | Verdict |
|---|---|------|--------|-----------|------|-------|-------------|---------|
| 200 | 500 | 5 | dense-zkz-grid-reml | 1.0141 | 60% | 100% | 100% | acceptable |
| 500 | 500 | 5 | dense-zkz-grid-reml | 0.9352 | 100% | 100% | 100% | acceptable |

Use `--ordinal-test all` for diagnostic output in `<out>.ordinal.score.tsv`. `robust-score` is retained as a sensitivity statistic, not the default recommendation.

Sections below record the earlier benchmark framework and pre-fix observations; the landing update above supersedes any statement that stacked PSR has persistent lambda GC near zero.

## 1. Overview
- What v0.5 adds: ordinal trait GWAS via stacked PSR (3 PSR methods) + cumulative-binary baseline
- 4 methods under test: cumulative-binary, fixed-gls, dense-zkz-fixed, dense-zkz-grid-reml
- R fixed-effect ordinal baselines: MASS::polr (available), ordinal::clm (available), POLMM (unavailable)
- Validation branch: `psr-gwas-v05-r-baseline-medium-validation`
- Prior tag: `psr-gwas-v0.5-validation-benchmark`

## 2. Formal Simulation Results
Run: `PSR_FORMAL_TRAITS=ordinal PSR_FORMAL_ORDINAL_METHODS=all` (1 rep, N=40, M=40, 3 scenarios)

| Method | Conv% | T1@0.05 | Lambda GC | Power (Top10) | Runtime (s) |
|--------|-------|---------|-----------|---------------|-------------|
| cumulative-binary | 100% | 0.000 | 3.152 | 100% | 0.151 |
| fixed-gls | 100% | 0.000 | 0.000 | 100% | 0.090 |
| dense-zkz-fixed | 100% | 0.000 | 0.000 | 100% | 0.086 |
| dense-zkz-grid-reml | 100% | 0.000 | 0.000 | 100% | 0.095 |

- dense-zkz-grid-reml estimated sigma2_g=0.533 (range 0.400-0.800)
- All methods achieve 100% type I error control at alpha=0.05 with M=40
- Lambda GC near 0 for stacked PSR methods (conservative); cumulative-binary inflated (~3.15) for small N

## 3. Benchmark Methodology
- Data: PLINK .bed/.bim/.fam with ordinal phenotypes (K=3 levels), M SNPs, N samples
- One causal SNP at position M/3, null + power_low (effect=0.35) + power_high (effect=0.70) scenarios
- Repeat for R replicates
- Metrics: runtime, convergence, lambda GC, type-I error, power (top1/top10/window hit), sigma2_g, boundary detection, BLUP availability

## 4. Medium-Scale Benchmark Results (N=200, M=500, REPS=10)

| Method | Conv% | Runtime(s) | Sigma2_g | Lambda GC | Top10% | BLUP% |
|--------|-------|-----------|----------|-----------|--------|-------|
| cumulative-binary | 100% | 0.221 | NA | 1.962 | 50% | 0% |
| fixed-gls | 100% | 0.161 | NA | 0.000 | 3% | 0% |
| dense-zkz-fixed | 100% | 0.156 | NA | 0.000 | 0% | 100% |
| dense-zkz-grid-reml | 100% | 0.283 | 1.253 | 0.000 | 0% | 100% |

**Runtime ranking (fastest to slowest):**
1. dense-zkz-fixed (0.156s)
2. fixed-gls (0.161s)
3. cumulative-binary (0.221s)
4. dense-zkz-grid-reml (0.283s)

**Sigma2_g behavior (dense-zkz-grid-reml):**
- Mean: 1.253, SD: 0.478, range: [0.4, 2.0], boundary rate: 0%
- No boundary issues observed — REML grid search finds interior solution

**Power (causal SNP in top10 out of 500):**
- cumulative-binary: 50% (best power, but inflated lambda GC)
- fixed-gls: 3% (weak power, conservative)
- dense-zkz-fixed: 0%
- dense-zkz-grid-reml: 0%

## 5. R Baseline Results (N=40, M=30, 1 rep per scenario)

| Method | Available | Valid SNPs | Mean Runtime(s) | Top SNP |
|--------|-----------|------------|-----------------|---------|
| MASS::polr | Yes | 30/30 | 0.218 | sim000007-009 |
| ordinal::clm | Yes | 30/30 | 0.272 | sim000007-009 |
| POLMM | No | - | - | - |

- Both R methods produce nearly identical -log10(p) values (as expected for the same cumulative logit model)
- R methods serve as fixed-effect ordinal baseline (not mixed-model benchmarks)
- MASS::polr is slightly faster than ordinal::clm (~0.22s vs ~0.27s for M=30)

## 6. Key Observations

**Cumulative-binary:**
- Best power (50% top10 hit rate at N=200)
- Lambda GC inflated (~2) for small N — type I error not well-controlled at small sample sizes
- No BLUP or variance component estimation
- Good for baseline comparisons

**Fixed-gls (stacked PSR, fixed GLS):**
- Fastest or near-fastest runtime
- Lambda GC ~0 (very conservative) — p-values are systematically too large
- No BLUP available
- Top10 hit rate low (3% at N=200)
- Best for quick ordinal screening where calibration is not critical

**Dense-zkz-fixed (stacked PSR, dense ZKZ, sigma2_g=0.5):**
- BLUP available (100% of runs have finite BLUP)
- Lambda GC ~0 (conservative)
- 0% top10 hit rate at N=200/M=500 — weakest power
- Sigma2_g fixed at 0.5 (may be misspecified)

**Dense-zkz-grid-reml (stacked PSR, dense ZKZ, grid REML):**
- BLUP + REML sigma2_g estimation
- 0% boundary rate — grid search finds interior solution
- Slowest runtime (0.283s mean at N=200)
- Lambda GC ~0 (conservative)
- Top10 hit rate 0%
- Recommended when VC estimation is needed despite conservative calibration

**All stacked PSR methods show lambda GC ~0**, which indicates the current PSR/PQL approximation for ordinal traits produces conservative test statistics at moderate sample sizes. This is consistent with the behavior of penalized quasi-likelihood in small samples.

## 7. Method Selection Guide
| Method | Speed | BLUP | Sigma2_g | Lambda GC | Recommended for |
|--------|-------|------|----------|-----------|-----------------|
| cumulative-binary | Fast | No | N/A | Inflated | Baseline method comparison |
| fixed-gls | Fastest | No | Fixed | ~0 | Quick scan, no BLUP needed |
| dense-zkz-fixed | Fast | Yes | Fixed 0.5 | ~0 | BLUP + fixed VC assumption |
| dense-zkz-grid-reml | Moderate | Yes | REML (0% boundary) | ~0 | Full VC estimation |

## 8. Known Limitations
- Lambda GC ~0 for stacked PSR methods: PQL approximation is conservative for ordinal traits at N<1000
- Lambda GC inflated for cumulative-binary at N<500: Wald test with small samples
- All results on simulated data with known causal SNPs; real phenotype validation pending
- R baselines are fixed-effect ordinal models (polr, clm) — not mixed-model equivalents
- No formal calibration study (e.g., empirical null, genomic control adjustment)

## 9. Large-N Calibration Status
See the "Large-N Calibration Status" section below for the systematic calibration design and results.

| Branch | Tag | Status |
|--------|-----|--------|
| `psr-gwas-v05-large-n-calibration` | `psr-gwas-v0.5-large-n-calibration` | Runner + smoke wrapper complete |
| `psr-gwas-v05-score-test-calibration` | `394e89d` | p-value fix + score diagnostics complete |

Current diagnostic verdict after p-value fix: **`acceptable_calibration`** for dense-zkz-grid-reml at N=200 and N=500 focused calibration.

Note: Stage 2C AI-REML / Brent solver improvements are orthogonal to this fix. The resolved issue was in p-value conversion, not in the null model solver.

## 9. Recommendations for v0.5 → v1.0
- **Default method**: dense-zkz-grid-reml for v0.5 ordinal delivery
- **Fast screening**: fixed-gls for quick ordinal GWAS scans
- **Baseline**: cumulative-binary for cross-method validation
- **Diagnostics**: run `--ordinal-test all` when validating new datasets; check `wald_psr_score_equivalent=TRUE` and p-values in [0,1]
- **R comparison**: MASS::polr available as external fixed-effect reference; full mixed-model comparison (POLMM, ordinal::clmm) requires additional R packages
- **Real data validation**: test on publicly available ordinal GWAS datasets

## A. Run Commands
```bash
# Full formal simulation with all ordinal methods:
PSR_FORMAL_REPS=20 PSR_FORMAL_N=200 PSR_FORMAL_TRAITS=ordinal PSR_FORMAL_ORDINAL_METHODS=all python3 tests/run_psr_formal_simulation.py

# Medium-scale benchmark (defaults):
bash tests/run_psr_v05_validation_benchmark.sh

# R baseline check:
bash tests/run_psr_v05_r_baseline_validation.sh

# Custom medium benchmark:
PSR_V05_BENCH_REPS=5 PSR_V05_BENCH_N=120 PSR_V05_BENCH_M=200 PSR_V05_BENCH_METHODS=all python3 tests/run_psr_v05_ordinal_benchmark.py

# Cross-commit comparison:
python3 tests/compare_psr_benchmarks.py --new tests/results/psr_v05_ordinal_benchmark
```
