# CosmicScanner Repeatability GWAS Methods and Results Draft

## Manuscript Positioning

CosmicScanner implements a scalable exact repeated-measure GWAS framework for related individuals. The method targets longitudinal or repeated-record phenotypes commonly observed in livestock, plant, and biomedical cohorts, where each genotyped subject can have multiple phenotypic records and record-level covariates such as batch, year, parity, or environment.

The main contribution is not a C++ reimplementation of RepeatABEL. The contribution is a structured record-level mixed model solver that preserves repeated-measure modeling while avoiding dense record-level covariance factorization through subject-level Woodbury algebra and incidence-aware block SNP scanning.

## Statistical Model

For record `j` of subject `i`, the repeatability GWAS null model is

```text
y = Xb + Zu + Zp + e
u ~ N(0, sigma_g^2 G)
p ~ N(0, sigma_pe^2 I)
e ~ N(0, sigma_e^2 W^-1)
```

where `G` is the subject-level genomic relationship matrix, `Z` maps records to subjects, and `W` is an optional diagonal matrix of record-specific residual weights. The marginal covariance is

```text
V = ZKZ' + R
K = sigma_g^2 G + sigma_pe^2 I
R = sigma_e^2 W^-1
```

CosmicScanner uses the Woodbury identity to avoid explicit construction and factorization of the dense record-level `V`:

```text
V^-1 = R^-1 - R^-1 Z (K^-1 + Z'R^-1Z)^-1 Z'R^-1
```

The core matrix requiring dense factorization is therefore subject-level:

```text
M = K^-1 + Z'R^-1Z
```

This changes the main covariance operation from record-level dense algebra to subject-level dense algebra with record-level aggregation through the incidence matrix.

## SNP Scan

For each SNP, CosmicScanner keeps genotype dosage at the subject level and never expands the SNP vector to a full record-level matrix block. SNPs are processed in blocks. The scan reuses the null-model Woodbury cache, computes projected score/Wald quantities, and writes association statistics for PLINK, PGEN, and BGEN inputs through a shared block scan path.

This design is incidence-aware: genotype blocks have dimension approximately `n_subjects x block_size`, while the repeated records contribute through subject-level sums and weighted residual projections.

## Variance Components

The default production setting is:

```bash
--repeatability-model --repeatability-mode exact --vce-mode auto
```

`--vce-mode auto` uses subject-level Woodbury VCE and skips full-record exact refinement by default. Full record-level exact refinement remains available for small-data validation:

```bash
--vce-mode full-record-exact
```

This policy is important for scalability. RepeatABEL constructs a dense record-level `V` and performs `eigen(V)`, whereas CosmicScanner keeps the exact repeated-measure covariance structure but solves through subject-level Woodbury matrices.

## Implementation Features

Implemented features include:

- Record-level repeatability LMM with additive genetic, permanent environmental, and residual variance components.
- Record-level quantitative and categorical covariates keyed by `FID IID TIME`.
- Optional record-specific residual weights through `--repeatability-weight`.
- Exact Woodbury scan for PLINK, PGEN, and BGEN.
- Global and LOCO scans.
- PLINK/PGEN/BGEN consistency tests for global, weighted global, and LOCO analyses.
- RepeatABEL-compatible benchmark and accuracy comparison scripts.
- Stage timing fields in `.repeatability.null.txt`: `time_grm_sec`, `time_vce_sec`, `time_scan_vbuild_sec`, and `time_snp_scan_sec`.

## Accuracy Against RepeatABEL

RepeatABEL 1.1 was installed from CRAN archive into a temporary R library and patched only for compatibility with R 4.3/modern compilers. Accuracy comparisons used identical simulated PLINK data, identical long-format phenotypes, and the same record-level batch covariate. RepeatABEL was run with `rGLS(y ~ batch, ...)`; CosmicScanner was run with exact repeatability GWAS.

### Accuracy Scenario 1

```text
N subjects = 80
records per subject = 3
records = 240
SNPs = 300
```

| Metric | Value |
|---|---:|
| SNP overlap | 300 / 300 |
| Signed beta correlation | -0.9932 |
| Best-orientation beta correlation | 0.9932 |
| Absolute beta correlation | 0.9846 |
| SE correlation | 0.9990 |
| log10P correlation | 0.9851 |
| Top10 overlap | 9 / 10 |
| Top20 overlap | 17 / 20 |
| Median absolute log10P difference | 0.0427 |
| Max absolute log10P difference | 0.2418 |

The negative signed beta correlation reflects opposite allele dosage orientation between the two tools. After orientation correction, effect estimates and association ranks are highly concordant.

### Accuracy Scenario 2

```text
N subjects = 80
records per subject = 5
records = 400
SNPs = 300
```

| Metric | Value |
|---|---:|
| SNP overlap | 300 / 300 |
| Signed beta correlation | -0.9996 |
| Best-orientation beta correlation | 0.9996 |
| Absolute beta correlation | 0.9992 |
| SE correlation | 0.99996 |
| log10P correlation | 0.9991 |
| Top10 overlap | 10 / 10 |
| Top20 overlap | 19 / 20 |
| Median absolute log10P difference | 0.0109 |
| Max absolute log10P difference | 0.1250 |

These results indicate that CosmicScanner Exact/Woodbury reproduces RepeatABEL association results with near-perfect effect-size, standard-error, and significance-rank concordance under the same model specification.

## Performance Results

Benchmarks were run with one BLAS/OpenMP thread unless otherwise stated. SNP count was fixed at 1000 for the expanded benchmark. Runtime includes data loading, GRM construction, VCE, and association scan.

| Subjects | Records/Subject | Records | SNPs | Method | Runtime (s) | Max RSS (MB) |
|---:|---:|---:|---:|---|---:|---:|
| 100 | 2 | 200 | 1000 | CS Exact/Woodbury | 0.290 | 29.04 |
| 100 | 2 | 200 | 1000 | CS Fast/Mean Approx. | 0.300 | 27.04 |
| 100 | 2 | 200 | 1000 | Subject Mean LMM | 0.080 | 26.05 |
| 100 | 2 | 200 | 1000 | RepeatABEL rGLS | 5.340 | 238.93 |
| 100 | 5 | 500 | 1000 | CS Exact/Woodbury | 0.310 | 29.03 |
| 100 | 5 | 500 | 1000 | CS Fast/Mean Approx. | 0.300 | 29.03 |
| 100 | 5 | 500 | 1000 | Subject Mean LMM | 0.080 | 26.06 |
| 100 | 5 | 500 | 1000 | RepeatABEL rGLS | 5.920 | 254.28 |
| 200 | 2 | 400 | 1000 | CS Exact/Woodbury | 1.060 | 35.80 |
| 200 | 2 | 400 | 1000 | CS Fast/Mean Approx. | 1.070 | 33.82 |
| 200 | 2 | 400 | 1000 | Subject Mean LMM | 0.100 | 30.82 |
| 200 | 2 | 400 | 1000 | RepeatABEL rGLS | 12.530 | 272.93 |
| 200 | 5 | 1000 | 1000 | CS Exact/Woodbury | 0.930 | 35.83 |
| 200 | 5 | 1000 | 1000 | CS Fast/Mean Approx. | 0.880 | 33.82 |
| 200 | 5 | 1000 | 1000 | Subject Mean LMM | 0.100 | 30.81 |
| 200 | 5 | 1000 | 1000 | RepeatABEL rGLS | 10.280 | 316.89 |
| 500 | 2 | 1000 | 1000 | CS Exact/Woodbury | 11.740 | 65.20 |
| 500 | 2 | 1000 | 1000 | CS Fast/Mean Approx. | 13.110 | 54.20 |
| 500 | 2 | 1000 | 1000 | Subject Mean LMM | 0.230 | 46.46 |
| 500 | 5 | 2500 | 1000 | CS Exact/Woodbury | 12.150 | 65.18 |
| 500 | 5 | 2500 | 1000 | CS Fast/Mean Approx. | 11.570 | 53.19 |
| 500 | 5 | 2500 | 1000 | Subject Mean LMM | 0.190 | 46.46 |

The main performance observation is that CosmicScanner Exact/Woodbury is largely insensitive to increasing repeated records from 2 to 5 per subject at fixed subject count, while RepeatABEL memory increases with record count because it constructs record-level covariance matrices. For `N=100-200`, CosmicScanner was approximately 10-20 times faster than RepeatABEL and used roughly one-tenth of the memory.

## Expanded Benchmark Round 2

Because the project filesystem `/disk228` was full during manuscript benchmarking, all second-round benchmark data and results were written to `/tmp/opencode` using the script-level overrides `BENCH_DATA_ROOT` and `BENCH_RESULT_ROOT`. This avoids writing large simulated genotype matrices and RepeatABEL intermediate files to the source tree.

The second-round CosmicScanner benchmark evaluated a larger subject count and genotype dimension:

```text
subjects = 1000
SNPs = 10000
records per subject = 2, 5, 10
methods = CS Exact/Woodbury, subject-mean LMM
```

| Subjects | Records/Subject | Records | SNPs | Method | Status | Runtime (s) | Max RSS (MB) |
|---:|---:|---:|---:|---|---|---:|---:|
| 1000 | 2 | 2000 | 10000 | CS Exact/Woodbury | ok | 98.150 | 154.46 |
| 1000 | 2 | 2000 | 10000 | Subject Mean LMM | ok | 1.780 | 109.23 |
| 1000 | 5 | 5000 | 10000 | CS Exact/Woodbury | ok | 92.470 | 157.17 |
| 1000 | 5 | 5000 | 10000 | Subject Mean LMM | ok | 1.760 | 109.24 |
| 1000 | 10 | 10000 | 10000 | CS Exact/Woodbury | ok | 94.890 | 157.48 |
| 1000 | 10 | 10000 | 10000 | Subject Mean LMM | ok | 1.740 | 108.24 |
| 2000 | 5 | 10000 | 10000 | CS Exact/Woodbury | ok | 728.820 | 451.04 |
| 2000 | 5 | 10000 | 10000 | Subject Mean LMM | ok | 6.710 | 349.82 |

The key scaling result is that increasing repeated records from 2 to 10 per subject did not materially increase CosmicScanner Exact/Woodbury runtime or memory. This supports the intended model design: repeated records are handled mainly through incidence-aware aggregation and Woodbury subject-level matrices rather than dense record-level matrix factorization. Increasing subject count from 1000 to 2000 substantially increased runtime, consistent with the expected dense subject-level GRM/K/M algebra becoming the main bottleneck.

RepeatABEL was benchmarked to document scaling limitations:

| Subjects | Records/Subject | Records | SNPs | Method | Status | Runtime (s) | Max RSS (MB) |
|---:|---:|---:|---:|---|---|---:|---:|
| 500 | 2 | 1000 | 1000 | RepeatABEL rGLS | ok | 144.280 | 503.97 |
| 500 | 5 | 2500 | 1000 | RepeatABEL rGLS | ok | 226.140 | 766.91 |
| 1000 | 2 | 2000 | 1000 | RepeatABEL rGLS | timeout | 300.010 | 76.05 |

RepeatABEL completed `N=500` but required minutes even for 1000 SNPs and used substantially more memory than CosmicScanner. At `N=1000, records=2, SNPs=1000`, RepeatABEL did not complete within the 300-second timeout. This is consistent with its dense record-level covariance construction and `eigen(V)` implementation.

## Interpretation

RepeatABEL is an important published reference because it implements repeated-measure GWAS with related individuals. However, its source code constructs the dense record-level covariance matrix `V` and performs `eigen(V)`. This gives memory complexity proportional to `n_records^2` and dense eigendecomposition cost proportional to `n_records^3`.

CosmicScanner preserves the same repeated-measure model but reduces the core covariance operations to subject-level matrices using Woodbury decomposition. This is the central scalability advantage. It also supports modern genotype formats and LOCO analysis, which are required for current GWAS workflows.

## Figure and Table Plan

Recommended manuscript figures:

- Figure 1: Model schematic showing records, subjects, incidence matrix `Z`, GRM `G`, and Woodbury reduction.
- Figure 2: Accuracy comparison with RepeatABEL: scatter plots of beta, SE, and `-log10(P)`.
- Figure 3: Runtime and memory benchmark: CS Exact/Woodbury vs RepeatABEL across subjects and records per subject.
- Figure 4: Exact vs Fast/Mean approximation simulation under balanced, unbalanced, and record-level batch/year effects.

Generated figure files from the current benchmark round are stored under `/tmp/opencode/repeatability_figures`:

| Figure file | Content |
|---|---|
| `repeatability_runtime_scaling.png/pdf` | Runtime scaling across CosmicScanner and RepeatABEL benchmark points |
| `repeatability_memory_scaling.png/pdf` | Peak memory scaling across CosmicScanner and RepeatABEL benchmark points |
| `repeatability_concordance_n80_r5_m300.png/pdf` | CS Exact/Woodbury vs RepeatABEL beta, SE, and `-log10(P)` concordance |
| `repeatability_benchmark_plot_data.tsv` | Combined plotting data used for runtime and memory figures |

Recommended tables:

- Table 1: Feature comparison between CosmicScanner and RepeatABEL.
- Table 2: Accuracy metrics against RepeatABEL.
- Table 3: Runtime and memory benchmark.
- Supplementary Table 1: PLINK/PGEN/BGEN consistency matrix.

## Reproducibility Commands

Accuracy comparison:

```bash
ACC_N=80 ACC_M=300 ACC_R=5 ACC_TIMEOUT_SEC=180 \
ACC_REPEATABL_R_LIB=/tmp/opencode/Rlib_repeatabel \
bash CosmicScanner/tests/run_repeatability_repeatabel_accuracy.sh
```

Expanded CS benchmark:

```bash
BENCH_DATA_ROOT=/tmp/opencode/repeatability_paper_bench_cs_data \
BENCH_RESULT_ROOT=/tmp/opencode/repeatability_paper_bench_cs_results \
BENCH_N_LIST=100,200,500 \
BENCH_M_LIST=1000 \
BENCH_R_LIST=2,5 \
BENCH_METHODS=exact,fast,mean \
BENCH_TIMEOUT_SEC=180 \
bash CosmicScanner/tests/run_repeatability_performance_benchmark.sh
```

RepeatABEL benchmark:

```bash
BENCH_DATA_ROOT=/tmp/opencode/repeatability_paper_bench_rb_data \
BENCH_RESULT_ROOT=/tmp/opencode/repeatability_paper_bench_rb_results \
BENCH_N_LIST=100,200 \
BENCH_M_LIST=1000 \
BENCH_R_LIST=2,5 \
BENCH_METHODS=repeatabel \
BENCH_TIMEOUT_SEC=180 \
BENCH_REPEATABL_R_LIB=/tmp/opencode/Rlib_repeatabel \
bash CosmicScanner/tests/run_repeatability_performance_benchmark.sh
```

Second-round CS benchmark on non-`/disk228` storage:

```bash
BENCH_DATA_ROOT=/tmp/opencode/repeatability_paper_bench2_cs_n1000_data \
BENCH_RESULT_ROOT=/tmp/opencode/repeatability_paper_bench2_cs_n1000_results \
BENCH_N_LIST=1000 \
BENCH_M_LIST=10000 \
BENCH_R_LIST=2,5,10 \
BENCH_METHODS=exact,mean \
BENCH_TIMEOUT_SEC=900 \
bash CosmicScanner/tests/run_repeatability_performance_benchmark.sh
```

Second-round CS subject-scaling point:

```bash
BENCH_DATA_ROOT=/tmp/opencode/repeatability_paper_bench2_cs_n2000_data \
BENCH_RESULT_ROOT=/tmp/opencode/repeatability_paper_bench2_cs_n2000_results \
BENCH_N_LIST=2000 \
BENCH_M_LIST=10000 \
BENCH_R_LIST=5 \
BENCH_METHODS=exact,mean \
BENCH_TIMEOUT_SEC=1200 \
bash CosmicScanner/tests/run_repeatability_performance_benchmark.sh
```

Second-round RepeatABEL scaling/timeout benchmark:

```bash
BENCH_DATA_ROOT=/tmp/opencode/repeatability_paper_bench2_rb_n500_data \
BENCH_RESULT_ROOT=/tmp/opencode/repeatability_paper_bench2_rb_n500_results \
BENCH_N_LIST=500 \
BENCH_M_LIST=1000 \
BENCH_R_LIST=2,5 \
BENCH_METHODS=repeatabel \
BENCH_TIMEOUT_SEC=300 \
BENCH_REPEATABL_R_LIB=/tmp/opencode/Rlib_repeatabel \
bash CosmicScanner/tests/run_repeatability_performance_benchmark.sh
```

For systems where the project disk is full, keep benchmark artifacts out of the source tree by setting:

```bash
BENCH_DATA_ROOT=/tmp/opencode/<benchmark_name>_data
BENCH_RESULT_ROOT=/tmp/opencode/<benchmark_name>_results
ACC_DATA_ROOT=/tmp/opencode/<accuracy_name>_data
ACC_RESULT_ROOT=/tmp/opencode/<accuracy_name>_results
```

Generate manuscript figures from existing benchmark and accuracy outputs:

```bash
python3 CosmicScanner/tools/plot_repeatability_benchmark.py \
  --benchmark /tmp/opencode/repeatability_paper_benchmark_combined.tsv \
  --benchmark /tmp/opencode/repeatability_paper_benchmark_round2_combined.tsv \
  --accuracy-cs /tmp/opencode/repeatability_repeatabel_accuracy_r5_results/n80_m300_r5.exact.assoc.txt \
  --accuracy-repeatabel /tmp/opencode/repeatability_repeatabel_accuracy_r5_results/n80_m300_r5.repeatabel.full.tsv \
  --accuracy-label n80_r5_m300 \
  --outdir /tmp/opencode/repeatability_figures
```

## Remaining Manuscript Work

Code development for the repeatability GWAS module is complete. Remaining work is manuscript and benchmark expansion:

- Run final large-scale CS-only benchmarks at `N=1k, 5k, 10k` with `M=10k-100k`.
- Run RepeatABEL until timeout/memory limit to document scaling limitations.
- Generate final plotting scripts for runtime, memory, beta concordance, and log10P concordance.
- Convert this draft into Methods, Results, and Supplementary Methods sections.
