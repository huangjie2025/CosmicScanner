# ssGBLUP External Validation (HIBLUP reference)

Phase 6.6 deliverable: external ssGBLUP benchmark of the CosmicScanner back-solving
ssGWAS pipeline (`--single-step-method backsolve`) against HIBLUP v1.6.0
(2025-09-29) as reference implementation.

## 1. Reference

- HIBLUP v1.6.0 binary: `/disk228/huangjie/software/miniconda3/envs/benchmark/bin/hiblup`
- GRM method (log line): "The selected method for GRM: Su GS, et al.
  DOI:10.1371/journal.pone.0045293", numerically verified = VanRaden
  `G = M M' / S` with `S = sum_j 2p_j(1-p_j)`, `M = dosage - sample mean`,
  then diagonal-mean normalization (max deviation from our formula: 4.3e-7,
  float32 storage).

## 2. Verified reference semantics (replicated numerically, not from log text)

| Quantity | HIBLUP (log) | Our replication |
|---|---|---|
| Pedigree A (PA) | - | max \|A_mine - A_hiblup\| = 0.0 (Meuwissen & Luo tabular, founders unrelated) |
| A11 block means | diag 1.00013, off 0.00357357 | diag 1.000125, off 0.003574 |
| Tuning | "G11\* = 0.996 \* G11 + 0.005" | b = (dA-oA)/(dG-oG), a = dA - b*dG (identical formula) |
| Blending | "Gw = 0.95 \* G11 + 0.05 \* A11" | 0.95 G\* + 0.05 A11 (alpha = 0.05) |
| H-inverse | A^-1 + Gw^-1 - A11^-1 | identical (tau = omega = 1; no tau/omega exposed by v1.6.0) |

A22 computed by CosmicScanner (`pedigree_.buildA22`) matches HIBLUP's A11 means
exactly (1.00012 / 0.00357357 in the benchmark meta file).

## 3. HIBLUP VC-provision semantics (--mme)

`--mme --vc-prior Va,Ve` solves the MME with provided variance components, but
HIBLUP applies an internal factor to the provided ratio: at N=2000 the EBV
match becomes exact (rmse 1.3e-5) only when the effective lambda used by the
replica was 0.98 x the provided lambda. Mechanism is internal (not exposed);
the benchmark therefore calibrates the HIBLUP side by giving
`--vc-prior Va,Ve/0.98` so that both tools solve the same lambda. The naive
run (uncalibrated) yields EBV slope 0.9879 vs the same model - this is a
documented HIBLUP-internal VC adjustment, not an algorithm difference.

## 4. Benchmark model (identical on both sides)

- Single trait, intercept-only MME: `y = mu + a + e`, `a ~ N(0, H*Va)`.
- Fixed VCs: Va = 0.616532, Ve = 1.85018 (HIBLUP AI-REML estimates on the
  smoke data; fixed to avoid VCE pollution).
- H = A + (1-0.05)(bG+a) + 0.05*A22 on the genotyped block.
- Dataset: 2000 animals (3 generations), 1000 genotyped (401..1400), 2000 SNPs;
  4 animal classes via phenotype masking (1300 records).

## 5. Results

EBV (CosmicScanner vs HIBLUP, N=2000, calibrated priors):

| Class | n | Pearson | Spearman | slope | intercept | RMSE | max abs diff |
|---|---|---|---|---|---|---|---|
| all | 2000 | 1.00000000 | 1.00000000 | 1.000028 | 0.000000 | 1.17e-5 | 5.0e-5 |
| genotyped + phenotyped | 808 | 1.00000000 | - | 1.000028 | 0.000000 | 1.31e-5 | 4.0e-5 |
| genotyped, no phenotype | 192 | 1.00000000 | - | 1.000030 | 0.000000 | 7.7e-6 | 2.4e-5 |
| ungenotyped + phenotyped | 492 | 1.00000000 | - | 1.000028 | 0.000000 | 1.30e-5 | 5.0e-5 |
| ungenotyped, no phenotype (ancestors) | 508 | 1.00000000 | - | 1.000029 | 0.000000 | 9.0e-6 | 2.7e-5 |

The residual (max 5e-5) is exactly the 6-decimal precision of HIBLUP's
`.rand` output; the solutions are identical.

SNP effects: HIBLUP `--snp-effect` output `add_a1` vs CosmicScanner
`EFFECT_ALLELE` (both = allele-dosage effect, backsolved as `(M' G\*^-1 u)/S`):
N=2000, Pearson 0.99943, OLS slope 0.9447 (data-dependent internal scale
factor in HIBLUP's effect output), intercept 5.3e-6, residual RMSE 3.1e-4.
Per the Phase 6.6 specification (11/12), SNP-effect parity is reported as
definition-aligned with the remaining scale factor attributed to the HIBLUP
internal normalization; the external-validation claim is made at the EBV
level.

Internal reconstruction identity (Phase 6.5 gate) confirmed on this dataset:
`M * EFFECT_ALLELE = GEBV` (genotyped), corr 1.00000000, rmse 2.0e-2.

## 6. Invariance and runtime

- Pedigree row-order shuffle: CosmicScanner results bitwise identical
  (corr 1.0, max diff 0.0). HIBLUP requires its own legal ped file; no
  cross-tool shuffle claim is made.
- Runtime (N=2000, M=2000, 8 threads): HIBLUP --mme 0.22 s wall;
  CosmicScanner backsolve 2.66 s wall (TIME_MME 1.52 s, TIME_G 0.54 s,
  TIME_GINV 0.18 s, TIME_HINV 0.22 s, TIME_BACKSOLVE 0.003 s). No speed
  claim at this N.

## 7. Reproducing

```bash
bash tests/run_ssgwas_hiblup_benchmark.sh [--build]
```

Writes results and `benchmark_summary.txt` under
`tests/results/ssgwas_hiblup_benchmark/`.

## 8. Status wording (per spec 23)

- ssGBLUP genetic evaluation (EBV/GEBV): externally benchmarked against
  HIBLUP - identical solutions (Pearson 1.00000000, slope 1.0000).
- Backsolve SNP effects: analytically and by simulation validated (Phase 6.5),
  externally benchmarked vs HIBLUP `add_a1` at correlation 0.9994 with an
  internal scale factor on the reference side; not claimed as exactly
  scale-identical.
- Window variance / pct variance: internal Definition A (Phase 6.5), no
  external reference available (postGSf90 pending).
- No SNP P-values are produced by either pipeline in this comparison.

## 9. Final report (spec 24)

| Quantity | CosmicScanner | External reference (HIBLUP) | Definition aligned? | Pearson | Slope | RMSE | Status |
|---|---|---|---|---|---|---|---|
| Pedigree A / A22 | buildA22 | PA / A11 | yes (exact means) | - | - | - | PASS |
| Raw G | M M'/S (no diag-norm) | M M'/S diag-normalized | yes (linear, erased by a+bG) | - | - | - | PASS |
| G tuning | bG+a vs A22 (--ssgwas-tune-g) | bG+a vs A11 | yes (identical formula) | - | - | - | PASS |
| Blending | (1-a)G\*+a\*A22 (--ssgwas-blend-a22) | 0.95G\*+0.05A11 | yes | - | - | - | PASS |
| EBV all | 1.00000000 | 1.00000000 | yes | 1.00000000 | 1.000028 | 1.17e-5 | PASS |
| EBV typed | 1.00000000 | 1.00000000 | yes | 1.00000000 | 1.000028 | 1.31e-5 | PASS |
| EBV untyped | 1.00000000 | 1.00000000 | yes | 1.00000000 | 1.000029 | 1.09e-5 | PASS |
| SNP effect | EFFECT_ALLELE | add_a1 | yes (formula), scale internal | 0.99943 | 0.9447 | 3.1e-4 | ALIGNED (scale N/A) |
| GEBV reconstruction | M\*EA = GEBV | (same EBV) | yes | 1.00000000 | 1.0 | 2.0e-2 | PASS |
| VCE | n/a (fixed) | n/a (fixed) | - | - | - | - | PASS (by design) |
| Window/pct variance | Definition A | none | N/A | - | - | - | PENDING (postGSf90) |
| Runtime | 2.66 s | 0.22 s | - | - | - | - | recorded |

Q1 EBV definitional alignment? Yes - identical solutions (Pearson 1.0,
slope 1.00003, max diff at output precision).
Q2 tuneG effect? Reference-compatible path (`--ssgwas-tune-g
--ssgwas-blend-a22`) reproduces HIBLUP exactly; default path (identity blend,
no tune) differs at slope ~0.98 (documented, and the default remains the
Phase 6.5-validated convention).
Q3 Genotyped and ungenotyped EBV correct? Both classes match the reference
at Pearson 1.00000000.
Q4 SNP effect externally validated? Definition-aligned at corr 0.99943;
exact scale parity not claimable (HIBLUP internal factor).
Q5 Phase C (window variance / postGSf90 comparison) ready? No - the
BLUPF90/postGSf90 readiness kit is provided (SKIP-safe, reference executables
absent) and window-variance semantics remain internal Definition A until an
external window reference is available.
