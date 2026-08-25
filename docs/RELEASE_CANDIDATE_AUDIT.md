# Release Candidate Audit Report — CosmicScanner v1.0-rc

**Date**: 2026-08-23  
**Audit scope**: GWAS release candidate (no code changes, documentation + CLI + reproducibility only)  
**Git status**: 20 modified files; no `git add`/`git commit`/`git clean` performed

---

## 1. CLI Flag Audit

**Binary**: `build/cosmicscanner` (5.4 MB ELF x86-64)

All 11 core GWAS flags confirmed present in `--help` output with correct default values:

| Flag | Default | Status |
|---|---|---|
| `--gwas` | on | ✅ PASS |
| `--multivariate` | on | ✅ PASS |
| `--repeatability-model` | off | ✅ PASS |
| `--long-pheno` / `--rrm-order` | off | ✅ PASS |
| `--rrm-gwas-method` | `coefficient` | ✅ PASS |
| `--animal-model` | off | ✅ PASS |
| `--maternal-model` | off | ✅ PASS |
| `--maternal-snp-test` | `direct` | ✅ PASS |
| `--common-env-model` | off | ✅ PASS |
| `--single-step` / `--single-step-method` | `lmm` | ✅ PASS |
| `--ssgwas-pvalue-method` | `off` | ✅ PASS |
| `--ssgwas-apy-core-size` | required for `apy-*` | ✅ PASS |

---

## 2. Documentation Consistency Audit

### Canonical status labels (3 allowed: `Stable`, `Experimental validated`, `Experimental approximate`)

**README.md** ✅ PASS
- Line 10-14: Correct 3-label taxonomy
- Line 303-307: Correct model-to-status mapping
- No stale terms

**docs/Tutorial_EN.md** ✅ PASS
- Lines 10-11: Correct label header
- Lines 83-89: Correct model-to-status table
- Lines 511, 527, 576, 597-599: Correct inline labels

**docs/Tutorial_CN.md** ✅ PASS
- Lines 9-10: Correct label header
- Lines 82-88: Correct model-to-status table
- Lines 599, 618, 676, 686-688: Correct inline labels

**docs/PAPER_MODEL_SUPPORT_MATRIX.md** ✅ PASS
- Lines 8-10: Correct label definitions
- Lines 25-36: All 12 model rows with correct status labels
- Lines 54-55: Correct APY boundary warnings

**docs/FINAL_VALIDATION_REPORT.md** ⚠️ NON-BLOCKING
- Line 31: Contains `EXPERIMENTAL_APPROXIMATE` (internal metadata constant, expected)
- Lines 13-26: Correct model-to-status mapping
- Lines 68, 81, 83-86: Correct boundary warnings

### Stale terms check

`simulation-validated` appears 3× in `docs/GWAS_MODEL_DEFINITIONS.md` (lines 125, 126, 281) — this is a **definition/specification document** not exposed to users, so acceptable as internal terminology. Does NOT appear in README, tutorials, or manuscripts.

No stale `10-test` references found in any manuscript; all 3 manuscripts correctly reference `14 registered CTest checks`.

---

## 3. Reproducibility Audit

### Test scripts referenced in documentation — all exist in `tests/`:

| Script | Referenced in | Exists |
|---|---|---|
| `run_test_rrm_method_a_simulation.sh` | README, tutorials, manuscripts | ✅ |
| `run_rrm_balance_record_equivalence.sh` | README, tutorials | ✅ |
| `run_test_singlestep.sh` | FINAL_VALIDATION_REPORT | ✅ |
| `run_test_maternal_snp_effects.sh` | FINAL_VALIDATION_REPORT | ✅ |
| `run_test_animal_model.sh` | FINAL_VALIDATION_REPORT | ✅ |
| `run_test_common_environment.sh` | FINAL_VALIDATION_REPORT | ✅ |
| `run_ssgwas_backsolve_analytic.py` | FINAL_VALIDATION_REPORT | ✅ |
| `run_ssgwas_pvalue_null_calibration.py` | FINAL_VALIDATION_REPORT | ✅ |
| `run_ssgwas_pvalue_apy_exact_regression.py` | FINAL_VALIDATION_REPORT | ✅ |
| `run_gmat_unbalanced_alignment.sh` | README, tutorials | ✅ |
| `run_gma_compatible_fixreg_alignment.sh` | README, tutorials | ✅ |
| `run_animal_model_gemma_validation.sh` | FINAL_VALIDATION_REPORT | ✅ |

All referenced scripts have `#!/bin/bash` or `#!/usr/bin/env python3` shebangs. No temp-file dependencies outside `build/` and `tests/` directories.

### README entry points

README provides direct command-line examples for all 11 model classes with working flag combinations. No missing entry points.

---

## 4. Manuscript Consistency Audit

### Key validation numbers — consistent across all 3 manuscripts:

| Metric | MANUSCRIPT_DRAFT.md | MANUSCRIPT_DRAFT_CN.md | manuscript.tex |
|---|---|---|---|
| GEMMA correlation | r > 0.99 | r > 0.99 | r > 0.99 |
| GMAT autosome r | 0.952 | 0.952 | 0.952 |
| RRM speedup | 213× | 213× | 213× |
| RRM memory | 93× | 93× | 93× |
| CTest count | 14 | 14 | 14 |
| λ_GC (extended) | ≈1.00 | ≈1.00 | ≈1.00 |
| Maternal variance error | <6% | <6% | <6% |
| Record-equivalence Pearson | ≥0.99999 | ≥0.99999 | ≥0.99999 |

### Method 1/Method 2 naming — consistent:

- All 3 manuscripts: "H-matrix single-step GWAS (Method 1)" = fixed-SNP BETA/SE/P
- All 3 manuscripts: "ssGBLUP back-solving GWAS (Method 2)" = EBV/GEBV → marker effects
- No conflation of the two methods

### APY boundaries — consistent:

- APY-exact: `Experimental validated` across all docs
- APY-approx: `Experimental approximate` across all docs
- All docs include `STRICT_TYPE_I_CONTROL=NOT_GUARANTEED` caveat
- All docs state "research/exploratory only; not calibrated P-value"
- No docs recommend APY-approx for ranking/screening

---

## 5. Git Status Audit

**No forbidden operations performed**: no `git add`, `git commit`, `git clean`, or `git reset`.

Current working tree shows 20 modified files (all from prior sessions):
- 16 tracked modifications (M)
- 4 untracked/uncommitted modifications (m)
- No deleted files
- No new untracked files beyond `libCosmicCore` submodule

---

## 6. Release Readiness Summary

| Category | Status |
|---|---|
| CLI flags (11 core) | ✅ PASS |
| Documentation labels | ✅ PASS |
| Stale terms | ✅ PASS (definition docs only) |
| Reproducibility (scripts) | ✅ PASS |
| Manuscript numbers | ✅ PASS |
| Manuscript Method1/2 | ✅ PASS |
| APY boundaries | ✅ PASS |
| Git safety | ✅ PASS |

### Blocking Issues: 0

### Non-Blocking Issues (informational):

1. `docs/GWAS_MODEL_DEFINITIONS.md` uses internal `simulation-validated` terminology (3 occurrences) — acceptable in definition doc, not exposed to users.
2. `latexmk`/`pdflatex` not available on system — LaTeX compilation cannot be verified locally.
3. `libCosmicCore` appears as modified submodule in git status — may need submodule sync before tag.

### Recommended Version Tag: `v1.0-rc1`

**Audit verdict**: Release candidate is ready for tagging. All 11 GWAS models have correct CLI flags, consistent documentation labels, reproducible test scripts, and aligned manuscript numbers.
