# Release Checklist — CosmicScanner v1.1.0

**Date**: 2026-08-23  
**Target tag**: `v1.1.0`  
**Scope**: Focused GWAS release — single-trait, multivariate, repeatability, A/G/H matrix

---

## 1. Version Consistency

| Source | Version | Status |
|---|---|---|
| `CMakeLists.txt` line 2 | `1.1.0` | ✅ |
| `VERSION` file | `v1.1.0` | ✅ |
| `CHANGELOG.md` line 30 | `v1.1.0` (2026-07-13) | ✅ |
| `docs/Tutorial_CN.md` line 4 | `1.1.0` | ✅ |
| `docs/Tutorial_EN.md` line 4 | `1.1.0` | ✅ |
| `dist/` artifacts | `CosmicScanner-v1.1.0-*` | ✅ |
| README.md header | **MISSING** | ⚠️ |

**Issue**: `README.md` does not display the version number in its header or metadata.

---

## 2. Build Status

| Item | Status |
|---|---|
| Binary exists | ✅ `build/cosmicscanner` (5.4 MB ELF x86-64) |
| CMake version injection | ✅ `COSMICSCANNER_VERSION="${PROJECT_VERSION}"` |
| CTest registered | ✅ 14 tests |
| `--help` output | ✅ All 11 core GWAS flags present |
| `--help` version display | ❌ No version string in `--help` header |

---

## 3. Supported Models

| Model | CLI Flag | Status |
|---|---|---|
| General LMM GWAS | `--gwas` | ✅ Stable |
| Multivariate GWAS | `--multivariate` | ✅ Stable |
| Repeatability GWAS | `--repeatability-model` | ✅ Experimental validated |
| RRM Method A | `--rrm-gwas-method record` | ✅ Experimental validated |
| RRM Method B | `--rrm-gwas-method coefficient` | ✅ Stable/compatibility |
| Animal Model GWAS | `--animal-model --ped` | ✅ Experimental validated |
| Maternal GWAS | `--maternal-model` | ✅ Experimental validated |
| Common Environment GWAS | `--common-env-model` | ✅ Experimental validated |
| H-matrix ssGWAS (Method 1) | `--single-step --single-step-method lmm` | ✅ Experimental validated |
| ssGBLUP back-solving (Method 2) | `--single-step-method backsolve` | ✅ Experimental validated |
| Exact SNP P-value | `--ssgwas-pvalue-method exact` | ✅ Stable |
| APY-exact SNP P-value | `--ssgwas-pvalue-method apy-exact` | ✅ Experimental validated |
| APY-approx SNP P-value | `--ssgwas-pvalue-method apy-approx` | ✅ Experimental approximate |

---

## 4. Validation Status

| Suite | Last Run | Status |
|---|---|---|
| GEMMA accuracy | beta_corr 0.9998708 | ✅ PASS |
| Multivariate validation | standalone PASS | ✅ PASS |
| Repeatability exact/fast/LOCO | smoke PASS | ✅ PASS |
| RRM Method A (scenarios A–G) | simulation PASS | ✅ PASS |
| Maternal SNP suite (A–H) | simulation PASS | ✅ PASS |
| Animal model (7 scenarios) | simulation PASS | ✅ PASS |
| Common environment | simulation+LOCO PASS | ✅ PASS |
| ssGWAS Method 1 smoke | PASS | ✅ PASS |
| ssGWAS Method 2 backsolve | analytic+simulation PASS | ✅ PASS |
| Exact SNP P-value (C4/C5) | PASS | ✅ PASS |
| APY-exact (C2-G3) | PASS | ✅ PASS |
| APY-approx boundary (D3-D7) | PASS (boundary) | ✅ PASS |
| CTest 14/14 | PASS | ✅ PASS |
| GMAT alignment (r=0.952) | PASS | ✅ PASS |
| GMA alignment (r=1.000) | PASS | ✅ PASS |

---

## 5. Documentation Status

| Document | Exists | Version Match | Status |
|---|---|---|---|
| `README.md` | ✅ | ⚠️ No version | ⚠️ |
| `CHANGELOG.md` | ✅ | ✅ v1.1.0 | ✅ |
| `docs/Tutorial_CN.md` | ✅ | ✅ 1.1.0 | ✅ |
| `docs/Tutorial_EN.md` | ✅ | ✅ 1.1.0 | ✅ |
| `docs/PAPER_MODEL_SUPPORT_MATRIX.md` | ✅ | N/A | ✅ |
| `docs/FINAL_VALIDATION_REPORT.md` | ✅ | N/A | ✅ |
| `docs/COSMICGWAS_V1.1.0_RELEASE_NOTES.md` | ✅ | ✅ v1.1.0 | ✅ |
| `docs/MANUSCRIPT_DRAFT.md` | ✅ | N/A | ✅ |
| `docs/MANUSCRIPT_DRAFT_CN.md` | ✅ | N/A | ✅ |
| `docs/manuscript.tex` | ✅ | N/A | ✅ |
| `docs/references.bib` | ✅ | N/A | ✅ |

---

## 6. Examples Status

| Item | Status |
|---|---|
| `examples/mouse_clean/` exists | ✅ |
| `.bed/.bim/.fam` PLINK data | ✅ |
| `.pheno` phenotype file | ✅ |
| Example commands in README | ❌ MISSING |
| Expected outputs documented | ❌ MISSING |
| Path validity (relative paths) | ✅ |

---

## 7. Release Artifacts

| Artifact | Exists | SHA-256 |
|---|---|---|
| `dist/CosmicScanner-v1.1.0-linux-x86_64.tar.gz` | ✅ | ✅ |
| `dist/CosmicScanner-v1.1.0-source.tar.gz` | ✅ | ✅ |
| `install.sh` | ✅ | N/A |
| `docs/COSMICGWAS_V1.1.0_RELEASE_NOTES.md` | ✅ | N/A |

---

## 8. Blocking Issues

**None.**

---

## 9. Non-Blocking Issues

1. **README missing version number**: Header does not display `v1.1.0`.
2. **README missing Installation section**: No dedicated build/install instructions in README (only in tutorials).
3. **README missing Citation section**: No citation or BibTeX entry in README.
4. **README missing Benchmark section**: GMAT/GMA comparison data exists inline but no dedicated section.
5. **README missing Example commands for `examples/`**: `mouse_clean` data exists but no commands reference it.
6. **`--help` missing version string**: Binary help output does not show version number.
7. **`latexmk`/`pdflatex` unavailable**: LaTeX compilation cannot be verified locally.
