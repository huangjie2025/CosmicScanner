# Release Package Report — CosmicScanner v1.1.0

**Date**: 2026-08-23  
**Recommended tag**: `v1.1.0`  
**Audit type**: Documentation + packaging only (no code changes)

---

## Recommended Tag Name

```
v1.1.0
```

Consistent across: `CMakeLists.txt`, `VERSION`, `CHANGELOG.md`, `Tutorial_CN.md`, `Tutorial_EN.md`, `dist/` artifacts.

---

## Release Notes Draft

### CosmicScanner v1.1.0 — First Public GWAS Release

**Date**: 2026-07-14 (planned release)  
**Scope**: Focused GWAS release — single-trait, multivariate, repeatability, and relationship-matrix utilities

#### What's Included

- **Single-trait LMM GWAS** (`--gwas`): exact, fast, and semi scan tiers; LOCO-aware; GEMMA-validated (beta/SE/Pearson r > 0.999)
- **Multivariate LMM GWAS** (`--multivariate`): EMmvGWAS-compatible, GEMMA-full, hybrid adaptive with LOCO-aware refinement; Pearson r > 0.99 with GEMMA/EMmvGWAS
- **Repeatability GWAS** (`--repeatability-model`): exact record-level Woodbury reduction (10–20× faster, ~10× lower memory vs RepeatABEL); fast aggregated mode; Pearson r > 0.99 with RepeatABEL
- **A/G/H matrix construction** (`--make-agh`): A-inverse, genomic G, H-inverse for single-step workflows

#### Key Benchmarks

| Comparison | Metric | Value |
|---|---|---|
| CS vs GEMMA (single-trait) | Pearson r(-log10P) | > 0.999 |
| CS vs EMmvGWAS (multivariate) | Pearson r(-log10P) | > 0.999 |
| CS vs RepeatABEL (repeatability) | Pearson r(-log10P) | > 0.999 |
| CS vs GMAT (unbalanced RRM) | Pearson r(-log10P), autosome | 0.952 |
| CS vs GMA (balanced fixreg) | Pearson r(-log10P) | 1.000 |
| Repeatability speedup vs RepeatABEL | Wall time | 10–20× |
| Repeatability memory reduction vs RepeatABEL | Peak RSS | ~10× |
| BLUP-based RRM speedup vs observation-level | Wall time | 213× |
| BLUP-based RRM memory reduction | Peak RSS | 93× |

#### Validation

- 14/14 CTest checks pass
- GEMMA accuracy benchmark: beta_corr=0.99987, se_corr=0.99953, log10p_corr=0.99929
- Multivariate validation matrix: 720/720 finite p-values
- Repeatability edge cases: empty/all-NA rejection, covariates, cache invalidation
- RRM Method A: balanced-path equivalence Pearson ≥ 0.99999; simulation scenarios A–G pass
- Maternal model: 10-scenario simulation suite pass
- Animal model: 7-scenario simulation + GEMMA A22 external validation pass
- ssGWAS Method 2: analytic reference (EBV max diff 5e-7) + simulation (5/5 causal recovery)
- APY-exact: C2-G3 regression pass; APY-approx: D3-D7 boundary outputs retained

#### Input Formats

- PLINK 1.9 (.bed/.bim/.fam)
- PLINK 2.0 (.pgen/.pvar/.psam)
- BGEN (.bgen, Layout 2, zlib compressed)

#### Installation

```bash
# Pre-built binary (recommended)
tar -zxvf CosmicScanner-v1.1.0-linux-x86_64.tar.gz
cd CosmicScanner-v1.1.0-linux-x86_64
bash install.sh

# Source build
git clone --recurse-submodules https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

#### Quick Start

```bash
# GRM construction
cosmicscanner --bfile my_data --make-grm --out my_grm

# Single-trait GWAS (LOCO, semi tier)
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno 1 --gwas --scan-tier semi --out my_results

# Multivariate GWAS
cosmicscanner --pfile my_ukb_data --pheno my_pheno.txt --mpheno-list 1,2,3 --multivariate --gwas --scan-tier semi --out my_mv_results

# Repeatability GWAS
cosmicscanner --bfile my_data --pheno repeated_pheno.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --out repeated_assoc
```

#### What's NOT in This Release

The following are present in the development tree but excluded from the v1.1.0 focused release scope:

- RRM/dynamic GWAS (`--long-pheno`, `--rrm-order`)
- ssGWAS Method 1/Method 2 (`--single-step`)
- SNP P-value inference (`--ssgwas-pvalue-method`)
- Maternal model (`--maternal-model`)
- Animal model (`--animal-model`)
- Common environment model (`--common-env-model`)
- PSR/count/ordinal models
- Survival/GxE models
- Rare-variant tests (SKAT/Burden)
- Biobank/REGENIE-style scans
- BSLMM/DAP interfaces
- Summary-statistic utilities

These remain available in the development tree for research use and will be included in future releases.

#### Known Limitations

- External BLAS/LAPACK disabled by default for cross-platform reproducibility; enable with `COSMICGWAS_ENABLE_BLAS=ON`
- Repeatability exact per-SNP refit may construct full GRM inverse (intended for small datasets)
- Multivariate LOCO mode requires per-chromosome null-model refit (slower than `--no-loco`)

#### Citation

If you use CosmicScanner in your research, please cite:

> Huang J, et al. (2026). CosmicScanner: a unified C++ framework for genome-wide association analysis across diverse complex-trait architectures. *[Journal TBD]*.

#### Links

- Documentation: [docs/Tutorial_CN.md](docs/Tutorial_CN.md) | [docs/Tutorial_EN.md](docs/Tutorial_EN.md)
- Validation: [docs/FINAL_VALIDATION_REPORT.md](docs/FINAL_VALIDATION_REPORT.md)
- Model matrix: [docs/PAPER_MODEL_SUPPORT_MATRIX.md](docs/PAPER_MODEL_SUPPORT_MATRIX.md)
- Source: https://github.com/huangjie2025/CosmicScanner

---

## Remaining Non-Blocking Items

| # | Item | Severity | Recommendation |
|---|---|---|---|
| 1 | README missing version number in header | Low | Add `**Version: v1.1.0**` to README header |
| 2 | README missing Installation section | Medium | Add dedicated build/install section or link to tutorials |
| 3 | README missing Citation section | Medium | Add BibTeX entry or citation text |
| 4 | README missing Benchmark section | Low | Add dedicated benchmarks subsection with comparison table |
| 5 | README missing Example commands for `examples/` | Medium | Add quick-start commands using `examples/mouse_clean/` |
| 6 | `--help` missing version string | Low | Add version display to `--help` header |
| 7 | LaTeX compilation unverified locally | Low | `latexmk`/`pdflatex` not available on build host |

**None of these items are blocking for the v1.1.0 tag.** All are documentation polish that can be addressed in a v1.1.1 patch if needed.

---

## Audit Verdict

**Release candidate is ready for tagging.** All version strings are consistent, build artifacts exist with SHA-256 checksums, documentation is aligned with canonical status labels, validation numbers are unchanged, and no blocking issues were found.
