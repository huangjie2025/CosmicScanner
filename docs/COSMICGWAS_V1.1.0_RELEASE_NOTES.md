# CosmicGWAS v1.1.0 Release Notes

Date: 2026-07-13

## Release Scope

CosmicGWAS v1.1.0 is a focused GWAS release for three workflows:

- Single-trait GWAS using linear mixed models.
- Multivariate GWAS using multivariate mixed models.
- Repeatability GWAS for repeated-measure phenotypes, with exact record-level analysis as the main route.

The v1.1.0 release package intentionally excludes broader experimental or future modules, including RRM/dynamic GWAS, PSR/count/ordinal models, ssGWAS, maternal models, survival models, GxE reaction norms, rare-variant tests, biobank/REGENIE-style scans, BSLMM/DAP interfaces, and summary-statistic utilities.

## Release Artifacts

The validated artifacts are under `dist/`:

- `CosmicGWAS-v1.1.0-linux-x86_64.tar.gz`: Linux x86_64 binary package.
- `CosmicGWAS-v1.1.0-linux-x86_64.tar.gz.sha256`: binary package checksum.
- `CosmicGWAS-v1.1.0-source.tar.gz`: focused source package for the manuscript-facing release.
- `CosmicGWAS-v1.1.0-source.tar.gz.sha256`: source package checksum.

## Source Package Layout

The focused source package uses a flat implementation layout:

- `src/`: one flat source directory, with no `libCosmicCore/`, `core/`, or `models/` subtree.
- `scripts/build.sh`: portable release build script.
- `scripts/validate_release.sh`: packaged validation script.
- `docs/`: manuscript-facing scope, methods, validation, reproducibility, and input/output notes.
- `examples/` and `test/`: minimal data and packaged smoke validation assets.

## Validation Status

The final source tarball was unpacked and validated from the archive. The validation covered:

- Single-trait LMM association output.
- Multivariate LMM joint association output.
- Repeatability exact association, null-model, and cache-summary outputs.
- Repeatability fast comparison output.
- SHA-256 checksum verification for both source and binary archives.

The source package defaults to Eigen-native solvers and disables external BLAS/LAPACK unless explicitly enabled with `COSMICGWAS_ENABLE_BLAS=ON bash scripts/build.sh`.

## Manuscript Status

The software release package is prepared and validated. Manuscript submission still requires additional evidence:

- Final type-I error and power simulation matrix.
- Exact repeatability versus subject-mean aggregation comparison.
- Record-imbalance and record-level covariate sensitivity analysis.
- Final RepeatABEL reference table with allele harmonization rules.
- Publishable real repeated-measure application.
- Final manuscript title, journal target, DOI, and citation metadata.
