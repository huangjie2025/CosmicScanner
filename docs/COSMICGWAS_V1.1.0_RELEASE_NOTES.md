# CosmicGWAS v1.1.0 Release Notes

Date: 2026-07-14

## Release Scope

CosmicGWAS v1.1.0 is a focused release for three Gaussian LMM GWAS workflows and one relationship-matrix utility:

- Single-trait GWAS using linear mixed models.
- Multivariate GWAS using multivariate mixed models.
- Repeatability GWAS for repeated-measure phenotypes, with exact record-level analysis as the main route.
- A-inverse, genomic G, and H-inverse matrix construction.

The v1.1.0 release package intentionally excludes broader experimental or future modules, including RRM/dynamic GWAS, PSR/count/ordinal models, ssGWAS association scans, maternal models, survival models, GxE reaction norms, rare-variant tests, biobank/REGENIE-style scans, BSLMM/DAP interfaces, and summary-statistic utilities.

## Version History

CosmicGWAS v1.1.0 is the first public, GWAS-focused release. Earlier numbers below are retrospective internal milestones reconstructed from repository history; they are not Git tags and must not be represented as prior public releases.

- `1.0.0` (internal milestone, 2026-04-07): independent GWAS baseline with single-trait LMM/LOCO, an initial multivariate LMM path, and an aggregated repeatability prototype.
- `1.0.5` (internal milestone, 2026-05-22 to 2026-06-24): dedicated repeatability model with exact and fast routes, expanded multivariate compatibility paths, and consolidated reference and regression validation.
- `v1.1.0` (public release, 2026-07-14): focused and validated source and Linux packages for the four supported functions.

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
- Repeatability exact with a within-subject centered record-level covariate, guarding the zero-right-hand-side VCE solve.
- ID-aligned A-inverse, G, and H-inverse matrix construction.
- SHA-256 checksum verification for both source and binary archives.

The source package defaults to Eigen-native solvers and disables external BLAS/LAPACK unless explicitly enabled with `COSMICGWAS_ENABLE_BLAS=ON bash scripts/build.sh`.

The default repeatability exact/fast Woodbury routes retain the subject GRM without forming its full inverse. Per-SNP exact-refit validation may construct a full inverse and is intended for small datasets or selected variants.

Repeatability null summaries report the number of unique subjects retained in the phenotype analysis rather than the total number of samples in the genotype file.

## Release Readiness

The v1.1.0 artifacts meet the software packaging gate for public release: package construction, archive checksums, clean source rebuild, and packaged workflow validation have passed. This status establishes operational release readiness within the documented four-function scope; it does not claim exhaustive statistical validation for every cohort design or data distribution.
