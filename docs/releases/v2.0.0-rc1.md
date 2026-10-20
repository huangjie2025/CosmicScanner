# CosmicScanner v2.0.0-rc1 — Release candidate

**Status: pre-release, not a stable release.** This candidate provides mixed-model association workflows with model-specific support boundaries. Longitudinal RRM and the documented extensions remain experimental.

## Package update

The 2026-10-19 distribution update standardizes the four platform archives:

- Consistent `CosmicScanner-v2.0.0-rc1-<platform>` filenames and top-level directories.
- A shared `PACKAGE_ID.json` schema (version 2) that separates candidate-reference and platform-provenance fields and distinguishes distribution version `2.0.0-rc1` from executable version `2.0.0`.
- Complete Chinese and English manuals in Markdown and Word, bilingual CLI guides, and installation instructions.
- The same six introductory dataset files: `synth.bed`, `synth.bim`, `synth.fam`, `quantitative.pheno`, `multivariate.pheno`, and `long.pheno`.
- Removal of extra Linux developer fixtures from the user distribution.
- Archive and per-file SHA-256 checksums, with original package identity retained under `provenance/`.

**Executables and runtime libraries are unchanged from the supplied 15 October candidates.** This update changes distribution materials, not statistical methods. Original macOS acceptance metadata are retained; a consistent RC1 label does not imply new scientific certification.

## Supported scope

| Family | RC1 scope |
|---|---|
| Single-trait Gaussian LMM | Supported within scope for documented additive LM/LMM/LOCO settings. |
| Multi-trait LMM | BED supported within scope; PGEN/BGEN experimental. The `exact` tier requires no-LOCO. |
| Repeatability LMM | BED, record-level `approx`/`semi` supported within scope; other documented extensions experimental. |
| Longitudinal RRM | Both record-level and coefficient methods remain experimental. |

Other models and formats retain their documented restrictions. File-reading capability does not imply support for every model.

## Checks performed for this packaging update

All four archives were checked for consistent metadata, matching common files, archive and member checksums, and unchanged executable/runtime bytes.

The Windows x86_64 executable was run with the RC1 synthetic fixture: **60 individuals, 180 variants, 180 association rows**, and finite valid additive P values. A separate Windows package self-test completed the single-trait, multivariate, repeatability and GRM examples. The website's new Manhattan and QQ figures use the actual Windows RC1 single-trait output.

**Linux and macOS executables were not rerun in this packaging pass.** Introductory workflow checks do not establish scientific validation of every model or platform. RRM is not included in the four-workflow self-test.

## Downloads

| Platform | Package and checksum | Requirement |
|---|---|---|
| Windows x86_64 | [CosmicScanner-v2.0.0-rc1-windows-x86_64.zip](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-windows-x86_64.zip) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-windows-x86_64.zip.sha256) | Windows 10 1903+ recommended |
| Linux x86_64 | [CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz.sha256) | glibc 2.34 or later |
| macOS Apple Silicon · arm64 | [CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz.sha256) | macOS 15 or later |
| macOS Intel · x86_64 | [CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz.sha256) | macOS 15 or later |

macOS requires version 15 or later and retains the original ad-hoc signature, without Apple Developer ID signing or notarization. Windows ARM64, Linux ARM64 and 32-bit packages are not provided.

`release.json` and `SHA256SUMS` describe these exact archive bytes. Verify the archive, keep the extracted directory intact, and begin with `README.txt` and the manuals.

```text
107fc1b0407a501ad86b1f803ae971b89cb6914d569b651ad346f073f8c06749  CosmicScanner-v2.0.0-rc1-windows-x86_64.zip
1c6206ba86b00205c16e8b13a84da58b83d4cfa52d0886f9e81816a14a2b1f4f  CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz
1281f6271e961e46b28adc207c45027fe238cd560647ce9e81cf0c5c356ab9df  CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz
5d8f41e48242a0c809445991295dd738073d2cec003a1dc67efa6560943891b6  CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz
```

## Source and version identity

- Distribution: `2.0.0-rc1`.
- Executable version: `2.0.0`.
- Candidate-manifest reference revision: `e80062ae31942dc518b8a13110c509470bfe7421`.
- Packaging revision: `rc1-consistent-distribution-20261019`.
- GitHub publication date: `2026-10-20`.
- The source tag includes the candidate source plus the updated manuals, examples, citation metadata and release documentation. No C++ source or build configuration was changed in this documentation update; the original candidate reference remains recorded separately in the distribution manifest and package identities.
- Original candidate date: `2026-10-15`; normalized distribution date: `2026-10-19`.

The reference revision comes from the supplied candidate manifest. Original platform build metadata are preserved separately in each archive. The macOS packages record workflow run `38678035126`, whose public head is `86dce9681ba31e22de72c2b771f1d38e1b5d2a31`; the reference revision follows it with documentation and packaging changes only. A public workflow record is not an independent byte-for-byte reconstruction of the binary, and this packaging pass did not reproduce every platform build.

The historical `v1.1.3` tag points to a different source snapshot. Its commands, source archives, synthetic fixtures and archived manuscript comparisons must not be relabeled as RC1 material. The source archive generated for this pre-release must match the verified RC1 source revision, not an automatically selected historical default-branch commit.

The manuscript remains in preparation. No publication, archive DOI, broad calibration result, or cross-platform performance claim is introduced by this release.

## Feedback

Report issues through the [public tracker](https://github.com/huangjie2025/CosmicScanner/issues) or [huangjie2025@zju.edu.cn](mailto:huangjie2025@zju.edu.cn). Include `PACKAGE_ID.json`, the package checksum, the command, relevant logs and a small shareable synthetic reproducer.

