# CosmicScanner

[![Build and regress](https://github.com/huangjie2025/CosmicScanner/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/huangjie2025/CosmicScanner/actions/workflows/ci.yml)

CosmicScanner is a command-line toolkit for mixed-model genome-wide association analysis of continuous traits, correlated traits, and repeated measurements.

**Current candidate: v2.0.0-rc1.** This is a release candidate for evaluation, not a stable release. Support depends on the model, input format, and analysis settings. Longitudinal random regression remains experimental.

**[Download V2.0 RC1 for Windows, Linux and macOS](https://github.com/huangjie2025/CosmicScanner/releases/tag/v2.0.0-rc1)** · [Supported configurations](docs/SUPPORTED_MODES.md) · [Automatic builds and verification](docs/CONTINUOUS_INTEGRATION.md)

## Model support

| Phenotype structure | Model | Current scope |
|---|---|---|
| One continuous measurement per individual | Single-trait LMM | Supported within scope for documented standard additive LM/LMM/LOCO configurations. |
| Several correlated continuous traits per individual | Multi-trait LMM | Supported within scope with BED input. PGEN/BGEN are experimental; the `exact` scan tier requires no-LOCO. |
| Repeated measurements of the same trait | Repeatability LMM | Supported within scope with BED input and record-level `approx`/`semi` scans. Fast, STCG and full-refit extensions are experimental. |
| A genetic trajectory over time | Longitudinal RRM | Experimental. Record-level and two-stage coefficient methods are separate analyses with limited validation. |

Input-format support is model-specific. An available command or a successful run does not establish statistical calibration for every study design. Read the model-specific manual sections before analyzing your own data.

## Download and verify

Choose a package from the [v2.0.0-rc1 pre-release](https://github.com/huangjie2025/CosmicScanner/releases/tag/v2.0.0-rc1). Each contains the executable and required bundled runtimes, complete Chinese and English manuals, bilingual CLI guides, and the same synthetic examples.

| Platform | Package and checksum | Requirement |
|---|---|---|
| Windows x86_64 | [CosmicScanner-v2.0.0-rc1-windows-x86_64.zip](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-windows-x86_64.zip) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-windows-x86_64.zip.sha256) | Windows 10 1903+ recommended |
| Linux x86_64 | [CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz.sha256) | glibc 2.34 or later |
| macOS Apple Silicon · arm64 | [CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz.sha256) | macOS 15 or later |
| macOS Intel · x86_64 | [CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz) · [SHA-256](https://github.com/huangjie2025/CosmicScanner/releases/download/v2.0.0-rc1/CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz.sha256) | macOS 15 or later |

Keep the entire extracted directory. Analysis does not require Python, Conda, Homebrew, or a compiler; Python 3 is used only for the optional package self-test. macOS packages are ad-hoc signed and are not Apple Developer ID signed or notarized.

Download the matching `.sha256` file and compare its hash before extraction. For example, on Linux:

```bash
sha256sum -c CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz.sha256
```

Use `shasum -a 256 -c` on macOS, or `Get-FileHash -Algorithm SHA256` in Windows PowerShell.

## Quick start

From the extracted package directory on Linux or macOS:

```bash
./bin/cosmicscanner --version
mkdir -p demo_results
./bin/cosmicscanner gwas --model lmm \
  --bfile examples/dataset/synth \
  --pheno examples/dataset/quantitative.pheno \
  --mpheno 1 --threads 1 --out demo_results/lmm
```

In Windows PowerShell:

```powershell
New-Item -ItemType Directory -Force -Path demo_results | Out-Null
.\bin\cosmicscanner.exe gwas --model lmm --bfile examples/dataset/synth --pheno examples/dataset/quantitative.pheno --mpheno 1 --threads 1 --out demo_results/lmm
```

The bundled fixture contains **60 individuals and 180 variants**. The association table `demo_results/lmm.assoc.txt` should contain a header and 180 variant rows. Read `BETA_A`, `SE_A`, and `P_A` by column name and inspect the model diagnostics and run manifest. This small synthetic workflow demonstrates execution; it is not scientific or performance validation.

## Manuals

Read the [complete English manual](docs/USER_MANUAL_EN.md), [Chinese manual](docs/USER_MANUAL_CN.md), [English CLI guide](docs/CLI_GUIDE_EN.md), or [Chinese CLI guide](docs/CLI_GUIDE_CN.md) directly in this repository. Word editions are in [docs/](docs/).

Inside each extracted package:

- `docs/USER_MANUAL_EN.md`: complete English manual.
- `docs/USER_MANUAL_CN.md`: complete Chinese manual.
- `docs/USER_MANUAL_EN_QUICK.md`: concise English quick start.
- `docs/CLI_GUIDE_EN.md` and `docs/CLI_GUIDE_CN.md`: command-line guides.
- `docs/INSTALL.md`: platform installation instructions.
- `docs/CosmicScanner_v2.0.0_User_Manual_EN.docx` and `docs/CosmicScanner_v2.0.0_User_Manual_CN.docx`: complete Word manuals.

Use `--help-all` for the executable's command reference. Trait indices start after FID and IID.

## Versions and reproducibility

The distribution version is **2.0.0-rc1**, while the unchanged executable reports **2.0.0**. Retain the archive checksum and `PACKAGE_ID.json` with each analysis.

The candidate manifest records reference revision [`e80062ae31942dc518b8a13110c509470bfe7421`](https://github.com/huangjie2025/CosmicScanner/commit/e80062ae31942dc518b8a13110c509470bfe7421). The 2026-10-19 packaging update standardizes names, documentation, examples and metadata without changing executables or runtime libraries. The package identity records this candidate reference separately from platform build provenance. Original platform metadata are retained under `provenance/`; this packaging update did not independently rebuild every platform binary. The macOS packages record a successful original workflow whose head is `86dce9681ba31e22de72c2b771f1d38e1b5d2a31`; the reference revision adds only documentation and packaging changes.

The default branch includes the RC1 source and updated distribution documentation. The `v2.0.0-rc1` tag includes that candidate source with the updated distribution documentation. Its C++ source and build configuration are unchanged from the candidate reference revision above; the tag also preserves the manuals, examples and release records for this publication. Prebuilt distributions retain their separately recorded packaging and platform-build identities. The [historical v1.1.3 tag](https://github.com/huangjie2025/CosmicScanner/tree/v1.1.3) and its examples belong to a different software snapshot and must not be used as RC1 source or evidence.

Archived manuscript comparisons are historical results, not validation of RC1. The manuscript is in preparation; no published paper, preprint or DOI is claimed here. Report the software version, package checksum, command and analysis settings in your methods.

## Release records

See the [release notes](docs/releases/v2.0.0-rc1.md), [distribution manifest](docs/releases/v2.0.0-rc1.json), and [changelog](CHANGELOG.md).

## Support and license

Use the [public issue tracker](https://github.com/huangjie2025/CosmicScanner/issues) or contact [huangjie2025@zju.edu.cn](mailto:huangjie2025@zju.edu.cn). Include the package version, command, relevant log and a small synthetic reproducer.

CosmicScanner is distributed under the [MIT license](LICENSE). Bundled third-party libraries retain their own licenses in the package.
