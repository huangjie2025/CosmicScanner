# Automatic builds and minimum regression checks

Every push to `main`, pull request targeting `main`, and manual run starts
[Build and regress](https://github.com/huangjie2025/CosmicScanner/actions/workflows/ci.yml).
The four native jobs cover Windows x86_64 (Server 2022/UCRT64), Linux x86_64
(Ubuntu 24.04), macOS 15 Apple Silicon, and macOS 15 Intel.

## Downloading a development build

Open a successful workflow run and select the artifact for your platform.
A GitHub login may be required to download Actions artifacts. Extract the artifact
ZIP, then extract the **inner** platform archive while preserving its directory
structure. Artifacts are retained for 30 days. For a published candidate with a
permanent public download, use [V2.0 RC1](https://github.com/huangjie2025/CosmicScanner/releases/tag/v2.0.0-rc1).

CI archives are development snapshots, named
`CosmicScanner-v2.0.0-dev.<source-sha-prefix>.<run-number>-<platform>-<arch>`.
They do not replace the frozen RC1 downloads. The executable's internal version
is the CMake project version; `PACKAGE_ID.json` identifies its exact build.
Linux CI snapshots require glibc 2.39 or later (Ubuntu 24.04 baseline); this
differs from the published RC1 Linux package's glibc 2.34 baseline. macOS builds
target macOS 15 and use ad-hoc signing, without Developer ID notarization.

## Source → environment → package → checksum

Each successful job supplies:

| File | Purpose |
|---|---|
| `BUILD_ENVIRONMENT.json` | Full source commit and Git tree, workflow and run URL, runner image, compiler/CMake/Python versions, installed dependency versions, configuration flags, downloaded dependency URLs and hashes. Also included inside the package. |
| `PACKAGE_ID.json` (inside archive) | Package version, exact source/build identity, bundled runtimes, and packaged executable SHA-256. |
| `CosmicScanner-….zip.sha256` or `….tar.gz.sha256` | SHA-256 of the downloadable inner archive. |
| `SHA256SUMS` (inside archive) | SHA-256 of each packaged file other than this checksum list itself. |
| `CI_RESULT.json` | PASS only after all checks below; archive name/hash, packaged executable hash, exact source commit and run URL. |
| Regression logs and `build.log` | Build and test evidence, retained even when a job fails. |

Compare the source commit in `CI_RESULT.json`, `BUILD_ENVIRONMENT.json` and
`PACKAGE_ID.json` with the commit shown on the workflow run. Verify the inner
archive before extraction, for example:

```bash
sha256sum -c CosmicScanner-<version>-linux-x86_64.tar.gz.sha256
# macOS: shasum -a 256 -c <archive>.sha256
```

Windows: `Get-FileHash -Algorithm SHA256 <archive>.zip` and compare against the
matching `.sha256`. After extraction, the optional package self-test is
`python test/verify_native_package.py`. Python is needed for testing, not analysis.

Compression sources are hash-verified and linked statically. Windows compiler
DLLs, macOS libomp and Linux non-system runtimes are bundled with license texts.
System libc/OS libraries remain host requirements. Package managers may advance
between runs; the environment manifest records their resolved versions rather
than claiming bit-for-bit reproducibility.

## Minimum checks

1. Phenotype record-field initialization regression.
2. Existing CLI compatibility and error-path regression suite.
3. Existing review regression suite for input/output and documented analyses.
4. Four packaged workflows: additive single-trait LMM, multivariate LMM,
   repeatability LMM and GRM. The synthetic fixture has 60 individuals and 180
   variants. Checks require the expected row/sample counts, finite P values in
   range, summary output and readable run manifests.
5. Archive extraction, all packaged file hashes, source identity, then the four
   workflows again with a clean executable search path and no loader overrides.

This is a minimum software regression gate. It does not establish power,
false-positive calibration, large-cohort scalability, or validation of every
model/format combination. Model support remains as documented in
[SUPPORTED_MODES.md](SUPPORTED_MODES.md); longitudinal RRM remains experimental.
The workflow has read-only repository access and does not publish releases or
move tags. A failed job can still have a logs artifact; check both the green job
status and `CI_RESULT.json` before using a development package.
