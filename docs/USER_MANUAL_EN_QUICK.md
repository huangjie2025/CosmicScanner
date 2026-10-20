# CosmicScanner v2.0.0 Quick Start

This concise guide covers the bundled example. For complete documentation, read USER_MANUAL_CN.md (Chinese) or USER_MANUAL_EN_FULL.md (English).
Choose Windows x86_64, Linux x86_64, macOS arm64 (Apple Silicon), or macOS x86_64 (Intel) as appropriate. macOS packages require macOS 15 or later. Keep the complete extracted folder, including all runtime libraries.

From the extracted package directory on macOS/Linux:
```sh
./bin/cosmicscanner --version
./bin/cosmicscanner --help
mkdir -p demo_results
./bin/cosmicscanner gwas --model lmm --bfile examples/dataset/synth --pheno examples/dataset/quantitative.pheno --mpheno 1 --threads 1 --out demo_results/lmm
```
The association file demo_results/lmm.assoc.txt contains a header and 180 SNP rows. Use the named BETA_A, SE_A and P_A columns. Trait indices start after FID/IID.

On Windows, use .\bin\cosmicscanner.exe instead. To run the optional four-workflow package self-test, use Python 3: `python3 test/verify_native_package.py` (Windows: `python test/verify_native_package.py`). Python is not required for analysis.

See NATIVE_INSTALL_CN.md for platform installation details and USER_MANUAL_CN.md for model-specific commands and output interpretation. The complete CLI reference is available with `--help-all`.

Verify the archive SHA-256 and retain PACKAGE_ID.json. Passing these package checks does not establish scientific validation of every model. macOS packages are ad-hoc signed, without Apple Developer ID notarization. Earlier macOS releases have not been verified.
