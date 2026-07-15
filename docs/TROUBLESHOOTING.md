# Troubleshooting Guide

## Build Issues

### CMake cannot find Eigen3
Set CMAKE_PREFIX_PATH to the Eigen3 prefix.

### Compiler does not support C++17
Use GCC 8+ or Clang 10+. Set CXX=g++-10 or equivalent.

### OpenMP not found
Install libomp-dev (Ubuntu) or libomp-devel (CentOS).

## Runtime Errors

### "This focused build requires --gwas or --make-agh"
Add --gwas or --make-agh to the command line.

### "Unsupported option in focused source build"
The parameter is not in v1.1.2 scope. Run cosmicscanner --help for supported
options. RRM, threshold, survival, GxE, SKAT, BSLMM, biobank mode are excluded.

### "Missing value for --<option>"
Provide the required argument.

### Genotype file errors
PLINK BED must be SNP-major (0x01). BGEN requires --sample. PGEN requires
.pgen/.pvar/.psam.

### Phenotype file errors
Minimum 3 columns (FID, IID, trait). Check --mpheno/--tpheno column numbers.

### Repeatability exact tier limitations
Requires --no-loco, PLINK input, and no --repeatability-weight.

### Multivariate exact requires --no-loco
Add --no-loco or use --scan-tier semi.

### GRM fails: "No valid SNPs found"
Loosen MAF/HWE thresholds. Check extract/exclude files.

### VCE does not converge
Increase --reml-maxit or loosen --reml-tol.

## Output Issues

### Empty output
All SNPs filtered out; check MAF/HWE/extract/exclude/chr settings.

### NaN p-values
Small sample size, monomorphic SNP, or collinearity with covariates.

## Performance

### Slow runtime
Increase --threads, use --scan-tier approx, reduce SNPs with --extract/--chr.

### Out of memory
Reduce samples/SNPs, use --chr to limit, avoid --scan-tier exact on large data.

## Platform

### Linux binary GLIBC error
Build from source on the target system.

### Windows/macOS
Use the matching v1.1.2 GitHub Release archive. Windows packages target x86_64
UCRT64. On macOS, approve the downloaded executable in System Settings if
Gatekeeper blocks an ad-hoc-signed archive. Source builds require Homebrew
`eigen`, `libomp`, and `zstd` on macOS or MSYS2 UCRT64 dependencies on Windows.

## Getting Help

Open an issue at https://github.com/huangjie2025/CosmicScanner/issues with the
full command line, error output, and system information.
