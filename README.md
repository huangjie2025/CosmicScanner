# CosmicScanner Source Package

CosmicScanner v1.1.0 is a focused, reproducible source release for the first
manuscript. It contains four supported functions:

- Single-trait GWAS
- Multivariate GWAS
- Repeatability GWAS
- A-inverse, G, and H-inverse construction

All retained C++ files are placed directly under a single flat `src/` directory.
The command-line interface exposes only the four supported functions. Shared
numerical implementation files required by those workflows remain internal and
do not define additional supported commands.

## Quick Start

```bash
bash scripts/build.sh
bash scripts/validate_release.sh
```

The release defaults to Eigen's native solvers. BLAS/LAPACK are opt-in because
results depend on the selected vendor runtime:

```bash
COSMICSCANNER_ENABLE_BLAS=ON bash scripts/build.sh
```

## Documentation

- `docs/USER_MANUAL_EN.md`: comprehensive English user manual.
- `docs/USER_MANUAL_CN.md`: comprehensive Chinese user manual.
- `docs/TROUBLESHOOTING.md`: solutions for common build, runtime, and output issues.

## Dependencies

- Linux x86_64
- CMake 3.15 or newer
- C++17 compiler
- Eigen3
- OpenMP
- zlib and zstd

If CMake cannot find Eigen3, set `CMAKE_PREFIX_PATH` to the prefix containing
`share/eigen3/cmake/Eigen3Config.cmake`.

## Citation

Citation metadata are provided in `CITATION.cff`.
