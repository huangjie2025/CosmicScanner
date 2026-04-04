# CosmicScanner

CosmicScanner is a high-performance Genome-Wide Association Study (GWAS) engine designed for large-scale genetic evaluation. It implements an exact Mixed Linear Model (LMM) with "on-the-fly" Genetic Relationship Matrix (GRM) construction, eliminating I/O bottlenecks and offering blazing fast performance.

## Features

- **Exact LMM & High Accuracy**: Perfect alignment with industry gold standards (e.g., GEMMA). Pearson correlation coefficient of -log10(P) > 0.999.
- **On-the-fly GRM**: Dynamically constructs VanRaden GRM in memory, significantly reducing disk I/O and runtime.
- **EBI GWAS-SSF Support**: One-click output to the modern GWAS Summary Statistics Format standard (`--out-ssf`), ensuring seamless downstream integration.
- **Modern I/O Compatibility**: Fully supports PLINK (`.bed`/`.bim`/`.fam`), PGEN, and BGEN binary formats.

## Getting Started

### Prerequisites
A Linux environment. (Windows executable `cosmicscanner.exe` is also provided in releases).

### Quick Run

We provide a small test dataset in the `example/` directory.

```bash
./cosmicscanner --bfile example/test_data --pheno example/test_data.pheno --gwas --no-loco --out results
```

### Outputs
The output will be saved as `results.assoc.txt` containing the standard GWAS summary statistics (CHR, SNP, BP, A1, A2, FREQ, N, BETA, SE, P, LOG10P).

## License
MIT License
