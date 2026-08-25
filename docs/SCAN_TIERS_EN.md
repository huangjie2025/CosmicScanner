# Gaussian Mixed-Model GWAS Scan Tiers

`--scan-tier approx|semi|exact` controls per-SNP nuisance variance handling for single-trait, multivariate, and repeatability Gaussian LMMs.

| Tier | Variance handling |
|---|---|
| `approx` | Reuse all null-model variance parameters |
| `semi` | Fix variance ratios/covariance shape and profile a SNP-specific residual scale/covariance |
| `exact` | Refit alternative variance parameters per SNP |

Canonical commands explicitly include `--gwas`, phenotype columns, the tier, and an output prefix:

```bash
# Single trait
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 \
  --gwas --scan-tier semi --out results/trait

# Multivariate
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier semi --out results/mv

# Repeatability
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas \
  --scan-tier semi --out results/repeatability
```

Multivariate and repeatability `exact` currently require `--no-loco`; repeatability `exact` is additionally limited to PLINK input without record weights. `--mv-method` remains an advanced/compatibility interface and must not conflict with `--scan-tier`.
