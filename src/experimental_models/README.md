# Experimental Models

These models are **not part of CosmicScanner v1.0** official scope.

## Models Included

| Model | Files | Description |
|-------|-------|-------------|
| Ordinal PSR-GWAS | `ordinal_model.cpp/h` | Cumulative-binary & stacked PSR for ordinal traits |
| Threshold/Probit | `threshold_model.cpp/h` | Binary trait GWAS (PQL + probit/logit) |
| Count/Binomial PSR-GWAS | `count_model.cpp/h` | Poisson/Negative Binomial/Binomial traits |
| Survival/Cox | `survival_model.cpp/h` | Time-to-event GWAS (Cox + martingale residual) |
| GxE Reaction Norm | `gxe_model.cpp/h` | SNP-by-environment interaction |
| BSLMM Compatibility | `gemma_bslmm_model.cpp/h` | GEMMA-style BSLMM |
| GEMMA VC Compatibility | `gemma_vc_model.cpp/h` | GEMMA-style variance components |
| Biobank REGENIE 2-step | `gwas_biobank_kernels.cpp`, `gwas_biobank_scan.cpp` | UK Biobank scale 2-step |
| Rare Variant Test | `rare_variant_test.cpp/h` | SKAT / Burden tests |

## Build

Experimental models are **disabled by default** for v1.0 release.

To enable:
```bash
cmake -DCOSMICSCANNER_ENABLE_EXPERIMENTAL=ON ..
```

## Status

- Source code preserved for research continuity
- Not tested in v1.0 release validation
- Not documented in official Tutorial/Release Notes
- May have incomplete features or known issues
- API may change in future versions

## Warning

**Do not use for production analysis.** These models are research prototypes.