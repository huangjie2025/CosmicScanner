# PSR-GLMM/R Reference Status

Official URL reported by the paper:

`http://106.13.144.206/psrglmm/index.html`

Current access status from this development host: timeout with zero bytes received using both header and page requests.

Commands attempted:

```bash
curl -I --max-time 20 "http://106.13.144.206/psrglmm/index.html"
curl -L --max-time 30 "http://106.13.144.206/psrglmm/index.html"
```

License status: not confirmed. Until an explicit license is available, PSR-GLMM/R should be treated as a read-only algorithm reference and benchmark target. Do not copy source code into CosmicScanner.

Validation hook:

`tests/run_psr_gwas_validation.sh` supports `PSR_GLMM_R_SCRIPT=/path/to/adapter.R`. The local `tests/psr_glmm_r_adapter.R` currently provides a deterministic R GLM baseline and a placeholder for sourcing the official PSR-GLMM/R code once its API and license are available.
