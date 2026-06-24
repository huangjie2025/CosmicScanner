# CosmicScanner Official Tutorial

**Author**: Jie Huang
**Version**: 1.1.0
**GitHub**: [https://github.com/huangjie2025/CosmicScanner](https://github.com/huangjie2025/CosmicScanner)

---

## 1. Introduction
**CosmicScanner** is a high-performance computational tool built for Genome-Wide Association Studies (GWAS) and is a core component of the **Cosmic** pan-omics analysis suite. It shatters traditional bottlenecks in computing speed, memory consumption, and model versatility, easily handling ultra-large cohort genotype data and complex phenotypes.

**Supported Platforms**:
- **Linux**: Native support for major Linux distributions (x86_64 architecture).
- **Windows**: Recommended to run via WSL2 (Windows Subsystem for Linux).

---

## 2. Core Features
CosmicScanner currently supports the following core features:
1. **Standard Single-Trait LMM GWAS**: Utilizes exact variance component estimation and Leave-One-Chromosome-Out (LOCO) to thoroughly account for population stratification and proximal contamination.
2. **Multivariate LMM GWAS**: Directly tests multiple phenotypes jointly, significantly increasing statistical power.
3. **Repeatability Model GWAS**: Analyzes long-format data where individuals have multiple measurements, offering both record-level exact GLS scanning and mean-aggregated fast scanning.
4. **Random Regression Model (RRM) GWAS**: Employs Legendre polynomials to perform joint trajectory association tests for longitudinal phenotypes across time or environmental gradients.
5. **Rare Variant Analysis**: Integrates SKAT and Burden tests, and supports Saddle Point Approximation (SPA) for binary traits.
6. **PSR-GWAS for non-normal traits**: Supports binary, binomial event/trial, and Poisson traits using PSR approximation plus weighted LMM scanning.
7. **Biobank-Scale Cohort Analysis**: Built-in REGENIE-style Ridge Regression (Biobank Mode) to tackle cohorts with millions of samples.
8. **Ultra-Fast GRM Construction**: Optimized Bit-Parallel algorithm with zero-copy `mmap` cache reuse.

---

## 3. Installation & Compilation

### 3.1 Pre-compiled Binaries (Recommended)
You can download the latest pre-compiled binaries from the GitHub Release page:
```bash
tar -zxvf CosmicScanner-1.1.0-Linux-x86_64.tar.gz
chmod +x cosmicscanner
./cosmicscanner --help
```

### 3.2 Compiling from Source
Ensure `CMake (>= 3.15)` and a C++17 compatible compiler (e.g., GCC 9+) are installed:
```bash
git clone --recurse-submodules https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
```
The executable will be located at `build/cosmicscanner`.

---

## 4. Data Preparation

CosmicScanner requires standard format inputs:
- **Genotypes**: PLINK 1.9 (`.bed/.bim/.fam`), PLINK 2.0 (`.pgen/.pvar/.psam`), or BGEN (`.bgen/.sample`).
- **Phenotypes (`--pheno`)**: Text files with 3+ columns (FID, IID, Pheno1...), missing values coded as `NA` or `-9`. Long-format data is supported for Repeatability and RRM models.
- **Covariates (`--covar` / `--qcovar`)**:
  - `--covar`: Categorical covariates (e.g., sex, batch).
  - `--qcovar`: Quantitative covariates (e.g., age, BMI).

---

## 5. Detailed Usage & Examples

### 5.1 Standard Single-Trait LMM GWAS
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --mpheno 1   --covar data/my_covar.txt   --out results/standard_gwas
```

### 5.2 Multivariate LMM GWAS
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2,3   --multivariate   --mv-method emmv   --out results/multi_gwas
```
The default `--mv-method emmv` uses an EMmvGWAS-style semi-exact multivariate Wald scan: the null model estimates trait covariance matrices after GRM eigendecomposition, then the SNP scan reuses the null variance structure and updates the per-SNP residual covariance for the multivariate Wald test. `--mv-method wald` keeps the fixed-null Wald route for debugging and comparison. `--mv-method gemma` provides a GEMMA-style per-SNP variance-component refit Wald route; it is slower and intended for small validation runs or top-SNP refinement. `--mv-method gemma-full --mv-test wald,lrt,score` outputs Wald, LRT, and Score p-values for the GEMMA-full MVP route.

`--mv-method hybrid` is the recommended performance-oriented route. It runs a genome-wide EMmvGWAS scan, optionally audits top Stage 1 SNPs with a Score-only prepass, and refines the union of top-K, p-threshold, and Score/Wald discrepancy candidates with GEMMA-full. Without `--no-loco`, hybrid uses LOCO-aware screening and per-chromosome LOCO GEMMA-full refinement.

```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2   --multivariate   --mv-method hybrid   --mv-test wald,lrt,score   --mv-refine-top 100   --mv-refine-p 1   --mv-refine-score-prepass 500   --mv-refine-discrepancy-top 50   --mv-refine-discrepancy 1   --gwas   --out results/mv_hybrid_loco
```

Hybrid outputs `*.assoc.mv.txt`, `*.refine.snplist`, `*.refined.assoc.mv.txt`, and `*.hybrid.summary.tsv`. If `--mv-refine-score-prepass` is used, it also writes `*.score_prepass.refined.assoc.mv.txt`.

For method-compatibility validation against the EMmvGWAS R package, use:
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2   --multivariate   --mv-method emmv   --mv-compatible emmvgwas   --no-loco   --out results/emmvgwas_compatible
```
`--mv-compatible emmvgwas` automatically enables the EMmvGWAS-compatible GRM, equivalent to `--make-grm-alg 6`: it uses `Z = dosage - 1` (`-1/0/1` coding) and normalizes by `mean(diag(Z Z'))`. This mode is intended for compatibility validation against EMmvGWAS R; it is not recommended as the default replacement for VanRaden or GEMMA-style GRMs. For routine GWAS, keep the default `--make-grm-alg 1`.

`--mv-method gemma-full` now defaults to a fast GEMMA-style path (null variance components), matching GEMMA mvLMM behavior. For higher-accuracy per-SNP variance refit, use `--mv-refit-p`:
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2,3   --multivariate   --mv-method gemma-full   --mv-test wald,lrt,score   --mv-refit-p 0.05   --no-loco   --out results/mv_gemma_full_refit
```
The default `--mv-refit-p 1.0` uses the fast null-component path (no per-SNP refit). Values below 1 first compute the null Score test and only fully refit SNPs with Score p-values at or below the threshold; other SNPs keep Score-only output. For production-scale genome-wide scans, `--mv-method hybrid` is still recommended.

### 5.3 Repeatability Model GWAS
For repeated measures data (long format):
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/repeat_pheno.txt   --repeatability-model   --repeatability-grm-cache grm_cache   --out results/repeat_gwas
```
The recommended production path is `--repeatability-mode exact`, which uses a record-level repeated-measure LMM with subject-level Woodbury solves. `--repeatability-mode fast` is a mean-aggregated approximation intended for quick screening only, especially when repeated records are balanced and record-level covariates are not important. `--vce-mode auto` uses scalable subject-level Woodbury VCE; use `--vce-mode full-record-exact` only for small-data validation.

### 5.4 PSR-GWAS For Non-Normal Traits

PSR-GWAS provides a fast scan for non-normal traits. The current implementation uses PSR/PQL linearization to produce a pseudo response and working residual variance, then runs a weighted LMM SNP scan. It is not an exact GLMM. Significant loci should be followed by candidate refinement or external exact GLMM/GMMAT/POLMM validation.

Binary traits can use the probit/logit threshold model:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/binary.pheno \
  --mpheno 1 \
  --gwas \
  --threshold-model \
  --link probit \
  --out results/psr_binary
```

Poisson count traits support log-offset/exposure. The `--offset` file contains one log-offset value per sample in genotype sample order:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/count.pheno \
  --mpheno 1 \
  --gwas \
  --count-model \
  --family poisson \
  --offset data/count.log_offset \
  --out results/psr_poisson
```

Binomial proportion traits should use event/trial columns instead of precomputed proportions. The file format is `FID IID EVENT TRIAL`:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --count-model \
  --family binomial \
  --binomial-pheno data/binomial_event_trial.pheno \
  --out results/psr_binomial
```

Ordinal traits (3+ ordered categories) should use the v0.5 stacked proportional-odds PSR path for production scans:

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --ordinal-model \
  --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --out results/psr_ordinal
```

The default `.gwas.tsv` statistic is the corrected Wald / pseudo-score result. Use `--ordinal-test all` to write `.ordinal.score.tsv` diagnostics; use `--ordinal-test robust-score` only for cluster-robust sensitivity analysis.

Current limitations and recommendations:

- `--family negbin` is still experimental; overdispersion estimation has not been fully validated and should not be used for production claims yet.
- The official PSR-GLMM/R website currently times out from this development host, so the official source/demo data are not integrated. See `tests/PSR_GLMM_R_REFERENCE_STATUS.md`.
- Candidate refinement MVP is available at `tests/psr_candidate_refine.R`; it keeps R `glm` as the default and also exposes optional `firth`, `glmer`, and `glmmTMB` hooks. These optional methods require `logistf`, `lme4`, and `glmmTMB`, respectively. Missing packages are reported per SNP through `converged=FALSE` and `error_message` instead of crashing the validation suite.

PSR/PQL tuning options are available when difficult traits need more iterations or less aggressive clipping: `--psr-maxit`, `--psr-tol`, `--psr-mu-eps`, `--psr-min-weight`, and `--psr-max-weight`. Defaults are `50`, `1e-6`, `1e-6`, `1e-8`, and `1e8`.

Convergence warnings should be treated as diagnostics, not automatic failures. A run with `status_class=CONVERGED_WITH_WARNING` met the convergence criterion but had boundary variance, non-monotone quasi-likelihood, or clipping diagnostics. A run with `status_class=NOT_CONVERGED` should be treated cautiously in formal analyses. `WARNING_PSR_MAXIT_REACHED` means the iteration limit was reached, `WARNING_PSR_CLIPPED` means `mu` or working weights hit clipping bounds, `WARNING_PSR_SIGMA_BOUNDARY` means the final genetic variance is close to zero, and `WARNING_PSR_QUASILIK_DECREASED` flags a decrease in the quasi-likelihood proxy.

Run the PSR-GWAS smoke, R baseline, and candidate refinement validation:

```bash
bash tests/run_psr_gwas_validation.sh
```

The script also generates:

- `*.psr.convergence.tsv`: PSR/PQL iteration diagnostics, including `mu` range, weight range, clipping count, and convergence flag.
- `psr_convergence_summary.tsv`: one-row-per-trait convergence summary with warning flags.
- `weighted_lmm_equivalence.tsv`: small fixed-variance checks that `V = K sigma_g2 + R` is equivalent to the weighted LMM transform.
- `psr_type1_summary.tsv`, `psr_power_summary.tsv`, and `psr_qq_points.tsv`: lightweight simulation smoke summaries generated by actual CosmicScanner PSR-GWAS runs. Publication-grade validation still requires larger simulations and external GLMM/official PSR-GLMM/R baselines.
- `method_metrics.tsv`: per-run simulation metrics with `engine=cosmicscanner_psr`, seed, trait type, scenario, runtime, convergence flag, and `status_class`.
- `psr_validation_report.md`: Markdown report combining convergence, weighted equivalence, simulation, candidate refinement, and known limitations.

See `docs/PSR_GWAS_ROADMAP.md` for the PSR-GWAS roadmap, `docs/PSR_GWAS_V0.2_VALIDATION_SUMMARY.md` for the v0.2 validation record, and `docs/PSR_GWAS_V0.3_STABLE_SUMMARY.md` for the stable sprint summary.

Include PSR-GWAS in the core regression suite:

```bash
bash tests/run_gwas_core_regression.sh --fast --with-psr
```

### 5.5 Random Regression Model (RRM) GWAS
For longitudinal/trajectory phenotypes:
```bash
./cosmicscanner   --bfile data/my_geno   --long-pheno data/long_pheno.txt   --mpheno 2   --tpheno 1   --rrm-order 2   --gwas   --out results/rrm_gwas
```

### 5.6 Rare Variant Analysis (SKAT / Burden)
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --set-file data/gene_sets.txt   --skat   --out results/skat_gwas
```

### 5.7 Biobank-Scale Analysis (REGENIE Mode)
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --biobank-mode   --block-size 1000   --out results/biobank_gwas
```

---

## 6. Parameter Reference

| Category | Argument | Description |
| --- | --- | --- |
| **I/O** | `--bfile <prefix>` | PLINK 1.9 binary prefix |
| | `--pfile <prefix>` | PLINK 2.0 binary prefix |
| | `--bgen <file>` | BGEN format genotype file |
| | `--out <prefix>` | Output prefix (default: cosmic_gwas) |
| **Pheno/Covar**| `--pheno <file>` | Phenotype file path |
| | `--mpheno <n>` | Target phenotype column index |
| | `--mpheno-list <l>`| Phenotype indices for multivariate (e.g. 1,2,3) |
| | `--covar <file>` | Categorical covariate file |
| | `--qcovar <file>` | Quantitative covariate file |
| **Models** | `--multivariate` | Enable Multivariate LMM |
| | `--mv-method <name>` | Multivariate scan method: `emmv`, `wald`, `gemma`, `gemma-full`, or `hybrid`; default `emmv` |
| | `--mv-compatible <name>` | Multivariate compatibility mode; `emmvgwas` enables the EMmvGWAS-compatible GRM for R-package benchmarking |
| | `--mv-test <list>` | GEMMA-full test list: `wald,lrt,score`; all enabled by default |
| | `--mv-refine-top <n>` | Hybrid/GEMMA refinement SNP limit; default 100 |
| | `--mv-refine-p <p>` | Hybrid/GEMMA refinement p-value threshold; default 1e-3 |
| | `--mv-refit-p <p>` | GEMMA-full full-refit threshold for fast mode; default 1 means full refit for all SNPs |
| | `--threshold-model` | Binary PSR/probit or PSR/logit GWAS |
| | `--link probit\|logit` | Link function for threshold model; default `probit` |
| | `--count-model` | Count/binomial PSR-GWAS entry point |
| | `--family poisson\|binomial\|negbin` | Count model family; `negbin` is still experimental |
| | `--offset <file>` | Poisson log-offset/exposure, one value per sample |
| | `--binomial-pheno <file>` | Binomial event/trial phenotype file: `FID IID EVENT TRIAL` |
| | `--ordinal-model` | Enable ordinal PSR-GWAS |
| | `--ordinal-pheno <file>` | Ordinal phenotype file: `FID IID ORDINAL_LEVEL` |
| | `--ordinal-mode <name>` | Ordinal algorithm: `experimental-cumulative-binary` or `stacked-psr` |
| | `--ordinal-stacked-engine <name>` | Stacked-PSR engine: `fixed-gls` or `dense-zkz` |
| | `--ordinal-vc-method <name>` | Ordinal VC estimation: `fixed` or `grid-reml` |
| | `--ordinal-test <name>` | Stacked-PSR association test / diagnostics: `wald`, `score`, `robust-score`, or `all` |
| | `--psr-maxit <n>` | PSR/PQL maximum iterations; default 50 |
| | `--psr-tol <f>` | PSR/PQL convergence tolerance; default 1e-6 |
| | `--psr-mu-eps <f>` | PSR mean clipping epsilon; default 1e-6 |
| | `--psr-min-weight <f>` | Minimum PSR working weight; default 1e-8 |
| | `--psr-max-weight <f>` | Maximum PSR working weight; default 1e8 |
| | `--repeatability-model` | Enable Repeatability Model for repeated measures |
| | `--long-pheno <file>` | Enable RRM and set the long-format phenotype file |
| | `--rrm-order <n>` | Set the Legendre polynomial order for RRM; requires `--long-pheno` |
| | `--fast-lmm` | Enable approximate LMM scanning |
| | `--spa` | Enable Saddle Point Approximation for binary traits |
| | `--skat` / `--burden`| Enable rare variant aggregation tests |
| | `--biobank-mode` | Enable REGENIE-style ridge regression mode |
| **Optimization**| `--make-grm` | Calculate GRM and exit |
| | `--make-grm-alg 6` | EMmvGWAS-compatible GRM: `Z=dosage-1`, normalized by `mean(diag(Z Z'))`; intended for compatibility validation only |
| | `--repeatability-mixed-precision`| Use Float32 to reduce memory footprint |
| | `--repeatability-grm-cache <pre>`| Enable `mmap` cache for extremely fast reruns |
