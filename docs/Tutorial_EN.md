# CosmicScanner Official Tutorials

CosmicScanner is a high-performance Genome-Wide Association Study (GWAS) scanning engine designed for **Biobank-scale** datasets within the Cosmic ecosystem. It not only benchmarks against the classic software GEMMA in computational precision but also natively supports modern high-compression genomic formats (PLINK 2.0 PGEN, BGEN) at the architectural level, easily handling association analyses for millions of samples.

This tutorial will detail CosmicScanner's input formats, data filtering, Genomic Relationship Matrix (GRM) construction, mixed-model association analyses (Exact LMM/Fast LMM), and multivariate phenotype analyses.

---

## 1. Genomic Data Input and Filtering

CosmicScanner's data I/O module is highly optimized to rapidly parse the following three mainstream formats without the need for intermediate conversions:

1.  **PLINK 1.9 Binary Format (`.bed`, `.bim`, `.fam`)**:
    Specify the file prefix using `--bfile <prefix>`.
2.  **PLINK 2.0 High-Compression Format (`.pgen`, `.pvar`, `.psam`)**:
    Specify the file prefix using `--pfile <prefix>` (recommended for datasets with hundreds of thousands of samples).
3.  **BGEN Format (Layout 2 zlib compressed)**:
    Specify using `--bgen <file>` along with `--sample <file>`.

### Quality Control (Data Filtering)
Before analysis, you can use CosmicScanner's built-in QC parameters to filter variants on the fly:
*   `--maf <float>`: Minimum Allele Frequency filter (default: 0.01).
*   `--hwe <float>`: Hardy-Weinberg Equilibrium P-value filter (default: 0.0, i.e., disabled).

```bash
# Example: Read BGEN data and filter low-frequency variants
./cosmic_scanner --bgen ukb_chr22.bgen \
                 --sample ukb_chr22.sample \
                 --maf 0.05 \
                 --hwe 1e-6
```

---

## 2. Genomic Relationship Matrix (GRM) Construction

Before performing Mixed Linear Model (LMM) analysis, a Genomic Relationship Matrix needs to be constructed to correct for population stratification and cryptic relatedness.

CosmicScanner provides an independent GRM construction mode:
```bash
./cosmic_scanner --bfile demo \
                 --make-grm \
                 --make-grm-alg 1 \
                 --out demo_grm
```

| Algorithm (`--make-grm-alg`) | Description |
| :--- | :--- |
| **0 (Yang)** | GCTA's default GRM construction algorithm. |
| **1 (VanRaden)** | Default additive model coding, more suitable for agricultural breeding and routine association analysis. |

> **Advanced Features**: If you are analyzing admixed populations or need a weighted GRM, you can also provide `--pop-class <file>` or `--snp-weight <file>` to construct a kinship matrix with specific structures.

---

## 3. Genome-Wide Association Scan (GWAS Scan)

CosmicScanner offers various modes to accommodate different sample sizes and precision requirements.

### Standard Exact Mixed Model (Exact LMM)
The exact model re-evaluates the fixed effects at every tested SNP locus while correcting for covariates and the polygenic background effects from the GRM. Its precision is completely consistent with GEMMA.

```bash
./cosmic_scanner --bfile demo \
                 --pheno phenotype.txt \
                 --mpheno 1 \
                 --covar discrete_covar.txt \
                 --qcovar quantitative_covar.txt \
                 --gwas \
                 --out my_exact_results
```
*   `--pheno`: Phenotype file (`FID` `IID` `Pheno`...).
*   `--mpheno <n>`: Specifies the column number in the phenotype file.
*   `--covar` / `--qcovar`: Discrete (categorical) covariate and quantitative (continuous) covariate files.

### Saddle Point Approximation (SPA) for Imbalanced Traits
In cases where Case-Control data is extremely imbalanced (e.g., prevalence is only 1%), traditional LMMs produce severe P-value inflation (false positives). By adding the `--spa` flag, CosmicScanner will automatically use Saddle Point Approximation to perform an ultra-high precision P-value recalibration on marginally significant loci, thoroughly resolving false positives (benchmarking against SAIGE).

### Accelerated Approximate Mixed Model (Fast LMM)
For Biobank data with hundreds of thousands of samples, Exact LMM can be extremely time-consuming. In such cases, `--fast-lmm` is recommended. In this mode, the program first fits a null model containing only covariates, fixes the ratio of genetic variance to residual variance ($\delta = \sigma_e^2 / \sigma_g^2$), and then reuses this variance ratio during the genome-wide scan, increasing speed by tens of times.

### Biobank-Scale Ridge Regression Mode (REGENIE-style)
When sample sizes reach millions, even calculating a GRM becomes a memory catastrophe. You can activate `--biobank-mode`.
In this mode, CosmicScanner completely adopts a Block-wise Cross-Validated Ridge Regression algorithm similar to **REGENIE**, bypassing the $N \times N$ kinship matrix to directly compute the polygenic background of each individual and generate Leave-One-Chromosome-Out (LOCO) predictions.

```bash
./cosmic_scanner --pfile ukb_data \
                 --pheno ukb_pheno.txt \
                 --gwas \
                 --biobank-mode \
                 --block-size 1000 \
                 --out my_biobank_results
```

### Multivariate Phenotype Joint Analysis (Multivariate GWAS)
When facing complex diseases or highly correlated physiological indicators, CosmicScanner supports direct multivariate joint GWAS scans.
```bash
./cosmic_scanner --bfile demo \
                 --pheno multi_pheno.txt \
                 --mpheno-list 1,2,3 \
                 --multivariate \
                 --gwas \
                 --out my_mv_results
```

---

## 4. Rare Variant Association Test

When the Minor Allele Frequency (MAF) of variants is extremely low, single-point association tests suffer from severely diminished power. CosmicScanner supports region-based (e.g., gene-level) aggregation tests for rare variants, directly benchmarking against SAIGE-GENE and the SKAT R package.

You need to provide a `.set` formatted gene region file (two columns: `SetID` and `SNP_ID`) and use it in conjunction with the `--skat` or `--burden` flags.

```bash
./cosmic_scanner --bfile demo \
                 --pheno phenotype.txt \
                 --set-file genes.set \
                 --skat \
                 --out my_rare_variants
```
*   `--skat`: Sequence Kernel Association Test. Assumes that different variants within the same gene may have opposite effect directions.
*   `--burden`: Burden Test. Assumes that all rare variants within the same gene have the same effect direction.

---

## 5. Output Result Format

By default, CosmicScanner outputs a tab-separated file (`.assoc.txt`) similar to GEMMA/GCTA, containing:
`CHR`, `SNP`, `BP`, `A1`, `A2`, `FREQ`, `N`, `BETA`, `SE`, `P`, `LOG10P`.

**EBI GWAS-SSF Format Support**:
To facilitate direct submission of your research findings to international public databases (such as the NHGRI-EBI GWAS Catalog), you can add the `--out-ssf` parameter. The program will automatically output association summary statistics files that comply with the official GWAS-SSF standards.
