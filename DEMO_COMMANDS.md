# CosmicScanner Demo Commands

Below are some common commands to get you started with CosmicScanner.

## 1. Basic Exact LMM GWAS (PLINK format)
This command runs a standard Linear Mixed Model using `test_data.bed/.bim/.fam` and outputs standard summary statistics.
```bash
./cosmicscanner --bfile example/test_data --pheno example/test_data.pheno --gwas --out results/basic_run
```

## 2. Disable LOCO (Global GRM)
By default, CosmicScanner uses the Leave-One-Chromosome-Out (LOCO) strategy to avoid proximal contamination. You can disable this to use a single global GRM for all chromosomes.
```bash
./cosmicscanner --bfile example/test_data --pheno example/test_data.pheno --gwas --no-loco --out results/no_loco_run
```

## 3. EBI GWAS-SSF Output Format
To generate results compatible with the EBI NHGRI GWAS Catalog Summary Statistics Format (SSF), append `--out-ssf`.
```bash
./cosmicscanner --bfile example/test_data --pheno example/test_data.pheno --gwas --out-ssf --out results/ssf_run
```

## 4. BGEN Format (Dosage Data)
CosmicScanner fully supports BGEN files (often used in Biobank imputations).
```bash
./cosmicscanner --bgen data.bgen --sample data.sample --pheno data.pheno --gwas --out results/bgen_run
```

## 5. PGEN Format (PLINK 2.0)
For ultra-large datasets, PGEN format is recommended.
```bash
./cosmicscanner --pgen data.pgen --pvar data.pvar --psam data.psam --pheno data.pheno --gwas --out results/pgen_run
```
