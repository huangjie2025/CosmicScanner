# CosmicScanner v2.0.0 User Manual

Software for genome wide association analysis

Developer: Huang Jie

Contact: huangjie2025@zju.edu.cn

Website: see the accompanying CosmicScanner website

Platforms: Windows / Linux / macOS

macOS packages: Apple Silicon arm64 and Intel x86_64; macOS 15 or later

This manual is for researchers using CosmicScanner for genetic breeding and association analysis. It explains how to obtain and install the software, prepare data, select models, run analyses and check results. Complete the bundled examples in Chapter 4 first, then use Chapter 5 to select the modules appropriate for your study design.

This manual covers CosmicScanner v2.0.0. The current distribution is **v2.0.0-rc1**, a release candidate. Command examples use Linux or macOS terminals and Windows PowerShell. Use the package matching your operating system and preserve its directory structure. Package checks and successful introductory examples do not establish scientific validation of every model or compatibility with every system configuration.

This complete English edition corresponds to the Chinese manual supplied on 15 October 2026. The original short English quick start remains available separately. Command names, option values, filenames and the program's English help output are retained; explanatory prose and comments are translated.

## Contents

Chapter 1 Software overview
Chapter 2 Installation and execution
Chapter 3 Input file formats
Chapter 4 Quick start
Chapter 5 Analysis modules
Chapter 6 Website documentation and viewing results
Chapter 7 Output files
Chapter 8 Troubleshooting and reproducibility
Appendix A Complete command line help
Appendix B Common command line options
Appendix C Command migration and analysis records

## Chapter 1 Software overview

### 1.1 Introduction

CosmicScanner is command-line software for genome-wide association analysis and relationship-matrix computation. A typical workflow reads genotypes, phenotypes and covariates, matches samples, filters variants, fits the selected model and writes variant effects, standard errors, test statistics and analysis records. Computation runs locally; research data remain on the user's computer or server.

The software provides entry points for single-trait, multi-trait and repeated-record analyses. Select a model according to phenotype type, repeated measurements, population structure and relatedness before choosing a computational strategy. Faster computation, a more complex model or the word `exact` in an option name does not by itself establish that a result is more appropriate for the research question.

### 1.2 Main capabilities

| Capability | Main use | Usage considerations |
|---|---|---|
| Linear and mixed models | Single-marker association for continuous traits | Complete the BED example first |
| Additive and dominance analysis | Separate additive, dominance and joint effects | Check coding, input format and random-effect structure together |
| Multi-trait analysis | Joint analysis of correlated traits | Specify phenotype columns and scan strategy |
| Repeatability and random regression | Repeated measurements and changes over time | Distinguish records, individuals and time |
| Relationship matrices and variance components | Construct GRMs and fit variance structures | Match matrix entries to sample IDs |
| Pedigree and single-step analysis | Combine pedigree and genotype information | Different Methods estimate different quantities |
| Binary and other distribution models | Binary, count and ordinal traits | Choose the link for the phenotype distribution and inspect fitting |
| Regional and post-association analysis | Set tests, conditional analysis, COJO, meta-analysis and fine-mapping | Check alleles, LD and statistical assumptions |
| Input adaptation and conversion | Organize tables and multiple genotype formats | Reader support does not imply support in every model |

### 1.3 Choosing an analysis module

Select a model from the phenotype and study design. Start with LM or LMM for a continuous trait with one record per individual. Use a multi-trait model for several correlated continuous traits and a repeatability model for repeated measurements of the same individual. Consider random regression when the objective is to describe genetic effects that vary over time.

| Data and research question | Suggested module | Check before analysis |
|---|---|---|
| One continuous trait | LM or LMM | Whether relatedness and population structure require adjustment |
| Several continuous traits | Multi-trait model | Trait columns, missingness and covariance fitting |
| Repeated measurements per individual | Repeatability model | Individual IDs, time and record-level covariates |
| Genetic effects varying over time | RRM | Time units, measurement coverage and basis order |
| Binary trait coded 0/1 | Logistic GLM | Case coding, number of cases and variant MAC |
| Joint use of pedigree and genotypes | Single-step analysis | Pedigree matching, Method and variance treatment |
| Follow-up of summary statistics | COJO, meta-analysis or fine-mapping | Alleles, statistical fields and LD source |

Choose low-rank, matrix-free, two-stage prediction and extended-distribution models only after completing the basic workflow, using the conditions in Chapter 5. For approximate algorithms, compare results with the baseline model on a computationally manageable subset and inspect fitting, numerical stability and fallback records.

### 1.4 Platforms and computing resources

| Platform | Current availability | Requirements and notes |
|---|---|---|
| Linux x86_64 | Platform package available | glibc 2.34 or later; retain bundled libraries and launcher |
| Windows x86_64 | Platform package available | Windows 10 1903 or later recommended for Chinese paths; retain all DLLs |
| macOS Apple Silicon | arm64 package | macOS 15 or later; retain `bin` and `lib` |
| macOS Intel | x86_64 package | macOS 15 or later; retain `bin` and `lib` |

Storage for a dense relationship matrix grows with the square of the sample count. Multiple matrices, factorizations and multi-trait structures require additional memory. Requirements depend on sample size, trait count and model. Run a small pilot, observe memory use and elapsed time, then choose the scale of the full analysis. `--memory-gib` is a diagnostic budget setting, not a hard limit guaranteeing that the process will stay below that amount.

## Chapter 2 Installation and execution

### 2.1 Selecting a package and verifying files

Choose the package for your operating system and processor from the table below and download its matching `.sha256` file. All four packages contain the program, runtime libraries, example data and instructions. Analysis itself does not require Python, Conda, Homebrew or a compiler. Only the optional self-test script requires Python 3.

| Platform | Package filename |
|---|---|
| Linux x86_64 | CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz |
| Windows x86_64 | CosmicScanner-v2.0.0-rc1-windows-x86_64.zip |
| macOS Apple Silicon | CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz |
| macOS Intel | CosmicScanner-v2.0.0-rc1-macos-x86_64.tar.gz |

Linux requires glibc 2.34 or later. A 64-bit Windows 10 1903 or later system is recommended. macOS requires version 15 or later. This release set provides no native Windows ARM64, Linux ARM64 or 32-bit packages. Apple Silicon users should choose arm64; the Intel package is unnecessary for those machines. Archive names and the distribution identity use `2.0.0-rc1`. The unchanged executable reports `2.0.0` when invoked with `--version`.

On Linux, verify from the directory containing the download:

```bash
sha256sum -c CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz.sha256
```

On macOS, use the system's `shasum`. This example is for Apple Silicon; Intel users replace `arm64` with `x86_64` in the filename:

```bash
shasum -a 256 -c CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz.sha256
```

`OK` means that the downloaded contents match the checksum file. On Windows, calculate the hash in PowerShell and compare it with the first field in the matching `.sha256` file:

```powershell
Get-FileHash -Algorithm SHA256 -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip
Get-Content -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip.sha256
```

Hash comparison is case-insensitive. If the values differ, obtain the archive and its matching checksum file again. Save `PACKAGE_ID.json` and the archive hash to identify the build actually used.

### 2.2 Linux installation and first run

Extract into your chosen software directory and enter the extracted package root:

```bash
tar -xzf CosmicScanner-v2.0.0-rc1-linux-x86_64.tar.gz
cd CosmicScanner-v2.0.0-rc1-linux-x86_64
./bin/cosmicscanner --version
./bin/cosmicscanner --help
```

Preserve the package directories, including `bin`, `lib`, `libexec` and the launcher scripts. Do not move only the main executable. The launcher loads the bundled dependencies. If glibc is too old, use a system meeting the requirements or compile the matching source according to its build instructions; do not replace core system libraries.

### 2.3 Windows installation and first run

Save and extract the ZIP into a directory of your choice, such as a software directory on drive E. Do not run the executable from a ZIP preview. Run these commands in PowerShell from the directory containing the ZIP:

```powershell
Expand-Archive -LiteralPath CosmicScanner-v2.0.0-rc1-windows-x86_64.zip -DestinationPath .
Set-Location .\CosmicScanner-v2.0.0-rc1-windows-x86_64
& '.\bin\cosmicscanner.exe' --version
& '.\bin\cosmicscanner.exe' --help
```

Keep the DLLs in `bin` with the executable. In PowerShell, use `&` to invoke a program path stored in a variable or enclosed in quotes, and quote paths containing spaces. Software and data can reside on your chosen drive; copying them to drive C is unnecessary. If the destination already exists, preserve the previous version and its records, and choose a separate directory for the new package.

### 2.4 macOS installation and first run

Check the chip type in About This Mac. Apple M-series computers use arm64; Intel computers use x86_64. This example is for Apple Silicon. Intel users replace both occurrences of `arm64` with `x86_64`.

```bash
tar -xzf CosmicScanner-v2.0.0-rc1-macos-arm64.tar.gz
cd CosmicScanner-v2.0.0-rc1-macos-arm64
./bin/cosmicscanner --version
./bin/cosmicscanner --help
```

Keep `bin` and `lib` intact; the program loads the bundled `libomp.dylib`. For a `bad CPU type` error, check that the archive matches the processor architecture. For `Permission denied`, confirm that the package was fully extracted to a local directory that permits executable files and that `bin/cosmicscanner` retains execute permission.

The macOS packages are ad-hoc signed and have not been signed with Apple Developer ID or notarized. After a browser download, macOS may display a developer-verification prompt. Verify the download source and SHA256, then consult Privacy & Security in System Settings and [Apple's guidance](https://support.apple.com/en-ie/102445). A successful native run in a cloud environment does not replace this post-download system check.

macOS can use the full terminal commands in Section 4.2, with the same analysis options as Linux.

### 2.5 Running the bundled self test

After a first installation, upgrade or change of computer, run the four examples in Chapter 4. If Python 3 is installed, you can instead run the bundled self-test from the package root. The script uses only the Python standard library.

Linux and macOS:

```bash
python3 test/verify_native_package.py
```

Windows:

```powershell
python test/verify_native_package.py
```

On success, the final record reports `status` as `PASS`, 180 variants for each of `lmm`, `mv` and `rep`, and `grm_samples` as 60. It also gives the self-test output directory. Outputs are placed in the system's temporary directory and do not overwrite your research results. Without Python, use the Chapter 4 commands directly; regular analysis does not require installing Python.

If the self-test fails, check the system version, processor architecture, file integrity and package directory. Save the full command and error messages and include the version information from `PACKAGE_ID.json` when contacting the developer. The self-test checks bundled synthetic examples. Your own data still require checks of ID matching, missingness, covariates, fitting and output counts.

### 2.6 Command organization and help

The recommended order is task, model, input, model-specific options, resources and output. `gwas` and `grm` make the task explicit. Existing flat commands remain supported.

```text
cosmicscanner gwas --model lmm --bfile PREFIX --pheno FILE --out PREFIX
cosmicscanner grm --bfile PREFIX --out PREFIX
cosmicscanner help models
cosmicscanner help input
cosmicscanner --help-all
```

In your terminal, replace `cosmicscanner` with the executable path for your platform. `--model` selects a model family; `--gwas-model` selects the effect structure within a single-trait model; `--scan-tier` selects the treatment of per-variant variance. Appendix A contains full help, and Appendix B indexes commonly used options.

### 2.7 Building from source

Source builds are intended for developers and users requiring custom builds. Obtain the matching `libCosmicCore`, third-party dependencies and full build instructions. Downloading only the main directory or mixing Core versions may cause build failures or behavioral differences.

Source builds commonly use `./build/cosmicscanner`, with examples in `examples/data/r33_synthetic/dataset`. Platform packages use `./bin/cosmicscanner` and `examples/dataset`. Choose paths appropriate to how the software was obtained, and run the Chapter 4 examples after building.

## Chapter 3 Input file formats

### 3.1 Genotype files

PLINK BED is recommended for a first analysis. `PREFIX.bed`, `PREFIX.bim` and `PREFIX.fam` must come from the same dataset, with matching sample and variant orders. Supply only the prefix to `--bfile`, without the `.bed` suffix.

| Input | Main options | Files and limitations |
|---|---|---|
| PLINK BED | `--bfile PREFIX` | BED, BIM and FAM; recommended introductory format |
| PLINK PGEN | `--pfile PREFIX` | PGEN, PVAR and PSAM; check support separately for each model |
| BGEN | `--bgen FILE --sample FILE` | Implemented Layout 2 zlib scope; probability and hard-call paths differ |
| VCF | `--vcf FILE` | Compressed indexed reading requires a matching index; see `--vcf-tbi` |
| HapMap | `--hmp FILE` | The specified 11 metadata columns; check allele coding |
| Text matrix | `--geno-matrix FILE` | Also supply `--geno-var` and `--geno-sample`; specify matrix orientation |
| NPY | `--npy FILE` | Little-endian float32, with variant and sample information |
| PED MAP | `--ped-genotypes FILE --map FILE` | This PED contains genotypes and is not a pedigree file |

BGEN dosage, hard calls and non-additive coding follow different data paths. Support for an additive model does not establish support for dominance, multi-trait, record-level or exact-scan models using that format.

### 3.2 Single trait phenotypes

Use space- or tab-delimited plain text. The first two columns are FID and IID; phenotype fields follow. This example uses `--pheno quantitative.pheno --mpheno 1`.

```text
S001 S001 0.246761
S002 S002 1.519086
S003 S003 0.831204
```

`--mpheno` is a 1-based index after FID/IID, not an absolute column number in the file. In `FID IID HEIGHT WEIGHT`, HEIGHT is 1 and WEIGHT is 2. Do not copy another program's phenotype-column numbering without checking.

Use missing-value conventions supported by the selected reader, such as `NA` or the default `-9`, and inspect the retained sample count. Zero is meaningful for binary and count traits and cannot serve as a universal missing-value placeholder. Never assume phenotypes match genotypes solely because they occupy corresponding row positions.

### 3.3 Multi trait phenotypes

Place several continuous traits after the two IDs. For `FID IID HEIGHT WEIGHT LENGTH`, use `--mpheno-list 1,2,3`. Check trait units, transformations and missingness before a full analysis. Highly collinear, duplicate or zero-variance traits can cause variance-estimation or factorization failures.

The multi-trait joint P value tests a joint hypothesis. A trait-specific BETA and SE describe the effect for that trait. A significant joint P value is not evidence that every trait is individually significant.

### 3.4 Repeated measurements and time columns

A long table can contain several rows for the same FID/IID, with separate time and measurement fields. This input uses `--tpheno 1 --mpheno 2`.

```text
FID IID TIME VALUE
S001 S001 1 1.947087
S001 S001 2 1.463607
S002 S002 1 0.723190
```

Multiple records for one individual do not increase the number of genotyped individuals. Check record count, independent individual count, each individual's time range and missingness. Do not interchange time labels, actual age and batch without explaining the design. RRM imposes additional requirements on the time scale and design; see Section 5.11.

### 3.5 Covariates

`--qcovar` supplies continuous covariates, and `--covar` supplies categorical covariates. Covariate files must use FID/IID consistent with the genotypes and phenotypes. Age and principal components are commonly continuous covariates; batch and group may be categorical according to the design.

```text
FID IID AGE PC1
S001 S001 120 0.013
S002 S002 125 -0.008
```

Avoid duplicate columns, linear combinations of other columns and constant columns duplicating the intercept. Categorical coding and reference levels affect coefficient interpretation. Record models differ in their support for individual-level and record-level covariates. Do not treat a one-row-per-individual file as a one-row-per-record file.

### 3.6 Relationship matrices and pedigrees

A binary GRM must be accompanied by its matching ID file. Do not reorder IDs independently after generating the matrix. Sample order in FAM, matrix IDs and the samples actually matched in the model must agree. GRM construction methods use different scaling and variance interpretations.

Use `--pedigree FILE` for a pedigree prepared as `ID SIRE DAM`. Follow the selected module's conventions for missing parents and non-genotyped individuals. Check duplicate individuals, self-parenting, pedigree cycles and unmatched IDs. The meaning of legacy `--ped` depends on the task; distinguish `--pedigree` from `--ped-genotypes` explicitly.

### 3.7 Set and summary statistic inputs

Set tests require a variant-grouping file. COJO, meta-analysis and fine-mapping require summary-statistic fields meeting each module's format requirements and the necessary LD information. Variant ID, chromosome, position, effect allele and other allele must be verifiable.

Renaming external columns does not complete harmonization. Address allele swaps, strand orientation, palindromic variants, genome-build differences, duplicate variants and frequency discrepancies. Check current help and module requirements for exact fields. Before a full analysis, inspect reading logs on a small set of variants.

### 3.8 Adapting raw tables

Raw CSV, TSV and TXT files can be organized into model inputs through `--input-file` and the configuration workflow. `--preview-schema` and `--infer-columns` preview fields and assist identification. `--prepare-input-only` prepares files without running analysis. This workflow is not an Excel workbook reader; export `.xlsx` to plain text first.

`--id-col`, `--iid-col`, `--time-col` and `--trait-cols` select fields in a raw table. `--tpheno` selects a time column in an already prepared phenotype file; the options are not aliases. Save the adaptation configuration, duplicate-record policy and output directory. Inspect date parsing, wide-to-long conversion and individual deduplication before fitting a model.

## Chapter 4 Quick start

### 4.1 Example data and working directory

Run this chapter's commands from the complete package root. `examples/dataset` contains synthetic data for 60 individuals and 180 SNPs across two chromosomes, with single-trait, multi-trait and repeated-measurement phenotypes. These examples check the program and workflow; they do not demonstrate real biological effects or performance.

Linux/macOS and Windows use different shell syntax. Do not copy Bash backslash continuations or variable assignments directly into PowerShell. Set BIN and DATA again whenever you open a new terminal.

### 4.2 Complete Linux and macOS workflow

`$BIN` and `$DATA` below are Bash variables. Initialize them in every new terminal.

```bash
BIN=./bin/cosmicscanner
DATA=examples/dataset
mkdir -p demo_results

# Single-trait LMM with default LOCO
"$BIN" gwas --model lmm --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" \
  --mpheno 1 --threads 1 --out demo_results/lmm

# Three-trait native fixed-null analysis
"$BIN" gwas --model multivariate --bfile "$DATA/synth" \
  --pheno "$DATA/multivariate.pheno" --mpheno-list 1,2,3 \
  --scan-tier approx --reml-maxit 1200 --reml-tol 1e-10 --no-loco \
  --threads 1 --out demo_results/mv

# Record-level repeatability; fixed variances only for this synthetic example
"$BIN" gwas --model repeatability --bfile "$DATA/synth" \
  --pheno "$DATA/long.pheno" --tpheno 1 --mpheno 2 \
  --repeatability-mode exact --scan-tier approx --no-loco \
  --reml-vg 0.5 --reml-vpe 0.2 --reml-ve 0.3 --threads 1 --out demo_results/rep

# Relationship matrix
"$BIN" grm --bfile "$DATA/synth" --threads 1 --out demo_results/grm
```

The fixed variances 0.5/0.2/0.3 are demonstration settings, not general estimates for real studies. For real data, omit the three options to estimate the components and inspect fitting diagnostics.

| Output | Checks |
|---|---|
| `lmm.assoc.txt` | 180 SNPs, chr1/chr2; current unified fields include BETA_A/SE_A/P_A |
| `lmm.cosmic.sumstats` | The same 180 SNPs, with the header and both chromosomes |
| `lmm.scan.meta.txt`, `lmm.run_manifest.json` | Resolved settings, samples and original command |
| `mv.assoc.mv.txt` | 180 rows, joint test and trait-specific effects/standard errors |
| `rep.assoc.txt`, `rep.repeatability.null.txt` | 180 rows and record-model diagnostics; schemas need not match LMM |
| `grm.grm.bin`, `grm.grm.id` | Matrix and matched sample IDs |

```bash
head -2 demo_results/lmm.assoc.txt
head -7 demo_results/lmm.cosmic.sumstats
cat demo_results/lmm.scan.meta.txt
```

File existence alone does not establish success. Check exit status, record and SNP counts, convergence diagnostics and non-finite values.

### 4.3 Complete Windows workflow

```powershell
$BIN = '.\bin\cosmicscanner.exe'
$DATA = 'examples/dataset'
New-Item -ItemType Directory -Force demo_results | Out-Null
& $BIN --version
& $BIN gwas --model lmm --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" --mpheno 1 --threads 1 --out demo_results/lmm
& $BIN gwas --model multivariate --bfile "$DATA/synth" --pheno "$DATA/multivariate.pheno" --mpheno-list 1,2,3 --scan-tier approx --reml-maxit 1200 --reml-tol 1e-10 --no-loco --threads 1 --out demo_results/mv
& $BIN gwas --model repeatability --bfile "$DATA/synth" --pheno "$DATA/long.pheno" --tpheno 1 --mpheno 2 --repeatability-mode exact --scan-tier approx --no-loco --reml-vg 0.5 --reml-vpe 0.2 --reml-ve 0.3 --threads 1 --out demo_results/rep
& $BIN grm --bfile "$DATA/synth" --threads 1 --out demo_results/grm
(Import-Csv -LiteralPath demo_results/lmm.assoc.txt -Delimiter ([char]9)).Count
(Import-Csv -LiteralPath demo_results/mv.assoc.mv.txt -Delimiter ([char]9)).Count
(Import-Csv -LiteralPath demo_results/rep.assoc.txt -Delimiter ([char]9)).Count
```

The final three counts should each be 180, and the GRM ID file should contain 60 samples. The fixed repeatability variances are only synthetic-example settings.

### 4.4 Success criteria and next steps

Check the program's exit status first. Immediately after the analysis, use `echo $?` on Linux/macOS or `$LASTEXITCODE` in PowerShell. Zero indicates normal completion; also check result files, effective sample size, variant count, convergence and finite values. Do not overwrite the exit status by running another command before inspecting it.

Each of the three example association tables should contain 180 data rows, and the GRM ID file should have 60 rows. The unified single-trait table uses P_A; multi-trait and repeatability tables have their own fields. Do not force all results through the same P-column parser. These counts apply only to the bundled synthetic fixture. Missingness and filtering affect counts in your data.

After completing the examples, substitute your own genotype prefix and phenotype path, choose a new output prefix and start with a small analysis. Do not reuse demonstration repeatability variances 0.5, 0.2 and 0.3 for real data without justification. Estimate them or use justified values and inspect diagnostics.

## Chapter 5 Analysis modules

This chapter explains each module's purpose, inputs, main options and results. Chapter 4 runs directly on bundled data. Paths such as `data/cohort` and `data/trait.txt` below are placeholders for user data; replace them with actual paths. Before running the Bash examples, set BIN to the executable path and create the `results` directory.

### 5.1 Genotype quality control

Quality control excludes samples and variants that do not meet study requirements. `--maf` sets the minimum minor allele frequency, and `--hwe` controls Hardy-Weinberg filtering. The BED path supports `--geno` and `--mind` for variant and individual missingness. `--keep`/`--remove` select samples; `--extract`/`--exclude` select variants.

```bash
"$BIN" gwas --model lmm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 \
  --threads 4 --out results/qc_lmm
```

These thresholds illustrate the options and are not universal standards for every species, population or design. HWE treatment in breeding populations, inbred populations and case-control studies needs a study-specific rationale. Check retained sample and variant counts at each step and save the filtering criteria.

### 5.2 Relationship matrix construction

Relationship matrices estimate genomic relationships between individuals for mixed models or subsequent matrix analysis. The recommended task is `grm`; `--make-grm` is also available. Construction algorithms use different scaling. Report the method and variant filtering.

```bash
"$BIN" grm --bfile data/cohort --threads 4 --out results/cohort
```

Main outputs are `.grm.bin`, `.grm.id` and `.grm.N.bin`. `--make-grm-sparse` produces supported sparse relationship results. Sparsification changes relationship information and is not automatically a lossless replacement for the full matrix. `--make-grm-alg` selects 1 VanRaden, 2 Zeng, 3 Yang, 4 Vitezica, 5 TrueVar or 6 the EMmvGWAS-compatible algorithm. See Appendix A for option values.

### 5.3 Linear model LM

LM performs fixed-effect association analysis for continuous traits. It does not automatically account for relatedness through a random genetic effect; choose covariates and models according to population structure and design. Select `--model lm`; `--lm-mode` controls supported test modes.

```bash
"$BIN" gwas --model lm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --qcovar data/qcov.txt --threads 4 --out results/lm
```

A common output is `.lm.assoc.txt`; verify columns from the actual header. Check phenotype variance, covariate rank and effect direction. Differences between LM and LMM P values may reflect relatedness. Do not select a result solely by its number of significant variants.

### 5.4 Linear mixed model LMM and LOCO

LMM models sample correlation through a random genetic effect and tests fixed effects at each variant. Select `--model lmm`. The standard multi-chromosome path uses LOCO by default: when scanning one chromosome, its information is excluded from the relationship structure. `--no-loco` explicitly disables LOCO.

```bash
"$BIN" gwas --model lmm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --qcovar data/qcov.txt --threads 4 --out results/lmm
```

The current unified `.assoc.txt` uses additive fields BETA_A, SE_A, P_A and LOG10P_A. LOCO changes the null model and computational cost; record its setting when comparing analyses. For a single chromosome or a restricted path, inspect logs for the actual behavior.

`--reml-maxit` and `--reml-tol` control the relevant fitting iterations. Raising the iteration limit addresses only some nonconvergence problems. First inspect inputs, covariates, matrices and variance boundaries. Wald, LRT and score modes selected by `--lmm-gemma-mode` have different hypotheses. Compatibility naming does not establish identical implementation details.

### 5.5 Additive dominance effects and variant coding

`--gwas-model additive`, `additive_dominance` and `dominance` select the effect structure within a single-trait family. Additive-dominance joint models and dominance models cannot be interpreted through a single phenotype-output field alone. The joint additive-dominance model uses BED/PGEN hard calls; dominance also has the corresponding BGEN probability path.

```bash
"$BIN" gwas --model lmm --gwas-model additive_dominance \
  --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --threads 4 --out results/additive_dominance
```

`--snp-model additive`, `dominant`, `recessive` and `heterozygote` specify the coding of the tested variant, distinct from the random-effect structure above. Inspect additive, dominance and joint fields separately. A joint test is not a single additive P value.

For BGEN dominance, also inspect missingness, INFO_D and probability-entropy diagnostics. Some compressed or additional unified output formats have restrictions on additive-dominance paths. Respect explicit program checks instead of bypassing them and combining incompatible outputs.

### 5.6 Low rank and matrix free LMM

`--fastlmm`, `--lmm-low-rank` and `--fastlmm-matrix-free` provide low-rank or matrix-free computation for analyses needing to manage matrix storage and cost. Approximation error depends on rank; compare with a computationally feasible full model. `--fastlmm-trace-threshold` controls a related approximation target. Minimum and maximum ranks, random seed and error tolerances also affect results.

`--fast-lmm` and `--fastlmm` are different switches. The former is an alias for a related fixed-null scan; the latter concerns low-rank computation. Do not substitute them because their names look similar. Compare effects, standard errors, P values and diagnostics against standard LMM on manageable data, and record any fallback.

### 5.7 FarmCPU and BLINK

FarmCPU and BLINK use iterative or variant-selection procedures for multi-locus association. Entry options include `--farmcpu` and `--blink`. FarmCPU's `--farmcpu-max-iterations` and `--farmcpu-pseudo-qtn` control iterations and the number of candidate pseudo-QTNs. Inspect stopping reasons and selected variants.

Save iteration logs, selected variants and stopping reasons. Changing the iteration budget, candidate count or filtering changes the final model. More significant results do not prove greater accuracy. Retain an LMM analysis on the same data to compare effects, significance and population-structure adjustment.

FarmCPU example; replace the data paths:

```bash
"$BIN" --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --farmcpu --farmcpu-max-iterations 5 --out results/farmcpu
```

### 5.8 Binary logistic analysis and Firth correction

A 0/1 logistic GLM is suitable for binary outcomes in independent samples. `--logistic score` or `--logistic wald` selects supported tests. Confirm that cases are coded 1 and controls 0, and that missingness and covariates are handled correctly.

```bash
"$BIN" --bfile data/cohort --pheno data/binary.txt --mpheno 1 \
  --logistic score --threads 4 --out results/logistic
```

Complete separation, very few cases and low MAC can destabilize ordinary estimates. `--firth-logistic` selects a Firth-corrected path; inspect estimates and convergence. `--logistic-glmm` uses a PQL approximation. `--saige-style` selects the corresponding mixed-model score and SPA calculations and requires `--biobank-mode`. These methods differ in fitting and testing; record the actual method. Inspect estimation and testing stability carefully for very few cases or low-MAC variants.

### 5.9 Multi trait mixed models

Multi-trait models jointly analyze correlated continuous traits. Use `--model multivariate` and `--mpheno-list` to select fields after the two IDs. Chapter 4 provides the checked three-trait BED fixed-null example.

`--mv-method` selects an implementation, `--mv-test` selects supported tests and `--scan-tier` controls per-variant variance treatment. Native `approx` denotes the corresponding fixed-null approximation. Native multi-trait computation supports `approx` and conditionally supported `exact`; `semi` is available only with `--mv-compatible emmvgwas`. Confirm that the compatibility model and its variance treatment match the intended analysis before using it.

The `exact` tier requires `--no-loco` and the path-specific input conditions. It re-estimates variance components for each variant; P values still use the reference distribution of the selected test. In `.assoc.mv.txt`, WALD_STAT and P represent the joint test, while BETA_1/SE_1 and corresponding columns hold trait effects. Inspect residual and genetic covariance fitting, numerical rank and convergence.

### 5.10 Repeatability models

Repeatability analysis uses repeated measurements of the same individual and separates genetic, permanent environmental and residual effects. Use `--model repeatability` with `--tpheno` and `--mpheno` for a long table. Chapter 4 demonstrates record-level computation with a fixed-null scan and deliberately fixed variances.

`--repeatability-mode auto`, `exact`, `fast` and `stcg` select record-model computation paths. `--scan-tier approx`, `semi` and `exact` select per-variant variance updates. `--repeatability-mode exact` is not `--scan-tier exact`. The current per-variant exact tier requires BED, no record weights and LOCO disabled.

`--repeatability-weight` supports specified individual- or record-weight files. Weights need a measurement-precision or design rationale. Cache contents must match genotypes, filters and configuration; do not force reuse across datasets. Save `.repeatability.null.txt` with `.assoc.txt`, and inspect individual count, record count, variances and fitting status.

### 5.11 Random regression model RRM

RRM models genetic effects that vary over time. Select `--model rrm`. `--rrm-order` sets a supported order from 1 to 4; `--rrm-basis` selects `legendre` or `bspline`. Balanced wide tables and unbalanced long tables require their corresponding inputs and settings. Changing only the model name is insufficient.

`--rrm-gwas-method record` or `coefficient` selects scanning at the record or coefficient level. `--rrm-curve-grid` and `--rrm-top` control curve-output settings, and `--rrm-unbalance` selects the corresponding unbalanced workflow. Compatibility options `gma-fixreg` and `gmat-unbalance` have specific input and variance-treatment conventions; match the input to those conventions.

Before analysis, establish time units, coverage, records per individual and basis order. High orders with sparse time points may be unidentifiable or unstable. Interpret curves and variant effects in relation to the basis and time scale. RRM remains within the documented experimental scope of RC1.

### 5.12 Pedigree animal maternal and common environment models

`--model animal` uses pedigree relationships and requires `--pedigree`. Its A-kernel structure differs from genomic LOCO. `--model maternal` adds maternal structure and requires maternal IDs according to the module. Do not apply the relative phenotype-column numbering of `--mpheno` to `--maternal-id`.

`--model common-env` models shared environmental groups. Supply the correct grouping file and distinguish individual-level and record-level organization. Maternal, pen or batch factors that completely coincide with other structures may prevent separate variance estimation.

Check pedigree and group matching first, inspect fitting diagnostics and compare against known-structure data or an independent implementation. Do not interpret unmeasured shared environment as genetic contribution.

Animal-model example; replace the data paths:

```bash
"$BIN" gwas --model animal --pedigree data/pedigree.txt \
  --bfile data/cohort --pheno data/trait.txt --mpheno 1 --out results/animal
```

### 5.13 Single step analysis Method 1

The unified entry is `--model ssgwas`, which defaults to Method 1. Single-step analysis combines pedigree and genomic relationships. `--single-step`, `--make-agh` and `--write-agh` support the corresponding analysis and matrix outputs. Tau, omega and blend change scaling or blending and require explicit justification.

Method 1 fits on the genotyped H subset and scans fixed SNP effects. It does not produce Method 1 EBVs. Use Method 2 for individual EBVs and variant backsolving. When preparing pedigree, genotyped-individual lists and phenotypes, check genotyped and non-genotyped individuals separately.

### 5.14 Single step backsolving Method 2

`--single-step-method backsolve` explicitly selects Method 2. The example below uses BED, a VanRaden relationship structure, LDLT and fixed Va/Ve, with individual EBVs, variant effects and window summaries. Methods 1 and 2 estimate different quantities and their results must not be mixed without explanation.

Method 2 does not output variant P values by default. Supported paths with explicit `--ssgwas-pvalue-method exact` add uncertainty and P-value outputs. This option uses a selected inverse under LDLT to calculate variant uncertainty; P values use the corresponding reference distribution. PCG and APY use different solving or approximation strategies. Check convergence, core-individual selection and uncertainty calculation separately.

Main files are `.ss_backsolve.ebv.txt`, `.ss_backsolve.snp.txt`, `.ss_backsolve.window.txt` and `.ss_backsolve.meta.txt`. Supported P-value workflows additionally produce `.pvalue.txt`. Retain allele definitions when reporting effects, explain window construction and interpret EBVs at the individual level.

Method 2 fixed-variance example. These values illustrate syntax; replace them with justified study values:

```bash
"$BIN" gwas --model ssgwas --single-step-method backsolve \
  --bfile data/cohort --pedigree data/pedigree.txt \
  --pheno data/trait.txt --mpheno 1 --reml-vg 0.5 --reml-ve 0.5 \
  --ssgwas-solver ldlt --out results/backsolve
```

### 5.15 Threshold count and ordinal models

`--model threshold` analyzes binary traits with a `probit` or `logit` link. Dichotomize a continuous trait only when the study design provides a clear rationale. `--threshold-dichotomize` is not generally recommended preprocessing.

`--model count` provides corresponding `poisson`, `negbin` and `binomial` paths. `--offset` expects the specified per-sample numeric log-exposure input, not a FID/IID table. Prepare `--binomial-pheno` as `FID IID EVENTS TRIALS`. Check nonnegative counts, events not exceeding trials and the definition of exposure.

`--model ordinal` analyzes ordered traits. `experimental-cumulative-binary` and `stacked-psr` are different implementations. Record `fixed-gls` versus `dense-zkz` and fixed versus grid-based variance strategies. Inspect logs for fallback permitted by `--ordinal-allow-fallback`; do not report fallback results as the originally requested model.

A threshold-model summary convergence flag alone does not establish successful fitting. Also inspect iterations, variance boundaries and numerical results. Do not assume ordinal categories have equal continuous spacing.

Count-model example; replace the data paths:

```bash
"$BIN" gwas --model count --bfile data/cohort \
  --pheno data/count.txt --mpheno 1 --family poisson --out results/count
```

### 5.16 Gene environment interaction

Select `--model gxe`; `--env-covar` identifies the environmental covariate according to that model's conventions. The relevant path includes a genetic main effect and gene-environment interaction, with a corresponding joint two-degree-of-freedom test.

Main fields include BETA_G, SE_G, BETA_GxE, SE_GxE, JOINT_CHISQ and JOINT_P. A significant joint P value does not establish a significant interaction term on its own. Units, centering and coding of the environmental variable affect main-effect interpretation. Inspect environmental coverage and population confounding.

### 5.17 Variance components and REML

`--reml` or `--vce` enters the corresponding variance-estimation workflow. Available `--reml-method` solvers include `ai`, `em`, `brent`, `idul`, `idul_plus` and `aidul`, with allowed combinations depending on the model. Identical option names do not establish identical behavior between internal scan null fitting and standalone variance analysis.

`--reml-maxit` and `--reml-tol` control relevant iterations. On corresponding paths, `--reml-vg`, `--reml-vpe` and `--reml-ve` specify variances. Fix components only with a study rationale and report that choice. Do not enter arbitrary values simply to remove nonconvergence.

Inspect nonnegative variance boundaries, likelihood or objective, gradient or iteration status, matrix factorization and sample size. Relationship-matrix scaling, fixed effects and sample subsets change variance interpretation. Heritabilities should not be compared across analyses without these conditions.

### 5.18 BSLMM BVSR and empirical Bayes extensions

`--bslmm`, `--bvsr` and `--bslmm-mcmc` select different approximation or sampling paths. Default approximate outputs belong to the corresponding approximate model. Posterior sampling requires explicit `--bslmm-mcmc` or `--bvsr`. The ridge-regression path `--bslmm 2` is also distinct from sampling models.

Record `--mcmc-burnin`, `--mcmc-samples`, `--mcmc-thin` and `--seed` for sampling. Inspect mixing, effective sample size and stability across seeds. Variant rankings or approximate evidence in some outputs are not necessarily calibrated posterior inclusion probabilities.

EB LMM and EB-Wald use empirical-Bayes-related estimation and testing. Retain standard LMM as a comparison when applying EB adjustment, and record original tests and adjusted metrics separately. Report them in separate columns. Assess significance according to the selected EB method rather than automatically applying ordinary Wald thresholds.

### 5.19 Burden SKAT and set tests

Set tests analyze predefined groups of variants. Use `--burden` or `--skat` with `--set-file`. `--skat-o` uses the currently implemented fixed grid, and `--acat-v` combines within-set evidence using a Cauchy procedure. Check weighting, dependence treatment and P-value calibration against the study design. Define sets, variant weights and allele coding before analysis.

The current default `--rare-max-maf` is 0.5, so it does not automatically restrict inputs to conventionally rare variants. Explicitly set the desired threshold if the design requires, for example, MAF at most 0.01. Choose `--rare-min-mac` according to sample size and design.

Burden has a specific interpretation for aggregated effects with consistent directions. SKAT tests different set hypotheses, and the relevant tail probabilities currently use an approximation. Save the variants actually included, filters and statistical method. Do not compare sets of different sizes solely by raw statistics.

Burden example; replace the data paths:

```bash
"$BIN" --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --burden --set-file data/sets.txt --rare-max-maf 0.01 --out results/burden
```

### 5.20 LD pruning clumping and conditional analysis

`--ld-prune` provides supported LD pruning. `--ld-window` is a physical-distance window; `--ld-step` and `--ld-r2` control step size and correlation threshold. `--ld-clump` uses a statistical table containing SNP, P, CHR and BP, with options such as `--clump-p` controlling selection. See Appendix A for units and defaults.

`--condition-snp` specifies conditioning variants, and `--conditional-joint` controls corresponding joint treatment. Missing conditioning variants, ID conflicts and high correlation may cause instability. Inspect which variants actually enter the model.

These extensions are not automatically equivalent to the full behavior of PLINK or other tools. For an external LD reference, check ancestry, variant matching and genome-coordinate build. Unchecked external LD affects subsequent interpretation.

### 5.21 COJO

`--cojo`, `--cojo-slct` and `--cojo-joint` select supported conditional and joint workflows. `--cojo-summary` supplies summary statistics and `--cojo-ld-ref` supplies the required reference LD. `--cojo-p`, `--cojo-r2` and `--cojo-wind` control significance selection, collinearity and range conditions.

Signed-r LD must agree with variant order and allele direction. Do not substitute an r-squared matrix for signed-r. Use `--cojo-cond` for explicit conditioning variants. Interpret outputs as conditional effects under the selected LD reference and selection conditions, not established independent causal variants.

COJO example; replace the data paths:

```bash
"$BIN" --cojo --cojo-slct --cojo-summary data/summary.txt \
  --cojo-ld-ref data/signed_ld.txt --out results/cojo
```

### 5.22 Meta analysis

`--meta-analysis` reads studies from the file list supplied to `--meta-studies`. `--meta-method` supports `fixed`, `random-dl` and `random-reml`, corresponding to the implemented fixed- and random-effect strategies.

```bash
"$BIN" --meta-analysis --meta-studies data/studies.list \
  --meta-method fixed --meta-min-studies 2 --out results/meta
```

Each listed file must satisfy this module's field requirements. Harmonize alleles, units, frequencies and genome builds before interpreting combined effects and heterogeneity. `--meta-gc` changes the results; record whether it was enabled and the lambda setting. Overlapping study samples require additional consideration and are not automatically resolved by software.

### 5.23 Fine mapping

`--fine-map` uses `--finemap-summary` and `--finemap-region` to specify statistics and region. `--credible-level`, `--finemap-prior-var` and `--finemap-max-signals` control supported credible-set and signal assumptions.

This path requires actual genotype-based LD. A summary table alone is insufficient to claim a fully informed multi-signal analysis. Credible sets depend on priors, LD, candidate region and model assumptions. They do not assign equal causal plausibility to every member or replace functional experiments.

Fine-mapping example; replace the data paths:

```bash
"$BIN" --fine-map --bfile data/cohort --finemap-summary data/summary.txt \
  --finemap-region 1:1000000-2000000 --credible-level 0.95 --out results/finemap
```

### 5.24 Two stage ridge prediction and association

`--biobank-mode`, `--regenie-model-snps`, `--write-pred`, `--use-pred` and `--regenie-step` control the corresponding prediction and scan workflows. Step 0 is the combined workflow; steps 1 and 2 are prediction and scanning. Consult actual help for specific inputs.

Stage 1 generates LOCO predictions using blockwise ridge regression; stage 2 reads those predictions for association scanning. Save training settings, selected variants, prediction files and sample correspondence. Regenerate matching predictions after changing samples or genotypes. Before using this approximate workflow, compare against standard LMM on a manageable subset.

### 5.25 PCA and format conversion

`--pca` and `--pca-k` provide BED/GRM principal-component analysis. Check component count, sample order and matrix construction. Principal components can describe population structure, but study design determines how many to include.

`--convert` with `--to` converts supported inputs to `text`, `npy`, `ped` or `vcf`. Supply the relevant genotype-input options. Check sample and variant counts, IDs, coordinates and alleles before and after conversion. Converting between probabilities, dosages and hard calls may lose information. Successful conversion does not establish support for the resulting format in every downstream model.

## Chapter 6 Website documentation and viewing results

### 6.1 Website entry points and navigation

The website introduces the software, lists platform downloads and provides online documentation. First read the installation instructions, then run the Chapter 4 package examples and check the results.

CosmicScanner executes on a local computer or Linux server. Run analysis commands in a terminal and choose the output location using `--out`. Result files are saved in that computing environment.

### 6.2 Reading result tables

Use a plain-text viewer, R, Python or spreadsheet software supporting tab-delimited files to open association tables. Select fields by header rather than fixed column position. Excel and similar tools may automatically transform variant IDs, date-like strings or scientific notation. Preserve original output and inspect a copy.

The Chapter 4 PowerShell `Import-Csv` examples use a tab delimiter. For very large files, prefer streaming or analysis tools rather than forcing the entire table into a spreadsheet and exhausting memory.

### 6.3 Manhattan and QQ plots

A Manhattan plot displays genomic position against negative-log P values. A QQ plot compares observed statistics with a reference distribution. For current single-trait LMM, use P_A; for multi-trait analysis, use joint P; for repeatability, use its P field. Filter or separately record non-finite and missing values. For P values underflowing to zero, consult the LOG10P field instead of applying a negative logarithm blindly.

Choose significance thresholds in advance according to study design and multiple-testing strategy. QQ deviations may arise from genuine polygenic effects, population structure, model misspecification or numerical problems. The plot alone cannot identify the cause. Plots are representations of statistical results; also retain original tables and plotting settings.

### 6.4 Matching downloads and documentation versions

Choose the package matching your operating system and check the version, SHA256 and `PACKAGE_ID.json`. Update manuals and example data when upgrading. On macOS, choose arm64 or x86_64 according to the chip. Section 2.4 explains installation and system security prompts.

## Chapter 7 Output files

### 7.1 Output prefixes and directories

`--out` specifies a prefix, not just a folder. `--out results/study` creates files beginning with `study` inside `results`. Create the parent directory before running. Use distinct prefixes for models, filters and repeated runs to avoid overwriting results.

### 7.2 Common files

| File suffix | Contents | Main checks |
|---|---|---|
| `.assoc.txt` | Unified single-trait or certain record-model association table | Interpret the header according to the model |
| `.lm.assoc.txt` | LM association results | Check actual test fields |
| `.cosmic.sumstats` | Summary statistics on supported paths | Variant completeness and allele definitions |
| `.assoc.mv.txt` | Multi-trait joint results | Joint statistics and trait effects |
| `.repeatability.null.txt` | Repeatability null diagnostics | Individuals, records, variances and status |
| `.grm.bin` | Binary relationship matrix | Matching ID and count files |
| `.grm.id` | Matrix sample IDs | Order and sample matching |
| `.grm.N.bin` | Matrix variant-count information | Match the matrix with the same prefix |
| `.scan.meta.txt` | Metadata on supported scan paths | Actual model, strategy and sample scale |
| `.run_manifest.json` | Run record on supported paths | Command, configuration and run identity |
| `.ss_backsolve.ebv.txt` | Method 2 individual EBVs | IDs and estimand |
| `.ss_backsolve.snp.txt` | Method 2 variant backsolving | Effects and allele definitions |
| `.ss_backsolve.window.txt` | Method 2 window summary | Window boundaries and aggregation |
| `.ss_backsolve.meta.txt` | Method 2 configuration | Variances, matrices and solver path |

Not every path produces every file. A script checking a universal list of suffixes may incorrectly label a module as failed. Check the outputs promised by the selected module.

### 7.3 Single trait fields

| Field | Meaning | Interpretation requirements |
|---|---|---|
| CHR SNP BP | Chromosome, variant ID and position | Confirm genome build and duplicate variants |
| A1 A2 | Output alleles | Confirm effect direction using model conventions |
| BETA_A | Additive effect estimate | Units depend on phenotype and coding |
| SE_A | Additive-effect standard error | Investigate non-finite or unusual values |
| P_A | Additive-test P value | Interpret with the test and multiple comparisons |
| LOG10P_A | Corresponding negative-log probability field | Check numerical records for extreme values |

Dominance and joint models add corresponding fields. Do not hard-code the first BETA or P column position or assume A1 follows one universal effect-allele convention across all outputs and programs.

### 7.4 Multi trait and repeatability fields

Example multi-trait fields are CHR SNP BP A1 A2 FREQ WALD_STAT P and BETA_1 SE_1, BETA_2 SE_2, BETA_3 SE_3. P belongs to the joint test. Trait effects follow the order in `--mpheno-list`.

The reference repeatability table uses CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P. Check path metadata for the counting level of N; do not equate it automatically with independent individual count. Report both individuals and measurement records in record models.

### 7.5 Logs and abnormal results

Check that output P values lie between 0 and 1 and are numerically valid. Investigate NA, infinities, negative standard errors, abnormal frequencies and extensive boundary variances. Identical results at all variants, zero effective samples or nonconvergence cannot be accepted merely because a file exists.

Saving standard output and error output helps diagnose problems, but logs do not replace original results. For chromosome-wise LOCO, confirm that aggregate files contain all expected chromosomes. Do not mistake the last block for the complete result.

## Chapter 8 Troubleshooting and reproducibility

### 8.1 Installation and path problems

| Symptom | Common cause | Action |
|---|---|---|
| Program or DLL not found | Incomplete extraction or only the executable copied | Extract the complete package again and run from its root |
| Linux glibc incompatibility | System older than package requirement | Use a compatible environment or matching source build |
| BED or phenotype not found | Wrong working directory or incorrect suffix | Check paths; omit `.bed` from `--bfile` |
| Path with spaces fails | Unquoted path | Quote input, output and executable paths |
| Results cannot be written | Missing parent directory or insufficient permissions | Create the directory and use a writable location |
| Mac package cannot run | System version, architecture or execution permissions mismatch | Select the matching macOS 15+ package and follow Section 2.4 |

### 8.2 Input and model problems

| Symptom | Check |
|---|---|
| Fewer effective samples than expected | FID/IID matching, missingness, sample filters and covariate completeness |
| Wrong trait or time column | `--mpheno` and `--tpheno` are counted after the two IDs |
| Singular covariate matrix | Constant, duplicate or collinear columns and categorical coding |
| Multi-trait nonconvergence | Trait variances, correlations, missingness and covariates before increasing iterations |
| Native semi rejected | Use native approx/exact; semi requires the corresponding EMmvGWAS-compatible model |
| Exact scan rejected | no-LOCO, BED, weights and module-specific conditions |
| P column not found | Select P_A, P or another explicitly named model field |
| No Method 2 P values | They are off by default; use explicit options on supported paths |

### 8.3 Performance and resource problems

Start at a manageable scale and determine the algorithm path and memory requirements. More threads do not guarantee linear speedup. Running several analyses with every thread occupied may reduce throughput. Relationship matrices and multi-trait models can consume considerable memory; compressed genotype-file size alone is insufficient to estimate resources.

A nonzero exit status requires inspecting error messages. Do not remove input checks, force arbitrary variances or ignore convergence failures to generate formal conclusions. Explain differences from reference computation when using approximate, low-rank or matrix-free strategies.

### 8.4 Reproducibility records

For each analysis, save at least software version and package identity, input hashes, full command, platform, thread count, model and effect structure, LOCO, scan tier, sample and variant filters, variance treatment, output tables, diagnostics and plotting settings. For sampling, iterative or approximate methods, additionally retain seeds, iterations, fallback records and approximation settings.

Data may be restricted by privacy rules or project agreements. When reporting problems, provide a shareable small reproducer, complete command, error log and dataset dimensions. Do not publish private raw data, passwords or internal server records on the website.

## Appendix A Complete command line help

This appendix reproduces the concise and full command-line help supplied in the Chinese v2.0.0 manual. Search by option name. Appendix B provides a categorized index, and Chapter 5 explains model inputs and combination requirements. Help visibility does not promote an experimental or internal option to supported use.

### A.1 Concise help

```text
CosmicScanner: genetic association analysis
Usage: ./cosmicscanner [gwas|grm] [options]

Common options (PLINK/GCTA-style long flags):
  --bfile PREFIX         PLINK BED/BIM/FAM input
  --pfile PREFIX         PLINK PGEN/PVAR/PSAM input
  --bgen FILE --sample FILE  BGEN input (model-specific support)
  --pheno FILE           FID IID followed by phenotype columns
  --mpheno N             1-based column AFTER FID/IID; default 1
  --mpheno-list 1,2      Phenotype columns for multivariate analysis
  --covar FILE           Categorical covariates (FID IID ...)
  --qcovar FILE          Quantitative covariates (FID IID ...)
  --gwas                 Association task; optional after 'gwas' command
  --model NAME           lmm, lm, multivariate, repeatability, ...
  --scan-tier TIER       approx|semi|exact; see help models for limits
  --no-loco              Global scan; otherwise LOCO where supported
  --maf F --hwe P        Variant filters
  --extract FILE --exclude FILE  Variant-ID inclusion/exclusion
  --threads N            Number of threads
  --out PREFIX           Output file prefix

Example:
  cosmicscanner --bfile data --pheno trait.txt --gwas --model lmm --out run
  cosmicscanner gwas --model multivariate --bfile data --pheno traits.txt --mpheno-list 1,2 --scan-tier approx --reml-maxit 1000 --out mv


Help: --help [topic], help [topic], --help-all (complete legacy/advanced reference)
Manual: docs/CLI_GUIDE_EN.md / docs/CLI_GUIDE_CN.md
Existing flat commands, aliases and statistical defaults are retained.
```

### A.2 Complete help

```text
CosmicScanner v2.0.0
  --version             Print software version and exit
Public scope: A=scoped official contract, B=experimental, C=internal/legacy.
See COSMICSCANNER_V2_SCOPE_FREEZE_PROPOSAL.md; final acceptance is pending.
Formats, engines and output schemas are model-specific; no universal equivalence.
Usage: ./cosmicscanner --bfile <prefix> [options]

=== INPUT/OUTPUT ===
  --pca --pca-k <n>     Experimental BED GRM PCA; default k=10
  --geno <f> --mind <f> Maximum variant/sample missing fraction [0,1]; BED QC
  --keep <f> --remove <f> BED sample lists: FID IID; remove takes precedence
  --ld-prune            Experimental BED window LD pruning
  --ld-clump <f>        Experimental clumping; summary header SNP P CHR BP
  --ld-window <bp> --ld-step <n> --ld-r2 <f> Defaults: 250000, 1, 0.2
  --clump-p <f>         Lead P threshold (default 5e-8)
  --skat-o              Fixed-grid dependence-adjusted SKAT-O
  --acat-v              Experimental set-wise Cauchy combination (--set-file)
  --rare-max-maf <f>    Maximum set-test MAF in (0,0.5] (default 0.5)
  --rare-min-mac <n>    Minimum set-test minor allele count (default 1)
  --snp-model <m>       additive (default), dominant, recessive, heterozygote
                        Fixed-effect marker coding for standard LM/LMM;
                        non-additive --snp-model requires hard-call genotypes.
                        --gwas-model dominance/additive_dominance requires
                        additive --snp-model.
  --gwas-model <m>      additive (default) | additive_dominance | dominance
                        additive:          a + A
                        additive_dominance: a + d + A + D LMM scan (PLINK BED or PGEN hardcall)
                        dominance:         d + D LMM scan (BED/PGEN hardcall or BGEN probability)
  --bgen-dom-max-missing <f>  BGEN dominance variant missingness filter (default 0.05)
  --bgen-dom-min-info <f>     Minimum dominance INFO_D (default 0.01)
  --bgen-dom-entropy-warning <f>  Normalized entropy warning threshold (default 0.50)
  --farmcpu             Experimental native FarmCPU multi-locus GWAS
  --farmcpu-max-iterations <n>  Iteration cap (default 5)
  --farmcpu-pseudo-qtn <n>      Pseudo-QTN cap (default 10)
  --blink               Experimental native BLINK multi-locus GWAS
  --condition-snp <f>   Experimental conditional GWAS; one lead SNP ID per line
  --conditional-joint  Experimental LD/BIC stepwise joint selection for --condition-snp
  --firth-logistic     Experimental Firth penalized binary GWAS
  --logistic <test>    Ordinary logistic GWAS: score or wald
  --logistic-glmm      B: PQL mixed null with full projected working-score variance
                        Not an exact mixed-likelihood fit.
  --saige-style        Experimental binary mixed score/SPA scan; requires --biobank-mode
  --cojo               COJO-style conditional analysis from summary stats + LD reference
  --cojo-slct          COJO stepwise conditional analysis (forward selection)
  --cojo-joint         COJO joint model fit of selected SNPs
  --cojo-summary <f>   GWAS summary statistics file
  --cojo-ld-ref <f>    Dense signed-r LD: SNP CHR BP A1 A2 <ordered SNP IDs>
  --cojo-p <f>         P-value threshold for COJO selection (default: 5e-8)
  --cojo-r2 <f>        LD r^2 threshold for collinearity (default: 0.9)
  --cojo-wind <bp>     Window size in bp for LD (default: 10000000)
  --cojo-cond <f>      Condition on these SNP IDs (one per line)
  --fine-map           Native summary-statistic fine-mapping with real genotype LD
  --finemap-summary <f> Canonical .cosmic.sumstats input
  --finemap-region <r> Locus in chr:start-end form
  --credible-level <f> Credible-set posterior level (default: 0.95)
  --finemap-prior-var <f> Prior effect variance W (default: 0.04)
  --finemap-max-signals <n> Maximum conditional signals (default: 5)
  --finemap-signal-p <f> Conditional signal threshold (default: 5e-8)
  --meta-analysis      Summary-statistics meta-analysis (fixed/random effects)
  --meta-studies <f>   File listing study summary stat files
  --meta-method <m>    Meta-analysis method: fixed, random-dl, random-reml (default: fixed)
  --meta-het-p <f>     Heterogeneity threshold p-value (default: 0.05)
  --meta-min-studies <n> Min studies required (default: 2)
  --meta-gc            Apply genomic control correction
  --meta-gc-lambda <f> Genomic control lambda (default: 1.0)
  --meta-out-format <f> Output format: metal, simple (default: metal)
  --bfile <prefix>       PLINK 1.x binary prefix (.bed/.bim/.fam)
  --bgen <file>          BGEN file (Layout 2, zlib compressed)
  --sample <file>        Sample file for BGEN (Oxford format)
  --pfile <prefix>       PLINK 2.0 prefix (.pgen/.pvar/.psam)
  --pgen <file>          PLINK 2.0 genotype file
  --pvar <file>          PLINK 2.0 variant file
  --psam <file>          PLINK 2.0 sample file
  --vcf <file>           VCF file (bgzip compressed, tabix indexed)
  --vcf-tbi <file>       Tabix index file for VCF (optional if .tbi alongside)
  --vcf-sample <file>    Sample mapping file for VCF (optional)
  --hmp <file>           Standard HapMap file (11 metadata columns + samples)
  --geno-matrix <file>   Text genotype matrix (default: variant rows, sample columns)
  --geno-var <file>      Text variant metadata: ID CHROM POS REF ALT
  --geno-sample <file>   Optional text sample metadata: IID or FID IID
  --geno-orientation <o> variants-by-samples (default) or samples-by-variants
  --npy <file>           NPY genotype matrix (float32, little-endian)
  --npy-var <file>       NPY variant metadata: ID CHROM POS REF ALT
  --npy-sample <file>    NPY sample metadata: FID IID
  --npy-orientation <o>  variants-by-samples (default) or samples-by-variants
  --ped <file>           PLINK PED file
  --map <file>           PLINK MAP file (must be used with --ped)
  --convert <format>     Convert genotype format (input: --bfile/--pfile/--bgen/--vcf/--hmp/--geno-matrix/--npy/--ped+--map)
  --to <format>          Target format for conversion: text, npy, ped, vcf (requires --convert)
  --out <prefix>         Output prefix for conversion
  --batch-models <list>  C: developer batch runner: lm,lmm,fastlmm,saige,farmcpu,blink
  --batch-out <prefix>   Output prefix for batch run (each model writes its own output)
  --out/-o <prefix>      Output prefix (default: cosmic_gwas)
  --out-ssf              EBI GWAS-SSF output on supporting paths only
  --out-gz               Gzip on supporting paths; AD/dominance guards reject it
                          No all-model output schema; sumstats INFO=1 is not measured INFO.
  --extract <file>       File listing SNP IDs to include in analysis
  --exclude <file>       File listing SNP IDs to exclude from analysis

=== PHENOTYPE / COVARIATES ===
  --pheno/-p <file>      Phenotype file: FID IID [TIME] VALUE
  --mpheno/-n <n>        Phenotype column number, 1-indexed (default: 1)
  --tpheno <n>           Time column number for repeatability/RRM
  --qcovar <file>        Quantitative covariates: FID IID VAL1 VAL2 ...
  --covar <file>         Categorical covariates: FID IID CAT1 CAT2 ...

=== INPUT ADAPTATION MVP ===
  --input-file <file>    Raw csv/tsv/txt table to prepare before modeling
  --input-config <json>  Load input adaptation settings from JSON config
  --write-template-config <json> Write an editable input config and exit
  --input-format <fmt>   auto|csv|tsv|txt (default: auto; xlsx not in MVP)
  --input-model <model>  repeatability|rrm|lmm|multivariate (default: repeatability)
  --prepared-input-dir <dir> Directory for prepared manifests and model files
  --prepare-input-only   Write prepared inputs and exit before GWAS
  --id-col <name>        Raw ID column; matched to genotype IID or FID:IID
  --iid-col <name>       Optional IID column when --id-col is FID
  --time-col <name>      Record time/date/batch column
  --trait-cols <list>    Comma-separated trait columns
  --input-trait <name>   Trait to run after preparation (default: first trait)
  --record-qcovars <l>   Record-level quantitative covariates
  --record-ccovars <l>   Record-level categorical covariates
  --subject-qcovars <l>  Subject-level quantitative covariates
  --subject-ccovars <l>  Subject-level categorical covariates
  --time-parser <mode>   auto|date|integer|batch (default: auto)
  --time-unique-policy <p> strict|within-time-index (default)
  --time-index-multiplier <n|auto> Multiplier for duplicate-time keys
  --preview-schema       Write schema_preview.tsv and exit before preparation
  --infer-columns        Suggest ID/time/trait/covariate columns in preview
  --wide-to-long         Treat --trait-cols as wide trait-time columns
  --maf <float>          Minor allele frequency filter (default: 0.01)
  --hwe <float>          HWE exact test p-value filter (default: 0.0, disabled)
  --gwas                 Run standard GWAS; some model flags enable scans implicitly
  --model <name>         Model family: lmm, lm, multivariate, repeatability, rrm,
                          animal, maternal, ssgwas, common-env, threshold, count, ordinal, gxe
                          Legacy --model lmm-eb/eb-lmm remains an EB modifier.
                          Documented analysis commands should include it explicitly.
  --no-loco              Disable leave-one-chromosome-out (run global)
  --scan-tier <tier>     Per-SNP VC handling: approx, semi, exact
                        Multivariate semi requires --mv-compatible emmvgwas; native uses approx
                          approx=fixed null; semi=fixed ratio/profile scale;
                          exact=full per-SNP variance-component refit
                          Not finite-sample exact P. mv/repeatability exact tier needs no-LOCO;
                          repeatability exact tier also requires BED and no record weights.
  --fast-lmm             Legacy alias for fixed-null single-trait scan
  --fastlmm              B: BED Gaussian trace-ranked randomized low-rank LMM
  --fastlmm-trace-threshold <f> Minimum retained GRM trace fraction (default: 0.95)
  --fastlmm-min-rank <n> Minimum automatic rank (default: 1)
  --fastlmm-max-rank <n> Maximum automatic rank (default: 256)
  --fastlmm-matrix-free  B: BED matrix-free low-rank path; rank adequacy is data-dependent
  --fastlmm-no-fallback  Don't fall back to full rank if trace threshold not met
  --lmm-low-rank <k>     Opt-in randomized rank-k LMM (0 keeps exact Eigen default)
  --lmm-low-rank-oversample <n> Randomized eigensolver oversampling (default: 10)
  --lmm-low-rank-seed <n> Deterministic randomized eigensolver seed (default: 42)
  --lmm-low-rank-error-tol <f> Relative Frobenius error limit; exact fallback above it
  --lmm-gemma-mode/-lmm <n> GEMMA LMM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --lmm-gemma-refit      Per-SNP alternative REML refit for GEMMA LMM Wald/LRT
  --lm                   Ordinary linear model flag; use --lm-mode <n>
  -lm <n>                GEMMA-style short LM alias
  --lm-mode <n>          GEMMA LM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --model lmm-eb         Empirical Bayes random-SNP LMM GWAS (experimental; report standard LMM too)
  --gwas-method eb-lmm   Alias for --model lmm-eb
  --test eb-wald         Alias enabling EB-LMM Wald test
  --adjust eb-mbonf      Experimental threshold in .eb.summary.txt; do not report alone
  --spa                  Use Saddle Point Approximation for binary traits
  --threads <n>          Number of OpenMP threads (default: max)
  --memory-gib <value>   SolverPlanner diagnostic memory budget in GiB
  --chr <chr>            Restrict scan to one chromosome
  --from-bp <pos>        Restrict scan start position
  --to-bp <pos>          Restrict scan end position

=== VCE OPTIONS (model-dependent) ===
  --lmm-null-solver <s>  Null model solver: brent, idul, idul_plus, aidul
                          (default: brent)
  --reml-maxit <n>       Max VCE iterations (default: 50)
  --reml-tol <f>         VCE convergence tolerance (default: 1e-6)
  --reml-vg <f>          Fix genetic variance (skip VCE)
  --reml-vpe <f>         Fix permanent environment variance (skip VCE)
  --reml-ve <f>          Fix residual variance (skip VCE)

=== MODEL SELECTION (choose exactly one) ===

--- Standard LMM (default, no flag needed) ---
  Exact single-trait mixed model GWAS with Eigen decomposition.
  Standard for quantitative traits with one record per individual.

--- Repeatability Model ---
  --repeatability-model  Multi-record per individual (long format).
  --pheno format: FID IID TIME VALUE  (--tpheno <col>)
  VCE: additive genetic + permanent environment + residual.
  --repeatability-mode <m> Strategy: auto, exact, fast, stcg
                          A: BED record-level approx/semi; B: fast/STCG/full refit.
                          Auto may switch strategy; exact mode differs from exact tier.
  --repeatability-grm-cache <p> Save/load subject GRM
  --repeatability-weight <file> Record weights: FID IID WEIGHT or FID IID TIME WEIGHT
  --vce-mode <m>         Repeatability VCE refinement: auto, fast,
                          subject-level-exact, full-record-exact
                          auto uses scalable subject-level Woodbury VCE;
                          full-record-exact is validation-only.
  --repeatability-refit-validation Per-SNP VCE refit (small data)
  --repeatability-refit-limit <n> Refit limit (0 = all)

--- Common-Environment Model (Experimental) ---
  --common-env-model     y = Xb + x_j*beta_j + Z_a*a + Z_c*c + e
                         a ~ N(0, K*sigma2_a) (genomic GRM), c ~ N(0, I*sigma2_c)
  --pheno format: FID IID [TIME] VALUE  (long format, --tpheno <col>)
  --common-env-file <f>  FID IID GROUP per subject, or per record with
  --common-env-record-level
  --reml-vg/-vc/-ve <v>  Fix variance components (skip VCE)

--- Random Regression Model (RRM, B experimental) ---
  --rrm-order <n>        Basis order/df-1 (Legendre default, 1-4)
  --rrm-basis <m>        Basis: legendre (default) or bspline
  --rrm-curve-grid <n>   Time grid points for dynamic curves (default: 21)
  --rrm-top <n>          Number of top SNPs for RRM curve outputs (default: 100)
  --rrm-ai-iter <n>      RRM null AI refinement iterations (-1=auto, 0=fast)
  --rrm-sig-threshold <p> P threshold for dynamic time windows (default: 0.05)
  --rrm-unbalance         Unbalanced longitudinal GWAS (long-format --pheno)
  --rrm-unbalance        Unbalanced longitudinal RRM scan (default: coefficient method)
  --rrm-gwas-method <m>  Unbalanced RRM-GWAS: record (record-level x*Phi Method A) | coefficient (two-stage Method B, default)
  --rrm-loco-global-null  Fast approximate RRM LOCO using global null model
  --tpoint <time_col>     Column name for time points (unbalanced mode)
  --tclass <cols>         Time-dependent class covariates (comma-sep column names)
  --tcovar <cols>         Time-dependent q-covariates (comma-sep column names)
  --rrm-compatible <m>   Compatibility mode: gma-fixreg (balanced PLINK) | gmat-unbalance
  --rrm-compatible-varcom <file> GMA .var file for gma-fixreg mode
  --long-pheno <file>    Long-format pheno: FID IID TIME VALUE
  Fits subject-level Legendre coefficients with GRM-based
  random regression. Outputs joint multi-df chi-square test.

--- Multivariate Model ---
  --multivariate         Requires --mpheno-list <c1,c2,...>
                          A: scoped BED reference; PGEN/BGEN extensions B; VCF unsupported.
                          Hybrid requires BED; exact scan tier requires no-LOCO.
  --mpheno-list <list>   Comma-separated phenotype column numbers
  --mv-method <m>        emmv, wald, eb, gemma, gemma-full, gemma-exact, gemma-gxe, gemma-gxe-full, hybrid
  --mv-test <list>       GEMMA-full tests: wald,lrt,score
  --mv-gemma-mode <n>    GEMMA mvLMM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --mv-pnr <p>           GEMMA mvLMM refinement threshold (default 0.001)
  --mv-emi <n>           GEMMA exact EM iterations (default 10000)
  --mv-nri <n>           GEMMA exact refinement iterations (default 100)
  --mv-emp <f>           GEMMA exact EM precision (default 0.0001)
  --mv-nrp <f>           GEMMA exact refinement precision (default 0.0001)
  --mv-crt              GEMMA-style Edgeworth/PCRT p-value correction after refinement
  --mv-refine-discrepancy-top <n> Hybrid: add top Score/Wald discrepancy SNPs
  --mv-refine-discrepancy <x>     Hybrid: min |log10P score-wald diff| (1.0)
  --mv-refine-score-prepass <n>   Hybrid: Score-only audit top n Stage1 SNPs
  --mv-vce-solver <s>    em (default) or aireml
  --mv-em-min-iter <n>   Min EM iterations before AI-REML (10)
  --mv-ai-max-iter <n>   Max AI-REML iterations (50)
  --mv-output-vbeta      Add upper-triangle beta covariance columns
  --mv-emmvgwas-output   Also write EMmvGWAS-style G/R/beta/association CSV files
  --mv-emmvgwas-G <csv>  EMmvGWAS semi_gwas mode: external G.csv from vc_estimation
  --mv-emmvgwas-R <csv>  EMmvGWAS semi_gwas mode: external R.csv from vc_estimation
  --emmvgwas             Alias: --multivariate --mv-method emmv --mv-compatible emmvgwas --no-loco
  --emmvgwas-vc-only     EMmvGWAS vc_estimation route: write G/R/beta and stop before scan
  --emmvgwas-output      Alias for --mv-emmvgwas-output
  --emmvgwas-G/R <csv>   Aliases for EMmvGWAS semi_gwas external G/R inputs
  --env-covar <col>      Environment column for --mv-method gemma-gxe or --gxe-model

--- Maternal Effects Model (B experimental) ---
  --maternal-model       Direct + maternal + PE multi-component model
  --pheno format: FID IID DAM_ID PHENO  (4+ columns)
  --maternal-id <col>    Maternal ID column (default: 3, 1-indexed)
  VCE: direct additive (ad), maternal additive (am),
       direct-maternal covariance (adm), PE (c), residual (e)
  Scan uses effective VCs: Vg=ad+am+(c) with LOCO GRM.

--- Single-Step GWAS (ssGWAS) Model ---
  --single-step          B Method 1 by default: genotyped H-subset fixed-SNP LMM
                          Not full-pedigree MME association; no Method-1 EBV writer.
  --make-agh             Build A-inverse, G, and H-inverse matrices and exit
  --write-agh            Write A-inverse, G, and H-inverse during ssGWAS
  --ped <file>           Pedigree: ID SIRE DAM (no header)
  --genotyped-list <f>   File listing genotyped individual IDs
  --write-ebv            C: reserved Method-1 output switch, writer not implemented
  --h2-tau <f>          Genomic inverse scaling (default: 1.0)
  --h2-omega <f>        A22 inverse correction scaling (default: 1.0)
  --blend-g <f>         G blend weight with I (default: 0.05)
  Method 1 estimates VCs on the genotyped H subset, then scans fixed SNPs.

--- Method 2: full-pedigree ssGBLUP back-solving ---
  --single-step-method lmm|backsolve
                        lmm       = Method 1 fixed-SNP association (default, backward compatible)
                        backsolve = Method 2: ssGBLUP genetic evaluation -> GEBV -> SNP
                                    back-solving -> SNP variance sums/window summaries
                        A reference: BED, fixed positive Va/Ve, VanRaden, LDLT; P default off.
                        PGEN/BGEN backsolve rejected; VC/PCG/operator extensions are B.
  --reml-vg <f> --reml-ve <f>
                        Fixed variance components for Method 2 (canonical/deterministic path)
  --ssgwas-vce <mode>   VCE estimation mode for Method 2 (ai|emai|em|hi|exact|...).
                        If omitted, --reml-vg/--reml-ve are required (fixed-VC path)
  --ssgwas-window-snps <n>
                        SNP-count windows (default: 10); denominator total SNP variance,
                        not model Va or LD-adjusted window variance.
  --ssgwas-tune-g       Align G to A22 (G = b*G + a, preGSf90-style) before inversion
  --ssgwas-blend-a22    Blend G* with A22 instead of I (HIBLUP/preGSf90 reference; backsolve only)
  --ssgwas-solver ldlt|pcg
                        MME solver (default: ldlt reference; pcg = matrix-free PCG)
  --ssgwas-pcg-tol <x>  Matrix-free MME PCG tolerance (default: 1e-8)
  --ssgwas-pcg-maxit <n> Matrix-free MME PCG maximum iterations (default: 2000)
  --ssgwas-pvalue-method off|exact|apy-exact|apy-approx
                        Frequentist P-values for back-solved SNP effects (Phases C1/C2/D).
                        Default off. exact = C1 reference (selected MME inverse via LDLT,
                        O(n_genotyped^2) memory). apy-exact = APY route (Leite et al.
                        2024 GSE 56:59): APY G^-1 in the MME + exact core-slice selected
                        inverse (requires --ssgwas-solver ldlt; O(n_core^2) uncertainty
                        block; choose an explicit estimator below).
  --ssgwas-apy-exact-estimator <core-v2|legacy-hybrid-v1>
                        core-v2: matched projected/core SNP effects and exact selected-core
                        covariance (Bermann 2022 Eq. 21-22), kappa=1/S. Initial scope:
                        fixed VC, tau=omega=1, no tuning, identity blend 0.05, dense Gcc
                        solver and fixed core selection. B: limited validation, no universal calibration.
                        Default legacy-hybrid-v1 preserves full-point/core-variance outputs;
                        C compatibility: unmatched pairing, not supported calibrated inference.
                        apy-approx = Phase D Approx_GinvAPY (Leite Eq. 8/9):
                        APY G^-1 + approximate core PEV C_core = (D_c + lambda*Gcc^-1 +
                        Q*S*Q')^-1 (Schur/Woodbury form of Eq. 8; D = Eq. 9 diagonal
                        effective records); NO selected MME inverse and NO dense n_g^2
                        uncertainty block; works with --ssgwas-solver pcg; status
                        EXPERIMENTAL_APPROXIMATE: raw P-values, no genomic control,
                        strict type-I control NOT guaranteed. Use for research/exploration,
                        not as a calibrated inferential P-value path.
  --ssgwas-apy-core-size <n>
                        Core animals for the APY G^-1 (apy-* modes; >= 1, <= n_genotyped;
                        n = n_genotyped reproduces the C1 exact path)
                        'auto' = C developer-only sizing; excluded from v2 supported inference
                        due to the trace defect; no C98/C99 coverage guarantee.
  --ssgwas-apy-core-file <file>
                        Core-animal ID file (one ID per line, pedigree IDs); overrides
                        the seeded random core selection
  --ssgwas-apy-seed <s>  Seed of the deterministic random core selection (default 12345)
  --ssgwas-apy-core-target-ratio <x>
                        Harness-assigned label core_size / K98 (recorded in meta only; 0 = unset)
  --ssgwas-apy-core-variance <x>
                        Target proportion of genomic variance for automatic core sizing
                        (C98=0.98, C99=0.99; used when --ssgwas-apy-core-size auto)
  --ssgwas-apy-core-eigensolver <method>
                        Eigensolver for automatic core sizing: 'auto' (dense for small,
                        Lanczos for large), 'exact' (dense always), 'lanczos' (matrix-free)
  --ssgwas-grm-mode <mode>
                        VanRaden G construction: 'auto' (dense for small n_g, block-streamed
                        for large), 'dense' (full n_g^2 G), 'streamed' (block-streamed,
                        never materializes the full G; APY methods only)
  --ssgwas-apy-inverse-mode <dense|operator>
                        APY effect inverse application backend (default dense); operator
                        streams Gcn blocks and does not materialize Q in the operator
  --ssgwas-mme-apy-mode <dense|operator>
                        APY precision backend used by the MME (default dense); operator
                        requires --ssgwas-solver pcg and fixed variance components
  --ssgwas-pev-apy-mode <dense|operator>
                        APY Eq. 8 core PEV backend (default dense); operator streams
                        Gcn and avoids Q, explicit Gcc^-1, Kc^-1, and core Vhat
  --ssgwas-hinv-mode <dense|operator>
                        Pedigree A22^-1 backend (default dense); operator uses the
                        exact sparse Schur action and requires PCG + fixed VC
  --ssgwas-pev-scan-mode <scalar|batched>
                        SNP uncertainty scan backend (default batched; scalar is reference)
                        Exact/core-v2 use batched only when explicitly requested; omission
                        preserves historical scalar execution. Legacy hybrid stays scalar.
  --ssgwas-pev-scan-block-size <n>
                        Marker block width for batched uncertainty scan (default 256)
  --allow-large-exact-ssgwas-pvalue
                        Bypass the 8 GiB selected-inverse memory guard of
                        --ssgwas-pvalue-method exact (and apy-exact core slice)
  Method 2 outputs: <out>.ss_backsolve.ebv.txt / .snp.txt / .window.txt / .meta.txt
  With --ssgwas-pvalue-method exact additionally:
    <out>.ss_backsolve.pvalue.txt (EFFECT_ALLELE / SE_ALLELE / Z / P / LOG10P)

--- Pedigree Animal Model (B experimental) ---
  --animal-model        Pedigree-based Simple Animal Model GWAS: y = Xb + x*beta + Za + e,
                        a ~ N(0, A*sigma2_a) with full-pedigree A22 for the analysis sample set.
                        Requires --ped. A is built from the pedigree only (no genomic GRM);
                        LOCO does not apply (pedigree kernel is the same for all chromosomes).
                        Outputs *.assoc.txt (standard format) and *.animal_model.null.txt.
  --write-ped-kernel    Write the analysis A22 kernel (*.animal_model.k.txt + .k.id) for
                        external benchmarking (e.g. GEMMA -k).

--- Threshold/Probit Model (B experimental) ---
  --threshold-model      Binary trait GWAS (PQL + probit/logit)
  --threshold-dichotomize  Allow --threshold-model to dichotomise a non-binary
                         response at 0.5 (off by default; binary input is the contract)
  --link probit|logit    Link function (default: probit)
  PQL with fixed Ve=1 (probit) or Ve=pi^2/3 (logit).
  SPA is not a universal calibration guarantee; null-summary convergence flag is unreliable.

--- Count/Binomial PSR-GWAS Model (B experimental) ---
  --count-model          Count/binomial trait GWAS (PSR/PQL + weighted LMM)
  --family poisson|negbin|binomial Distribution (default: poisson)
  --offset <file>        One numeric log-offset/exposure value per sample for count traits
  --binomial-pheno <file> Binomial FID IID EVENT TRIAL file for event/trial traits

--- Ordinal PSR-GWAS Model (B experimental) ---
  --ordinal-model        Ordinal trait GWAS (cumulative-binary or stacked PSR)
  --ordinal-pheno <file> Ordinal FID IID ORDINAL phenotype file
  --ordinal-mode <name>  experimental-cumulative-binary (default) or stacked-psr
  --ordinal-stacked-engine <name> fixed-gls (default) or dense-zkz
  --ordinal-vc-method <name> fixed (default) or grid-reml
  --ordinal-test <name>  wald (default), score, robust-score, or all
  --ordinal-allow-fallback Fall back from stacked-psr to cumulative-binary
  --psr-null-only        Output null model (Vg, Ve, fixed effects, BLUP, quasi-likelihood)
                        then exit without SNP scan (available for threshold/count/ordinal)
  --psr-maxit <n>       Max PSR/PQL iterations (default: 50)
  --psr-tol <f>         PSR/PQL convergence tolerance (default: 1e-6)
  --psr-mu-eps <f>      Mu clipping epsilon for PSR (default: 1e-6)
  --psr-min-weight <f>  Minimum PSR working weight (default: 1e-8)
  --psr-max-weight <f>  Maximum PSR working weight (default: 1e8)

--- Survival prototype (C developer-only, not supported inference) ---
  --survival-model       Incomplete residual-like LMM prototype; no fitted Cox null
  --pheno format: FID IID TIME STATUS  (2 columns after FID/IID)
  runCoxNull only prints a message; do not use for screening or inference.

--- GxE Reaction Norm Model (B experimental) ---
  --gxe-model            SNP-by-environment interaction GWAS
  --env-covar <col>      Environment variable column in --pheno
  VCE: Vg0 (intercept) + Vg1 (slope on env) + residual.
  Scan: 2-df Wald test for beta_g + beta_gxe jointly.
  Output: BETA_G, SE_G, BETA_GxE, SE_GxE, JOINT_CHISQ, JOINT_P

=== INFRASTRUCTURE ===
  --make-grm             Calculate GRM and exit
  -gk [n]                GEMMA-compatible GRM alias; optional algorithm
  --make-grm-sparse      PLINK-only sparse VanRaden GRM and exit
  --sparse-cutoff <f>    Sparse GRM absolute-value cutoff (default: 0.05)
  --gs/-gs <n>           GEMMA-style S matrix for summary-stat VC (mode 1)
  --gq/-gq <n>           GEMMA-style q vector from --beta summary stats (mode 1)
  --gemma-eigen/-eigen   Calculate GRM eigen decomposition and write eigenD/eigenU
  --calccor/-calccor/-gw GEMMA-style SNP correlation windows; writes <out>.cor.txt
  --windowbp <n>         BP window for --calccor (default: 1000000; 0=no BP limit)
  --windowcm <f>         cM window for --calccor (default: 0; PLINK BIM only)
  --windowns <n>         SNP-count window for --calccor (default: 0=no SNP-count limit)
  --vc/-vc <n>           GEMMA-style VC estimate: 1 HE, 2 REML, 3 AI/IDUL REML
  --reml, --vce          Standalone production REML/VCE estimation
  --reml-method <m>      ai, em, brent, idul, idul_plus, or aidul
  --ci/-ci <n>           GEMMA-style summary-stat PVE confidence interval
  --beta/-beta <file>    Summary statistics for --gq or direct --vc with --cor
  --cor/-cor <file>      GEMMA .cor.txt for direct summary-stat --vc with --beta
  --cat/-cat <file>      SNP category file for summary-stat routes
  --mcat/-mcat <file>    List of SNP category files for summary-stat routes
  --wcat/-wcat <file>    SNP-by-category weights for summary-stat routes
  --wsnp/-wsnp <file>    SNP weights for summary-stat routes
  --study/-study <p>     Study prefix for summary-stat --vc (<p>.Vq/q/size.txt)
  --ref/-ref <p>         Reference prefix for summary-stat --vc (<p>.S/size.txt)
  --mstudy/-mstudy <f>   File listing study prefixes for summary-stat --vc/--ci
  --mref/-mref <f>       File listing reference prefixes for summary-stat --vc/--ci
  --make-grm-alg <n>     Method: 1=VanRaden, 2=Zeng, 3=Yang, 4=Vitezica,
                          5=TrueVar, 6=EMmvGWAS-compatible (default: 1)
                          Note: alg 6 is auto-set by --mv-compatible emmvgwas
  --mass-regress <file>  Multi-phenotype mass regressor (testing)
  --biobank-mode         B: native REGENIE-style two-step research path, not reference equivalence
  --regenie-model-snps <f> Model SNP list for Step 1
  --write-pred <file>    Write step 1 LOCO predictions
  --use-pred <file>      Use existing LOCO predictions
  --regenie-step <n>     0=both, 1=predictions only, 2=association only
  --predict-effect/-epm <f> GEMMA-style prediction from SNP effect table
  --predict-score <f>    Output prediction file (default: <out>.predict.txt)
  --predict-snp-col <c>  SNP ID column in effect table (default: SNP/rs/id)
  --predict-effect-col <c> Effect column (default: BETA/beta/alpha/effect)
  --bslmm/-bslmm <n>     B: native/compatibility modes: sparse ridge, ridge, probit, DAP-style
                          Native sampling is separate; ridge parity does not validate posterior parity.
  --bslmm-top <n>        Sparse top-SNP count for --bslmm 1/3 approximation (default: 100)
  --bslmm-h2 <f>         Ridge/GBLUP shrinkage h2 in (0,1), default 0.5
  --bslmm-mcmc           Native posterior-sampling BSLMM
  --bvsr                 Native posterior-sampling BVSR
  --mcmc-burnin <n>      MCMC burn-in iterations
  --mcmc-samples <n>     MCMC post-burn-in iterations
  --mcmc-thin <n>        MCMC thinning interval
  --seed <n>             Deterministic BSLMM/BVSR seed
  --bf/-bf <file>        DAP BF file for --bslmm 5
  --hyp/-hyp <file>      DAP hyperparameter file for --bslmm 5
  --bslmm-hgrid <n>      DAP h grid size for --bslmm 4 (default: 5)
  --bslmm-rho-grid <n>   DAP rho grid size for --bslmm 4 (default: 5)
  -hmin/-hmax <f>        GEMMA BSLMM/DAP h prior range
  -rmin/-rmax <f>        GEMMA BSLMM/DAP rho prior range
  -pmin/-pmax <f>        GEMMA BSLMM/DAP log10(pi) prior range
  -gmean <f>             GEMMA BSLMM gamma-count prior mean hint
  -hscale/-rscale/-pscale <f> GEMMA BSLMM proposal scale hints
  -mh <n>                GEMMA BSLMM MH step-count hint
  -smin/-smax <n>        GEMMA sparse gamma size bounds (-smax aliases --bslmm-top)
  --block-size <n>       Ridge block size (default: 1000)
  --skat                 Scoped set test with tail approximation, not exact-tail equivalence
  --burden               Run Burden rare-variant test (--set-file)
  --set-file <file>      Variant set definition for SKAT/Burden

Model notes:
  Models are mutually exclusive. Use exactly one --*-model flag.
  Default (no flag): standard univariate LMM (single-trait).
  Report bugs at: https://github.com/anomalyco/PanOmics/issues
```

## Appendix B Common command line options

This appendix groups frequently used options by task. PREFIX means a file prefix, FILE a full file path, N an integer and F a numeric value. Input formats and combinations must still satisfy the selected module. Run `help models`, `help input` or `--help-all` using the executable for further options.

### B.1 Input output and resources

| Option | Meaning and usage |
|---|---|
| `--bfile PREFIX` | BED BIM FAM prefix, without suffix |
| `--pfile PREFIX` | PGEN PVAR PSAM prefix |
| `--bgen FILE --sample FILE` | BGEN and matching sample file |
| `--pheno FILE` | Phenotype file with FID IID first |
| `--mpheno N` | 1-based phenotype index after both IDs |
| `--mpheno-list 1,2,3` | Trait columns; effects follow the supplied order |
| `--tpheno N` | 1-based time index after both IDs |
| `--qcovar FILE` | Continuous covariates |
| `--covar FILE` | Categorical covariates |
| `--pedigree FILE` | Pedigree as ID SIRE DAM |
| `--ped-genotypes FILE --map FILE` | PED MAP genotype input |
| `--threads N` | Threads; defaults to the maximum available when omitted |
| `--memory-gib F` | Memory-budget diagnostic, not an enforced process limit |
| `--out PREFIX` | Output prefix; create the parent directory first |

### B.2 Filtering and relationship matrices

| Option | Meaning and default |
|---|---|
| `--maf F` | Minimum MAF; default 0.01 |
| `--hwe F` | HWE P-value filter; default 0 disables it |
| `--geno F --mind F` | Maximum BED variant/individual missing fractions, from 0 to 1 |
| `--keep FILE --remove FILE` | BED sample-ID lists; remove takes precedence |
| `--extract FILE --exclude FILE` | Include/exclude variant-ID lists |
| `--chr CHR` | Chromosome to scan |
| `--from-bp N --to-bp N` | Physical-position range |
| `--make-grm-alg N` | Default 1; 1 VanRaden, 2 Zeng, 3 Yang, 4 Vitezica, 5 TrueVar, 6 EMmvGWAS-compatible |
| `--make-grm-sparse` | Sparse BED VanRaden GRM |
| `--sparse-cutoff F` | Sparse-matrix cutoff; default 0.05 |

### B.3 Models effects and variance estimation

| Option | Meaning and choices |
|---|---|
| `--model NAME` | lmm, lm, multivariate, repeatability, rrm, animal, maternal, ssgwas, common-env, threshold, count, ordinal, gxe |
| `--gwas-model NAME` | additive, additive_dominance, dominance; default additive |
| `--snp-model NAME` | additive, dominant, recessive, heterozygote; non-additive coding requires hard calls; dominance structure requires additive variant coding |
| `--no-loco` | Disable LOCO and use the genome-wide model |
| `--scan-tier NAME` | approx fixes the null; semi fixes variance ratios while handling scale; exact refits per variant; module conditions apply |
| `--lmm-null-solver NAME` | brent, idul, idul_plus, aidul; default brent |
| `--reml-maxit N` | Relevant VCE iteration limit; default 50 |
| `--reml-tol F` | Relevant VCE convergence tolerance; default 1e-6 |
| `--reml-vg F --reml-vpe F --reml-ve F` | Fixed genetic, permanent environmental and residual variances on corresponding models |
| `--lm-mode N` | 1 Wald, 2 LRT, 3 score, 4 all |
| `--lmm-gemma-mode N` | Corresponding LMM tests: 1 Wald, 2 LRT, 3 score, 4 all |
| `--lmm-gemma-refit` | Per-variant refitting for the corresponding test |

### B.4 Multi trait repeated measurement and random regression

| Option | Meaning and choices |
|---|---|
| `--mv-method NAME` | Multi-trait method; see Section 5.9 and program help |
| `--mv-test LIST` | wald, lrt, score for the corresponding method |
| `--mv-compatible emmvgwas` | EMmvGWAS-compatible model; required for semi |
| `--mv-vce-solver NAME` | em or aireml; default em |
| `--mv-output-vbeta` | Additional effect-covariance output |
| `--repeatability-mode NAME` | auto, exact, fast, stcg; default auto; distinct from scan-tier |
| `--repeatability-grm-cache PREFIX` | Write/read individual GRM cache matching data and settings |
| `--repeatability-weight FILE` | FID IID WEIGHT or FID IID TIME WEIGHT |
| `--rrm-order N` | Basis order 1 to 4, defined as the corresponding basis degrees of freedom minus 1 |
| `--rrm-basis NAME` | legendre or bspline; default legendre |
| `--rrm-unbalance` | Unbalanced long-table workflow |
| `--rrm-gwas-method NAME` | record or coefficient; default coefficient |
| `--rrm-curve-grid N` | Curve time-grid points; default 21 |
| `--rrm-top N` | Variants for curve output; default 100 |

### B.5 Single step and other phenotype models

| Option | Meaning and choices |
|---|---|
| `--single-step-method NAME` | lmm is Method 1; backsolve is Method 2 |
| `--genotyped-list FILE` | Genotyped-individual IDs in the pedigree |
| `--ssgwas-solver NAME` | ldlt or pcg; default ldlt |
| `--ssgwas-pvalue-method off` | Default Method 2 setting: no variant P values |
| `--ssgwas-pvalue-method exact` | Selected-inverse uncertainty under LDLT and corresponding fixed-variance conditions |
| `--ssgwas-window-snps N` | Variants per Method 2 window; default 10 |
| `--logistic NAME` | score or wald; phenotype coded 0/1 |
| `--firth-logistic` | Firth correction; inspect estimates and convergence |
| `--family NAME` | Count model: poisson, negbin or binomial; default poisson |
| `--offset FILE` | Per-sample numeric log-exposure for count models |
| `--binomial-pheno FILE` | FID IID EVENT TRIAL binomial phenotype |
| `--ordinal-pheno FILE` | FID IID ORDINAL ordered phenotype |
| `--link NAME` | Threshold model: probit or logit; default probit |
| `--env-covar COL` | Environmental phenotype column for the corresponding GxE model |

### B.6 Set tests and post association analysis

| Option | Meaning and default |
|---|---|
| `--burden --set-file FILE` | Burden set test |
| `--skat --set-file FILE` | SKAT set test with the corresponding tail approximation |
| `--rare-max-maf F` | Maximum set-variant MAF; default 0.5; explicitly set for rare-variant studies |
| `--rare-min-mac N` | Minimum set-variant MAC; default 1 |
| `--ld-prune --ld-window N --ld-step N --ld-r2 F` | BED LD pruning; defaults 250000 bp, step 1, correlation threshold 0.2 |
| `--ld-clump FILE --clump-p F` | Clump a SNP P CHR BP table; default P threshold 5e-8 |
| `--condition-snp FILE` | One conditioning SNP ID per line |
| `--cojo-summary FILE --cojo-ld-ref FILE` | COJO summary statistics and signed-r LD |
| `--cojo-p F --cojo-r2 F --cojo-wind N` | Defaults 5e-8, 0.9 and 10000000 bp |
| `--meta-studies FILE --meta-method NAME` | Study-file list; fixed, random-dl or random-reml; default fixed |
| `--meta-min-studies N` | Minimum study count; default 2 |
| `--finemap-summary FILE --finemap-region REGION` | Summary statistics and chr:start-end; requires genotype-based LD |
| `--credible-level F` | Credible-set probability level; default 0.95 |
| `--finemap-prior-var F` | Fine-mapping effect prior variance; default 0.04 |
| `--finemap-max-signals N` | Maximum conditional signals; default 5 |

### B.7 Input preparation and help

| Option or command | Meaning |
|---|---|
| `--input-file FILE` | Raw CSV TSV TXT table |
| `--input-config FILE` | Input-preparation JSON configuration |
| `--write-template-config FILE` | Write an editable configuration |
| `--preview-schema --infer-columns` | Preview fields and assist identification |
| `--id-col NAME --iid-col NAME` | Individual fields in the raw table |
| `--time-col NAME --trait-cols LIST` | Time and trait fields in the raw table |
| `--prepare-input-only` | Prepare inputs only; inspect before analysis |
| `--prepared-input-dir DIR` | Prepared-file output directory |
| `--version` | Program version |
| `--help` | Common help |
| `help models` | Model help |
| `help input` | Input help |
| `--help-all` | Complete current program options |

## Appendix C Command migration and analysis records

### C.1 Common entry point migration

| Existing entry | Recommended form | Preserve |
|---|---|---|
| `--lm` | `gwas --model lm` | Input and test mode |
| `--multivariate` | `gwas --model multivariate` | Trait columns, method and scan tier |
| `--repeatability-model` | `gwas --model repeatability` | Time column, record strategy and variances |
| `--single-step` | `gwas --model ssgwas` | Explicitly confirm Method 1 or backsolve |
| `--make-grm` | `grm` | Matrix algorithm and filters |
| `--ped` for pedigrees | `--pedigree FILE` | Do not confuse with genotype PED |
| `--ped` for genotypes | `--ped-genotypes FILE --map FILE` | Do not confuse with pedigrees |

Existing flat commands remain supported. Preserve the model definition, phenotype columns, LOCO and test strategy when migrating entry points. Update result-reading scripts to follow current headers, and save software identities and results separately for runs before and after upgrading.

### C.2 Suggested analysis record template

```text
Study and analysis identifier
Program version and PACKAGE_ID
Downloaded archive SHA256
Operating system and hardware resources
Genotype and phenotype files with SHA256
Sample counts before and after matching
Variant counts before and after filtering
Model and effect structure
LOCO and scan-tier
Covariates and time coding
Estimated or fixed variances and their rationale
Complete execution command
Exit status and convergence diagnostics
Main output files
Multiple-testing and plotting settings
Approximation settings seeds and fallback records
```
