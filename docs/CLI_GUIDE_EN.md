# CosmicScanner Unified CLI User Guide

Updated 2026-10-14. Applies to v2.0 development sources containing this CLI update;
it does not certify a public v2.0 release package. The scope contract remains
`COSMICSCANNER_V2_SCOPE_FREEZE_PROPOSAL.md`. Algorithms, encoding, schemas and
numerical defaults are unchanged.

## Recommended structure

Keep familiar PLINK/GCTA-style input/output flags and select a model explicitly:

```bash
cosmicscanner --bfile data --pheno trait.txt --gwas --model lmm --threads 4 --out results/run
cosmicscanner gwas --model lmm --bfile data --pheno trait.txt --threads 4 --out results/run
cosmicscanner grm --bfile data --threads 4 --out results/grm
```

`gwas` adds only `--gwas`; `grm` adds only `--make-grm`. Flat legacy commands remain
supported. Always include `--gwas` or the `gwas` command in recommended association
examples rather than relying on implicit behavior of some old model flags.
Historical flat `--model lmm/fixed-lmm/lmm-eb/eb-lmm` modifiers retain their legacy
EB-switch behavior, including existing scripts with another model flag. The new
task entry checks conflicting model selections.

## Build and help

Start with an obtained source tree and its matching Core submodule. Do not infer
that a public download exists from these source-build instructions. Requirements:
C++17, CMake >=3.15, Eigen, OpenMP, zlib, zstd and the runtime numerical libraries
selected by the build. Current verification is on the project Linux host.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cosmicscanner -j2
BIN=./build/cosmicscanner
"$BIN" --version
"$BIN" --help
"$BIN" help input
"$BIN" --help models
"$BIN" --help repeatability
"$BIN" --help-all
```

Default help is concise. Topics: gwas, models, input, grm, repeatability, rrm,
ssgwas. `--help-all`, `--help all`, or `--help advanced` retains the complete
legacy/advanced reference. An exposed internal option is not a supported feature.

## Common options and indices

| Option | Meaning |
|---|---|
| `--bfile PREFIX` | BED/BIM/FAM input |
| `--pfile PREFIX` | PGEN/PVAR/PSAM input |
| `--bgen FILE --sample FILE` | BGEN input; model-specific boundaries apply |
| `--pheno FILE` | FID IID followed by phenotype fields |
| `--mpheno N` | 1-based index **after FID/IID** |
| `--mpheno-list 1,2` | Comma-separated trait indices under the same rule |
| `--tpheno N` | Time index after FID/IID |
| `--covar FILE` / `--qcovar FILE` | Categorical / quantitative covariates |
| `--maf F --hwe P` | Variant filters |
| `--extract FILE --exclude FILE` | Included/excluded variant IDs |
| `--threads N --out PREFIX` | Threads/output prefix; create parent directory first |
| `--pedigree FILE` | ID SIRE DAM pedigree for animal/ssGWAS/AGH |
| `--ped-genotypes FILE --map FILE` | PED/MAP genotype input |

`FID IID VALUE` uses `--mpheno 1`. `FID IID TIME VALUE` uses `--tpheno 1 --mpheno 2`.
These are not PLINK2 absolute `--pheno-col-nums`. Legacy `--ped` remains
context-dependent. Raw-table adaptation `--time-col` is not an alias for `--tpheno`.

## Models and computation

| Unified model | Existing route | Boundary |
|---|---|---|
| lmm | Default Gaussian LMM | Scoped reference |
| lm | `--lm` | Linear regression |
| multivariate | `--multivariate` | Requires `--mpheno-list`; BED reference |
| repeatability | `--repeatability-model` | Explicit record strategy recommended |
| rrm | RRM flags | Experimental; retain balanced/unbalanced input settings |
| animal, maternal, common-env | Corresponding model flags | Experimental |
| threshold, count, ordinal, gxe | Corresponding model flags | Experimental |
| ssgwas | `--single-step` | Method 1 by default; Method 2 requires backsolve |

`--gwas-model additive|additive_dominance|dominance` selects effect structure within
a single-trait family; it is distinct from `--model`.
`--scan-tier approx|semi|exact` controls SNP variance treatment, not the model.
Native multivariate semi is not implemented and is rejected. Explicit
`--mv-compatible emmvgwas --scan-tier semi` uses a separate coding/null contract.
Multivariate/repeatability tier exact requires no-LOCO; repeatability additionally
requires BED and no record weights. Record mode exact and scan tier exact are
different dimensions. Omitted tiers preserve historical behavior; native legacy
emmv uses fixed-null calculation and is recorded as approx. No experimental scope
is promoted by this interface update; survival/internal options remain outside
the recommended model selector.

## Complete synthetic examples

Run from the source root. The included R33 fixture contains 60 samples and 180
variants across two chromosomes. Reinitialize these Bash variables in each shell.

```bash
BIN=./build/cosmicscanner
DATA=examples/data/r33_synthetic/dataset
mkdir -p demo_results
"$BIN" --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" \
  --mpheno 1 --gwas --model lmm --threads 1 --out demo_results/lmm
"$BIN" gwas --model multivariate --bfile "$DATA/synth" \
  --pheno "$DATA/multivariate.pheno" --mpheno-list 1,2,3 \
  --scan-tier approx --reml-maxit 1200 --reml-tol 1e-10 --no-loco \
  --threads 1 --out demo_results/mv
"$BIN" gwas --model repeatability --bfile "$DATA/synth" \
  --pheno "$DATA/long.pheno" --tpheno 1 --mpheno 2 \
  --repeatability-mode exact --scan-tier approx --no-loco \
  --reml-vg 0.5 --reml-vpe 0.2 --reml-ve 0.3 --threads 1 --out demo_results/rep
"$BIN" grm --bfile "$DATA/synth" --threads 1 --out demo_results/grm
```

The fixed repeatability variances are synthetic demonstration values, not
estimates to reuse in real studies. Omit those three options to estimate them
for actual data and inspect the fitting diagnostics.

Expected: lmm.assoc.txt has 180 variants and both chromosomes, with current
BETA_A/SE_A/P_A columns; lmm.cosmic.sumstats has the same variant set and a header.
mv.assoc.mv.txt has 180 rows and joint/trait-specific statistics. rep.assoc.txt has
180 rows and rep.repeatability.null.txt records its model diagnostics. GRM writes
grm.grm.bin and grm.grm.id. Output schemas remain model-specific.

```bash
head -2 demo_results/lmm.assoc.txt
head -7 demo_results/lmm.cosmic.sumstats
cat demo_results/lmm.scan.meta.txt
```

Check exit status, marker/sample counts, finite values and diagnostics, not file
existence alone. The run manifest retains the command originally typed by the user.

## Errors and migration

Conflicting tasks (`grm --model lm`, `gwas --make-grm`) fail. Conflicting new model
selectors fail. Multivariate requires explicit columns. Nonconverged multivariate
null fitting fails; inspect data/model settings and increase the iteration budget
when appropriate. Reader support is not an all-model format guarantee.

Migration replaces only model selectors: `--lm` to `--model lm`, `--multivariate`
to `--model multivariate`, `--repeatability-model` to `--model repeatability`.
Retain input, indices, tier, LOCO, variance, thread and output settings. Flat old
scripts still run; historical result/version provenance is not rewritten.

## Design references and validation

Command organization draws on [PLINK2](https://www.cog-genomics.org/plink/2.0/general_usage),
[GCTA](https://yanglab.westlake.edu.cn/software/gcta/),
[GEMMA](https://github.com/genetics-statistics/GEMMA),
[REGENIE](https://rgcgithub.github.io/regenie/options/) and
[HIBLUP](https://www.hiblup.com/tutorials). This is not a claim of universal option
or statistical equivalence.

Run `python3 tests/run_cli_regressions.py` or CTest `cli_regressions` for help,
invalid combinations, new/old numerical parity, GRM parity, pedigree aliases and
typed-command provenance.
