# Synthetic Example Data

All files in this directory are deterministic synthetic fixtures generated
with Python's random seed `20260714` and committed as fixed regression fixtures. They contain
no human, animal, cohort, or real-world identifiers or measurements.

- `test_data.*`: equivalent PLINK 1, PGEN, and BGEN data with 120 synthetic
  samples and 240 SNPs.
- `demo_multivariate.pheno`: two correlated, non-collinear continuous traits.
- `demo_repeatability.*`: three noisy records for each of 60 subjects.
- `demo.ped` and `demo.genotyped.txt`: a synthetic ordered pedigree and subset.

These compact files are intended only for build and smoke validation; they are
not evidence for statistical power or biological conclusions.
