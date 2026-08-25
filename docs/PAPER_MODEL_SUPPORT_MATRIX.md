# CosmicScanner Paper Model Support Matrix

Date: 2026-08-24

This matrix separates **paper inclusion** from **software maturity**. A model
may be recommended for the main manuscript because its formulation and
validation evidence are publishable while still carrying an
`Experimental validated` software label. The only allowed software maturity
labels in this matrix are `Stable`, `Experimental validated`, and
`Experimental approximate`.

## Terminology

- **Stable**: externally benchmarked core route with release regression coverage.
- **Experimental validated**: formulation implemented and supported by
  model-specific analytic, simulation, self-consistency, or external evidence,
  but not claimed as a production-stable interface across every input backend.
- **Experimental approximate**: explicitly approximate route with quantified
  statistical boundaries; strict inferential calibration is not guaranteed.

## Final Support Matrix

| Model | Mathematical formulation | CLI | Validation evidence | External comparison | Paper status | Software status |
|---|---|---|---|---|---|---|
| General LMM GWAS | `y = Xb + x_j beta_j + g + e`, `g ~ N(0, G sigma_g^2)` | `--gwas` | core regression, null/VCE and scan tests | GEMMA beta/SE/log10P correlations >0.999 | Main manuscript recommended | Stable |
| Multivariate GWAS | `Y = XB + x_j beta_j + U + E`, `U ~ MN(0,G,Sigma_g)` | `--multivariate` | multivariate validation matrix, mode/filter/LOCO tests | GEMMA and EMmvGWAS | Main manuscript recommended | Stable |
| Repeatability GWAS | `y = Xb + x_j beta_j + Za + Wc + e`; additive, permanent-environment, residual components | `--repeatability-model` | exact/cache/LOCO/fast validation and edge cases | RepeatABEL-aligned historical comparison | Main manuscript recommended | Experimental validated |
| RRM Method A | record-level `y_r = X_r b + x_subject(r),j Phi(t_r)' gamma_j + u_r + pe_r + e_r`; joint `H0: gamma_j=0` | `--long-pheno --rrm-order ... --rrm-gwas-method record` | scenarios A-G, balanced record-level equivalence, rank diagnostics | no direct external record-level `x×Phi` implementation | Main manuscript recommended | Experimental validated |
| Animal Model GWAS | `y = Xb + x_j beta_j + Za + e`, `a ~ N(0,A sigma_a^2)` | `--animal-model --ped` | simulation, pedigree-order invariance, h2 boundaries, incomplete pedigree | GEMMA `-k A22` | Main manuscript recommended | Experimental validated |
| Maternal GWAS | direct: offspring genotype fixed effect with direct/maternal random background; maternal/joint: dam genotype or 2-df direct+maternal fixed test | `--maternal-model --maternal-snp-test direct|maternal|joint` | scenarios A-H, null/truth recovery, missing dam, collinearity | sommer background comparison for direct path; no external maternal SNP reference | Main manuscript recommended | Experimental validated |
| Common Environment GWAS | `y = Xb + x_j beta_j + Za + Wc + e`, common-environment group effect `c` | `--common-env-model` | simulation, LOCO, record mapping, fixed-VC checks | simulation only | Main manuscript recommended | Experimental validated |
| H-matrix single-step GWAS (Method 1) | fixed-SNP H-kernel LMM: `y = Xb + x_j beta_j + u + e`, `u ~ N(0,H sigma_a^2)` | `--single-step --single-step-method lmm` | H/A/G construction and single-step smoke/self-consistency | not equivalent to HIBLUP/postGSF90 back-solving | Main manuscript recommended, named **H-matrix GWAS** | Experimental validated |
| ssGBLUP back-solving GWAS (Method 2) | ssGBLUP MME on full pedigree, then `ahat_j = M_j' G^-1 uhat_g / S`; SNP/window variance | `--single-step --single-step-method backsolve` | analytic reference, causal recovery, shuffle invariance, LDLT/PCG consistency | HIBLUP EBV and definition-aligned SNP-effect comparison | Main manuscript recommended, named **ssGBLUP back-solving GWAS** | Experimental validated |
| Exact SNP P-value | exact selected MME inverse and exact back-solved SNP-effect uncertainty | `--ssgwas-pvalue-method exact` | C1 C4/C5, null calibration, causal simulation, frozen regression | Aguilar/Leite exact formula lineage | Main manuscript recommended | Stable |
| APY-exact SNP P-value | APY `G_APY^-1` in MME/backsolve plus exact core selected-inverse uncertainty | `--ssgwas-pvalue-method apy-exact --ssgwas-apy-core-size N` | C2-G3, D3/D4, D5, paired D6.5, D7 comparison | Leite Exact_GinvAPY semantics | Main manuscript recommended as scalable APY comparison path | Experimental validated |
| APY-approx SNP P-value | APY/core-only effect plus Eq.8/Eq.9 approximate core uncertainty | `--ssgwas-pvalue-method apy-approx --ssgwas-apy-core-size N` | D3-D7 approximation, calibration, cancellation, and ranking boundary audits | Leite Approx_GinvAPY; literature reports inflation and discusses GC | Experimental / limitations section only | Experimental approximate |

## Required ssGWAS Naming

The generic label `ssGWAS` must not be used without a method qualifier in the
paper when describing a test statistic or output:

1. **H-matrix single-step GWAS (Method 1)** is a fixed-SNP association scan
   with H as the random-effect covariance background. Its output is
   fixed-effect `BETA/SE/P`.
2. **ssGBLUP back-solving GWAS (Method 2)** solves breeding values on the full
   pedigree and transforms genotyped GEBV into marker effects and window
   variance. Optional exact/APY P-values concern the back-solved marker effect,
   not the Method-1 fixed SNP coefficient.

## APY Publication Boundary

- APY-exact may appear in the main Methods and Results as an
  `Experimental validated` scalable APY uncertainty comparison path.
- APY-approx must remain `Experimental approximate`, with
  `STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`, raw uncorrected P-values, and no
  claim of calibrated inference.
- D6/D6.5 calibration limitations and D7 ranking limitations must not be
  hidden by causal-recovery results.
