# CosmicScanner GWAS Release Readiness Report

Date: 2026-08-24

Scope: release finalization audit for GWAS/ssGWAS modules. No statistical model
or algorithm was changed in this phase; changes were limited to status wording,
CLI help text, metadata labels, and documentation.

## Final Model Matrix

| Model | Formula / background | CLI | Status | Validation | External benchmark | Recommended usage |
|---|---|---|---|---|---|---|
| General LMM GWAS | `y=Xb+x_j beta_j+g+e`, `g~N(0,G)` | `--gwas` | Stable | core regression + GEMMA accuracy | GEMMA beta/SE/log10P corr >0.999 | primary single-trait GWAS |
| Multivariate GWAS | `Y=XB+x_j B_j+U+E`, `U~G⊗Σg` | `--multivariate` | Stable | multivariate validation matrix + mode checks | GEMMA / EMmvGWAS | primary multi-trait GWAS |
| Repeatability GWAS | record model with additive + permanent environment + residual | `--repeatability-model` | Experimental validated | exact/cache/LOCO/fast + edge cases | RepeatABEL-aligned historical gates | repeated-record GWAS |
| Random Regression GWAS Method A | record-level `x_subject,j × Phi(t)` joint test | `--long-pheno --rrm-order ... --rrm-gwas-method record` | Experimental validated | A-G simulation suite + balanced equivalence | no direct external record-level x⊗Phi reference | time-varying SNP research |
| RRM coefficient Method B | SNP scan on BLUP coefficient phenotypes | `--rrm-gwas-method coefficient` | Stable / compatibility path | smoke + GMAT-compatible checks | GMAT-compatible coefficient route | backward-compatible unbalanced route |
| Common Environment GWAS | record model with additive + group common environment + residual | `--common-env-model` | Experimental validated | simulation, LOCO, record-level, fixed-VC checks | simulation only | shared-environment-adjusted GWAS |
| Maternal GWAS direct | offspring genotype fixed effect, maternal effects in background | `--maternal-model --maternal-snp-test direct` | Experimental validated | maternal suite A-H | sommer VCE background for direct path | maternal-adjusted GWAS |
| Maternal SNP tests | dam genotype or direct+maternal joint fixed test | `--maternal-snp-test maternal|joint` | Experimental validated | 10 scenarios incl. null/collinarity/missing dam | simulation only | model-specific maternal SNP tests (PLINK only) |
| Animal Model GWAS | `a~N(0,A)` pedigree kernel | `--animal-model --ped` | Experimental validated | simulation + order invariance + h2 cases | GEMMA `-k A22` | pedigree-kinship GWAS |
| ssGWAS Method 1 | H-kernel fixed SNP scan | `--single-step --single-step-method lmm` | Experimental validated | smoke + matrix construction checks | HIBLUP not directly equivalent | single-step-informed fixed SNP scan |
| ssGWAS Method 2 backsolve | ssGBLUP MME -> GEBV -> `Z'G^-1 u_g` | `--single-step-method backsolve` | Experimental validated | analytic + causal + shuffle + solver consistency | HIBLUP EBV external validation | classic marker-effect ssGWAS |
| Exact SNP P-value | selected MME inverse, dense exact uncertainty | `--ssgwas-pvalue-method exact` | Stable reference within ssGWAS P-value path | C4/C5, null calibration, C1 frozen gates | Aguilar-style equations | reference / small-to-medium exact inference |
| APY-exact SNP P-value | APY `G^-1` + exact core selected inverse | `--ssgwas-pvalue-method apy-exact` | Experimental validated | C2-G3, D3/D4, D5/D6.5/D7 comparisons | Leite 2024 route semantics | validated scalable APY uncertainty comparison path |
| APY-approx SNP P-value | APY `G^-1` + Eq.8/Eq.9 approximate core uncertainty | `--ssgwas-pvalue-method apy-approx` | Experimental approximate | D6/D6.5 calibration boundary + D7 ranking boundary | Leite 2024 reports inflation / GC discussion | research/exploratory only; not calibrated P-value |

## APY Release Semantics

- `apy-exact`: validated scalable APY uncertainty path; raw P-values; exact core selected inverse; core is required.
- `apy-approx`: `APPROXIMATION_STATUS=EXPERIMENTAL_APPROXIMATE`, `CALIBRATION_STATUS=LIMITED / REGIME_DEPENDENT`, `STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`, `SCALABILITY_SCOPE=UNCERTAINTY_LAYER_ONLY`, `GENOMIC_CONTROL_APPLIED=false`.
- `apy-approx` must not be described as simulation-validated, recommended for ranking/screening, or calibrated P-value inference.
- CosmicScanner does not implement or apply genomic control/rescale for APY-approx; Leite 2024 inflation and GC discussion are documented as literature context only.

## CLI Safety Audit

- `--single-step-method lmm|backsolve`: default remains `lmm`; `backsolve` requires pedigree and fixed/estimated variance components; no output behavior changed.
- `--ssgwas-pvalue-method off|exact|apy-exact|apy-approx`: default remains `off`; help now identifies `apy-exact` as validated scalable APY uncertainty path and `apy-approx` as experimental approximate with no strict type-I guarantee.
- `--ssgwas-apy-core-size`: still required for `apy-*`; default behavior and core selection are unchanged.
- `--maternal-snp-test direct|maternal|joint`: defaults to `direct`; invalid values are rejected; PLINK-only guard remains for maternal/joint.
- `--rrm-gwas-method record|coefficient`: default remains coefficient; record is explicitly record-level x×basis Method A and experimental validated.

## Final Regression Results

| Suite | Result | Notes |
|---|---|---|
| Build | PASS | `make -C build cosmicscanner -j8` |
| ctest | PASS | 14/14 PASS after rerunning serially; earlier parallel run conflicted with standalone multivariate validation output directory |
| GEMMA accuracy | PASS | beta_corr 0.9998708, se_corr 0.9995331, log10p_corr 0.9992931 |
| Multivariate validation | PASS | standalone PASS; ctest rerun PASS |
| Repeatability validation | PASS | exact/cache/LOCO/fast smoke |
| Repeatability edge cases | PASS | empty/all-NA reject, covariates, cache invalidation |
| RRM suite | PASS | base RRM smoke |
| RRM Method A simulation | PASS | scenarios A-G PASS |
| Maternal SNP suite | PASS | scenarios A-H PASS |
| Animal model suite | PASS | all validations PASS |
| Common environment suite | PASS | simulation, LOCO, fixed VC PASS |
| ssGWAS Method 1 | PASS | `run_test_singlestep.sh` smoke |
| ssGWAS Method 2 backsolve | PASS | full backsolve validation + analytic + simulation |
| Exact SNP P-value | PASS | C4/C5 tiny reference |
| APY exact | PASS | C2-G3 apy-exact regression |
| APY approximate | PASS | D3/D4/D5 smoke; D6/D6.5/D7 boundary outputs retained |

## Known Non-Blocking Limitations

- APY-approx is not a calibrated inferential P-value path; D6/D6.5 showed regime-dependent raw-P inflation/deflation and error cancellation, especially outside livestock-like redundancy.
- APY-approx scalability claim is limited to the uncertainty layer; upstream dense G/A22 construction remains a separate scalability backlog.
- RRM Method A, common environment, maternal/joint SNP tests, animal model, and ssGWAS routes remain experimental validated rather than Stable.
- Maternal/joint SNP tests are PLINK-only.
- D7 showed APY-approx@K98 can alter top-k prioritization; APY-approx@full preserved ranking in the tested causal scenarios.

## Release Readiness

Release readiness: **conditional yes** for a software release candidate if experimental labels and APY caveats remain visible in README, tutorials, metadata, and manuscript text.

Modules suitable for manuscript main methods: General LMM, multivariate GWAS,
repeatability GWAS, RRM Method A, animal model GWAS, maternal GWAS, common
environment GWAS, H-matrix single-step GWAS (Method 1), ssGBLUP back-solving
GWAS (Method 2), exact back-solved SNP P-value, and APY-exact as a scalable APY
uncertainty comparison path. Main-manuscript inclusion does not promote an
`Experimental validated` route to `Stable` software status.

Modules requiring `Experimental validated` labeling: repeatability, RRM Method
A, common environment, maternal/joint SNP tests, animal model, H-matrix
single-step GWAS, ssGBLUP back-solving GWAS, and APY-exact. APY-approx remains
`Experimental approximate` and belongs only in experimental/limitations text,
not the recommended main-method set. The canonical paper/software split is in
`PAPER_MODEL_SUPPORT_MATRIX.md`.

Blocking release issues: none after serial regression rerun. The only transient failure was caused by concurrently running ctest and the standalone multivariate validation suite against the same output directory; rerunning serially passed.
