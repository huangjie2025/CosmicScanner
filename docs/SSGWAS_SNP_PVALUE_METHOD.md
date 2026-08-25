# ssGWAS SNP P-value Method (Phase C1)

Frequentist inference for back-solved ssGBLUP SNP effects (Method 2), following
Aguilar et al. (2019), GSE 51:28, DOI 10.1186/s12711-019-0469-3 - same
PEV/PEC route as BLUPF90/postGSf90 SSGWAS.

Status: **Phase C1 complete. Implementation approved and all gates C1-C7
PASS on commit date (2026-08-17).** See the gate-status table in §14 and the
verification harnesses listed there.

## 1. Why Method-1 P != Method-2 P

- Method 1 (`--single-step-method lmm`): SNP fitted as a **fixed effect** in an
  H-kernel mixed model; output BETA/SE/Wald-P/Score-P. Pests the SNP as a
  fixed-effect term against the H-kernel background (type-C LMM association).
- Method 2 (backsolve): SNP effects are **linear functions of the EBV**
  (random-effect back-solving, Equation 3 below); the inferential object is
  the breeding-value prediction error, not a fixed-effect SE.
- Copying Method-1 SE to Method-2 effects is forbidden (project rules 8/9):
  the covariance structure of `\hat a` differs from the fixed-effect SE
  formula (it depends on `Var(\hat u_g)` and on the back-solving operator B).

## 2. Reference formulas (verified 2026-08, full text: PMC6584984)

Conventions in Aguilar et al. (2019):
- MME written in inverse-residual-scaled form
  `[X'X σe⁻², X'W σe⁻²; W'X σe⁻², W'W σe⁻² + H⁻¹ σu⁻²] [β̂; û] = [X'y σe⁻²; W'y σe⁻²]`,
  so the coefficient matrix C already contains `σe⁻²` and `Var(u − û) = C^{uu}`
  **directly** ("sparse inverse" - selected elements of `C⁻¹` corresponding to
  nonzero entries of the Cholesky factor `LL'`, supernodal factorization,
  e.g. YAMS; factorization cost ~ cubic in number of genotyped animals).
- Back-solve (their Eq. for effect estimates):

  `â = (1−α)b · Z'/(2Σ pᵢqᵢ) · G⁻¹ · û₂`                                    (R1)

  with `û₂` = genotyped EBV, `Z` = gene content matrix, and
  `G = (1−α)(a + b·ZZ'/2Σpq) + α A22` their tuned+blended G.
- SNP-effect prediction error variance (their Eq.):

  `Var(âᵢ) = k² · zᵢ' G⁻¹ (G σu² − C^{u₂u₂}) G⁻¹ zᵢ`,  k = (1−α)b/(2Σpᵢqᵢ)   (R2)

  where `C^{u₂u₂}` is the genotyped block of `C⁻¹`. Mathematically:
  `Cov(u−û, û) = 0 ⇒ Var(û₂) = Var(u₂) − Var(u₂−û₂) = Gσu² − C^{u₂u₂}`,
  and since the marker effect is treated as a fixed parameter,
  `Var(âᵢ − aᵢ) = Var(âᵢ)` (sampling variance of the estimate).
- Test statistic and P (their Eq.):

  `z = âᵢ/sd(âᵢ)`,  `pvalᵢ = 2(1 − Φ(|z|))`                                 (R3)

  with `Φ` the standard normal CDF. **No t-statistic, no empirical df**:
  the reference uses the normal approximation (DF = ∞). Only the fixed-effect
  SE differs between EMMAX-style tests and (R3); the paper shows both are
  z-tests on algebraically equivalent information for the SNP-BLUP model.
- Scaling invariance is used in the paper (standardized effect b̂ᵢ with
  `Var(b̂ᵢ) = (2pᵢqᵢ/σu²)·Var(âᵢ)`; statistic unchanged).

## 3. MME inverse convention in CosmicScanner (mapped, not assumed)

Our MME (single_step_model.cpp `solveGeneticEvaluation`, λ-form):

```
C_λ · θ̂ = rhs,   C_λ = [X'X  X'Z        ],  rhs = [X'y]
                       [Z'X  Z'Z + λ·H⁻¹]        [Z'y]
λ = Ve/Va,  H⁻¹ = A⁻¹ + [τ·Gw⁻¹ − ω·A22⁻¹] on the genotyped block
```

- The coefficient matrix has **no** σe⁻² factor (MME_SCALE = λ-form). Therefore
  `Var(θ − θ̂) = σe² · C_λ⁻¹` (PEV_SCALE = σe² · C_λ⁻¹; this is the standard
  Henderson result: PEV = last-block of scaled inverse in the R⁻¹-form MME).
- `C_λ⁻¹` is the full inverse incl. the intercept column; PEV of animals =
  `σe²·(C_λ⁻¹)_{aa}` (the intercept is coupled through `C_λ⁻¹`, no absorption
  subtlety because we take the full inverse).
- Consistency: `C_λ⁻¹` must be computed from the **same** matrix that was
  factored for the solve. In the LDLT path we retain/sparse-refactorize the
  built LHS. The PCG path (implicit matrix) is **not** supported by the exact
  PEV mode (error if requested together; `--ssgwas-solver ldlt` is the
  default).
- Genotyped block selection: MME animal rows = pedigree order
  (`recs[i].aid-1`), so `(C_λ⁻¹)_{gg}` columns/rows = `genotyped_indices_ + 1`
  (0 = intercept).

## 4. PEV derivation (GATE C2)

BLUP property: `Cov(u − û, û) = 0`. Hence

```
Var(û) = Var(u) − Var(u − û) = H σa² − σe²·C_λ⁻¹
```

- `Var(u) = Hσa²` with H the joint relationship matrix of the model
  (the H whose inverse enters the MME).
- Genotyped block: `H_gg = Gw` (the blended matrix `(1−α)G*+α·A22` in
  `--ssgwas-blend-a22` mode; `(1−α)G_raw+α·I` by default; coherent because
  the MME uses exactly `tau·Gw⁻¹ − omega·A22⁻¹`).
- Therefore the sampling covariance of the genotyped EBV is

```
Var(û_g) = Gw σa² − σe²·(C_λ⁻¹)_{gg}          (D1)
```

The negative term is the PEV of the genotyped animals. (D1) must be PSD;
any negative eigenvalue is clamped only within tolerance (see §8).

## 5. Back-solve operator in CosmicScanner (GATE C3)

Implemented (`backsolveSnpEffects`):

```
q = G*⁻¹·â_g                  (G* = tuned-unblended matrix when
                               --ssgwas-blend-a22, else Gw⁻¹)
zt_j = z_j′ q,  z_j = (x_j − mean)/sd_j       (standardized genotype)
EA_j = (zt_j / sd_j) / S = (M_j′ G*⁻¹ â_g)/S   (M_j = dosage − mean, allele
                                                count deviations; S = Σ_j 2p_j(1−p_j))
```

So the linear map from `â_g` to `EA` is the row vector

```
B_j = M_j′·G*⁻¹ / S,   EA = B · â_g           (D2)
```

which equals Aguilar's (R1) with `k = 1/S`, `α=0, b=1`, Z = M (Centered
gene content; identical estimator class - empirically the HIBLUP benchmark
showed this backsolve matches reference `add_a1` at corr 0.99943 on the
allele scale). IMPORTANT: `B` is OUR implemented operator (our G*),
not Aguilar's blended-G operator; PEV propagation must use (D2) exactly.

## 6. SNP PEC derivation (GATE C1/C3)

Pushing (D1) through (D2) (linear transformation of a BLUP):

```
Var(EA_j) = B_j · Var(û_g) · B_j′
          = [ M_j′ G*⁻¹ (Gw σa² − σe² (C_λ⁻¹)_{gg}) G*⁻¹ M_j ] / S²        (D3)
```

= varA_j − varP_j with

```
varA_j = σa² · w_j′ Gw w_j / S²,   varP_j = σe² · w_j′ (C_λ⁻¹)_{gg} w_j / S²,
w_j = G*⁻¹ M_j column j.
```

- `Var(û_g)` from (D1) is computed once (dense n_g×n_g) and reused for all
  SNPs; per-SNP cost = two quadratic forms.
- Aguilar's (R2) is the special case G* = Gw = G.

## 7. SE / statistic / df (GATE C1)

```
SE_j = sqrt(max(Var(EA_j), 0))     (tolerance-clamped, see §8)
z_j  = EA_j / SE_j
p_j  = 2·(1 − Φ(|z_j|))            Φ via erfc; df = ∞ (normal, per (R3))
```

No t-distribution: the reference formula uses the standard normal (verified
against PMC6584984 full text). ddf is not part of the reference method. The
metadata emitted is `TEST_STATISTIC=Z`, `REFERENCE_DISTRIBUTION=STANDARD_NORMAL`,
`SIDES=2` and **no `DF` field is emitted** (no DF=inf anywhere in the exact
ssGWAS metadata; the only `DF` columns in the program belong to the RRM /
maternal joint F-tests, which are genuinely different statistics).

## 8. Numerical safety (project §21)

- `PEC_jj < 0` with `|PEC_jj| ≤ 1e-10 · max(1, max_diag(Var(û_g))) · (|EA_j|²+1)`
  → clamp to 0, count in metadata, SNP STATUS="NEGATIVE_PEC_CLAMPED".
- `PEC_jj < 0` beyond tolerance → SNP STATUS="NUMERICAL_FAILURE" (no abs(),
  no silent z). Program continues; metadata reports the count.
- `SE_j ≈ 0` (≤ machine-epsilon scale) → STATUS="ZERO_SE" (z/p = NaN)
- `Gw·σa² − σe²Cinv_gg` non-PSD is detected via the clamp counter and
  reported (not silently accepted). PSD tolerance used for the (D1) check is
  `1e-6 · max(1, max|eig|)` (float32 GRM storage noise ~1.2e-7 relative times
  a factor-10 safety margin; the earlier 1e-8 relative bound fired false
  `INVALID_COVARIANCE` on marginally-non-PSD draws at n_g=20, e.g. min eig
  -1.88e-8 vs tol 1.26e-8).

## 9. Scale invariance (GATE C4)

For any constant `c` per SNP: EA′ = c·EA, PEC′ = c²·PEC, SE′ = |c|·SE ⇒
`z′ = z`, `p′ = p` (exact in exact arithmetic). Implemented as a **test**:
- Representation A: allele-count effect EA (canonical, documented).
- Representation B: EFFECT_STD = EA·sd_j (reconstruction-consistent
  standardized scale).
Gate: max |z_A − z_B| / max|z_A| ≤ 1e-9 and identical p (reported in the
tiny-reference test).

## 10. Exact implementation design (Phase C1 only)

- CLI: `--ssgwas-pvalue-method exact` (default: **off**; the only mode in C1,
  "exact/reference"). Not default because the PEV step is O(n_g·nnz(L)) +
  O(n_g²) memory (see §12).
- Flow after `solveGeneticEvaluation` (LDLT required):
  1. Reuse/refactorize the sparse `C_λ` (SimplicialLDLT).
  2. Selected inverse: for each genotyped column `j`: solve `C_λ x = e_j`
     (LDLT.sparseSolve), keep the n_g-length slice → dense `(C_λ⁻¹)_{gg}`
     (n_g × n_g). Symmetrize by averaging.
  3. `V = Var(û_g) = Gw·σa² − σe²·(C_λ⁻¹)_{gg}` (dense n_g×n_g; PSD check).
  4. `W = G*⁻¹·M′` (n_g × m) via existing dense Ginv.
  5. Per SNP: `EA_j` (existing formula), `Var(EA_j)` = quadratic forms vs V,
     `SE, z, p, -log10p`; write `*.ss_backsolve.pvalue.txt` ONLY when the
     flag is on (**no change** to `*.ss_backsolve.snp.txt` schema).
- Metadata (always added fields, value "off" when disabled):
  `P_VALUE_METHOD`, `MME_DIMENSION`, `MME_INVERSE_METHOD` (selected-inverse
  LDLT, n_g RHS), `MME_SCALE` (lambda-form, no sigma_e^-2), `PEV_SCALE`
  (sigma_e2*C_lambda_inv), `EFFECT_SCALE` (allele-dosage, EFFECT_ALLELE),
  `SNP_PEC_SCALE` (Var(estimator), Aguilar), `TEST_STATISTIC` (Z),
  `REFERENCE_DISTRIBUTION` (STANDARD_NORMAL), `SIDES` (2),
  `REFERENCE_FORMULA` (Aguilar2019-GSE-51-28 R1-R3),
  `EXACT_OR_APPROX` (exact), `PEC_NEGATIVE_CLAMPED`, `PEC_NUMERICAL_FAILURE`.

## 11. Independent tiny reference (GATE C5) - DONE

`tests/run_ssgwas_pvalue_tiny_reference.py` - pure Python/numpy, no
CosmicScanner helpers: build A (Meuwissen-Luo), G, H, H⁻¹, full MME C_λ,
full dense C_λ⁻¹, EBV via MME solve, PEV = σe²·(C_λ⁻¹)_{aa}, backsolve EA,
PEC via (D3) with dense matrices, SE/z/p. Compare vs program output on the
same seed-fixed tiny data (n_ped=30, n_geno=20, M=30, 6 tags = a22/identity
x h² 0.25/0.60/0.80). Also gates C4 (scale invariance: EFFECT_STD/allele
representations at 4.4e-16 / 1.1e-16) on the program's own outputs.

Result (2026-08-17): all 6 tags PASS at tolerance 1e-5; measured max diffs
EBV/EA/SE/Z/LOG10P ≈ 1.4e-6..5.0e-6. The floor is the C++ GRM storage
precision: `GrmBuilder` accumulates `G += Z·Z′` from float32 geno blocks, so
the fitted Gw differs from the exact-double reference by ~eps32·‖G‖ per
entry; the inverse amplifies that by κ(Gw) (~30-40 at n_g=20) into H⁻¹/G⁻¹
entry noise up to ~1e-6 and EBV/Z effects ~5e-6. Plain double arithmetic
(with exact inverses, after removing the former unconditional
ridge=1e-8 first-attempt in `invertSpdWithRidge`) reproduces the reference
to ~1e-9. This is a documented storage-noise floor, not a modeling error.

Null draws in this harness use `H = A` with the genotyped block replaced by
`Gw` for the ANALYSIS EQUALITY check; this is fine for C5 (C++ vs numpy
comparison) but NOT for calibration - see §11a for the model-consistent
draw used by the null gates.

## 11a. Model-consistent null draws (GATE C6/C7 basis)

The ssGBLUP precision the program fits is `Hinv = A⁻¹ + [0 0; 0 τGw⁻¹ − ωA22⁻¹]`,
which is the inverse of `H = (A⁻¹ + D)⁻¹`, NOT the naive "A with the
genotyped block replaced by Gw" (the two coincide only when
`A21·A11⁻¹·A12 = 0`, which fails for any non-trivial pedigree; on the n=600
test pedigree the two H differ at O(1) in the precision!). Drawing null
phenotypes from the naive block-replaced H puts the data OUTSIDE the fitted
model and inflates every SNP test (measured: lambda_GC 1.07-1.13, FPR@0.05
up to 6.4% at n=600/m=8000). The null gates therefore draw
`a ~ N(0, (A⁻¹ + D)⁻¹ · Va)`, y = 1 + a + e - data and model coincide and
the z statistic is exactly chi2-1 distributed (verified).

## 11b. Null calibration (GATE C6) - DONE

`tests/run_ssgwas_pvalue_null_calibration.py`: n_ped=600 (120 founders
untyped), n_g=480, n_pheno=300, M=8000 spread over 20 chromosomes (single-
chromosome bim would correlate the tests and void the binomial QQ/FPR
bounds), 10 replicates per tag, 4 tags = {a22, identity} x {h² 0.25, 0.60},
fixed VCs = data-generating. Gates (LD-robust: the per-replicate lambda SD
is used to inflate all bounds by the observed correlation factor):
- lambda_GC within 3 sigma of 1; measured 0.9934-1.0276 (nominal SE 0.005)
- FPR at 0.05/0.01/0.001 within 3 binomial MC-SE; measured
  0.0478-0.0520 / 0.0096-0.0109 / 0.00086-0.00126
- QQ decile max deviation within 3 sigma; measured 0.0018-0.0036
All 4 tags PASS. Summary: `tests/results/ssgwas_pvalue_null/null_calibration_summary.json`.

## 11c. Causal recovery (GATE C7) - DONE

`tests/run_ssgwas_pvalue_causal_sim.py`: same population, M=8000, 10 causal
SNPs with fixed alternating-sign effects (total SNP variance 50% of Va;
random per-SNP draws of N(0,σ²) were rejected after measuring |Z|≈1-3 and
<25% of causals above the top-50 threshold at this n - a design-power issue,
not an inference one), y = 1 + a + Xc·b + e; 3 replicates. Both inference
routes on identical inputs: M1 = `--single-step-method lmm` (fixed-SNP LMM,
BETA/SE/Wald-P) and M2 = backsolve + `--ssgwas-pvalue-method exact`.
Result: M2 recovers 17/30 causal SNPs in top-50 (56.7%; chance 2.5%), >=50%
of causal tests at P<0.01; M1 similar or slightly higher (consistent with
the fixed-effect SE being free of the BLUP shrink; documented, not a bug).
M1-vs-M2 concordance: Spearman rho(-log10P) 0.92-0.93, top-50 overlap
0.74-0.80. Summary:
`tests/results/ssgwas_pvalue_causal/causal_sim_summary.json`.

## 12. Runtime / memory (exact mode)

- Factorization: as the LDLT solve already does (n = 1+n_ped).
- Selected inverse: n_g sparse solves ≈ O(n_g·nnz(L)); memory n×n_g ws +
  n_g² block.
- `V = Var(û_g)`: O(n_g²); W build via dense G*⁻¹ M′ O(n_g²·m) (existing
  factors reused).
- Per-SNP O(n_g). Measured budget at n_g=2000/m=2000: seconds (recorded in
  meta TIME_PEV; no speed claims).
- PEV dominates for large n_g (O(n_g·nnz)) - the exact mode is a reference
  tool, not the production scale path.

## 13. Scalable design (design only, NOT implemented in C1)

- Sparse selected inverse of `LL′` (YAMS-style supernodal; the 2019 paper's
  production approach) - only the genotyped columns needed.
- APY-based scalable route per Leite et al. (2024), GSE 56:59,
  DOI 10.1186/s12711-024-00925-3, PMC11340074 ("Marker effect p-values for
  single-step GWAS with the algorithm for proven and young in large genotyped
  populations"): G_APY⁻¹ + core-animal GEBV PEV approximation (Bermann 2022a/b
  lineage) - Phase C2 in progress (see SSGWAS_APY_PVALUE_METHOD.md).
- PCG diagonal extraction / Hutchinson-probing (no explicit inverse at all).
- All deferred to Phase C2.

## 14. Gates C1-C7 (status)

| Gate | Content | Status |
|---|---|---|
| C1 | Aguilar/reference formulas independently verified (full text PMC6584984, R1-R3 quote) | DONE (this doc) |
| C2 | MME inverse -> PEV scaling verified (λ-form vs σe⁻²-form; PEV_SCALE=σe²·C_λ⁻¹) | DONE (derivation + tiny reference numeric) |
| C3 | Back-solve PEC independently derived (D2/D3 from implemented EFFECT_ALLELE) | DONE (derivation + tiny reference numeric) |
| C4 | Effect/SE scale invariance | DONE (tiny reference; 4.4e-16 / 1.1e-16, all 6 tags) |
| C5 | Analytic tiny reference | DONE (6 tags PASS @ 1e-5; float32-GRM storage floor documented, §11) |
| C6 | Null calibration (lambda_GC, QQ, FPR; no threshold tuning) | DONE (4 tags PASS, §11b) |
| C7 | Causal recovery (direction, rank, top-k, power; shrinkage != beta documented) | DONE (M2 top-50 17/30, power@0.01 >= 50%, rho 0.92-0.93, §11c) |

Phase C2 (APY + core-PEV scalable route per Leite et al. 2024) started
2026-08-17 — Phase A audit complete, see docs/SSGWAS_APY_PVALUE_METHOD.md.

## 15. BLUPF90 readiness (C1 scope)

`tests/ssgwas_external/blupf90/` extended (pending implementation): parser
for postGSf90-style SNP effect output + protocol to compare EBV (solutions),
SNP effect, SE and P once reference executables exist. SKIP-safe (exits 0
when binaries absent). HIBLUP: v1.6.0 has no definition-consistent SNP
SE/P output; it remains the genetic-evaluation (EBV) external reference
(Phase 6.6), not a P-value reference.

## 16. External validation status

- EBV: externally benchmarked (HIBLUP, Pearson 1.0, Phase 6.6) - the PEV
  inputs (σe², C_λ⁻¹) therefore inherit the validated evaluation.
- SNP effect scale: definition-aligned (corr 0.99943); the P-value path
  does not depend on HIBLUP's internal scaling (its own PEC is computed
  from (D3), not from a reference-side SE).
- No external P-value reference yet (postGSf90 absent; pending installation).

## 17. Limitations (C1)

- Exact mode: memory O(n·n_g + n_g²), time O(n_g·nnz(L) + n_g²·m); for
  large n_ped/n_g use Phase C2 designs (not implemented here).
- VCs are fixed (or VCE-provided, then treated as known): PEV treats VCs as
  known; standard in the reference method (Aguilar 2019, "variances assumed
  known").
- z/P assume the normal approximation per (R3); no t/df.
- No window P-values (prohibited this round; window variance output
  unchanged).
- No weighted/iterative ssGWAS, no APY (prohibited this round).

## 18. Bibliography (as verified)

- Aguilar, I., Legarra, A., Cardoso, F.F., Masuda, Y., Lourenco, D., Misztal,
  I. (2019). Frequentist p-values for large-scale-single step genome-wide
  association, with an application to birth weight in American Angus cattle.
  Genet. Sel. Evol. 51:28. DOI 10.1186/s12711-019-0469-3. PMC6584984.
- Leite, N.G., Bermann, M., Tsuruta, S., Misztal, I., Lourenco, D. (2024).
  Marker effect p-values for single-step GWAS with the algorithm for proven
  and young in large genotyped populations. Genet. Sel. Evol. 56:59.
  DOI 10.1186/s12711-024-00925-3. PMC11340074 (Phase C2 reference).
- Bermann, M., Lourenco, D., Forneris, N.S., Legarra, A., Misztal, I. (2022a).
  On the equivalence between marker effect models and breeding value models
  and direct genomic values with the Algorithm for Proven and Young.
  Genet. Sel. Evol. 54:52. PMC9288049.
- Bermann, M., Lourenco, D., Misztal, I. (2022b). Efficient approximation of
  reliabilities for single-step genomic best linear unbiased predictor models
  with the Algorithm for Proven and Young. J. Anim. Sci. 100:skab353.
  PMC8827023.
- VanRaden (2008) JDS; Stranden & Garrick (2009) GSE; Wang et al. (2012)
  JDS - back-solving lineage (Phase 6/6.5 verified).

## 19. C0 glossary (C1 closeout)

| Quantity | Mathematical definition | Code variable | Reference terminology |
|---|---|---|---|
| `PEV(a)` | prediction error variance of the genetic value: `Var(u − û) = σe²·(C_λ⁻¹)`; genotyped block `σe²·(C_λ⁻¹)_{gg}` | `mme_lhs_` (retained C_λ); selected inverse `cinv_gg_`; metadata `PEV_SCALE=SIGMA_E2*C_LAMBDA_INV` | Aguilar 2019/Leite 2024: `C^{u₂u₂}` (breeding-value prediction error (co)variance, genotyped animals) |
| `Var(ahat)` = `Var(uhat_g)` | sampling (co)variance of the genotyped GEBV: `Var(û_g) = Gw·σa² − σe²·(C_λ⁻¹)_{gg}` | `vhat_g_` (single_step_model.cpp:1149), computed as `sigma2_a_backsolve_*G_bs_ − sigma2_e_backsolve_*cinv_gg_` | Leite 2024 Eq. 4: `G·σu² − C^{u₂u₂}` |
| `SNP_EFFECT_ALLELE` | back-solved SNP effect on the EFFECT_ALLELE (dosage-deviation) scale: `EA_j = (M_j′·G*⁻¹·û_g)/S`, `S = Σ_j 2p_j(1−p_j)` | per-SNP `snp_effect_allele_(j)`; output column `EFFECT_ALLELE`/`EA` | Leite 2024 Eq. 1: `â = (1−β)b(σ²_u/σ²_a)Z′G⁻¹û` (same estimator class; their (1−β)b/S-scale factor, see §20-Q1) |
| `VAR_EFFECT` | sampling variance of the SNP-effect estimate: `w_j′(Gwσa² − σe²Cinv_gg)w_j/S²`, `w_j = G*⁻¹M_j` | per-SNP `var_effect` (single_step_model.cpp:1221) | Leite 2024 Eq. 4: `k²·z_i′G⁻¹(Gσu² − C^{u₂u₂})G⁻¹z_i`, `k = (1−β)b/(2Σpq)` |
| `SE_EFFECT` | `sqrt(max(VAR_EFFECT, 0))`, tolerance-clamped (see §8) | per-SNP `se`, stored `snp_se_allele_` | Leite 2024 Eq. 3: `sd(â_i)` |
| `Z` | Wald-like standardized effect: `EA/SE` | `snp_z_(j)` | Leite 2024 Eq. 3: `â_i/sd(â_i)` |
| `P` | two-sided standard-normal tail: `2(1 − Φ(|z|))` (erfc) | `snp_p_(j)`, `snp_log10p_(j)` | Leite 2024 Eq. 3: `p-value_i = 2(1−Φ(|\u00e2_i|/sd(â_i)))` |

## 20. C1 closeout: Q1–Q8 final report (2026-08-17)

- **Q1. Exact P 是否与 Aguilar 2019 全文公式一致？** 是。R1–R3 (§2) 逐条对
  PMC6584984 全文核对；CosmicScanner 的 λ-form 取舍已在 §3 映射（PEV =
  σe²·C_λ⁻¹ 与论文 R⁻¹-form 的 C^{-1} 同构）。常量侧记：论文 backsolve 带
  (1−β)b·σu²/σa²（blended-tuned G），CosmicScanner 用 G*-backsolve + 1/S；
  z/P 对全局常数不变（§9, C4），数值上等价。
- **Q2. 与独立参考的数值一致性？** C5 tiny-reference 6 个 tag 全部 PASS
  @tol 1e-5；实测 EBV/EA/SE/Z/LOG10P 最大差 1.4e-6..5.1e-6。纯 double
  运算（移除无条件 ridge 后）与 numpy 参考一致到 ~1e-9；剩余下限是
  float32 GRM 存储噪声（§11）。
- **Q3. Scale invariance？** C4：scale-invariance 测试在全部 6 tag 上
  Z 差 4.4e-16、P 差 ~1e-16。
- **Q4. Null calibration？** C6：λ_GC 0.9934–1.0276（nominal 1）、
  FPR@0.05 = 4.78–5.20%、@0.01 = 0.96–1.09%、@0.001 = 0.86–1.26e-3、
  QQ 分位偏离 ≤0.0036；4/4 tag PASS（LD-robust gate）。关键前置条件：
  null 必须从模型一致精度 (A⁻¹+D)⁻¹ 抽样（§11a），naive block-replaced
  H 会系统性通胀。
- **Q5. Causal recovery？** C7：M2（backsolve+exact）top-50 恢复 17/30
  (56.7%，chance 2.5%)，P<0.01 的 causal 比例 ≥50%；M1 (LMM fixed) 略高
  （固定效应 SE 无 BLUP shrink；预期行为）。M1 vs M2：rho(-log10P)
  0.92–0.93，top-50 overlap 0.74–0.80。
- **Q6. 数值底线？** float32 GRM eps32≈1.2e-7 相对噪声经 κ(Gw)（n_g=20 时
  ~30–40）放大到 EBV/Z ~5e-6 的确定性下限；PSD 容差 1e-6·max(1,max|eig|)
  可安全区分"真实非正定"与"存储噪声"。非建模误差，已文档化（§8, §11）。
- **Q7. Metadata/schema？** exact 模式输出 `P_VALUE_METHOD=AGUILAR2019_EXACT`、
  `TEST_STATISTIC=Z`、`REFERENCE_DISTRIBUTION=STANDARD_NORMAL`、`SIDES=2`，
  **无 DF 字段**（不存在 "DF=inf"，已核实代码与产物文件）；`DF` 仅存在于
  RRM/maternal joint F 检验（独立统计量）。`.ss_backsolve.pvalue.txt` 列
  schema 未变动（CHR SNP BP A1 A2 FREQ EFFECT_ALLELE SE_ALLELE VAR_EFFECT
  Z P LOG10P STATUS）。
- **Q8. 剩余限制与 C2 结论？** exact 模式 memory O(n·n_g + n_g²)、time
  O(n_g·nnz(L) + n_g²·m)；瓶颈 = dense A22⁻¹、dense G⁻¹、dense Cinv_gg。
  作为 reference 工具保留；规模化即 Phase C2（Leite 2024 APY 路线：
  G_APY⁻¹ + core PEV 近似），见 docs/SSGWAS_APY_PVALUE_METHOD.md。