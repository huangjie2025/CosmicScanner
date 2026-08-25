# ssGWAS SNP P-value — Scalable APY Route (Phase C2)

Reference: Leite, N.G., Bermann, M., Tsuruta, S., Misztal, I., Lourenco, D.
(2024). *Marker effect p-values for single-step GWAS with the algorithm for
proven and young in large genotyped populations.* Genet. Sel. Evol. 56:59.
DOI 10.1186/s12711-024-00925-3. PMC11340074.

Status: **Phase A (full-text formula audit + C2 architecture map) COMPLETE
(2026-08-17).** Phases B–E not started (implementation deferred pending user
confirmation). The audit below is against the **published** full text
(PMC11340074, MathML equations extracted directly; not the abstract).

Companion docs: `SSGWAS_SNP_PVALUE_METHOD.md` (C1 exact route, frozen
reference: `--single-step-method backsolve --ssgwas-pvalue-method exact`).

---

## 1. What the paper actually does (per the full text)

Two approximations, **both** required for scalability (the paper is explicit:
replacing only `G⁻¹` by `G_APY⁻¹` while still inverting the LHS gives no
computational win — Exact_GinvAPY measured 110.98 h / 178.30 GB vs
Exact_Ginv 106.76 h / 159.66 GB):

- **A. `G⁻¹ → G_APY⁻¹`** (Eq. 5–6): sparse inverse via the APY recursion.
- **B. core-animal GEBV prediction-error (co)variance without LHS inversion**
  (Eq. 7–10, lineage Bermann 2022a/b): the SNP-effect variance needs only the
  **core block** `C_{u2Cu2C}`; that block is approximated by `C_{u2Cu2C,appx}`
  (Eq. 8), a block-sparse inversion of `D + λG_APY⁻¹`-type LHS contributions
  (D = diagonal effective record contributions, Eq. 9).
- p-values then use `z = â/sd(â)`, `p = 2(1−Φ(|z|))` (Eq. 3) exactly as C1.

Paper pipeline: PREGSF90 (build/save `G` pieces + `G_APY⁻¹`), BLUP90IOD3
(APY + PCG GEBV), ACCF90GS2 (PEVCCappx block-sparse inversion per
Bermann 2022b), POSTGSF90 (backsolve + var(â) + p-values).

Paper numbers: 50K genotyped — Exact_Ginv 106.76 h/159.66 GB,
Exact_GinvAPY 110.98 h/178.30 GB, Approx_GinvAPY 2.83 h/16.62 GB
(≈38× wall-clock, ≈10× memory). 450K genotyped — Approx_GinvAPY
24.47 h/87.64 GB, **linear** in the number of noncore animals.
Correlations vs Exact (50K, 3 replicates): SNP effects 0.87–0.89,
var(â) 0.91–0.92; between the two APY methods ≈ unit.

---

## 2. Reference equations (published text, verbatim reconstruction)

Paper conventions: `β` = blending parameter (5%); `b` = tuning parameter;
`σ²_u` = additive genetic variance; `σ²_a` = SNP variance;
`Z` = gene content centered by 2p; `α = λ = σ²_e/σ²_u`; `S = 2Σp_i(1−p_i)`;
`τ = ω = 1` in their H⁻¹. MME is in inverse-residual form
(`X′X`, `W′W + α·T⁻¹`, no σe⁻² factor) — the same λ-form family as C1.

| Eq | Formula (paper) | Role |
|---|---|---|
| (1) | `â = (1−β)b·(σ²_u/σ²_a)·Z′G⁻¹û` | SNP backsolve (VanRaden 2008; Strandén & Garrick 2009; Wang 2012) with tuned/blended-G correction factor |
| (2) | `G = (1−β)[11′·α + b·ZZ′/S] + β·A22` | G construction: VanRaden type-I `ZZ′/S`, tune `b·(·)+α·11′` (Vitezica 2011), blend 5% with A22 |
| (3) | `p-value_i = 2(1−Φ(|\u00e2_i|/sd(â_i)))` | p-value, standard normal, two-sided |
| (4) | `var(â_i) = k²·z_i′G⁻¹(Gσ²_u − C_{u₂u₂})G⁻¹z_i`, `k = (1−β)b/S` | **Exact** SNP-effect variance; `C_{u₂u₂}` = GEBV PEC from full inverse of LHS |
| (5) | `G_APY⁻¹ = [Gcc⁻¹ 0; 0 0] + [−Gcc⁻¹Gcn; I]·Mnn⁻¹·[−GncGcc⁻¹ I]` | APY inverse (Misztal 2014/2016); identical to Bermann 2022a Eq. (3) |
| (6) | `m_nn,j = g_jj − g_jc′Gcc⁻¹g_cj` | noncore conditional variance (diagonal Mnn) |
| (7) | `var(â_i) = k²·z_cj′Gcc⁻¹(Gccσ²_u − C_{u₂Cu₂C})Gcc⁻¹z_cj` | SNP variance under APY, **exact** given the core GEBV PEC `C_{u₂Cu₂C}`; core-only (Bermann 2022a Eq. 21–22) |
| (8) | `C_{u₂Cu₂C}appx = (D_c + λG_cc − λG_cn(λM_nn⁻¹ + D_n)⁻¹λG_nc)⁻¹` | **Approx** core GEBV PEC; no LHS inversion (Bermann 2022b lineage) |
| (9) | `D ≈ W′(I − X(X′X)⁻¹X′)W` | diagonal effective record contributions after fixed-effect absorption (VanRaden & Freeman 1985; Misztal & Wiggans 1988) |
| (10) | `var(â_i) ≈ k²·z_cj′Gcc⁻¹(Gccσ²_u − C_{u₂Cu₂C}appx)Gcc⁻¹z_cj` | **Approx** SNP variance (combines APY + PEVCCappx) |
| (11) | `y = Xβ + Wu + e` | single-trait animal model (contemporary-group fixed effects in the paper) |
| (12) | `H⁻¹ = A⁻¹ + [0 0; 0 T⁻¹ − A22⁻¹]`, `T⁻¹ = G_APY⁻¹` | ssGBLUP joint inverse (Aguilar 2010) |

Equations (1)–(12) cross-checked against the primary sources:
Bermann 2022a (PMC9288049, its Eq. 1–3, 21–22, 31–32, 34), Bermann 2022b
(PMC8827023, its Eq. 1–4 and block-sparse-inversion description), Misztal
2014/2016 (G_APY⁻¹ recursion), Pocrnic 2016 (core-size rule).

---

## 3. Mapping table (paper → CosmicScanner)

| Ref. eq. | Paper notation | CosmicScanner notation | Implementation object | Exact/Approx |
|---|---|---|---|---|
| (1) | `â_j = (1−β)b(σ²_u/σ²_a)z_j′G⁻¹û` | `EA_j = (M_j′·G*⁻¹·û_g)/S`, `S = 2Σp(1−p)` | `backsolveSnpEffects`, per-SNP `snp_effect_allele_(j)` (C1, frozen); C2: `M_j′G_APY⁻¹û_g/S` (exact-apy) or `M_cj′Gcc⁻¹û_c/S` (approx-apy, Bermann 2022a Eq. 21) | exact (C1) / exact & approx (C2) |
| (2) | `G = (1−β)(α·11′ + bZZ′/S) + βA22` | `G_raw = MM′/S`; tune `a,b`; blend `Gw = (1−α)G*+α·A22` (or `+α·I`) | `GrmBuilder`, `G_bs_tuned_`/`Gw` (`--ssgwas-blend-a22`) | exact |
| (3) | `p_i = 2(1−Φ(|\u00e2_i|/sd(â_i)))` | `P = erfc(\|Z\|/√2) = 2Φ(−\|Z\|)` | `snp_p_`, `snp_log10p_` | exact & approx (same in both) |
| (4) | `var(â_i) = k²z_i′G⁻¹(Gσ²_u−C_{u₂u₂})G⁻¹z_i` | `Var(EA_j) = w_j′·vhat_g_·w_j/S²`, `w_j = G*⁻¹M_j`, `vhat_g_ = Gwσa² − σe²·cinv_gg_` | C1 `cinv_gg_`/`vhat_g_` (frozen reference) | exact |
| (5)+(6) | `G_APY⁻¹` = 2×2 block formula; `m_jj = g_jj − g_jc′Gcc⁻¹g_cj` | `G_APY⁻¹` block-sparse: `Gcc⁻¹` (n_c), `Q = Gcc⁻¹Gcn`, `Mnn⁻¹` (diag) | **C2 new**: `ApyGinv` (no dense n_g×n_g) | approx (exact in the core_size=n_g limit) |
| (7) | `var(â_i) = k²z_cj′Gcc⁻¹(Gccσ²_u − C_{cc})Gcc⁻¹z_cj` | `Var(EA_j)ᵃᵖʸ = w_cj′·vhat_c·w_cj/S²`, `w_cj = Gcc⁻¹M_cj`, `vhat_c = Gccσa² − σe²·(C_λ⁻¹)_cc` | **C2 new**: `cinv_cc_` (exact core slice of selected inverse) | exact (given LHS inverse) |
| (8) | `C_cc^appx = (D_c + λG^{cc} − λG^{cn}(λM_nn⁻¹+D_n)⁻¹λG^{nc})⁻¹` (superscript blocks of Eq. 5) `≡ (D_c + λGcc⁻¹ + Q·S·Q′)⁻¹`, `S_j = λd_j/(λ+d_jm_jj)` | `pec_c_appx_` (Eq. 8 Schur/Woodbury, single production form; gate D1.5: equivalence verified to 8.9e-16) | **Phase D new**: `computeSnpPvaluesApyApprox`; **no n_g² allocation** | approx |
| (9) | `D = diag(W′(I−X(X′X)⁻¹X′)W)` | `D_animals = diag(W′W − W′X(X′X)⁻¹X′W)` per genotyped animal (fixed design known at MME build) | **C2 new**: from design matrices (reuse MME builders) | approx |
| (10) | `var(â_i) ≈ k²z_cj′Gcc⁻¹(Gccσ²_u − C_cc^appx)Gcc⁻¹z_cj` | per-SNP quadratic form on n_c only; batched via `U = Gcc⁻¹M_c` | **C2 new** | approx |
| (12) | `H⁻¹ = A⁻¹ + [0 0; 0 T⁻¹ − A22⁻¹]` | `Hinv = A⁻¹ + [τ·T⁻¹ − ω·A22⁻¹]` on genotyped block (`T⁻¹ = G_APY⁻¹` in C2) | exists (C1); C2 replaces dense `G⁻¹` adapter with APY matvec | exact / approx |

Consistency rule (project §14): in every mode the **effect operator and the
uncertainty operator use the same `G^-1` object** — exact-ginv uses `G*⁻¹`
for both; exact-apy uses `G_APY⁻¹` for the effect and the cc-slice of the
exact LHS inverse for the variance; approx-apy uses `G_APY⁻¹` for the effect
and Eq. (8)/(10) for the variance. Mixing (effect with `G⁻¹`, SE with
`G_APY⁻¹`) is only allowed as an explicitly labeled benchmark option.

### 3.1. Key derivation facts behind the mapping

- **Bermann 2022a Eq. 21–22** (verified in PMC9288049): under the APY model
  the backsolved SNP effect and its variance depend **only on core animals**:
  `â|û = k·Z_c′·Gcc⁻¹·û_c` and `Var(â|û) = k²·Z_c′Gcc⁻¹(Gccσ²_u − PEV_core)Gcc⁻¹Z_c`.
  This is the mathematical basis of Eqs. (7)/(10): only `C_{cc}` (n_c×n_c)
  is needed, never the full `C_{u₂u₂}`.
- **Bermann 2022b Eq. (4)** (verified): the PEV approximation is
  `diag((D + G_APY⁻¹·α)⁻¹)` — block-sparse inversion (Henderson & Searle
  1981 formula 8) of `(D + λG_APY⁻¹)`. Empirical anchor: correlations
  between approximated and exact reliabilities 0.97–0.99 (their Data1).

---

## 4. Audit findings that must be handled in the implementation

### 4.1. Eq. (1) vs (4) scale constants

The paper's backsolve carries `(1−β)b·σ²_u/σ²_a`; their var(â) carries
`((1−β)b/S)²`. In CosmicScanner C1 the effect uses `1/S` with the
unblended tuned `G*`. The two conventions differ by a global constant
(`C4` scale invariance ⇒ z/P identical), but **EFFECT/VAR_EFFECT numeric
values will differ from a BLUPF90 run by that constant**. Cosmetic only for
z/P; recorded in metadata (`EFFECT_SCALE`), not changed in C1.

### 4.2. Eq. (7) exact-vs-APY effect nuance

Bermann 2022a Eq. 21 (`â = kZ_c′Gcc⁻¹û_c`) is exact when the APY recursion
`û_n = P_nc·û_c` holds exactly in the solution. In a fitted ssGBLUP with
noncore phenotypes, the full product `M′G_APY⁻¹û`
(= `M_c′Gcc⁻¹û_c + (M_n′ − M_c′Gcc⁻¹G_cn)·Mnn⁻¹·ξ̂`) differs from the
core-only form by the noncore residual `ξ̂`. Consequences:
- **Exact_GinvAPY effect** = full APY product (like POSTGSF90).
- **Approx_GinvAPY effect** = core-only form (Bermann Eq. 21).
- The difference is part of the "APY error" decomposition measured in
  Phases C/E (`Exact_Ginv vs Exact_GinvAPY` isolates the G-approximation
  error; `Exact_GinvAPY vs Approx_GinvAPY` isolates the PEV-approximation
  error).

### 4.3. Eq. (8): resolved by gate D1.5 (Phase D, 2026-08-17)

1. **Notation resolution (NOT a Leite equation failure).** Eq. (8)'s blocks
   `G^{cc}, G^{cn}, G^{nc}` are the **APY block-inverse pieces** of
   `G_APY⁻¹` (superscript blocks; Eq. 5), **not** the ordinary relationship
   blocks `G_cc, G_cn, G_nc`:
   `G^{cc} = G_cc⁻¹ + Q R Q'`,  `G^{cn} = −Q R`,  `G^{nc} = −R Q'`,
   `R = M_nn⁻¹`, `Q = G_cc⁻¹ G_cn`.
   Substituting into Eq. (8) and simplifying (Woodbury on the diagonal
   `(λR + D_n)`) gives the single Schur/Woodbury form:
   `C_core⁻¹ = D_c + λG_cc⁻¹ + Q·S·Q'`,   `S_j = λ·d_j/(λ + d_j·m_jj)`
   with `m_jj` the noncore conditional variance (`M_nn` diagonal).
   **Gate D1.5**: 100 random SPD partitions (n_c 3–12, n_n 5–18, random
   λ, D with zeros), direct Eq. (8)-superscript vs Schur — worst `max|Δ| =
   8.9e-16` (< 1e-12). The earlier Phase A note of an "O(1) deviation of the
   printed Eq. (8)" was traced to substituting the **ordinary** blocks
   `G_cc/G_cn` (that mis-notated form deviates O(1)–2.4e5 on the D1.5 grid);
   it is recorded as a **notation/implementation misinterpretation on our
   side, NOT a paper equation error**, and the O(1)-deviation conclusion is
   deleted. Production implements the single, algebraically equivalent Schur
   form (`CORE_UNCERTAINTY_METHOD=LEITE2024_EQ8_BERMANN2022_SCHUR`); no
   separate CLI variants; the direct Eq. (8) reference remains only inside
   the equivalence unit test (`test_apy_ginv_d15`).
   D=0 gate (no information): `D_c = D_n = 0 ⇒ S = 0 ⇒ C_core = (1/λ)G_cc`,
   i.e. `σ_e2·C_core = σ_a2·G_cc`, the exact no-data PEC limit (verified
   bit-wise on the D1.5 grid).
2. **σ²_e placement.** Eq. (8) as printed has no `σ²_e` factor; the exact
   counterpart `C_{u₂u₂} = σ²_e·(C_λ⁻¹)_{gg}` does. The λ-form MME conventions
   of C1 (`MME_SCALE`, `PEV_SCALE`) fix the correct placement:
   `Var(u_hat_core) = (G_w)_cc·σ_a2 − σ_e2·C_core_approx`, with the exact
   core slice `σ_e2·(C_λ⁻¹)_cc` used as benchmark (gate D3c quantification;
   z/P are **not** scale-invariant under an *additive* distortion of
   `G_ccσ²_u − C_cc`, so this matters for calibration, gate C2-G7/D6).

### 4.4. D is diagonal by construction (paper-stated noise source #2)

Eq. (9) drops off-diagonal elements of the fixed-effect absorption
(`I − X(X′X)⁻¹X′` term); the paper explicitly attributes part of
Approx_GinvAPY's extra noise to this ("off-diagonal elements are not
considered during the absorption of environmental effects into the mixed
model equations for constructing D"). No diagonal-only invention; this IS
the reference formula.

### 4.5. Core-size rule (paper-used; Pocrnic 2016)

The paper's core = **number of eigenvalues of G explaining 98% of the
genetic variance** (13,030 of 450,673 genotyped in the Angus data; most
livestock ≤15K per Pocrnic 2016/Cesarani 2022). Complexity is cubic in the
core size, linear in noncore.

Leite et al. (2024) used 13,030 core genotypes because this count explained
approximately 98% of the genetic variance in their ~450K-genotype matrix; core
animals were then randomly selected. The practical APY advantage therefore
depends on `n_g` being substantially larger than the genomic dimensionality
captured by the leading eigenvalues. A high `K98/n_g` in a small/medium D5
pilot (for example, `n_g=500`) is not, by itself, an algorithmic anomaly: it
means the benchmark regime has limited genomic redundancy at that sample size.
D5 records such cases as high-dimensional stress rather than tuning the data
or eigenvalues to make APY look better.

---

## 5. Current architecture map (inputs for the A22⁻¹/G⁻¹ audit)

| Object | Current construction | Storage | Where |
|---|---|---|---|
| `A⁻¹` | Meuwissen–Luo recursion | sparse, O(n) | `pedigree_relationship.cpp` |
| `A22` | Schur complement of sparse `A⁻¹` blocks: `(Ainv22 − Ainv21·Ainv11⁻¹·Ainv12)⁻¹`, dense inverse | dense n_g² (transient) | `PedigreeRelationship::buildA22` (pedigree_relationship.cpp:147) |
| `A22_inv_` | `invertSpdWithRidge(A22)` — second dense LDLT inverse of n_g×n_g | dense n_g² (persistent) | `computeA22` (single_step_model.cpp:244-245) → MME adapter (SplitMatrixAdapter), VCE H-sub |
| `G_` / `Ginv_bs_` / `Ginv_bs_tuned_` | GRM via `GrmBuilder` (float32 accumulation); dense LDLT inverses | dense n_g² each (persistent) | single_step_model.cpp:256, 283, 635-637; Ginv_bs_ feeds the MME adapter (line 926) |
| `H_inv_` | `A_inv_` + genotyped-block `τGinv − ωA22inv` | sparse with ~n_g² entries in the block (VCE path only) | `computeHInverse` (252-276) |
| MME (`C_λ`) | λ-form, `λ = Ve/Va`; PCG: ImplicitMME + SplitMatrixAdapter (matrix-free); LDLT: sparse lhs with dense genotyped block | factor fill grows with the dense block | `solveGeneticEvaluation` (885-1007) |
| `cinv_gg_` (C1 PEV) | n_g selected solves with SimplicialLDLT factor | dense n_g² | p-value exact step (1104-1140) |
| `vhat_g_` | `σa²·G_bs_ − σe²·cinv_gg_` | dense n_g² | line 1149 |
| SNP scan | `w_j = Ginv_bs_·M_j`; two quadratic forms vs `vhat_g_` | O(n_g) per SNP | lines 1219-1250 |

**Bottleneck accounting (exact C1, simultaneous peak):**
`A22_inv_ + G_ + Ginv_bs_ + Ginv_bs_tuned_ + cinv_gg_ + vhat_g_` ≈
**6·n_g² doubles = 48·n_g² bytes** + LDLT fill with dense-block nnz.
n_g = 50K ⇒ ≈120 GB (matches the paper's Exact_Ginv 159.66 GB scale,
their extra GBs from the full sparse-inverse machinery).
The two dense-inverse chains (`A22⁻¹` and `G⁻¹`) are each O(n_g³) time.
`A22⁻¹` **is** the other dense bottleneck the paper warns about (§23 of the
brief): it is currently materialized as a dense n_g×n_g matrix twice over
(A22 + its inverse) and kept, even though the MME only ever *applies* it
(matrix-vector products via the adapters, and, in exact mode, solves with it
inside `C_λ`).

**C2 consequences (mapped, not yet implemented):**
- approx-apy removes every dense n_g² object from the p-value path
  (`cinv_gg_`, `vhat_g_`, dense `G⁻¹`): nothing n_g×n_g remains — runtime
  guard required (gate C2-G5).
- `A22⁻¹` in the EBV solve: keep the sparse adapters; for PCG the required
  operation is `(A22⁻¹)·x`, available as forward/back-solves on the sparse
  `A⁻¹` Schur factors (reuse `buildA22` internals) — no dense `A22⁻¹`
  necessary. If that reuse proves too invasive this round, the limitation is
  stated explicitly (brief §24) rather than inventing an A22 approximation.
- EBV solution for approx-apy: PCG (already matrix-free) — never a hard
  LDLT fallback force because of p-values (brief §25).

---

## 6. The three modes (design, CLI proposed)

| Mode | G inverse | GEBV uncertainty | LHS inverse | Memory |
|---|---|---|---|---|
| `exact-ginv` (C1 frozen) | `G*⁻¹` dense | `σe²·cinv_gg_` (full n_g block) | n_g selected solves | O(n_g²) |
| `exact-apy` | `G_APY⁻¹` block-sparse | exact `σe²·(C_λ⁻¹)_cc` (core slice; n_g selected solves, kept n_c×n_c) | same | O(n_c² + n_g·nnz factor) — no dense n_g² G-inverse; Cinv slice n_c² only |
| `approx-apy` | `G_APY⁻¹` block-sparse | Eq. (8) single Schur form (`LEITE2024_EQ8_BERMANN2022_SCHUR`, gate D1.5) | **none** | O(n_c² + n_c·m + n_g) |

CLI (names per project style; final alpha after Phase B smoke):
- `--ssgwas-pvalue-method exact` (unchanged, backward-compatible)
- `--ssgwas-pvalue-method apy-exact`  (= Exact_GinvAPY)
- `--ssgwas-pvalue-method apy-approx` (= Approx_GinvAPY)
- `--ssgwas-apy-core-size N` | `--ssgwas-apy-core-file FILE` (core ID list)
- `--ssgwas-apy-seed S` (deterministic random core, baseline; never the
  only option)
- `--ssgwas-apy-core-auto` (eigenvalue-based target, §4.5; randomized
  eigensolver so initialization stays far below O(n_g³))
- `--ssgwas-apy-genomic-control` (explicit opt-in only; **off by default**;
  first release reports raw approximate P)

Metadata additions (all modes): `P_VALUE_METHOD` (`AGUILAR2019_APY_EXACT` /
`AGUILAR2019_APY_APPROX`), `G_INVERSE_METHOD` (`GINV_APY`), `APY_CORE_SIZE`,
`APY_NONCORE_SIZE`, `APY_CORE_SELECTION`, `APY_SEED`, `APY_CORE_VAR_FRACTION`,
`GEBV_UNCERTAINTY_METHOD` (`EXACT_LHS_CORE_SLICE` | `LEITE2024_EQ8` |
`BERMANN2022B_SCHUR`), `EXACT_OR_APPROX` (`APROX_GINV_APY` labels the whole
run; never per-SNP), `DENSE_NG2_ALLOCATED`, `ESTIMATED_MEMORY`,
`GENOMIC_CONTROL_APPLIED=false`, plus existing TIME_* series.

---

## 7. Complexity targets (brief §13)

| Quantity | Exact C1 | Exact GinvAPY | Approx GinvAPY |
|---|---|---|---|
| core-dependent | — (n_g is all) | n_c³ (Gcc⁻¹) | n_c³ (Gcc⁻¹) + n_c²·m (batched scan) |
| noncore-dependent | — | n_g·nnz(L) solves (Cinv) | linear in n_g (G_APY⁻¹ build, GEBV PCG) |
| marker-scan | O(n_g·m) | O(n_c·m) | O(n_c·m) |
| memory | 6·n_g²·8 B | ~n_g²·8 B (Cinv slice) + n_c² | n_c²·8 B + n_c·m·4 B + O(n_g) |

Paper anchors: Approx = 38× faster / 10× less memory than Exact at 50K;
24.47 h / 87.64 GB at 450K, linear in noncore. Memory and runtime must be
*measured* on CosmicScanner's own path in Phase E (gate C2-G10).

---

## 8. Gate plan for Phase C2 (mapped to Phases B–E; NOT yet executed)

| Gate | Content | Phase | Target |
|---|---|---|---|
| C2-G1 | APY G⁻¹ formula/reference verified | A **DONE** (§2, §3): Eq. 5-6 cross-checked vs Misztal 2014/2016 + Bermann 2022a Eq. 3; symmetry = 2×2 block formula; PD of Mnn diag by construction | keep open in B (code review) |
| C2-G2 | core_size = n_g ⇒ degenerate to exact G⁻¹ | B/C | n_g 30–100, exact-G⁻¹ vs APY with all-core: max|diff| ≤ ~1e-6; G_APY⁻¹ = G⁻¹ formal (empty noncore ⇒ G_APY⁻¹ = Gcc⁻¹) |
| C2-G3 | Exact_GinvAPY validated | C | tiny reference: z/P vs exact C1 ≤ C5 tolerance on 6 tags; Pearson(EA), Pearson(SE), Pearson(−log10P) ≥ 0.999 on a22/identity |
| C2-G4 | Approx core uncertainty matches reference equation | C/D **DONE (Phase D)**: single production Schur form == Leite Eq. 8 (superscript blocks) to 8.9e-16 (gate D1.5, 100 cases); D=0 limit exact; vs exact σe²(C_λ⁻¹)_cc quantified (gate D3c: SE Pearson ≈ 0.98, Z slope ≈ 1.05–1.07 at n_g=20/ml=21 records; remaining error = Eq. 9-diagonal + A22-structure, not APY) |
| D4 | uncertainty-error decomposition (tiny; diagnostic-only, not in CLI) | **DONE (2026-08-18)**: R0 exact == R1 apy-exact == R2 (core=a) block-identical; R1-vs-R2 structural (MME vs D-model, full-D) SE slope 0.9996 with Z 1.0545; R2-vs-R3 (full-D vs Eq.9 diag-D) SE slope 0.9640 with Z 1.0008; Eq.8/Eq.9 total (R1-vs-R3) SE 0.9616 Z 1.0554 LOG10P 1.0645 (== D3 corrected for ddof); R1≡R2 by construction (core slice of the actual-MEE inverse IS the exact core Schur complement); B==C gate: production Eq.9 D == diag(D_full) maxdiff 0; D_full off-diag = −1/n (intercept-only, analytic); intercept-only N-dependence: Z slope 1.077 (n=20) → 1.024 (n=1000), LOG10P 1.125 → 1.037 — converges toward 1.00, NOT a constant ≈1.10 mismatch; tau/omega four cases (tau=1): omega=0 (no A22⁻¹ correction) gives SE slope ≈0.52–0.58 / Z ≈1.76–1.87 / LOG10P ≈2.35–2.55; with omega=1 the pedigree-derived contribution on the typed block cancels and the genomic term aligns with G⁻¹ — mechanics: when the ungenotyped animals are eliminated from the pedigree precision A⁻¹ by Schur complementation, the effective pedigree precision of the typed animals is exactly A22⁻¹ (this is the definition of the A22⁻¹ block in single-step, not an approximation Ainv-GG ≈ A22⁻¹); in H⁻¹ = A⁻¹ + [0 0; 0 τG⁻¹ − ωA22⁻¹], the ωA22⁻¹ term subtracts that reduced typed-animal pedigree precision, so at ω=1 the pedigree contribution to the typed-animal reduced precision cancels and the genomic contribution is aligned with G⁻¹; the Eq.8 D-model (D + λG⁻¹, which carries no explicit A22⁻¹ term) is consistent with the MME only at ω=1; production omega=1 recovers SE ≈0.96 / Z ≈1.05–1.07; fixed-effect designs (intercept/+binary/+continuous/+categorical4) keep Z slope 1.042–1.055; operator full-vs-core (same PEV, Bermann Eq. 21) exactly 1.0 at core=n_g, 1.0284 at core=15 (SE); Python replica gates vs C++: exact/apy-exact/apy-approx SE match to <1e-6 for identity AND a22 at cores 20 and 15 |
| C2-G5 | no dense O(n_g²) uncertainty block in approx mode | B–E | runtime guard: `DENSE_NG2_ALLOCATED=false`, memory audit; test fails on dense Cinv_gg n_g² allocation |
| C2-G6 | Approx vs exact concordance quantified | D | n_g 1k–5k, M 5k–50k: Pearson/Spearman/slope/RMSE of EFFECT, SE, −log10P + top-k overlap (brief §17) |
| C2-G7 | null calibration quantified | E | C6-style per mode; raw λ_GC/FPR/QQ; NO auto-GC (report only) |
| C2-G8 | causal recovery quantified | E | C7-style per mode: rank, top-k, power, region (brief §19) |
| C2-G9 | core-size / core-selection sensitivity | E | ≥3 seeds × {25%, 50%, 100%-of-target} core sizes (brief §20-21) |
| C2-G10 | memory/runtime scaling demonstrated | E | small → medium → large synthetic; memory break of the n_g² curve; runtimes recorded |

Status label: **"Experimental scalable APY ssGWAS p-value"** (never
"Stable"). Current Phase D status after strict D6 null calibration is:
`APPROXIMATION_STATUS=EXPERIMENTAL_APPROXIMATE`,
`CALIBRATION_STATUS=LIMITED / REGIME_DEPENDENT`,
`STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`,
`SCALABILITY_SCOPE=UNCERTAINTY_LAYER_ONLY`, and
`GENOMIC_CONTROL_APPLIED=false`. The path must not be described as
simulation-validated or recommended for ranking/screening until D6.5 has
paired calibration boundaries and D7 has causal-ranking evidence. The exact C1
path status is unchanged.

Leite et al. (2024) reported Approx_GinvAPY QQ inflation (approximately
1.25±0.03) and discussed genomic control as a practical correction. CosmicScanner
currently reports raw uncorrected approximate P-values: no genomic control,
rescale, empirical correction, or result-dependent core selection is applied.

---

## 9. Frozen baseline and prohibitions (Phase C2)

Frozen exact baseline (must be bit-stable): `--single-step-method backsolve
--ssgwas-pvalue-method exact` (C1 gates C4–C7 re-verified during closeout).
Prohibited this round: weighted/iterative SNP ssGWAS; window P-values;
maternal/RRM modifications; automatic genomic control; core tuning to
maximize test correlations; any change to the exact C1 numbers; dense
n_g² blocks labeled "scalable"; commits; release-version edits.

---

## 10. Final report template (filled after Phases C/E)

- Q1. Exact C1 bottleneck: dense `G⁻¹`, `A22⁻¹`, `cinv_gg_`/`vhat_g_` —
       ≈6·n_g²·8 B peak (§5).
- Q2. APY G⁻¹ correctness: pending Phase B/C (G2-G3).
- Q3. core=n_g degeneration: pending (G2).
- Q4. Exact_GinvAPY − Exact_Ginv error source: APY G approximation
       (+ noncore ξ̂ term in the effect, §4.2).
- Q5. Approx_GinvAPY − Exact_GinvAPY error source: Eq. (8) PEV
       approximation (D diagonal absorption; §4.3/§4.4).
- Q6. Approx P null calibration: pending (G7, raw values).
- Q7. Genomic control needed?: report-only; off by default.
- Q8. Memory truly off O(n_g²)?: pending (G5/G10).
- Q9. A22⁻¹ new bottleneck?: reuse sparse Schur solves for PCG (limit if
       not feasible this round, §5).
- Q10. Recommended n_g: pending Phase E measurements.

Status: Phase A complete (this document). Next: **Phase B** on user
confirmation.
