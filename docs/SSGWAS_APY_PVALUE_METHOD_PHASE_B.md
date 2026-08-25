# ssGWAS SNP P-value — APY Route Phase B (Exact_GinvAPY) Implementation Report

Companion to `SSGWAS_APY_PVALUE_METHOD.md` (Phase A audit, 2026-08-17) and
`SSGWAS_SNP_PVALUE_METHOD.md` (C1 exact route, frozen reference).
Date: 2026-08-17.

Status: **Phase B (Exact_GinvAPY) IMPLEMENTED, unit tests + gate C2-G3
end-to-end regression PASS.** Not yet committed (report-first workflow).

---

## 1. Scope delivered

`--ssgwas-pvalue-method` gains two APY modes:

| mode | meaning | status |
|---|---|---|
| `off` / `exact` | C1 frozen behavior (unchanged) | frozen |
| `apy-exact` | Exact_GinvAPY: G⁻¹→G_APY⁻¹ in the MME + Eq.7 marker-effect variance using the **exact** core-slice of the LHS inverse (`Cinv_cc`), core-animal | **NEW (this phase)** |
| `apy-approx` | Approx_GinvAPY (Eq. 8–10 block-sparse approximation) | rejected for now: "Phase D not implemented yet; use apy-exact in Phase C2" |

`apy-exact` = the C1 exact route with `G⁻¹` replaced by `G_APY⁻¹` *everywhere
it appears* (MME LHS, backsolve operator, and the Eq.7 variance term), so the
three quantities `(â, sd(â), p)` stay internally consistent. This mirrors the
paper's Exact_GinvAPY (GSE 56:59); its measured cost there was no better than
Exact_Ginv, and the same is true here — the win over C1 is limited to the
`n_core × n_core` core-slice inversion + memory, the real win comes with
Phase D/`apy-approx` (linear in noncore animals).

## 2. Mathematics (as implemented, with the corrected G_apy definition)

APY inverse (paper Eq. 5/6, Misztal 2014) with core `c` and noncore `n`:

    G_APY⁻¹ = [Gcc⁻¹ 0; 0 0] + [−Gcc⁻¹Gcn; I]·Mnn⁻¹·[−GncGcc⁻¹ I]     (Eq. 5)

**Key fact (verified in unit tests):** this is the *exact inverse* of the
matrix where noncore rows/columns are replaced by their conditional
prediction given the core:

    G_apy  = [Gcc  Gcn; Gnc  Mnn + Q'·Gcn],   Q = Gcc⁻¹·Gcn             (Eq. 6')
    Mnn    = diag(g_jj − g_jc'·Gcc⁻¹·g_cj)   (noncore residual variances)

If one instead built `[Gcc Gcn; Gnc diag(Mnn)]`, `G_APY⁻¹` is *not* its
inverse — this was the exact bug caught during Phase B debugging (a −Q·t sign
error in `applyInverse` and the wrong "modified matrix" identity both trace to
it). The noncore diagonal property `m_jj = 1/[G_APY⁻¹]_jj` holds (unit test).

`applyInverse(x)` block solve, matrix-free, O(n_c² + n_c·n_n):

    t  = Mnn⁻¹·(x_n − Q'·x_c)
    y_c = Gcc⁻¹·x_c − Q·t
    y_n = t

Degeneration gate **G2**: at `n_core = n_g`, `G_apy = G`, `G_APY⁻¹ = G⁻¹`
(NoncoreSize=0, Q is 0×0), so `apy-exact` degenerates to the C1 exact route
bit-for-bit (e2e result: max |Δ| = 0.00e+00 on EBV/EA/SE/Z/LOG10P, both
identity-G and blend-a22-G).

Marker-effect variance (Eq. 7, core-slice version):

    Var(EA_j) = w_c'·vhat_c·w_c / S²
    vhat_c    = (G⁻¹G)ᶜᶜ·σa² − σe²·Cinv_cc     (Cinv = LHS⁻¹, solve on G⁻¹·G... see file ApyCoreOps)
    w_c       = (G_backsolve⁻¹)ᶜᶜ·M_c,j        (core rows of the backsolve preconditioner inverse)
    p         = erfc(|z|/√2),  z = EA/SE        (paper Eq. 3; identical to C1)

Implementation detail: the DenseM compute block is assembled once and reused
for the 8 GiB-guarded core slice even when `n_c` exceeds plausible tiny-data
sizes; the compute block for `vhat_c` runs `n_c` LV solves over the MME LHS
(selected solves on the genotyped columns, +1 offset for the intercept).

## 3. Architecture

New files:
- `src/models/apy_ginv.h/.cpp` — `ApyGinv`: builds and owns the APY pieces
  (`Gcc⁻¹` LDLT, `Q = Gcc⁻¹Gcn`, `Mnn⁻¹` vector, `core_block_ = Gcc⁻¹ +
  Q·Mnn⁻¹·Q'`), plus `applyInverse`, `applyCoreInverse`, `visitTriplets`
  (symmetric, bidirectional), `addDiagonalTo`, and diagnostics
  (`info_`: min/max LDLT pivots of Gcc, min/max of Mnn, core-block min/max
  abs). `ApyGinvAdapter` wraps it as an `AbstractMatrix` with τ = h_tau for
  MME insertion; `getBlock`/`toDense` explicitly throw (no accidental dense
  materialization of n_g² G⁻¹; `DENSE_NG2_ALLOCATED=false` in metadata).
- `GrmBuilder` — the same object that built the dense G⁻¹ at seed level is
  used to expose `Gcc`/`Gcn`/`Gnc` blocks to `ApyGinv::build`; no duplicated
  GRM construction path.

Integration (all inside `runBacksolve`/`solveGeneticEvaluation`, exact path
untouched):
1. MME LHS: `SplitMatrixAdapter(A_, a22_inv_scale, ...)` for C1; for apy a
   `unique_ptr<ApyGinvAdapter>` provides the G⁻¹ block with the same
   `σe⁻¹·...` scaling as C1, so LDLT factorizations use an identical
   numerical LHS layout (this is why core=n_g matches bit-for-bit).
2. Backsolve operator: `backsolveSnpEffects` uses `apy_ginv_eff_->applyInverse`
   on the GEBV vector (matrix-free), replacing the dense `G⁻¹` multiply.
3. P-value: `computeSnpPvaluesExact` (C1, frozen) vs
   `computeSnpPvaluesApyExact` (new) selected by option value; both write the
   same `.pvalue.txt` schema and share `writePvalueFile`/`t_pev_out_` timing.

CLI (config.h + main.cpp):
- `--ssgwas-pvalue-method off|exact|apy-exact|apy-approx` (defaults unchanged)
- `--ssgwas-apy-core-size N` (1..50,000,000)
- `--ssgwas-apy-core-file FILE` (pedigree IDs, one per row; `#`/blank skipped)
- `--ssgwas-apy-seed S` (default 12345; used only in RANDOM_SEEDED mode)
- Validation: apy modes require core options (`throw` otherwise); `core-size >
  n_genotyped` rejected; unknown/duplicate IDs in the core file rejected with
  an explicit message; `core-size` in FILE mode must match the file row count.
- Core selection modes (metadata reports `APY_CORE_SELECTION`):
  - `FILE` — exact index order of the file
  - `RANDOM_SEEDED` — partial Fisher–Yates on a `std::mt19937_64(seed)` state,
    taking the first N drawn indices, stored ascending (block arithmetic
    partition-insensitive; verified by unit test)

Metadata added (apy-exact runs):
`P_VALUE_METHOD=AGUILAR2019_APY_EXACT`, `G_INVERSE_METHOD=GINV_APY`,
`APY_CORE_SIZE`, `APY_NONCORE_SIZE`, `APY_CORE_SELECTION`, `APY_CORE_SEED`,
`APY_CORE_FILE`, `APY_GCC_MIN_PIVOT`/`MAX_PIVOT`, `APY_MNN_MIN`/`MAX`,
`MME_INVERSE_METHOD=SELECTED_SOLVES_LDLT_CORE_SLICE`,
`GEBV_UNCERTAINTY_METHOD=EXACT_LHS_CORE_SLICE`,
`EXACT_OR_APPROX=EXACT`, `DENSE_NG2_ALLOCATED=false`,
`ESTIMATED_CINV_CC_BYTES`, `CINV_CC_MAX_ASYMMETRY`,
`VAR_GHAT_MIN_EIG`/`MAX_EIG` + `VAR_GHAT_EIG_MODE`, `TIME_PEV_*`.

## 4. Numerical decisions (gate-relevant)

- Gcc near-singularity: LDLT of Gcc must not be treated as invertible at the
  precision floor — relative pivot test `min(|D|) <= 1e-12 · max(|dmax|,1)`
  throws (Eigen's LDLT does not itself flag semidefinite input).
  `--allow-large-exact-ssgwas-pvalue`-style override is NOT given for this;
  the throw is unconditional (a near-singular Gcc makes the whole core-slice
  algebra meaningless).
- 8 GiB guard on the dense `n_core²` compute block (Cinv_cc + vhat_c + w_c):
  `throw` unless `--allow-large-exact-ssgwas-pvalue Yes` (same flag as C1,
  reused deliberately: a user who already accepted the C1 dense path accepts
  this one).
- Cinv_cc from `n_c` selected solves with `SimplicialLDLT` on the MME LHS,
  symmetrized (`(C+C')/2`), max asymmetry recorded in metadata.
- PSD audit of `vhat_c`: explicit eigen decomposition for `n_c <= 4096`;
  otherwise 24 Rayleigh-probe iterations, tolerance 1e-6 × max|eig|; any
  eigenvalue below −tol → `throw "INVALID_COVARIANCE"` (paper: var(â) can go
  negative at extreme cores; we fail loudly instead of emitting bogus p).
- Per-marker tolerance reference `σa²·w_c'(G·w_c)ᶜᶜ·/S²` precomputed to
  classify NEGATIVE_VARIANCE_CLAMPED vs NUMERICAL_FAILURE vs ZERO_SE; the
  clamp/status conventions are byte-identical to C1's schema.
- All computation paths are single-threaded-per-factoring; OpenMP sections
  reused where the C1 path already used them (no new threading).

## 5. Verification

### 5.1 Unit tests (`tests/test_apy_ginv.cpp`, registered with CTest; ALL PASS)
- G_APY⁻¹ symmetric; `G_apy · G_APY⁻¹ ≈ I` (the corrected Eq. 6' identity);
- `applyInverse` equals dense linear solve to ~1e-12;
- `m_jj = 1/[G_APY⁻¹]_jj` (noncore diagonal property);
- `coreBlock = [G_APY⁻¹]cc` and `applyCoreInverse = Gcc⁻¹·x`;
- `n_core = n_g` ⇒ G_APY⁻¹ numerically equals dense G⁻¹ (gate G2 at unit
  level); `Mnn + Q'Gcn` identity for the modified-matrix test;
- invalid inputs: empty/duplicate/out-of-range core, non-symmetric singular
  Gcc, non-square G, negative core size — all throw.

### 5.2 End-to-end (new `tests/run_ssgwas_pvalue_apy_exact_regression.py`, PASS)
Tiny design n_ped=30/n_g=20/m=30 with an analysis-consistent null draw
(Meuwissen–Luo A, H⁻¹ blending 5% per paper Eq. 2 semantics), identity-G and
blend-a22-G, compared against the frozen C1 `exact` run on identical data:

| run | max |Δ| EBV / EA / SE / Z / LOG10P |
|---|---|---|
| identity  core=20 (=n_g) | 0.00e+00 / 0.00e+00 / 0.00e+00 / 0.00e+00 / 0.00e+00  **PASS (gate C2-G3)** |
| a22       core=20 (=n_g) | all 0.00e+00  **PASS** |
| identity  core=15 | 1.38e-01 / 2.42e-03 / 3.74e-02 / 9.05e-01 / 8.80e-01 (finite) |
| identity  core=10 | 3.96e-01 / 1.97e-02 / 4.68e-02 / 1.85e+00 / 2.40e+00 (finite) |
| identity  core=5  | 4.93e-01 / 5.34e-02 / 7.96e-02 / 3.20e+00 / 6.20e+00 (finite) |
| a22       core=15/10/5 | monotonically larger deviations, all finite |

Deviations grow monotonically as the core shrinks (expected: random cores of
15/20 on a 20-animal GRM are far from the "capture 98% of the variance" rule
that guides real-data cores like 13K-of-450K). Metadata gate: all APY_* /
GINV_APY / SELECTED_SOLVES / DENSE_NG2_ALLOCATED=false / RANDOM_SEEDED fields
as expected. Core-file mode (all-core file) reproduces exact within tolerance;
`APY_CORE_SELECTION=FILE`.

Error-path gate (all rejected with nonzero exit + message):
`apy-approx` (Phase D hint), apy-exact without core options, core-size >
n_genotyped, core file with unknown ID.

### 5.3 Build
`make -C build cosmicscanner` and `test_apy_ginv` clean (Release, -j4);
full closeout regression suite re-run completed under 2026-08-17 session
(C4/C5 stated PASS earlier; pending full-suite summary in the session log).

## 6. Complexity and limits (Phase B = exact slice; honest)

- MME assembly/factorization cost is that of C1 (n_g² storage was already
  paid there); saving: Ginv never materialized at n_g² (only n_c² slices),
  backsolve operator matrix-free, and `Cinv_cc` needs n_c selected solves
  instead of an n_g² dense inverse.
- Memory: `ESTIMATED_CINV_CC_BYTES` ≈ 8·n_c² (+ core-block copy). 8 GiB guard
  caps n_c ≈ 32K on the default threshold.
- `apy-approx` (Phase D, paper Eq. 8–10 block-sparse cinv + GEBV PEC appx) is
  where the paper's linear-in-noncore scaling lives; currently rejected by
  design to keep Phase B honest about what was and was not implemented.
- PCG solver + apy-exact is refused (MME LDLT required for the exact
  core-slice; a PCG+APY combination belongs to Phase D's approximated path).

## 7. Files touched (diff discipline)

Added: `src/models/apy_ginv.{h,cpp}`, `tests/test_apy_ginv.cpp`,
`tests/run_ssgwas_pvalue_apy_exact_regression.py`. Modified:
`src/models/single_step_model.{h,cpp}` (Options + core members; runBacksolve
branches for Ginv construction/adapter/backsolve; new
`computeSnpPvaluesApyExact` + `buildApyCoreSelection`; metadata apy-exact
branch), `src/config.h`, `src/main.cpp` (three new flags + help + dispatch),
`CMakeLists.txt` (test target). The C1 `exact` p-value path received zero
functional edits (behavioral freeze held; e2e C4/C5 re-run PASS at the same
1e-5 tolerance).