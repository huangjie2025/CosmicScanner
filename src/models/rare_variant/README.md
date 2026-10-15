# Rare-Variant Extensions (R2)

## Integration API

Include `models/rare_variant/rare_variant_extensions.h` and link
`src/models/rare_variant/rare_variant_extensions.cpp`. The production build
links this source explicitly and exposes it through `--skat-o` and `--acat-v`.

All pure functions are in `cosmic::rare_variant`:

```cpp
MixtureTailResult mixtureChiSquareUpperTail(
    double q, const Eigen::VectorXd& lambda);

SkatOGridResult skatOConservativeGrid(
    const Eigen::VectorXd& scores,
    const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights,
    const std::vector<double>& rho_grid =
        {0.0, 0.01, 0.04, 0.09, 0.16, 0.25, 0.5, 1.0});

SkatOResult skatOCanonical(
    const Eigen::VectorXd& scores,
    const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights);

AcatVResult acatV(
    const Eigen::VectorXd& p_values, const Eigen::VectorXd& weights);

AcatVResult acatVFromScores(
    const Eigen::VectorXd& scores,
    const Eigen::MatrixXd& covariance,
    const Eigen::VectorXd& weights);
```

The existing `RareVariantTest::runTest(set_id, X_set, scores, P_matrix, MAF)`
now dispatches `Method::SKAT_O` and `Method::ACAT_V` to these kernels, computing
`covariance = X_set.transpose() * P_matrix * X_set`. Existing `BURDEN` and `SKAT`
algorithms are unchanged. `Method::SKAT_O` no longer silently falls back to SKAT.

`RareVariantTest::Result` adds `method_name`, `selected_rho`, and `min_p_value`.
Use `p_value`, NOT `min_p_value`, as the set-level p-value. The wrapper uses the
default fixed rho grid and existing Beta(MAF; a,b) score weights. Its ACAT-V
combination weights are Beta(MAF; a,b)^2 * MAF * (1-MAF). Direct pure-function
calls support explicitly supplied weights and grids without CLI changes here.

## Statistical Contract

`Method::SKAT_O` uses `skatOCanonical` and reports `SKAT_O_CANONICAL`. It evaluates
the fixed rho grid under the shared Gaussian score null with deterministic Halton
quasi-Monte Carlo, then integrates the union event at the observed minimum
marginal p-value. The output records `P_METHOD=DEPENDENCE_ADJUSTED_GAUSSIAN_QMC`
and whether a fallback was used. `skatOConservativeGrid` remains available as an
explicit legacy helper and retains `SKAT_O_EXPERIMENTAL_CONSERVATIVE_GRID`.

For U approximately distributed as N(0,C), and W=diag(w), the kernel is

```
Q_rho = (1-rho) * sum((w_i U_i)^2) + rho * (sum(w_i U_i))^2
A_rho = (1-rho) * W^2 + rho * w w'
null(Q_rho) = sum(lambda_i(rho) * chi-square_1)
```

Eigenvalues come from a covariance-root transform, not from marginal variances
alone. All rho values must be fixed independently of the observed scores.
The grid is sorted and exact duplicate rho values removed. Canonical output is
dependence-adjusted rather than `m * min(p)`. Single-variant and singleton-grid
tests reduce exactly to their marginal score test. The legacy conservative helper
continues to report `min(1, m * min_rho(p_rho_upper))`.

The mixture evaluator uses exact chi-square survival recurrences for rank-one
and equal spectra. For unequal positive eigenvalues, beta=min(lambda) yields
Ruben's positive gamma-mixture expansion:

```
c_0 = product(sqrt(beta / lambda_i))
d_k = 0.5 * sum((1 - beta / lambda_i)^k)
c_j = (1/j) * sum_{k=1}^j(d_k * c_{j-k})
P(Q >= q) = sum_{j>=0}(c_j * P(chi-square_(rank+2j) >= q/beta))
```

The coefficients are nonnegative and sum to one. The finite partial sum plus
remaining coefficient mass bounds the omitted contribution from above. The
implementation uses at most 512 terms and takes the smaller of that upper bound
and an optimized Chernoff/MGF upper bound. It does not reuse the legacy
Satterthwaite approximation, which is not generally a conservative tail bound.
Highly unequal spectra and extreme tails may yield loose bounds and lower power.
These are mathematical bounds evaluated with floating-point arithmetic and
roundoff padding, not certified interval arithmetic. Eigenvalue computations
also have ordinary floating-point precision limits. The pure helper exposes
`conservative_bound=true` for unequal spectra.

For score APIs, weights must be finite and nonnegative with positive total weight
on nonempty inputs. Covariance must be finite, symmetric and PSD. Roundoff-sized
negative eigenvalues are clipped; positive modes are retained even when tiny.
Scores outside the numerical covariance support are rejected.
For the grid statistic, tolerated nullspace rounding residuals are projected out
so that the statistic and its null kernel use the same covariance support.
Empty sets and zero-score/zero-covariance sets are neutral. Validation errors throw
`std::invalid_argument`; numerical failures/overflow throw runtime/overflow errors.

For grid tests, input weights are score multipliers, not their squares. They are
normalized to maximum one, leaving p-values unchanged. Returned Q statistics use
those normalized weights in the original score units. `n_variants` counts positive
weight, positive-variance variants. Covariance must correspond to the supplied
scores under the same fitted null model and include all relevant residual-scale
factors. Estimated-null Gaussian score calibration is asymptotic, not an exact
small-sample or rare-binary-trait guarantee.

ACAT-V uses `sum(normalized_weight_i * cot(pi * p_i))` and the standard Cauchy
reference tail. Its input weights are direct combination weights, not squared
score multipliers. The implementation scales cotangents before division and uses
`atan2` for the final tail, preserving tiny and subnormal p-values without
evaluating `tan(pi*(0.5-p))` or subtracting a tail from one.

ACAT's Cauchy reference is exact for independent continuous uniform p-values;
under general LD it is a small-tail approximation, not a universal conservative
bound. This implementation combines individual-variant p-values only and does
not implement an ultra-rare-variant burden-pooling rule. Single-variant score
p-values use `erfc(abs(U_i)/sqrt(2*C_ii))`, omitting zero-variance variants.

Direct ACAT accepts p in [0,1]. A positive-weight exact zero produces p=0, an
exact one produces p=1, and simultaneous active zeros and ones are undefined and
throw. Zero-weight entries are ignored after validation. Computed positive tails
are floored only at the smallest positive double on numerical underflow, never
at a conventional threshold such as 1e-15. The Cauchy statistic can be infinite
at endpoints or overflow when the wrapper converts it to double; the p-value is
still valid under the stated boundary convention. Do not interpret or label that
statistic as chi-square. In particular, an exact p=1 is not a neutral ACAT input.

## Scanner Output Handoff

`gwas_scanner_rv.cpp` writes `<out>.rare.txt` directly from each computed set
result with this schema:

```text
CHR START END VARIANT_SET_ID METHOD STATISTIC P_VALUE
```

No fabricated effect sizes or standard errors are written. SKAT-O and ACAT-V
retain the additional `<out>.rare.meta.tsv` audit fields, including rho-grid and
P-value method diagnostics where applicable.

Filter monomorphic variants before wrapper ACAT calls, because all-MAF-zero sets
have zero total ACAT weight and are intentionally rejected. Handle empty sets
explicitly and do not convert invalid covariance errors into significant output.

## Focused Verification

The standalone test requires only Eigen and the old/new rare-variant translation
units; no whole-project build is needed:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -I src \
  -isystem /disk228/huangjie/software/miniconda3/include/eigen3 \
  tests/test_rare_variant_extensions.cpp \
  src/models/rare_variant/rare_variant_extensions.cpp src/rare_variant_test.cpp \
  -o /tmp/opencode/test_rare_variant_extensions
/tmp/opencode/test_rare_variant_extensions
```

Tests cover analytic rank-one/equal-spectrum and unequal-spectrum convolution
tails, ill-conditioned bounds, rho endpoints, correction and deduplication,
weight/score scaling, LD and singular covariance, ACAT analytic identities,
subnormal tails, invalid and empty input, wrapper dispatch, and a fixed-seed
Gaussian-null regression guard. They do not establish equivalence to canonical
SKAT-O, validate scanner I/O, or replace binary-trait/LD/power benchmarking.
