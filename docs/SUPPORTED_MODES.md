# CosmicScanner v1.1.3 Supported Modes

This focused release freezes the public interface by three criteria: the mode is
callable from the CLI, the code path is complete, and the bundled examples can
exercise the path successfully.

## Support Matrix

| Workflow | Mode | Status | Public command entry |
|---|---:|---|---|
| Single-trait LM | Wald, LRT, Score, all | Stable | `--lm 1/2/3/4` |
| Single-trait LMM approx | Fixed-null LMM scan | Stable | `--scan-tier approx` |
| Single-trait LMM semi | Fixed ratio, profile-scale scan | Stable | `--scan-tier semi` |
| Single-trait LMM exact | Per-SNP variance-component refit | Stable | `--scan-tier exact` |
| Single-trait LMM tests | Wald, LRT, Score, all | Stable | `--lmm-gemma-mode 1/2/3/4` |
| Null-model solver | Brent, IDUL, IDUL+, AIDUL | Stable | `--lmm-null-solver brent/idul/idul_plus/aidul` |
| Multivariate fixed-null | Wald | Stable | `--mv-method wald` or `--scan-tier approx` |
| Multivariate semi | EMMV/EMmvGWAS-style Wald | Stable | `--mv-method emmv` or `--scan-tier semi` |
| Multivariate GEMMA-style | Refit covariance, Wald primary; GEMMA equivalence validation pending | Experimental | `--mv-method gemma` |
| Multivariate GEMMA full | Refit covariance, Wald/LRT/Score columns; GEMMA equivalence validation pending | Experimental | `--mv-method gemma-full` |
| Multivariate exact full | GEMMA-style exact full refit, no-LOCO; GEMMA equivalence validation pending | Experimental | `--mv-method gemma-exact-full --no-loco` or `--scan-tier exact --no-loco` |
| Multivariate EB | Random-SNP EB Wald path | Experimental | `--mv-method eb` |
| EMmvGWAS compatibility | EMmvGWAS coding/no-LOCO compatibility | Experimental | `--emmvgwas` |
| Repeatability exact mode | Record-level GLS with approx/semi/exact scan tiers | Stable | `--repeatability-mode exact --scan-tier approx/semi/exact` |
| Repeatability fast | Woodbury-style accelerated scan | Stable | `--repeatability-mode fast` |
| Repeatability STCG | Matrix-free repeatability scan | Stable | `--repeatability-mode stcg` |
| A/G/H matrices | A-inverse, G, H-inverse | Stable | `--make-agh` |
| Multivariate Hybrid | Hybrid refinement route | Unsupported | rejected by `--mv-method hybrid` |
| Multivariate GxE | GxE routes | Unsupported | rejected by `--mv-method gemma-gxe` |
| Rare variant, biobank/REGENIE, BSLMM, DAP, summary-stat tools | Historical/internal fields only | Unsupported | no public focused CLI |

## CLI Semantics

- `--scan-tier approx`: fixed null variance components are estimated once and
  reused for all SNP tests.
- `--scan-tier semi`: null variance-component ratio is fixed; per-SNP residual
  scale is profiled.
- `--scan-tier exact`: per-SNP alternative variance components are refit.
- For repeatability, `exact`, `fast`, and `stcg` are `--repeatability-mode`
  values. `approx`, `semi`, and `exact` are scan tiers used only with
  `--repeatability-mode exact`; `semi` is not a repeatability mode.
- `--lmm-gemma-mode 1/2/3/4`: selects Wald, LRT, Score, or all three LMM test
  outputs. Mode 4 writes `P_WALD`, `P_LRT`, and `P_SCORE`; the primary `P`
  column remains Wald-oriented.
- `--lmm-null-solver` affects the null model only. Brent is the default. IDUL,
  IDUL+, and AIDUL can fall back to Brent if the null solver fails to converge.
- Default LOCO uses per-chromosome leave-one-chromosome-out GRMs when multiple
  chromosomes are present. `--no-loco` uses one global GRM for all scanned SNPs.
- Multivariate `--scan-tier approx/semi/exact` maps to `wald/emmv/gemma-exact`
  internally. For publication validation, use explicit `--mv-method` commands.

## Release Notes v1.1.3

- Public help now lists legal values, defaults, and scope for scan tiers,
  association tests, null-model solvers, multivariate methods, repeatability
  modes, and LOCO.
- Unsupported multivariate Hybrid and GxE modes now fail during argument parsing
  with a clear focused-release error.
- Scan metadata now records model, scan tier, statistical test, null solver,
  LOCO/no-LOCO policy, fallback policy, and multivariate method.
- README, user manuals, and validation scripts were aligned to this frozen
  support matrix.
