# PSR-GWAS v0.4 Parity Summary

PSR-GWAS v0.4 extends the v0.3 stable binary/binomial/Poisson implementation into a parity-ready validation surface covering the main trait classes claimed by PSR-GLMM/R: normal, binary, binomial, Poisson, and ordinal.

## Status

| Trait | Status | Entrypoint | Notes |
| --- | --- | --- | --- |
| normal | stable special case | standard `--gwas --no-loco` LMM | Gaussian identity-link traits reuse existing LMM GWAS. |
| binary | v0.3 stable | `--threshold-model --link probit|logit` | PSR/PQL pseudo-response plus weighted LMM scan. |
| binomial | v0.3 stable | `--count-model --family binomial --binomial-pheno` | Event/trial PSR/PQL. |
| poisson | v0.3 stable | `--count-model --family poisson` | Count PSR/PQL with optional offset. |
| ordinal | experimental v0.4 MVP | `--ordinal-model --ordinal-pheno` | Experimental cumulative-binary approximation; not full stacked proportional-odds PSR-GLMM parity. |

## Validation

Run:

```bash
bash tests/run_psr_full_parity_validation.sh
```

The wrapper runs stable CI, binary/binomial/Poisson validation, normal special-case equivalence, ordinal smoke validation, PSR-GLMM/R availability tracking, formal simulation, accuracy/speed benchmark, and benchmark summarization.

## Benchmark Outputs

- `tests/results/psr_gwas_smoke/psr_accuracy_speed_summary.tsv`
- `tests/results/psr_gwas_smoke/psr_accuracy_detail.tsv`
- `tests/results/psr_gwas_smoke/psr_runtime_summary.tsv`
- `tests/results/psr_gwas_smoke/psr_method_availability.tsv`
- `tests/results/psr_gwas_smoke/psr_gwas_parity_matrix.tsv`
- `tests/results/psr_gwas_smoke/psr_accuracy_speed_report.md`
- `tests/results/psr_formal_simulation/psr_formal_simulation_report.md`
- `tests/results/psr_ordinal_full_design/full_ordinal_design_check.tsv`

## v0.4.1 Additions

- Formal simulation (`tests/run_psr_formal_simulation.py`) expanded with scenarios for binary null_rare, link_logit, poisson null_offset, binomial null_rare_prop/null_balanced_prop, and generates `psr_formal_simulation_report.md`.
- Accuracy/speed benchmark (`tests/run_psr_accuracy_speed_benchmark.py`) now generates `psr_accuracy_speed_report.md` with per-trait accuracy, runtime, and method availability.
- CLI: `--ordinal-mode experimental-cumulative-binary|stacked-psr` and `--ordinal-allow-fallback` added.
- Full ordinal PSR design document: `docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md` (status: `design_only`).
- Full ordinal design check: `tests/run_psr_full_ordinal_design_check.sh` reports `full_ordinal_psr_status`.

## External References

The official PSR-GLMM/R URL is tracked by `tests/check_psr_glmm_r_availability.sh`. Unavailability is recorded but does not fail CI. Source code is not embedded in CosmicScanner.

## Scope Limits

PSR-GWAS remains a PSR/PQL approximation, not exact GLMM inference. Repeated-measure PSR-GLMM and random-regression PSR-GWAS remain roadmap items. Negative binomial remains experimental and should not be treated as stable parity output. Stacked proportional-odds PSR-GWAS is design-only.
