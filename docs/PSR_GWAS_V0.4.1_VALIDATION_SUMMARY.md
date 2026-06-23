# PSR-GWAS v0.4.1 Validation Summary

## Status: Reportable Validation Infrastructure

v0.4.1 extends the v0.4 parity-ready surface with reportable formal simulation, accuracy/speed benchmarks, ordinal design documentation, and CLI scaffolding.

## What's New

| Area | Description |
| --- | --- |
| Formal simulation | `tests/run_psr_formal_simulation.py` expanded with null_rare, link_logit, null_offset, null_rare_prop, null_balanced_prop scenarios; outputs `link`, `ordinal_mode`, `has_offset`, `error_message`; generates `psr_formal_simulation_report.md`. |
| Accuracy/speed benchmark | `tests/run_psr_accuracy_speed_benchmark.py` now generates `psr_accuracy_speed_report.md` with per-trait accuracy, runtime, method availability. |
| Ordinal CLI | `--ordinal-mode` (`experimental-cumulative-binary`/`stacked-psr`) and `--ordinal-allow-fallback` added to `config.h`, `main.cpp`, `ordinal_model.h`, `ordinal_model.cpp`. |
| Ordinal design doc | `docs/PSR_GWAS_FULL_ORDINAL_DESIGN.md` documents stacked proportional-odds PSR with C-1 thresholds, shared beta, stacked pseudo-response, W, K_stack, V_stack. Status: `design_only`. |
| Ordinal design check | `tests/run_psr_full_ordinal_design_check.sh` confirms design doc presence and reports `full_ordinal_psr_status=design_only`. |
| Full parity wrapper | `tests/run_psr_full_parity_validation.sh` includes formal simulation, accuracy/speed report, ordinal design check. |
| `.gitignore` | Added to ignore `build/`, `tests/__pycache__/`, `tests/results/`. |
| Tag | `psr-gwas-v0.4-parity-ready` frozen at v0.4 stable point. |

## CLI Changes

```
--ordinal-mode experimental-cumulative-binary|stacked-psr
--ordinal-allow-fallback
```

- `experimental-cumulative-binary` (default): runs C-1 binary PSR scans per threshold; picks best p-value per SNP.
- `stacked-psr`: design-only; raises `std::runtime_error` unless `--ordinal-allow-fallback` is set, which falls back to experimental cumulative-binary with a warning.

## Key Results

- Normal LMM log10p_corr = 1.0 (special case).
- Binary PSR log10p_corr ≥ 0.9996 vs R glm.
- Binomial PSR log10p_corr ≥ 0.9996 vs R glm cbind.
- Poisson PSR log10p_corr ≥ 0.971 vs R glm (PSR/PQL approximation).
- Ordinal: cumulative-binary path validated; stacked-psr design-only.
- R availability: `lm`, `glm`, `MASS::polr` available; `logistf`, `lme4`, `glmmTMB`, `ordinal`, `POLMM` unavailable (noted but non-blocking).

## Validation Commands

```bash
bash tests/run_psr_full_parity_validation.sh
PSR_FORMAL_REPS=1 python3 tests/run_psr_formal_simulation.py   # sanity
python3 tests/run_psr_accuracy_speed_benchmark.py              # full benchmark
bash tests/run_psr_full_ordinal_design_check.sh                # design check
```

## Outputs

- `tests/results/psr_formal_simulation/psr_formal_simulation_report.md`
- `tests/results/psr_gwas_smoke/psr_accuracy_speed_report.md`
- `tests/results/psr_gwas_smoke/psr_accuracy_speed_summary.tsv`
- `tests/results/psr_gwas_smoke/psr_method_availability.tsv`
- `tests/results/psr_gwas_smoke/psr_gwas_parity_matrix.tsv`
- `tests/results/psr_ordinal_full_design/full_ordinal_design_check.tsv`

## Known Limitations

- Stacked proportional-odds PSR is Stage 1 (fixed-VC weighted GLS first pass; see `docs/PSR_GWAS_V0.5_STACKED_ORDINAL_SUMMARY.md`). Stage 2 (dense ZKZ' random effects) is not yet started.
- `ordinal::clm` and `POLMM` are unavailable; no ordinal R baseline for accuracy comparison.
- Negative binomial remains experimental.
- Exact GLMM, repeated-measure, and random-regression PSR remain future roadmap items.
- Ordinal cumulative-binary path uses experiment-wise MIN p-value (best across C-1 thresholds) without multiplicity adjustment. Full stacked model would provide a single joint test.
