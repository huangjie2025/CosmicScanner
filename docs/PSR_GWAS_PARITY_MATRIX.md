# PSR-GWAS Parity Matrix

CosmicScanner PSR-GWAS v0.4 tracks parity against the trait classes claimed by the PSR-GLMM/R paper tool while keeping PSR-GLMM/R itself as a non-blocking external reference. The official PSR-GLMM/R site is availability-tracked but not required for CI because it has not been reliably reachable.

| Trait type | PSR-GLMM/R claimed support | CosmicScanner status | Entrypoint | Validation | Benchmark | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| normal | yes | implemented special case | standard LMM GWAS via `--gwas --no-loco` | v0.4 normal special-case validation | included | Gaussian identity-link case reuses existing LMM GWAS; no separate PSR iteration is needed. |
| binary | yes | implemented v0.3 stable | `--threshold-model --link probit|logit` | v0.3 stable CI | included | PSR/PQL pseudo-response plus weighted LMM scan. |
| binomial | yes | implemented v0.3 stable | `--count-model --family binomial --binomial-pheno` | v0.3 stable CI | included | Event/trial PSR/PQL with bounded means and working weights. |
| poisson | yes | implemented v0.3 stable | `--count-model --family poisson` | v0.3 stable CI | included | Poisson PSR/PQL with optional offsets. |
| ordinal | yes | experimental v0.4 MVP | `--ordinal-model --ordinal-pheno` | v0.4 smoke validation | included | Current implementation is an experimental cumulative-binary approximation, not full stacked proportional-odds PSR-GLMM parity. |

Machine-readable copy: `tests/psr_gwas_parity_matrix.tsv`.

Repeated-measure PSR-GLMM and random-regression PSR-GWAS remain roadmap items. Negative binomial remains experimental and should not be reported as a stable parity trait.
