# CosmicScanner v1.1.3 Release Test Coverage

This matrix binds each public focused-release mode to at least one smoke test.
`semi` appears only as a scan tier, not as a `--repeatability-mode` value.

| Model | CLI command pattern | Stability | Default | Smoke test name | Expected output | Required columns | Metadata checks |
|---|---|---|---|---|---|---|---|
| Single-trait LM | `--gwas --lm 4` | Stable | off | `single-trait LM all tests` | `.lm.assoc.txt` | `P_WALD,P_LRT,P_SCORE` | n/a |
| Single-trait LMM approx | `--gwas --scan-tier approx --no-loco` | Stable | approx | `single-trait LMM approx` | `.assoc.txt`, `.scan.meta.txt` | `P` | `release_status=stable; scan_tier=approx` |
| Single-trait LMM semi | `--gwas --scan-tier semi --no-loco` | Stable | approx | `single-trait LMM` | `.assoc.txt`, `.scan.meta.txt` | `P` | `release_status=stable; scan_tier=semi` |
| Single-trait LMM exact all tests | `--gwas --scan-tier exact --lmm-gemma-mode 4 --no-loco` | Stable | approx, legacy P if no mode | `README single-trait exact Wald/LRT/Score` | `.assoc.txt`, `.scan.meta.txt` | `P_WALD,P_LRT,P_SCORE` | `release_status=stable; scan_tier=exact; statistical_test=wald,lrt,score; loco=disabled_no_loco` |
| Multivariate EMMV | `--multivariate --mv-method emmv --no-loco` or `--scan-tier semi` | Stable | emmv | `multivariate LMM` | `.assoc.mv.txt`, `.scan.meta.txt` | `WALD_STAT,P,BETA_1,SE_1` | `release_status=stable; multivariate_method=emmv; scan_tier=semi` |
| Multivariate fixed-null Wald | `--multivariate --mv-method wald --no-loco` | Stable | emmv | `multivariate fixed-null Wald` | `.assoc.mv.txt`, `.scan.meta.txt` | `WALD_STAT,P` | `release_status=stable; multivariate_method=wald` |
| Multivariate GEMMA Wald | `--multivariate --mv-method gemma --no-loco` | Stable | emmv | `multivariate GEMMA Wald` | `.assoc.mv.txt`, `.scan.meta.txt` | `WALD_STAT,P` | `release_status=stable; multivariate_method=gemma` |
| Multivariate GEMMA full | `--multivariate --mv-method gemma-full --no-loco` | Stable | emmv | `multivariate GEMMA full` | `.assoc.mv.txt`, `.scan.meta.txt` | `P_WALD,P_LRT,P_SCORE` | `release_status=stable; multivariate_method=gemma-full; statistical_test=wald,lrt,score` |
| Multivariate GEMMA exact full | `--multivariate --mv-method gemma-exact-full --no-loco` | Stable | emmv | `README multivariate GEMMA exact full` | `.assoc.mv.txt`, `.scan.meta.txt` | `P_WALD,P_LRT,P_SCORE` | `release_status=stable; multivariate_method=gemma-exact-full; scan_tier=exact` |
| Repeatability exact mode | `--repeatability-mode exact --scan-tier semi --no-loco` | Stable | auto mode; approx tier if omitted | `repeatability exact mode, scan-tier semi` | `.assoc.txt`, `.repeatability.null.txt`, `.scan.meta.txt` | `P` | `release_status=stable; scan_tier=semi` |
| Repeatability fast mode | `--repeatability-mode fast --no-loco` | Stable | auto | `repeatability fast comparison` | `.assoc.txt` | `P` | output exists |
| Repeatability STCG mode | `--repeatability-mode stcg --no-loco` | Stable | auto | `repeatability STCG filtered smoke`; `critical_regressions` | `.stcg.assoc.txt` | `P` | output exists |
| A/G/H matrices | `--make-agh --ped <file> [--genotyped-list <file>]` | Stable | n/a | `A/G/H matrix construction` | `.Ainv.bin/.id`, `.G.bin/.id`, `.Hinv.bin/.id`, `.agh.meta.txt` | n/a | files exist |
| Multivariate EB | `--multivariate --mv-method eb --no-loco` | Experimental | off | `multivariate EB experimental` | `.assoc.mv.txt`, `.scan.meta.txt` | `WALD_STAT,P` | `release_status=experimental; multivariate_method=eb` |
| EMmvGWAS compatibility | `--emmvgwas --gwas` | Experimental | off | `multivariate EMmvGWAS compatibility experimental` | `.assoc.mv.txt`, `.scan.meta.txt` | `WALD_STAT,P` | `release_status=experimental; multivariate_method=emmv; compatibility_mode=emmvgwas` |
| Multivariate Hybrid | `--mv-method hybrid` | Unsupported | n/a | `bad_hybrid` | non-zero exit | n/a | error contains `unavailable in this focused release` |
| Multivariate GxE | `--mv-method gemma-gxe` | Unsupported | n/a | `bad_gxe`; `critical_regressions` | non-zero exit | n/a | error contains `unavailable in this focused release` |
