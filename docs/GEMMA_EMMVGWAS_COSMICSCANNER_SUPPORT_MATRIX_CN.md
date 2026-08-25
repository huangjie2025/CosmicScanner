# GEMMA / EMmvGWAS / CosmicScanner 功能支撑矩阵

本文档用于验收 CosmicScanner 对 GEMMA 与 EMmvGWAS 功能路线的覆盖情况。状态定义如下：

- 完全支持：CosmicScanner 有原生实现，并已纳入本地或 reference 验证。
- 兼容支持：命令、输入输出和主要统计路线已覆盖，可用于 GEMMA/EMmvGWAS 风格工作流。
- 近似支持：提供 GEMMA-style 文件接口和确定性近似实现，但不等价于 GEMMA 原始采样器或内部优化器。
- 条件验证：功能已实现，但外部 reference 验证依赖 GEMMA、R、EMmvGWAS 包或可开 socket 的运行环境。
- 不适用：该软件本身不提供该类模型或工具路线。

## 总体结论

CosmicScanner 当前覆盖 GEMMA 的主干 GWAS、mvLMM、summary-stat、VC、GRM/eigen、calccor、prediction 和 BSLMM/DAP 文件接口；同时覆盖 EMmvGWAS 的 `vc_estimation`、`full_gwas` 和 `semi_gwas` 三条路线。除 GEMMA BSLMM MCMC/DAP 内核属于确定性近似外，GEMMA 和 EMmvGWAS 的核心 GWAS 路线已经可以在 CosmicScanner 中直接执行。

## 三方功能矩阵

| 功能路线 | GEMMA | EMmvGWAS | CosmicScanner | CosmicScanner 验证状态 |
|---|---|---|---|---|
| PLINK 1 binary 输入 | 支持 | 不直接主打 | 完全支持 | smoke/reference 覆盖 |
| PLINK 2 PGEN 输入 | 不主打 | 不适用 | 完全支持 | 核心回归覆盖 |
| BGEN 输入 | 不主打 | 不适用 | 完全支持 | 核心回归覆盖 |
| GRM 构建 `-gk` | 支持 | 内部使用 K | 完全支持 | smoke 覆盖 |
| GEMMA eigen `-eigen` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| 单性状 LM `-lm 1` Wald | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LM `-lm 2` LRT | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LM `-lm 3` Score | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LM `-lm 4` all | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LMM `-lmm 1` Wald | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LMM `-lmm 2` LRT | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LMM `-lmm 3` Score | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LMM `-lmm 4` all | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| 单性状 LMM per-SNP alternative refit | 支持 | 不适用 | 兼容支持 `--lmm-gemma-refit` | smoke 覆盖 |
| GEMMA `-lmm 5` null/VC route | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| GEMMA `-lmm 9` LRT alias | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| VC `-vc 1` HE | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| VC `-vc 2` REML | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| VC `-vc 3` AI/IDUL REML | 支持 | 不适用 | 完全支持 | smoke 覆盖 |
| SNP correlation `-calccor` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Summary reference `-gs` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Summary q vector `-gq` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Summary VC/CI `-study/-ref -vc/-ci` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Summary category `-cat/-mcat/-wcat/-wsnp` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Direct summary `--beta --cor --vc` | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| Prediction `-epm` 风格 | 支持 | 不适用 | 兼容支持 | smoke 覆盖 |
| GEMMA mvLMM mode 1 Wald | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA mvLMM mode 2 LRT | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA mvLMM mode 3 Score | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA mvLMM mode 4 all | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA full alternative mvLMM | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA exact-style mvLMM | 支持 | 不适用 | 完全支持 | reference 硬阈值覆盖 |
| GEMMA CRT/Edgeworth/PCRT | 支持 | 不适用 | 兼容支持 | reference 输出结构覆盖 |
| GEMMA multivariate GxE | 支持相关路线 | 不适用 | 兼容支持 | reference 输出结构覆盖 |
| EMmvGWAS `vc_estimation` | 不适用 | 支持 | 完全支持 `--emmvgwas-vc-only` | smoke 覆盖 |
| EMmvGWAS `full_gwas` | 不适用 | 支持 | 完全支持 `--emmvgwas` | smoke/reference 覆盖 |
| EMmvGWAS `semi_gwas` with G/R | 不适用 | 支持 | 完全支持 `--emmvgwas-G/R` | reference 硬阈值覆盖 |
| EMmvGWAS G/R/beta/association CSV | 不适用 | 支持 | 完全支持 `--emmvgwas-output` | smoke 覆盖 |
| BSLMM mode 1 | MCMC BSLMM | 不适用 | 近似支持 | smoke 覆盖接口 |
| BSLMM mode 2 | Ridge/GBLUP | 不适用 | 完全支持 | GEMMA reference 覆盖 |
| BSLMM mode 3 | Probit BSLMM | 不适用 | 近似支持 | smoke 覆盖接口 |
| BSLMM/DAP mode 4 | DAP BF/HYP | 不适用 | 近似支持 | smoke 覆盖接口 |
| BSLMM/DAP mode 5 | DAP fine mapping | 不适用 | 近似支持 | smoke 覆盖接口 |
| BSLMM/DAP prior 参数 | 支持 | 不适用 | 兼容支持 | smoke 覆盖，含 `log10(pi)` |
| LOCO | 有相关能力 | 不主打 | 完全支持 | 核心回归覆盖 |
| Hybrid EMmvGWAS + GEMMA refinement | 不适用 | 不适用 | CosmicScanner 扩展 | multivariate validation 覆盖 |
| Repeatability GWAS | 不适用 | 不适用 | CosmicScanner 扩展 | 独立验证覆盖 |
| RRM/dynamic GWAS | 不适用 | 不适用 | CosmicScanner 扩展 | 独立验证覆盖 |
| Maternal / ssGWAS / threshold / count / ordinal / survival / GxE reaction norm | 不适用 | 不适用 | CosmicScanner 扩展 | 独立 smoke/calibration 覆盖 |

## 当前硬性验收命令

```bash
cmake --build build --target check_gemma_compat -j 8
cmake --build build --target check_gemma_reference -j 8
cmake --build build --target check_gemma_mv_reference -j 8
cmake --build build --target check_gemma_all -j 8
```

`check_gemma_all` 是聚合入口，覆盖本地兼容 smoke、外部 GEMMA BSLMM2 reference、多性状 GEMMA/EMmvGWAS reference validation。外部 GEMMA 缺失时，reference 测试按 CTest `SKIP_RETURN_CODE=2` 跳过；EMmvGWAS R 包或 R socket 不可用时，多性状脚本会保留 GEMMA reference 验证，并把 EMmvGWAS R reference 标记为 blocked。

## 已记录的 reference 指标

最近一次多性状 GEMMA reference validation 指标：

| 指标 | 数值 | 验收阈值 |
|---|---:|---:|
| EMmvGWAS-compatible EMMV vs GEMMA `log10P` corr | 0.989535528325 | >= 0.98 |
| GEMMA-full vs GEMMA `log10P` corr | 0.998486743763 | >= 0.995 |
| GEMMA-full beta1 corr | 0.998216150373 | >= 0.995 |
| GEMMA-full beta2 corr | 0.999581254704 | >= 0.995 |
| GEMMA-exact vs GEMMA `log10P` corr | 0.997284819487 | >= 0.995 |
| GEMMA-exact beta1 corr | 0.999735913436 | >= 0.995 |
| GEMMA-exact beta2 corr | 0.999609615151 | >= 0.995 |
| EMmvGWAS full vs semi fixed G/R corr | 1.0 | >= 0.999999 |

最近一次 BSLMM2 GEMMA reference validation 指标：

| 指标 | 验收阈值 |
|---|---:|
| alpha corr | >= 0.85 |
| absolute alpha corr | >= 0.70 |
| breeding value corr | >= 0.95 |

## BSLMM/DAP 边界说明

CosmicScanner 的 `--bslmm 2` 是 deterministic ridge/GBLUP 路线，已与 GEMMA `-bslmm 2` 做 reference validation。`--bslmm 1`、`--bslmm 3`、`--bslmm 4` 和 `--bslmm 5` 目前是 GEMMA-style deterministic approximation：

- `--bslmm 1`：sparse+ridge approximation，不运行 GEMMA 原始 MCMC sampler。
- `--bslmm 3`：binary/probit sparse+ridge approximation，不运行 GEMMA 原始 probit MCMC sampler。
- `--bslmm 4`：Wakefield-style ABF grid，输出 GEMMA-style BF/HYP 文件，不等价于 GEMMA DAP exact engine。
- `--bslmm 5`：基于 BF/HYP 的 posterior-weight summary，输出 gamma/coef/hyp，不等价于 GEMMA full DAP EM annotation model。

因此，对外表述应使用“GEMMA-style BSLMM/DAP interface and deterministic approximation”，不能把 mode 1/3/4/5 描述成 GEMMA MCMC/DAP 内核逐行复刻。

## EMmvGWAS R reference 运行条件

外部 EMmvGWAS R reference 需要：

1. 可执行 `Rscript`；
2. R 包 `EMmvGWAS` 已安装；
3. 当前运行环境允许 R `parallel::makeCluster()` 开启本地 PSOCK socket。

当前受限 CTest 环境不能开启 socket，因此脚本会记录：

```text
emmvgwas_preflight    blocked_no_Rscript_EMmvGWAS_package_or_R_socket
emmvgwas_status       blocked_no_Rscript_or_output
```

这不影响 CosmicScanner 自身 EMmvGWAS-compatible full/semi/G/R 路线，也不影响 GEMMA reference validation。在允许 R socket 的普通 shell 环境中，直接运行：

```bash
tests/run_gemma_multivariate_validation.sh
```

如果 EMmvGWAS R 包可用，脚本会自动执行 `gwas_all` reference comparison，并要求 CS vs EMmvGWAS R 的 `log10P` 与 beta 相关性均不低于 `0.999999`。
