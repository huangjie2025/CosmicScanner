# GWAS 板块收口验收记录

日期：2026-07-07

## 1. 收口范围

本记录用于 GWAS model suite 集成后的发布前验收。覆盖范围包括：

- 单性状 LM/LMM 与 GEMMA-style `-lm/-lmm` 路线。
- 多性状 GWAS、GEMMA full/exact/GxE、EMmvGWAS-compatible full/semi/vc 路线。
- 重复力 GWAS、RRM/dynamic GWAS。
- threshold/count/binomial/ordinal PSR-GWAS。
- maternal、ssGWAS、survival、GxE reaction norm。
- rare variant SKAT/Burden、SPA、Biobank/REGENIE-style、EB-GWAS。
- GEMMA-style summary-stat、VC、eigen、calccor、prediction、BSLMM/DAP interface。

## 2. 本地验收命令

最小必跑：

```bash
cmake --build build --target cosmicscanner -j4
ctest --test-dir build -R gemma_emmvgwas_compat_smoke --output-on-failure
bash tests/run_gwas_core_regression.sh --compat-smoke-only
```

扩展本地验收：

```bash
bash tests/run_gwas_core_regression.sh --fast --with-psr
bash tests/run_repeatability_release_validation.sh
bash tests/run_rrm_dynamic_simulation.sh --loco-smoke
bash tests/run_eb_gwas_suite.sh --quick
```

条件 reference 验收：

```bash
cmake --build build --target check_gemma_reference -j8
cmake --build build --target check_gemma_mv_reference -j8
cmake --build build --target check_gemma_all -j8
```

这些 reference 命令依赖外部 GEMMA、R、EMmvGWAS R 包以及允许 R PSOCK socket 的运行环境。依赖缺失时，脚本应跳过或标记 blocked，不应解读为 CosmicScanner 原生路径失败。

## 3. 本轮实际验收结果

2026-07-07 在 `/disk228/huangjie/PanOmics/CosmicScanner` 完成以下验收：

| 命令 | 结果 | 说明 |
|---|---|---|
| `cmake --build build --target cosmicscanner -j4` | PASS | 主程序构建通过 |
| `ctest --test-dir build -R gemma_emmvgwas_compat_smoke --output-on-failure` | PASS | 本地 GEMMA/EMmvGWAS 兼容 smoke 通过 |
| `bash tests/run_gwas_core_regression.sh --compat-smoke-only` | PASS | 脚本层兼容 smoke 通过 |
| `bash tests/run_test_multivariate.sh` | PASS | 多性状 EMMV/GEMMA/GEMMA-full/hybrid smoke 通过 |
| `bash tests/run_test_repeatability.sh` | PASS | repeatability exact/cache/LOCO/fast smoke 通过 |
| `bash tests/run_repeatability_release_validation.sh` | PASS | repeatability smoke、edge cases、exact-vs-fast simulation smoke 通过 |
| `bash tests/run_gwas_core_regression.sh --fast --with-psr` | PASS | EB smoke、多性状 smoke、GEMMA/EMmvGWAS validation、PSR smoke、ordinal dense-zkz validation 通过 |
| `bash tests/run_eb_gwas_suite.sh --full` | PASS | EB 单性状、重复力、多性状、validation matrix、LOCO、calibration、calibration grid、PGEN/BGEN consistency、performance 全部通过 |
| `cmake --build build --target check_eb_gwas` | PASS | EB-GWAS CTest quick suite 通过 |
| `cmake --build build --target check_eb_gwas_full` | PASS | EB-GWAS CMake full suite 通过 |
| `ctest --test-dir build --output-on-failure` | PASS | CTest 4/4 通过：compat smoke、BSLMM2 reference、多性状 reference、EB-GWAS quick suite |

本轮同时修复了两个测试脚本的发布前可复现性问题：

- `tests/run_test_multivariate.sh`：缺少 `examples/mouse_clean` 时自动使用 `tests/data/mv_eb_smoke` 内置 PLINK 数据。
- `tests/run_test_repeatability.sh`：缺少 `examples/mouse_clean` 时自动使用 `tests/data/repeatability_eb_smoke` 内置 PLINK/long phenotype 数据。
- 两个脚本均动态读取 `.fam` 样本数，并兼容当前测试数据的 5 列 `.bim` 格式。
- repeatability smoke 的输出检查改为按输出前缀追加 `.assoc.txt`，避免 `Path.with_suffix()` 错误替换前缀后缀。

条件 reference 说明：

- 本机 GEMMA 可用，BSLMM2 reference 和多性状 GEMMA reference 已通过 CTest。
- EMmvGWAS R reference 预检记录为 `blocked_no_Rscript_EMmvGWAS_package_or_R_socket`；CosmicScanner 原生 EMmvGWAS-compatible full/semi/G/R 路线已在 smoke/reference 框架中通过。
- PSR R/GLM adapter 在受限环境中会打印 `yp_bind_client_create_v3` RPC 警告，但候选精修输出正常生成，聚合验收命令最终 PASS。

## 4. 发布表述边界

- `--bslmm 2` 是 deterministic ridge/GBLUP 路线，可作为 GEMMA-style BSLMM2 reference 覆盖功能表述。
- `--bslmm 1/3/4/5` 是 GEMMA-style deterministic approximation/interface，不是 GEMMA 原始 MCMC BSLMM、probit MCMC 或 DAP exact engine。
- EB-GWAS 默认路线已通过稳定验收，可作为标准 LMM/mvLMM/重复力 GWAS 的辅助信号增强路线；正式报告仍必须同时保留标准主结果。当前稳定范围包括单性状 EB-LMM、重复力 fast aggregated subject-level EB-LMM、多性状 summary-level EB shrinkage。exact record-level repeatability EB、full rank-q per-SNP multivariate EB refit 和 RRM EB-GWAS 仍属于后续增强。
- PSR-GWAS 是 PSR/PQL 近似路线。二元、计数、二项和序数模型应报告 convergence、warning、candidate refinement 或外部 GLMM/ordinal baseline 复核。
- `--repeatability-mode fast` 是 subject-level 均值聚合近似，不替代 `--repeatability-mode exact` 正式主分析。
- `--repeatability-mode stcg` 是大规模矩阵自由原型，默认不作为课题组生产路线。
- RRM dynamic GWAS 和 repeatability GWAS 是不同模型；RRM 用于动态轨迹和时间窗口，repeatability 用于重复记录总体效应。
- `--rrm-loco-global-null` 是快速近似 LOCO，不替代默认逐染色体重估 null model 的 LOCO。
- EMmvGWAS-compatible 路线可直接运行 CosmicScanner 原生实现；与 EMmvGWAS R 包的逐项 reference 需要满足 R 包和 socket 条件。

## 5. 文档收口

本轮统一后的主文档入口：

- `README.md`
- `docs/Tutorial_CN.md`
- `docs/Tutorial_EN.md`
- `docs/GEMMA_EMMVGWAS_COSMICSCANNER_SUPPORT_MATRIX_CN.md`
- `docs/PSR_GWAS_ROADMAP.md`
- `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md`
- `docs/REPEATABILITY_GWAS_RELEASE_CN.md`
- `docs/RRM_GWAS_VALIDATION_CN.md`
- `docs/EB_GWAS_MODULE_CN.md`

旧阶段性 PSR 总结文档已合并到 `docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md` 和路线图；旧 repeatability/RRM 文档名已统一为全大写命名。
