# GWAS 板块收口验收记录

日期：2026-07-08

状态：发布前阻断项已清理。当前代码可以构建，核心 smoke/reference/unit 测试已纳入 CTest 门禁，CosmicGWAS 可以按稳定板块发布；外部 reference 和近似/interface 路线仍按各自方法边界说明。

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
ctest --test-dir build --output-on-failure
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

2026-07-08 在 `/disk228/huangjie/PanOmics/CosmicScanner` 完成以下发布前代码修复：

| 文件 | 处理结果 | 发布影响 |
|---|---|---|
| `src/models/repeatability_model.cpp` | 移除 `curl`/`std::system` 调试上报，保留本地 `std::cerr` 诊断 | 消除发布包中的外部调试上报和 shell 调用风险 |
| `src/models/repeatability_scan_exact.cpp` | 移除 `curl`/`std::system` 调试上报，失败时输出本地诊断 | 消除发布包中的外部调试上报和 shell 调用风险 |
| `src/gwas_app.cpp`、`CMakeLists.txt` | 通过 `COSMICSCANNER_VERSION` 输出 CMake 项目版本 | 程序日志版本与构建版本统一 |
| `CMakeLists.txt` | 将 `test_sparse_reml`、`test_multi_vc`、`test_spa_accuracy`、`test_regenie`、`test_input_adaptation_units` 加入 `add_test` | 默认 CTest 门禁从 5 项扩展到 10 项 |
| `src/lmm_engine.cpp` | PGEN biobank scan 增加 keep-index 子集到全样本向量的转换包装 | 修复过滤样本后 PGEN 矩阵乘法的维度/样本错配风险 |

2026-07-07 至 2026-07-08 在 `/disk228/huangjie/PanOmics/CosmicScanner` 完成以下发布前代码阅读和最小验证：

| 命令 | 结果 | 说明 |
|---|---|---|
| `git status --short` | PASS | 工作区进入评估前干净 |
| `cmake --build build -j2` | PASS | 主程序和测试可执行文件构建通过 |
| `ctest --test-dir build --output-on-failure` | PASS | 2026-07-07 CTest 5/5 通过：compat smoke、BSLMM reference、多性状 reference、EB-GWAS quick suite、ssGWAS A22 smoke |
| `./build/test_sparse_reml` | PASS | Sparse REML 单元测试通过 |
| `./build/test_multi_vc` | PASS | 多方差组分 null model 单元测试通过 |
| `./build/test_spa_accuracy` | PASS | SPA 数值测试通过 |
| `./build/test_regenie` | PASS | REGENIE Step 1 单元测试通过 |
| `./build/test_input_adaptation_units` | PASS | 输入适配单元测试通过 |

CTest 门禁说明：

- 2026-07-08 起，默认 CTest 包含 10 项：5 个 smoke/reference 测试 + 5 个核心单元测试。
- `tests` 顶层仍保留大量扩展 benchmark/reference 脚本；这些脚本用于发布增强验收或条件 reference，不全部纳入默认快速门禁。

条件 reference 说明：

- 本机 GEMMA 可用，BSLMM2 reference 和多性状 GEMMA reference 已通过 CTest。
- EMmvGWAS R reference 预检记录为 `blocked_no_Rscript_EMmvGWAS_package_or_R_socket`；CosmicScanner 原生 EMmvGWAS-compatible full/semi/G/R 路线已在 smoke/reference 框架中通过。
- PSR R/GLM adapter 在受限环境中会打印 `yp_bind_client_create_v3` RPC 警告，但候选精修输出正常生成，聚合验收命令最终 PASS。

## 4. 发布表述边界

- 标准 LM/LMM、GRM、GEMMA/EMmvGWAS 兼容 smoke 覆盖路径、EB-GWAS quick、ssGWAS smoke、SPA/REGENIE/Input adaptation 单元测试覆盖路径可以进入稳定发布说明。
- `--bslmm 2` 是 deterministic ridge/GBLUP 路线，可作为 GEMMA-style BSLMM2 reference 覆盖功能表述。
- `--bslmm 1/3/4/5` 是 GEMMA-style deterministic approximation/interface，不是 GEMMA 原始 MCMC BSLMM、probit MCMC 或 DAP exact engine。
- EB-GWAS 默认路线已通过稳定验收，可作为标准 LMM/mvLMM/重复力 GWAS 的辅助信号增强路线；正式报告仍必须同时保留标准主结果。当前稳定范围包括单性状 EB-LMM、重复力 fast aggregated subject-level EB-LMM、多性状 summary-level EB shrinkage。exact record-level repeatability EB、full rank-q per-SNP multivariate EB refit 和 RRM EB-GWAS 仍属于后续增强。
- PSR-GWAS 是 PSR/PQL 近似路线。二元、计数、二项和序数模型应报告 convergence、warning、candidate refinement 或外部 GLMM/ordinal baseline 复核。
- `--repeatability-mode fast` 是 subject-level 均值聚合近似，不替代 `--repeatability-mode exact` 正式主分析。
- `--repeatability-mode stcg` 是大规模矩阵自由原型，默认不作为课题组生产路线。
- RRM dynamic GWAS 和 repeatability GWAS 是不同模型；RRM 用于动态轨迹和时间窗口，repeatability 用于重复记录总体效应。
- `--rrm-loco-global-null` 是快速近似 LOCO，不替代默认逐染色体重估 null model 的 LOCO。
- EMmvGWAS-compatible 路线可直接运行 CosmicScanner 原生实现；与 EMmvGWAS R 包的逐项 reference 需要满足 R 包和 socket 条件。
- 输入适配层仍是 MVP。帮助信息中明确 `xlsx not in MVP`，发布说明应写成 csv/tsv/txt prepared-input 工具，不应写成完整 ETL/表格适配平台。
- `--ordinal-model` 默认 cumulative-binary 路线仍带 experimental 标识，发布说明应要求报告 convergence、warning、fallback 和外部 baseline 复核。
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
