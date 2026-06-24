# EB GWAS 发布收口清单

本文档定义 CosmicScanner EB GWAS experimental module 的阶段性结束标准。当前目标不是把 EB 作为正式主检验替代模型，而是作为 validated experimental module 随主线 GWAS 一起发布。

## 发布定位

EB GWAS 当前定位：

1. 单性状 EB-LMM：Wang 2015 风格 random-SNP effect，逐 SNP marker variance ratio。
2. 重复力 EB-GWAS：fast aggregated subject-level EB-LMM。
3. 多性状 EB-GWAS：summary-level EB shrinkage，`Phi_k = lambda_k Ve`。

EB GWAS 不替代标准 fixed-effect LMM、mvLMM、重复力 GWAS。正式分析必须同时报告标准模型结果、普通 Bonferroni/FDR 和 EB 结果。

## 明确不做

当前发布周期不做以下高风险扩展：

1. exact record-level repeatability EB。
2. full REML rank-q per-SNP multivariate EB refit。
3. EB-mBonf 作为唯一正式显著性阈值。

这些内容只有在后续论文/审稿明确需要，并完成大规模 Type I error 验证后再进入新周期。

## Release Gate

快速 gate：

```bash
bash CosmicScanner/tests/run_eb_gwas_suite.sh --quick
```

完整 EB gate：

```bash
bash CosmicScanner/tests/run_eb_gwas_suite.sh --full
```

核心回归集成 gate：

```bash
bash CosmicScanner/tests/run_gwas_core_regression.sh --fast --with-eb-suite
```

发布前至少要求完整 EB gate 和核心回归集成 gate 均通过。

## 必须通过的验证

| 验证 | 脚本 | 目的 |
|---|---|---|
| 单性状 EB smoke | `run_eb_lmm_smoke.sh` | 检查 EB-LMM 输出有限、causal SNP 可恢复 |
| 重复力 EB smoke | `run_repeatability_eb_smoke.sh` | 检查 fast aggregated repeatability EB 路径 |
| 多性状 EB smoke | `run_mv_eb_smoke.sh` | 检查 summary-level mvEB 输出 |
| EB validation matrix | `run_eb_gwas_validation.sh` | null/power 覆盖单性状、重复力、多性状 |
| standard-vs-EB report | `run_eb_comparison_report.sh` | 强制生成标准 LMM 与 EB 对照 |
| EB LOCO validation | `run_eb_loco_validation.sh` | 检查 LOCO 下 EB summary 和输出完整性 |
| EB calibration benchmark | `run_eb_calibration_benchmark.sh` | 快速 null/power calibration |
| EB calibration grid | `run_eb_calibration_grid.sh` | 小型 N/M/effect/noise/pop-structure 网格 |
| PGEN/BGEN consistency | `run_eb_format_consistency.sh` | 检查 PLINK/PGEN/BGEN EB 一致性 |
| EB performance benchmark | `run_eb_performance_benchmark.sh` | 检查 EB-LMM runtime/memory |

## 发布产物

必须生成：

```text
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.status.tsv
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.dashboard.tsv
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.dashboard.md
```

Dashboard 需要至少包含：

1. suite 通过步骤数。
2. validation matrix 的 null/power 指标。
3. standard-vs-EB top SNP 和 top10 overlap。
4. calibration/grid summary。
5. PGEN/BGEN 最大差异。
6. EB performance 最优 rows/sec。

## 结束标准

当以下条件同时满足，即可结束当前 EB GWAS 开发周期：

1. `run_eb_gwas_suite.sh --full` 通过。
2. `run_gwas_core_regression.sh --fast --with-eb-suite` 通过。
3. CLI help 明确标注 EB 为 experimental，并提示不要单独报告 EB-mBonf。
4. `docs/EB_GWAS_Experimental_Module_CN.md` 和本 checklist 均已更新。
5. 不再新增 EB 模型复杂度，只保留 bug fix、验证增强和文档维护。

## 后续论文级模拟

论文前可放大 calibration grid，例如：

```bash
REPS=50 GRID_SPEC="n200_m1000_e0.40_noise0.90_pop0.04 n500_m5000_e0.30_noise1.00_pop0.08" \
  bash CosmicScanner/tests/run_eb_calibration_grid.sh
```

论文级模拟属于后续研究验证，不阻塞当前 experimental module 发布。
