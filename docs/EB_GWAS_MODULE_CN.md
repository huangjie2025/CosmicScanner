# CosmicScanner EB-GWAS 模块

本文档记录 CosmicScanner 当前经验贝叶斯 GWAS 模块的实现范围、推荐用法和验证策略。EB-GWAS 默认路线已经完成稳定落地，可作为标准 fixed-effect LMM/mvLMM/重复力 GWAS 的辅助信号增强路线；正式分析仍应同时报告标准主模型结果和普通 Bonferroni/FDR。

## 当前落地结论

已稳定落地的默认路线：

- 单性状 EB-LMM：`--model lmm-eb`。
- 重复力 EB-GWAS：`--repeatability-model --model lmm-eb`，采用 fast aggregated subject-level EB-LMM。
- 多性状 EB-GWAS：`--multivariate --mv-method eb`，采用 summary-level EB shrinkage。

本模块已经通过 smoke、validation matrix、standard-vs-EB report、LOCO validation、calibration benchmark、calibration grid、PLINK/PGEN/BGEN format consistency 和 performance benchmark。总入口为：

```bash
bash tests/run_eb_gwas_suite.sh --full
```

CMake/CTest 入口：

```bash
cmake --build build --target check_eb_gwas
cmake --build build --target check_eb_gwas_full
ctest --test-dir build -R eb_gwas_suite_quick --output-on-failure
```

## 已支持模型

| 模块 | 命令 | 当前实现 | 定位 |
|---|---|---|---|
| 单性状 EB-LMM | `--model lmm-eb` 或 `--test eb-wald` | Wang 2015 风格 random-SNP effect，逐 SNP 估计 marker variance ratio | 稳定辅助路线 |
| 重复力 EB-GWAS | `--repeatability-model --model lmm-eb` | 加权 subject-level 聚合均值 EB-LMM | 稳定辅助路线 |
| 多性状 EB-GWAS | `--multivariate --mv-method eb` | summary-level EB shrinkage，`Phi_k = lambda_k Ve` | 稳定辅助路线 |

## 方法边界

当前重复力 EB 不是 exact record-level EB。它先将重复记录聚合为加权 subject-level 均值，再复用单性状 EB-LMM 引擎。该路径稳定、快，已经作为默认落地路线；exact record-level EB 需要额外处理记录层协方差、LOCO、Woodbury 和 Type I error，暂不作为默认功能。

当前多性状 EB 不是完整 REML rank-q per-SNP variance-component refit。它先用 mvLMM fixed-null 得到 SNP 效应向量和协方差，再用 `c_k ~ N(0, lambda_k Ve)` 做 summary-level shrinkage。这样每 SNP 只估计一个 `lambda_k`，避免完整 `q x q` SNP 协方差在小样本下不稳定。该 summary-level EB 已作为稳定辅助路线落地；full rank-q per-SNP EB refit 仍属于后续增强。

## 推荐用法

单性状：

```bash
CosmicScanner/build/cosmicscanner \
  --bfile data_prefix \
  --pheno pheno.txt \
  --mpheno 1 \
  --gwas \
  --model lmm-eb \
  --adjust eb-mbonf \
  --out result/single_eb
```

重复力：

```bash
CosmicScanner/build/cosmicscanner \
  --bfile data_prefix \
  --pheno long_pheno.txt \
  --mpheno 1 \
  --repeatability-model \
  --gwas \
  --model lmm-eb \
  --adjust eb-mbonf \
  --out result/repeatability_eb
```

多性状：

```bash
CosmicScanner/build/cosmicscanner \
  --bfile data_prefix \
  --pheno pheno.txt \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method eb \
  --gwas \
  --out result/mv_eb
```

## 输出文件

单性状 EB 输出：

```text
<out>.assoc.txt
<out>.eb.summary.txt
```

重复力 EB 输出：

```text
<out>.assoc.txt
<out>.repeatability.eb.summary.txt
```

多性状 EB 输出：

```text
<out>.assoc.mv.txt
<out>.eb.mv.summary.txt
```

核心字段包括：

| 字段 | 含义 |
|---|---|
| `EB_BETA` / `EB_BETA_i` | shrink 后的 SNP 效应 |
| `EB_SE` / `EB_SE_i` | posterior SE |
| `EB_WALD` | EB-Wald statistic |
| `MARKER_LAMBDA` | marker-specific prior variance ratio |
| `PRIOR_VAR` | 单性状/重复力 EB 的 prior variance |
| `POST_VAR` | 单性状/重复力 EB 的 posterior variance |
| `EB_CONFIDENCE` | effective-test contribution |
| `effective_tests` | `sum(EB_CONFIDENCE)` |
| `effective_tests_per_marker` | `effective_tests / tested_markers`，用于快速审查 shrinkage 强度 |
| `eb_mbonf_0.05` | `0.05 / effective_tests` |

说明：单性状和重复力 EB 的 `EB_CONFIDENCE` 范围是 `[0, 1]`；多性状 EB 的 `EB_CONFIDENCE` 范围是 `[0, trait_count]`，表示多性状联合 shrinkage 的 effective contribution。

## 验证脚本

EB GWAS 路线板块总入口：

```bash
bash CosmicScanner/tests/run_eb_gwas_suite.sh --quick
bash CosmicScanner/tests/run_eb_gwas_suite.sh --full
```

输出 dashboard：

```text
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.status.tsv
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.dashboard.tsv
CosmicScanner/tests/results/eb_gwas_suite/eb_gwas_suite.dashboard.md
```

该入口串联单性状 EB、重复力 EB、多性状 EB、validation matrix、standard-vs-EB report、LOCO validation、calibration benchmark、PGEN/BGEN format consistency 和 performance benchmark，是当前 EB GWAS 模块的 board-level 验证入口。

发布收口清单见：

```text
CosmicScanner/docs/EB_GWAS_MODULE_CN.md (this document)
```

快速 smoke：

```bash
bash CosmicScanner/tests/run_eb_lmm_smoke.sh
bash CosmicScanner/tests/run_repeatability_eb_smoke.sh
bash CosmicScanner/tests/run_mv_eb_smoke.sh
```

EB 验证矩阵：

```bash
bash CosmicScanner/tests/run_eb_gwas_validation.sh
```

标准 LMM 与 EB-LMM 对照报告：

```bash
bash CosmicScanner/tests/run_eb_comparison_report.sh
```

输出：

```text
CosmicScanner/tests/results/eb_comparison_report/standard_vs_eb.tsv
CosmicScanner/tests/results/eb_comparison_report/standard_vs_eb.summary.tsv
CosmicScanner/tests/results/eb_comparison_report/standard_vs_eb.summary.md
```

该报告把标准 fixed-effect LMM 与 EB-LMM 的 `P/BETA/SE`、`MARKER_LAMBDA`、`EB_CONFIDENCE` 放在同一个表中，便于正式分析时同时审查标准模型和 EB 模型。

EB LOCO 专项验证：

```bash
bash CosmicScanner/tests/run_eb_loco_validation.sh
```

该脚本构造 5 条染色体模拟数据，默认启用 LOCO，检查所有染色体输出完整、`EB_CONFIDENCE` 范围合理、`effective_tests` 与逐 SNP 求和一致，并确认 power 场景 causal SNP 排名。

EB calibration benchmark：

```bash
bash CosmicScanner/tests/run_eb_calibration_benchmark.sh
```

可通过 `REPS=<n>` 控制重复次数。输出：

```text
CosmicScanner/tests/results/eb_calibration_benchmark/eb_calibration.summary.tsv
CosmicScanner/tests/results/eb_calibration_benchmark/eb_calibration.aggregate.tsv
CosmicScanner/tests/results/eb_calibration_benchmark/eb_calibration.summary.md
```

该 benchmark 同时运行标准 LMM 与 EB-LMM 的 null/power 重复模拟，汇总普通 Bonferroni、EB-mBonf、`effective_tests` 和 causal top10 命中情况。它是默认工程校准入口；正式论文如需对某个新物种/群体结构/LD 场景做强声明，仍建议放大 Type I error 模拟。

EB calibration grid：

```bash
bash CosmicScanner/tests/run_eb_calibration_grid.sh
```

默认运行一个小型 N/M/effect/noise/pop-structure 网格，并汇总：

```text
CosmicScanner/tests/results/eb_calibration_grid/eb_calibration_grid.summary.tsv
CosmicScanner/tests/results/eb_calibration_grid/eb_calibration_grid.summary.md
```

可通过环境变量放大模拟：

```bash
REPS=50 GRID_SPEC="n200_m1000_e0.40_noise0.90_pop0.04 n500_m5000_e0.30_noise1.00_pop0.08" \
  bash CosmicScanner/tests/run_eb_calibration_grid.sh
```

`GRID_SPEC` 每个单元格式为 `n<N>_m<M>_e<EFFECT>_noise<SD>_pop<STRUCTURE>`。

PGEN/BGEN EB 一致性验证：

```bash
bash CosmicScanner/tests/run_eb_format_consistency.sh
```

该脚本使用 `plink2` 将同一 PLINK hard-call 数据导出为 PGEN 和 BGEN，分别运行 EB-LMM，并逐列比较 `EB_BETA`、`EB_SE`、`EB_WALD`、`P`、`MARKER_LAMBDA`、`POST_VAR`、`EB_CONFIDENCE` 等核心字段。

EB-LMM 性能 benchmark：

```bash
bash CosmicScanner/tests/run_eb_performance_benchmark.sh
```

可通过环境变量配置规模和线程：

```bash
N=500 M=5000 THREADS_LIST="1 4 8" bash CosmicScanner/tests/run_eb_performance_benchmark.sh
```

输出：

```text
CosmicScanner/tests/results/eb_performance_benchmark/eb_performance.summary.tsv
CosmicScanner/tests/results/eb_performance_benchmark/eb_performance.summary.md
```

核心回归中追加 EB 验证矩阵：

```bash
bash CosmicScanner/tests/run_gwas_core_regression.sh --fast --with-eb-validation
```

追加所有 EB 扩展验证：

```bash
bash CosmicScanner/tests/run_gwas_core_regression.sh --fast \
  --with-eb-validation \
  --with-eb-report \
  --with-eb-loco \
  --with-eb-calibration \
  --with-eb-formats \
  --with-eb-perf
```

或者直接运行 board-level EB suite：

```bash
bash CosmicScanner/tests/run_gwas_core_regression.sh --fast --with-eb-suite
```

验证矩阵覆盖：

| 场景 | 目标 |
|---|---|
| 单性状 null | 检查 EB 输出有限、effective tests 合理、null 下 EB-mBonf 命中不过度膨胀 |
| 单性状 power | 检查 causal SNP rank |
| 重复力 null | 检查加权聚合 EB 的 null 行为 |
| 重复力 power | 检查重复记录场景 causal SNP rank |
| 多性状 null | 检查 summary-level EB shrinkage 的 null 行为 |
| 多性状 power | 检查联合效应 causal SNP rank |
| EB LOCO | 检查 LOCO 下 EB 输出完整性、summary 一致性和 causal SNP rank |
| EB calibration | 快速重复模拟，审查 null/power、Bonferroni/EB-mBonf 和 effective tests |
| standard-vs-EB report | 强制生成标准 LMM 与 EB-LMM 对照表，避免只报告 EB 结果 |

## 性能说明

单性状 PLINK EB-LMM 扫描已对 SNP 循环做 OpenMP 并行化；线程数遵循 `--threads` 或运行环境中的 OpenMP 设置。BGEN/PGEN EB 路径仍保持顺序 streaming 读取，避免破坏 reader 状态。

## 报告建议

正式分析不要只报告 EB-mBonf。建议同时输出并报告：

1. 标准 fixed-effect LMM/mvLMM/重复力 GWAS 结果。
2. EB-GWAS 结果。
3. 普通 Bonferroni/FDR。
4. EB-mBonf 作为辅助阈值。
5. 如研究对象的群体结构、LD、重复记录不平衡或性状相关性明显不同于内置验证场景，追加目标数据的 Type I error 模拟。

## 后续可选增强

以下路线未作为当前稳定默认功能承诺，只有在论文/审稿确实需要时再考虑：

1. exact record-level repeatability EB。
2. full REML rank-q per-SNP multivariate EB refit。
3. 随机回归 EB-GWAS。

这些增强不建议作为当前主线开发目标。
