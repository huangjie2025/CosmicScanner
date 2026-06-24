# CosmicScanner GWAS 1.0 开发、测试与验证方案

## 目标定位

1.0 版本优先完成 GWAS 板块，定位为“主要 GWAS 模型可运行、可测试、可汇报、可打包”的功能全面版。当前版本以已实现并可验证的模型为发布核心，实验性模型仅作为后续扩展储备，不作为 1.0 主打承诺。

## 1.0 发布核心模型

| 模型 | 命令入口 | 1.0 状态 | 发布要求 |
|---|---|---|---|
| 单性状 LMM GWAS | 默认 `--gwas` | 核心发布功能 | 与 GEMMA 对比，检查 beta、SE、P 值相关性 |
| 多性状 GWAS | `--multivariate --mv-method emmv` | 核心发布功能 | 接入 EMmvGWAS 半精准路线，小样本运行通过，输出联合检验结果；`--mv-method gemma` 和 `--mv-method gemma-full` 作为精准 refit 对照路线 |
| 重复力 GWAS | `--repeatability-model` | 核心发布功能 | exact、LOCO、cache、fast 路线验证通过 |
| 随机回归 GWAS | `--long-pheno` + `--rrm-order` | 核心发布功能 | 长格式数据运行通过，输出联合检验与曲线结果 |

## 1.0 暂不主打模型

以下模型源码中已有入口或原型，但 1.0 不作为正式发布主打功能，避免对老师和用户过度承诺：

| 模型 | 原因 | 后续版本 |
|---|---|---|
| 母体效应 GWAS | 需要独立数据格式、方差结构和外部软件对标 | 1.1/2.0 |
| ssGWAS / single-step | 需要系谱、基因组、H 矩阵流程系统验证 | 1.1/2.0 |
| 阈值/二分类模型 | 需要二元性状校准、SPA/GLMM 对标 | 2.0 |
| 计数性状模型 | 需要 Poisson/NegBin 模型验证 | 2.0 |
| 生存模型 | 需要 Cox/混合生存模型验证 | 2.0 |
| GxE 反应规范模型 | 需要环境梯度与交互效应基准验证 | 2.0 |
| Biobank mode | 可作为技术储备，不作为 1.0 生产主线 | 2.0/3.0 |

## 多性状 GWAS 路线

| 路线 | 命令 | 定位 | 推荐用途 |
|---|---|---|---|
| EMmvGWAS 半精准路线 | `--mv-method emmv` | 默认主线 | 常规多性状 GWAS，全基因组扫描 |
| GEMMA 风格精准 refit 路线 | `--mv-method gemma` | 精准对照/精修 | 小样本验证、top SNP 精修、与 GEMMA 对标 |
| GEMMA-full MVP 路线 | `--mv-method gemma-full --mv-test wald,lrt,score` | 完整检验输出 | 小样本验证、Wald/LRT/Score 对照 |
| Hybrid 优化路线 | `--mv-method hybrid --mv-refine-p 1e-3 --mv-refine-top 100` | 性能优化主线 | 全基因组快扫 + 候选 SNP 精修 |
| 固定 null Wald 路线 | `--mv-method wald` | 调试对照 | 方法开发、排查差异 |

GEMMA 风格路线当前实现为逐 SNP 重新估计多性状方差组分并进行 Wald 检验，属于比 null-reuse 更精细但更慢的路线。GEMMA-full MVP 当前输出 Wald、LRT、Score 三类 p 值，其中 LRT 基于逐 SNP refit 后的似然差，Score 基于 null 方差结构的一步检验。Hybrid 路线复用 EMmvGWAS 快速全基因组扫描，仅对候选 SNP 做 GEMMA-full 精修，是当前最接近“功能全面且速度优先”的路线。完整 GEMMA 的 PX-EM + Newton-Raphson 深度精修和 p 值校正仍作为后续增强目标。

当前 GEMMA/EMmvGWAS 对标脚本会自动生成 96 个体、360 SNP、2 性状模拟数据，并输出 `tests/results/multivariate_gemma_validation/multivariate_validation.summary.txt`。R 环境使用 `/disk228/huangjie/software/miniconda3/envs/R_GS/bin/Rscript`，已安装 `doParallel`、`foreach`、`RSpectra` 和本地 `EMmvGWAS` 包。

## 优先开发清单

1. 统一 1.0 文档口径：README、教程、发布说明只强调单性状、多性状、重复力、随机回归四类已验证模型。
2. 检查 CLI 与文档一致性：重点确认 `--long-pheno`、`--rrm-order`、`--repeatability-model`、`--multivariate` 示例命令可直接运行。
3. 补齐每类模型的小样本示例：每类模型至少保留一个可自动生成或可复用的数据集。
4. 补齐每类模型的验证报告：记录输入规模、运行命令、输出文件、关键结果、耗时和内存。
5. 打包发布前运行核心回归脚本，形成 1.0 验证记录。

## 核心测试矩阵

| 测试目标 | 脚本 | 验证内容 |
|---|---|---|
| 单性状准确性 | `tests/run_gemma_accuracy.sh` | 与 GEMMA 对比 beta、SE、P 值 |
| 单性状速度 | `tests/run_gemma_speed_benchmark.sh` | 与 GEMMA 或 Cosmic-only 速度对比 |
| 多性状功能 | `tests/run_test_multivariate.sh` | 多表型输入、联合结果输出、无 NaN/崩溃 |
| 多性状 GEMMA 对标 | `tests/run_gemma_multivariate_validation.sh` | 模拟双性状数据，比较 CosmicScanner 与 GEMMA 的 beta、p 值相关性和耗时 |
| 重复力功能 | `tests/run_test_repeatability.sh` | exact、cache、LOCO、fast 路线一致性 |
| 重复力边界 | `tests/run_test_repeatability_edge_cases.sh` | 单记录、非平衡记录、缺失等边界情况 |
| 随机回归功能 | `tests/run_test_rrm.sh` | RRM 元数据、联合检验、曲线输出 |
| PGEN 一致性 | `tests/run_pgen_consistency.sh` | PGEN 输入与 PLINK 路线一致性 |

## 发布前验收标准

1. 四类核心模型均有可运行命令和小样本结果。
2. 单性状 GWAS 与 GEMMA 对比达到报告中设定的相关性阈值。
3. 重复力 GWAS exact/cache/LOCO/fast 测试通过。
4. RRM GWAS 能输出 `.rrm_joint.assoc.txt`、`.rrm_joint.null.txt`、`.rrmcoef.txt`、`.rrm.top_curve.txt`。
5. 教程中的 1.0 核心命令与实际 CLI 一致。
6. 发布包包含可执行文件、中文教程、小样本数据、测试报告和已知问题说明。

## 给老师汇报口径

当前开发策略已调整为优先完成 GWAS 板块。1.0 版本将以单性状、多性状、重复力和随机回归 GWAS 为核心功能，通过小样本数据完成自动化测试和外部工具对标验证，形成可打包发布版本。GS 板块保留基础能力，后续版本再系统补齐多性状、重复力和随机回归等育种值估计模型。
