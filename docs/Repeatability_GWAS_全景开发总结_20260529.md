# CosmicScanner 重复力 GWAS 全景开发总结

> 生成日期：2026-05-29
> 覆盖范围：原理、实现流程路线、代码、测试数据集、测试方法、验证方法、验证结果、准确性、时间内存消耗、文件废弃分析

---

## 一、统计模型与原理

### 1.1 模型定义

重复力 GWAS 面向同一个体存在多次重复观测的表型数据。设共有 N 个基因分型个体，第 i 个个体对应 n_i 条表型记录，总记录数为 n_obs = Σn_i。第 i 个个体的第 j 次记录模型为：

```
y_ij = x_ij' β + g_i + pe_i + m_i α + e_ij
```

其中：
- x_ij：记录级固定效应与协变量（如胎次 PARITY）
- β：固定效应系数
- g_i：个体加性遗传随机效应，g ~ N(0, Vg·G)
- pe_i：个体永久环境随机效应，pe ~ N(0, Vpe·I)
- m_i：SNP 剂量（个体内所有重复记录共享同一基因型）
- α：待检验 SNP 效应
- e_ij：记录级残差，e ~ N(0, Ve·I_obs)

### 1.2 协方差结构

记录级总协方差矩阵：

```
Var(y) = Z(G·Vg + I·Vpe)Z' + I·Ve
```

其中 Z 为记录到个体的 incidence matrix。对角块结构：

- 同一个体不同记录间：Cov(y_ij, y_ik) = Vg·G_ii + Vpe
- 不同个体记录间：Cov(y_ij, y_i'k) = Vg·G_ii'

### 1.3 关键参数

| 参数 | 定义 | 生物学含义 |
|---|---|---|
| Vp = Vg + Vpe + Ve | 总表型方差 | — |
| h2 = Vg / Vp | 狭义遗传力 | 加性遗传方差占比 |
| r2 = (Vg + Vpe) / Vp | 重复力 | 个体稳定差异占比 |

### 1.4 与 RRM 的区别

| 特征 | 纯重复力模型 | RRM 随机回归模型 |
|---|---|---|
| 触发参数 | `--repeatability-model` | `--long-pheno` |
| 随机效应 | g_i + pe_i（时不变） | 遗传效应随时间/阶数变化 |
| 适用场景 | 重复力 GWAS | 纵向轨迹建模 |
| 结果目录 | `results_pure_rep_*` | `results_rrm_model/` |

**重要**：RRM 结果不能当作纯重复力结果汇报，二者模型假设不同。

---

## 二、实现流程路线

### 2.1 整体架构

```
main.cpp → GwasApp::run() → RepeatabilityModel
                                ├── loadPhenotype()    长格式表型读入
                                ├── loadCovariates()   记录级协变量对齐
                                ├── buildGRM()         个体间 GRM 构建
                                ├── buildMMEComponents() 组装随机效应
                                ├── runVCE()           方差组分估计
                                └── run() → scanExact() / scanFast() / scanSTCG()
```

### 2.2 三路线扫描策略

| 路线 | 原理 | 复杂度 | 成熟度 | 推荐场景 |
|---|---|---|---|---|
| **Exact (R1)** | 记录层精确 GLS，Woodbury 降维到个体层 O(N³) | O(n_subj³ + n_obs·n_subj) | ⭐⭐⭐ 最成熟 | 正式结果、论文主分析 |
| **Fast (R2)** | 个体均值聚合 + 加权近似 | O(n_subj³) | ⭐⭐ 可用 | 大样本预筛 |
| **STCG (R3)** | 矩阵自由 PCG/score 型 | O(n_iter·nnz) | ⭐ 原型 | Biobank 级长期储备 |

Auto 模式选择逻辑：
- n_records > 80000 → STCG
- n_records > 20000 或 `--fast-lmm` → Fast
- 其他 → Exact

### 2.3 VCE 路线演进（四代）

| 代次 | 方法 | 特点 | 时间 |
|---|---|---|---|
| 第一代 | multi-start EM-REML + 长 EM refinement (300次) | 稳定但 Vg/Vpe 拆分慢 | 2026-04 |
| 第二代 | exact trace：tr(Q_k⁻¹ C_kk) 精确计算 | Vpe 稳定性提升 | 2026-05-26 |
| 第三代 | GLS fix：SNP 扫描公式修正为 full record-level GLS | 非平衡重复记录精度修复 | 2026-05-27 |
| 第四代 | fast subject-level REML refinement（当前推荐） | 速度提升 1.35x–1.59x | 2026-05-28 |

### 2.4 Exact 扫描核心公式

利用 Woodbury 恒等式将 record-level V⁻¹ 降维到 subject-level：

1. 构建 K = Vg·G + Vpe·I（subject-level）
2. H = K + Ve·diag(1/D_i)，其中 D_i = 个体 i 的记录数
3. V⁻¹y = (y_i - u_i)/Ve，其中 u = H⁻¹·(K·y_sub) ./ D
4. SNP score = g_sub'·Py_sub_sum
5. SNP var = g'V⁻¹g - (X'V⁻¹g)'(X'V⁻¹X)⁻¹(X'V⁻¹g)

### 2.5 LOCO 实现

两种路径：
1. **单次全基因组拆分 GRM**（默认，无 cache 时）：一次遍历基因型数据，同时计算 global GRM 和 per-chromosome GRM，然后 G_loco = (S_global·G_global - S_chr·G_chr) / S_loco
2. **cache 路径**：先计算并缓存 global + per-chr GRM 到 `.f32grm.bin` + `.meta.txt`，后续运行通过 mmap 加载

---

## 三、代码文件清单与修改时间

### 3.1 核心源代码

| 文件 | 行数 | 最后修改 | 职责 |
|---|---|---|---|
| [repeatability_model.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/src/models/repeatability_model.cpp) | 1733 | 2026-05-28 22:14 | 重复力 GWAS 主模型：数据读入、GRM、VCE、三路线扫描、LOCO、cache、refit |
| [repeatability_model.h](file:///disk228/huangjie/PanOmics/CosmicScanner/src/models/repeatability_model.h) | 113 | 2026-05-27 18:13 | 类定义、RepeatabilityMode 枚举、Options 结构体 |
| [gwas_app.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/src/gwas_app.cpp) | 101 | 2026-05-19 18:25 | 应用入口，`--repeatability-model` 分支路由 |
| [main.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/src/main.cpp) | 133 | 2026-05-27 18:22 | CLI 解析，所有重复力相关参数 |
| [config.h](file:///disk228/huangjie/PanOmics/CosmicScanner/src/config.h) | — | 2026-05-27 18:14 | GwasConfig 定义，含所有重复力 CLI 参数 |

### 3.2 VCE 底座（libCosmicCore）

| 文件 | 行数 | 最后修改 | 职责 |
|---|---|---|---|
| [vce.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/vce.cpp) | 3587 | 2026-05-27 22:34 | AI/EM-REML 核心，多组分支持，exact trace |
| [vce.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/vce.h) | 606 | 2026-05-22 10:14 | VCESolver 基类、AI_REML、策略枚举、SolverReport |
| [stcg_vce.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/stcg_vce.cpp) | 414 | 2026-05-27 21:44 | STCG 矩阵自由 VCE |
| [stcg_vce.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/stcg_vce.h) | 56 | 2026-05-12 18:41 | STCG_VCE 类定义 |
| [eigen_vce.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/eigen_vce.h) | 94 | 2026-05-19 21:03 | EigenVce（特征分解 REML，Fast 路线使用） |
| [eigen_vce.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/eigen_vce.cpp) | 231 | 2026-05-19 21:03 | EigenVce 实现 |
| [mme_builder.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/mme_builder.h) | 96 | 2026-04-07 20:01 | MMELHSBuilder、RandomComponent 定义 |
| [mme_builder.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/mme_builder.cpp) | 425 | 2026-04-07 20:01 | MME LHS 构建 |
| [matrix_free_pcg.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/matrix_free_pcg.h) | 75 | 2026-04-07 20:01 | PCG 求解器（STCG 路线使用） |

### 3.3 支撑模块（libCosmicCore）

| 文件 | 行数 | 最后修改 | 职责 |
|---|---|---|---|
| [covar_reader.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/covar_reader.h) | 72 | 2026-05-27 18:13 | 协变量读取，支持 FID:IID:TIME 记录级对齐 |
| [covar_reader.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/covar_reader.cpp) | 534 | 2026-05-27 18:13 | 协变量实现，含 alignRecords() |
| [pheno_reader.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/pheno_reader.h) | 76 | 2026-04-07 20:01 | 表型读取，含 loadLongitudinalRecords() |
| [pheno_reader.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/pheno_reader.cpp) | 286 | 2026-05-28 15:57 | 表型实现 |
| [design.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/design.h) | 61 | 2026-05-27 18:14 | GenRecord、FixedDesignG 定义 |
| [grm_builder.h](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/grm_builder.h) | 122 | 2026-05-20 20:07 | GRM 构建，含 computeSplit() |
| [grm_builder.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/libCosmicCore/src/grm_builder.cpp) | 2174 | 2026-05-20 20:07 | GRM 实现 |

**代码总量**：约 10,234 行（仅统计核心+VCE+支撑模块）

### 3.4 修改时间线

| 日期 | 修改内容 |
|---|---|
| 2026-04-07 | 基础架构：mme_builder、pheno_reader、matrix_free_pcg |
| 2026-04-20 | CosmicSolver 重复力 GS 测试脚本；CosmicScanner smoke test |
| 2026-04-23 | Repeatability_Closeout_CN.md |
| 2026-04-29 | edge case 测试脚本 |
| 2026-05-12 | stcg_vce.h |
| 2026-05-19 | gwas_app.cpp、eigen_vce、grm_builder |
| 2026-05-20 | grm_builder（computeSplit） |
| 2026-05-22 | vce.h（策略枚举、SolverReport）、开发总结文档 |
| 2026-05-24 | Methods 论文稿 |
| 2026-05-26 | 测试数据集（precision、snp_gradient、indiv_gradient） |
| 2026-05-27 | **关键修正日**：vce.cpp exact trace、stcg_vce PE 识别、covar_reader 记录级对齐、design.h GenRecord.time、main.cpp/config.h CLI 参数、repeatability_model.h |
| 2026-05-28 | **VCE 优化日**：repeatability_model.cpp fast subject-level REML、pheno_reader.cpp、RESULTS_SUMMARY.md、导师汇报文档、验证脚本 |

---

## 四、CLI 参数完整清单

| 参数 | 功能 | 默认值 | 状态 |
|---|---|---|---|
| `--repeatability-model` | 启用重复力 GWAS | off | ✅ |
| `--repeatability-mode auto/exact/fast/stcg` | 扫描策略 | auto | ✅ |
| `--gwas` | 执行 GWAS 扫描 | off | ✅ |
| `--no-loco` | 禁用 LOCO | off(默认 LOCO) | ✅ |
| `--fast-lmm` | Fast 均值聚合模式 | off | ✅ |
| `--repeatability-mixed-precision` | Float32 mixed precision | off | ✅ |
| `--repeatability-grm-cache <prefix>` | GRM 缓存 + mmap 复用 | "" | ✅ |
| `--repeatability-refit-validation` | per-SNP refit gold standard | off | ✅ |
| `--repeatability-refit-limit <N>` | refit SNP 数上限 | 0(全部) | ✅ |
| `--tpheno <col>` | 记录级表型 TIME 列号 | 0(禁用) | ✅ |
| `--pheno <file>` | 长格式表型文件 | — | ✅ |
| `--mpheno <col>` | 表型列号 | 1 | ✅ |
| `--covar <file>` | 记录级 categorical 协变量 | — | ✅ |
| `--qcovar <file>` | 记录级 quantitative 协变量 | — | ✅ |
| `--reml-vg/vpe/ve <float>` | VCE 初始值覆盖 | -1(自动) | ✅ |
| `--reml-maxit <n>` | VCE 最大迭代 | 50 | ✅ |
| `--reml-tol <float>` | VCE 收敛容差 | 1e-6 | ✅ |
| `--maf <float>` | MAF 过滤 | 0.01 | ✅ |
| `--hwe <float>` | HWE 过滤 | 0.0(禁用) | ✅ |
| `--make-grm-alg <n>` | GRM 算法 | 1(VanRaden) | ✅ |

---

## 五、测试数据集

### 5.1 Smoke Test 数据集

| 数据集 | 来源 | 个体数 | 记录数 | SNP 数 | 染色体 | 用途 |
|---|---|---|---|---|---|---|
| mouse_repeatability_subset | mouse_clean 子集 | 20 | 40 | 128 | 2 | CI smoke test |

生成方式：从 `examples/mouse_clean/` 中抽取 20 个有表型的个体，chr1 和 chr2 各 64 个 SNP，每个个体生成 2 条重复记录（原始值 + 偏移值）。

### 5.2 精度验证数据集

| 数据集 | 个体数 | 记录数 | SNP 数 | 用途 |
|---|---|---|---|---|
| precision_s100 | 100 | 385 | 500 | 精度验证、独立 R GLS 对照、sommer 对照 |
| snp_10k | 1092 | 4143 | 10,000 | SNP 梯度 benchmark |
| snp_100k | 1092 | 4143 | 100,000 | SNP 梯度 benchmark |

### 5.3 个体梯度数据集

| 数据集 | 个体数 | 用途 |
|---|---|---|
| s100 | 100 | 小规模基准 |
| s300 | 300 | 中规模 |
| s600 | 600 | 中大规模 |
| s1092 | 1092 | 全量 |

### 5.4 模拟回收数据集

使用 precision_s100 真实基因型，已知参数 Vg=1, Vpe=1, Ve=8, h2=0.1, r2=0.2，5 次独立模拟。

---

## 六、测试方法与覆盖

### 6.1 Smoke Test

**脚本**：[run_test_repeatability.sh](file:///disk228/huangjie/PanOmics/CosmicScanner/tests/run_test_repeatability.sh)

5 种模式覆盖：

| 模式 | 命令关键参数 | 验证内容 |
|---|---|---|
| exact | `--repeatability-model --repeatability-mixed-precision --repeatability-grm-cache --no-loco` | assoc.txt + null.txt 存在、n_records=40、n_subjects=20、scan_mode=exact_repeatability_global、Vg/Vpe/Ve/h2/r2 可解析、kernel_precision=float32_mixed、grm_cache_status=computed_save |
| exact_cached | 同 exact（cache 已存在） | grm_cache_status=mmap_load |
| loco | `--repeatability-model --repeatability-mixed-precision --repeatability-grm-cache` | scan_mode=exact_repeatability_loco、global+chr cache 文件完整 |
| loco_cached | 同 loco（cache 已存在） | chr cache mmap_load |
| fast | `--repeatability-model --fast-lmm --no-loco` | scan_mode=fast_repeatability_mean_model、agg.pheno 存在、mean_records>1 |

**状态**：✅ 通过

### 6.2 Edge Case Test

**脚本**：[run_test_repeatability_edge_cases.sh](file:///disk228/huangjie/PanOmics/CosmicScanner/tests/run_test_repeatability_edge_cases.sh)

8 个边界场景：

| 测试 | 预期行为 | 状态 |
|---|---|---|
| 空表型文件 | 正确拒绝（非零退出码） | ✅ |
| 全 NA 表型 | 正确拒绝 | ✅ |
| 单观测/个体 | 正常运行 | ✅ |
| MAF=0 边界 | 正常运行 | ✅ |
| MAF=0.5 边界 | 正常运行 | ✅ |
| categorical covariate | 正常运行 | ✅ |
| 缓存损坏 | 优雅重算 | ✅ |
| quantitative covariate | 正常运行 | ✅ |

**状态**：✅ 全部通过

### 6.3 CosmicSolver 侧重复力 GS 测试

**脚本**：[run_test_repeatability.sh](file:///disk228/huangjie/PanOmics/CosmicSolver/tests/run_test_repeatability.sh)

验证内容：
- `.vars` 含 Pe 和 Repeatability 行
- `.rand` 含 Pe 和 n_obs 列
- `.obsrand` 含 Fixed/Fitted/Residual
- 每个个体恰好 2 条观测
- 日志含 `AI-REML applied Pe absorption` 和 `Applied Pe absorption in final MME solve`

**状态**：✅ 通过

---

## 七、验证方法与验证结果

### 7.1 七线验证闭环

| 验证线 | 方法 | 数据 | 核心结果 | 状态 |
|---|---|---|---|---|
| 1. 协方差模型正确性 | 数学推导 | — | Var(y)=Z(G·Vg+I·Vpe)Z'+I·Ve 与经典重复力模型一致 | ✅ |
| 2. 记录级表型/协变量对齐 | 日志确认 | precision_s100, snp_10k | "Loaded 385/4143 record-level categorical records" | ✅ |
| 3. 独立 R full record-level GLS 逐 SNP 验证 | 独立 R 脚本从头解码 PLINK BED、重建 GRM、构建完整 record-level V、逐 SNP full GLS | precision_s100, 160 SNP | BETA corr=0.999787, SE corr=0.999961, P corr=0.998348 | ✅ |
| 4. 10k vs 100k 共有 SNP 一致性 | 不同 marker 集构建 GRM | 4195 共有 SNP | BETA corr=0.993280, SE corr=0.999974, FREQ corr=1.000000 | ✅ |
| 5. 可控模拟回收 | 已知 Vg=1,Vpe=1,Ve=8,h2=0.1,r2=0.2 | 5 次模拟 | h2 均值=0.0998（几乎无偏） | ✅ |
| 6. 独立 C++/Eigen REML + 外部 sommer::mmer 对照 | 独立 REML 实现 + R sommer | precision_s100 null model | Ve 和 r2 几乎一致，Vg/Vpe 进入外部 REML 支持区间 | ✅ |
| 7. 外部 sommer 10-SNP 逐 SNP refit | sommer::mmer 逐 SNP refit | 10 个 SNP | BETA corr=0.999281, SE corr=0.995166 | ✅ |

### 7.2 验证线 3 详细结果（最重要）

独立 R full record-level GLS 验证：

| 指标 | SNP 数 | 相关系数 | 平均绝对差 | 最大绝对差 |
|---|---:|---:|---:|---:|
| BETA | 160 | **0.999787** | 0.018414 | 0.133908 |
| SE | 160 | **0.999961** | 0.002585 | 0.034599 |
| P | 160 | **0.998348** | 0.012457 | 0.056428 |
| -log10(P) | 160 | **0.999168** | 0.015708 | 0.122489 |

验证脚本：`validate_precision_independent_gwas.R`

### 7.3 验证线 6 详细结果

precision_s100 null model 三方对照：

| 来源 | Vg | Vpe | Ve | h2 | r2 |
|---|---:|---:|---:|---:|---:|
| CosmicScanner fast VCE | 0.984842 | 2.39909 | 8.43302 | 0.0833415 | 0.286363 |
| 独立 C++/Eigen REML | 0.955067 | 2.42695 | 8.43288 | 0.0808358 | 0.286250 |
| 外部 sommer::mmer | 0.961298 | 2.42159 | 8.43262 | 0.0813590 | 0.286309 |

验证脚本：`validate_precision_reml.cpp`、`validate_precision_sommer.R`

### 7.4 验证线 5 详细结果

模拟回收（true: Vg=1, Vpe=1, Ve=8, h2=0.1, r2=0.2）：

| 指标 | 真值 | 5 次均值 | 5 次标准差 |
|---|---:|---:|---:|
| h2 | 0.100 | 0.0998 | 0.0257 |
| Vg | 1.000 | 0.9896 | 0.2934 |
| Ve | 8.000 | 8.1431 | 0.5630 |
| r2 | 0.200 | 0.1722 | 0.0381 |

### 7.5 VCE 四代迭代精度对比

| 数据集 | 版本 | Vg | Vpe | Ve | h2 | r2 |
|---|---|---:|---:|---:|---:|---:|
| precision_s100 | 旧版 | 1.14181 | 2.27241 | 8.43494 | 0.09636 | 0.28814 |
| precision_s100 | exact trace | 1.10853 | 2.30765 | 8.43854 | 0.09351 | 0.28817 |
| precision_s100 | exact trace + GLS fix | 1.14375 | 2.27574 | 8.43826 | 0.09646 | 0.28838 |
| precision_s100 | + VCE refinement | 1.00886 | 2.37434 | 8.43302 | 0.08538 | 0.28632 |
| precision_s100 | **+ fast VCE（推荐）** | **0.98484** | **2.39909** | **8.43302** | **0.08334** | **0.28636** |
| snp_10k | 旧版 | 1.25010 | 0.97056 | 11.0177 | 0.09443 | 0.16774 |
| snp_10k | **fast VCE（推荐）** | **1.26625** | **0.98911** | **11.0080** | **0.09547** | **0.17004** |
| snp_100k | 旧版 | 1.54698 | 0.76655 | 10.9842 | 0.11633 | 0.17398 |
| snp_100k | **fast VCE（推荐）** | **1.59556** | **0.69259** | **11.0123** | **0.11996** | **0.17204** |

---

## 八、时间与内存消耗

### 8.1 推荐结果（vcefast 路线）

| 数据集 | 记录数 | 个体数 | SNP 数 | Vg | Vpe | Ve | h2 | r2 | 耗时(s) | 峰值内存(MB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| precision_s100 | 385 | 100 | 500 | 0.985 | 2.399 | 8.433 | 0.083 | 0.286 | 3.6 | 27 |
| snp_10k | 4143 | 1092 | 10,000 | 1.266 | 0.989 | 11.008 | 0.095 | 0.170 | 377 | 654 |
| snp_100k | 4143 | 1092 | 100,000 | 1.596 | 0.693 | 11.012 | 0.120 | 0.172 | 398 | 688 |

### 8.2 VCE 路线性能对比

| 数据集 | 长 EM refinement | Fast subject-level REML | 加速比 |
|---|---:|---:|---:|
| precision_s100 | 4.917s | 3.639s | 1.35x |
| snp_10k | 573.3s | 377.3s | 1.52x |
| snp_100k | 632.7s | 397.9s | 1.59x |

### 8.3 内存对比

| 数据集 | 长 EM (MB) | Fast VCE (MB) | 差异 | 说明 |
|---|---:|---:|---:|---|
| precision_s100 | 24.0 | 27.2 | +3.2 | 小样本额外矩阵开销很小 |
| snp_10k | 656.5 | 654.3 | -2.1 | 大样本内存持平 |
| snp_100k | 679.5 | 687.8 | +8.3 | 100k SNP 下略增，同量级 |

### 8.4 个体梯度性能

| 样本量 | 记录数 | 耗时 |
|---|---|---|
| s100 | 385 | 1.76s |
| s300 | ~1200 | 11.64s |
| s600 | ~2400 | 53.71s |
| s1092 | 4143 | 199s |

### 8.5 每千 SNP 耗时

| 数据集 | SNP 数 | 总耗时 | 每千 SNP |
|---|---:|---:|---:|
| precision_s100 | 500 | 3.6s | 7.3s |
| snp_10k | 10,000 | 377s | 37.7s |
| snp_100k | 100,000 | 398s | 4.0s |

snp_100k 每千 SNP 耗时远低于 snp_10k，因为固定成本（GRM、VCE、预处理）被摊薄。

---

## 九、输出文件格式

### 9.1 assoc.txt

标准 GWAS 输出（默认格式或 `--out-ssf` EBI GWAS-SSF 格式）：

```
CHR  SNP  BP  A1  A2  FREQ  N  BETA  SE  P  LOG10P
```

### 9.2 repeatability.null.txt

```
Item    Value
n_records       4143
n_subjects      1092
Vg      1.26625
Vpe     0.989107
Ve      11.008
Vp      13.2634
h2      0.0954694
r2      0.170044
kernel_precision        float32_mixed
grm_cache_prefix       /path/to/cache/grm
grm_cache_status       mmap_load
scan_mode      exact_repeatability_loco
```

### 9.3 repeatability.cache.txt

```
Label   Status
global  mmap_load
chr_1   mmap_load
chr_2   mmap_load
```

### 9.4 Fast 模式额外输出

- `*.agg.pheno`：聚合后个体均值表型（FID IID PHENO_MEAN N_RECS）

### 9.5 Refit 模式额外输出

- `*.refit.assoc.txt`：per-SNP refit 结果
- `*.refit.meta.txt`：refit 元信息

---

## 十、文件废弃分析与清理建议

### 10.1 可清理的结果目录

| 目录 | 大小 | 最后修改 | 状态 | 建议 |
|---|---|---|---|---|
| `results_pure_rep_model/` | **2.1 GB** | 2026-05-27 | 旧版（个体级协变量 fallback） | ⚠️ **可清理**，已被 vcefix 替代 |
| `results_rrm_model/` | **929 MB** | 2026-05-27 | RRM 模型结果，非纯重复力 | ⚠️ **可归档**，不能当纯重复力汇报 |
| `results_pure_rep_record_time_smoke_20260527/` | 24 KB | 2026-05-27 | 中间 smoke | 可清理 |
| `tests/results/pig_repro_v2_gradient_opt/` | 7.9 MB | 2026-05-07 | 旧梯度优化 | 可归档 |
| `tests/results/pig_repro_v2_gradient/` | 56 KB | 2026-04-30 | 旧梯度 | 可清理 |
| `tests/results/pig_repro_gradient/` | 164 KB | 2026-05-17 | 旧梯度 | 可清理 |
| `tests/results/tmp_repeatability_fix3/` | 0 | 2026-05-07 | 空目录 | 可删除 |
| `tests/results/pig_snp_gradient/` | 13 MB | 2026-05-18 | 旧 SNP 梯度 | 可归档 |
| `tests/results/final_eval/` | 172 KB | 2026-05-24 | 旧评估 | 可归档 |

**潜在可回收空间**：约 3.0 GB

### 10.2 vcefix 目录内的中间版本文件

`results_pure_rep_record_time_vcefix/` 内包含四代 VCE 迭代的中间结果：

| 前缀 | 版本 | 是否推荐 | 建议 |
|---|---|---|---|
| `precision_s100` / `snp_10k` / `snp_100k` | 旧版（无 exact trace） | ❌ | 可清理 |
| `*_exacttrace` | 第二代 | ❌ | 可清理 |
| `*_exacttrace_glsfix` | 第三代 | ❌ | 可清理 |
| `*_exacttrace_glsfix_vceopt2` | 第三代+VCE refinement | ❌ | 可清理 |
| `*_exacttrace_glsfix_vcefast` | **第四代（推荐）** | ✅ | **保留** |

### 10.3 仍需保留的文件

| 文件/目录 | 原因 |
|---|---|
| `results_pure_rep_record_time_vcefix/*_vcefast.*` | 推荐最终结果 |
| `validate_precision_independent_gwas.R` | 独立 R GLS 验证脚本 |
| `validate_precision_reml.cpp` | 独立 C++ REML 验证 |
| `validate_precision_sommer.R` | sommer null model 对照 |
| `validate_precision_sommer_gwas.R` | sommer SNP refit 对照 |
| `run_record_time_validation.sh` | 验证运行脚本 |
| `RESULTS_SUMMARY.md` | 结果汇总 |
| `VCE_ROUTE_AUDIT.md` | VCE 路线审计 |
| `REPEATABILITY_GWAS_MENTOR_BRIEF.md` | 导师汇报文档 |
| `common_snp_consistency*.tsv` | 共有 SNP 一致性 |
| `sim_recovery/` | 模拟回收 |

### 10.4 文档重复度分析

| 文档 | 位置 | 内容 | 重复度 | 建议 |
|---|---|---|---|---|
| Repeatability_GWAS_开发总结与Methods初稿_CN.md | docs/ | 开发总结 + Methods 初稿 | 基础版 | 保留，作为论文底稿 |
| Repeatability_GWAS_Methods_论文稿_CN.md | docs/ | Methods 论文稿 | 进阶版 | 保留 |
| REPEATABILITY_GWAS_MENTOR_BRIEF.md | tests/data/.../vcefix/ | 导师汇报 | 最详细 | 保留 |
| Repeatability_Closeout_CN.md | tests/ | 收官说明 | 与 MENTOR_BRIEF 部分重叠 | 保留（GS 侧重） |
| RESULTS_SUMMARY.md | tests/data/.../vcefix/ | 结果汇总 | 数据层 | 保留 |
| VCE_ROUTE_AUDIT.md | tests/data/.../vcefix/ | VCE 审计 | 修正层 | 保留 |
| GWAS_AUDIT_20260527.md | tests/results/ | 全 GWAS 审计 | 含重复力段落 | 保留 |

---

## 十一、关键模型修正历史

| 修正项 | 旧问题 | 修正方案 | 修正日期 | 影响 |
|---|---|---|---|---|
| 记录级表型对齐 | 协变量可能退化为个体级 fallback | `--tpheno` 支持 FID:IID:TIME | 2026-05-27 | 修正固定效应校正 |
| 记录级协变量对齐 | categorical covariate 按 FID:IID 对齐 | 支持 FID:IID:TIME key | 2026-05-27 | 修正记录级协变量 |
| SNP 扫描公式 | subject-level shortcut 仅在平衡重复时等价 | 修正为 full record-level GLS 等价 | 2026-05-27 | 非平衡记录精度修复 |
| VCE trace 计算 | 随机迹估计对 Vpe 不稳定 | direct solver 下精确 trace | 2026-05-27 | Vpe 稳定性提升 |
| STCG PE 识别 | 大矩阵自由 VCE 可能遗漏 Vpe | 显式识别 sparse identity + record-level id_map | 2026-05-27 | STCG 路线修复 |
| 多组分 AI 更新 | 旧代码可能只更新 component 0 | 修复为正确更新所有组分 | 2026-05-22 | 多组分 VCE 正确性 |
| VCE refinement | 长 EM 300 次迭代慢 | subject-level exact REML refinement | 2026-05-28 | 速度提升 1.35x–1.59x |

---

## 十二、CosmicSolver 侧重复力 GS 闭环

| 能力 | 状态 |
|---|---|
| 标准重复力模型 `--pe` | ✅ |
| AI-REML / EM-REML / Exact / STCG VCE | ✅ |
| Pe absorption（AI-REML 迭代 + 最终 MME 求解） | ✅ |
| STCG 原型（V-implicit + stochastic trace + CG + EM warmup + damped AI） | ✅ |
| `.beta / .rand / .obsrand / .vars` 输出 | ✅ |

---

## 十三、当前限制与风险

| 风险项 | 严重程度 | 说明 |
|---|---|---|
| ASReml 未安装 | 中 | 尚未做 ASReml 对照，已有 sommer 对照 |
| sommer 全量逐 SNP refit 很慢 | 低 | 目前只做了 10-SNP 子集 |
| STCG 大规模 matrix-free 路线 | 中 | PE component 识别已修复，但需专门 benchmark |
| `--threads` 未由 main.cpp 解析 | 低 | 线程依赖 OPENBLAS/OMP 环境变量 |
| 时变协变量 | 中 | 当前稳定支持时不变协变量，时变尚未纳入正式主扫描 |
| 多组分标准误 | 低 | 通用 MME 路线尚未实现 |
| AI-REML 日志标签 | 低 | EM 更新迭代仍打印 AI 标签，属 cosmetic |
| BGEN/PGEN refit validation | 低 | 当前 refit 仅支持 PLINK |

---

## 十四、下一步落地优先级

### 第一优先级（尽快落地）

1. 固化 Exact + LOCO + cache + mixed precision 为正式主线
2. 在真实项目数据上形成一版可复现分析流程
3. 补充导师汇报所需图表和 benchmark 表
4. 清理废弃结果目录（可回收约 3 GB）

### 第二优先级

5. 扩展更大样本量与更多 SNP 数的 benchmark
6. 对 Fast 与 Exact 的误差做系统量化
7. 梳理 BGEN/PGEN 路径下的 I/O 性能差异
8. 补充 ASReml 对照（如可获取许可）

### 第三优先级

9. 推进 STCG 路线从原型走向正式 benchmark
10. 评估时变协变量纳入主扫描的改造方案
11. 准备英文稿件 Methods / Results / Supplementary
12. 实现 `--threads` CLI 解析

---

## 十五、验证脚本与结果文件索引

### 验证脚本

| 脚本 | 位置 | 用途 |
|---|---|---|
| validate_precision_independent_gwas.R | vcefix/ | 独立 R full record-level GLS 逐 SNP 验证 |
| validate_precision_reml.cpp | vcefix/ | 独立 C++/Eigen REML null model 验证 |
| validate_precision_sommer.R | vcefix/ | sommer::mmer null model 对照 |
| validate_precision_sommer_gwas.R | vcefix/ | sommer 10-SNP 逐 SNP refit |
| run_record_time_validation.sh | vcefix/ | 完整验证运行脚本 |

### 推荐最终结果文件

| 文件 | 数据集 |
|---|---|
| `precision_s100_exacttrace_glsfix_vcefast.assoc.txt` | 精度验证集 |
| `precision_s100_exacttrace_glsfix_vcefast.repeatability.null.txt` | 精度验证集 null |
| `snp_10k_exacttrace_glsfix_vcefast.assoc.txt` | 10k SNP 梯度 |
| `snp_10k_exacttrace_glsfix_vcefast.repeatability.null.txt` | 10k SNP null |
| `snp_100k_exacttrace_glsfix_vcefast.assoc.txt` | 100k SNP 梯度 |
| `snp_100k_exacttrace_glsfix_vcefast.repeatability.null.txt` | 100k SNP null |

### 验证结果文件

| 文件 | 用途 |
|---|---|
| `precision_s100_exacttrace_glsfix_vcefast.independent_gwas.tsv` | 独立 R GLS 结果 |
| `precision_s100_exacttrace_glsfix_vcefast.independent_gwas_compare.tsv` | CosmicScanner vs R 对比 |
| `precision_s100_exacttrace_glsfix_vcefast.independent_gwas_correlation.tsv` | 相关系数 |
| `precision_independent_reml.tsv` | 独立 C++ REML 结果 |
| `precision_sommer_reml.tsv` | sommer REML 结果 |
| `precision_s100_exacttrace_glsfix_vcefast.sommer_gwas_subset10.tsv` | sommer 10-SNP refit |
| `precision_s100_exacttrace_glsfix_vcefast.sommer_gwas_subset10.compare.tsv` | sommer 对比 |
| `common_snp_consistency_vcefast.tsv` | 10k vs 100k 共有 SNP |
| `sim_recovery/sim_recovery_summary.tsv` | 模拟回收汇总 |
