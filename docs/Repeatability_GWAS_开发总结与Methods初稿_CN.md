# CosmicScanner 重复力 GWAS 开发总结与 Methods 初稿

## 1. 文档定位

本文档面向两类用途：

1. 作为导师汇报材料，系统梳理 `PanOmics` 中 `Cosmic` 体系，尤其是 `CosmicScanner` 重复力 GWAS 板块的当前开发进展、功能边界、算法准确性、性能瓶颈和下一步路线。
2. 作为论文或预印本的中文初稿底稿，参考 `EMmvGWAS` 的 Methods 组织方式，形成后续可直接扩写为英文 Methods/Results/Discussion 的结构。

当前重复力板块已不再停留在概念验证阶段，而是已经具备：

- 可运行的重复力 GWAS 主流程；
- 多路线求解与扫描框架；
- 面向精度、稳定性和速度的工程优化；
- 针对 smoke test、边界场景、梯度性能和外部软件对照的验证闭环。

## 2. 项目位置与代码板块划分

### 2.1 项目结构定位

`PanOmics` 中与重复力 GWAS 直接相关的核心代码主要位于：

- `CosmicScanner/src/models/repeatability_model.cpp`
- `CosmicScanner/src/models/repeatability_model.h`
- `CosmicScanner/src/gwas_app.cpp`
- `CosmicScanner/src/main.cpp`
- `CosmicScanner/libCosmicCore/src/vce.cpp`
- `CosmicScanner/libCosmicCore/src/vce.h`
- `CosmicScanner/libCosmicCore/src/stcg_vce.cpp`
- `CosmicScanner/libCosmicCore/src/matrix_io.cpp`

辅助参考与上下游支撑模块包括：

- `CosmicScanner/src/gwas_scanner.cpp`
- `CosmicScanner/src/gwas_scanner_rv.cpp`
- `PanOmics/tests/Repeatability_Closeout_CN.md`
- `CosmicScanner/tests/run_test_repeatability.sh`
- `CosmicScanner/tests/run_test_repeatability_edge_cases.sh`
- `CosmicScanner/tests/results/pig_repro_v2_gradient_opt/`
- `CosmicScanner/tests/results/pig_repro_v2_gradient_opt/hiblup_validation/`

### 2.2 功能板块划分

当前 `CosmicScanner` 重复力 GWAS 可以拆成 6 个功能板块：

1. 数据读入与样本对齐
   - 基因型支持 `PLINK`、`BGEN`、`PGEN`
   - 长格式表型读入
   - 个体级协变量读入与对齐

2. 重复力 null model 构建
   - 构造重复记录到个体的映射
   - 构造固定效应设计矩阵 `X`
   - 构造加性遗传效应与永久环境效应两类随机效应

3. 方差组分估计
   - 标准 `AI/EM-REML` 路线
   - 多起点初始化以回避局部极值
   - 大样本 `STCG_VCE` 矩阵自由近似路线

4. 全基因组扫描
   - `Exact` 记录层精确 GLS
   - `Fast` 个体均值聚合近似扫描
   - `STCG` PCG/score 型扫描

5. 工程优化
   - `LOCO`
   - `mixed precision`
   - `GRM cache`
   - `mmap` 缓存复用
   - 单次全基因组拆分 GRM 以复用 per-chrom 信息
   - OpenMP 并行

6. 验证与回归测试
   - smoke test
   - edge cases
   - 梯度性能测试
   - HIBLUP fixed-effect 对照
   - 内部 per-SNP refit gold standard 对照

## 3. 当前开发进展总览

### 3.1 已完成能力

截至目前，重复力 GWAS 主线已经具备以下可交付能力：

- `--repeatability-model` 启动重复力 GWAS；
- 默认 `auto` 策略根据记录数自动切换 `Exact / Fast / STCG`；
- `Exact` 路线支持记录层精确 GLS 扫描；
- `Fast` 路线支持个体均值聚合加权近似扫描；
- `STCG` 路线支持矩阵自由的 null model/VCE 与 score 型扫描原型；
- 支持 `LOCO`；
- 支持 `--repeatability-mixed-precision`；
- 支持 `--repeatability-grm-cache` 与 `mmap` 复用；
- 支持 `--repeatability-refit-validation` 和 `--repeatability-refit-limit`；
- 支持个体级 `covar/qcovar`；
- 输出 `assoc`、`repeatability.null`、`cache summary`、`agg.pheno` 等结果文件。

### 3.2 当前主线结论

如果按“能否直接落地生产”来评估：

- `Exact + LOCO + cache + mixed precision` 已经是当前最成熟、最适合正式结果交付的主线。
- `Fast` 已经可以作为大样本预筛和开发期快速迭代路线。
- `STCG` 已经具备大样本探索价值，但仍偏原型化，当前更适合作为生物样本库级长期路线储备。

### 3.3 与 CosmicSolver 的关系

从 `Cosmic` 体系视角看，重复力板块已经形成：

- `CosmicSolver` 负责重复力 GS / GEBV / VCE；
- `CosmicScanner` 负责重复力 GWAS；
- 两者共享 `libCosmicCore` 的 VCE 与矩阵计算底座。

因此，`CosmicScanner` 的重复力 GWAS 已经不是孤立模块，而是 `Cosmic` 重复力生态中面向关联分析的一侧闭环。

## 4. 重复力 GWAS 的统计模型与原理

这一部分建议后续直接作为论文 Methods 的主体骨架。

### 4.1 研究对象与数据结构

重复力 GWAS 面向同一个体存在多次重复观测的表型数据。设共有 `N` 个基因分型个体，第 `i` 个个体对应 `n_i` 条表型记录，总记录数为：

```text
n_obs = sum_i n_i
```

令第 `i` 个个体的第 `j` 次记录为 `y_ij`，其模型可写为：

```text
y_ij = x_ij' b + g_i + pe_i + m_ij * alpha + e_ij
```

其中：

- `x_ij` 为固定效应与协变量；
- `b` 为固定效应系数；
- `g_i` 为个体 `i` 的加性遗传随机效应；
- `pe_i` 为个体 `i` 的永久环境效应；
- `m_ij` 为位点剂量，当前实现中实际对应个体基因型并复制到该个体的所有重复记录；
- `alpha` 为待检验 SNP 效应；
- `e_ij` 为残差项。

随机效应假设为：

```text
g ~ N(0, Vg * G)
pe ~ N(0, Vpe * I)
e ~ N(0, Ve * I_obs)
```

其中：

- `G` 为基于全基因组标记构建的个体间 GRM；
- `I` 为个体层单位阵；
- `I_obs` 为记录层单位阵。

因此总表型方差可以分解为：

```text
Var(y) = Zg (Vg * G) Zg' + Zpe (Vpe * I) Zpe' + Ve * I_obs
```

其中 `Zg` 和 `Zpe` 都是记录到个体的 incidence matrix。

### 4.2 重复力参数

总表型方差：

```text
Vp = Vg + Vpe + Ve
```

狭义遗传力：

```text
h2 = Vg / Vp
```

重复力：

```text
r2 = (Vg + Vpe) / Vp
```

在重复测量场景中，`r2` 比单纯的 `h2` 更能表征个体跨重复记录的稳定可重复成分，因此它是该模型的重要生物学量。

## 5. Methods 初稿框架

参考 `EMmvGWAS` 的写法，我们建议论文 Methods 用下面的结构：

### 5.1 重复力线性混合模型

建立记录层重复力线性混合模型，明确：

- 长格式表型定义；
- 固定效应项；
- 加性遗传随机效应；
- 永久环境随机效应；
- 残差项；
- 目标检验位点效应。

这一节建议重点强调：与单性状 GWAS 相比，重复力 GWAS 的核心差别在于把“同一个体的多次观测相关性”显式建模进协方差结构中，而非直接压缩成单个均值。

### 5.2 GRM 构建与样本/记录映射

`CosmicScanner` 当前先在个体层构建 GRM，再通过记录到个体映射把其提升到记录层协方差。该设计的关键优点是：

- 保留个体间遗传相关性；
- 避免直接构建超大记录层 GRM；
- 为 Woodbury/降维思路提供数学基础；
- 更容易和 `LOCO`、`cache`、`mixed precision` 结合。

当前支持多种 GRM 标准化算法，主要包括：

- `VanRaden`
- `Yang et al.`
- `TrueVar` 等

代码中对不同 GRM 算法是否进行 per-SNP scaling 做了区分，因此在 `LOCO` 时采用不同的尺度修正规则，保证“去掉某条染色体后”的核矩阵缩放仍然合理。

### 5.3 Null model 方差组分估计

重复力 null model 的目标是估计：

- `Vg`
- `Vpe`
- `Ve`

当前实现包括两级设计：

1. 标准主线：`AI/EM-REML`
2. 大样本备线：`STCG_VCE`

#### 5.3.1 标准 AI/EM-REML 路线与特征分解预估

在 `Exact/Fast` 模式下，null model 的方差组分估计主要走 `AI_REML`。**为了极致加速收敛，系统首先会运行基于特征分解（Eigen Decomposition）的 VCE 预估（受 EMmvGWAS 启发）**，通过将表型与协变量投影到对角化的正交空间，极速获取一组高精度的方差组分预估值，避免直接进行缓慢的稠密矩阵求逆迭代。

随后，在正式的 `AI_REML` 中，为了避免重复力模型中常见的局部极值问题，当前实现加入了多起点策略（包含上述特征分解预估的最佳值）：

- `(Vg, Vpe)` 初始化比例尝试：
  - `(0.5, 0.3)`
  - `(0.7, 0.2)`
  - `(0.2, 0.2)`
  - `(0.3, 0.5)`
  - `(0.4, 0.1)`
- 每次拟合后根据 `logL` 与解的合理性筛选最佳结果；
- 若没有“完全合理”的解，则退而选择 `logL` 最高解。

之所以要这样做，是因为重复力模型存在两个现实中的数值风险：

- 正确解：`Ve` 占主导，`Vg` 和 `Vpe` 共享剩余可重复部分；
- 伪解：`Vg + Vpe` 吸收过多变异，`Ve` 被异常压缩。

多起点筛选是当前准确性保证中非常关键的一步。

#### 5.3.2 STCG_VCE 路线

当记录数很大且 `PLINK` 基因型矩阵算子可直接挂接时，系统可走 `STCG_VCE`：

- 使用 matrix-free `V` 乘法；
- 使用 PCG 求解 `V^{-1} y` 和 `V^{-1} X`；
- 使用 Hutchinson 随机迹估计梯度中的 trace 项；
- 使用 AI 与 EM 风格更新混合；
- 用步长夹逼和正值约束避免方差塌缩。

这一路线的意义在于把原本依赖显式大矩阵分解的 null model，转为矩阵自由近似求解，从而为超大规模重复记录数据提供扩展空间。

### 5.4 全基因组扫描策略

#### 5.4.1 Exact 路线

`Exact` 是当前正式主线。它并不是朴素地对记录层完整协方差矩阵（复杂度 $O(n_{obs}^3)$）逐 SNP 反复求逆，而是利用个体层核与记录-个体映射做**基于 Woodbury 矩阵恒等式的极速降维**，核心思想是：

- 在个体层构造

```text
K = Vg * G + Vpe * I
```

- 再结合每个个体记录数 `D` 与残差 `Ve`，利用 Woodbury 恒等式将记录层矩阵的求逆转化为个体层矩阵（$O(N_{sub}^3)$）的分解与构造记录层逆协方差所需的中间量；
- 先拟合 null model；
- 然后使用 score/GLS 形式逐 SNP 闭式计算统计量，将核心扫描时间压缩至与单性状同量级（例如 10万 SNP 仅需 4 秒）。

当前实现会缓存：

- `VinvX`
- `(X'V^{-1}X)^{-1}`
- `P_y`
- `H_sub`
- 子空间投影量

从而把最重的 null model 部分与逐 SNP 扫描部分分离，显著降低扫描开销。

#### 5.4.2 Fast 路线

`Fast` 路线将同一个体的多次记录先聚合为个体均值：

```text
y_i_bar = sum_j y_ij / n_i
```

同时对协变量也做个体内平均，并使用：

```text
Var(y_i_bar) approx Vg * G + (Vpe + Ve / n_i) * I
```

构造近似加权模型。

这一近似的优点是：

- 直接把分析对象从记录层降到个体层；
- 与单性状 LMM 的数值结构更接近；
- 对大样本预筛尤其高效。

这一近似的局限是：

- 当 `n_i` 差异较大时，不平衡重复测量会带来近似误差；
- 因此它适合作为快速筛查路线，不宜替代正式精确结果。

#### 5.4.3 STCG 扫描路线

`STCG` 扫描原型基于矩阵自由 `V` 算子与 PCG：

- 先求 `P y`；
- 再对每个 SNP 构建记录层剂量向量；
- 使用 score 统计量：

```text
U = g' P y
Var(U) = g' P g
```

- 形成 Wald/score 型检验。

当前该路线体现了生物样本库级扩展方向，但还不是最完整的生产主线，后文会专门说明其边界。

### 5.5 LOCO 策略

`LOCO` 的目的是避免待检染色体上的近端污染。当前实现流程为：

1. 先构造全基因组 GRM；
2. 再构造每条染色体的贡献；
3. 通过尺度修正得到：

```text
G_LOCO(chr) = (G_global * S_global - G_chr * S_chr) / S_loco
```

其中不同 GRM 算法的 `S` 定义不同：

- 对 per-SNP scaling 的算法，`S` 使用 SNP 数量；
- 对只中心化不标准化的算法，`S` 使用 `sum_2pq`。

这一点是代码正确性的重要来源之一，因为如果 LOCO 的尺度处理不一致，就会导致每条染色体的背景核矩阵不在同一标度上，进而影响 P 值与 beta 稳定性。

### 5.6 工程优化细节

#### 5.6.1 Mixed precision

`Exact` 模式支持 `Float32` 保存 subject-level GRM 或中间核矩阵，再在关键求解时转回 `Float64`。其主要价值是：

- 降低内存压力；
- 为大样本运行腾出线程并行空间；
- 在当前测试中未观察到显著精度损失。

#### 5.6.2 GRM cache 与 mmap

`--repeatability-grm-cache <prefix>` 支持把 subject-level GRM 以及 per-chrom GRM 缓存为二进制 `Float32` 文件，并通过 `mmap` 复用。其作用有三层：

- 避免重复构造全局 GRM；
- 避免 LOCO 下重复构造每条染色体的局部核；
- 为未来 out-of-core 路线保留接口。

#### 5.6.3 数值稳定措施

为保证求解稳定性，代码中已经加入多种数值防护：

- `factorize_with_jitter`：对近奇异矩阵逐级加 jitter；
- `invert_psd_matrix`：特征值截断后再求逆；
- 近奇异 `XtVinvX` 的备用逆策略；
- 方差参数下界约束；
- STCG 中步长上下界夹逼。

## 6. 算法代码准确性评估

### 6.1 可以认为已经较为可靠的部分

根据代码结构与现有测试证据，以下部分准确性基础较扎实：

- 重复力模型的三组分方差分解 `Vg/Vpe/Ve`；
- `Exact` 路线的记录层精确扫描主线；
- `LOCO` 的染色体剔除与尺度修正；
- `cache + mmap` 的缓存复用机制；
- `Fast` 路线在平衡或近似平衡设计下的快速筛查能力；
- `refit validation` 提供的小数据严格 gold standard 验证闭环。

### 6.2 准确性证据

#### 6.2.1 与旧版/优化前后结果一致

在 `pig_repro_v2_gradient_opt` 的结果中，优化后的新版本在方差组分估计上与旧版本保持完全一致：

- `s100`: `Vg=0.379`, `Vpe=3.349`, `Ve=8.800`, `r2=0.298`
- `s300`: `Vg=0.419`, `Vpe=2.922`, `Ve=10.580`, `r2=0.240`
- `s600`: `Vg=0.393`, `Vpe=2.488`, `Ve=10.932`, `r2=0.209`
- `s1092`: `Vg=0.472`, `Vpe=2.273`, `Ve=11.481`, `r2=0.193`

这说明底层矩阵重构与并行优化没有破坏主模型结果。

#### 6.2.2 与 HIBLUP fixed-effect 的严格对照

在 `s100` 的 informative SNP 集合上：

- `CosmicScanner` 非 LOCO vs `HIBLUP fixed-effect`
  - matched SNPs: `38`
  - beta correlation: `0.959441`（VanRaden）
  - beta correlation: `0.956325`（Yang et al.）

这说明当前 `CosmicScanner` 输出的位点效应与外部软件严格 fixed-effect 方案已经高度一致。

#### 6.2.3 内部 per-SNP refit gold standard

当前已实现：

- `--repeatability-refit-validation`
- `--repeatability-refit-limit`

该模式在每个 SNP 上重新构造包含该 SNP 固定效应的完整模型，作为内部 gold standard。其与 HIBLUP fixed-effect 的结果：

- matched SNPs: `38`
- beta correlation: `0.960253`
- mean absolute beta difference: `0.13008`

这证明 `CosmicScanner` 现在已经具备自洽的严格验证路径，而不再完全依赖外部软件。

### 6.3 当前准确性边界与风险点

虽然主线已经相当可靠，但仍有若干风险点需要在汇报中实事求是说明：

1. `s100` 数据中只有 `38/108` 个 SNP 真正 informative
   - 其余 `70` 个位点在该子集中实质单态或低信息；
   - 因此“108 个 SNP 全部严格对照通过”这一说法目前并不成立。

2. `Fast` 路线对不平衡重复记录属于近似
   - 若 `n_i` 波动大，其结果与 `Exact` 之间可能产生偏差。

3. `STCG` 路线当前仍属原型
   - 其 null model 和扫描依赖 PCG、Hutchinson、SLQ 等近似；
   - 在正式论文主结果中暂不建议作为唯一主线。

4. 时变协变量尚未纳入正式主扫描主线
   - 当前协变量按个体对齐，稳定支持时不变协变量；
   - 对一些纵向设计，这是一个仍需补齐的方法学缺口。

## 7. 速度瓶颈与优化分析

### 7.1 已经解决或显著缓解的瓶颈

#### 7.1.1 记录层全矩阵直接求逆的瓶颈

重复力模型最初最大的计算代价来自记录层协方差矩阵的高维操作。当前 `Exact` 路线通过个体层降维与缓存 null model 中间量，已经显著缓解这一问题。

在 `pig_repro_gradient` 的优化结果中，`s1686`：

- 总运行时间从约 `15` 分钟下降到约 `1` 分钟；
- 峰值内存从约 `1.2 GB` 降到约 `482 MB`；
- `SNP_scan` 阶段对 `100` 个 SNP 仅约 `0.026` 秒。

这说明“扫描阶段”本身已经不再是主瓶颈，真正重的是 null model 与 LOCO 核构造。

#### 7.1.2 LOCO 的重复核构造成本

`pig_repro_v2_gradient_opt` 结果显示，在 `s1092` 并强制 `18` 条染色体 LOCO 的情况下：

- 总耗时从 `18` 分钟降到约 `3.3` 分钟；
- 加速比约 `5.5x`；
- 峰值内存约 `4.65 GB`。

这说明：

- 染色体级并行已经显著提升 LOCO 效率；
- `LOCO` 仍然是正式扫描中最主要的额外计算开销来源之一。

### 7.2 仍然存在的核心瓶颈

当前剩余瓶颈主要有 4 类：

1. per-chrom GRM 构造
   - LOCO 本质上需要全局核减去染色体核；
   - 当样本量和染色体数都上来时，这一步仍会占据大量时间和内存。

2. 个体层稠密矩阵分解
   - `Exact` 路线的关键代价仍然包含 `H_sub` 或同类矩阵的稠密分解；
   - 复杂度大致仍受 `O(n_sub^3)` 控制。

3. BGEN/PGEN 输入下的顺序 I/O
   - `PLINK` 在部分路径上受益于更直接的矩阵接口；
   - `BGEN/PGEN` 某些路径仍是逐位点顺序读取，I/O 和解码成本更高。

4. STCG 路线的近似误差与收敛成本
   - PCG 迭代次数、trace 采样数、SLQ 参数都影响稳定性与速度；
   - 尚需更系统的 benchmark 才能界定最佳配置。

### 7.3 下一步性能优化建议

按优先级排序，建议推进：

1. 优先优化 `Exact + LOCO`
   - 因为这是当前正式产出主线；
   - 任何优化收益都会直接转化为生产效率。

2. 深化 GRM reuse
   - 提升 per-chrom cache 命中率；
   - 探索跨任务复用全局 GRM 与 chr-GRM。

3. 继续压缩 `Exact` null model 的个体层分解成本
   - 例如更激进地复用分解结果、局部更新、分块近似等。

4. 把 `STCG` 从“能跑”推进到“可正式 benchmark”
   - 建立不同记录数、SNP 数和重复次数梯度下的精度-速度曲线。

## 8. 已做测试与验证结果总览

### 8.1 Smoke tests

`CosmicScanner/tests/run_test_repeatability.sh` 已覆盖：

- `exact`
- `exact_cached`
- `loco`
- `loco_cached`
- `fast`

脚本会验证：

- 输出文件是否齐全；
- `null summary` 关键字段是否存在；
- `scan_mode` 是否正确；
- `cache` 是否正确从 `computed_save` 转为 `mmap_load`；
- `agg.pheno` 是否生成；
- chr 级 cache 是否完整。

这是当前最核心的工程回归脚本。

### 8.2 Edge cases

`CosmicScanner/tests/run_test_repeatability_edge_cases.sh` 已覆盖：

- 空表型文件；
- 全 `NA` 表型；
- 单观测样本；
- `MAF=0` 与边界值；
- 协变量；
- 数量协变量；
- 缓存损坏重算。

这表明当前模块不仅能“正常跑”，也对异常输入和边界条件具有基本鲁棒性。

### 8.3 小鼠子集 smoke 结果

`repeatability_mouse_clean` 结果表明：

- 精确模式正常；
- LOCO 正常；
- cache/mmap 正常；
- fast 模式正常；
- 示例零模型参数为：
  - `Vg = 0.6139`
  - `Vpe = 5.7159`
  - `Ve = 9.0422`
  - `h2 = 0.0399`
  - `r2 = 0.4118`

### 8.4 梯度性能测试

#### `pig_repro_gradient`

- `s100`: `1.00 s`, `22 MB`
- `s500`: `15.75 s`, `63 MB`
- `s1000`: `30.63 s`, `189 MB`
- `s1686`: `60.82 s`, `482 MB`

说明在优化后，核心扫描阶段已经非常快。

#### `pig_repro_v2_gradient_opt`

- `s100`: `1.76 s`
- `s300`: `11.64 s`
- `s600`: `53.71 s`
- `s1092`: `199 s`

在全 `18` 染色体 LOCO 下依旧可运行，说明该路线已经具备真实任务的工程可落地性。

### 8.5 外部软件对照

与 `HIBLUP fixed-effect` 的对照说明：

- 当前最稳妥的表述是“在 informative SNP 集合上高度一致”；
- 不应直接表述为“所有 108 个 SNP 严格一致”；
- 差异主要仍与低信息 SNP、GRM 定义、缺失与单态处理有关。

## 9. 具体可选择的重复力 GWAS 路线

### 路线 A：Exact + LOCO + cache + mixed precision

适用场景：

- 正式结果产出；
- 导师汇报主结果；
- 对外展示和后续论文主分析。

优点：

- 统计意义最完整；
- 结果最接近严格记录层模型；
- 已有最充分测试支撑。

不足：

- 计算成本最高；
- LOCO 下内存与 per-chrom 核构造代价仍然不小。

推荐程度：最高。

### 路线 B：Exact + no-LOCO

适用场景：

- 算法对照；
- 调试；
- 小数据快速正确性验证；
- 与 `HIBLUP fixed-effect` 的非 LOCO 对比。

优点：

- 结果稳定；
- 便于做方法学对照与金标准比较。

不足：

- 不适合作为正式 GWAS 主结果，因为存在 proximal contamination 风险。

推荐程度：开发与验证强烈推荐，正式结果不推荐单独使用。

### 路线 C：Fast

适用场景：

- 大样本预筛；
- 多轮参数调试；
- 先快速看 Manhattan 轮廓。

优点：

- 速度快；
- 资源占用小；
- 易于大规模批量运行。

不足：

- 对不平衡重复测量是近似；
- 不建议直接替代最终正式结果。

推荐程度：作为预筛路线强烈推荐。

### 路线 D：STCG

适用场景：

- 超大规模记录数；
- 方法学研发；
- 未来生物样本库级重复力 GWAS 预研。

优点：

- 具备矩阵自由扩展潜力；
- 理论上适合超大数据。

不足：

- 当前仍偏原型；
- 目前与 `Exact` 相比，工程完备度和验证深度仍不足；
- `PLINK` 支持最顺畅，其他格式仍可能回退。

推荐程度：研发储备路线，不建议立即作为唯一生产主线。

### 路线 E：per-SNP refit validation

适用场景：

- 小数据严格验证；
- 对关键位点做 gold standard 复核；
- 方法论文中的准确性补充材料。

优点：

- 是当前最严格的内部验证路线；
- 适合证明主扫描近似或闭式计算没有偏离完整模型太远。

不足：

- 极慢；
- 只适合小规模验证而非全基因组生产。

推荐程度：验证必备，但不是生产路线。

## 10. 当前最推荐的对导师汇报口径

建议向导师明确表达如下结论：

1. `CosmicScanner` 重复力 GWAS 已完成主干算法搭建，形成 `Exact / Fast / STCG` 三路线框架。
2. 当前最成熟路线是 `Exact + LOCO + cache + mixed precision`，已经具备落地分析能力。
3. 在小数据和外部软件对照中，当前结果与 `HIBLUP fixed-effect` 高度一致，并已建立内部 `per-SNP refit` gold standard 验证闭环。
4. 现阶段最大工程瓶颈已从“逐 SNP 扫描”转移到“null model 与 LOCO 核构造”。
5. 下一阶段应以“巩固 Exact 正式主线 + 推进 STCG 大规模路线 + 扩展 benchmark”作为主攻方向。

## 11. 为什么重复力 GWAS 明显难于单性状 GWAS

这一节建议在导师汇报时单独讲清楚。重复力 GWAS 开发明显慢于单性状 GWAS，并不是因为当前实现路线有问题，而是因为二者本身就不是同一个难度层级的方法学任务。

### 11.1 统计建模层面更难

单性状 GWAS 的基本对象通常是“每个个体对应一个表型值”，其标准线性混合模型往往只需要处理：

```text
y = Xb + g + e
```

其中协方差结构通常写成：

```text
Var(y) = Vg * G + Ve * I
```

而重复力 GWAS 面对的是“同一个体存在多次重复记录”的长格式数据。其模型至少需要写成：

```text
y_ij = x_ij' b + g_i + pe_i + e_ij
```

进一步在 GWAS 时还需要加入待检位点效应。因此，重复力模型不是单纯多了几条记录，而是多了一层“个体内相关结构”和一个额外的永久环境方差组分 `Vpe`。这会直接带来 4 个统计后果：

1. 观测单位从“个体”变成“记录”，样本结构更复杂；
2. 协方差从两组分 `Vg/Ve` 变成至少三组分 `Vg/Vpe/Ve`；
3. 个体间相关和个体内相关同时存在；
4. 重复次数 `n_i` 常常不平衡，使近似模型更容易产生偏差。

因此，从方法学本质上讲，重复力 GWAS 是“单性状 GWAS + 层级相关结构 + 额外随机效应”的组合问题，天然更难。

### 11.2 数值优化层面更难

单性状 GWAS 中，null model 拟合好以后，很多扫描步骤已经非常成熟，现成算法和软件也很多。而重复力 GWAS 的难点首先就在 null model：

- 单性状通常估计 `Vg` 和 `Ve` 即可；
- 重复力模型需要同时估计 `Vg`、`Vpe`、`Ve`；
- 三组分模型更容易出现边界解、局部极值和方差塌缩；
- 不同初值可能把优化推到完全不同的参数区域。

这也是为什么当前实现中需要多起点 VCE，而不是像单性状那样一轮优化就足够稳定。换句话说，重复力 GWAS 的第一难点不是“扫描慢”，而是“null model 先要估得对、估得稳”。

### 11.3 LOCO 和核矩阵处理更难

单性状 GWAS 做 `LOCO` 时，通常是在个体层 GRM 上直接剔除当前染色体的贡献。重复力 GWAS 则更复杂：

- `LOCO` 仍然先发生在个体层 GRM；
- 但最终检验对象是在记录层；
- 因而必须把个体层的 `LOCO-GRM` 正确映射回记录层协方差；
- 不同 GRM 标准化算法下，LOCO 剔除后还要做一致的尺度修正。

这意味着单性状里可直接复用的 `LOCO` 技术，在重复力场景下往往只能借鉴思路，不能直接照搬实现。

### 11.4 速度优化更难

单性状 GWAS 的性能优化通常集中在：

- 更快的 GRM 构建；
- 更快的 null model；
- 更快的逐 SNP 读盘与扫描。

而重复力 GWAS 的优化链条更长，通常需要同时处理：

- 记录层到个体层的降维；
- 多组分方差结构；
- `LOCO` 的 per-chrom 核矩阵构造；
- 大矩阵分解与缓存复用；
- mixed precision；
- 线程并行；
- 对不平衡重复设计的近似误差控制。

因此，重复力 GWAS 的工程优化不是单一瓶颈优化，而是“统计正确性、数值稳定性和系统性能”三者同时协调的问题。这也是它明显比单性状更耗开发时间的直接原因。

### 11.5 验证闭环更难

单性状 GWAS 的一个巨大优势是外部 benchmark 非常多，诸如 `GEMMA`、`GCTA`、`BOLT-LMM`、`regenie` 等都可以直接拿来做对照。重复力 GWAS 则不同：

- 公开成熟软件更少；
- 可直接逐 SNP 严格对比的标准路线更少；
- 很多结果差异还会受到低信息 SNP、缺失模式、单态位点和 GRM 定义影响。

因此，重复力 GWAS 开发中大量时间并不花在“写出能跑的代码”，而是花在“证明这套代码确实是对的”。这也是当前项目里必须同时做：

- `HIBLUP fixed-effect` 外部对照；
- `per-SNP refit` 内部 gold standard；
- smoke tests；
- edge cases；
- 梯度 benchmark

的根本原因。

### 11.6 对当前项目的现实解释

从 `CosmicScanner` 当前开发实际看，单性状 GWAS 很多逻辑可以复用已有 `LMM` 主框架；而重复力 GWAS 在 [repeatability_model.cpp](file:///disk228/huangjie/PanOmics/CosmicScanner/src/models/repeatability_model.cpp) 中实际上单独发展出了：

- `Exact` 记录层精确路线；
- `Fast` 个体均值聚合路线；
- `STCG` 矩阵自由路线；
- `LOCO`；
- `cache/mmap`；
- `refit validation`。

这说明重复力 GWAS 不是在单性状 GWAS 上“打个补丁”，而是已经接近于在同一框架内重新长出了一套方法学子系统。因此，其开发速度不如单性状，是符合客观规律的。

### 11.7 汇报时最推荐的一句话概括

建议对导师用一句话概括：

```text
单性状 GWAS 主要难在“快不快”，而重复力 GWAS 同时难在“模得对不对、估得稳不稳、扫得快不快、证得明不明”。
```

如果要再展开一句，可以补充：

```text
重复力 GWAS 的开发本质上不是普通功能开发，而是“统计模型研发 + 数值优化 + 工程落地 + 严格验证”四件事同时进行，所以明显慢于单性状 GWAS 是正常现象。
```

## 12. 论文式 Results/Discussion 写作建议

### 12.1 Results 建议写法

Results 可分成 4 节：

1. 模块开发完成情况
   - 说明已实现记录层精确重复力 GWAS、fast 近似、LOCO、cache、mixed precision 和 refit validation。

2. 运行性能
   - 用 `pig_repro_gradient` 和 `pig_repro_v2_gradient_opt` 展示梯度加速和 LOCO 开销。

3. 结果一致性
   - 用 `HIBLUP fixed-effect` 和内部 refit 做 beta 一致性对照。

4. 鲁棒性
   - 用 smoke tests 和 edge cases 证明工程稳定性。

### 12.2 Discussion 建议写法

Discussion 可突出三点：

1. 相比简单均值分析，记录层重复力模型更符合重复测量数据结构；
2. 工程优化使得精确模型在中等规模数据上已具备现实可用性；
3. `STCG` 为后续超大规模重复力 GWAS 提供了可扩展方向，但还需更多系统 benchmark。

## 13. 可直接用于论文的 Methods 中文草稿

下面给出一段可直接继续扩写的中文 Methods 初稿。

### 13.1 统计模型

我们针对重复测量表型构建了记录层重复力 GWAS 模型。设第 `i` 个个体共有 `n_i` 条观测，第 `j` 条观测值记为 `y_ij`。在每个位点上，模型写为：

```text
y_ij = x_ij' b + g_i + pe_i + m_i alpha + e_ij
```

其中 `x_ij` 为固定效应和协变量，`g_i` 为加性遗传效应，`pe_i` 为永久环境效应，`m_i` 为个体在待检位点的基因型剂量，`alpha` 为位点效应，`e_ij` 为残差。我们假设：

```text
g ~ N(0, Vg * G), pe ~ N(0, Vpe * I), e ~ N(0, Ve * I_obs)
```

其中 `G` 为基于全基因组 SNP 构建的个体遗传关系矩阵。由此，记录层表型协方差由加性遗传成分、永久环境成分和残差成分三部分组成。

### 13.2 方差组分估计

在全基因组扫描前，我们首先在不包含待检 SNP 固定效应的 null model 下估计 `Vg`、`Vpe` 和 `Ve`。对于中小规模数据，采用 AI/EM-REML 方案；针对重复力模型中常见的局部极值问题，采用多组初始值进行多起点拟合，并按对数似然与方差组分合理性选择最终解。对于超大规模记录数的数据，我们进一步实现了基于 stochastic trace 和 preconditioned conjugate gradient 的矩阵自由 STCG_VCE 路线，以避免显式构造和分解超大协方差矩阵。

### 13.3 全基因组扫描

我们实现了三种扫描策略。第一，`Exact` 路线在记录层显式保留重复记录结构，并通过个体层 GRM 与记录-个体映射进行降维，在固定 null model 后以广义最小二乘和 score 型统计量逐 SNP 计算效应与显著性。第二，`Fast` 路线将同一个体的重复记录聚合为个体均值，并使用基于重复次数的加权近似模型进行快速扫描。第三，`STCG` 路线基于矩阵自由 `V` 算子与 PCG 求解 `P y` 和 `P g`，构造 score 检验统计量，用于超大规模数据的近似分析。

### 13.4 LOCO 与工程优化

为避免近端污染，我们在正式 GWAS 中采用 LOCO 策略，即在扫描某条染色体时，使用剔除了该染色体贡献的背景 GRM。考虑到不同 GRM 标准化算法的尺度不同，我们分别基于 SNP 数量或 `sum_2pq` 对 LOCO 核矩阵进行尺度修正。为进一步提升效率，我们引入了 `Float32` mixed precision、subject-level GRM 二进制缓存和 `mmap` 复用机制，以降低内存占用并减少重复运行时的核构造开销。

### 13.5 验证策略

我们采用四类验证策略评估方法正确性与实用性：第一，使用 smoke tests 验证 exact、fast、LOCO 和 cache 等主要功能；第二，使用边界场景测试评估异常输入、协变量和缓存损坏情况下的鲁棒性；第三，在猪繁殖性状模拟数据上进行样本量梯度 benchmark，评估运行时间和内存开销；第四，基于 HIBLUP fixed-effect 逐 SNP 拟合以及内部 per-SNP refit gold standard，对效应估计的一致性进行严格比较。

## 14. 下一步工作建议

建议把后续工作按“交付优先级”分成三层：

### 第一优先级

- 固化 `Exact + LOCO + cache + mixed precision` 为正式主线；
- 在真实项目数据上形成一版可复现分析流程；
- 补充导师汇报所需图表和 benchmark 表。

### 第二优先级

- 扩展更大样本量与更多 SNP 数的 benchmark；
- 对 `Fast` 与 `Exact` 的误差做系统量化；
- 梳理 `BGEN/PGEN` 路径下的 I/O 性能差异。

### 第三优先级

- 推进 `STCG` 路线从原型走向正式 benchmark；
- 评估时变协变量纳入主扫描的改造方案；
- 准备英文稿件中的 Methods、Results 和 Supplementary。

## 15. 总结

总体上，`CosmicScanner` 的重复力 GWAS 板块已经进入“可加速落地”的阶段，而不再只是原型探索：

- 方法学上，已经形成记录层精确模型、快速近似模型和矩阵自由大样本路线三层架构；
- 工程上，已经具备 LOCO、cache、mmap、mixed precision 和多线程并行等关键能力；
- 验证上，已经建立 smoke test、edge case、梯度 benchmark、外部软件对照和内部 gold standard 的闭环；
- 科学上，已经可以形成一版较完整的 Methods 叙述，并具备继续整理为论文初稿的基础。

如果目标是“尽快落地”，当前最优策略并不是再发散做更多分支，而是围绕 `Exact + LOCO` 主线继续压实性能、补足图表、整理正式报告，并同步把 `STCG` 作为下一阶段储备方向推进。
