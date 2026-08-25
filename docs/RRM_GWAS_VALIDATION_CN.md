# CosmicScanner 纵向/动态性状 GWAS 方法与验证说明

## 1. 模块定位

纵向/动态性状 GWAS 面向同一个体在多个时间点、发育阶段或环境梯度下观测到的表型。传统 GWAS 通常估计 SNP 的平均效应，而本模块将 SNP 效应建模为时间或环境的函数，从而识别稳定型、阶段特异型和动态变化型遗传位点。

命令入口：

```bash
cosmicscanner \
  --bfile my_data \
  --long-pheno longitudinal_pheno.txt \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-order 2 \
  --rrm-curve-grid 21 \
  --rrm-top 100 \
  --rrm-ai-iter -1 \
  --rrm-sig-threshold 0.05 \
  --gwas \
  --out dynamic_assoc
```

## 2. 输入数据

纵向表型文件使用 long format：

```text
FID IID TIME Y
F001 I001 0.00 1.23
F001 I001 0.25 1.45
F001 I001 0.50 1.60
```

参数含义：

- `--long-pheno`：纵向表型文件。
- `--pheno`：非平衡 RRM 使用的 long-format 表型文件；与 `--rrm-unbalance` 搭配使用。
- `--tpheno`：时间或环境梯度列，按 FID/IID 之后的列编号计数。
- `--mpheno`：表型列，按 FID/IID 之后的列编号计数。
- `--rrm-order`：RRM basis 阶数；`legendre` 下为多项式阶数，`bspline` 下为 basis 数减 1。
- `--rrm-basis`：RRM 基函数类型；支持 `legendre` 和 `bspline`，默认 `legendre`。
- `--rrm-ai-iter`：RRM null 方差 AI refinement 迭代次数；`-1` 为自动模式，`0` 为 fast null，正整数为强制 AI 迭代次数。
- `--rrm-unbalance`：启用非平衡 longitudinal RRM；适用于个体间观测时间点或记录数不一致的 long-format 数据。
- `--rrm-loco-global-null`：快速近似 LOCO；复用 global null model，只按染色体分组扫描 SNP。默认不开启，默认 LOCO 为逐染色体重估 null model。
- `--tclass`：时间依赖分类协变量列名，内部自动编码。
- `--tcovar`：时间依赖连续协变量列名。

## 3. 随机回归模型

设第 `i` 个个体在时间或环境梯度 `t_ij` 下的表型为 `y_ij`。模型为：

```text
y_ij = x_ij' beta + sum_k phi_k(t_ij) a_ik + e_ij
```

其中：

- `x_ij` 为固定效应和协变量。
- `phi_k(t)` 为标准化时间上的 RRM 基函数，默认是 Legendre 多项式。
- `a_ik` 为个体 `i` 的第 `k` 个随机回归系数。
- `e_ij` 为残差。

随机回归系数的遗传协方差结构为：

```text
Cov(a_k, a_l) = K_kl G
```

其中 `G` 为个体基因组关系矩阵，`K` 为随机回归系数之间的遗传协方差矩阵。

默认基函数为 Legendre 多项式。可通过 `--rrm-basis bspline` 使用 open-uniform B-spline basis；此时 `--rrm-order n` 表示构建 `n+1` 个 basis 系数，spline degree 为 `min(3, n)`。GMA/GMAT 对齐建议继续使用默认 `legendre`。

实现中默认使用 `--rrm-ai-iter -1` 自动策略：小样本使用 AI refinement 细化 RRM null 方差分量；大样本使用 subject-level 随机回归系数协方差和残差方差初始化的 fast null，以避免对多个 observation-level 导数核进行高成本 AI 迭代。`*.rrm_joint.null.txt` 会记录 `AIIter` 和 `NullMode`，用于复现和审计。

### 3.1 非平衡 longitudinal RRM

非平衡模式用于每个个体观测时间点不完全一致、缺失时间点不同或重复记录数不同的 longitudinal 数据。输入仍为 long format，但命令使用 `--pheno` 指向表型文件，并显式加入 `--rrm-unbalance`：

```bash
cosmicscanner \
  --bfile my_data \
  --pheno longitudinal_pheno.txt \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-unbalance \
  --rrm-order 3 \
  --tclass Sex_class \
  --rrm-ai-iter 20 \
  --no-loco \
  --out unbalanced_rrm
```

模型包含三类协方差来源：

```text
V = Z (K0 ⊗ A) Z' + Z (P0 ⊗ I) Z' + I sigma2_e
```

其中 `K0` 是 Legendre 系数遗传协方差矩阵，`A` 是 GRM，`P0` 是永久环境或个体内非遗传协方差矩阵，`sigma2_e` 是观测层残差。早期只使用 `K0 + residual` 的模型会系统偏离 GMAT 非平衡模型；加入永久环境效应后，P 值排序和 top signal 明显改善。

生产实现不直接构建 observation-level 大矩阵。流程为：

- 对每个个体的非平衡记录拟合 Legendre BLUP/系数，得到 `N × (order+1)` 的个体层系数矩阵。
- 在个体层系数上运行 Kronecker multivariate LMM：`Y_coef = X B + G_snp Gamma + U + E`。
- 使用 `MvEigenVce` 做 multivariate AI-REML null 估计，使用 `MvLmmEngine` 做多自由度 Wald SNP scan。
- `--rrm-ai-iter <n>` 控制 BLUP coefficient LMM 的 AI-REML 迭代次数；推荐非平衡对齐验证使用 `20`。

该 BLUP-based 路径将计算规模从 `n_obs × n_obs` 的 observation-level 协方差矩阵降为 `N × (order+1)` 的个体层多性状 LMM，因此更适合大规模非平衡 longitudinal GWAS。LOCO 模式会按染色体重建 LOCO GRM 并重新估计 multivariate null model，支持 PLINK/PGEN/BGEN，但运行时间会随染色体数增加。

若需要快速 exploratory LOCO，可加入 `--rrm-loco-global-null`。该模式保留按染色体扫描 SNP，但复用 global BLUP null model，不再逐染色体重建 LOCO GRM 和重估 null model，因此速度接近多次 filtered scan；该模式是近似校正，不替代默认精确 LOCO。

## 3.2 GMA-compatible balanced fixreg 模式

除默认 dynamic RRM 外，CosmicScanner 提供显式兼容验证入口：

```bash
cosmicscanner \
  --bfile plink_subset \
  --pheno phe.balance.cosmic.txt \
  --mpheno-list 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 \
  --rrm-compatible gma-fixreg \
  --rrm-compatible-varcom gma_reference_varcom.var \
  --rrm-order 3 \
  --gwas \
  --out cs_compatible
```

该模式目标是复现 GMA `balance_longwas_fixreg` 的 balanced fixed-regression longitudinal GWAS 路径，而不是替代默认 dynamic RRM。输入要求为 PLINK `--bfile`、宽格式 balanced 表型 `--pheno` + `--mpheno-list`，并要求同一 PLINK 前缀下存在 GMA `cal_kin()` 生成的 `*.grm` 文件；`--rrm-compatible-varcom` 指向 GMA `balance_varcom()` 输出的 `.var` 文件。内部复现 GMA Legendre scaling、GRM eigen rotation、per-SNP fixed regression design、10 次 per-SNP AI/EM variance update 和多自由度 Wald 检验。

输出文件为 `*.rrm.gma_fixreg.assoc.txt`，主要字段包括 `CHR`、`SNP`、`BP`、`A1`、`A2`、`CHISQ`、`P`、`chi_val`、`p_val` 和 `BETA_L0...BETA_Lp`。其中 `chi_val/p_val` 用于和 GMA reference 输出直接合并比较。

## 4. SNP 动态效应

对每个 SNP，模型估计其在 Legendre 基上的效应向量：

```text
b = (b_0, b_1, ..., b_p)'
```

在任意时间点 `t`，SNP 的动态效应重构为：

```text
b(t) = phi(t)' b
```

其方差使用完整协方差矩阵计算：

```text
Var[b(t)] = phi(t)' Var(b) phi(t)
```

时间点特异检验为：

```text
z(t) = b(t) / SE[b(t)]
P(t) = P(ChiSq_1 >= z(t)^2)
```

全局动态效应检验使用多自由度 Wald 检验：

```text
W = b' Var(b)^-1 b
P_joint = P(ChiSq_df >= W)
```

## 5. 输出文件

主要输出：

- `*.rrm_joint.assoc.txt`：observation-level RRM 联合检验结果。
- `*.rrm_joint.null.txt`：null model 方差组分和随机回归系数协方差矩阵。
- `*.rrm.dynamic_curve.txt`：top SNP 在时间网格上的动态效应曲线。
- `*.rrm.dynamic_summary.txt`：top SNP 的关键时间、显著窗口和轨迹分类。
- `*.rrmcoef.txt`：个体层随机回归系数。
- `*.rrm.top_curve.txt`：基于系数辅助扫描的近似曲线输出，主要用于调试和兼容。
- `*.rrm_unbal.assoc.txt`：非平衡 RRM BLUP-based 联合检验结果，包含 `JOINT_CHISQ`、`JOINT_P`、`BETA_L0...BETA_Lp` 和 `SE_L0...SE_Lp`。
- `*.rrm_unbal.null.txt`：非平衡 RRM null model 估计；`K_*` 字段为 Legendre 系数遗传协方差矩阵，`PE_*` 字段为 BLUP coefficient LMM 的非遗传/永久环境合并协方差矩阵。
- `*.rrm_unbal.dynamic_curve.txt`：非平衡模式 top SNP 的动态效应曲线。
- `*.rrm_unbal.dynamic_summary.txt`：非平衡模式 top SNP 的关键时间、显著窗口和轨迹分类。
- `*.rrm.gma_fixreg.assoc.txt`：`--rrm-compatible gma-fixreg` 输出的 GMA-compatible balanced fixed-regression Wald scan 结果。

`*.rrm.dynamic_summary.txt` 关键字段：

| 字段 | 含义 |
|---|---|
| `JointP` | 多自由度动态效应联合检验 P 值 |
| `BestTime` | 时间点特异 P 值最小的位置 |
| `BestTimeP` | 最小时间点特异 P 值 |
| `BetaStart` / `BetaEnd` | 曲线起点和终点效应 |
| `DeltaBeta` | 终点效应减起点效应 |
| `BetaRange` | 曲线最大效应减最小效应 |
| `MaxAbsBeta` | 最大绝对效应 |
| `EffectRatio` | `BetaRange / MaxAbsBeta` |
| `PeakTime` | 最大绝对效应对应时间 |
| `Pattern` | 轨迹分类 |
| `SigThreshold` | 关键窗口使用的 P 值阈值 |
| `SigWindowStart` / `SigWindowEnd` | 显著时间窗口起止 |
| `SigGridPoints` | 显著网格点数 |

## 6. 轨迹分类

当前实现输出 9 类轨迹：

- `stable`
- `increasing`
- `decreasing`
- `stage_specific`
- `nonlinear_dynamic`
- `early_specific`
- `late_specific`
- `transient`
- `sign_switch`

分类依据包括起止效应、峰值时间、效应范围和效应变化比例。分类结果用于辅助解释，不替代正式统计检验；正式显著性仍以 `JointP` 和时间点特异 `P(t)` 为准。

## 7. 传统方法对比验证

动态模拟脚本：

```bash
tests/run_rrm_dynamic_simulation.sh
```

可选 LOCO smoke 验证：

```bash
tests/run_rrm_dynamic_simulation.sh --loco-smoke
```

该模式在 3 条染色体的模拟数据上额外运行一次 RRM dynamic LOCO 扫描，并检查 `LOCO_1`、`LOCO_2`、`LOCO_3` null-model scope 以及动态曲线和动态 summary 输出是否生成。

可选 PGEN 输入验证：

```bash
tests/run_rrm_dynamic_simulation.sh --pgen-smoke
```

该模式使用 `plink2` 将模拟 PLINK 数据转换为 PGEN，并通过 `--pfile` 运行一次 RRM dynamic smoke test，检查 joint association、动态曲线和动态 summary 输出，并确认动态 causal SNP 可检出。

可选 BGEN 输入验证：

```bash
tests/run_rrm_dynamic_simulation.sh --bgen-smoke
```

该模式使用 `plink2` 将模拟 PLINK 数据转换为 BGEN，并通过 `--bgen` + `--sample` 运行一次 RRM dynamic smoke test，检查 joint association、动态曲线和动态 summary 输出，并确认动态 causal SNP 可检出。

可选阶数敏感性验证：

```bash
tests/run_rrm_dynamic_simulation.sh --order-sweep
```

该模式额外运行 `--rrm-order 1`、`--rrm-order 2` 和 `--rrm-order 3`，检查动态输出完整性，并生成 `rrm_order_sweep.tsv` 汇总 5 个 causal SNP 在不同阶数下的排名、`JointP` 和轨迹分类。

该脚本模拟 5 类 causal QTL：

- stable
- sign-switch
- early-specific
- late-specific
- transient

并同时运行：

- RRM dynamic GWAS
- 平均表型 GWAS
- time0 GWAS
- time1 GWAS
- 多性状 GWAS

脚本生成 `comparison.tsv`，记录 causal SNP 在不同方法中的排名。当前验证结果显示，RRM dynamic GWAS 能同时恢复五类 causal SNP，并输出连续效应曲线、关键窗口和轨迹分类；平均表型 GWAS 会削弱 sign-switch 类型动态信号；端点 GWAS 只能反映离散时间点效应；多性状 GWAS 可检出部分离散多时间点效应，但不提供连续轨迹和关键窗口。

多 replicate 方法比较脚本：

```bash
tests/run_rrm_method_comparison.sh
```

默认配置为 `N=120`、`M=300`、`TIMES=5`、`REPS=5`，比较以下方法：

- RRM fast：`--rrm-ai-iter 0`
- RRM exact：`--rrm-ai-iter 8`
- 平均表型 GWAS
- time0 GWAS
- time1 GWAS
- 多性状 GWAS

输出文件包括：

- `method_metrics.tsv`：每个 replicate、每个方法、每个目标 SNP 的 rank、P 值、runtime、max RSS、动态曲线 RMSE 和 RRM 轨迹分类。
- `summary.tsv`：按方法汇总 dynamic power、stable power、null FPR、median dynamic rank、平均 runtime、平均内存和动态曲线 RMSE。
- `summary.md`：同一汇总结果的 Markdown 表格。

当前默认结果目录为 `tests/results/rrm_method_comparison_n120_m300_r5`。在该小规模设置下，RRM fast 与 RRM exact 均达到 dynamic power 1.0、stable power 1.0、null FPR 0、median dynamic rank 1；RRM fast 平均运行约 0.29 秒，RRM exact 约 2.71 秒，动态曲线 RMSE 均约 0.085。平均表型 GWAS 对 sign-switch dynamic SNP 的 power 为 0，说明均值汇总会抵消方向相反的动态效应。多性状 GWAS 可检出 dynamic SNP，但不输出连续效应曲线；若某些 SNP 被多性状扫描路径跳过，汇总中按 `NA` 记录。

GMA-compatible fixreg 对齐脚本：

```bash
tests/run_gma_compatible_fixreg_alignment.sh
```

该脚本使用 GMA 官方 balanced case，抽取默认 `N=120, M=20, TIME=16` 子集，先运行 GMA `cal_kin()`、`balance_varcom()` 和 `balance_longwas_fixreg()` 作为 reference，再运行 CosmicScanner `--rrm-compatible gma-fixreg`，最后合并逐 SNP `chi_val/p_val`。当前默认对齐结果为 `Pearson_log10P=0.9999999999999999`、`Spearman_P=1.0`、最大 P 值绝对差约 `3e-9`、最大 chi-square 绝对差约 `1.4e-7`，top SNP 均为 `JAX00241881`。

GMAT `longwas --unbalance` 对齐验证：

- 数据：GMAT mouse_long 非平衡 longitudinal example，1212 个体、11833 SNP、19392 条观测、16 个时间点，`--rrm-order 3`，包含时间依赖分类协变量 `Sex_class`。
- 命令：`--rrm-unbalance --tpheno 1 --mpheno 2 --tclass Sex_class --rrm-ai-iter 20 --no-loco`。
- 性能：CosmicScanner BLUP-based RRM 运行约 7.9 秒，峰值内存约 195 MB；早期 observation-level moment 路径约 1682 秒和 17.7 GB。
- 对齐：与 GMAT `longwas --unbalance` 合并 11833 个 SNP 后，autosome `Pearson r(-log10P)=0.9521`，`Spearman rho(P)=0.9401`（全序含 X 染色体为 0.9470 / 0.9393；差异与 X 染色体离群 SNP 的方法学来源见 `docs/RRM_GWAS_METHOD_AUDIT.md` §4.3）。
- top signal：CosmicScanner 和 GMAT 的最强 SNP 均为 `UNC18848064`，GMAT P 值约 `4.2e-6`，CosmicScanner P 值约 `8.1e-6`。
- 方差组分：BLUP-based AI-REML 估计的 `K0` 对角线约为 `[1.63, 0.24, 0.04, 0.02]`，接近 GMAT 参考 `[1.34, 0.21, 0.04, 0.02]`。

可复现脚本：

```bash
tests/run_gmat_unbalanced_alignment.sh
```

脚本默认读取 `/tmp/opencode/gmat_data` 下的 GMAT mouse_long 数据和 GMAT reference 输出，运行 CosmicScanner `--rrm-unbalance --no-loco`，并生成 `gmat_unbalanced_alignment.tsv`、`gmat_unbalanced_alignment_summary.tsv` 和 `gmat_unbalanced_alignment_summary.md`。脚本同时记录全序与 autosome 两组指标；默认为在 autosome SNP 上评估 top-SNP 一致性（通过阈值为 autosome `Pearson_log10P >= 0.94`、`Spearman_P >= 0.93` 且 autosome top SNP 一致）。可通过环境变量 `GMAT_DATA`、`RESULT_DIR`、`RRM_AI_ITER`、`MIN_PEARSON` 和 `MIN_SPEARMAN` 调整。

### 7.1 论文 Methods/Results 摘要表

RRM GWAS 模型和实现路径：

| 场景 | 输入格式 | Null model | SNP test | 主要输出 | 验证对象 |
|---|---|---|---|---|---|
| Balanced dynamic RRM | wide 或 long-format balanced records | observation-level RRM，`K0 ⊗ A + P0 ⊗ I + I sigma2_e` | basis 多自由度 Wald | `*.rrm_joint.assoc.txt`、动态曲线、动态 summary | 模拟 dynamic QTL、GMA-compatible fixreg |
| GMA-compatible balanced fixreg | wide balanced phenotype + GMA `.var` | 复现 GMA fixed-regression variance component path | GMA-compatible Wald | `*.rrm.gma_fixreg.assoc.txt` | GMA `balance_longwas_fixreg` |
| Unbalanced longitudinal RRM | long-format unbalanced records | BLUP coefficient Kronecker mvLMM，`Vg ⊗ A + Ve ⊗ I` | coefficient-level 多自由度 Wald | `*.rrm_unbal.assoc.txt`、动态曲线、动态 summary | GMAT `longwas --unbalance` |

对齐和性能结果：

| 比较 | 数据规模 | 准确性 | Top SNP | CosmicScanner 性能 | 参考/旧路径性能 |
|---|---:|---:|---|---:|---:|
| GMA balanced fixreg | N=120, SNP=20, TIME=16 | `Pearson_log10P=1.0`, `Spearman_P=1.0` | `JAX00241881` 一致 | 秒级 | GMA Python wrapper 秒级 |
| GMAT unbalanced | N=1212, SNP=11833, obs=19392, TIME=16 | autosome `Pearson_log10P=0.9521`, `Spearman_P=0.9401`（全序 0.9470/0.9393） | `UNC18848064` 一致 | 约 7.9 秒，195 MB | observation-level 旧路径约 1682 秒，17.7 GB |
| Dynamic QTL simulation | N=140, SNP=360, TIME=5 | stable/sign-switch/early/late/transient 均可检出或解释 | causal SNP 排名靠前 | 秒级 | 平均表型 GWAS 对 sign-switch 信号失效 |

新增可选优化验证：GMAT mouse_long 上 `--rrm-loco-global-null --rrm-ai-iter 3` smoke 输出 11833/11833 SNP 且 P 值均 finite，用时约 4.3 秒、峰值内存约 146 MB；默认精确 LOCO 仍保留逐染色体 null model 重估。`--rrm-basis bspline --rrm-order 3` smoke 输出 11833/11833 SNP 且 P 值均 finite，用时约 2.1 秒、峰值内存约 188 MB。

## 8. 可视化

绘图工具：

```bash
python3 tools/plot_rrm_dynamic.py \
  --curve dynamic_assoc.rrm.dynamic_curve.txt \
  --summary dynamic_assoc.rrm.dynamic_summary.txt \
  --comparison comparison.tsv \
  --top 5 \
  --outdir figures
```

输出：

- `rrm_dynamic_effect_curves.png/pdf`
- `rrm_dynamic_effect_heatmap.png/pdf`
- `rrm_dynamic_method_comparison.png/pdf`
- `rrm_dynamic_plot_summary.tsv`

## 9. 性能 benchmark

性能 benchmark 脚本：

```bash
tests/run_rrm_dynamic_benchmark.sh
```

默认配置为 `N=500`、`M=3000`、`TIMES=5`、`--rrm-order 2`、`--rrm-top 50`。脚本生成 PLINK 纵向模拟数据，运行 RRM dynamic GWAS，并输出：

- `summary.tsv`
- `summary.md`
- `runtime.tsv`
- `logs/cosmicscanner.log`

默认 benchmark 当前结果：500 个体、3000 SNP、5 个时间点、2500 条观测，自动 fast null 模式下用时 10.70 秒，峰值内存 338088 KB；5 类 causal SNP 均进入 top5，动态 sign-switch causal SNP 排名第 1。

非平衡 BLUP-based RRM 的 GMAT mouse_long 验证显示，`N=1212`、`M=11833`、`n_obs=19392`、`order=3` 时，`--no-loco` 运行约 7.9 秒、峰值内存约 195 MB。相同数据上 observation-level covariance 路径需要约 1682 秒和 17.7 GB；BLUP-based 路径约 213 倍加速、约 93 倍内存降低。

可通过环境变量调整规模，例如：

```bash
N=1000 M=10000 TIMES=5 tests/run_rrm_dynamic_benchmark.sh
```

## 10. 回归测试

RRM dynamic 相关测试：

```bash
tests/run_test_rrm.sh
tests/run_rrm_dynamic_simulation.sh
tests/run_rrm_dynamic_simulation.sh --loco-smoke
tests/run_rrm_dynamic_simulation.sh --pgen-smoke
tests/run_rrm_dynamic_simulation.sh --bgen-smoke
tests/run_rrm_dynamic_simulation.sh --order-sweep
tests/run_rrm_dynamic_benchmark.sh
tests/run_rrm_bspline_benchmark.sh
tests/run_rrm_method_comparison.sh
tests/run_gma_compatible_fixreg_alignment.sh
tests/run_gmat_unbalanced_alignment.sh
tests/run_gwas_core_regression.sh --rrm-only
```

其中 `--rrm-only` 会运行 RRM smoke test 和多类型动态 QTL 模拟验证。

## 11. 当前限制

- 当前默认基函数为 Legendre；B-spline 已作为可选 basis (`--rrm-basis bspline`) 提供 smoke 验证，但论文对齐结果仍基于 Legendre。
- 当前显著窗口按时间网格计算，暂未输出多个不连续显著窗口。
- 当前 RRM 模拟主验证使用 PLINK 输入，PGEN/BGEN 路径已有 smoke 验证。
- `--rrm-compatible gma-fixreg` 当前为 PLINK balanced validation path，依赖 GMA `.grm` 和 `.var` 文件；PGEN/BGEN 和内置 GMA null VCE 可后续扩展。
- 非平衡 `--rrm-unbalance` 默认推荐 BLUP-based coefficient LMM；LOCO 模式会对每条染色体重新构建 GRM 并重估 null model，因此比 `--no-loco` 明显更慢。
- `--rrm-loco-global-null` 是快速近似 LOCO，适合 exploratory scan；正式 LOCO 结果建议使用默认逐染色体 null 重估。
- 非平衡 BLUP null 输出中的 `PE_*` 表示 coefficient-level 非遗传/永久环境合并协方差，不是 observation-level 的单独残差方差矩阵。
- 轨迹分类为解释性规则，正式显著性以统计检验结果为准。
