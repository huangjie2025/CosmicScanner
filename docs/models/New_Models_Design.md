# CosmicScanner 新模型设计文档

## 多性状 GWAS Methods 草稿（已完成模块）

CosmicScanner 的多性状 GWAS 模块实现了面向多个相关性状的 multivariate linear mixed model。对 `n` 个个体和 `d` 个性状，模型写作 `Y = XB + g + e`，其中 `g ~ MN(0, K, Vg)`，`e ~ MN(0, I, Ve)`。`K` 为基因组关系矩阵，`Vg` 和 `Ve` 分别为遗传和残差性状协方差矩阵。程序对 `K` 做特征分解，将表型和协变量旋转到特征空间，并通过 EM 或 AI-REML 估计 null model 方差组分。

快速 `emmv` 路径复用 null model 方差结构进行 genome-wide SNP scan。每个位点对 genotype dosage 做样本子集化和缺失均值填补，在旋转混合模型空间估计跨性状 SNP 效应，并使用多自由度 Wald statistic 检验所有性状 SNP 效应同时为零的联合原假设。`--mv-compatible emmvgwas` 提供与 EMmvGWAS R 包一致的验证路径，包括 `Z = dosage - 1` 编码、trace-normalized GRM、`lambda = G R^{-1}` profiling、uncentered genotype scan 以及官方公式中的 general matrix inverse。

`gemma` 和 `gemma-full` 路径实现 GEMMA-style multivariate mixed-model validation。`gemma-full` 输出 Wald、LRT 和 Score p 值；默认使用 null-component fast route，也可通过 `--mv-refit-p` 对部分位点执行逐 SNP variance-component refit。

Adaptive hybrid 路径用于在接近 GEMMA-full 精度的同时降低计算量。Stage 1 使用 EMmvGWAS-style scan 做全基因组筛选；启用 LOCO 时按染色体构建 LOCO GRM，估计 chromosome-specific multivariate null model，并只扫描当前染色体。Stage 1.5 可选对 Stage 1 top SNP 做 Score-only audit，产生 null-based `P_SCORE`。Stage 2 对三类候选的并集执行 GEMMA-full refinement：Stage 1 p-threshold、Stage 1 topK、以及 `abs(-log10(P_SCORE) - -log10(P_WALD))` 最大的 Score/Wald discrepancy SNP。LOCO-hybrid Stage 2 按候选 SNP 所在染色体分组，复用对应 LOCO GRM 做 per-chromosome GEMMA-full refinement，并 append 到统一 refined association file。

实现上，快速多性状 scan 对 PLINK/PGEN/BGEN 使用 block-wise genotype streaming；GRM 构建和 SNP scan 共享 `extract`、`exclude`、染色体和 bp range filters，避免 GRM SNP set 与 scan SNP set 不一致。Hybrid refinement 当前限定 PLINK 输入，因为它需要 repeated selective SNP access 和 per-chromosome LOCO refinement。

验证包括 EMmvGWAS R 官方 RIL example、GEMMA 0.98.5 multivariate LMM 对标、多性状 LOCO/no-LOCO validation matrix、adaptive hybrid smoke，以及系统 benchmark。默认系统 benchmark 比较 EMmvGWAS R、GEMMA、CS-fast、CS-hybrid 和 CS-LOCO-hybrid，并输出 runtime、memory、p-value correlation、top overlap、causal recovery 以及图表。


## 前言

本文档记录 CosmicScanner GWAS 模块的新模型开发设计。现有 7 种模型/模式
（单性状 LMM、重复力、RRM、多性状、罕见变异、Biobank、SPA）
已覆盖生猪育种约 60-70% 的应用场景。

为覆盖 95% 以上生产需求，新增以下 6 个模型：

| # | 模型 | 优先级 | 主要应用 |
|---|---|---|---|
| 1 | MaternalModel | P0 | 出生重/断奶重等母体效应性状 |
| 2 | SingleStepModel | P0 | 系谱+基因组联合 ssGBLUP/ssGWAS |
| 3 | ThresholdModel | P1 | 死亡率/抗病力等二分类性状 |
| 4 | CountModel | P2 | 产仔数/乳头数等计数性状 |
| 5 | SurvivalModel | P2 | 母猪使用年限/淘汰风险 |
| 6 | GxEModel | P3 | 跨环境反应范/基因型-环境互作 |

---

## 1. MaternalModel — 母体效应模型

### 1.1 统计模型

```
y = Xb + Z_d a_d + Z_m a_m + Z_c c + e

a_d ~ N(0, K_d · σ²_ad)        直接加性遗传效应
a_m ~ N(0, K_m · σ²_am)        母体遗传效应
cov(a_d, a_m) = K_dm · σ_adm   直接-母体遗传协方差
c   ~ N(0, I · σ²_c)           母体永久环境效应
e   ~ N(0, I · σ²_e)
```

方差组分：`σ²_ad, σ²_am, σ_adm, σ²_c, σ²_e`

### 1.2 数据结构映射

每条记录映射到：

```
record → animal_id (FID:IID)
record → dam_id    (母体ID, 来自--maternal-id列)
```

### 1.3 CLI 接口

```
--maternal-model             启用母体效应模型
--maternal-id <col/name>     母体ID列 (默认: 从pheno第3列)
--maternal-pe                启用母体永久环境效应
--no-direct-maternal-cov     固定直接-母体协方差为0
--maternal-grm <file>        预计算母体GRM (可选)
--direct-grm <file>          预计算直接GRM (可选)
--maternal-comp-maxit <n>    多组分VCE最大迭代 (默认: 100)
```

### 1.4 实现路线

```
v1.0: EM-REML 多组分VCE + SNP固定效应score test
v2.0: AI-REML加速 + IDUL扩展
v3.0: 大规模 PCG + average trace approximation
```

### 1.5 关键代码模块

| 模块 | 位置 | 说明 |
|---|---|---|
| MaternalModel | `models/maternal_model.h/.cpp` | 顶层面板 |
| MultiComponentVCE | `libCosmicCore/src/eigen_vce.h` (扩展) | 多组分AI-EM REML |
| GRM builder | `libCosmicCore/src/grm_builder.cpp` | 两组GRM |

---

## 2. SingleStepModel — ssGBLUP / ssGWAS

### 2.1 统计模型

单步关系矩阵：

```
H⁻¹ = A⁻¹ + [0    0
              0  G⁻¹ - A₂₂⁻¹]
```

```
A  = 系谱关系矩阵 (pedigree)
A₂₂ = 有基因型个体的系谱子矩阵
G  = 基因组关系矩阵 (VanRaden法)
H  = 单步关系矩阵
```

GWAS模型：

```
y = Xb + g_j β_j + a + e
a ~ N(0, H · σ²_g)
```

### 2.2 系谱文件格式

标准三列格式：

```
ID SIRE DAM
1   0    0
2   0    0
3   1    2
...
```

支持缺失值：`0` / `NA` / `.`

### 2.3 CLI 接口

```
--single-step               启用ssGBLUP/ssGWAS
--ped <file>                系谱文件
--genotyped-list <file>     有基因型个体列表
--h-tau <float>             τ参数 (默认: 1.0)
--h-omega <float>           ω参数 (默认: 1.0)
--blend-g <float>           G = (1-α)G + αA₂₂ (默认: 0.05)
--ssf-ebv                   输出GEBV (默认: 只输出GWAS)
```

### 2.4 A⁻¹ 构建 (Henderson方法)

```
规则:
- 父母已知: 对(i,i),(i,sire),(i,dam),(sire,sire),(sire,dam),(dam,dam)贡献
- 父母未知: 对(i,i)贡献1
- 一个父母已知: 贡献4/3因子
```

### 2.5 实现路线

```
v1.0: A⁻¹ sparse构建 + G + A₂₂ + H⁻¹ matrix-free operator
v1.1: PCG求解 + 单步GWAS marker fixed-effect scan
v2.0: SNP效应 backsolve from GEBV
v3.0: 大规模系谱(>100万个体)支持
```

### 2.6 关键代码模块

| 模块 | 位置 | 说明 |
|---|---|---|
| SingleStepModel | `models/single_step_model.h/.cpp` | 顶层面板 |
| Pedigree | `libCosmicCore/src/pedigree.h/.cpp` | 系谱解析+A⁻¹ |
| H-inverse operator | `libCosmicCore/src/single_step_matrix.h/.cpp` | H⁻¹x 算子 |

---

## 3. ThresholdModel — 阈值模型

### 3.1 统计模型

Probit link 阈值模型：

```
y_i ∈ {0, 1}

l_i = Xb + Za + e     (潜变量 liability)
y_i = 1 if l_i > 0

P(y_i = 1 | u) = Φ(Xb + Za)

a ~ N(0, K · σ²_g)
e ~ N(0, I)
```

### 3.2 PQL 算法

每轮：

```
1. η = Xb + Zu                      (线性预测子)
2. μ = Φ(η)                          (拟合概率)
3. dμ/dη = φ(η)                     (导数)
4. z = η + (y - μ) / (dμ/dη)       (工作变量)
5. w = (dμ/dη)² / (μ(1-μ))         (权重)
6. 解加权LMM: z = Xb + Za + e*     (复用现有LMM)
7. 更新方差组分                      (复用EigenVce)
8. 收敛检查
```

### 3.3 CLI 接口

```
--threshold-model           启用阈值模型
--link probit               链接函数 (默认: probit)
--case-control              病例-对照设计
--prevalence <float>        群体患病率 (用于liability尺度转换)
--liability-scale           输出liability尺度方差组分
--pql-maxit <n>             PQL最大迭代 (默认: 20)
--pql-tol <float>           PQL收敛容忍度 (默认: 1e-6)
```

### 3.4 实现路线

```
v1.0: PQL + 加权LMM + score test GWAS
v1.1: SPA校正完整接入
v2.0: Laplace近似 + AGHQ小样本精确似然
v3.0: 阈值-线性多性状联合
```

### 3.5 关键代码模块

| 模块 | 位置 | 说明 |
|---|---|---|
| ThresholdModel | `models/threshold_model.h/.cpp` | 顶层面板 |
| PQL solver | `models/threshold_model.cpp` | IRLS迭代 |
| LmmEngine | `lmm_engine.h/.cpp` | 加权LMM扫描复用 |

---

## 4. CountModel — 计数性状模型

### 4.1 统计模型

```
y_i ~ Poisson(λ_i)
log(λ_i) = Xb + Za

a ~ N(0, K · σ²_g)
```

或负二项扩展：

```
y_i ~ NB(μ_i, φ)
Var(y_i) = μ_i + φ · μ²_i
log(μ_i) = Xb + Za
```

### 4.2 PQL 算法 (Poisson)

```
1. η = Xb + Za
2. μ = exp(η)
3. z = η + (y - μ) / μ
4. w = μ
5. 解加权LMM (复用EigenVce)
6. 收敛检查
```

### 4.3 CLI 接口

```
--count-model               启用计数模型
--family poisson            分布族 (默认: poisson)
--family negbin             负二项
--offset <column>           偏移量列 (如log(有效时间))
--overdispersion            估计过度离散参数
```

### 4.4 实现路线

```
v1.0: Poisson PQL + 加权LMM
v1.1: 负二项 + 过度离散估计
v2.0: 零膨胀模型
```

---

## 5. SurvivalModel — 生存分析模型

### 5.1 统计模型

Cox比例风险 mixed model：

```
h(t|X, Z, g) = h₀(t) · exp(Xb + gβ + Zu)

u ~ N(0, Kσ²_u)
```

或 Weibull 动物模型：

```
h(t) = λ · p · (λt)^(p-1) · exp(Xb + Za)

a ~ N(0, Kσ²_a)
```

### 5.2 算法路线

MVP：Cox score test with kinship correction

```
1. 拟合Cox null模型 → 得到鞅残差 M = δ - H₀(t)exp(Xb)
2. 用鞅残差作为"表型"拟合LMM:
   M = Xb + Za + e, a ~ N(0, Kσ²_g)
3. 对每个SNP: score test of β=0
```

### 5.3 CLI 接口

```
--survival-model            启用生存模型
--time-col <n/name>         时间列
--event-col <n/name>        事件指示列 (0/1)
--model cox                 选用Cox (默认)
--model weibull             选用Weibull
--cox-efron                 使用Efron近似处理结
```

### 5.4 实现路线

```
v1.0: Cox score test + 鞅残差LMM
v1.1: frailty Cox mixed model
v2.0: Weibull animal model
```

---

## 6. GxEModel — 反应范/基因型-环境互作模型

### 6.1 统计模型

反应范模型 (扩展RRM)：

```
y_ij = Xb + Σₖ a_0k·Φ_k(time) + Σₖ a_1k·Φ_k(time)·env_j + e

[a_0, a_1] ~ N(0, K ⊗ G_env)
```

其中 `env_j` 是环境协变量(温度、群体、场址等)。

### 6.2 CLI 接口

```
--gxe-model                 启用GxE模型
--env-covar <col/name>      环境协变量列
--reaction-norm-order <n>   反应范阶数 (默认: 1)
```

### 6.3 实现路线

```
v1.0: 在RrmModel基础上增加env协变量交互
v1.1: 多环境方差组分估计
v2.0: 环境协变量连续GxE
```

---

## 开发顺序与里程碑

### Phase 1 — 核心模型 (当前)

```
1. MaternalModel      [开发中]
2. SingleStepModel    [开发中]
3. ThresholdModel     [开发中]
```

### Phase 2 — 扩展模型

```
4. CountModel
5. SurvivalModel
```

### Phase 3 — 高级生产模型

```
6. GxEModel
7. 多性状阈值-线性模型
8. 贝叶斯GWAS (BayesR/BayesCπ)
```

---

## 文件结构

```
CosmicScanner/docs/models/
├── New_Models_Design.md          ← 本文档

CosmicScanner/src/models/
├── univariate_model.h/.cpp       (现有)
├── repeatability_model.h/.cpp    (现有)
├── rrm_model.h/.cpp              (现有)
├── multivariate_model.h/.cpp     (现有)
├── maternal_model.h/.cpp         ← 新增
├── single_step_model.h/.cpp      ← 新增
├── threshold_model.h/.cpp        ← 新增
├── count_model.h/.cpp            ← 新增
├── survival_model.h/.cpp         ← 新增
└── gxe_model.h/.cpp              ← 新增
```
