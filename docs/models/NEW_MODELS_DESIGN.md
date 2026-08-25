# CosmicScanner 新模型设计文档

## GEMMA compatibility layer

GEMMA-compatible standalone routes now share one dispatch rule: only one of `--make-grm/-gk`、`--gs/-gs`、`--gq/-gq`、`--gemma-eigen/-eigen`、`--vc/-vc`、`--ci/-ci`、`--calccor/-calccor`、`--predict-effect/-epm`、`--bslmm/-bslmm` can run in a single command. Common GEMMA short aliases are accepted for validation convenience: `-o` for output prefix, `-p` for phenotype, `-n` for phenotype column, `-lm` for ordinary linear model, and `-lmm` for LMM mode. GEMMA `-lmm 5` maps to the null-model-only VC route (`--vc 2` style REML output), and `-lmm 9` is accepted as a GEMMA LRT-style alias over the existing LMM scan. This keeps CS long options available while allowing GEMMA-style test commands to map to the same implementation paths.

## GEMMA-style 单性状 LM

`--lm` 实现 GEMMA `-lm` 的普通线性模型 GWAS 路线，不估计随机效应和 LMM null model。`--lm-mode 1/2/3/4` 分别对应 Wald、LRT、Score 和 all；输出文件为 `.lm.assoc.txt`。实现使用和 GEMMA 源码一致的协变量投影统计量：`x'Px`、`x'Py`、`y'Py` 和 `y'P_x y`，Wald/Score 使用 F 分布尾概率，LRT 使用 1 自由度卡方尾概率。该路线支持 PLINK/PGEN/BGEN 输入，以及 `--maf`、`--extract`、`--exclude`、`--chr`、`--from-bp`、`--to-bp` 等扫描过滤。

## GEMMA-style 单性状 LMM mode

`--lmm-gemma-mode 1/2/3/4` 对应 GEMMA `-lmm` 的 Wald、LRT、Score 和 all 输出。默认路线复用单性状 LMM null model 方差结构，逐 SNP 在 GLS/LMM 旋转空间计算 Score/Wald 统计量，并用 alternative 固定方差 GLS 残差二次型下降计算 LRT；mode 4 输出 `P_WALD`、`P_LRT`、`P_SCORE`，mode 1/2/3 将主 `P` 列分别切换为 Wald、LRT 或 Score。

`--lmm-gemma-refit` 启用更慢的逐 SNP alternative profile refit：Wald 使用 alternative REML profile 的 SNP beta/SE，LRT 使用 null/alternative ML profile log-likelihood，Score 保持 null-model score test。该路线更接近 GEMMA `-lmm 4` 的 Wald/LRT/Score 数值，但计算量高于默认 null-structure scan。

## GEMMA-style prediction

`--predict-effect <file>` 提供 GEMMA `-predict` 风格的基因组预测路线：从 SNP effect 表读取 `SNP/rs/id/variant_id` 与 `BETA/beta/alpha/effect/EB_BETA` 列，按 genotype dosage 与 effect 加权求和，输出 `FID IID PRED N_VARIANTS_USED`。可通过 `--predict-snp-col` 和 `--predict-effect-col` 指定非标准列名，通过 `--predict-score <file>` 指定输出文件。该路线服务于 GEMMA/CS association、后续 BSLMM/ridge effect 文件和外部 effect 表的预测复用；它不是 BSLMM/BVSR MCMC 采样器本体。

## GEMMA-style eigen utility

`--gemma-eigen` 复用 GRM 构建器并对 dense GRM 做 self-adjoint eigen decomposition，写出 `<out>.eigenD.txt` 和 `<out>.eigenU.txt`。该路线覆盖 GEMMA `-eigen` 的常见复用场景：预先保存 GRM 特征值和特征向量，供 LMM/BSLMM 或外部脚本读取。其 GRM 标准化由 `--make-grm-alg`、`--maf`、`--hwe`、`--extract`、`--exclude` 等 CS 过滤和算法参数控制。

## GEMMA-style calccor utility

`--calccor`/`--gemma-calccor` 实现 GEMMA `-calccor` 风格的 SNP correlation window 输出，文件为 `<out>.cor.txt`，列为 `chr rs ps n_mis n_obs allele1 allele0 af window_size var cor`。`--windowbp` 默认 1000000，`0` 表示不设 bp 距离上限；`--windowcm` 默认 0，当前仅 PLINK BIM 的 cM 坐标有意义；`--windowns`/`-windowns` 默认 0，表示不设右侧 SNP 个数上限。输出中的 `cor` 为当前 SNP 与右侧同染色体窗口内 SNP 的逗号分隔相关系数；没有右侧邻居时写 `NA`。该路线支持 PLINK/PGEN/BGEN 输入，以及 `--maf`、`--hwe`、`--extract`、`--exclude`、`--chr`、`--from-bp`、`--to-bp` 过滤。与 GEMMA 源码一致，窗口扫描依赖输入变异顺序；实际 LD/相关分析建议输入按 `chr, bp` 排序。

## GEMMA-style variance component utility

`--vc 1/2/3` 提供 GEMMA `-vc` 风格的单性状 variance component 独立估计入口，不执行 SNP association scan。`--vc 1` 使用 Haseman-Elston/MinQue 风格矩估计；`--vc 2` 使用 REML；`--vc 3` 使用 AI/IDUL-style REML 求解。输出 `<out>.vc.txt` 包含 `Vg`、`Ve`、PVE 和占位 SE 列，`<out>.vc.summary.txt` 记录 method、样本数、协变量数、SNP 数、`sigma2_g`、`sigma2_e`、`pve` 和 `delta`。该路线复用 CS 的 GRM builder、phenotype/covariate 对齐、缺失样本过滤和 `--maf/--hwe/--extract/--exclude` 过滤。

## GEMMA-style summary-stat reference utilities

`--gs`/`-gs` 生成 GEMMA summary-stat VC 使用的 reference S matrix，输出 `<out>.S.txt`、`<out>.size.txt` 和 `<out>.snps.txt`。单 component 模式对 centered/scaled K 按 GEMMA `compAKtoS` 形式计算，并用 leave-one-sample jackknife 近似输出 S variance 行；加入 `--cat/-cat`、`--mcat/-mcat`、`--wcat/-wcat` 或 `--wsnp/-wsnp` 时，按 category 构建 component GRM，输出多 component S/Svar 矩阵和 component size。`--gq`/`-gq --beta <file>` 从 summary association 表生成 `<out>.Vq.txt`、`<out>.q.txt` 和 `<out>.size.txt`；beta 文件需包含 `N` 或 `N_MIS+N_OBS`，并包含 `Z`、或 `BETA+SE`、或 `CHISQ`、或 `P`。多组件输出使用 GEMMA 风格矩阵/向量文本。`--study/-study <prefix> --ref/-ref <prefix> --vc/-vc 1` 或 `--mstudy/-mstudy <list> --mref/-mref <list>` 读取这些中间文件输出 `<out>.summary.vc.txt` 和 metadata；`--ci/-ci` 输出 `<out>.summary.ci.txt` 的 PVE 95% CI，单组件、多组件和多前缀合并均支持。也支持 GEMMA 直接 summary 输入组合：`--vc/-vc 1 --beta/-beta <assoc> --cor/-cor <ref.cor.txt>`，按 GEMMA `ReadFile_beta`/`ReadFile_cor` 的权重、LD pair 和 `CalcVCss` 形式直接输出 summary VC；加入 `--cat/--mcat/--wcat/--wsnp` 时输出多 component PVE、SE、sigma2 和 enrichment。

## GEMMA-style multi-locus ridge/BSLMM outputs

`--bslmm <mode>` 提供 GEMMA `-bslmm` 风格的 multi-locus 输出。`--bslmm 2` 实现 standard ridge regression/GBLUP no-MCMC 路线：对过滤后的 SNP dosage 做均值填补、中心化和标准化，使用 `--bslmm-h2` 转换为 ridge penalty，输出 `.param.txt`、`.bv.txt` 和 `.hyp.txt`。`.param.txt` 使用 GEMMA 风格列 `chr rs ps n_miss alpha beta gamma`，`.bv.txt` 输出每个个体的 breeding value / genomic score，`.hyp.txt` 输出 `h pve rho pge pi n_gamma`。

`--bslmm 1` 和 `--bslmm 3` 当前实现为 deterministic sparse+ridge approximation：先做全量 ridge shrinkage，再按 ridge effect 绝对值选 `--bslmm-top`/`-smax` 个 sparse gamma 位点并做稀疏 ridge refit；mode 3 要求 0/1 表型并把 `.bv.txt` 输出转换到概率尺度。`--bslmm 4` 支持 GEMMA DAP 风格 BF/HYP 生成，输出 `<out>.bf.txt` 和 `<out>.hyp.txt`；`--bslmm 5 --bf <file> --hyp <file>` 支持 DAP fine-mapping 汇总，输出 `<out>.gamma.txt`、`<out>.coef.txt` 和 `<out>.hyp.txt`。DAP grid 支持 GEMMA 风格 `-hmin/-hmax/-rmin/-rmax/-pmin/-pmax/-gmean/-smin/-smax`，其中 `-pmin/-pmax` 与 GEMMA 一致表示 `log10(pi)`；并提供 `--bslmm-hgrid/--bslmm-rho-grid` 控制网格密度。`-w/-s/-rpace/-wpace/-hscale/-rscale/-pscale/-mh/-seed` 可解析以兼容 GEMMA 命令，但 deterministic 近似不执行采样循环。mode 4/5 使用 deterministic Wakefield-style ABF / posterior-weight summary 近似，覆盖 GEMMA 文件接口和路线入口，但不是 GEMMA 源码中的 full DAP EM annotation engine；mode 1/3 也不是 full MCMC sampler。

## 多性状 GWAS Methods 草稿（已完成模块）

CosmicScanner 的多性状 GWAS 模块实现了面向多个相关性状的 multivariate linear mixed model。对 `n` 个个体和 `d` 个性状，模型写作 `Y = XB + g + e`，其中 `g ~ MN(0, K, Vg)`，`e ~ MN(0, I, Ve)`。`K` 为基因组关系矩阵，`Vg` 和 `Ve` 分别为遗传和残差性状协方差矩阵。程序对 `K` 做特征分解，将表型和协变量旋转到特征空间，并通过 EM 或 AI-REML 估计 null model 方差组分。

快速 `emmv` 路径复用 null model 方差结构进行 genome-wide SNP scan。每个位点对 genotype dosage 做样本子集化和缺失均值填补，在旋转混合模型空间估计跨性状 SNP 效应，并使用多自由度 Wald statistic 检验所有性状 SNP 效应同时为零的联合原假设。`--mv-compatible emmvgwas`/`--emmvgwas` 提供与 EMmvGWAS R 包一致的验证路径，包括 `Z = dosage - 1` 编码、trace-normalized GRM、`lambda = G R^{-1}` profiling、uncentered genotype scan 以及官方公式中的 general matrix inverse。`--mv-emmvgwas-output`/`--emmvgwas-output` 输出 EMmvGWAS 风格 `G.csv`、`R.csv`、`beta.csv` 和 association CSV；`--emmvgwas-vc-only` 覆盖 `vc_estimation` 路线，只估计并输出 G/R/beta 后停止；`--mv-emmvgwas-G/--mv-emmvgwas-R` 或 `--emmvgwas-G/--emmvgwas-R` 支持读入外部方差组分执行 `semi_gwas` 等价扫描。

`gemma`、`gemma-full` 和 `gemma-exact` 路径实现 GEMMA-style multivariate mixed-model validation。`--mv-gemma-mode 1/2/3/4` 对应 GEMMA mvLMM 的 Wald、LRT、Score 和 all；`gemma-full` 输出 Wald、LRT 和 Score p 值，并按 mode 控制主 `P` 列。`gemma-exact` 支持 GEMMA EM/NR 控制参数和逐 SNP full alternative refit；`--mv-crt` 启用 GEMMA-style Edgeworth/PCRT 小样本校正，输出 `P_WALD_CRT`、`P_LRT_CRT`、`P_SCORE_CRT`，并按 mode 将主 `P` 列切换为校正后 p 值。该实现为 CosmicScanner/Eigen 原生矩阵实现，没有直接链接 GEMMA 的 GSL 代码。

`gemma-gxe` 和 `gemma-gxe-full` 路径覆盖 GEMMA 多性状 GxE 思路：null model 包含协变量和环境主效应，每个 SNP 条件在 SNP 主效应下检验 `SNP x environment` 跨性状交互；full 路线输出 GxE Wald、LRT 和 Score p 值。该路线支持 PLINK/PGEN/BGEN 输入，当前要求 `--no-loco`。

Adaptive hybrid 路径用于在接近 GEMMA-full 精度的同时降低计算量。Stage 1 使用 EMmvGWAS-style scan 做全基因组筛选；启用 LOCO 时按染色体构建 LOCO GRM，估计 chromosome-specific multivariate null model，并只扫描当前染色体。Stage 1.5 可选对 Stage 1 top SNP 做 Score-only audit，产生 null-based `P_SCORE`。Stage 2 对三类候选的并集执行 GEMMA-full refinement：Stage 1 p-threshold、Stage 1 topK、以及 `abs(-log10(P_SCORE) - -log10(P_WALD))` 最大的 Score/Wald discrepancy SNP。LOCO-hybrid Stage 2 按候选 SNP 所在染色体分组，复用对应 LOCO GRM 做 per-chromosome GEMMA-full refinement，并 append 到统一 refined association file。

实现上，快速多性状 scan、GEMMA-style refit 和 GEMMA-style GxE 对 PLINK/PGEN/BGEN 使用 block-wise genotype streaming；GRM 构建和 SNP scan 共享 `extract`、`exclude`、染色体和 bp range filters，避免 GRM SNP set 与 scan SNP set 不一致。Hybrid refinement 当前限定 PLINK 输入，因为它需要 repeated selective SNP access 和 per-chromosome LOCO refinement。

验证包括 EMmvGWAS R 官方 RIL example、GEMMA 0.98.5 multivariate LMM 对标、多性状 LOCO/no-LOCO validation matrix、adaptive hybrid smoke，以及系统 benchmark。默认系统 benchmark 比较 EMmvGWAS R、GEMMA、CS-fast、CS-hybrid 和 CS-LOCO-hybrid，并输出 runtime、memory、p-value correlation、top overlap、causal recovery 以及图表。


## 前言

本文档记录 CosmicScanner GWAS 模块的新模型开发设计。现有 7 种模型/模式
（单性状 LMM、重复力、RRM、多性状、罕见变异、Biobank、SPA）
已覆盖生猪育种约 60-70% 的应用场景。

为覆盖 95% 以上生产需求，新增以下 6 个模型：

| # | 模型 | 优先级 | 主要应用 |
|---|---|---|---|
| 1 | MaternalModel | P0 | 出生重/断奶重等母体效应性状 |
| 2 | SingleStepModel | P0 | 系谱+基因组联合 ssGWAS |
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

## 2. SingleStepModel — ssGWAS

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
--single-step               启用 ssGWAS
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

## 7. CommonEnvironmentModel — 公共环境模型（已完成，实验性）

### 7.1 统计模型

同一受试对象存在多条记录、且记录来自不同环境（组）时，按组引入随机截距：

```
y_r = X_r b + x_jr·beta_j + (Z_a a)_r + (Z_c c)_r + e_r

a ~ N(0, K·sigma2_a)      # 跨受试对象的基因组 GRM 遗传效应
c ~ N(0, I·sigma2_c)      # 按组共享的环境随机截距（组间独立）
e ~ N(0, I·sigma2_e)
```

- `Z_a`：记录→受试对象（subject）映射；`K` 为 subjects×subjects 的 GRM（VanRaden，默认算法 0）。
- `Z_c`：记录→组（group）映射；同组记录共享 `sigma2_c` 相关，组间独立。
- 方差分量 `(sigma2_a, sigma2_c, sigma2_e)` 通过 REML 估计（可 `--reml-vg/-vc/-ve` 固定）。

### 7.2 数据结构映射

| 概念 | 实现 |
|------|------|
| 记录 | `GenRecord`（按 `FID IID` 对齐表型长格式记录） |
| 受试对象 | `obs_to_sample_`：记录→subject 秩 |
| 组 | `obs_to_group_`：记录→组 level（subject 级或 record 级分组文件） |
| GRM | `G_subject_`（subjects×subjects），扫描时按保留 subject 排序 `subsetByRank` |
| 固定效应 | `FixedDesignG`：intercept + 记录级协变量（`--covar/--qcovar`） |

### 7.3 CLI 接口

```
--common-env-model            启用公共环境模型
--common-env-file <f>         分组文件（FID IID GROUP，或 record 级 FID IID [TIME] GROUP）
--common-env-record-level     分组文件按记录解析
--reml-vg/-vc/-ve <v>         固定方差分量（跳过 VCE）
--no-loco                     关闭逐染色体 LOCO
```

### 7.4 扫描算法（记录层分数检验）

对每个 SNP 块（block_size=256）构造 `G_block`（n_records × bsz，均值中心化），
基于已估 null 模型做 GLS score 检验：

```
scores   = G'·P·y            P = V^-1 - V^-1 X (X'V^-1 X)^-1 X'V^-1
gVg      = colsum(G ⊙ V^-1 G)
xvX      = G'·V^-1 X
term2    = rowwise( (xvX · (X'V^-1X)^-1) ⊙ xvX )
var_score = gVg - term2
beta = score / var_score,  se = sqrt(1/var_score),  p = erfc(|beta/se|/√2)
```

- `V = vg·Z_a·K·Z_a' + vc·Z_c·Z_c' + ve·I`，由 `multi_random_effect_lmm`（LDLT + jitter 因子化）求解。
- MAF/HWE QC（`--min-maf/--min-hwe`）与频率报告复用 `calculateHWE`；`n_obs` 为有效记录数。
- LOCO：`computeSplit` 得到全局与分染色体标准化基因型交叉积，`G_loco = G*(S_global/S_loco) - raw_chr/S_loco`（`raw_chr` 为未归一化 `Z_chr·Z_chr'`，公式与 `gwas_scanner.cpp` 的 PLINK-style LOCO 一致），并按染色体重估 VCE。

### 7.5 关键代码模块

| 文件 | 职责 |
|------|------|
| `src/models/common_environment_model.{h,cpp}` | 模型入口、分组解析、VCE、扫描、输出 |
| `src/models/multi_random_effect_lmm.{h,cpp}` | 通用 multi-random-effect LMM 层（prepare/Vinv/Py/scoreTest） |
| `src/gwas_app.cpp` / `src/main.cpp` | 互斥检查、CLI 参数、dispatch |
| `tests/run_test_common_environment.sh` | 模拟验证脚本 |

### 7.6 验证结果

模拟：N=120、M=160（chr1=100, chr2=60）、12 组×10 个体、360 条记录（2+i%3/受试对象）、真值 `sigma2a/c/e = 0.45/0.25/0.30`、因果 SNP rs8（β=0.7）。

| 检查项 | 结果 |
|--------|------|
| VCE 回收（global） | 0.4877 / 0.3475 / 0.3227（误差 0.0377 / 0.0975 / 0.0227，全部在容差内） |
| 因果 SNP 排名 | top hit，p=6.5e-07（global）、5.7e-10（LOCO） |
| 关联行数 | 160/160（global/LOCO/record-level/fixed-VC 全部） |
| 固定 VC 路径 | 精确使用 `--reml-vg/-vc/-ve` 指定值 |
| 分组语义 | subject 级按 GROUP LABEL 建 level（120 个体→12 组），record 级按行对齐 |
| 回归 | `run_test_repeatability.sh`、`run_rrm_dynamic_simulation.sh` 通过 |

---

## 8. AnimalModel — 系谱动物模型（已完成，与 ssGWAS 共享 A/A22）

### 8.1 统计模型

Simple Animal Model 用系谱亲缘矩阵 A 替代基因组 GRM 描述多基因背景：

```
y = Xb + x_j·beta_j + Z·a + e

a ~ N(0, A·sigma2_a)      # 系谱加性遗传效应
e ~ N(0, I·sigma2_e)
```

方差分量 `(sigma2_a, sigma2_e)` 通过 REML 估计（复用 `EigenVce` 单 GRM 路径）。

### 8.2 与 ssGWAS 的代码共享

- 新增共享类 `PedigreeRelationship`（`src/models/pedigree_relationship.{h,cpp}`）：系谱解析（ID SIRE DAM，`NA`/`.`/`0` 表示缺失）、去重/自亲/环检测、重亲序（parents-before-offspring，乱序不变性由 min-heap 拓扑序保证）。
- `buildAInverse()`：Henderson 方法构建稀疏 A⁻¹（`invertA_henderson_fixed`），置换回原系谱顺序。
- `buildA22(indices)`：对 A⁻¹ 的 Schur complement `A22 = A22 − A21 A11⁻¹ A12` 做 block-inverse 求解（`SimplicialLDLT`），**不物化稠密 A**。
- `invertSpdWithRidge()` 迁入 cosmic 命名空间（LDLT + ridge 1e-8×10 递增 6 次 + FullPivLU 兜底）。
- ssGWAS 已重构为复用上述 helper（删除本地 PedRecord/ped_id_to_idx_/本地 invertSpdWithRidge）；单一实现保证两个模型 A/A22 数值完全一致。

### 8.3 CLI 接口

```
--animal-model            启用系谱动物模型（隐式 --gwas）
--ped <f>                 系谱文件：ID SIRE DAM（必需）
--write-ped-kernel        导出分析子集精确 A22 到 *.animal_model.k.txt + .k.id
--no-loco                 关闭 LOCO（多染色体时消除 LOCO-not-applicable 提示）
```

ID 匹配：先 `FID:IID`，回退裸 `IID`；分析样本必须都在系谱中。

### 8.4 扫描算法

- A22 按 keep 样本（去缺失表型后）构建；`EigenVce::prepare(y, X, A22) → runNullModel()` 估计 null 方差分量。
- `LmmEngine::prepareFromVce` 复用 null 统计量，单次全扫描输出标准 `*.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`）。
- Null summary：`*.animal_model.null.txt`（`sigma2_a/sigma2_e/h2` + 系谱解析/Ainv/A22/VCE/scan 5 项 timing）。
- LOCO 不适用：只有一个关系矩阵；多染色体且未 `--no-loco` 时打印提示（仪器输出警告，非错误）。

### 8.5 关键代码模块

| 文件 | 职责 |
|------|------|
| `src/models/pedigree_relationship.{h,cpp}` | 共享系谱解析、A⁻¹、A22（Schur complement）、`invertSpdWithRidge`（新增） |
| `src/models/animal_model.{h,cpp}` | 模型入口：prepare、VCE、扫描、null summary、kernel 导出（新增） |
| `src/models/single_step_model.{h,cpp}` | 重构为复用 `PedigreeRelationship`（删除本地重复实现） |
| `libCosmicCore/src/eigen_vce.{h,cpp}` | 新增 `converged()`/`lastIterations()` getter（Brent 路径报 iterations=-1） |
| `src/config.h` / `src/main.cpp` / `src/gwas_app.cpp` | `--animal-model`/`--write-ped-kernel`、校验、dispatch |
| `tests/run_animal_model_simulation.py`、`tests/run_test_animal_model.sh`、`tests/run_animal_model_gemma_validation.sh` | 模拟器与验证脚本（新增） |

### 8.6 验证结果

| 检查项 | 结果 |
|--------|------|
| 默认模拟 causal 检出 | rs1 进入 top10，h2 VCE 回收在容差内 |
| ssGWAS tau=0 等价 | `--h2-tau 0 --h2-omega 0 --no-loco` 下 β/SE/-log10P Pearson = 1.0（maxabs 1e-6/0/1e-7），null VCE 完全一致 |
| 系谱乱序不变性 | 打乱系谱行序后 assoc 与 null VCE 一致到 <1e-12 |
| σ²a=0 | h2 < 0.1 |
| 高遗传力 | h2 > 0.5（输入 0.8 场景） |
| 不完整系谱 | 缺失单亲/双亲仍能构建，causal top10 |
| 无祖先系谱 | typed 全独立（隐式 I 背景），600 行有限值、converged=yes |
| A⁻¹ / A22 数值 | 对照 Python 暴力 inv(A) maxabs=0.0；Schur A22 vs 暴力 A22 maxabs≈5e-8 |
| GEMMA -k 外部基准 | GEMMA 0.98.5 读 `--write-ped-kernel` 导出的 A22 跑 `-lmm 4`：β Pearson≈0.99996、SE≈0.9992、-log10P≈0.9984；top5 完全重合、因果 SNP 排名 1/1；GEMMA null VCE 与 CS 精确一致 |

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
├── NEW_MODELS_DESIGN.md          ← 本文档

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
├── gxe_model.h/.cpp              ← 新增
├── multi_random_effect_lmm.h/.cpp ← 新增（通用多随机效应 LMM 层）
└── common_environment_model.h/.cpp ← 新增（公共环境模型）
```
