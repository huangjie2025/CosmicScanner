# CosmicScanner：面向复杂性状结构的统一、可扩展全基因组关联分析框架

**Jie Huang**

浙江大学动物科学学院，杭州 310058，中国

通讯作者：huangjie2025@zju.edu.cn

---

## 亮点

- 单一 C++ 引擎统一支持 11 类 GWAS 模型，包括标准 LMM、多性状 LMM、重复力模型、随机回归模型、母体效应、ssGWAS、阈值/序数模型、生存模型、计数/二项 PSR-GWAS 和 GxE 反应规范模型。
- 针对重复测定数据，Woodbury 主体水平降维在保持重复力混合模型协方差结构精确性的同时，将主要线性代数从记录水平降至个体水平，相比 RepeatABEL 实现 10-20 倍加速和约 10 倍内存降低。
- 针对非平衡纵向数据，BLUP-based Kronecker RRM 先估计个体 Legendre 系数，再在系数空间扫描 SNP；在 GMAT mouse_long benchmark 上实现 213 倍加速和 93 倍内存降低。
- 与 GEMMA、EMmvGWAS、RepeatABEL 和 GMA 的主结果达到近乎完全一致（Pearson r > 0.99）；与 GMAT 非平衡纵向 RRM benchmark 达到 close concordance（Pearson r = 0.951，Spearman rho = 0.940，top SNP 一致）。
- 六类扩展模型（母体效应、threshold、GxE、生存、Poisson PSR、binomial PSR）通过 type-I error calibration（lambda_GC 约为 1.00）；母体效应模型进一步通过模拟真值恢复验证。
- v1.1.0 release candidate 默认 CTest 包含 10 个测试，覆盖兼容性 smoke、外部 reference 路线、sparse REML、multi-VC null model、SPA、REGENIE 和输入适配层。
- 原生支持 PLINK 1.9、PLINK 2.0 PGEN 和 BGEN，减少格式转换带来的存储、时间和误差成本。

---

## 摘要

现代遗传研究中的表型数据已经远超标准单性状连续性状：重复测定、纵向轨迹、非高斯结局、多性状相关结构、母体效应和环境互作在农业、医学和进化遗传学中都十分常见。然而，现有 GWAS 软件生态高度碎片化，不同模型类别通常依赖不同工具、不同输入格式、不同 GRM 构建方式和不同输出 schema。重复测定和纵向模型尤其受限，因为许多工具需要构建记录或观测水平的稠密协方差矩阵，其 O(n^3) 分解复杂度使样本规模很快触及计算瓶颈。

本文介绍 CosmicScanner，一个统一的 C++ GWAS 框架，在同一引擎中实现 11 类 GWAS 模型，共享基因型 I/O、GRM 构建、null model 缓存和 block-wise SNP 扫描。针对最难扩展的模型类别，CosmicScanner 引入两项核心计算设计：第一，重复力 GWAS 使用精确 Woodbury 主体水平降维，通过 incidence matrix Z 重写重复测定协方差，把稠密计算限制在个体水平；第二，纵向 RRM 使用 BLUP-based Kronecker multivariate LMM，先估计 Legendre 随机回归系数，再在系数空间扫描 SNP，避免构造观测水平协方差矩阵。

CosmicScanner 原生读取 PLINK 1.9、PLINK 2.0 和 BGEN，并支持 memory-mapped、多线程 I/O。与 GEMMA、EMmvGWAS、RepeatABEL、GMA 和 GMAT 的验证显示：单性状、多性状、重复力和 balanced RRM 结果达到近乎完全一致（Pearson r > 0.99），非平衡 longitudinal RRM 与 GMAT 达到 close concordance（r = 0.951）。重复力 GWAS 相比 RepeatABEL 加速 10-20 倍、内存降低约 10 倍；纵向 RRM 相比 observation-level 路线加速 213 倍、内存降低 93 倍。六类扩展模型通过 type-I error calibration，母体效应模型可恢复模拟残差方差（误差 < 6%）。Release build 通过 10-test CTest suite。总体而言，CosmicScanner 为复杂表型结构 GWAS 提供了一个统一、可验证、可扩展的分析框架。

**关键词：** GWAS；线性混合模型；纵向性状；重复测定；多性状分析；随机回归；Woodbury identity；基因组关系矩阵

---

## 重要性

复杂表型数据迫使研究者在多个专用工具之间切换：重复记录用 RepeatABEL，纵向模型用 GMA 或 GMAT，多性状用 GEMMA 或 EMmvGWAS，序数或非高斯性状又需要其他模型工具。这种流程不仅繁琐，也容易在 GRM、协变量、缺失值处理和输出解释上产生不可见的不一致。CosmicScanner 的贡献在于：用一个统一引擎覆盖 11 类 GWAS 模型，同时对重复测定和纵向模型引入模型感知的降维计算。重复力模型中的 Woodbury 变换保持精确协方差结构；RRM 中的 BLUP coefficient-space 策略在较小统计偏差下换取数量级级别的计算收益。该框架使复杂性状 GWAS 从“多工具拼接”转向“统一、可审计、可复现”的单一流程。

---

## 1. 引言

自 2005 年以来，GWAS 已成为解析复杂性状遗传基础的核心方法。标准单性状 LMM 在控制群体结构和亲缘关系方面非常有效，并通过 GEMMA、GCTA、BOLT-LMM、REGENIE 等工具得到广泛应用。然而，真实研究中的表型结构往往更复杂：畜禽和植物育种常见重复测定和不同环境下的纵向记录；生物医学 cohort 同时包含连续、二分类、计数、序数和生存结局；多性状联合分析利用性状间遗传相关提高检出力。

现有工具虽然各自在特定模型上成熟，但软件生态存在三个结构性缺口。

第一，模型覆盖碎片化。没有一个工具能够完整覆盖现代研究中常见的模型谱系。研究者通常需要组合 GEMMA、RepeatABEL、POLMM、GMA、GMAT 或其他工具，每个工具有自己的输入格式、GRM 计算方式、方差组分估计策略和输出格式。这会使同一研究中的多个模型结果难以严格比较。

第二，重复测定和纵向 GWAS 难以扩展。RepeatABEL 构造记录水平协方差矩阵 V 并进行稠密特征分解，复杂度随记录数立方增长。若 1,000 个个体每个有 10 条记录，则需要处理 10,000 x 10,000 的矩阵。Observation-level RRM 也面临同样问题，协方差维度等于总观测数。

第三，现代基因型格式支持不完整。PLINK 2.0 PGEN 和 BGEN 已成为大规模队列数据的常用格式，但许多模型工具仍只接受 PLINK 1.9 BED，导致额外格式转换和存储负担。

CosmicScanner 的目标是同时解决这三个问题：在统一 C++ 引擎中实现 11 类模型；共享 genotype I/O、GRM、null model 和 SNP scanning 基础设施；对重复测定和纵向模型使用模型结构本身进行降维，而不是简单牺牲统计模型。

本文主要贡献包括：

1. **统一模型覆盖。** 支持单性状 LMM、多性状 LMM、重复力模型、RRM、母体效应、ssGWAS、threshold/probit、count/binomial PSR、ordinal PSR、生存 Cox 和 GxE reaction-norm。
2. **重复力 GWAS 的 Woodbury 主体水平降维。** 保持 record-level mixed model 的精确协方差结构，但将主要矩阵运算从记录水平降到个体水平。
3. **非平衡纵向 GWAS 的 BLUP-based Kronecker RRM。** 先估计个体随机回归系数，再在系数空间进行 SNP 联合检验。
4. **系统验证。** 与 GEMMA、EMmvGWAS、RepeatABEL、GMA 达到 r > 0.99；与 GMAT 非平衡 RRM 达到 r = 0.951、Spearman rho = 0.940 且 top SNP 一致；六类扩展模型通过 type-I error calibration。

---

## 2. 材料与方法

### 2.1 统一输入和 GRM 构建

CosmicScanner 原生支持 PLINK 1.9 BED/BIM/FAM、PLINK 2.0 PGEN/PVAR/PSAM 和 BGEN。三类格式通过共同 reader interface 向 GWAS scanner 提供 genotype dosage vector。对 PGEN 和 BGEN，程序使用 memory-mapped I/O 和多线程解压，按 SNP block 加载数据，避免全基因型矩阵常驻内存。

GRM 使用 VanRaden-style centered genotype matrix 构建：

```text
G = (M - P)(M - P)' / [2 sum_i p_i(1-p_i)]
```

其中 M 为 0/1/2 基因型计数矩阵，P 为对应等位基因频率的列均值矩阵。GRM 构建支持 PLINK、PGEN 和 BGEN 输入，也支持 chromosome-partitioned GRM 和 LOCO 分析。`--extract`、`--exclude`、`--chr`、`--from-bp`、`--to-bp` 等过滤条件同时用于 GRM 和 SNP scanning，以保证 null model 与扫描集合一致。

表型和协变量使用 PLINK-style 文件，同时为重复测定和纵向模型支持 long-format phenotype，其中每一行代表一个 observation record，通过 `FID IID TIME` 或相关字段索引。输入适配层可将 CSV/TSV/TXT 业务表转换为可审计的中间格式。

### 2.2 单性状 LMM GWAS

标准模型为：

```text
y = X beta + x_j alpha_j + g + e
g ~ N(0, sigma_g^2 G)
e ~ N(0, sigma_e^2 I)
```

CosmicScanner 对 GRM 做特征分解，将表型、协变量和 SNP dosage 旋转到特征向量空间，使协方差矩阵对角化，从而实现每 SNP O(n) 的 weighted least squares 检验。程序同时支持 fast-LMM approximation 和 LOCO，用于减少 proximal contamination。

### 2.3 IDUL/AIDUL null-model 加速

CosmicScanner 提供三类 null-model solver：基础 IDUL、带 monotone likelihood safeguard 的 IDUL+，以及加入 Anderson acceleration 的 AIDUL。默认仍保留 Brent solver 作为稳健后备；当迭代 solver 未收敛时，程序自动回退，保证输出有效。

### 2.4 多性状 LMM GWAS

多性状模型写作：

```text
Y = XB + x_j a_j' + U + E
U ~ MN(0, G, Sigma_g)
E ~ MN(0, I, Sigma_e)
```

CosmicScanner 实现五条多性状路线：EMmvGWAS-compatible、GEMMA-full、Hybrid、Empirical Bayes 和 GEMMA-compatible Wald/LRT/Score。EMmvGWAS-compatible 路线用于对齐 EMmvGWAS；GEMMA-full 进行 per-SNP variance-component refit；Hybrid 通过快速筛选和选择性精修平衡速度与精度。

### 2.5 重复力模型 GWAS

重复力模型针对一个个体多条记录：

```text
y = X beta + Z a + Z p + e
a ~ N(0, sigma_a^2 G)
p ~ N(0, sigma_p^2 I)
e ~ N(0, sigma_e^2 I)
```

其中 Z 是 record-to-subject incidence matrix。传统方法直接构造记录水平 V 并分解。CosmicScanner 利用 Woodbury identity 将求逆和 log determinant 重写到 subject-level core matrix，使稠密计算主要依赖 subject count，而非 record count。程序提供 exact、fast、STCG 和 mixed precision 四种策略，其中 exact Woodbury 是正式主结果路线。

### 2.6 随机回归模型 RRM GWAS

RRM 使用 Legendre polynomial 描述 SNP 效应随时间变化：

```text
alpha_j(t) = sum_k alpha_jk phi_k(t)
```

对非平衡纵向数据，CosmicScanner 首先对每个个体估计 Legendre coefficient BLUP，然后在 coefficient space 中进行 Kronecker covariance multivariate scan。该设计避免 observation-level covariance，同时保留对动态效应曲线的解释能力。

### 2.7 非高斯性状 GWAS

CosmicScanner 的 PSR-GWAS 框架支持 binary、count、binomial 和 ordinal 等非高斯结局。方法基于 penalized quasi-likelihood linearization 构造 working variate，并将 working residual variance 作为权重传入 LMM scanner。二分类路径结合 saddlepoint approximation 以改善 case-control imbalance 下的尾部 p 值校准。

### 2.8 母体效应、ssGWAS、GxE 和 rare-variant tests

母体效应模型包含 direct additive、maternal additive、direct-maternal covariance、maternal permanent environmental 和 residual variance components。ssGWAS 通过 H inverse 融合 pedigree A 和 genomic G 信息。GxE reaction-norm 模型将 SNP 效应建模为环境协变量 E 的线性函数，并进行 main effect 与 interaction effect 的 joint Wald test。Rare-variant analysis 支持 SKAT 和 burden tests。

### 2.9 实现

CosmicScanner 使用 C++17 和 Eigen3 实现，核心并行化依赖 OpenMP，GRM 和 SNP block scanning 的内层循环使用 AVX2/AVX-512 SIMD 加速。软件模块包括 genotype readers、GRM builder、LMM engine、model dispatch 和 shared core library。代码以 MIT License 发布，主仓库为 https://github.com/huangjie2025/CosmicScanner。

---

## 3. 结果

验证围绕四条主线展开：与 reference tools 的统计一致性、PLINK/PGEN/BGEN 格式一致性、重复测定和纵向模型的可扩展性、非高斯扩展模型的校准。整体结果显示，CosmicScanner 在外部 reference 覆盖的模型上达到 high concordance，并在重复力和纵向模型上获得一到两个数量级的运行时间和内存改善。

**图 1. CosmicScanner 架构和模型覆盖。** 图件见 `docs/figures/figure1_architecture.pdf`。该图展示 data loaders、GRM builder、LMM engine、model dispatch、11 类模型和端到端数据流。

### 3.1 单性状 LMM：准确性和格式一致性

与 GEMMA 的 benchmark 中，CosmicScanner exact LMM 在 effect size、standard error 和 p value 上高度一致。代表性结果如下：

| 指标 | 与 GEMMA 的一致性 |
|---|---:|
| beta correlation | 0.99987 |
| SE correlation | 0.99953 |
| -log10(P) correlation | 0.99929 |

PLINK BED、PGEN 和 BGEN 三类输入在同一数据集上产生一致的 beta、SE 和 p value，最大 p 值差异低于 1e-12，说明格式 reader 向 scanner 提供了等价 dosage interface。

### 3.2 多性状 GWAS：准确性和性能

EMmvGWAS-compatible 路线在官方 RIL example dataset 上与 EMmvGWAS 对齐。CosmicScanner 在 0.32-0.43 秒内完成，而 EMmvGWAS R 实现用时约 24.7 秒，约 60 倍加速。

GEMMA-full、Hybrid、EMMV 和 EB 路线在 N=240、M=1,200、2 traits 的系统 benchmark 中均能恢复主要 causal SNP。快速路线相比 GEMMA 实现 60-350 倍加速；LOCO-aware hybrid 在保持较高 accuracy concordance 的同时仍明显快于 GEMMA。

**图 5. 多性状 GWAS benchmark。** 图件见 `docs/figures/figure5_multivariate_benchmark.pdf`。该图展示 runtime、memory、accuracy concordance 和 causal SNP recovery。

### 3.3 重复力 GWAS：准确性和可扩展性

CosmicScanner exact Woodbury 路线与 RepeatABEL 的重复测定 GWAS 结果高度一致。虽然原始 beta 方向可能因 allele dosage coding 方向不同出现符号相反，但校正方向后，effect size、standard error 和 association rank 均达到近乎完全一致。

| 比较 | 结论 |
|---|---|
| CosmicScanner exact Woodbury vs RepeatABEL | Pearson r 约 0.999 |
| 运行时间 | CosmicScanner 快 10-20 倍 |
| 内存 | CosmicScanner 约为 RepeatABEL 的十分之一 |
| record count 敏感性 | CosmicScanner 对每个个体记录数增加基本不敏感 |

**图 2. Woodbury 主体水平降维。** 图件见 `docs/figures/figure2_woodbury.pdf`。图中对比 record-level covariance 和 subject-level core matrix。

**图 6. 重复力 GWAS scaling。** 图件见 `docs/figures/figure6_repeatability_scaling.pdf`。增加 repeated records 不显著增加运行时间和内存；增加 subject count 则符合 subject-level dense factorization 的复杂度预期。

### 3.4 纵向 RRM GWAS：准确性和可扩展性

在 GMAT mouse_long benchmark（1,212 个体、19,392 observations、11,833 SNPs）中，CosmicScanner BLUP-based RRM 与 GMAT 的 unbalanced longitudinal result 达到：

| 指标 | 结果 |
|---|---:|
| Pearson r(-log10 P) | 0.951 |
| Spearman rho(P) | 0.940 |
| top SNP | 一致 |
| CosmicScanner runtime | 约 7.9 秒 |
| CosmicScanner peak memory | 约 195 MB |
| observation-level 路线 runtime | 约 1,682 秒 |
| observation-level 路线 memory | 约 17.7 GB |

因此，CosmicScanner 在该 benchmark 上实现 213 倍加速和 93 倍内存降低。对 balanced longitudinal data，GMA-compatible path 与 GMA 的 fixed-regression scan 达到 Pearson r = 1.000。

**图 3. BLUP-based Kronecker RRM。** 图件见 `docs/figures/figure3_blup_rrm.pdf`。

**图 7. RRM 动态效应曲线。** 图件见 `docs/figures/figure7_rrm_dynamic.pdf`。模拟 stable、sign-switch、early-specific、late-specific 和 transient 五类 QTL 轨迹，CosmicScanner 均能恢复并分类。

### 3.5 PSR-GWAS 校准

Ordinal PSR-GWAS 在 null simulation 中保持良好 type-I error calibration。所有扩展 PSR 类模型的 lambda_GC 接近 1.00，经验 type-I error 位于 binomial 95% CI 内。二分类和计数路径仍需在更极端 case-control imbalance 或低事件数情形下进一步扩大验证。

### 3.6 IDUL/AIDUL null-model solver

IDUL、IDUL+ 和 AIDUL 均可减少部分场景下的 null-model 迭代次数。IDUL+ 的主要价值在于收敛稳健性：当基础 IDUL 出现 likelihood decrease 或震荡时，monotone safeguard 能保证更稳定的求解。GEMMA accuracy benchmark 在不同 solver 下保持不变，说明 solver 加速不影响 association scan 的统计结果。

### 3.7 11 类模型覆盖总结

| 模型 | 目标性状结构 | 主要入口 | 验证状态 |
|---|---|---|---|
| Single-trait LMM | 连续性状 | `--gwas` | GEMMA reference |
| Multivariate LMM | 多个相关性状 | `--multivariate` | GEMMA / EMmvGWAS |
| Repeatability | 重复测定 | `--repeatability-model` | RepeatABEL |
| Random regression | 纵向/动态性状 | `--rrm-order` | GMA / GMAT |
| Maternal effect | 母体效应 | `--maternal-model` | type-I + truth recovery |
| ssGWAS | pedigree + genotype | `--ssgwas` | smoke/self-consistency |
| Threshold/probit | 二分类/阈值性状 | model flag | type-I calibration |
| Ordinal PSR | 序数性状 | `--ordinal-model` | calibration |
| Poisson/Binomial PSR | 计数/比例 | `--psr` | type-I calibration |
| Survival Cox | time-to-event | `--survival-model` | type-I calibration |
| GxE reaction norm | 环境互作 | `--gxe-model` | type-I calibration |

外部 reference 覆盖的模型中，single-trait LMM、multivariate LMM、repeatability GWAS 和 balanced RRM 均达到 r > 0.99；unbalanced GMAT longitudinal benchmark 达到 r = 0.951 且 top-hit concordance。ssGWAS 已实现 A22-corrected H inverse block 并通过 smoke/self-consistency validation；BLUPF90/POSTGSF90 外部 benchmark 可作为后续工作。

**图 4. CosmicScanner 与 reference tools 的 accuracy concordance。** 图件见 `docs/figures/figure4_accuracy.pdf`。

### 3.8 扩展模型 type-I error calibration

六类扩展模型在 null simulation 中均表现出接近 1.00 的 genomic inflation factor，经验 type-I error 落在 alpha=0.05 的二项分布 95% CI 内。这说明这些模型作为 release-level 功能已经具备基础统计校准，但部分模型仍需针对具体应用场景做更大规模 external validation。

### 3.9 母体效应模型真值恢复

在 N=500、M=2,000 的模拟中，母体效应模型能将 residual variance 恢复到 <6% 相对误差，并在 null-maternal scenario 中正确估计 maternal variance 接近 0。这说明 direct-maternal covariance 没有简单吸收残差方差。需要注意的是，在 direct 和 maternal GRM 高度重叠时，direct/maternal genetic variance 的精确拆分存在结构性可辨识度限制，这是 GRM-based maternal model 的共同统计问题。

### 3.10 Release-level reproducibility checks

v1.1.0 release candidate 的默认 CTest suite 包含 10 个测试：GEMMA/EMmvGWAS compatibility smoke、GEMMA BSLMM reference validation、GEMMA multivariate reference validation、EB-GWAS quick validation、ssGWAS A22 smoke、sparse REML unit test、multi-variance-component null-model unit test、SPA numerical test、REGENIE unit test 和 input-adaptation unit test。当前 release workstation 上 10/10 通过。该 suite 不能替代完整 benchmark，但可作为发布回归门禁，覆盖主可执行文件、兼容路径、数值核心和输入准备层。

---

## 4. 讨论

### 4.1 主要发现

CosmicScanner 表明，GWAS 软件生态中的三个长期问题可以在一个框架中同时改善：模型碎片化、重复测定/纵向模型扩展性差、现代 genotype format 支持不完整。验证结果支持两个核心观点。

第一，统计一致性和可扩展性并不必然冲突。重复力模型中，Woodbury subject-level reduction 保持精确协方差结构，同时将主要稠密计算从 record level 移至 subject level。RRM 中，BLUP-based coefficient-space strategy 与 GMAT 存在小幅统计差异（r = 0.951），但换来 213 倍运行时间和 93 倍内存收益。

第二，统一框架不等于牺牲模型特异性。CosmicScanner 共享 genotype I/O、GRM 和 scanner 基础设施，但每类模型保留自身统计假设和检验方式。

### 4.2 与现有工具比较

CosmicScanner 与现有工具的主要差别有两点。第一，它在一个框架中统一多个模型类别，提供一致 interface、GRM 和输出 schema。第二，对重复测定和纵向模型，它不是依赖通用代码优化，而是利用模型结构进行计算降维。CosmicScanner 并不是为了取代 BOLT-LMM 或 REGENIE 这类百万样本单性状 GWAS 工具；它的主要价值是复杂性状模型覆盖和重复/纵向模型的可扩展实现。

### 4.3 算法创新意义

Woodbury reduction 对畜牧和植物育种数据尤其重要。若 50,000 个 genotyped subjects 平均每个有 5 条记录，record-level covariance 将达到 250,000 x 250,000；subject-level reduction 则将核心矩阵保持在 50,000 x 50,000。更重要的是，增加重复记录数不会线性推高主要稠密矩阵维度。

BLUP-based RRM 同样改变了纵向 GWAS 的计算可行性。GMAT mouse_long benchmark 如果走 observation-level direct route，需要处理 19,392 x 19,392 的观测水平协方差；CosmicScanner 将核心维度降至 1,212 个 subject 的 coefficient-space scan。

### 4.4 多性状 GWAS 策略

多性状模块提供从严格 reference 对齐到快速筛选的多个层级。EMmvGWAS-compatible route 适合与 EMmvGWAS 对齐；GEMMA-full 适合高精度比较；Hybrid route 适合全基因组扫描中兼顾速度、LOCO 和 selective refinement。

### 4.5 局限性

第一，exact LMM、repeatability 和 RRM 仍依赖 dense factorization，在 5,000-10,000 个体以上会出现计算压力。Biobank two-step mode 和 matrix-free PCG 是更大规模样本的替代路线。

第二，PSR-GWAS 使用 PQL linearization。PQL 在低 case count 或极端 case-control imbalance 下可能有偏；SPA 可部分改善尾部校准，但完整 penalized-likelihood 或 Bayesian approach 是后续方向。

第三，母体效应、ssGWAS、生存和 GxE 当前主要通过 type-I calibration 或 self-consistency 检查验证，尚未与 dedicated external reference tools 做全面 benchmark。本文的外部 reference 重点覆盖 single-trait LMM、multivariate LMM、repeatability、RRM 和 ordinal/PSR calibration。

第四，当前实现主要面向 common variants。虽然已有 SKAT 和 burden tests，针对 millions of rare variants 的 sequence-scale 优化仍需增强。

### 4.6 后续方向

后续开发重点包括：随机化 eigendecomposition 和 stochastic trace estimation；非高斯模型中的 full penalized-likelihood 或 variational Bayes；多组学整合；cloud-native chunked SNP processing；以及更大规模 external benchmark。

---

## 5. 结论

CosmicScanner 证明，GWAS 中模型覆盖广度和计算可扩展性之间的矛盾并非不可避免。通过在一个 C++ 引擎中统一 11 类 GWAS 模型，并对重复测定和纵向模型进行 subject-level 或 coefficient-space 降维，它能够在保持与成熟 reference tools 高一致性的同时，将关键模型的运行时间和内存成本降低一到两个数量级。原生 PLINK/PGEN/BGEN 支持进一步降低了数据转换负担。该框架适用于连续、纵向、重复测定、非高斯和多性状表型的可复现 GWAS 分析。

---

## 数据和代码可用性

- 源代码：https://github.com/huangjie2025/CosmicScanner
- License：MIT License
- 文档：仓库 `docs/` 目录
- 验证脚本：仓库 `tests/` 目录
- 预编译二进制：GitHub Release 页面提供 Linux x86_64 版本

## 作者贡献

JH 构思并设计研究，完成软件实现，开展验证实验，并撰写稿件。

## 致谢

感谢 GEMMA、EMmvGWAS、RepeatABEL、GMA、GMAT 和 PLINK 开发者公开软件和 benchmark 数据。

## 利益冲突

作者声明无利益冲突。

---

## 图表文件索引

- Figure 1：`docs/figures/figure1_architecture.pdf`
- Figure 2：`docs/figures/figure2_woodbury.pdf`
- Figure 3：`docs/figures/figure3_blup_rrm.pdf`
- Figure 4：`docs/figures/figure4_accuracy.pdf`
- Figure 5：`docs/figures/figure5_multivariate_benchmark.pdf`
- Figure 6：`docs/figures/figure6_repeatability_scaling.pdf`
- Figure 7：`docs/figures/figure7_rrm_dynamic.pdf`

