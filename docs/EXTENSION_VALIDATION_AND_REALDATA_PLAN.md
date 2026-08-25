# 扩展模型验证 + 真实数据 GWAS 分析方案

> 本方案对接 `MANUSCRIPT_DRAFT.md` 的 Phase 1，目标是把 Table 7 中 6 个 "—" 缺失项至少补齐 3 项，并新增 1 个真实数据应用章节，使手稿达到 Bioinformatics 投稿门槛。

---

## 一、现状盘点（基于代码核查）

| 模型 | 实现文件 | 参考工具（仓库内） | 现状 |
|------|----------|----------|------|
| 母体效应 | `maternal_model.cpp/.h` | 无 WOMBAT/DMU/MTG2 | 实现完整，未验证 |
| ssGWAS | `single_step_model.cpp` | 无 BLUPF90 | 实现完整，未验证 |
| Survival Cox | `survival_model.cpp/.h` | 无coxme/survival | 实现完整，未验证 |
| GxE reaction-norm | `gxe_model.cpp/.h` | 无 | 实现完整，未验证 |
| Threshold (binary) | `threshold_model.cpp/.h` | 无 SAIGE/POLMM | 已有 SPA，未正式对比 |
| Count (Poisson/Binomial) | `count_model.cpp` | 无 | 实现完整，未验证 |

**关键约束**：`References_and_Tools` 仓库不包含 WOMBAT/DMU/BLUPF90/coxme 等专用参考工具。因此验证策略以 **R/Python 自参考实现 + 模拟已知真值 + 跨方法一致性** 为主，避免引入无法本地运行的外部工具。

---

## 二、扩展模型验证方案

### 2.1 母体效应模型验证（优先级 P0，必做）

**目标**：证明 CosmicScanner 母体模型的方差分量与 SNP 效应与已知方法一致。

**参考实现路径**：用 R `ASReml-R`（若有许可证）或 `lme4qtl`/`sommer` 免费替代。若都不可用，退回到 **解析验证**：构造一个母体效应已知的解析案例，验证估计无偏。

**数据设计**：
- 模拟 $n = 1{,}000$ 个幼畜，500 头母畜，每头母畜平均 2 个后代
- 5 个方差分量配置：
  - $V_a = 0.30$, $V_m = 0.10$, $\sigma_{am} = -0.05$, $V_c = 0.05$, $V_e = 0.55$（典型畜群）
  - 4 个额外配置覆盖高母体/低直接/负相关等场景
- 1,000 个 SNP，其中 5 个因果（直接效应、母体效应、直接+母体）

**验证指标**：
| 指标 | 通过门槛 |
|------|---------|
| 方差分量相对偏差 | < 5%（10 replicates 均值） |
| $\sigma_{am}$ 相对偏差 | < 10% |
| 因果 SNP top-10 命中 | ≥ 4/5 |
| $\beta$ 与真值相关 | > 0.95 |
| 与 R 参考实现 $\beta$ 相关 | > 0.99（若 R 参考可用） |

**输出**：新增 Results 3.7.1，Table 7 中母体行从 "—" 改为 "R sommer / 模拟真值"。

---

### 2.2 ssGWAS 验证（优先级 P1，强烈建议）

**目标**：证明 CosmicScanner ssGWAS 的 H 矩阵构造与 SNP 效应估计正确。

**参考实现路径**：`References_and_Tools` 无 BLUPF90。采用 **两阶段自洽验证**：

1. **H 矩阵正确性**：构造小规模 pedigree（30 个体，3 代），其中 10 个基因型个体。手工构造 $\mathbf{A}$、$\mathbf{A}_{22}$、$\mathbf{G}$，验证 $\mathbf{H}^{-1} = \mathbf{A}^{-1} + \text{diag}(\mathbf{G}^{-1} - \mathbf{A}_{22}^{-1})$ 与 CosmicScanner 输出一致（最大相对误差 < 1e-8）。
2. **GEBV 一致性**：构造全基因型场景（pedigree 全部基因型），此时 ssGWAS 应退化为标准 GBLUP。比较 CosmicScanner ssGWAS 与 CosmicScanner LMM 在全基因型数据上的 GEBV，相关应 = 1.000。
3. **SNP 效应**：用 GEBV 反推 SNP 效应（`ssGWAS --snp-effect`），与直接 GBLUP SNP 效应比较，相关 > 0.999。

**数据设计**：
- Pedigree：500 个体，3 代，100 个基因型个体（20% 基因型比例）
- 1,000 SNP，3 个因果

**验证指标**：
| 指标 | 通过门槛 |
|------|---------|
| H 矩阵最大相对误差 | < 1e-8 |
| 全基因型退化 GEBV 相关 | > 0.9999 |
| 部分基因型 GEBV 与 GBLUP-only 相关 | > 0.95 |
| SNP 效应相关 | > 0.999 |

**输出**：Results 3.7.2。

---

### 2.3 Threshold / Binary 验证（优先级 P1，与 SAIGE 风格对比）

**目标**：证明 binary GWAS 在 case-control 不平衡下校准。

**参考实现路径**：`References_and_Tools` 无 SAIGE。采用 **R `coxph`/`glm` + SPA 自实现** 对照：
- 平衡场景（case:control = 1:1）：与 R `glm` + LMM（`coxmeg` 或 `lme4qtl`）对比
- 不平衡场景（1:10, 1:100）：仅与 SPA 校准理论曲线对比

**数据设计**：
- $n = 2{,}000$，$M = 5{,}000$
- 3 个因果 SNP（MAF=0.1/0.3/0.05），OR=1.3/1.5/2.0
- 4 种 case 比例：0.5, 0.2, 0.1, 0.05
- 每场景 10 replicates

**验证指标**：
| 指标 | 通过门槛 |
|------|---------|
| $\lambda_{\text{GC}}$（空模型） | 0.95–1.05 |
| 与 R `glm` $\beta$ 相关（平衡） | > 0.99 |
| SPA 校准后 $P$ 值分布 KS test | $p > 0.05$ |
| 因果 SNP power（$\alpha = 5\times10^{-8}$） | ≥ 0.8（OR≥1.5） |

**输出**：Results 3.7.3。

---

### 2.4 GxE Reaction-Norm 验证（优先级 P2，可选项）

**目标**：证明 GxE 模型能正确识别交互效应。

**参考实现路径**：与 R `lme4` 的 `(snp + env + snp:env)` LMM 对比。

**数据设计**：
- $n = 1{,}500$，环境变量为连续值（如温度 5–35°C，标准化）
- 2 个因果：主效应 + 交互；纯交互；纯主效应
- 1,000 SNP

**验证指标**：$\beta_0$、$\beta_1$ 与 R 相关 > 0.99；2-df 联合 $P$ top-5 命中 ≥ 4/5。

---

### 2.5 Survival Cox 验证（优先级 P2，可选项）

**参考实现路径**：与 R `coxme`（frailty 模型）对比。

**数据设计**：
- $n = 1{,}500$，50% 删失
- 模拟 Weibull 基线风险，3 个因果 SNP（HR=1.2/1.5/2.0）

**验证指标**：$\beta$ 与 `coxme` 相关 > 0.95；top-10 命中 ≥ 3/3。

---

### 2.6 Count Model 验证（优先级 P3，可选项）

**与 Poisson `glmer` 对比，模拟计数性状（子女数、窝产仔数）。门槛：$\beta$ 相关 > 0.95。**

---

## 三、真实数据 GWAS 应用方案

### 3.1 推荐数据集（按可得性与适配度排序）

| # | 数据集 | 样本量 | 性状类型 | 验证模型 | 获取难度 |
|---|-------|-------|---------|---------|---------|
| **1** | **UKB 子集**（身高/BMI） | 50K | 连续 | LMM | 低（公开申请） |
| **2** | **1000G + HGDP**（已知身高位点） | ~3K | 连续 | LMM/mvLMM | 低（公开） |
| **3** | **Cattle GT30K**（产奶量/体高） | 30K | 重复/纵向 | Repeatability/RRM | 中（合作） |
| **4** | **Arabidopsis 1001G**（开花时间） | ~1K | 纵向/环境 | RRM/GxE | 低（公开） |
| **5** | **Mouse HS**（体重/血糖） | ~2K | 多性状 | mvLMM | 中（合作） |

**推荐组合**：
- **必做**：UKB 子集（50K）做单性状 LMM——证明 biobank 可用性 + 找回身高已知位点
- **加分**：Arabidopsis 1001G 做 RRM/GxE——证明纵向+环境交互在真实植物数据上有效
- **保底**：1000G+HGDP 做身高——即使 UKB 申请受阻也能完成

### 3.2 UKB 50K 子集分析方案（主推）

**样本筛选**：
- 从 UKB 英格兰队列随机抽 50K 欧洲血统个体
- 排除性染色体异常、亲缘密切个体（`--king-cutoff 0.044`）

**基因型预处理**：
- 染色体 1–22，MAF ≥ 0.01，HWE $P > 1\times10^{-6}$，缺失率 < 0.02
- 用 PLINK 2.0 `--make-pgen` 输出，测试 PGEN 原生读取
- 预计 ~600K SNP

**性状**：身高（field 50）、BMI（field 21001）

**分析流程**：

```
# Step 1: 基因型 QC + PGEN 转换
plink2 --bfile ukb_cal --keep ukb50k.sample --maf 0.01 \
       --hwe 1e-6 --geno 0.02 --make-pgen --out ukb50k_qc

# Step 2: GRM 构造（LOCO，chrom 1-22）
cosmicscanner --pgen ukb50k_qc --grm-build --loco --out ukb50k_grm

# Step 3: 单性状 LMM（身高），LOCO，fast-LMM
cosmicscanner --pgen ukb50k_qc --pheno ukb50k_pheno.txt --mpheno 1 \
              --qcovar age.txt --covar sex.txt --gwas --loco --fast-lmm \
              --out ukb50k_height

# Step 4: 多性状 LMM（身高+BMI）
cosmicscanner --pgen ukb50k_qc --pheno ukb50k_pheno.txt \
              --multivariate --mv-traits 1,2 --mv-method hybrid \
              --qcovar age.txt --covar sex.txt --out ukb50k_height_bmi

# Step 5: biobank 两步法（验证可扩展性，全 50K）
cosmicscanner --pgen ukb50k_qc --pheno ukb50k_pheno.txt --mpheno 1 \
              --biobank-step1 --step1-blocksize 1000 --out ukb50k_step1
cosmicscanner --pgen ukb50k_qc --pheno ukb50k_pheno.txt --mpheno 1 \
              --biobank-step2 --use-pred ukb50k_step1.pred \
              --qcovar age.txt --covar sex.txt --out ukb50k_height_2step
```

**验证指标**：
| 指标 | 通过门槛 |
|------|---------|
| 已知身高位点恢复 | HMGA2、GH1、GDF5 等 top-50 至少命中 5 个 |
| $\lambda_{\text{GC}}$ | 1.0–1.10 |
| 与 REGENIE 公开 summary 比对（同 UKB 子集） | $\beta$ 相关 > 0.98 |
| Biobank 两步法 vs exact LMM $\beta$ 相关 | > 0.99 |
| 50K 全量运行时间 | < 30 min（单机） |

**预期 figure**：曼哈顿图 + QQ plot（身高 vs BMI）。

### 3.3 Arabidopsis 1001G RRM/GxE 方案（加分项）

**数据**：1001G 公开 VCF + 开花时间表型（10 个环境/时间点）

**分析**：
```
# RRM（纵向开花时间）
cosmicscanner --pgen ara1001 --pheno flowering_long.txt --tpheno 2 \
              --rrm-order 2 --tclass env --loco --out ara_flowering_rrm

# GxE（温度交互）
cosmicscanner --pgen ara1001 --pheno flowering.txt --mpheno 1 \
              --gxe-model --env-covar-col 3 --out ara_flowering_gxe
```

**验证指标**：识别已知开花时间基因（FLC、FRIGIDA）；动态效应曲线合理。

### 3.4 1000G 保底方案

若 UKB 申请受阻，用 1000G 欧洲子集（~500）+ HGDP 表型做身高 GWAS，方法同 3.2 但样本量小一个量级，仅作为 "real-data proof-of-concept"。

---

## 四、新增 Results 章节规划

| 新章节 | 内容 | 对应 Phase |
|--------|------|----------|
| 3.7 扩展模型验证 | 6 个子节（母体/ssGWAS/threshold/GxE/survival/count） | Phase 1.1 |
| 3.8 真实数据应用 | UKB 50K 身高/BMI 曼哈顿、QQ、已知位点恢复 | Phase 1.2 |
| 3.9 Type-I error 与 power | 4 个核心模型（LMM/mvLMM/Repeatability/RRM）的空模型 λ_GC、power 曲线 | Phase 1.3 |

**手稿整体长度预估**：从现 ~10 页扩到 ~18–22 页（含图表），适合 Bioinformatics Article 体例。

---

## 五、时间与执行优先级

| 优先级 | 任务 | 工作量 | 收益 |
|-------|------|-------|------|
| **P0** | 母体模型验证（模拟+R参考） | 2–3 天 | Table 7 关键红杆消除 |
| **P0** | type-I error / power 校准（LMM/mvLMM/Repeatability/RRM） | 2–3 天 | 审稿人必问，影响所有声明 |
| **P0** | UKB 50K 真实数据 GWAS | 5–7 天（含数据申请） | 决定能否投 Bioinformatics |
| **P1** | ssGWAS 自洽验证 | 1–2 天 | 扩展模型覆盖度 |
| **P1** | threshold/SPA 不平衡校准 | 2 天 | binary 声明可信度 |
| **P2** | GxE / survival / count | 各 1–2 天 | 锦上添花 |

**最小可投 Bioinformatics 集合**：P0 三项 + P1 一项 = 约 2 周完成。

---

## 六、风险与对策

| 风险 | 对策 |
|------|------|
| UKB 数据申请周期长（4–8 周） | 立即启动申请；同时用 1000G 做 PoC |
| R 参考工具（ASReml/coxme）许可证 | 优先 `sommer`/`lme4qtl`/`coxme`（CRAN 免费） |
| 50K UKB 内存不够 | 用 biobank 两步法 + 稀疏 GRM（`--grm-sparse`） |
| 母体模型 R 参考与 CosmicScanner 实现有差异 | 以模拟真值为基准，R 参考作为辅助证据 |
| 性能 benchmark 硬件环境未规范 | 所有实验统一在 `c6a.4xlarge` 或本地 32-core 机器上跑，记录 CPU/RAM/OS/编译选项 |

---

## 七、下一步行动

1. **立即启动**：UKB 数据申请（若未申请），并行用 1000G 跑 PoC
2. **本周内**：完成母体模型模拟验证脚本 + R 参考实现
3. **下周内**：完成 type-I error / power 模拟框架
4. **2 周内**：补齐 ssGWAS + threshold，UKB 数据到位后跑真实数据
5. **3 周内**：完成全部实验，更新手稿 Results，进入润色阶段

需要我现在就开始起草母体模型验证脚本（模拟+R参考），还是先做 type-I error 校准框架？
