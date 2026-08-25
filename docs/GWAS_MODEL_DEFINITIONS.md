# GWAS 模型统计定义 Audit（模型语义定义）

本文档把 CosmicScanner 的每个扩展 GWAS 模型严格拆分为五个层次：

1. **Genetic evaluation model**（若无 → 不适用）
2. **GWAS background covariance model**
3. **SNP-effect design**
4. **Association test**
5. **Primary output**

并给出 **Association type** 分类（A/B/C）与 **Naming accurate?** 判定。

> Audit 方法：沿 CLI → model dispatch → null model → covariance construction → SNP design matrix →
> association statistic → output 逐层追踪源码（`src/models/*.cpp`, `src/gwas_rrm_scan.cpp`,
> `src/gwas_app.cpp`, README/Tutorial 仅作旁证）。本轮为纯审计 + 文档；**未改任何算法代码**。

---

## 统一分类框架

| 类型 | 定义 | 本仓库成员 |
|------|------|-----------|
| **Type A** | complex-covariance-adjusted GWAS：`y = Xb + x_j β_j + random effects + e`，SNP 是普通 direct fixed SNP effect，复杂模型只控制背景相关性 | Standard LMM、Animal Model、Repeatability、Common-environment、Maternal（当前实现） |
| **Type B** | model-specific SNP-effect GWAS：SNP 设计本身随模型改变（母亲基因型、SNP×时间基、SNP×环境、多性状联合 SNP 效应） | RRM joint scan（平衡路径 + GMA-fixreg 路径）、Multivariate、GxE、**maternal SNP（`--maternal-snp-test maternal|joint`，已实现）** |
| **Type C** | breeding-value/coefficient derived GWAS：先 BLUP/REML 得到系数/GEBV，再对系数或 GEBV 扫描或 backsolve | RRM 系数空间辅助扫描（非平衡主导路径）、（未来）classic ssGWAS backsolve |

---

## 1. Standard LMM (UnivariateModel / LmmEngine)

- Background：`y = Xb + u + e`, `u ~ N(0, G·σ²g)`（VanRaden/Zeng/Yang/Vitezica/TrueVar GRM）
- SNP design：个体自身 genotype `x_i,j`（中心化后）
- Test：`H0: β_j = 0`（score/Wald，LOCO 可用）
- Primary output：`BETA SE P LOG10P`
- Association type：**ordinary fixed SNP**（Type A）
- Naming：**标准 LMM GWAS** ✓

## 2. Animal Model（--animal-model）

- Background：`y = Xb + x_j β_j + a + e`, `a ~ N(0, A22·σ²a)`（A22 为分析样本系谱关系子矩阵，Schur-complement 求解，不物化稠密 A）
- SNP design：个体自身 genotype（LmmEngine 标准扫描）
- Test：`H0: β_j = 0`；LOCO 不适用（A 与染色体无关），多染色体未 `--no-loco` 时 warning
- Primary output：`*.assoc.txt` + `*.animal_model.null.txt`（σ²a/σ²e/h2）+ 可选 `*.animal_model.k.txt`/`.k.id`（A22 核，GEMMA `-k` 可用）
- Association type：**ordinary fixed SNP**（Type A）
- Naming：**pedigree-kinship LMM GWAS / Animal Model-adjusted GWAS** ✓（准确）
- 验证：模拟套件 7 段 + GEMMA -k A22 核外部基准（β/SE/-log10P corr ≈ 0.99996/0.9992/0.9984）

## 3. Repeatability GWAS（--repeatability-model）

- Background：`y = Xb + x_j β_j + Z_a a + Z_pe pe + e`, `a ~ N(0, G·σ²g)`, `pe ~ N(0, I·σ²pe)`（记录级 Woodbury 精确分解 / exact GLS）
- SNP design：个体自身 genotype，按记录重复展开（一次记录一份 genotype）
- Test：`H0: β_j = 0`（记录级 score/Wald；exact tier 逐 SNP 全 VC refit）
- Primary output：`*.assoc.txt` + null
- Association type：**ordinary fixed SNP**（Type A）
- Naming：**repeatability-adjusted GWAS / GWAS with permanent-environment random effect** ✓
- 外部基准：RepeatABEL 对齐（r>0.99）

## 4. Common-environment GWAS（--common-env-model）

- Background：`y = Xb + x_j β_j + Z_a a + Z_c c + e`, `a ~ N(0, K·σ²a)`, `c ~ N(0, I_g·σ²c)`（组共享随机截距）
- SNP design：个体自身 genotype（记录级展开）
- Test：`H0: β_j = 0`（`RecordLevelNullModel::testSnp` score/Wald，`src/models/multi_random_effect_lmm.cpp:152`）
- Primary output：`*.assoc.txt` + null
- Association type：**ordinary fixed SNP**（Type A）
- Naming：**shared-environment-adjusted GWAS** ✓

## 5. Maternal Model（--maternal-model）【审计重点】

### 5.1 当前实现（只读代码确认，`src/models/maternal_model.cpp`）

- 随机效应模型（VCE 阶段）：`y = Xb + Z_d a_d + Z_m a_m + (Z_c c) + e`
  - `a_d ~ N(0, G_d·σ²ad)`，`G_d` = 子代个体基因组 GRM
  - `a_m ~ N(0, G_m·σ²am)`，**当前 `G_m = G_d`**（`buildGRMs()` 明确打印 "Maternal GRM uses direct GRM (full implementation pending)"，行 201-204）
  - `Cov(a_d, a_m)` 参数 `σ_adm` 存在但**未被估计**：AI-REML 是 2 分量（合并）或 3 分量（σ²ad/σ²am/σ²e，无协方差块），σ_adm 恒为初值 0（行 252-254, 312-446）
  - `merged_grm` 阈值（`‖K−K_m‖/‖K‖<1e-4`）在 `G_m = G_d` 时必然触发 → 实际合并为单遗传分量、`σ²am→0`

- **SNP 扫描阶段（`scanGwas()` / `runLocoScan()`）**：
  - 先把记录级表型**按个体平均**为 `y_agg`、协变量同样聚合（行 637-660）
  - 用**个体自身 genotype**（`LmmEngine` 全 SNP 扫描 `all_indices = 0..n_animal-1`）对 `y_agg` 做标准 LMM 扫描
  - 方差分量坍缩为单核：`vg_eff = σ²ad + σ²am + 2σ_adm (+σ²c)`，`ve_eff = σ²e`，扫描协方差 `vg_eff·G_sub + ve_eff·I`
  - 输出：`*.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`）+ `*_maternal_null.txt`（σ²ad/σ²am/σ_adm/σ²c/σ²e）

### 5.2 严格判定

| 问题 | 答案 |
|------|------|
| A. 扫描的是 `x_i,j·β_direct,j`（子代自身 genotype, direct SNP association）？ | **是**。SNP 设计是个体自身 genotype，对个体级平均表型做标准 LMM 扫描。 |
| B. 是否使用 `x_dam,j·β_m,j`（母亲 genotype 的 maternal SNP effect）？ | **否**。母亲 genotype 从未进入 SNP 设计；dam 映射只进入 VCE 的 `K_m`（且当前 `G_m=G_d`）。 |
| C. direct + maternal 双参数 SNP 模型？ | **否**。不存在 `x_iβ_d + x_damβ_m` 双列设计（任意 SNP 均单列个体 genotype）。 |
| D. 2-df joint test `H0: β_d=β_m=0`？ | 默认 `direct` 模式**否**（单列）；`--maternal-snp-test joint` 下**是**（见 5.3）。 |
| E. 最准确名称？ | **"maternal-adjusted GWAS" / "GWAS with maternal random effects"**。直接声称 "Maternal-Effect GWAS"（可检验母亲基因型 SNP 效应）**不准确**。 |

**结论：默认路径（`--maternal-model` = `--maternal-snp-test direct`）是 Type A（背景协方差含 direct+maternal 随机效应的 GWAS），其 `*.assoc.txt` 的 BETA/SE/P 语义 = 子代个体的 direct fixed SNP 效应，与标准 LMM 同级。** 自 `--maternal-snp-test maternal|joint` 实现后，maternal_model 同时覆盖 Type B（母亲基因型单列检验 / direct+maternal 2-df joint），见 5.3。

### 5.3 实现：model-specific SNP tests（`--maternal-snp-test direct|maternal|joint`，已实现）

```
--maternal-snp-test direct     y = Xb + x_offspring,j β_d,j + Z_d a_d + Z_m a_m + e      → 现状 scanGwas（等价回归已验证 byte-identical）
--maternal-snp-test maternal   y = Xb + x_dam,j β_m,j + random effects + e               → scanMaternalSnp（`*.maternal.assoc.txt`）
--maternal-snp-test joint      y = Xb + x_off,j β_d,j + x_dam,j β_m,j + random effects + e
                                H0: β_d,j = β_m,j = 0（2-df Wald/Score，`*.joint.assoc.txt`）
                                输出 BETA_DIRECT SE_DIRECT BETA_MATERNAL SE_MATERNAL STAT_JOINT DF P_JOINT COND RANK STATUS
```

**实现要点（`src/models/maternal_model.cpp` `scanMaternalSnp()`）**：
- **eligible subset**：只保留"有表型记录 **且** dam 为已基因型个体（subject index < n_samples）"的动物；`N_TOTAL / N_WITH_DAM_GENOTYPE / N_DROPPED_MISSING_DAM` 三列写入 assoc 与 `*.maternal.snp.meta.txt`（缺失策略 = `filter_to_eligible_subset`）。BGEN/PGEN 暂未接线（PLINK-only，入口 throw）。
- **SNP 设计**：同一染色体位置 `x_off`（子代个体 ALT dosage）与 `x_dam`（dam 的 ALT dosage，BED code 00=2.0/10=1.0/11=0.0，01=缺失）；逐 SNP 均值中心化、缺失按 0 填补（与 LMM engine 一致）。
- **Null 协方差**：与 maternal VCE 一致的全分量协方差 `V = σ²ad·K(offspring) + σ²am·K(dams) + σ²e·I`（dam 行取同一 subject-indexed keep-GRM；每染色体 LDLT 分解一次复用）。`σ_adm` 未被 VCE 估计（恒定 0）；`σ²c` 仅当 `--maternal-pe` 且 >0 时加对角。
- **检验**：score = S'P_reml y；M = S'V⁻¹S（maternal 1 列 / joint 2 列，pseudo-inverse + eigenvalue 截断 + rank/cond 诊断）；maternal 单列 Wald `beta²/var`；joint 2-df Wald `stat = sc' M⁻¹ sc`，`df=rank`。`rank<2` 或 `cond>1e6` → `STATUS=COLLINEAR`。maternal 模式 `M(0,0)≤1e-12` 的 SNP 跳过。
- **LOCO**：与 direct 相同策略（`G_loco = G_global·S_global/S_loco − G_chr/S_loco`，subject-indexed），母/子两列用同一 LOCO 校正 GRM。
- **VCE 注意**：existing scanGwas（direct 路径）以"记录位置==subject 顺序"约定运行；`--maternal-snp-test` 的 eligible 逻辑按 subject id 索引 GRM，不依赖行序假设。

**验证（`tests/run_maternal_snp_effect_simulation.py` + `tests/run_test_maternal_snp_effects.sh`，seed=7, n=250, m=1200, nf=120 → 10/10 PASS）**：
| 场景 | 真值 | 结果 | 判定 |
|------|------|------|------|
| A_null | β_d=β_m=0（含 σ²am 母体背景） | λ_maternal=1.061, λ_joint=1.092 | PASS（maternal/joint 均无膨胀） |
| B_direct_only | (0.8, 0) | direct p=7.6e-4 / maternal 弱 / joint p=4.3e-11 | PASS |
| C_maternal_only | (0, 0.8) | maternal p=4.3e-8 / direct 弱 / joint p=1.7e-7 | PASS |
| D_same_dir | (0.8, 0.5) | direct/maternal/joint 全部显著 | PASS |
| E_opposite | (0.8, −1.2) | joint 最强（−log10 17.8 > marginals） | PASS |
| F_dm_cov | β=0, σ_adm=−0.06 | λ_direct/maternal/joint 均 <1.45（协方差不产生虚假信号） | PASS |
| G_missing10/30/50 | (0,0.8), dam 缺失 10/30/50% | meta 正确报告 dropped；maternal p<0.05 且可用 N 递减 | PASS |
| H_self_dam | 子代即母体 | joint 100% COLLINEAR（rank<2），诊断如预期 | PASS |

- Status：**simulation-validated（无外部参考）；不声称 production 级验证**。
- 外部参考：sommer `mmer(u=~animal+dam)` 为**背景 VCE/育种值**参考（`tests/run_maternal_validation.py` 只对 direct 效应 β 对齐）；真正 maternal SNP effect 外部参考（Ragsdale-Gravel 型双基因型 GWAS 或等效 two-genotype LMM）不存在于本地工具集 → 标记 "simulation-validated only"，**不伪造外部验证**。

## 6. RRM-GWAS（--long-pheno / --rrm-order）【审计重点】

### 6.1 平衡路径（`runRrmScan`, `src/gwas_rrm_scan.cpp:576`）——两条并存输出线

**线 1（主）：观测水平联合扫描 `.rrm_joint.assoc.txt` —— 真 · time-varying SNP 固定效应模型**

- 模型：`y_it = Xb + φ(t)' γ_j · x_ij + (random genetic trajectory) + (permanent env trajectory) + e_it`
  - SNP 固定设计：`x_i,j × φ(t)`（Legendre/B-spline 基函数向量），系数向量 `γ_j`（维数 = order+1）——**验证：`scan_plink_rrm_joint`（行 446-506）中 `score(a) = score_weights[a]·subject_g`，`score_weights[a]` 由 `Φ_obs(:,a)·P_y` 逐观测聚合到个体（行 205-210），info `M(a,b)` 由 `P` 矩阵与 `Φ_obs` 的外积聚合（行 212-233）→ 这正是 `x⊗Φ` 设计下的 score/info 积分**
  - Test：`stat = score'M⁻¹score`, `p = χ²(df = order+1)` → **`H0: γ_j = 0`（2/3/…-df 联合 Score/Wald）**
  - β(t) 重构：`.rrm.dynamic_curve.txt` 用 `β(t)=φ(t)'γ`、`SE(t)=√(φ'M⁻¹φ)`（行 1169-1178，**严格二次型，正确**）
  - 背景协方差：`V = (G ⊗ ΦKΦ') + (I ⊗ ΦPΦ') + σ²e I`（`build_rrm_observation_kernel`，行 151-183）→ **genetic = G⊗D、PE = I⊗P，代码等价形式成立**
- 输出语义：`BETA_Lk/SE_Lk` = SNP×基 γ_k 与其 SE；`JOINT_CHISQ/JOINT_P` = γ joint 检验

**线 2（辅助）：系数空间/逐基扫描**
- 每个个体做**个体内最小二乘** Legendre 拟合得到系数 `coef_mat(:,k)`（类 BLUP 系数，行 664-693），随后用标准 GRM-LMM（`runLocoScan`）对每个系数作表型扫描 `.Lk.assoc.txt`；`.rrm.assoc.summary.txt` 做 best-basis + 逐基 z² 合并；`.rrm.top_curve.txt` 用逐基 β/SE **对角线近似**重构曲线（行 1250-1278，忽略基间协方差，是近似曲线）
- 该方法 = "two-step coefficient GWAS"（Type C/derived-phenotype），**不是** `x⊗basis` 固定效应检验

**线 3：GMA-compatible fixreg（`--rrm-compatible gma-fixreg`, `rrm_model.cpp:454`）**
- 逐 SNP 建 `X[i] = [1·L | g_i·L]`（个体 i 每时间点的 design = [固定项基 | SNP×基]），对 d 维 `β_snp = tail(L 系数)` 做 **Wald χ²(df=d)**，输出 `chi_val/p_val/BETA_Lk`
- 这是**真 · SNP×时间基固定效应联合检验**，与 GMA `balance_longwas_fixreg` 作**外部对齐**（`run_gma_compatible_fixreg_alignment.sh` 对 chi_val/p_val 做 Pearson 比较）——即该外部 benchmark 对应 **SNP association statistic（Type B joint），不是 BLUP 系数**。

### 6.2 非平衡路径（`--rrm-unbalance`, `runRrmUnbalancedScan`, 行 1288）——两条并存扫描线（`--rrm-gwas-method record|coefficient`）

- **Method A（`--rrm-gwas-method record`，新增主扫描线）**：记录级固定效应联合检验。
  - SNP 固定设计：`D(r,k) = x_{subject(r),j} · φ_k(t_r)`（record×basis），**真 · time-varying SNP effect（Type B）**，与平衡 joint 线相同语义、但直接处理非平衡记录与个体缺失时间点。
  - Test：`score = D'P_reml y`、`M = D'V⁻¹D`，eigen 截断（df=rank≥1）后 `stat = score'M⁺score`，`p = χ²_df`；`STATUS ∈ {OK, RANK_DEFICIENT, ILL_CONDITIONED, NUMERICAL_FAILURE}`（不静默丢 SNP）。
  - 验证：平衡显式 joint 路径为该法的隐式等价形式（`run_rrm_balance_record_equivalence.sh` 上 BETA/SE/STAT/log10P Pearson ≥ 0.99999）；A–G 模拟套件（causal rank=1、轨迹相关 0.995–0.999、λ_GC≈0.92）与 A-vs-B 对比（`run_rrm_method_ab_comparison.sh`）见 `docs/RRM_GWAS_METHOD_AUDIT.md` §5。
- **Method B（`--rrm-gwas-method coefficient`，默认，backward compatible）**：BLUP 系数空间多性状扫描。
  - `Y_blup = coef_mat − coef_mean`（个体 Legendre 系数矩阵作为**多性状表型**）
  - `MvEigenVce` AI-REML（Vg, Ve，核 = 个体 GRM）→ `MvLmmEngine` 多性状 Wald（β 为逐基多性状 SNP 效应向量，`H0: β_L0=…=β_Lq=0`）
  - 动态曲线文件对 top SNP 用 `β(t)=φ(t)'β_snp`, `SE(t)=√(φ'Cov_snp φ)`（行 1822-1831，基于多性状 Wald 协方差，合理代理）
- 这是**"先 BLUP 系数、后对系数扫描"**的两阶段方法（Type C）。其 SNP 固定设计是**普通 genotype 对多性状系数表型**，**不是** `x⊗Φ` 固定效应设计。
- **GMAT 外部基准（`run_gmat_unbalanced_alignment.sh`）比较的是 `JOINT_P`（BLUP 系数空间多性状 Wald P）与 GMAT 输出的 `p`** → 该 benchmark 对应的是 **"coefficient/BLUP-derived 多性状 SNP 检验"**，不是观察水平 `x⊗Φ` 固定效应检验；常染色体 Pearson=0.9521/Spearman=0.9401、top SNP `UNC18848064` 双一致（全序含 X 染色体为 0.9470/0.9393，略降亦正体现两套方法在非平衡下并非同一种检验）。X 染色体（chr23）因 GMAT 排除 X 出 GRM 而存在方法学差异（审计注记见 `docs/RRM_GWAS_METHOD_AUDIT.md` §4.3）。
- **非平衡方法选择**：需严格 time-varying SNP 效应 → `--rrm-gwas-method record`；需与 GMAT/历史 `--rrm-unbalance` 衔接 → `--rrm-gwas-method coefficient`（默认）。两者均合法，但统计语义不同，文档须分别解释。

### 6.3 严格判定

| 问题 | 答案 |
|------|------|
| A. SNP 固定设计是 `x_j × basis(t)`？ | 平衡 joint 线：**是**；GMA-fixreg 线：**是**；非平衡 Method A（`--rrm-gwas-method record`）：**是**（记录级 `x·φ`）；非平衡 Method B（默认 coefficient）：**否**（普通 genotype × 系数表型）。 |
| B. 估计 γ_j0…γ_jq？ | 平衡 joint 线输出 `BETA_Lk`（=γ_k）与 SE：**是**；GMA-fixreg 输出同名列：**是**。 |
| C. β_j(t) 由 SNP-specific γ 重构？ | 平衡 joint 线 `.rrm.dynamic_curve.txt`：**是**（二次型 SE）；`.rrm.top_curve.txt`（辅助系数线）：对角线近似。 |
| D. joint test 检验 `H0: γ_j=0` 而非 BLUP 系数？ | 平衡 joint 线、GMA-fixreg：**是**（score/Wald on γ）。 |
| E. 每时间点 SNP 效应/SE？ | 平衡 joint：β(t)=φ'γ、SE(t)=√(φ'Covφ)；系数辅助线：对角线近似。 |
| F. genetic = G⊗D、PE = I⊗P？ | 平衡观察水平 null：**是**（`Φ·K·Φ'` 元素级 × G、`Φ·P·Φ'` 同个体）。 |
| G. GMA/GMAT benchmark 对应 SNP 关联检验还是 BLUP 系数？ | GMA-fixreg：SNP `x⊗Φ` 联合检验（真关联统计）。GMAT：BLUP 系数空间多性状检验（derived/Type C）。**二者不混** —— 需在文档中明确。 |

### 6.4 命名

- 平衡路径主输出：**longitudinal dynamic GWAS / RRM-GWAS（Type B）** ✓
- 非平衡路径与系数辅助线：**coefficient-based (BLUP) two-step GWAS** —— 合法但必须按此名称解释，不能伪装成 `x⊗Φ` 固定效应动态 GWAS。
- 现状文档（Manuscript_CN）已如实写 "先估计个体 Legendre 系数，再在系数空间扫描 SNP"——与实现一致。

## 7. Multivariate GWAS（--multivariate）

- Background：`Y = XB + U + E`, `Vec(U) ~ N(0, G⊗Σg)`
- SNP design：普通 genotype，多性状效应向量 `β_j ∈ R^t`（EMmvGWAS-compatible / GEMMA-full / GEMMA-VC / hybrid）
- Test：多性状 Wald（`H0: β_j = 0` over t traits）
- Association type：**multi-trait joint SNP**（Type B）
- Naming：**multivariate LMM joint SNP GWAS** ✓；外部基准 GEMMA/EMmvGWAS r>0.99

## 8. ssGWAS（--single-step）【最高优先级审计】

> **完整方法/公式核查见 `docs/SSGWAS_METHOD_AUDIT.md`**：CURRENT（H-kernel association，Method 1）vs CANONICAL（经典 backsolve ssGWAS，Method 2）逐特性对照、H⁻¹ 构造审计、back-solving 公式（Wang 2012; Aguilar 2019 已联网核实）、SNP/window 方差、iterative weighted（Method 3，仅设计）。本节为结论摘要。

### 8.1 H 构造（代码确认，`single_step_model.cpp`）

```
H⁻¹ = A⁻¹ + [ 0            0
               0   τ·G⁻¹ − ω·A22⁻¹ ]
```
- `A⁻¹`：Henderson 系谱求逆（`pedigree_relationship.h`）
- `A22`：分析基因型子块的系谱关系（Schur-complement，不物化稠密 A），再岭求逆
- `G`：基因组 GRM（可选 `--blend-g` α·I 混合）
- 对称性 max-abs 检查写入 `*.agh.meta.txt`

### 8.2 之后实际做了什么（逐行追踪）

1. **VCE（`runVCE`）**：`H_sub = (H⁻¹ 的基因型子块求逆)` → `EigenVce` REML 在 (y_sub, X_sub, H_sub) 上估计 `σ²g, σ²e`。**只使用基因型样本**（`ng = keep.n_keep` 必须相等）；未基因型动物不进入 VCE。
2. **扫描（`scanGwas`）**：`EigenVce` 固定 (σ²g, σ²e) → `LmmEngine` 对**个体自身 genotype** 做标准 LMM score/Wald 扫描，核 = H-块（LOCO 时按染色体重建 `buildHSubForG(G_loco)`）。
3. **输出**：`*.assoc.txt`（`BETA SE P LOG10P`——与标准 LMM 完全同 schema）+ `*_ssgwas_null.txt` + 可选 `*.Ainv/G/Hinv.bin`.

### 8.3 与经典 ssGWAS 的二分法

| | Method 1：H-kernel mixed-model association | Method 2：classic breeding-value / backsolve ssGWAS |
|---|---|---|
| 输入→ | `y = Xb + x_j β_j + u + e`, `u~N(0,Hσ²)` | pedigree+genotype+phenotype → ssGBLUP → EBV/GEBV |
| 检验 | 逐 SNP `H0: β_j=0`（fixed SNP effect, Wald/score） | 无逐 SNP 固定效应检验；或对 GEBV 反解 SNP 效应/窗口方差 |
| 产出 | BETA/SE/P | SNP effects/weights、SNP variance、sliding-window %var、迭代加权 |

### 8.4 严格判定

| 问题 | 答案 |
|------|------|
| A. 属于哪种？ | **Method 1（H-kernel fixed-SNP association / single-step-informed LMM）**。 |
| B. 计算 EBV/GEBV？ | **否**。`writeEbv()` 为 stub（打印 "GEBV output will be implemented in v1.1."）。 |
| C. SNP back-solving？ | **否**。无 `ĝ = M'(MM')⁻¹ û_g` 类代码。 |
| D. SNP weights？ | **否**。 |
| E. iterative weighted ssGBLUP？ | **否**。 |
| F. sliding-window explained variance？ | **否**。`window` 关键字仅存在于 GEMMA `-calccor` 与 RRM 动态时间窗，非窗口遗传方差。 |
| G. 当前输出？ | **fixed-effect BETA/SE/P**（与 LMM 相同语义）。 |

**结论：CosmicScanner `--single-step` 是 H 核背景下的固定 SNP 效应 LMM 扫描（Method 1），与 BLUPF90/POSTGSF90 生态的 classic backsolve ssGWAS（Method 2）不是同一种方法。**

### 8.5 命名评估

最准确（按优先级）：
1. **H-matrix GWAS**（H-kernel LMM association）
2. **single-step LMM GWAS**
3. **single-step association scan**（保留 `--single-step` CLI 时代替不当语境下）
4. ~~ssGWAS~~：仅当文档明确其 = H 核扫描而非 backsolve 时可用，否则与经典名冲突

**保留 `--single-step` CLI（不做 breaking change）**；文档必须加算法定义句。拆分接口（**已实现**）：

```
--single-step-method lmm        （默认 == 现状，backward compatible）
--single-step-method backsolve  （Phase B experimental，已实现并通过解析/模拟验证）
```

## 9. Classic backsolve ssGWAS 数学设计（**Phase B 已实现，2026-08**；本节保留完整设计记录）

> 全部公式已联网核实（2026-08）：H⁻¹ 结构（Aguilar 2010; Legarra 2009; Christensen-Lund 2010）、backsolve（VanRaden 2008; Wang 2012; Strandén-Garrick 2009）、SNP P 值（Aguilar et al. 2019, GSE 51:28）、窗口方差与 BLUPF90 教程流程。细节与出处见 `docs/SSGWAS_METHOD_AUDIT.md` §15。

1. **ssGBLUP**：`y = Xb + Za + e`，`a ~ N(0, Hσ²a)`，MME：
   ```
   [X'X      X'Z    ] [b]   [X'y ]
   [Z'X  Z'H⁻¹Z + A₁  ] [û] = [Z'y]
   ```
2. **H⁻¹**：如 8.1（τ/ω 参数可含 λ 与 A22 刻度校正）。
3. **EBV/GEBV**：`û`（全体含未基因型动物）；基因型动物子集的 `û_g`。
4. **genotyped breeding values**：`û_g = Z ĝ + ε`（Z 为中心化剂量矩阵）。
5. **SNP back-solving**：`û = D·M′·G*⁻¹·â_g`（Wang et al. 2012 一般形式，G*=MDM′）；D=I 时 `ĝ = Z'(ZZ')⁻¹ û_g`（VanRaden 2008）；代码最廉价等价式 **`ĝ = Z′·G⁻¹·â_g`**（Strandén-Garrick 2009; Leite et al. 2024 Eq.1）。当 n<m 用 ridge/pseudo-inverse。
6. **allele-frequency scaling**：对每个 SNP 以 `M = X − 2p` 构建；G 与 backsolve 的 Z 定义必须**成对一致**（标准刻度 `diag(√{2p(1-p)})`）。CosmicScanner `libCosmicCore` 已有原语：`PlinkReader::multiply_Zt_v`（`plink_reader.cpp:355,411`）+ `snp_means`/`snp_inv_stds=1/√var`——注意其 Z=(x−mean)/std 是标准化刻度，输出需换算（audit §10.3）。
7. **SNP effect**：`ĝ`（加性替代数目刻度）。
8. **SNP weight**：`w_j = |ĝ_j|×√(2p_j(1-p_j))`（或 `ĝ_j²`、`|ĝ_j|` 变体）。
9. **window genetic variance**：`%Var_window = Σ_{j∈window} 2p_j(1-p_j) ĝ_j² / Σ_all 2p_j(1-p_j) ĝ_j²`（BLUPF90 教程默认窗口 = 10 个连续 SNP）。
10. **iterative weighted ssGBLUP（twssGBLUP）**：第 r 轮用 `w^{(r-1)}`（须 **`w ← w/mean(w)` 归一化**以保持 G* 方差尺度）重标度 G → 重建 G* → 再跑 ssGBLUP → 更新 ĝ/w，通常 1-3 轮（设计规范见 audit §9——存档不实现）。
11. **SNP P 值（可选）**：`p_i = P_t( û_i / √(σ_i²/n), n−1 )`（Aguilar et al. 2019, GSE 51:28；σ_i² = SNP 遗传方差，n = 有该 SNP 基因型的动物数；需 MME 系数阵逆的基因型块 C^{u₂u₂}，大群体用 APY 近似 Leite et al. 2024）。**t 型检验，非固定效应 Wald P**——backsolve 输出的 P 语义与此不同。

**可复用组件**：`A⁻¹`（buildAInverse）、`A22`（buildA22）、`G`（GrmBuilder，**含 `snp_weight_file` → D 加权零新增**）、`H⁻¹`（computeHInverse）、genotype matrix（PlinkReader/BgenReader/PgenReader）含 `multiply_Zt_v`/`snp_means`/`snp_inv_stds`。
**需新增**：MME 联合求解器（含未基因型动物）——**CosmicSolver 已有可移植实现**（`src/pipeline/solve/gs_solve.cpp`：SimplicialLDLT/MatrixFreePCG；`src/pipeline/io/matrix_builder.cpp:500-633`：Hinv/tuneG/blend/ridge/τω）；backsolve + 刻度换算；weights/window 归并；迭代循环（仅设计）。实施/复用表见 `docs/SSGWAS_METHOD_AUDIT.md` §11。

## 10. 模型语义矩阵（Section 15 交付）

| Model | BLUP/GS source | Random-effect background | SNP design | SNP parameter | Null hypothesis | Primary output | Association type | External reference | Status | Naming accurate? |
|---|---|---|---|---|---|---|---|---|---|---|
| Standard LMM | GBLUP 核（VCE only） | G·σ²g | x_i（个体） | β_j (scalar) | β_j=0 | BETA/SE/P | ordinary fixed SNP | GEMMA r>0.99 | 稳定/发布 | ✓ |
| Animal Model | A-kernel VCE | A22·σ²a | x_i | β_j | β_j=0 | BETA/SE/P + A22 核 | ordinary fixed SNP | GEMMA -k A22 核 | 已交付 | ✓（pedigree-kinship LMM GWAS） |
| Repeatability | 记录级 VCE | G (个体)+I(PE)+I(e) | x_i（按记录展开） | β_j | β_j=0 | BETA/SE/P | ordinary fixed SNP | RepeatABEL r>0.99 | 已交付 | ✓ |
| Common env | 记录级 VCE | G+I_g 组 | x_i | β_j | β_j=0 | BETA/SE/P | ordinary fixed SNP | 模拟 | 已交付 | ✓（shared-env adjusted） |
| **Maternal** | 多分量 VCE（direct+maternal） | G_d + G_m(=G_d)+残差 | direct 档：x_i（个体）；maternal/joint 档：x_dam / (x_i,x_dam)（`--maternal-snp-test`，PLINK-only, 已实现） | β_j (direct) / β_m (dam) / 2-df joint | β_j=0 / β_m=0 / β_d=β_m=0 | BETA/SE/P；`*.maternal.assoc.txt`/`*.joint.assoc.txt` | direct 档：ordinary fixed SNP（背景含 maternal RE）；maternal/joint 档：model-specific maternal SNP | sommer VCE（β 仅 direct） | 已交付 (direct/maternal/joint) | direct 档应称 maternal-adjusted GWAS；maternal/joint 档为 maternal/genotype SNP tests（simulation-validated） |
| **RRM balanced joint** | 观察级 RRM VCE | G⊗D + I⊗P + Ie | **x_i × Φ(t)** | γ_j (q+1 向量) | γ_j=0（q+1-df） | JOINT_P + BETA_Lk + β(t) 曲线 | **time-varying SNP** | GMA-fixreg chi/p 对齐 | 已交付 | ✓（RRM-GWAS/dynamic） |
| RRM coefficient lines / unbalanced Method B | BLUP 个体 Legendre 系数 | G（系数表性） | x_i（普通，无 ×basis） | β_j（多性状/逐基） | β_j=0 | JOINT_P/BEST 或逐基 | **derived/BLUP-coefficient SNP** | GMAT 常染色体 r=0.952/top 一致 | 已交付 | ✓（须注明 coefficient-based two-step） |
| **RRM unbalanced Method A (record)** | 记录级 RRM VCE | G⊗D + I⊗P + Ie | **x_{subject(r)} × Φ(t_r)**（记录级） | γ_j (q+1 向量) | γ_j=0（df=rank） | `*.rrm_unbal.assoc.txt`（JOINT_P + BETA_Lk + β(t) 曲线） | **time-varying SNP** | 模拟 A–G + 平衡等价（Pearson≥0.99999 级） | 已交付（simulation-validated，`--rrm-gwas-method record`） | ✓（record-level RRM-GWAS，须注明 experimental） |
| Multivariate | mv VCE | G⊗Σg | x_i（多性状） | β_j∈R^t | β_j=0 over t | mvWald/P | multi-trait joint SNP | GEMMA/EMmvGWAS r>0.99 | 已交付 | ✓ |
| **ssGWAS（Method 1, `--single-step`）** | H·σ²g VCE（仅基因型样本） | **H 子块**（H-kernel） | x_i | β_j | β_j=0 | BETA/SE/P | **H-kernel fixed SNP** | BLUPF90/POSTGSF90 **未比对** | 已交付 | **✗** — 与经典 backsolve ssGWAS 同名不同法；文档改称 H-matrix GWAS / single-step LMM GWAS（audit 见 `SSGWAS_METHOD_AUDIT.md`） |
| **ssGWAS backsolve（Method 2, 经典）** | ssGBLUP MME（**含未基因型动物**） | H·σ²a（全模型） | back-solve `ĝ = Z′G⁻¹â_g`（= D=M′G*⁻¹â, D=I） | û_j（SNP 效应） | —（可选 Aguilar 2019 t-P，未实现） | SNP effects / Var_SNP / %Var_window | **marker-effect（back-solving）** | Wang 2012; Aguilar 2019; Strandén-Garrick 2009（**已联网核实**） | **已交付（`--single-step-method backsolve`；解析门 EBV diff 5e-7 + 模拟门 Top10=5/5 + 系谱洗牌逐位一致；见 SSGWAS_METHOD_AUDIT §16）** | —（Marker-effect 路径，与 Method 1 并列、互不混淆） |
| **ssGWAS APY-exact P-value** | ssGBLUP MME + APY `G_APY⁻¹` | H·σ²a（APY G inverse） | backsolve + exact core selected inverse | SNP effect / SE / Z / raw P | Aguilar-style SNP effect P | `*.ss_backsolve.pvalue.txt` | scalable APY uncertainty path | D gates C2-G3/D3/D4/D5/D6.5/D7 | **Experimental validated** | validated scalable uncertainty path; raw P, no GC |
| **ssGWAS APY-approx P-value** | ssGBLUP MME + APY `G_APY⁻¹` | Eq.8/Eq.9 approximate core uncertainty | core-only APY effect + approximate uncertainty | SNP effect / approximate SE / Z / raw P | approximate SNP effect P | `*.ss_backsolve.pvalue.txt` | approximate scalable uncertainty layer | D6/D6.5 calibration boundary + D7 ranking boundary | **Experimental approximate** | `STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`; not calibrated P-value |

## 11. 结论摘要（三个 Q 的直接回答）

- **Q1.** `maternal_model.cpp` 默认 direct 路径是 **maternal-adjusted GWAS（Type A）**：SNP 设计为子代个体 genotype、检验 `β_direct,j=0`、输出 BETA/SE/P；`σ_adm`/`G_m`（=G_d）饱和度低。**自 `--maternal-snp-test maternal|joint` 起**（§5.3），可检验 maternal SNP effect（`x_dam·β_m,j`）与 direct+maternal 2-df joint，属 Model-specific GWAS（Type B），simulation-validated（10 场景 PASS，λ_GC≈1）。文献层面引用时须区分 direct（Type A）与 maternal/joint（Type B）两个语义。
- **Q2.** **平衡 RRM 主路径（观察级 joint 扫描）与 GMA-fixreg 路径是真正的 time-varying SNP effect GWAS**：SNP 设计为 `x_j×Φ(t)`、联合检验 `H0: γ_j=0 (df=q+1)`、由 γ 重构 `β_j(t)` 并输出二次型 SE。**非平衡路径是 BLUP 系数空间多性状扫描**（`x_j` 对 Legendre 系数表型的普通联合 SNP 检验），不是 `x⊗Φ` 固定效应检验——两者都合法，但方法语义不同，文档须分别说明。
- **Q3.** **不同，且现在两种方法都可用：** CosmicScanner `--single-step`（默认）= Method 1（H 核 LMM 固定 SNP 扫描，输出 fixed-effect BETA/SE/P）；**`--single-step-method backsolve` = Method 2**（ssGBLUP MME → GEBV → `ĝ=Z′G⁻¹â_g` → SNP effects/variances/window %var）。差异要点：①Method 2 对随机遗传效应求评价解并反解（含未基因型动物）②产出 marker effect/SNP variance 型输出 ③迭代加权（twssGBLUP）与 Aguilar 2019 t-P **未实现**。两路径互不混淆（输出文件前缀 `.ss_backsolve.` 与 `.assoc.` 区分）。

## 12. Audit 边界与待办

- 本轮 **Phase 3 已实现** maternal `--maternal-snp-test direct|maternal|joint`（`scanMaternalSnp()`；PLINK-only；10 场景模拟验证全 PASS；默认 `direct` 与旧 `--maternal-model` 输出 byte-identical）。**仍待办**：maternal/joint 在 BGEN/PGEN 下的扫描。
- **Phase B backsolve 已实现（2026-08）**：`--single-step-method backsolve` 全管线（MME LDLT/PCG + G⁻¹ backsolve + SNP/window 方差 + `*.ss_backsolve.*` 输出），解析门 + 模拟门 + 回归网全 PASS（`SSGWAS_METHOD_AUDIT.md` §16）；顺带修复 `libCosmicCore` GRM 构建的 `std::vector<bool>` OpenMP 竞态（G 确定性）。**外部基准（2026-08-17，HIBLUP v1.6.0）已完成**：EBV 逐类一致（Pearson 1.00000000 / slope 1.000028，`SSGWAS_METHOD_AUDIT.md` §18 + `docs/SSGBLUP_EXTERNAL_VALIDATION.md`）；参考对齐通道 = `--ssgwas-tune-g --ssgwas-blend-a22`（tuneG 先行、混合目标 A22、SNP 反解用调优未混合 G*⁻¹，与 HIBLUP 语义一致；默认路径不变仍为单位阵混合）。**仍待办**：BLUPF90/postGSf90 外部基准（就绪包 `tests/ssgwas_external/blupf90/`，SKIP-safe）、Aguilar 2019 t-P（Phase C）、weighted/iterative（Phase D）。
- **RRM 非平衡 Method A（`--rrm-gwas-method record`）已实现**（记录级 x⊗Φ 联合线，§6.2），A–G 模拟套件 + 平衡等价 + A/B 对比全 PASS；Method B（默认 coefficient）GMAT 常染色体对齐通过（§4.2）。仍未实现/验证：BGEN/PGEN 下 record Method A；X 染色体默认排除选项（Method B 与 GMAT 全 GC 对等）。
- **联网恢复（2026-08）**：ssGWAS canonical 定义已复核——H⁻¹（Aguilar 2010; Legarra 2009; Christensen-Lund 2010）、backsolve（VanRaden 2008; Wang 2012; Strandén-Garrick 2009; Leite 2024）、SNP P（Aguilar 2019, GSE 51:28）、BLUPF90 教程流程（Masuda），见 `docs/SSGWAS_METHOD_AUDIT.md` §15。GMA/GMAT 原文仍未复核（RRM 语境）。
- 下一阶段建议（按风险排序）：
  1. 文档/命名修复（本文档已含；把 Method A/B 措辞同步到 README/Tutorial/Manuscript——低风险）。
  2. Method A record 的正式外部参考（GMAT 无观察级 x⊗Φ 输出线；可做 GMA unbalance 对齐复核）——中风险。
  3. X 染色体排除选项（`--exclude X` / GRM 自动排除性染色体）使 Method B/GMAT 全 SNP 集对等——中风险。
  4. classic backsolve ssGWAS（Method 2，§9 + `docs/SSGWAS_METHOD_AUDIT.md` 完整设计）作为 experimental 新路径——**已完成（Phase B）**：无权重 backsolve + SNP/窗口方差（不默认 weighted/iterative D，无 Aguilar 2019 P）。采用 CosmicScanner 原生 MME 求解（SimplicialLDLT / MatrixFreePCG）而非移植；验证 = 解析门（EBV 5e-7）+ 模拟门（Top10=5/5、系谱洗牌逐位一致）+ 回归网 + **HIBLUP 外部基准（EBV 逐类 Pearson 1.00000000，2026-08-17）**。仍待办：BLUPF90/postGSf90 基准（SKIP-safe 就绪包已建）、窗口方差外部对齐（当前为内部定义 A）。
