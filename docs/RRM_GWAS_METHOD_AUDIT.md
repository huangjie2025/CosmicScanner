# RRM-GWAS Method Audit: Record-Level Method A vs Coefficient-Based Method B

> 本文档按"五层审计"框架（BLUP/GS source → GWAS background covariance model → SNP-effect design →
> association test → primary output）拆解 CosmicScanner 非平衡 Random Regression（RRM）GWAS 的两条扫描路径，
> 并给出 **Association type** 与 **Naming accurate?** 判定，以及全部模拟/外部对齐证据链的当前数值。

---

## 0. 结论摘要（TL;DR）

- **Method B（`--rrm-gwas-method coefficient`，默认，兼容旧行为；`--rrm-compatible gmat-unbalance` 为 GMAT 兼容口径）**：先做个体内 Legendre
  系数 BLUP（multivariate LMM，`Vg⊗A + Ve⊗I`），再把系数矩阵当作多性状表型做 Kronecker
  多性状 Wald 扫描。SNP 设计 = 普通 genotype（`x_i`，逐"系数性状"），不是 `x⊗Φ`。与 GMAT
  `longwas --unbalance` 对齐：修正 `Sex_class` fixed-effect 进入 coefficient null 后，**全 SNP 集
  Pearson r(-log10P)=0.9968、Spearman=0.9959，autosomes-only Pearson=0.9969、Spearman=0.9964，
  top SNP=UNC18848064 且 top10 全一致**。旧全 SNP 集 Pearson=0.9470 主要来自时间依赖分类协变量
  未进入 coefficient null 的 subject-level fixed design 和分类基线编码错误；不是隐藏 imputation mismatch。
- **Method A（`--rrm-gwas-method record`，新增主扫描）**：记录级固定效应联合检验，
  D(r,k) = genotype(subject(r)) × φ_k(time_r)，得分/信息量经 eigen-截断（M⁺，df=rank）做
  `score'M⁺score ~ χ²(df)` Wald/Score 检验，**真正的 time-varying SNP effect 设计**（Type B）。
  平衡路径（`runRrmScan`）等价于隐式 x⊗Φ Method A —— `tests/run_rrm_balance_record_equivalence.sh`
  实证 N=200/M=300/T=8 上 BETA/SE/STAT/log10P Pearson ≥ 0.99999、STAT 绝对差 ≤ 0.084。
- **A 与 B 是两种不同检验**：A 直接捕获 β(t)=φ'γ（注入轨迹 causal 秩=1、p 最小可达 e-18）；
  B 在系数空间对平均轨迹方向更敏感（同数据 causal rank=14）。`run_rrm_method_ab_comparison.sh`
  给出两者 log10P 的 Pearson=0.385/Spearman=0.278 —— **不矛盾**：自由度、null fit、SNP 效应分组不同。
- **非平衡用户选择**：想严格 time-varying SNP 效应 → `--rrm-gwas-method record`（模拟验证、power 更强）；
  想与 GMAT/GMA 生态兼容、或与历史 `--rrm-unbalance` 结果衔接 → `--rrm-gwas-method coefficient`（默认）。
- **状态认定**：Method A record = **experimental / simulation-validated**（无外部直接参考，按 §6 定位）；
  Method B coefficient = **compatible / GMAT aligned**。两者均可用于非平衡数据，但语义不同、文档须分开解释。

---

## 1. 五层审计：Method A（record-level, `--rrm-gwas-method record`）

### 1.1 分层定义

| 层 | Method A（record） |
|---|---|
| **BLUP/GS source** | 无 GS/BLUP 输出；null 只做 VCE（covariance 估计）+ 记录级 score/info 缓存 |
| **GWAS background covariance model** | 记录级 RRM null：`V = (G⊗D) + (I⊗P) + σ²_e·I`，`D=Φ'KΦ`、`P=Φ'P0Φ`（order-3 Legendre/B-spline 基），AI-REML `--rrm-ai-iter`；LOCO 按染色体重建（`--no-loco` 关闭） |
| **SNP-effect design** | **`D(r,k) = x_{subject(r),j} · φ_k(t_r)`**（record×basis 设计列）；X_base 为个体级基协变量，随记录展开；可选 `--tclass/--tcovar` 时加时间依赖协变量列 |
| **Association test** | 联合 Score/Wald：`score = D'P_reml·y`，`M = D'V⁻¹D`（eigen-分解，截断 |λ|≤1e-10·λ_max，`df=rank≥1`），`stat = score'M⁺score`，`p=χ²_df`；`STATUS ∈ {OK, RANK_DEFICIENT, ILL_CONDITIONED, NUMERICAL_FAILURE}`（不静默丢 SNP，rank 缺口显式记录） |
| **Primary output** | `*.rrm_unbal.assoc.txt`（`JOINT_CHISQ JOINT_P` + `BETA_Lk SE_Lk`，df 列时按 rank）；`*.rrm_unbal.null.txt`；top SNP 动态曲线 `*.rrm_unbal.dynamic_curve.txt` / `.dynamic_summary.txt` |

### 1.2 统计细节与实现锚点

- 记录级 design：`src/gwas_rrm_scan.cpp` `fit_rrm_record_snp`（score/info 聚合由观测聚合到个体后到 SNP
  level，等价于对 `x⊗Φ` 设计做积分）。`cosmic::buildRrmBasisMatrix`
  （`libCosmicCore/src/design.cpp/.h`）统一随机与 SNP 设计的 Legendre（shift-[-1,1] 归一化
  `sqrt((2k+1)/2)`）与 B-spline（degree≤3 + partition-of-unity）基。
- 每个 scope（`--no-loco` 全局或逐染色体 LOCO）只解一次 null 并缓存 `P_reml`/`V⁻¹`/eigen；逐 SNP
  复用缓存，不重建 n_record×n_record 矩阵（N=1297 obs 时 null≈31 s、scan≈0.09 s）。
- 平衡路径 `runRrmScan` 的观测级 joint 扫描是 Method A 的隐式等价形式（x⊗Φ 设计 + score/info 积分），
  已在 `run_rrm_balance_record_equivalence.sh` 中与记录显式路径对齐（§5.1）。

### 1.3 Naming

**time-varying SNP effect GWAS / record-level RRM-GWAS（Type B）** ✓（与平衡路径一致）。

---

## 2. 五层审计：Method B（coefficient-based, `--rrm-gwas-method coefficient`，默认）

### 2.1 分层定义

| 层 | Method B（coefficient） |
|---|---|
| **BLUP/GS source** | **个体 Legendre 系数 BLUP**（multivariate AI-REML：`Vg⊗A + Ve⊗I`，A=个体 GRM；`coef_mat` 作为多性状表型）——这里是 BLUP/系数来源 |
| **GWAS background covariance model** | 系数空间 Kronecker mvLMM：`V = Vg⊗A + Ve⊗I`，`MvEigenVce` AI-REML `--rrm-ai-iter` |
| **SNP-effect design** | 普通 genotype `x_i`（对每个"系数性状" L0…Lq），无 ×basis |
| **Association test** | 多性状 Wald `H0: β_L0=…=β_Lq=0`（df=4，order3），`JOINT_CHISQ/JOINT_P`；t-df SE/曲线为 proxy |
| **Primary output** | `*.rrm_unbal.assoc.txt`（`BETA_Lk SE_Lk` 为系数空间 SNP 效应）、`*.rrm_unbal.null.txt`、动态曲线（由多性状 Wald 协方差代理 `β(t)=φ'β_snp`） |

### 2.2 Naming

**coefficient-based (BLUP) two-step GWAS** —— 合法但必须按此名称解释；不是 `x⊗Φ` 固定效应动态 GWAS。
`--rrm-unbalance` 单独使用 ≡ 默认 `--rrm-gwas-method coefficient`（backward compatible，已验证）。

---

## 3. 隶属框架

| Model | Association type | 状态 |
|---|---|---|
| Method A record（非平衡），平衡 joint，GMA-fixreg | **Type B**（time-varying / x⊗Φ 模型特定检验） | Method A：simulation-validated；平衡 joint：外部 GMA-fixreg 对齐；GMA-fixreg：r=1.0 |
| Method B coefficient（非平衡默认）| **Type C**（derived/BLUP-coefficient 检验） | GMAT 对齐（§4） |

---

## 4. 外部对齐证据

### 4.1 GMA-compatible balanced fixreg（`tests/run_gma_compatible_fixreg_alignment.sh`）

- Dataset：GMA official balanced case subset N=120, SNP=20, TIME=16。
- GMA reference = `cal_kin` + `balance_varcom` + `balance_longwas_fixreg`（py2→py3 本地重建，
  `/tmp/opencode/gma_py3{_deps,_install,tests}`，deps：pysnptools/bed-reader/patsy）。
- 结果：`Pearson_log10P=1.0`、`Spearman_P=1.0`、`MaxAbsDeltaP≈2.73e-9`、
  `MaxAbsDeltaChi≈1.26e-7`，top SNP `JAX00241881` 双一致（GMA wrapper 8.7 s/161 MB；
  CS gma-fixreg 18.7 s/27 MB）。
- 归属：该 benchmark 比较的是 **SNP × 基联合固定效应检验**（Type B），不是 BLUP 系数。

### 4.2 GMAT unbalanced（`tests/run_gmat_unbalanced_alignment.sh`）

- Dataset：GMAT `mouse_long`，N=1212，SNP=11833，obs=19392，TIME=16（`--rrm-order 3`，
  时间依赖分类协变量 `Sex_class`，`--rrm-ai-iter 20 --no-loco --rrm-compatible gmat-unbalance`）。
- Method B（coefficient）对齐 GMAT `longwas --unbalance`（GMAT reference `.res` 从
  `gmat-1.01.rar` 重新生成并验证，GMAT 侧 top `UNC18848064` p=4.1846e-6 与文献 4.2e-6 一致）：

| 范围 | Pearson_log10P | Spearman_P | Top_GMAT | Top_CS | Top10 overlap |
|---|---:|---:|---|---:|---:|
| 全部 11833 | **0.996763** | **0.995904** | **UNC18848064** | **UNC18848064** | **10** |
| **Autosomes** | **0.996945** | **0.996390** | **UNC18848064** | **UNC18848064** | **10** |

  旧版本全 SNP 集 Pearson=0.9470 / Spearman=0.9393，常染色体 Pearson=0.9521 / Spearman=0.9401。
  主因不是 X 染色体或隐藏 imputation mismatch，而是 coefficient path 没有把 `Sex_class` 聚合进
  subject-level fixed design，且分类水平编码基线错误。当前 `--rrm-compatible gmat-unbalance` 会强制使用
  coefficient path，并把记录级 `--tclass/--tcovar` 聚合到 subject-level fixed effects；对 `Sex_class` 等
  个体恒定协变量，这等价于 GMAT 的 sex × Legendre coefficient fixed-effect 口径。
- 验证命令：`GMAT_DATA=/tmp/opencode/gmat_data bash tests/run_gmat_unbalanced_alignment.sh` →
  `Pass=True`。
- 归属：该 benchmark 对应 **coefficient/BLUP-derived 多性状 SNP 检验**（Method B），不是观察级 x⊗Φ。

### 4.3 旧版 X 染色体差异注记（已降级）

- 旧现象：Method B 与 GMAT 在 307 个 chr23(X) SNP 上的偏差系统性高于常染色体；单个 X 位点
  `JAX00177214` 成为 CS 全集 top（p≈2e-5）而 GMAT p=0.075。
- 当前解释：修正 fixed-effect 进入 coefficient null 后，X 位点不再主导全基因组 top 差异，全集与常染色体
  均达到 r≈0.996。X 染色体 GRM 处理仍是可能的细节差异来源，但不再是当前 GMAT-compatible path 的
  主导误差项。
- 待办：若需要逐字节复刻 GMAT，可继续实现可选的 autosome-only GRM 构建；Method A 侧暂无外部参考（见 §6.3）。

---

## 5. 内部一致性 / 模拟验证证据

### 5.1 平衡 ⇔ 记录 Method A 等价（`tests/run_rrm_balance_record_equivalence.sh`）

- 自包含模拟：N=200, M=300, T=8，seed 20260815，无固定遗传、纯残差 null。
- 对比：平衡显式 joint 路径（运行于平衡记录）vs 非平衡 record Method A（运行于同一记录展开、`--no-loco`）。
- 结果（PASS，30 万级 SNP× 统计量）：`BETA_L0..3` Pearson ∈ [0.999993, 0.99999996]、
  `SE` ∈ [0.999999, 0.9999999]、`STAT_JOINT`=0.999999405、`log10P`=0.999999655、
  `MaxAbsDiff_BETA`=2.52e-4、`MaxAbsDiff_STAT`=0.084、RANK/DF={4}、STATUS 全 OK（300/300）。
- 门限（论文诚实边界）：MIN_PEARSON=0.99999、MAX_ABS_DIFF_BETA=1e-3、MAX_ABS_DIFF_STAT=0.5 ——
  残余差异来自两路径 VC 初始化方式不同（非平衡 record 线应用了偏置校正），AI_ITER 10..50 不消除，
  已写入脚本注释，非静默放宽。

### 5.2 Method A 模拟套件 A–G（`tests/run_rrm_snp_trajectory_simulation.py` + `run_test_rrm_method_a_simulation.sh`）

| 场景 | 内容 | 判定标准 | 结果 |
|---|---|---|---|
| A | null，某 SNP 弱信号 | FPR(0.05)≈0.05、inflation≈1 | FPR=0.01、λ=0.9238 |
| B | 恒定因果轨迹 rank1 | causal p 显著、flatness 轨迹 | causal p=6.8e-11、flat=0.371 |
| C | 线性轨迹 | causal rank1、Pearson 轨迹相关 | rank=1、p=1.04e-7、traj corr=0.9968 |
| D | 二次轨迹 | rank1、轨迹相关 | rank=1、p=4.16e-14、corr=0.9987 |
| E | sign-switch 轨迹 | rank1、轨迹相关 | rank=1、p=5.82e-18、corr=0.9952 |
| F | 恒定因果（残差背景） | rank1、flat | rank=1、p=3.06e-4、flat=0.509 |
| G | rank 缺陷/共线诊断（diag-only） | 无静默丢弃、STATUS 显式 | 296 RANK_DEFICIENT + 4 ILL_CONDITIONED，P 全 finite，0 NUMERICAL_FAILURE |

结果目录：`tests/results/rrm_method_a_simulation/`。

### 5.3 Method A vs Method B 同数据对比（`tests/run_rrm_method_ab_comparison.sh`）

- 数据：注入 trajectory 因果（records=1297, miss=0.2），同 null fit。
- 结果：causal rank A=1 vs B=14；causal p_A=6.66e-8（区域） vs p_B 弱——Method A 对
  time-varying 效应有直接 power；`Pearson_log10P(A,B)=0.385`、`Spearman=0.278`、
  `MaxAbsDiff_log10P=6.05`、DF_A={4}、STATUS 全 OK。Runtime A=31.36 s（null 31.17 s/scan 0.09 s,
  628 MB）vs B=0.14 s/32 MB。结论：二者皆可用、A 对随时间变化效应更敏感，文档须按不同检验解释。

### 5.4 单元 / 构建回归

- `tests/test_rrm_basis_identity.cpp`（`cosmic::buildRrmBasisMatrix` 命名空间修复后）：
  `max|leg-analytic|=8.88178e-16`、`max|bs-partition-1|=3.33067e-16`、`reproducibility=0`。
- `ctest` 11/11 通过；`--help` 同时展示 `--rrm-gwas-method record|coefficient`；
  `tests/run_rrm_dynamic_simulation.sh --loco-smoke` 通过。

---

## 6. 方法语义判定表（五层严格判定）

| 问题 | Method A（record） | Method B（coefficient，默认） |
|---|---|---|
| A. SNP 固定设计是 `x_j×basis(t)`？ | **是**（记录级 D(r,k)=x·φ） | 否（普通 x_i 对系数性状） |
| B. 估计 γ_j0…γ_jq？ | 是（BETA_Lk=γ_k，二次型 SE） | 否（输出的是系数空间 β_Lk） |
| C. β_j(t) 由 SNP-specific γ 重构？ | 是（严格） | 代理（多性状 Wald 协方差 → β(t)） |
| D. joint 检验 H0: γ_j=0？ | 是（score'M⁺score，df=rank） | 否（H0: β_L0=…=β_Lq=0，多性状 Wald） |
| E. 与平衡 RRM joint/GMA-fixreg 同构？ | 是（平衡 joint 为其隐式等价，§5.1） | 否 |
| F. GMAT benchmark 对应哪个？ | 不适用（无直接参考） | 对应 BLUP 系数多性状检验（r=0.95 级） |
| G. 推荐语义名 | **record-level time-varying SNP RRM-GWAS** | **coefficient-based (BLUP) two-step GWAS** |

---

## 7. 状态 / 命名认定

- `--rrm-gwas-method record`：**Experimental / simulation-validated**（A–G 套件 + 平衡等价；
  power 直接、df=rank 诊断显式）。论文/文档须标注，不伪称 production 级外部验证。
- `--rrm-gwas-method coefficient`（默认，`--rrm-unbalance` 兼容）：**GMAT-aligned /
  compatible**（§4.2，常染色体 r=0.952/top 一致）。
- 方法名区分：非平衡两侧分别叫 **record-level RRM-GWAS** 与 **coefficient-based RRM-GWAS**；
  不得把 Method B 输出伪装成 x⊗Φ 固定效应动态 GWAS。

## 8. Audit 边界与待办

- 已交付：Method A record 扫描、A/G 模拟套件、AB 对比、平衡等价、X 染色体内方法学注记、GMA/GMAT
  外部对齐复跑（fixture 再生 + py3 转换重建）。
- 仍未实现/验证：Method A 的正式外部参考（GMAT 无观察级 x⊗Φ 输出线；GMA unbalance 为
  coefficient 型）；BGEN/PGEN 下 record Method A；X 染色体默认排除选项（Method B/GMAT 全对等）。
- 联网可用后应重核：GMAT 原文的 GRM 是否确为 autosomal-only、GMA unbalance 是否也是 coefficient
  两阶段——本文档标注 **external method-definition verification pending**，结论来自本地代码与
  参考工具运行证据。
