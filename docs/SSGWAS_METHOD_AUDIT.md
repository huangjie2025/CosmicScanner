# ssGWAS Method Audit: H-Kernel Association vs Canonical ssGBLUP Back-solving

> Status: **设计审计 + 已实施（2026-08）**。Phase B（unweighted backsolve + SNP/window variance）已落地并通过解析/模拟验证（§16）；本文档保留完整审计与设计记录。
> Scope: CosmicScanner `--single-step` 现状 vs BLUPF90/postGSf90 生态的经典 back-solving ssGWAS（Wang et al. 2012; Aguilar et al. 2019）。
> 证据来源：代码逐行追踪（`single_step_model.cpp` / `pedigree_relationship.cpp` / `grm_builder.cpp` / `libCosmicCore/src/plink_reader.cpp`）+ 联网核实的文献公式（§15）。

---

## 0. 结论摘要（TL;DR）

1. **`--single-step` = Method 1（H-kernel fixed-SNP LMM association）**，不是经典 backsolve ssGWAS（Method 2）。全部严格判定见 `docs/GWAS_MODEL_DEFINITIONS.md` §8.4。
2. **Method 2 的全部数学公式已联网核实**（不再标注 external verification pending）：H⁻¹ 结构（Aguilar 2010 / Legarra 2009 / Christensen-Lund 2010）、backsolve（VanRaden 2008; Strandén-Garrick 2009; Wang 2012）、SNP 效应方差与 P 值（Aguilar et al. 2019，GSE 51:28）、窗口方差（Fragomeni 2014; Wang 2014）。见 §15。
3. **审计时的差距（已全部由 Phase B 实施关闭，见 §16）**：无 MME 联合求解（含未基因型动物）、无 GEBV 输出（`writeEbv()` 是 stub）、无 SNP back-solving、无窗口方差。§1 的对照表保留历史状态。
4. **CosmicSolver（同仓库族）提供 Method 2 骨架**：ssGBLUP MME 求解（SimplicialLDLT / MatrixFreePCG）+ `--snp-effect` back-solving + HIBLUP 基准（`tests/run_test_ssgblup_smoke.sh`）。Phase B 的主要移植来源（§11）。
5. **Phase B 建议（见 §11.4）**：B（unweighted backsolve + SNP/window variance）或 C（+ Aguilar 2019 P 值）。**已按 B 完整实施并通过验证**（§16）；C（P 值）与 D（weighted/iterative）未实现。

---

## 1. 术语与命名（Terminology & Naming）

| 名称 | 定义 | 示例/引用 |
|---|---|---|
| **Method 1：H-kernel mixed-model association**（= 现状） | 以 H（或其基因型子块）为随机核的 LMM，逐 SNP 对 `x_j` 做固定效应检验（Wald/score），输出 BETA/SE/P | CosmicScanner `--single-step`；GEMMA `-lmm` 的 H 核变体 |
| **Method 2：classic backsolve ssGWAS**（= 经典） | ssGBLUP（MME 求解含未基因型动物）→ GEBV → 线性变换反解 SNP 效应 →（可选）SNP P / 窗口方差 | BLUPF90 + postGSf90；Wang et al. 2012; Aguilar et al. 2019 |
| **Method 3：weighted / iterative ssGWAS**（= wssGBLUP） | Method 2 + SNP 权重 D（重标度 G* = MDM′）+ 迭代更新 | Wang et al. 2012; ZHANG et al. 2016 综述 |

命名结论（与 `GWAS_MODEL_DEFINITIONS.md` §8.5 一致）：
- 保留 `--single-step` CLI（breaking change 不值得）；文档层面本方法写作 **"single-step LMM GWAS（H-matrix GWAS / H-kernel association）"**。
- backsolve 路径经 **`--single-step-method lmm|backsolve`** 提供（**已实现**，lmm=现状默认，backsolve=Phase B；见 §16）。

### 1.1 CURRENT（CosmicScanner）vs CANONICAL（BLUPF90/postGSf90 生态）对照表

| 特性 | CosmicScanner `--single-step`（现状） | 经典 ssGWAS 流程（canonical） | 相同？ | 需要的改动 |
|---|---|---|---|---|
| 亲缘矩阵 | H⁻¹ = A⁻¹ + diag(τG⁻¹ − ωA22⁻¹)，G 默认 VanRaden 2p(1−p)，可选 α·I 混合 | 同结构；preGSf90 另有 G−A22 对齐校正 `G=bG+a`、α 混合、τ≥1 | **≈（缺对齐）** | （可选）加 `tuneG` 对齐，跨软件对比前必做 |
| VCE | 只对基因型子块 (y_sub, X_sub, H_sub) REML（`ng == keep.n_keep`） | AI-REML/EM 全模型（含未基因型动物对 Vg/Ve 的信息） | **✗（仅基因型）** | Method 2 需要 MME 全模型或等价 H-BLUP |
| EBV/GEBV 求解 | 无（`writeEbv()` stub；无 MME） | MME `[X′X X′Z; Z′X Z′Z+H⁻¹λ]` 解得全部 EBV；基因型动物 GEBV = â_g | **✗** | **需新增 MME 求解器**（可移植 CosmicSolver） |
| SNP back-solving | 无 | â = D·M′·G*⁻¹·â_g（D=I 时 â = M′(MM′)⁻¹â_g）；等价 â = Z′G⁻¹â_g | **✗** | **需新增**（Z′v 原语已在 libCosmicCore） |
| SNP 效应 P 值 | 无（backsolve 不适用） | Aguilar 2019：p_i = P_t(û_i/√(σ_i²/n), n−1)（t 分布，n−1 df） | **✗** | 需新增（O(n³) PEV；可先跳过） |
| SNP / 窗口遗传方差 | 无 | %Var_window = Σ 2p(1−p)û² / total（滑动窗口，如 10 SNP） | **✗** | 需新增（廉价，低风险） |
| SNP weights（D） | GRM builder 已支持 `snp_weight_file`（w=√w 作用于 Z 列） | G* = MDM′ 重标度 | **≈（原语已有）** | 需组装（Phase B/C） |
| 迭代加权（twssGBLUP） | 无 | 1–3 轮 G*←MDM′→refit→更新 û/w（权重均值归一化） | **✗** | **设计完成，不实现（用户限制）** |
| 随机效应检验语义 | 逐 SNP 固定效应 Wald/score，BETA/SE/P | 无固定效应检验；输出 marker effect +（可选）t 型 P | **不同语义（两者都合法）** | 文档明确区分 |

---

## 2. 当前实现逐行追踪（CURRENT，只读代码确认）

源码：`src/models/single_step_model.cpp`（535 行）、`single_step_model.h`、`pedigree_relationship.{h,cpp}`、`libCosmicCore/src/grm_builder.cpp`。

### 2.1 数据流

```
输入：PLINK/BGEN/PGEN + 系谱 + --genotyped-list [--write-agh --no-loco ...]
  │
  ├─ PedigreeRelationship（ID SIRE DAM；NA/. → 0）→ 拓扑排序 → A⁻¹（Henderson 算法）
  │     └─ A22：由 A⁻¹ 的未基因型块做 Schur 补（不物化稠密 A）→ 岭求逆 A22⁻¹
  │
  ├─ GrmBuilder（算法 1=VanRaden 默认，2p(1−p) 缩放；--grm-alg 1..6）
  │     └─ G 可选混合：G ← (1−α)G + αI（α=--blend-g，默认 0.05）
  │
  ├─ H⁻¹ ← A⁻¹ + 基因型块 [τ·G⁻¹ − ω·A22⁻¹]；对称化 + max-abs 对称性检查（写 *.agh.meta.txt）
  │
  ├─ VCE（runVCE）：H_sub = buildHSubForG(G_)（重建 Hinv→全 H 求逆→取基因型子块）
  │     → EigenVce REML 于 (y_sub, X_sub, H_sub)；σ²g, σ²e —— 只用基因型样本（ng == keep.n_keep）
  │
  └─ 扫描（scanGwas）：σ²g/σ²e 固定 → EigenVce.fixedComponents → LmmEngine.prepareFromVce
        → 逐 SNP（个体自身 genotype）：Wald/score → *.assoc.txt（BETA SE P LOG10P）
        LOCO：computeSplit 分染色体 GRM → 每染色体 G_loco → buildHSubForG(blend(G_loco))
```

### 2.2 代码锚点

| 步骤 | 锚点 |
|---|---|
| G 混合 | `single_step_model.cpp:189-193`（`(1−α)G + αI`） |
| A22⁻¹ | `single_step_model.cpp:196-212`；Schur 补 `pedigree_relationship.cpp:193-207`；岭求逆 `pedigree_relationship.cpp:219-236`（ridge 从 1e-8 几何增长） |
| H⁻¹ | `single_step_model.cpp:214-238`（`H_inv_ = A_inv_; += τG⁻¹ − ωA22⁻¹`；对称化 231） |
| H 基因型子块 | `buildHSubForG` `single_step_model.cpp:240-263`（重建 Hinv → 全 H 岭求逆 → 取基因型子块） |
| VCE | `runVCE` `single_step_model.cpp:265-299`（`ng != keep_.n_keep` 直接 throw） |
| 扫描 | `scanGwas` `single_step_model.cpp:301-449`；LOCO 块 424-449（`S_loco = S_global − S_chr` 重标定） |
| GEBV 输出 | `writeEbv()` stub（打印 "GEBV output will be implemented in v1.1."） |
| A⁻¹ 构造 | `libCosmicCore/src/a_matrix.cpp:446` `invertA_henderson_fixed`（含 phantom parents、ml/dii 校正） |

### 2.3 输出物

- `*.assoc.txt`：`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`（与标准 LMM **完全同 schema**，fixed-effect 语义）
- `*_ssgwas_null.txt`；可选 `*.Ainv.bin` / `*.G.bin` / `*.Hinv.bin`
- 无 EBV 文件、无 SNP effect 文件、无权重文件、无窗口文件

---

## 3. H⁻¹ 构造审计（当前 vs 经典）

### 3.1 公式

```
H⁻¹ = A⁻¹ + [ 0     0
              0  τ·G⁻¹ − ω·A22⁻¹ ]
```

- 代码默认：τ = `h_tau` = 1.0，ω = `h_omega` = 1.0，α = `blend_g` = 0.05（与 BLUPF90 家族默认 τ=1、ω=1、α_blend≈0.05 一致）。
- BLUPF90 家族：τ ≥ 1 常被建议（抵消"基因型已知"造成的双重计数高估，Christensen-Lund 2010）；preGSf90 默认 G 用 VanRaden 2p(1−p)，并可 `OPTION tuneG` 做 G=b·G+a 对齐（目标 A22），`OPTION alpha_blend` 设 α。
- **CosmicScanner 缺失**：G−A22 对齐校正（`tuneG`）。在 α、τ=1、同一频率源下差距通常很小，但**跨软件数值对比前必须统一**（用户 §21 前提：不得在未对齐时指控 BLUPF90 不一致）。

### 3.2 构造参数对照表

| 参数 | CosmicScanner | BLUPF90 家族 | 默认一致？ | 影响 |
|---|---|---|---|---|
| τ (tau) | 1.0 | `OPTION tau`，默认 1.0 | ✓ | 大：双重计数校正 |
| ω (omega) | 1.0 | `OPTION omega`，默认 1.0 | ✓ | 中 |
| α blend (G+αI / (1−α)G+αA22) | 0.05 与 **I** 混合 | 0.05 与 **A22** 混合（alpha_blend） | 部分（混合对象不同） | 中；与 A22 混合更贴近经典 |
| G 缩放 | 算法 1：VanRaden 2p(1−p)（--grm-alg 1..6） | 默认 2p(1−p)（method 可改） | ✓ | 大（频率源须一致：基因型样本频率） |
| G−A22 对齐 | **无** | `OPTION tuneG`（G=bG+a） | ✗ | 中；对比前必须启用 |
| 求逆稳定性 | 岭（ridge 1e-8 几何增长直至对称正定可逆） | 通常不岭（G 直接可逆；A22⁻¹ 由稀疏分解） | 方法不同 | 低；数值影响小 |
| A22⁻¹ | 岭求逆 | 稀疏 Cholesky 求解（Henderson 快速法） | 结果等价，实现不同 | 低 |
| 对称化 | 强制 `0.5(H+H′)` + max-abs 检查 | 构造即对称 | ✓ | 低 |
| 基因型顺序 | 必须等于过滤后样本集（`ng==n_keep`，否则 throw） | 灵活性更高 | — | 用法限制（§14） |

---

## 4. Method 1 数学（现状）

模型（逐 SNP）：

```
y = Xb + x_j·β_j + u + e,   u ~ N(0, H_sub·σ²g),  e ~ N(0, I·σ²e)
```

1. **VCE**（仅基因型样本）：`H_sub = { (H⁻¹)^{-1} }[g,g]`，REML 估计 σ²g, σ²e（EigenVce，`use_covariate_projection`）。
2. **扫描**：σ²g/σ²e 固定 → 对每个 SNP j（个体自身的 genotype 列 x_j）做 score/Wald 检验 H₀: β_j = 0，输出 BETA/SE/P。
3. LOCO：染色体 c 用 `G_loco ∝ (S_global·G_global − S_chr·G_chr)` 重建并重做 H⁻¹→H_sub（代价 O(n³) 全 H 求逆 × 染色体数）。

性质：u 是随机核（kernel），**不做 EBV 点估计求解**；不产出 marker-effect 型输出。这正是与 Method 2 的本质差异（§6）。

---

## 5. Canonical ssGBLUP（Method 2 核心，已核实）

模型（Wang et al. 2012; Aguilar et al. 2010）：

```
y = Xb + Zu + e,   u ~ N(0, Hσ²u),  e ~ N(0, Iσ²e),  λ = σ²e/σ²u
```

MME：

```
[X′X      X′Z     ] [b]   [X′y]
[Z′X  Z′Z + H⁻¹·λ ] [û] = [Z′y]
```

- û 覆盖**全部动物**（基因型 + 未基因型），基因型动物子集记 â_g（GEBV）。
- H⁻¹ 如 §3.1；G 可被 α-mixed、tuneG 对齐、τ/ω 缩放。
- VCE 常用 AI-REML / EM-REML（blupf90 生态：airemlf90、`OPTION` 参数化）。
- 求解：PCG 或直接分解（CosmicSolver 用 SimplicialLDLT / MatrixFreePCG）。

---

## 6. SNP Back-solving 公式（已核实）

### 6.1 三种等价表达（都出自文献）

| 形式 | 出处 | 说明 |
|---|---|---|
| **a) û = D·M′·(MDM′)⁻¹·â_g** | Wang et al. 2012（wssGBLUP 一般形式；G* = MDM′） | 加权反解；D=I 退化为 b) |
| **b) â = Z′(ZZ′)⁻¹·â_g**（D=I） | VanRaden 2008；Wang et al. 2012（特例）；常写作 `Z′(ZZ′)⁻¹û_g` | 最小二乘反解；标准 ssGWAS 教科书式 |
| **c) â = Z′·G⁻¹·â_g** | Strandén & Garrick 2009；Leite et al. 2024 Eq.(1) | 当 G = ZDZ′/S（S=Σ2p(1−p) 或 M）时与 a) 代数等价；**这是代码实现最廉价的形式** |

其中 M = X−2p（等位基因计数中心化），Z 为 M 的缩放版（缩放方式须与 G 构造一致，见 §10）。

### 6.2 CosmicScanner 可用原语（libCosmicCore，已确认在子模块内）

- `PlinkReader::multiply_Zt_v(v, y, n_threads)` / `(v, y, keep_indices, n_threads)` — `libCosmicCore/src/plink_reader.cpp:355,411`（y_j = Σ (x_ij−mean_j)·inv_std_j·v_i，即**标准化的 Z′·v**）
- `snp_means` / `snp_inv_stds = 1/√var`（`plink_reader.cpp:190-234`；`getMeans/getInvStds` `plink_reader.h:107-108`）
- 同样实现于 `PgenReader`（`pgen_reader.cpp:1005`）与 `BgenReader`（`bgen_reader.cpp:544`）
- 因此 **e_effects = Z′·(G⁻¹·â_g) 直接用现有 `multiply_Zt_v` 实现**（注意 §10.3 的刻度换算：`multiply_Zt_v` 是标准化 Z，back-solve 后需还原到替换刻度）。

> **注意刻度**：`multiply_Zt_v` 的 Z = (x−mean)/std（标准化），而文献 a)/b) 的 M/Z = x−2p 或 (x−mean)/√(2p(1−p))。**输出 SNP 效应定义域（allowance scale）必须与文档声明一致**：COSMIC 路径产出"标准化剂量刻度"效应，需乘 1/√var_j？——不：若 G 用 2p(1−p) 缩放，则 b) 的 (ZZ′)⁻¹ 恰抵消 G 的 1/S 因子，效应定义域由"G 的缩放 + Z 的定义"共同决定。实现时以"与 BLUPF90 对齐的对照实验"定标（§10.3 有完整换算表）。

### 6.3 SNP 效应的 P 值（已核实，Aguilar et al. 2019）

postGSf90（`OPTION snp_p_value`）的公式（GSE 51:28; 又见 fgene.2020.00692 使用例）：

```
p_i = P_t( û_i / √(σ_i²/n) ,  n − 1 )
```

- û_i：第 i 个 SNP 的 backsolve 效应；σ_i²：该 SNP 遗传方差；n：有该 SNP 基因型的动物数；P_t：t 分布 CDF。
- σ_i² 与 var(û_i) 需要 GEBV 的预测误差（协）方差 C^{u₂u₂}（= MME 系数阵逆的基因型块，O(n_g³)）；Leite et al. 2024 给出 APY 近似（GSE 56:53）。
- **语义**：这是"估计效应 / 其估计标准差"的 t 型检验，**不是固定效应 Wald P**。因此用户 §16/Q16 的答案：经典 backsolve ssGWAS **可以**有 P（Aguilar 2019），但①是 t 型、②需 MME 系数阵求逆（或 APY 近似）、③不是固定效应假设检验。**Phase B 若只做效应 + 方差，可暂不输出 P**（marker/window variance 排序是更传统的做法）。

---

## 7. SNP Weights（D 矩阵，Method 3 用）

文献变体（Wang et al. 2012; Zhang et al. 2016 综述）：

| 权重 | 公式 | 用途 |
|---|---|---|
| 经典 wssGBLUP | w_j = |û_j|·√(2p_j(1−p_j)) | 缩放效应权重（Wang 常用） |
| 方差占比形式 | w_j = û_j²·2p_j(1−p_j) | 方差贡献型权重 |
| 效应平方 | w_j = û_j² | 简化变体 |
| 迭代归一 | w ← w / mean(w) | 保持 G* 的方差尺度（twssGBLUP 必做） |

实现注记：`libCosmicCore/src/grm_builder.cpp` **已支持 `snp_weight_file`**（对 Z 列乘 √w，`w = sqrt(snp_weights_block[j])`），正好实现 `G* = M·diag(w)·M′/S`。所以 Method 3 的 G 构造端 0 新增代码（仅接口组装）。

---

## 8. SNP / Window 遗传方差（Method 2 标准产出）

```
Var_SNP_j     = 2p_j(1−p_j) · û_j²
%Var_window   = Σ_{j∈win} 2p_j(1−p_j)·û_j² / Σ_all 2p_j(1−p_j)·û_j²
```

- 窗口：滑动窗口，常见 10 个连续 SNP（BLUPF90 教程默认窗口=10 SNP；也可物理窗口如 1 Mb）。
- 传统显著性经验法则：窗口解释遗传方差 >1% 即视为显著区域（Fragomeni et al. 2014; Wang et al. 2014）。
- CosmicScanner 现状：无任何窗口归并（§8.4-F 已判定 `window` 关键字仅存在于 GEMMA/RRM 语境）。

---

## 9. 迭代加权 ssGBLUP（twssGBLUP）——设计规范（本轮**不实现**）

算法（Wang et al. 2012; Zhang et al. 2016）：

```
r = 0:  G* ← G（无权重）
repeat:
  r += 1
  1. 求解 ssGBLUP（MME，H*⁻¹ = A⁻¹ + diag(τ·G*⁻¹ − ω·A22⁻¹)）
  2. â_g^{(r)} = GEBV（基因型子集）
  3. û^{(r)} = D·M′·G*⁻¹·â_g^{(r)}（D=I 时 û = Z′G*⁻¹â_g）
  4. w_j^{(r)} = |û_j|·√(2p_j(1−p_j))（或变体）
  5. w^{(r)} ← w^{(r)} / mean(w^{(r)})            （归一化，保持 G* 方差尺度）
  6. G* ← M·diag(w^{(r)})·M′ / S
  until r == R_max（典型 1–3）或 ||w^{(r)} − w^{(r−1)}|| < tol
```

工程注意：归一化（§7 + 步骤 5）是避免 G* 尺度漂移、进而污染 VCE 的**必要措施**；每轮需重做 H*⁻¹ → MME → backsolve（代价 ≈ Method 2 × R）。此设计仅存档，**禁止本轮实现**（用户 §24 约束）。

---

## 10. 标准化 / 缩放变体与换算（对齐 BLUPF90 的关键）

### 10.1 G 构造变体（CosmicScanner `--grm-alg 1..6`）

| # | 算法 | 缩放 | BLUPF90 对应 |
|---|---|---|---|
| 1 | VanRaden（默认） | 1/Σ2p(1−p) | preGSf90 默认 |
| 2 | Zeng | ΣM² 归一 | — |
| 3 | Yang | 1/√(2p(1−p)) 每列 | GCTA/GEMMA 默认 |
| 4 | Vitezica | 纯合矫正 | — |
| 5 | TrueVar | 经验方差 | — |
| 6 | EMmvGWAS | 兼容路径 | — |

**对比基准时：频率需来自同一基因型样本、同一 MAF 过滤、同一算法（默认都用 VanRaden 2p(1−p)）。**

### 10.2 尺度必要对齐项（跨软件对比 checklist）

1. 频率源（基因型样本 vs 全群加权）；
2. MAF/HWE 过滤 前后一致；
3. G 算法与缩放一致；
4. α 混合对象（I vs A22）与 α 值；
5. τ、ω；
6. `tuneG`（G=bG+a → A22）开/关——**BLUPF90 默认对 ssGBLUP 有影响力，对比前必须两端一致或明确其效应**；
7. backsolve 用的 Z 定义（x−2p / (x−mean)/√var / (x−mean)/std）——直接决定效应的数值刻度；
8. ridge 是否启用（CosmicScanner 用岭求逆，BLUPF90 一般不；对效应数值影响二阶，但应记录）。

### 10.3 Back-solve 刻度换算表（设计值）

| 输入 Z 定义 | G 缩放 | backsolve 形式 | SNP 效应定义域 |
|---|---|---|---|
| M = X − 2p | 1/Σ2p(1−p) | û = M′G⁻¹â_g | 等位替代刻度（天然） |
| Z_std = (X−mean)/√var | （同上，但与 Z_std 不匹配） | 避免混用；换算 û_allow = û_std·(√var/√... ) | 需显式换算 |
| `multiply_Zt_v` 的 Z = (X−mean)/std | 仅当 G 也构建于该 Z | 直接 | 标准化刻度 |

**结论：实现原则 = "back-solve 的 Z 定义必须与 G 的构造定义成对一致"**，并在输出 header 声明刻度。若要与 POSTGSF90 数值对齐，应采用 ×10.3 第 1 行（M=X−2p 或等价），必要时用 `snp_means` 把 `multiply_Zt_v` 的标准化输出换算回来（`snp_means`/`snp_inv_stds` 都已提供）。

---

## 11. 计划实现（Phase B 设计——已于 2026-08 按 B 完整实施，实施细节见 §16）

### 11.1 目标架构：`SingleStepEngine`（现状 = Method 1；扩展 = Method 1 + Method 2 并存）

```
SingleStepEngine
├── 1. H 构造（现状，已实现）
│     A⁻¹（Henderson） · A22⁻¹（Schur+岭） · G（GrmBuilder 1..6, --snp-weight-file）· H⁻¹（τ/ω/α)
├── 2. 遗传评价（Method 2 新增）
│     MME：X′X, X′Z, Z′Z+H⁻¹λ   求解：移植 CosmicSolver `gs_solve.cpp`
│     （SimplicialLDLT 直接 / MatrixFreePCG+VCE warm start；λ 来自 §2.1 VCE 或外部 σ²）
│     产出：EBV（全群）→ GEBV（基因型子集）→ *.ebv.txt（替代 writeEbv stub）
├── 3. 固定 SNP 关联（Method 1，现状）
│     EigenVce 固定 VC → LmmEngine 扫描 → *.assoc.txt（keep）
└── 4. Back-solving（Method 2 新增）
      4a. q = G⁻¹·â_g（稠密；n_g 大时 APY 渐进——设计备注）
      4b. effects = Z′q（multiply_Zt_v，现有原语）+ 刻度换算（§10.3）
      4c. Var_SNP_j = 2p_j(1−p_j)·û_j²；%Var_window（滑动窗口，默认 10 SNP）
      4d.（可选，Phase C）PEV → Aguilar 2019 P 值
      4e.（存档，不做）wssGBLUP 迭代
```

### 11.2 复用表（Component | 位置 | 现状 | 复用 | 改动）

| 组件 | 位置 | 现状 | 复用？ | 改动 |
|---|---|---|---|---|
| 系谱解析 + A⁻¹ | `pedigree_relationship.{h,cpp}` / `a_matrix.cpp:446` | ✓ 稳定 | ✓ | 无 |
| A22⁻¹ | `pedigree_relationship.cpp:219-236` + `single_step_model.cpp:196-212` | ✓ | ✓ | 无 |
| G（6 算法 + blend + **snp_weight_file**） | `libCosmicCore/src/grm_builder.cpp` | ✓ | ✓ | 无（D 加权零新增） |
| G⁻¹（岭） | `invertSpdWithRidge`（`single_step_model.cpp` / `pedigree_relationship.h`） | ✓ | ✓ | 无 |
| H⁻¹ | `single_step_model.cpp:214-238` | ✓ | ✓ | 无 |
| VCE（REML，基因型子块） | `runVCE` | ✓（但仅基因型） | ✓（作 MME λ 来源或对比） | 记录差异 |
| **MME 联合求解（含未基因型）** | **CosmicScanner：无** | — | **从 CosmicSolver 移植**：`src/pipeline/solve/gs_solve.cpp`（LDLT/PCG）、`src/pipeline/io/matrix_builder.cpp:500-633`（Hinv/tuneG/blend/ridge/τω 组装逻辑） | 新模块 `src/models/single_step_gs.{h,cpp}` |
| Z′·v / 刻度 | `libCosmicCore/src/plink_reader.cpp:355,411`（PGEN:1005, BGEN:544） | ✓ | ✓ | 刻度换算层（§10.3） |
| Back-solve + window 归并 | **无** | — | 新写（轻量，~200–400 行） | 新模块 |
| SNP P（Aguilar 2019） | **无** | — | 设计：C^{u₂u₂} 需要 LHS 求逆或 APY；Phase C | 新模块 |
| 外部基准 | **无（本地）** | CosmicSolver `tests/run_test_ssgblup_smoke.sh`（vs HIBLUP，EBV 级） | ✓ | 扩展为 SNP-effect 级对比（§13） |

### 11.3 CLI 设计（不改现状）

```
--single-step                   现状（Method 1，backward compatible）
--single-step-method lmm        （默认 == 现状）
--single-step-method backsolve  （Phase B experimental；需系谱+基因型+表型）
  --ssgwas-window 10            （SNP 窗口大小）
  --ssgwas-write-ebv            （*.ebv.txt）
  --ssgwas-p-values              （Phase C；默认关）
--snp-weight-file <wgt>         （GrmBuilder 已支持，Method 3 预留）
```

### 11.4 Phase B 推荐（风险排序回答）

- **B（推荐首选）**：MME 求解 + 无权重 backsolve（û = Z′G⁻¹â_g）+ SNP/窗口方差。全部公式已核实；复用率最高（MME 有 CosmicSolver 现成实现、Z′v 在库）；基准用 HIBLUP/（未来）BLUPF90。
- **C（次选，+P）**：在 B 基础上做 Aguilar 2019 P 值——需要 MME 系数阵求逆（小数据可行）；大 n_g 须 APY 近似（Leite 2024），工程量大。
- **D（weighted/iterative，不建议本轮）**：GrmBuilder 的 `snp_weight_file` 使 G* 构造廉价，但迭代 + 归一化 + 数值验证成本高、且无本地外部基准需求。**明确不默认 D**（与用户 §23 一致）。

---

## 12. 验证计划（Phase B 落地的验收）

| 层级 | 验证 | 判据 |
|---|---|---|
| 内部数学 | A=I 退化（全部基因型 且 无亲缘 → H⁻¹ = I + G⁻¹ 语义自洽） | 已有 `tests/run_ssgwas_validation.py` 精神扩展 |
| 内部数学 | n_g=1 / 无系谱等退化路径 | 不崩溃、数值有限 |
| backsolve 一致性 | û = Z′G⁻¹â_g vs û = Z′(ZZ′)⁻¹â_g（同数据） | 至少 1e-8 级一致（G 用同一 Z/S） |
| backsolve 刻度 | 模拟 QTL：一个 SNP 效应 δ；GEBV 反解恢复 δ（rank 满时） | 单 SNP 模拟恢复 ±1% |
| 窗口方差 | 已知 QTL 窗口 | top 窗口包含真值 |
| 与 CosmicSolver | 同数据 ssGBLUP → GEBV 全等 | EBV corr ≥ 0.9999（在两端同 H 参数下） |
| 与 HIBLUP | 现有 `run_test_ssgblup_smoke.sh` 扩展 SNP-effect 对照 | effect corr、top-SNP overlap |
| 与 BLUPF90/postGSf90 | 需软件；先做 §10.2 对齐 checklist | 对齐后 effect corr / window 排序可比 |

---

## 13. 外部对比计划（BLUPF90 / HIBLUP）

1. **HIBLUP（本机可用）**：`/disk228/huangjie/software/miniconda3/envs/benchmark/bin/hiblup`。已用于 CosmicSolver ssGBLUP VCE/EBV 基准（`CosmicSolver/tests/run_test_ssgblup_smoke.sh`，N=2000/1000 基因型/2000 SNP）。扩展：HIBLUP 的 GBLUP/ssGBLUP → SNP effect 与 CosmicScanner backsolve 对照。
2. **BLUPF90 家族（本机无）**：blupf90 + postGSf90（`nce.ads.uga.edu/software`）。若后续可用：按 §10.2 checklist 对齐（频率源、算法、α 混合对象、τ/ω、tuneG、Z 定义）。**未对齐前不得声称"不一致"**（用户 §21）。
3. 数值协议：single-thread（`OMP_NUM_THREADS=1`）→ 确定性对比；报告参数记录表。

---

## 14. 当前限制（Limitations）

1. **VCE 仅基因型样本**（`single_step_model.cpp:270-272` 强制 `ng==keep.n_keep`）；未基因型动物不参与 σ² 估计（Method 1 语义，非 bug）。Method 2 的 MME 用同一 σ²（或 `--h2-tau/--h2-omega` 外部比率），全模型 VCE 未实现。
2. ~~`writeEbv()` stub~~ → 已由 backsolve 的 `*.ss_backsolve.ebv.txt` 替代（全群 EBV + 基因型 GEBV，§16）。
3. ~~无 MME ⇒ 无法 backsolve~~ → MME（LDLT/PCG）已实现；无 P 值（Phase C 预留）、无 SNP weights（Phase D 预留）。
4. H 子块重建是全 H 稠密求逆（O(n³)）；LOCO 每染色体一次 → 大 n 下昂贵（设计备注：可改为 H 子块的直接求逆技巧）。
5. `--genotyped-list` 与过滤后样本集的严格顺序绑定（`ng==n_keep`），灵活性受限。
6. 岭求逆（ridge 1e-8→几何增长）与 BLUPF90 的无岭实现有二阶数值差。
7. 无 G−A22 对齐（tuneG）；跨软件数值对比前必须处理。
8. `multiply_Zt_v` 输出标准化刻度效果——输出前必须换算并声明刻度（见 §16 的 EFFECT_ALLELE/EFFECT_STD 定义）。
9. Aguilar 2019 P 值需要 LHS 逆（O(n³)）或 APY 近似（Leite 2024）——Phase C。
10. 本机无 BLUPF90/postGSf90；外部数值基准以 HIBLUP 为当前现实下限。

---

## 15. 参考文献（联网核实，2026-08）

1. **Aguilar I, Misztal I, Johnson DL, et al. (2010)**. Hot topic: A unified approach to utilize phenotypic, full pedigree, and genomic information for genetic evaluation of Holstein final score. *J Dairy Sci* 93:743–752. https://doi.org/10.3168/jds.2009-2730 — H⁻¹ 的 τG⁻¹−ωA22⁻¹ 形式。
2. **Legarra A, Aguilar I, Misztal I (2009)**. A relationship matrix including full pedigree and genomic information. *J Dairy Sci* 92:4656–4663. https://doi.org/10.3168/jds.2009-2412 — H 的原始定义。
3. **Christensen OF, Lund MS (2010)**. Genomic prediction when some animals are not genotyped. *Genet Sel Evol* 42:2. https://doi.org/10.1186/1297-9686-42-2 — H 等价性；τ>1 讨论。
4. **Wang H, Misztal I, Aguilar I, Legarra A, Muir WM (2012)**. Genome-wide association mapping including phenotypes from relatives without genotypes. *Genet Res* 94:73–83. https://doi.org/10.1017/S0016672312000081 — **wssGBLUP 与 û = DM′G*⁻¹â**。
5. **VanRaden PM (2008)**. Efficient methods to compute genomic predictions. *J Dairy Sci* 91:4414–4423. https://doi.org/10.3168/jds.2007-0980 — backsolve 线性变换（Z′(ZZ′)⁻¹â_g 系）。
6. **Strandén I, Garrick DJ (2009)**. Technical note: Derivation of equivalent computing algorithms for genomic prediction and association of traits. *J Dairy Sci* 92:2971–2975. https://doi.org/10.3168/jds.2008-1294 — â = Z′G⁻¹â_g 等价式。
7. **Aguilar I, Legarra A, Cardoso FF, et al. (2019)**. Frequentist p-values for large-scale-single step genome-wide association, with an application to birth weight in American Angus cattle. *Genet Sel Evol* 51:28. https://link.springer.com/article/10.1186/s12711-019-0469-3 — **SNP P 值公式 p_i = P_t(û_i/√(σ_i²/n), n−1)**；`OPTION snp_p_value`。
8. **Leite NG, Bermann M, Tsuruta S, et al. (2024)**. Marker effect p-values for single-step GWAS with the algorithm for proven and young in large genotyped populations. *Genet Sel Evol* 56:53. https://link.springer.com/article/10.1186/s12711-024-00925-3 — backsolve（VanRaden/Strandén-Garrick/Wang 变换）+ PEV/APY 近似 P 值。
9. **Fragomeni BO, Misztal I, Lourenco DAL, et al. (2014)**; **Wang H, et al. (2014)** — SNP/窗口遗传方差与 >1% 经验显著性（见 [8] 引用综述）。
10. **Masuda Y.** BLUPF90 教程：ssGWAS 标准工作流（blupf90 → postGSf90 → SNP effects → 窗口）。 https://github.com/masuday/Blupf90TutorialStandard/blob/master/genomic_gwas.md
11. **Front. Genet. 11:692 (2020)**（Landrace 猪饲料效率 ssGWAS 应用例，采用 Wang 2012 + Aguilar 2019 公式）。 https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00692/full
12. **Misztal I, Tsuruta S, Lourenco D, et al.** BLUPF90 家族手册（MME、H⁻¹、tuneG/alpha_blend/tau/omega OPTION）。 http://nce.ads.uga.edu/wiki/doku.php?id=documentation

---

*本文档为审计/实施文档。Phase B 实施与验证见 §16；配套 `docs/GWAS_MODEL_DEFINITIONS.md` §8/§9/§10/§12 同步更新。*

---

## 16. Phase B 实施状态（2026-08）

### 16.1 已交付（`--single-step-method backsolve`）

| 组件 | 实现位置 | 状态 |
|---|---|---|
| CLI（`--single-step-method`, `--h2-tau`, `--h2-omega`, `--blend-g`, `--ssgwas-window-snps`, 求解器选择） | `src/config.h`, `src/main.cpp` | ✓ |
| 系谱 + H⁻¹ = A⁻¹ + diag(τG⁻¹ − ωA22⁻¹)，G = MM′/S（S = Σ2p(1−p)）。**τ、ω 是用户设置的调谐权重（默认 1.0），不是 h² 的派生量**（`--h2-tau/--h2-omega` 命名有历史误导性；Phase 6.5 审计 GATE2 已独立证实） | `src/models/single_step_model.cpp`（`runBacksolve`） | ✓ |
| MME `[X′X X′Z; Z′X Z′Z + λH⁻¹]`，**λ = σe²/σa²（Henderson 方向，代码 `single_step_model.cpp:861`）**；实现用 LDLT（拦截列最左消除 = mu-absorbed 精确等价，Phase 6.5 GATE3）或 PCG | 同上（SimplicialLDLT 直接 / MatrixFreePCG，`--ssgwas-solver pcg`） | ✓ |
| G⁻¹ 岭逆 ＋ backsolve q = G⁻¹·â_g；effects = Z′q（`multiply_Zt_v`） | `invertSpdWithRidge` + `libCosmicCore/src/plink_reader.cpp` | ✓ |
| 输出：`*.ss_backsolve.ebv.txt`（全部动物 EBV/GEBV + 统计）、`*.snp.txt`、`*.window.txt`（滑动窗口，不跨染色体）、`*.meta.txt` | `runBacksolve` | ✓ |
| 刻度：`EFFECT_STD` = `EFFECT_ALLELE`×σ_z（重建一致刻度，Z·EFFECT_STD = M·EFFECT_ALLELE = GEBV）；`EFFECT_ALLELE = (z′q/σ_z)/S`（等位替换效应）；SNP 方差 = 2p(1−p)·EFFECT_ALLELE²；窗口 %var（定义 A，见 §17） | §16.3 | ✓ |

### 16.2 验证结果（全部 PASS，2026-08-16）

1. **解析门**（`tests/run_ssgwas_backsolve_analytic.py`）：106 只两代系谱、3000 SNP 的封闭公式对拍（G2/O3 骨架，P1/P2 路径 + pcg 变体，τ=ω=0 的 pedigree-BLUP 退化）。max|diff|：EBV/GEBV 5e-7、EFFECT_ALLELE 4.9e-7、EFFECT_STD 1.7e-7、SNP_VAR 4.7e-8、PCT 4.7e-8、WIN/WIN_PCT 3.8e-7（阈值 REL=1e-4/ABS=1e-8）。Phase 6.5 审计更正：**拦截列最左消除与 mu-absorbed 型在 1e-15 级精确一致，系数阵非奇异（min_sv≈1，cond≈32，Hinv·1≠0），不存在 (μ, a→a−c) 自由规范自由度**（§17 GATE3）。
2. **模拟门**（`tests/run_ssgwas_backsolve_simulation.py`）：900 动物三代系谱、1000 SNP、5 因果位点、300 基因型 + 650 表型（含未基因型个体）→ Pearson(top200)=0.81、Top10=5/5、Top20=5/5；**系谱行序洗牌后 EBV/效应逐位一致（max diff = 0）**；LDLT vs PCG 一致（1e-8 级）。
3. **回归网**：Method-1 单步回归（`tests/run_test_singlestep.sh`）、GEMMA 精度（beta/SE/log10p corr ≥ 0.99）、共享 RRM 全套（`tests/run_test_rrm_regression_suite_ci.sh --build`）全部 PASS。
4. **GRM 确定性修复**：`libCosmicCore/src/grm_builder.cpp` 的 `std::vector<bool> snp_valid`（compute/computeSplit 两处）在 OpenMP 下是位压缩数组，不同元素的并发写共享字节 → G/Σ2p(1−p) 跨线程数漂移（实测 35.9–37.8 抖动）。改为 `std::vector<char>` 后完全确定（37.7769 = Python Σ2p(1−p) 参考 37.776857 一致，max freq diff 4.5e-7）。三份 libCosmicCore（scanner/solver/root 子模块）已同步修改。

### 16.3 Method 2 数学约定（实现用，与 §10.3 一致）

- 基因型 x ∈ {0,1,2}，M = x − mean，S = Σ_j 2p_j(1−p_j)，G = MM′/S；G⁻¹ 为岭逆（ridge 1e-8→几何增长，同 `invertSpdWithRidge`）。
- backsolve：q = G⁻¹·â_g；`EFFECT_ALLELE_j = (z′q/σ_z)_j / S`（allele-substitution 刻度，σ_z = √(2p(1−p))）；`EFFECT_STD_j = EFFECT_ALLELE_j × σ_z`（重建一致刻度：Z·EFFECT_STD = M·EFFECT_ALLELE = GEBV(基因型块)）。
- SNP 遗传方差 = 2p_j(1−p_j)·EFFECT_ALLELE_j²；%Var_window = Σ_window Var_SNP / Σ_all Var_SNP（滑动窗口默认 10 SNP，不跨染色体）。
- P 值（Aguilar 2019）**未实现**（Phase C 预留）；`BACKSOLVE` 输出表因此无 p 列，靠效应/方差排序。
---

## 17. Phase 6.5 独立数学红队审计（2026-08-16）

> 审计方式：**不信任** Phase 6 的全部实现/测试/文档结论，一切从第一性原理重新推导，并与外部权威文献逐条核对；实现代码在审计结论报告前未被改动。审计脚本 `tests/run_ssgwas_math_audit.py`（全部 7 门 PASS）。

### 17.1 审计结论汇总（Component | Phase 6 表述 | 独立理论/文献 | 判定 | 变更 | 验证）

| 组件 | Phase 6 表述 | 独立推导/文献 | 判定 | 变更 | 验证 |
|---|---|---|---|---|---|
| λ（MME 方差比） | λ = h²/(1−h²) | Henderson MME：R=σe²I, CoVar(a)=Hσa² ⇒ **λ = σe²/σa²** | **表述错误，代码正确**（`single_step_model.cpp:861`） | 仅文档 | GATE1：显式方差分量法（不用 λ/h² 捷径）下 0.4/0.6→λ=1.5、0.8/0.2→λ=0.25 两方向区分案例 max\|diff\|≈4.8e-6 |
| τ/ω 语义 | τ=h²，ω=(1−h²)τ | H⁻¹ = A⁻¹+τG⁻¹−ωA22⁻¹ 中 τ、ω 仅 H 调谐权重 | **表述错误，代码正确**（`main.cpp:536-537` 仅校验非负有限） | 仅文档 | GATE2：τ=1.5 / ω=1.2 各与对应显式参考一致且与原解显著不同（敏感度 0.26/10.2） |
| MME 结构 | 截距+动物系统存在 (μ, a→a−c) 规范自由度 | 拦截列最左消除 = mu-absorbed 型**恒等式**；Hinv·1≠0（实测 max=20）、min_sv≈1.0、cond≈32 ⇒ 非奇异 | **"规格自由度"表述不成立**；实现无碍 | 仅文档 | GATE3：5 数据集 canonical ≡ mu-absorbed（≤2.3e-15） |
| backsolve 重建 | a_g ≈ M·α | α=(M′G⁻¹â_g)/S ⇒ M·α = â_g（精确）；Z·γ=â_g ⇒ γ=α·σ | **正确** | 无 | GATE4：Pearson 0.9995、slope 0.991、max_abs 0.045（岭/混合同源噪声量级） |
| SNP_EFFECT_STD 刻度 | = z′q（standardized） | 重建一致刻度应为 α_j·σ_z_j（= z′q·σ²/S） | **代码刻度缺陷（F1）**：原列 ≠ 任何可解释"效应"，约为 γ 的 S/σ²≈90–200 倍 | **列重定义**：EFFECT_STD = EFFECT_ALLELE×σ_z | Gate 4 复跑：Z·EFFECT_STD ≡ M·EFFECT_ALLELE（RMSE 0.0187 逐位一致） |
| SNP 方差 | 2p(1−p)·EA² | 每位点边际加性方差的样本解释 Var_s(M_j·EA_j) 与 2pq·EA² 一致（加权均值比 0.965，受噪声 EA_j 影响） | **正确** | 无 | GATE5 |
| 窗口方差 | Σ_window 2p(1−p)û² | 定义 A（边际方差求和） vs 定义 B（Var(M_win·α_win)）**在 LD 下不可互换**（示例 B/A=7.97） | 实现=定义 A，正确但必须声明 | **元数据新增声明字段**（SNP_VARIABLE_DEFINITION 等） | GATE6：逐窗口 max\|diff\| 1.4e-7 |
| PCT 分母 | PCT=SNP_var/TOTAL | 分母=Σ_j 2p_j(1−p_j)EA_j²（**非 σa²**）；PCT 是"边际方差占比份额" | 正确，需声明 | 元数据新增 PCT_VARIANCE_DENOMINATOR | GATE7 |

### 17.2 GATE 细节（seed 11 分析数据集；`tests/results/ssgwas_math_audit/`）

- **GATE1**：显式方差分量参考（R=Ve·I、CoVar(a)=H·Va 不变式 + 拦截保留的完整规范解，全程无 λ/h² 捷径）。两个方向区分案例（σa²=0.4/σe²=0.6→λ=1.5 与 0.8/0.2→λ=0.25；**刻意不用 h²=0.5 的歧义案例**）max\|diff\| 4.6e-6/4.8e-6（浮点路径级）。
- **GATE2**：τ=1.5、ω=1.2 独立运行 vs 对应显式参考（τ 进 Hinv 的位置）max\|diff\| 4.7e-6/3.8e-5；相对 τ=ω=1 基线的敏感度 0.258/10.17 ⇒ τ/ω 是独立调谐权重。
- **GATE3**：5 个随机数据集的完整规范解 vs mu-absorbed：max\|diff\|=2.3e-15；min_sv(C_mu)∈[0.99,1.05]、cond∈[31,33]；max\|Hinv·1\|=20（=τ·Ginv·1，因 G·1=0.05·1 混同岭 ⇒ 非零）。**不存在真零方向**。
- **GATE4**：重建统计（dosage 刻度 M·EA vs GEBV：Pearson 0.9995/Spearman 0.9963/RMSE 0.0187/max_abs 0.045/slope 0.9908/intercept ≈0；OLS 一致 1.000000；EA vs OLS β Pearson 0.9985）。F1 修复后标准化刻度 Z·STD 与 dosage 刻度**逐位同值**（marker-scale 不变性达成）。
- **GATE5**：SNP_VARIANCE 列 vs 2p(1−p)EA² max\|diff\|≤1.4e-7；样本方差解释加权均值比 0.965/0.969（期望 1.0，噪声来自 EA_j 估计）；不涉及自由度修正误差。
- **GATE6**：独立 LD 数据集（n=200，m=10，隐性相关 0.95，MAF 0.3，等效应 α=1）：定义 A=3.93，定义 B=31.28，**B/A=7.97** ⇒ LD 下两定义不可互换；C++ 逐窗口 = 定义 A（82/82 窗口 max\|diff\|=1.4e-7）。
- **GATE7**：PCT=SNP_var/TOTAL=Σ2p(1−p)EA² 份额；与 σa² 无关，已在元数据声明。

### 17.3 外部文献核对（Reference | Equation | CosmicScanner | 判定）

| 文献 | 方程 | CosmicScanner 对应 | 判定 |
|---|---|---|---|
| Henderson (1984) MME | [X′R⁻¹X X′R⁻¹Z; Z′R⁻¹X Z′R⁻¹Z+G⁻¹] | ×σe² 后 λ=σe²/σa² 施加于动物块 | 一致 |
| VanRaden (2008) | G=MM′/Σ2p(1−p) | grm_builder VanRaden 算法（mean=2p、权重 2p(1−p)） | 一致（S=36.7035 对拍） |
| Strandén–Garrick (2009) | û = Z′G⁻¹â（G 由对应 Z 定义） | EFFECT_ALLELE=(Z′q)/S（M 刻度）、EFFECT_STD（Z 刻度）均重建 â_g | 一致 |
| Wang et al. (2012) | û = D M′G*⁻¹â（D=I 时 M′G⁻¹â） | q=G⁻¹â_g ⇒ EA=(M′q)/S（S 归一化） | 一致 |
| Legarra et al. (2009)/Aguilar et al. (2010)/Christensen–Lund (2010) | H⁻¹=A⁻¹+diag(G⁻¹−A22⁻¹) | SplitMatrixAdapter：A⁻¹+τG⁻¹−ωA22⁻¹（τ=ω=1 时为文献形） | 一致 |
| postGSf90/BLUPF90 文档 | backsolve 效应 + SNP/窗口方差报告 | 窗口=定义 A（边际方差求和）；无 P 值（Phase C 预留） | 一致（P 值缺失按设计） |

### 17.4 元数据新增字段（`.ss_backsolve.meta.txt`）

`SNP_VARIANCE_DEFINITION`（2p(1−p)EA²，边际加性方差）、`WINDOW_VARIANCE_DEFINITION`（定义 A：窗口内边际方差求和，非 Var(M_win·α_win)，LD 下不可互换）、`PCT_VARIANCE_DENOMINATOR`（Σ_j 2p_j(1−p_j)EA_j²，非 σa²）。

### 17.5 结论

完成 7/7 门通过。**无阻断性数学缺陷**；实施代码（MME 方向、H⁻¹ 结构、backsolve 重建、SNP/窗口方差、PCT 分母）经验证正确。修复两类非阻断问题：F1（`SNP_EFFECT_STD` 列重定义为重建一致的 EA×σ_z，原列 ~S/σ²≈90–200× 超出可解释刻度，无下游消费方）与 F2（"λ=h²/(1−h²)"、"(μ,a→a−c) 规范自由度"、τ/ω=h² 派生 三处文档表述更正 + 元数据定义字段）。Phase 6 的 EBV/GEBV/EFFECT_ALLELE/SNP_VARIANCE/WINDOW 全部数值结论在 F1 后不变（仅 SNP_EFFECT_STD 列值按新刻度重算）。

## 18. Phase 6.6 外部基准（HIBLUP 参考实现，2026-08-17）

> 完整报告与复现脚本见 `docs/SSGBLUP_EXTERNAL_VALIDATION.md` 与
> `tests/run_ssgwas_hiblup_benchmark.sh`（结果在
> `tests/results/ssgwas_hiblup_benchmark/`）。

### 18.1 已验证的外部语义（数值级，非日志文字）

| 项目 | HIBLUP | CosmicScanner | 判定 |
|---|---|---|---|
| A（PA） | Meuwissen-Luo | `buildA22`（同公式） | max diff = 0.0；A11 均值 1.000125/0.003574 ≡ HIBLUP 日志 1.00013/0.00357357 |
| 原始 G | M M′/Σ2pq + 对角线均值归一（Su GS/plink 诊断已证） | 公式同（不做归一；被 tuneG 线性项吸收） | 公式一致 |
| tuneG | b=(dA−oA)/(dG−oG), a=dA−b·dG（打印 0.996/0.005，确值 0.9956/0.0046） | `--ssgwas-tune-g` 同公式 | 一致 |
| 混合 | Gw = 0.95·G* + 0.05·A11 | `--ssgwas-blend-a22`（α=0.05） | 一致 |
| H⁻¹ | A⁻¹ + Gw⁻¹ − A11⁻¹（τ=ω=1） | 同 | 一致 |
| backsolve | SNP 效应用调优未混合的 G*⁻¹ | 同（`EFFECT_GINV=G_bs_tuned`） | 一致 |

### 18.2 关键发现：HIBLUP `--vc-prior` 内部缩放

HIBLUP `--mme --vc-prior Va,Ve` 会对其提供 VC 施加内部缩放（N=2000 下有效
λ = 0.98×提供值），表现为"假定完全对齐的实现只看似 slope≈0.988"。
通过 λ 扫描定位（精确匹配点 sf=0.98），对齐协议为：给 HIBLUP `Ve/0.98`。
校准后 CosmicScanner EBV 与 HIBLUP EBV **逐类完全一致**：
Pearson 1.00000000 / Spearman 1.00000000 / slope 1.000028 / intercept 0 /
RMSE 1.17e-5 / max|Δ| 5e-5（恰为 HIBLUP 输出 6 位小数精度；4 类动物
均同）。未校准（naive priors）时 slope 0.9879–0.9887，属参考侧内部行为，
非算法差异。

### 18.3 SNP 效应侧面

HIBLUP `--snp-effect` 输出 `add_a1`（原始刻度）vs `EFFECT_ALLELE`：
corr 0.99943、OLS slope 0.9447（数据相关内部刻度因子，与频率无关）、
intercept 5.3e-6。按规范 §11/12 报告为准定义对齐（corr 层面），
绝对刻度不做跨工具对等声明。GEBV 重建恒等式在本数据集复核：
M·EA vs 基因型 EBV corr 1.00000000、RMSE 2.0e-2。谱系行序洗牌：
结果逐位一致（max|Δ|=0）。

### 18.4 回归状态

Phase 6.5 全部门保持绿色：backsolve 全流程（analytic+simulation）、
数学红队 gates 1–7、SingleStep smoke、CTest 12/12、GEMMA accuracy、
RRM CI。运行时（N=2000/M=2000/8 线程）：HIBLUP mme 0.22s vs
CosmicScanner 2.66s（TIME_MME 1.52s / TIME_G 0.54s / TIME_GINV 0.18s /
TIME_HINV 0.22s / 其余 <0.01s）；仅记录，不做小 N 速度声明。
BLUPF90/postGSf90 就绪包在 `tests/ssgwas_external/blupf90/`（SKIP-safe）。
