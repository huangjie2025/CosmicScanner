# 变更记录

## 版本历史说明

- `v1.1.0` 是首个计划对外发布的 CosmicGWAS focused 版本，只正式交付单性状、多性状和重复力 GWAS。
- `1.0.0` 和 `1.0.5` 是根据 Git 历史回溯整理的内部开发里程碑，不是 Git 标签，也不代表曾经公开发布。
- 仓库开发分支还包含其他实验模块；这些模块不属于 CosmicGWAS v1.1.0 的支持范围。

## 开发分支（未发布实验功能）

- 新增 **ssGBLUP 育种值 back-solving ssGWAS（Method 2）**：`--single-step-method backsolve`（默认 `lmm` = 既有 H 核 LMM 关联 Method 1，backward compatible）。全系谱（含未基因型动物）求解单步 MME `[X'X X'Z; Z'X Z'Z + λH⁻¹]`（**λ = σe²/σa²（Henderson 方向）**；`--h2-tau/--h2-omega` 为 H⁻¹= A⁻¹+τG⁻¹−ωA22⁻¹ 的**调谐权重（默认 1.0，非 h² 派生）**；名称具有历史误导性，已在 Phase 6.5 审计 §17 澄清），SimplicialLDLT 直接求解（拦截列最左消除 = mu-absorbed 精确等价）或 `--ssgwas-solver pcg`（MatrixFreePCG + 初次 LDLT 预热），反解 `ĝ = Z'G⁻¹â_g`（VanRaden 2008 / Strandén-Garrick 2009 / Wang 2012；G 为 VanRaden 2p(1−p) 刻度 + `--blend-g` 默认 0.05）。
- 新增 ssGWAS SNP P-value APY 路径：`apy-exact` 标记为 validated scalable uncertainty path；`apy-approx` 标记为 `EXPERIMENTAL_APPROXIMATE`，`CALIBRATION_STATUS=LIMITED / REGIME_DEPENDENT`，`STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`，raw P-values only（`GENOMIC_CONTROL_APPLIED=false`），不得称为 calibrated P-value。
- Method 2 输出：`*.ss_backsolve.ebv.txt`（全部动物 EBV 与基因型子集 GEBV + Σκ/Rel/秩统计）、`*.ss_backsolve.snp.txt`（**`EFFECT_STD = EFFECT_ALLELE×σ_z`（重建一致刻度，Z·STD = M·ALLELE = GEBV）**、`EFFECT_ALLELE = (z'q/σ_z)/S`、`SNP_VAR = 2p(1−p)·EFFECT_ALLELE²`）、`*.ss_backsolve.window.txt`（滑动窗口 %Var_window，默认 10 SNP、不跨染色体；定义 A = 窗口内 2p(1−p)EA² 求和，非 Var(M_win·α_win)，LD 下不可互换）、`*.ss_backsolve.meta.txt`（含 `SNP_VARIANCE_DEFINITION`/`WINDOW_VARIANCE_DEFINITION`/`PCT_VARIANCE_DENOMINATOR`）。无 P 值（Aguilar 2019 需 LHS 逆，Phase C 预留）、无迭代加权（Phase D 预留）。
- Method 2 验证全 PASS：解析门 `tests/run_ssgwas_backsolve_analytic.py`（106 动物封闭公式对拍：EBV/GEBV max diff 5e-7、EFFECT_ALLELE 4.9e-7、SNP_VAR 4.7e-8、PCT/WIN 4.7e-8–3.8e-7；ldlt 与 pcg 一致；τ=ω=0 退化为 pedigree BLUP）；模拟门 `tests/run_ssgwas_backsolve_simulation.py`（900 动物三代系谱 1000 SNP：Pearson(top200)=0.81、Top10=5/5、Top20=5/5 因果恢复、系谱行序洗牌 EBV/效应逐位一致 = 0、ldlt vs pcg 1e-8）；回归网（`run_test_singlestep.sh`、GEMMA 精度、共享 RRM 套件 `run_test_rrm_regression_suite_ci.sh --build`）全 PASS。CTest 注册 `ssgwas_backsolve_validation`。
- **修复 libCosmicCore GRM 多线程不确定性问题**：`grm_builder.cpp` `compute()/computeSplit()` 使用 `std::vector<bool>`（位压缩）标记有效 SNP，OpenMP 并发写不同元素共享字节 → G 与 Σ2p(1−p) 随线程数漂移（实测 S 在 35.9–37.8 抖动）。改为 `std::vector<char>` 后完全确定（S=37.7769 = Python Σ2p(1−p) 参考 37.776857 一致）。三份 libCosmicCore（scanner/solver/root 子模块）同步修改。
- 文档同步：`SSGWAS_METHOD_AUDIT.md` 新增 §16 实施状态；`GWAS_MODEL_DEFINITIONS.md` §8/§9/§10/§12 更新；`MANUSCRIPT_DRAFT.md / _CN.md`、`manuscript.tex` + `references.bib`（新增 aguilar2010/legarra2009/christensen2010/stranden2009/wang2012/aguilar2019）写入 Method 2 公式与验证；`Tutorial_EN.md / Tutorial_CN.md` 新增 backsolve 用法与参数表。

- 新增 RRM 非平衡 record-level 时间变异 SNP 检验 `--rrm-gwas-method record`（Method A）与既有系数空间路径（`--rrm-gwas-method coefficient`，Method B，默认，保持 `--rrm-unbalance` 旧行为完全不变）。Method A 直接在记录上构建固定设计 `D(r,k)=x_subj(r),j·φ_k(t_r)`（精确 `x⊗Φ` Type B 模型），以 score/Wald 联合检验 `T=s'M⁺s`、`df=rank≥1`（1e-10·e_max 特征值截断伪逆），并输出显式诊断状态 `OK/RANK_DEFICIENT/ILL_CONDITIONED/NUMERICAL_FAILURE`，不再静默丢弃 SNP。null VCE 与 BLUP 阶段完全一致，扫描复用缓存分解。
- Method A 验证：`tests/run_rrm_balance_record_equivalence.sh` 平衡等价（N=200/M=300/T=8，BETA/SE/STAT/log10P Pearson ≥ 0.99999，max|ΔBETA|=2.5e-4，max|ΔSTAT|=0.084）；`tests/run_test_rrm_method_a_simulation.sh` 七情景仿真套件（常量/线性/二次/符号翻转轨迹因果恢复 + type-I error 校准）；`tests/run_rrm_method_ab_comparison.sh` 双路径对比（Method A 因果秩 1、Method B 秩 14，A 对时变效应功效更高）。设计矩阵/基函数单元测试 `tests/test_rrm_basis_identity.cpp`（Legendre 解析对照 max err 8.88e-16、B-spline 划分单元 3.33e-16）。
- RRM 非平衡外部对齐回归更新：`tests/run_gmat_unbalanced_alignment.sh` 改为在 autosome 上评估 top-SNP 一致性并同时记录全序/autosome 指标（GMAT 与 CS 的 GRM X 染色体覆盖范围差异已文档化为 docs/RRM_GWAS_METHOD_AUDIT.md §4.3；autosome Pearson r=0.9521 / Spearman 0.9401 / 全序 0.9470 / 0.9393 / autosome top UNC18848064 一致）；`tests/run_gma_compatible_fixreg_alignment.sh` 在 py3 重建的 GMA 上固定回归对齐 Pearson/Spearman=1.0、max|ΔP|=2.7e-9。

- 新增母体模型专属 SNP 检验 `--maternal-snp-test direct|maternal|joint`（Experimental）：`direct`（默认，与旧 `--maternal-model` 输出 byte-identical）保持子代自身基因型 single 固定效应；`maternal` 将**母亲基因型**（dam ALT dosage）作为单列固定效应检验 `H0: β_m=0`；`joint` 做 direct+maternal 2-df Wald `H0: β_d=β_m=0`（含 `COND/RANK/STATUS` 秩与条件数诊断，信息矩阵秩亏标记 `COLLINEAR`）。null 协方差与该模型 VCE 全分量一致：`σ²ad·K(子代)+σ²am·K(母亲)+σ²e·I`（subject-indexed，每染色体 LDLT 一次复用），LOCO 策略同 direct；无基因型 dam 显式过滤并在 assoc/`*.maternal.snp.meta.txt` 报告 `N_TOTAL/N_WITH_DAM_GENOTYPE/N_DROPPED_MISSING_DAM`。输出：`*.maternal.assoc.txt`、`*.joint.assoc.txt`、`*.maternal.snp.meta.txt`。PLINK-only（BGEN/PGEN 未接线）。模拟验证 10 场景全 PASS（A null λ_GC≈1 / B direct-only / C maternal-only / D 同向 / E 反向 joint 增强 / F dm-协方差无虚假信号 / G 缺失 dam 10/30/50% / H self-dam COLLINEAR）；`tests/run_test_maternal_snp_effects.sh` 为回归入口。simulation-validated only，不声称外部参考验证。
- 新增系谱动物模型 GWAS（`--animal-model`）：`y = Xb + x_j*beta_j + Z*a + e`，`a ~ N(0, A*sigma2_a)`，A 为系谱亲缘矩阵。共享 ssGWAS 的 `PedigreeRelationship` 系谱解析与 `A⁻¹`/`A22` 构建（`buildA22` 经 A⁻¹ 的 Schur complement 求解，不物化稠密 A）；`--write-ped-kernel` 导出分析子集精确 A22 供 GEMMA `-k` 使用；ssGWAS `--h2-tau 0 --h2-omega 0 --no-loco` 时两者逐 SNP 完全一致。
- `tests/run_test_animal_model.sh`：7 段模拟验证（causal 检出 + LOCO 提示、ssGWAS tau=0 等价、系谱乱序不变性、σ²a=0、高遗传力、不完整系谱、无祖先系谱）全部通过。
- `tests/run_animal_model_gemma_validation.sh`：外部 GEMMA 0.98.5 `-k` A22 核基准，β/SE/-log10P Pearson 通过阈值（beta≈0.99996、se≈0.9992、log10p≈0.9984），top5 完全重合、因果 SNP 排名 1/1，GEMMA null VCE 与 CS 精确一致。
- 新增公共环境模型 GWAS（`--common-env-model`，实验性）：`y = Xb + x_j*beta_j + Z_a*a + Z_c*c + e`，基因组 GRM 遗传效应 + 分组共享环境随机截距；支持 subject 级与 record 级（`--common-env-record-level`）分组文件、PLINK/BGEN/PGEN 输入、MAF/HWE QC、LOCO（含按染色体 VCE 重估）与固定方差分量（`--reml-vg/-vc/-ve`）两种扫描。
- 新增通用 multi-random-effect LMM 层（`multi_random_effect_lmm.{h,cpp}`）与矩阵适配器（GRM/Identity/Custom 三种随机效应核），供公共环境模型复用记录层分数检验。
- 模拟验证（N=120, M=160, 12 组, 360 记录, sigma2a/c/e=0.45/0.25/0.30）：VCE 回收 0.4877/0.3475/0.3227，因果 SNP 为 top hit（p=6.5e-07），LOCO 因果 p=5.7e-10，记录级分组与固定 VC 路径一致；重复力与 RRM 动态回归验证通过。

## v1.1.0 - 首个公开发布候选版（2026-07-13）→ Release Polish（2026-08-23）

- 发布范围收敛为单性状 LMM GWAS、多性状 LMM GWAS 和重复力 GWAS。
- 提供 Linux x86_64 二进制包和可独立构建的 focused 源码包。
- focused 源码包采用单一扁平 `src/`，排除 RRM、PSR、ssGWAS、survival、rare-variant、biobank 等非本版模型源码。
- 发布验证覆盖单性状、多性状、重复力 `exact` 以及 `fast` 对照流程。
- 移除 repeatability 路径中的 `curl`/`std::system` 调试上报，改为本地诊断输出。
- 程序版本由 CMake `PROJECT_VERSION` 统一注入。
- 默认关闭外部 BLAS/LAPACK，以提高跨环境数值复现性。
- 增加源码和二进制归档 SHA-256 校验、引用元数据及发布范围文档。
- 修复重复力 VCE 的 PCG 零右端项处理，支持个体内中心化后求和为零的记录级协变量；所有 VCE 初值失败时不再接受零方差结果。
- 默认 subject-level Woodbury `exact/fast` 路径不再无条件计算完整 GRM 逆矩阵；仅 STCG、full-record VCE 和 refit 路径按需求逆。
- 重复力 null summary 的 `n_subjects` 改为实际进入表型分析的唯一受试对象数，不再报告基因型文件总样本数。

## 1.0.5 - 回溯内部里程碑（2026-05-22 至 2026-06-24）

- 将早期个体均值聚合方案发展为独立重复力模型。
- 增加重复力 `exact`、`fast` 和实验性 `stcg` 路线，并明确以记录层 `exact` 作为正式分析路线。
- 完善多性状 EMmvGWAS-compatible、GEMMA-style 及 LOCO/过滤输入路径。
- 增加单性状 GEMMA accuracy、多性状 reference matrix、重复力 simulation/edge-case/format consistency 等验证脚本。
- 该阶段是开发功能和验证基础设施的集中整合，并非公开版本。

## 1.0.0 - 回溯内部里程碑（2026-04-07）

- 建立独立 CosmicScanner GWAS 代码基线。
- 支持单性状 LMM、GRM 和 LOCO 扫描。
- 提供初始多性状 LMM 联合检验路径。
- 提供基于个体加权均值聚合的早期重复测量分析原型；尚不是后续记录层精确重复力模型。
- 该里程碑仅用于描述代码演进，并非公开版本。

## 开发分支（Phase 6.6，2026-08-17）

- **外部基准（HIBLUP v1.6.0 参考实现）**：`docs/SSGBLUP_EXTERNAL_VALIDATION.md`
  发布 Phase 6.6 外部验证结论，复现脚本 `tests/run_ssgwas_hiblup_benchmark.sh`
  （结果 `tests/results/ssgwas_hiblup_benchmark/`）。
  - 逐一数值核实 HIBLUP 语义：A（Meuwissen-Luo，max diff 0.0）、原始 G
    （M M′/Σ2pq，Su GS）、tuneG（b=(dA−oA)/(dG−oG)，a=dA−b·dG，先行）、
    混合（0.95G*+0.05A11）、H⁻¹=A⁻¹+Gw⁻¹−A11⁻¹、backsolve 用调优未混合 G*⁻¹。
  - **新发现 HIBLUP `--vc-prior` 内部缩放**（N=2000 下有效 λ=0.98×提供值）：
    经 λ 扫描定位；对齐协议为给 HIBLUP 提供 Ve/0.98。校准后 EBV 逐类
    **完全一致**：Pearson 1.00000000 / slope 1.000028 / RMSE 1.17e-5 /
    max|Δ|=5e-5（HIBLUP 输出 6 位小数精度）；四类动物（基因型×表型 2×2）
    全部 Pearson 1.00000000。
  - SNP 效应（`add_a1` vs `EFFECT_ALLELE`）corr 0.99943，OLS slope 0.9447
    （HIBLUP 侧数据相关内部刻度，按规范 §11/12 仅声明定义对齐）；
    GEBV 重建恒等式复核 corr 1.00000000。谱系行序洗牌逐位一致。
- **新 CLI 选项（Argument 兼容，默认行为不变）**：
  - `--ssgwas-tune-g`：backsolve 路径 G 调优（bG+a，以 A22 均值刻度对齐；
    原为 H 核 lmm 路径专用，现共享）。默认关。
  - `--ssgwas-blend-a22`：backsolve 路径混合目标由单位阵改为 A22
    （Gw=(1−α)G*+α·A22，α 取 `--blend-g`=0.05），MME 用 Gw⁻¹、
    **SNP 效果反解用调优未混合 G*⁻¹**（HIBLUP 参考语义，meta `EFFECT_GINV`）。
    默认关（保持 Phase 6.5 已验证单位阵混合约定）。
  - meta 新增 `G_RAW_METHOD/G_RAW_DENOMINATOR`、`G_MEAN_DIAG/OFFDIAG_BEFORE`
    `A22_MEAN_DIAG/OFFDIAG` `G_MEAN_DIAG/OFFDIAG_AFTER` `BLEND_TARGET`
    `EFFECT_GINV`（与既有 `G_TUNED/G_SCALING_A/G_SCALING_B/BLEND_ALPHA` 并列）。
- **就绪包**：`tests/ssgwas_external/blupf90/`（引用描述 + 输入转换 + 运行器 +
  postGSf90 对比协议；本环境无 BLUPF90 可执行，SKIP-safe）。
- **回归**：Phase 6.5 全部门保持绿色（backsolve 全流程、数学红队 gates 1–7、
  SingleStep smoke、CTest 12/12、GEMMA accuracy、RRM CI）。
- 运行时（N=2000/M=2000）：CS backsolve 2.66s wall（MME 1.52s）vs
  HIBLUP mme 0.22s；仅记录不做速度声明。
