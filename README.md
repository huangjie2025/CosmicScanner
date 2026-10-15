# CosmicScanner

**Version: v2.0.0 — finalization in progress; not frozen or published**

Version governance: [R29 decision](docs/R29_VERSION_GOVERNANCE.md).
The default output schema changed since tagged v1.1.3; see the migration notes.

当前公开范围以 [v2.0 scope contract](COSMICSCANNER_V2_SCOPE_FREEZE_PROPOSAL.md)
为唯一入口：A = 限定正式支持契约，B = 实验，C = 内部/开发用途。A 不代表
当前 dirty 工作树已通过最终 build/package 验收。架构与来源见
[architecture baseline](COSMICSCANNER_V2_ARCHITECTURE_BASELINE.md)。

历史候选范围、基准、安装与复现记录（不是当前验收证明）：
[v1 能力矩阵](docs/COSMICSCANNER_V1_CAPABILITY_MATRIX.md)、
[v1 基准摘要](docs/COSMICSCANNER_V1_BENCHMARK_SUMMARY.md)、
[安装](INSTALL.md)、[复现](docs/COSMICSCANNER_REPRODUCIBILITY.md)。
历史 R19–R27 文档中的 1.1.3 标识保留原始二进制来源，不代表当前候选版本。

发布证据（R23）：[最终能力矩阵](docs/COSMICSCANNER_FINAL_CAPABILITY_MATRIX.md)、
[实测扩展性](docs/COSMICSCANNER_SCALABILITY_SUMMARY.md)、
[发布 README](docs/COSMICSCANNER_RELEASE_README.md)。这三个文档是下述
22 个 GWAS / post-GWAS 方法的历史证据入口，不替代 v2.0 的 A/B/C 范围。

CosmicScanner 是以 Gaussian 线性模型和混合模型为核心的 C++ 命令行遗传关联分析软件，提供具有明确支持边界的多性状、重复记录、系谱感知反解及关联后汇总统计分析。模型共享部分基础设施，但没有全模型统一 numerical engine、输入格式或输出 schema；style/compatibility 名称不表示与参考软件完全等价。

详细使用说明请参阅：**[CosmicScanner 中文教程](docs/Tutorial_CN.md)** | **[English Tutorial](docs/Tutorial_EN.md)**

## 统一命令行入口

推荐保留领域通用长参数，使用 `--gwas --model <模型>`；可选 `gwas`/`grm` 入口与旧平铺命令兼容。
`--help` 显示简明帮助，`--help models` / `--help input` 分组说明，`--help-all` 保留完整高级参考。
详见 [中文命令行指南](docs/CLI_GUIDE_CN.md) / [English CLI guide](docs/CLI_GUIDE_EN.md)，包含可运行合成数据、列号规则、预期结果和旧命令迁移。

```bash
cosmicscanner --bfile data --pheno trait.txt --gwas --model lmm --out run
cosmicscanner grm --bfile data --out data_grm
```

## 安装

### 预编译版本（待最终确认与发布）

下列为候选归档命名约定，不表示归档已经发布。旧 focused 打包脚本不适用于
本候选的完整方法范围；详见 R28 文档审计。当前可验证的安装路线为源码构建。

```bash
tar -zxvf CosmicScanner-v2.0.0-linux-x86_64.tar.gz
cd CosmicScanner-v2.0.0-linux-x86_64
bash install.sh
bin/cosmicscanner --help
```

### 源码编译

```bash
git clone --recurse-submodules https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

编译成功后，可执行文件位于 `build/cosmicscanner`。

### 依赖

- CMake >= 3.15
- 支持 C++17 的编译器 (GCC 9+ / Clang 10+)
- OpenMP（编译时必需）、Eigen、zlib 和 zstd
- 外部 BLAS/LAPACK 可由 CMake 自动检测；R23 审计的本地构建使用 Intel MKL，
  不保证不同数值库或线程配置之间逐位一致。

## 快速使用示例

```bash
# 1. 计算 GRM
cosmicscanner --bfile my_data --make-grm --out my_grm

# 2. 运行单性状 GWAS（默认 LOCO，推荐 semi）
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno 1 --gwas --scan-tier semi --out my_results

# 3. 运行多性状 GWAS
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno-list 1,2,3 --gwas --model multivariate --scan-tier approx --reml-maxit 1000 --out my_mv_results
```

## 发布状态

论文纳入范围、支持等级和验证结果是不同维度。A 范围包括限定的 LM/LMM/LOCO、
BED mvLMM、BED repeatability record-level approx/semi、单性状 AD/dominance、
ordinary logistic、burden/SKAT、native post-GWAS，以及 BED/fixed-VC/LDLT
Method 2 reference effects/EBV/window 和显式 exact P。RRM、animal、maternal、
common environment、PSR/GxE、EB、style-methods、APY core-v2/approx 等为 B。
Survival 原型、Method 1 reserved EBV、legacy APY inference、automatic C98/C99
及 BGEN AD research 为 C。完整格式和参数边界见 v2 scope contract。

支持等级不替代 R23 的 `PASS` / `PASS_WITH_LIMITATIONS` /
`NOT_VALIDATED` 历史验证状态。[R23 矩阵](docs/COSMICSCANNER_FINAL_CAPABILITY_MATRIX.md)
保留当时的证据身份；更广泛的
育种模型论文范围与历史证据见
[docs/PAPER_MODEL_SUPPORT_MATRIX.md](docs/PAPER_MODEL_SUPPORT_MATRIX.md) 和
[docs/FINAL_VALIDATION_REPORT.md](docs/FINAL_VALIDATION_REPORT.md)。

## 核心功能

*   **高性能数据读取**：
    *   原生读取 **PLINK 2.0 (`.pgen`)**；支持范围依模型而异，未据此证明百万级样本运行。
    *   支持 **BGEN (`.bgen`)** 和传统 **PLINK 1.9 (`.bed`)** 格式。
    *   采用内存映射（Memory Mapping）与多线程解压技术。

*   **先进的 LMM 引擎**：
    *   **LOCO (Leave-One-Chromosome-Out)**：在支持路径中排除当前染色体的背景 SNP；不保证所有数据都改善功效，也不表示所有模型逐染色体重估全部 VC。
    *   **Fast-LMM 近似**：通过固定方差比（Delta）显著加速全基因组扫描。
    *   **多性状分析 (MvLMM)**：支持多表型联合分析，利用性状间相关性提高统计效力。
    *   **重复测量模型**：支持同一特定个体的多次观测数据分析。

*   **并行计算**：
    *   基于 OpenMP 的多线程扫描与 GRM 构建。
    *   针对现代 CPU 指令集（AVX2/AVX-512）优化的矩阵运算。

### 实测规模证据（MEASURED，NOT ESTIMATED）

R19-E/R20 冻结二进制的 BED matrix-free LRLMM（rank=256）实测：
20k×1M 用时 1,907.261 s、峰值 RSS 6.034 GB；50k×1M 用时 4,655.433 s、
峰值 RSS 13.938 GB，两者均未构造 dense GRM。50k 全量重复结果逐位一致。
这些结果不代表所有模型或当前二进制的性能；没有 dense 50k 精度对照，
也不保证 rank=256 对任意数据足够。详见实测扩展性文档与二进制哈希。

## 快速使用示例

```bash
# 1. 计算 GRM
cosmicscanner --bfile my_data --make-grm --out my_grm

# 2. 运行单性状 GWAS（默认 LOCO，推荐 semi）
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno 1 --gwas --scan-tier semi --out my_results

# 3. 运行多性状 GWAS
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno-list 1,2,3 --gwas --model multivariate --scan-tier approx --reml-maxit 1000 --out my_mv_results
```

## Gaussian LMM 扫描档位

单性状、多性状和重复力 Gaussian LMM 支持统一的 `--scan-tier approx|semi|exact`。三档分别对应 fixed-null、fixed-ratio/profile-scale 和 per-SNP full variance-component refit。详细命令及当前边界见 [docs/SCAN_TIERS_CN.md](docs/SCAN_TIERS_CN.md)。

## 重复力 GWAS

`--repeatability-model` 适用于同一受试对象存在多次观测的长格式表型。正式分析应显式选择记录层引擎和扫描档位：

```bash
cosmicscanner \
  --bfile my_data \
  --pheno repeated_pheno.txt \
  --tpheno 1 --mpheno 2 \
  --qcovar covar_q.txt \
  --covar covar_c.txt \
  --repeatability-model \
  --repeatability-mode exact \
  --gwas \
  --scan-tier semi \
  --out repeated_assoc
```

输入要求：
- `--pheno` 使用长格式，前两列必须是 `FID IID`，后续为表型列。
- 同一 `FID IID` 可出现多行，每一行代表一条重复记录。
- 带时间列时，记录级协变量可按 `FID IID TIME` 对齐；无时间列时按 `FID IID` 对齐。

输出文件：
- `*.assoc.txt`：重复力 GWAS 结果。
- `*.repeatability.null.txt`：null model 的 `Vg/Vpe/Ve/h2/r2` 摘要。

快速模式：
- 使用 `--repeatability-mode fast` 切换到“个体均值聚合 + 加权近似”；`--fast-lmm` 仅保留旧命令兼容。
- fast 属于 B 级聚合近似。auto 可按规模切换 exact/fast/STCG；正式记录层示例须显式指定 exact mode，不能声称默认始终精确。per-SNP exact tier 为 B，限 BED、no-LOCO、无记录权重。
- fast mode 会额外输出 `*.agg.pheno`，并在 `*.repeatability.null.txt` 中将 `scan_mode` 标记为 `fast_repeatability_mean_model`。

工程优化选项：
- `--repeatability-mixed-precision`：在 exact 模式下使用 `Float32` 构建记录层核矩阵，并在最终求解前提升到 `Float64`，以降低中间内存占用。
- `--repeatability-grm-cache <prefix>`：将重复力 exact 模式的 subject-level GRM 缓存为 `Float32` 二进制文件；再次运行时会优先通过 `mmap` 复用缓存，作为后续 `out-of-core` 路线的基础接口。

## 公共环境模型 GWAS（实验性）

`--common-env-model` 适用于同一受试对象的多条记录来自不同环境、需要按环境（而非受试对象）引入随机截距的场景：

```bash
cosmicscanner \
  --bfile my_data \
  --pheno long_pheno.txt \
  --tpheno 2 \
  --common-env-model \
  --common-env-file group.txt \
  --gwas \
  --out ce_gwas
```

模型：`y = Xb + x_j*beta_j + Z_a*a + Z_c*c + e`，其中 `a ~ N(0, K*sigma2_a)`（基因组 GRM，跨受试对象）、`c ~ N(0, I*sigma2_c)`（按组共享的环境随机截距），同组记录共享 `sigma2_c` 相关。

- 表型为长格式 `FID IID [TIME] VALUE`；`--tpheno <col>` 指定值列。
- `--common-env-file <f>`：默认按受试对象分组（每行 `FID IID GROUP`）；加 `--common-env-record-level` 后按记录分组（`FID IID [TIME] GROUP`，需与表型记录顺序一致）。
- `--reml-vg/-vc/-ve <v>`：固定方差分量、跳过 VCE；`--no-loco` 关闭逐染色体 LOCO。
- 输出：`*.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`）与 `*.common_env.null.txt`（`sigma2_a/sigma2_common/sigma2_e/h2/common_environment_fraction`）。

> v2.0 Level B：研究路径；group common environment 与 subject permanent environment 不是同一效应。

## 系谱动物模型 GWAS（Animal Model）

`--animal-model` 使用系谱亲缘关系矩阵 A 描述多基因背景遗传效应：

```
y = Xb + x_j*beta_j + Z*a + e,   a ~ N(0, A*sigma2_a)
```

```bash
cosmicscanner \
  --bfile my_data \
  --pheno my_pheno.txt \
  --animal-model \
  --ped pedigree.txt \
  --gwas \
  --out animal_assoc
```

- 系谱文件列格式：`ID SIRE DAM`，`ID` 为个体标识（需包含全部分析样本），`SIRE`/`DAM` 无亲本时写 `NA`/`.`/`0`。
- `--write-ped-kernel`：将分析子集（去缺失表型后）上的精确 A22 子矩阵导出为 `*.animal_model.k.txt`（空格分隔，全矩阵）与 `*.animal_model.k.id`，可供 GEMMA `-k` 使用。
- `--no-loco`：默认 LOCO 不适用（A 只是环境/近交结构 G 的系谱替代），多染色体数据会打印提示；单染色体可忽略。
- 输出：`*.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`）与 `*.animal_model.null.txt`（`sigma2_a/sigma2_e/h2`）。

该模型与 ssGWAS 共享同一 `PedigreeRelationship` 系谱解析与 `A⁻¹`/`A22` 构建代码（`buildA22` 经由 A⁻¹ 的 Schur complement 求解，不物化稠密 A）。两种模型给出完全等价的 null 方差分量（ssGWAS `--h2-tau 0 --h2-omega 0 --no-loco` 时 β/SE/-log10P 与 animal model 逐 SNP 一致）。详情见 [docs/ANIMAL_MODEL_GWAS.md](docs/ANIMAL_MODEL_GWAS.md)。

## A/G/H 矩阵

```bash
cosmicscanner --bfile my_data --make-agh \
  --ped pedigree.txt --genotyped-list genotyped_ids.txt --out results/agh
```

输出 `Ainv`、`G`、`Hinv` 的 lower-triangle Float32 binary 与 ID 文件，以及 `*.agh.meta.txt`。ssGWAS 中加入 `--write-agh` 可同时保存这些矩阵。

## 随机回归 / 动态 GWAS

`--long-pheno` 与 `--rrm-order` 适用于纵向表型或时间梯度表型。模型使用 Legendre 基函数表示动态 SNP 效应，并输出联合检验与时间点效应曲线：

```bash
cosmicscanner \
  --bfile my_data \
  --long-pheno longitudinal_pheno.txt \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-order 1 \
  --rrm-curve-grid 21 \
  --rrm-top 100 \
  --rrm-ai-iter -1 \
  --rrm-sig-threshold 0.05 \
  --gwas \
  --out dynamic_assoc
```

主要输出文件：
- `*.rrm_joint.assoc.txt`：基于 observation-level RRM 的多自由度联合检验结果。
- `*.rrm.dynamic_curve.txt`：top SNP 在时间网格上的动态效应、标准误和时间点 P 值，使用联合模型的完整协方差矩阵计算。
- `*.rrm.dynamic_summary.txt`：top SNP 的最佳作用时间、起止效应、最大效应、轨迹类型和基于 `--rrm-sig-threshold` 的关键时间窗口。
- `*.rrmcoef.txt`：个体层 Legendre 系数。

性能选项：
- `--rrm-ai-iter -1`：默认自动模式；小样本使用 AI refinement，大样本使用 moment-initialized fast null。
- `--rrm-ai-iter 0`：强制 fast null。
- `--rrm-ai-iter <n>`：强制执行 `n` 次 AI refinement。

非平衡 longitudinal RRM：

```bash
cosmicscanner \
  --bfile my_data \
  --pheno longitudinal_pheno.txt \
  --tpheno 1 \
  --mpheno 2 \
  --rrm-unbalance \
  --rrm-gwas-method record \
  --rrm-basis legendre \
  --rrm-order 3 \
  --tclass Sex_class \
  --rrm-ai-iter 20 \
  --no-loco \
  --out unbalanced_rrm
```

该模式面向每个个体观测时间点不完全一致或记录数不同的 long-format 表型。`--rrm-gwas-method` 选择两条具有不同统计语义的扫描路径（文档详见 `docs/RRM_GWAS_METHOD_AUDIT.md`）：

- `--rrm-gwas-method record`：**记录级 time-varying SNP 效应**（Method A）——SNP 固定设计为 `x_{subject(r)} × φ_k(t_r)`（record×basis），联合 Score/Wald `score'M⁺score ~ χ²(df=rank)` 显式诊断 `STATUS ∈ {OK, RANK_DEFICIENT, ILL_CONDITIONED, NUMERICAL_FAILURE}`。真正的 `x⊗Φ` 固定效应检验，与平衡 joint 路径同构（模拟等价 Pearson ≥ 0.99999；A–G 模拟套件与 A/B 对比见 `tests/run_test_rrm_method_a_simulation.sh`、`tests/run_rrm_method_ab_comparison.sh`）。软件状态为 **Experimental validated**。
- `--rrm-gwas-method coefficient`（**默认**，`--rrm-unbalance` 单独使用与之等价，backward compatible）：内部先按个体拟合 Legendre 随机回归系数，再用 BLUP 系数上的 Kronecker multivariate LMM 扫描 SNP（Method B，coefficient-based two-step），避免构建大规模 observation-level 协方差矩阵；`--rrm-compatible gmat-unbalance` 启用 GMAT unbalanced longwas 兼容口径（全集 Pearson=0.9968、Spearman=0.9959、top10 全一致）。

`--tclass` 指定时间依赖分类协变量，`--tcovar` 指定时间依赖连续协变量；LOCO 模式支持 PLINK/PGEN/BGEN，但每条染色体会重新估计 multivariate null model，速度慢于 `--no-loco`。

RRM 基函数和 LOCO 选项：
- `--rrm-basis legendre`：默认 Legendre 多项式基函数，用于和 GMA/GMAT 对齐。
- `--rrm-basis bspline`：open-uniform B-spline 基函数；`--rrm-order n` 表示输出 `n+1` 个 basis 系数。
- `--rrm-loco-global-null`：快速近似 LOCO，只按染色体分组扫描 SNP，复用 global BLUP null model；默认不开启，默认 LOCO 仍逐染色体重估 null model。

非平衡模式输出：
- `*.rrm_unbal.assoc.txt`：非平衡 RRM 联合检验结果，包含 `JOINT_CHISQ`、`JOINT_P`、`BETA_L0...` 和对应标准误（Method A 为 `x⊗Φ` 基系数与 `STATUS/DF` 诊断；Method B 为 BLUP 系数空间多性状 Wald）。
- `*.rrm_unbal.null.txt`：BLUP coefficient LMM 的 null 方差组分矩阵，记录遗传协方差 `K_*` 和残差/永久环境合并协方差 `PE_*`。
- `*.rrm_unbal.dynamic_curve.txt`：top SNP 在时间网格上的动态效应曲线。
- `*.rrm_unbal.dynamic_summary.txt`：top SNP 的关键时间、显著窗口和轨迹分类。

GMAT mouse_long 非平衡数据验证结果（Method B，coefficient，`--rrm-compatible gmat-unbalance`）：1212 个体、11833 SNP、19392 条纵向观测、16 个时间点、`--rrm-order 3`。CosmicScanner BLUP-based RRM 在 `--no-loco` 下运行约 3.6 秒、峰值内存约 164 MB；与 GMAT `longwas --unbalance` 对齐达到全集 `Pearson r(-log10P)=0.9968`、`Spearman rho(P)=0.9959`，常染色体 `Pearson=0.9969`、`Spearman=0.9964`，top SNP 和 top10 均一致（`UNC18848064`）。旧结果全集 `Pearson=0.9470` 的主因是 `Sex_class` 时间依赖分类协变量未进入 coefficient null 的 subject-level fixed design，且分类编码存在基线错误；修正后不需要把差异归因于隐藏 imputation mismatch。

Method A（record）的历史文档报告平衡 joint 路径等价和 A–G 模拟结果。
但 R23 checkout 缺少 wrapper 引用的 `tests/run_rrm_snp_trajectory_simulation.py`，
因此这些 A–G 数值不能由当前 wrapper 复现，也不作为 R23 的新增验证证据。
R23 将 RRM 标为 `PASS_WITH_LIMITATIONS`；历史外部 GMAT 比较也不能视作
所有 RRM 路径的通用外部验证。

GMA-compatible balanced fixreg 验证模式：

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

该模式用于复现 GMA `balance_longwas_fixreg` 的 balanced fixed-regression Wald scan，当前限定 PLINK `--bfile` 输入，并读取同前缀 `*.grm` 以及 `--rrm-compatible-varcom` 指定的 GMA `.var` 方差分量文件。输出为 `*.rrm.gma_fixreg.assoc.txt`。该兼容模式不替代默认 RRM dynamic GWAS；默认路径仍负责动态效应曲线、关键窗口和轨迹分类。

可视化：

```bash
python3 tools/plot_rrm_dynamic.py \
  --curve dynamic_assoc.rrm.dynamic_curve.txt \
  --summary dynamic_assoc.rrm.dynamic_summary.txt \
  --comparison comparison.tsv \
  --top 5 \
  --outdir figures
```

该脚本会输出动态效应曲线图、效应热图、方法排名对比图和绘图用 summary 表。

验证脚本：
- `tests/run_test_rrm.sh`：小样本 RRM smoke/regression test。
- `tests/run_rrm_dynamic_simulation.sh`：模拟 stable、sign-switch、early-specific、late-specific 和 transient 五类 QTL，验证 causal SNP 检出、效应曲线和轨迹分类输出；同时运行平均表型 GWAS、端点逐时间 GWAS 和多性状 GWAS，输出 causal SNP 排名对比表 `comparison.tsv`。
- `tests/run_rrm_dynamic_simulation.sh --loco-smoke`：在多染色体模拟数据上额外运行 RRM dynamic LOCO smoke test。
- `tests/run_rrm_dynamic_simulation.sh --pgen-smoke`：将模拟 PLINK 数据转换为 PGEN 后运行 RRM dynamic smoke test，验证 PGEN 输入路径。
- `tests/run_rrm_dynamic_simulation.sh --bgen-smoke`：将模拟 PLINK 数据转换为 BGEN 后运行 RRM dynamic smoke test，验证 BGEN 输入路径。
- `tests/run_rrm_dynamic_simulation.sh --order-sweep`：额外比较 `--rrm-order 1/2/3`，输出 `rrm_order_sweep.tsv`，用于阶数敏感性检查。
- `tests/run_rrm_dynamic_benchmark.sh`：中等规模 RRM dynamic 性能 benchmark，默认 `N=500, M=3000, TIMES=5`，输出 `summary.tsv`、`summary.md` 和 runtime/RSS 记录。
- `tests/run_rrm_bspline_benchmark.sh`：同一模拟数据上比较 `--rrm-basis legendre` 与 `--rrm-basis bspline` 的 runtime、内存、causal SNP 排名和轨迹分类。
- `tests/run_rrm_method_comparison.sh`：多 replicate 方法比较，默认 `N=120, M=300, TIMES=5, REPS=5`，比较 RRM fast、RRM exact、平均表型 GWAS、端点 GWAS 和多性状 GWAS 的 power、null FPR、runtime 和动态曲线 RMSE。
- `tests/run_gma_compatible_fixreg_alignment.sh`：使用 GMA 官方 balanced case 子集，比较 GMA `balance_longwas_fixreg` 与 CosmicScanner `--rrm-compatible gma-fixreg`；默认 `N=120, M=20` 当前 `Pearson_log10P=1.0`、`Spearman_P=1.0`、最大 P 值差约 `3e-9`。
- `tests/run_gmat_unbalanced_alignment.sh`：使用 GMAT mouse_long 非平衡数据，比较 GMAT `longwas --unbalance` 与 CosmicScanner `--rrm-compatible gmat-unbalance`（Method B）；要求全集和常染色体 `Pearson_log10P >= 0.99`、`Spearman_P >= 0.99` 且常染色体 top SNP 一致。
- `tests/run_rrm_balance_record_equivalence.sh`：平衡显式 joint 路径 vs 非平衡 record Method A 的等价回归（N=200/M=300/T=8：BETA/SE/STAT/log10P Pearson ≥ 0.99999）。
- `tests/run_test_rrm_method_a_simulation.sh`：Method A 模拟套件 A–G（null FPR/膨胀、恒定/线性/二次/sign-switch 轨迹 power、rank 缺陷诊断）。
- `tests/run_rrm_method_ab_comparison.sh`：同数据上 Method A（record）vs Method B（coefficient）对比，验证 time-varying power 与二者统计差异。
- `tests/run_gemma_emmvgwas_compat_smoke.sh`：本地 GEMMA/EMmvGWAS 兼容 smoke，不依赖外部 GEMMA/R；覆盖 `-lm/-lmm` GEMMA mode、EMmvGWAS `vc_estimation/full/semi` alias、summary-stat 多 category `--gs/--gq/--vc/--ci`、BSLMM/DAP 参数兼容。
- `tests/run_gwas_core_regression.sh --with-rrm`：在 GWAS 核心回归中追加 RRM dynamic 验证；`--rrm-only` 可只运行动态 GWAS smoke 和模拟验证。
- `tests/run_gwas_core_regression.sh --compat-smoke-only`：只运行本地 GEMMA/EMmvGWAS 兼容 smoke，适合快速验收新增兼容层。
- `ctest --test-dir build -R gemma_emmvgwas_compat_smoke --output-on-failure`：同一兼容 smoke 的 CTest 入口，输出 manifest 位于 `build/test-results/gemma_emmvgwas_compat_smoke/compat_smoke_manifest.tsv`。
- `cmake --build build --target check_gemma_compat`：构建 `cosmicscanner` 后立即运行 GEMMA/EMmvGWAS 兼容 CTest，是当前最短落地验收命令。
- `tests/run_gemma_bslmm_reference_validation.sh` 或 `cmake --build build --target check_gemma_reference`：在本机存在 GEMMA 时运行外部 GEMMA `-bslmm 2` reference validation，比较 CS/GEMMA 的 `.param.txt` alpha 和 `.bv.txt` breeding value 相关性。
- `tests/run_gemma_multivariate_validation.sh` 或 `cmake --build build --target check_gemma_mv_reference`：运行外部 GEMMA 多性状 mvLMM reference validation，并在安装 EMmvGWAS R 包时同步比较 EMmvGWAS `gwas_all`。该测试对 GEMMA full/exact、GEMMA mode 1/2/3、exact CRT、GxE/full GxE、EMmvGWAS full/semi G/R 路线设置硬性相关性和输出形状阈值。
- `cmake --build build --target check_gemma_all`：聚合运行 GEMMA/EMmvGWAS 兼容 smoke、GEMMA BSLMM2 reference 和多性状 GEMMA/EMmvGWAS reference validation。
- `cmake --build build --target check_eb_gwas`：运行 EB-GWAS quick suite，覆盖单性状、重复力、多性状 EB smoke、validation matrix 和 standard-vs-EB report。
- `cmake --build build --target check_eb_gwas_full`：运行 EB-GWAS full suite，额外覆盖 LOCO、calibration、calibration grid、PGEN/BGEN consistency 和 performance benchmark。
- `cmake --build build --target check_ssgwas`：运行 ssGWAS A22-corrected H⁻¹ smoke，检查 `--single-step --ped --genotyped-list`、null summary、A22 correction 标记和 H⁻¹ 对称性诊断。
- `tests/run_test_animal_model.sh`：系谱动物模型模拟验证套件（causal 检出/LOCO 提示、ssGWAS tau=0 等价、系谱乱序不变性、σ²a=0、高遗传力、不完整系谱、无祖先系谱）。
- `tests/run_animal_model_gemma_validation.sh`：外部 GEMMA `-k` A22 核基准——将 `--write-ped-kernel` 导出的精确 A22 交给 GEMMA `-lmm 4`，比较 β/SE/-log10P（无 GEMMA 时报告 SKIP）。
- `docs/GEMMA_EMMVGWAS_COSMICSCANNER_SUPPORT_MATRIX_CN.md`：GEMMA、EMmvGWAS 与 CosmicScanner 三方功能支撑矩阵，明确完全支持、兼容支持、近似支持和条件验证边界。

方法说明：
- `docs/RRM_GWAS_VALIDATION_CN.md`：纵向/动态性状 GWAS 的模型、统计量、输出、模拟验证和可视化说明。

## 开发状态

### 数据与引擎基础设施
- [x] PLINK 1.9/2.0 & BGEN 读取支持
- [x] 高效 GRM 构建 (Dense & Sparse)
- [x] Exact LMM & Fast-LMM 扫描引擎
- [x] LOCO 策略实现
- [x] IDUL/AIDUL 加速 null 求解器 (`--lmm-null-solver brent|idul|idul_plus|aidul`)
- [x] REGENIE 风格 Biobank 两步法 (`--biobank-mode`)：PGEN keep-index 路径已增加子集适配。
- [x] GEMMA-style eigen 工具（`--gemma-eigen`，输出 `.eigenD.txt/.eigenU.txt`）
- [x] GEMMA-style summary-stat 工具（`--gs/-gs` 输出 `.S.txt/.size.txt/.snps.txt`，含 `--cat/--mcat/--wcat/--wsnp` 多 component S；`--gq/-gq --beta` 支持 `Z`、`BETA+SE`、`CHISQ` 或 `P`，输出 `.Vq.txt/.q.txt/.size.txt`；`--study/--ref` 和 `--mstudy/--mref --vc/--ci` 输出 summary VC/CI；`--beta/--cor --vc` 支持 GEMMA `.cor.txt` 直接 summary VC）
- [x] GEMMA-style SNP correlation 工具（`--calccor --windowbp/--windowcm`，输出 `.cor.txt`）
- [x] GEMMA-style variance component 工具（`--vc 1|2|3`，输出 `.vc.txt/.vc.summary.txt`）
- [~] 输入适配层 MVP（csv/tsv/txt → prepared input，支持 repeatability/rrm/lmm/multivariate）

### GWAS 模型套件（含 GEMMA/EMmvGWAS 兼容路线）
- [x] GEMMA-style 普通线性模型 GWAS（`--lm --lm-mode 1|2|3|4`，对应 GEMMA `-lm` Wald/LRT/Score/all；输出 `.lm.assoc.txt`）
- [x] 单性状 LMM GWAS（默认 `--gwas`；`--lmm-gemma-mode 1|2|3|4` 对应 GEMMA `-lmm` Wald/LRT/Score/all 输出；`--lmm-gemma-refit` 启用逐 SNP alternative profile refit）
- [x] GEMMA-style prediction（`--predict-effect <assoc/effect.tsv> --predict-score <out>`，基于 SNP effect 表输出个体预测值）
- [x] GEMMA-style multi-locus 输出（`--bslmm 2` ridge/GBLUP；`--bslmm 1/3` deterministic sparse+ridge/probit approximation；`--bslmm 4/5` DAP BF/HYP 与 fine-mapping summary 近似；兼容 `-hmin/-hmax/-rmin/-rmax/-pmin/-pmax/-gmean/-hscale/-rscale/-pscale/-mh/-smin/-smax/-bf/-hyp`；输出 `.param.txt/.bv.txt/.hyp.txt` 或 `.bf.txt/.gamma.txt/.coef.txt`）
- [x] 多性状 GWAS（`--multivariate`，EMMV、EMmvGWAS-compatible `vc_estimation/full_gwas/semi_gwas` alias、GEMMA mode 1/2/3/4、GEMMA-full、GEMMA-exact-style、GEMMA-GxE、hybrid、EB）
- [x] 重复力模型 GWAS（`--repeatability-model`，exact/fast/STCG/refit-validation 四种策略）
- [x] 随机回归 RRM GWAS（`--long-pheno` + `--rrm-order`，Legendre/B-spline 基函数，平衡/非平衡纵向）
- [~] 母体校正 GWAS / maternal-adjusted GWAS（`--maternal-model`，ad/am/adm/c/e 多组分 VCE；SNP 扫描为个体自身基因型的 direct 固定效应，maternal 通路属背景随机效应）；并支持模型专属 SNP 检验 `--maternal-snp-test direct|maternal|joint`（direct=现状；maternal=母亲基因型单列检验；joint=direct+maternal 2-df Wald，含 rank/条件数 COLLINEAR 诊断；缺失 dam 基因型记录显式过滤并在输出标注 N）。软件状态为 `Experimental validated`，证据为 10 场景 simulation suite。
- [~] H-matrix single-step GWAS（Method 1，`--single-step --single-step-method lmm`）：以 H 为背景协方差的 fixed-SNP LMM 关联，输出 `BETA/SE/P`；`Experimental validated`。
- [~] ssGBLUP back-solving GWAS（Method 2，`--single-step --single-step-method backsolve`）：全系谱 MME → EBV/GEBV → `ĝ=Z'G⁻¹â_g` → SNP effects/SNP variance/滑动窗口 %Var；`Experimental validated`。该标记效应与 Method 1 的固定 SNP 系数不是同一估计量。
- [x] exact back-solved SNP P-value（`--ssgwas-pvalue-method exact`）：exact selected-inverse uncertainty；`Stable`。
- [~] APY SNP P-values（`--ssgwas-pvalue-method apy-exact|apy-approx`）：`apy-exact` 为 `Experimental validated` scalable APY uncertainty path；`apy-approx` 为 `Experimental approximate`，`STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`，raw P-values only（no genomic control/rescale），且 scalability 仅指 uncertainty layer，不称 calibrated P-value。
- [x] 系谱动物模型 GWAS（`--animal-model`，A 矩阵，共享 ssGWAS 的 `A⁻¹`/`A22` 构建，GEMMA `-k` A22 核外部验证通过）
- [x] 阈值/Probit 模型（`--threshold-model --link probit|logit`）
- [x] 计数/二项 PSR-GWAS（`--count-model --family poisson|binomial`，PSR/PQL + weighted LMM）
- [~] 序数 PSR-GWAS（`--ordinal-model`，cumulative-binary + stacked-psr，dense-zkz + grid-reml VC 估计；当前默认 cumulative-binary 仍按 experimental 口径发布）
- [x] 生存 Cox 模型（`--survival-model`，Cox PH + 鞅残差 LMM）
- [x] GxE 反应范数模型（`--gxe-model`，Vg0+Vg1 双 VC，2-df 联合 Wald）
- [~] 公共环境模型 GWAS（`--common-env-model`，subject/record 两级分组，VCE+固定 VC 两种扫描；实验性）

### 检验与扩展
- [x] 罕见变异检验（SKAT + Burden，beta 权重，混合卡方 p 值）
- [x] SPA 鞍点近似（二元性状，CGF 求解）
- [~] 经验贝叶斯 EB-GWAS（单性状/重复力/多性状辅助路线；正式报告应同时保留标准 LMM/mvLMM/repeatability 主结果）
- [~] PSR-GWAS v0.5 序数 stacked-psr 落地（Wald/pseudo-score 校准通过；仍需按近似模型报告 convergence/warning）

### 当前开发焦点
- PSR-GWAS v0.5 多分支整合收尾（大 N 校准、score test 诊断已完成）
- Stage 2C AI-REML 优化与大规模 N 支持规划中
- 重复力/RRM PSR-GWAS 与 exact GLMM 为后续路线图

## Benchmark

CosmicScanner 与外部参考工具的验证结果：

| Comparison | Metric | Value |
|---|---|---|
| CS vs GEMMA (single-trait LMM) | Pearson r(-log10P) | > 0.999 |
| CS vs EMmvGWAS (multivariate) | Pearson r(-log10P) | > 0.999 |
| CS vs RepeatABEL (repeatability) | Pearson r(-log10P) | > 0.999 |
| CS vs GMAT (unbalanced RRM, autosome) | Pearson r(-log10P) | 0.952 |
| CS vs GMA (balanced fixreg) | Pearson r(-log10P) | 1.000 |
| Repeatability speedup vs RepeatABEL | Wall time | 10–20× |
| Repeatability memory reduction vs RepeatABEL | Peak RSS | ~10× |
| BLUP-based RRM speedup vs observation-level | Wall time | 213× |
| BLUP-based RRM memory reduction | Peak RSS | 93× |

详细验证报告见 [docs/FINAL_VALIDATION_REPORT.md](docs/FINAL_VALIDATION_REPORT.md)。

## 示例数据

`examples/mouse_clean/` 提供最小化示例数据集（PLINK 格式，1197 个体），可用于快速验证安装：

```bash
# 单性状 GWAS
cosmicscanner --bfile examples/mouse_clean/mouse_clean \
  --pheno examples/mouse_clean/mouse_clean.pheno \
  --mpheno 1 --gwas --no-loco --out example_results

# 查看结果
head example_results.assoc.txt
```

详见 [examples/mouse_clean/README.md](examples/mouse_clean/README.md)。

## Citation

如果 CosmicScanner 对您的研究有帮助，请引用：

> Huang J, et al. (2026). CosmicScanner: a unified C++ framework for genome-wide association analysis across diverse complex-trait architectures. *[Journal TBD]*.

BibTeX:

```bibtex
@article{huang2026cosmicscanner,
  title={CosmicScanner: a unified C++ framework for genome-wide association analysis across diverse complex-trait architectures},
  author={Huang, Jie and others},
  journal={Journal TBD},
  year={2026}
}
```

## 许可证

本项目采用 [MIT License](LICENSE) 许可证。

## 相关文档
- 中文教程：[docs/Tutorial_CN.md](docs/Tutorial_CN.md)
- PSR-GWAS 路线图：[docs/PSR_GWAS_ROADMAP.md](docs/PSR_GWAS_ROADMAP.md)
- PSR-GWAS 验证报告：[docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md](docs/PSR_GWAS_VALIDATION_REPORT_V0.5.md)
- 序数 PSR 教程：[docs/tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md](docs/tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md)
- 重复力 GWAS 发布说明：[docs/REPEATABILITY_GWAS_RELEASE_CN.md](docs/REPEATABILITY_GWAS_RELEASE_CN.md)
- RRM 方法与验证：[docs/RRM_GWAS_VALIDATION_CN.md](docs/RRM_GWAS_VALIDATION_CN.md)
- 输入适配层规范：[docs/INPUT_ADAPTATION_SPEC_CN.md](docs/INPUT_ADAPTATION_SPEC_CN.md)
- IDUL/AIDUL 加速路线图：[docs/IDUL_AIDUL_ROADMAP.md](docs/IDUL_AIDUL_ROADMAP.md)
- EB-GWAS 模块：[docs/EB_GWAS_MODULE_CN.md](docs/EB_GWAS_MODULE_CN.md)
- 新模型设计：[docs/models/NEW_MODELS_DESIGN.md](docs/models/NEW_MODELS_DESIGN.md)
- GEMMA/EMmvGWAS/CS 支撑矩阵：[docs/GEMMA_EMMVGWAS_COSMICSCANNER_SUPPORT_MATRIX_CN.md](docs/GEMMA_EMMVGWAS_COSMICSCANNER_SUPPORT_MATRIX_CN.md)
- GWAS 发布收口验收记录：[docs/GWAS_RELEASE_CLOSEOUT_CN.md](docs/GWAS_RELEASE_CLOSEOUT_CN.md)
- 公共环境模型开发总结：[docs/COMMON_ENVIRONMENT_GWAS_DEV_SUMMARY_CN.md](docs/COMMON_ENVIRONMENT_GWAS_DEV_SUMMARY_CN.md)
- 系谱动物模型 GWAS：[docs/ANIMAL_MODEL_GWAS.md](docs/ANIMAL_MODEL_GWAS.md)
- 支持模式总表：[docs/SUPPORTED_MODES_CN.md](docs/SUPPORTED_MODES_CN.md)

## Association R2 (Experimental)

[Association R2 usage and API guide](docs/ASSOCIATION_R2.md) documents the
implemented opt-in BED PCA, physical-subset QC, window LD pruning/clumping,
nonadditive BED LM (`--lm --lm-mode 1`), and experimental rare-variant extensions.
All new R2 features are experimental; existing legacy model maturity labels,
including the Stable models listed above, are unchanged. This is not a claim of
full IASQTL coverage or that all regression suites pass.

- PCA/QC/LD requires `--bfile` and saves a persistent `<out>.qc.*` subset without
  overwriting the input or existing genetics outputs. Integrated QC downstream
  support is Gaussian-only and rejects single-step and unsafe auxiliary inputs.
- PCA outputs unit GRM eigenvectors, not projected genotype scores. The default
  GRM is VanRaden; degenerate eigenspaces require subspace-based comparisons.
- LD retains dense sample-by-variant storage with bounded genomic comparisons
  and a conservative 512 MiB workspace guard, not a biobank-ready memory model.
- `--skat-o` is a conservative rho-grid omnibus, not canonical SKAT-O;
  `--acat-v` uses a Cauchy tail approximation under general LD, with no ultra-rare
  collapsing. Read `.rare.txt` as set statistics, not SNP effects, and retain the
  new `.rare.meta.tsv` method metadata; exact schemas are in the guide.

Focused verification and timing entrypoints (not a statement of test results):

```bash
ctest --test-dir build -L association-r2 --output-on-failure
cmake --build build --target benchmark_association_r2
```

R2 did not implement Conditional GWAS, BLINK, or FarmCPU. Their experimental R3
implementations and narrower support contracts are documented below.

## Association R3: Multi-Locus and Conditional GWAS (Experimental)

R3 adds native C++ FarmCPU-style, BLINK-style, and conditional Gaussian
association routes. These implementations are opt-in and do not change legacy
LMM, mvLMM, RRM, repeatability, ssGWAS, threshold, or survival defaults.

```bash
# FarmCPU-style: global LMM ranking, LD/BIC pseudo-QTN selection, fixed-effect iterations
build/cosmicscanner --bfile data --pheno phenotype.txt --farmcpu \
  --farmcpu-max-iterations 5 --farmcpu-pseudo-qtn 10 --out farmcpu

# BLINK-style: LD filtering, BIC pseudo-QTN selection, fixed-effect iterations
build/cosmicscanner --bfile data --pheno phenotype.txt --blink --out blink

# Conditional scan; condition.txt contains one SNP ID per line
build/cosmicscanner --bfile data --pheno phenotype.txt \
  --condition-snp condition.txt --out conditional
```

FarmCPU and BLINK write `<out>.cosmic.assoc`; conditional analysis writes
`<out>.conditional.assoc`. All routes also write `<out>.multilocus.meta.tsv`.
The `MODEL` column is `farmcpu`, `blink`, or `conditional`.

All three routes are **Experimental**. They currently materialize an additive
sample-by-marker matrix, use Gaussian fixed-effect iterations, and are guarded by
a 1 GiB genotype estimate. They are not binary-trait, LOCO, or biobank-scale
implementations and are not claimed to be numerically identical to GAPIT or other
reference packages. See [docs/ASSOCIATION_R3.md](docs/ASSOCIATION_R3.md).

```bash
ctest --test-dir build -L association-r3 --output-on-failure
cmake --build build --target benchmark_association_r3
```

## Association R4-A: Firth Logistic and Conditional Joint (Experimental)

R4-A adds a separation-safe binary-trait scan and an individual-level,
LD/BIC-guided extension of the R3 conditional scan:

```bash
# Phenotype must be coded exactly 0/1 and contain cases and controls
build/cosmicscanner --bfile data --pheno binary.txt \
  --firth-logistic --out firth

# leads.txt contains candidate lead SNP IDs, one per line
build/cosmicscanner --bfile data --pheno phenotype.txt \
  --condition-snp leads.txt --conditional-joint --out joint
```

Firth logistic writes `<out>.cosmic.assoc` and
`<out>.generalized.meta.tsv`. Conditional joint writes
`<out>.conditional_joint.assoc` and `<out>.multilocus.meta.tsv`. Both use the
common association schema shown above. `PSEUDO_QTN=1` marks LD/BIC-selected
joint leads.

Both routes are native C++ and **Experimental**. Firth reports a Wald P-value
from the penalized fit; conditional joint uses individual-level dosages rather
than external summary statistics and is not GCTA-COJO equivalence. Both
materialize a guarded sample-by-marker matrix and are intended for small to
medium datasets. See [docs/ASSOCIATION_R4A.md](docs/ASSOCIATION_R4A.md).

```bash
ctest --test-dir build -L association-r4a --output-on-failure
cmake --build build --target benchmark_association_r4a
```

## Biobank Binary GWAS: SAIGE-Style (Experimental)

R5-A adds a native binary mixed-model score workflow for streaming large marker
sets. It is selected independently from the existing linear REGENIE-style
biobank route:

```bash
build/cosmicscanner --bfile data --pheno binary.txt \
  --saige-style --biobank-mode --sparse-cutoff 0.05 --out binary_gwas
```

The phenotype must be coded exactly `0/1`. The workflow fits one logistic mixed
null model with covariates and a thresholded sparse GRM, estimates a genetic
variance component, calibrates an exact-to-naive score variance ratio on a
small marker subset, then streams BED/PGEN/BGEN variants. Markers passing the
normal trigger are calibrated with the existing SPA implementation.

Results use the common `<out>.cosmic.assoc` schema with
`MODEL=SAIGE_STYLE`. Null-model, sparse-solver, variance-ratio, case/control,
and SPA metadata are written to `<out>.saige_style.meta.tsv`.

This implementation is **SAIGE-style**, not SAIGE and not numerical or
algorithmic equivalence. R5-A avoids per-marker mixed-model solves and streams
the association pass, but its first release still constructs a dense GRM before
thresholding and uses sparse direct LDLT for the null model. It is therefore
marker-scale capable but not yet validated for biobank-scale sample counts. See
[docs/ASSOCIATION_R5A.md](docs/ASSOCIATION_R5A.md).

```bash
ctest --test-dir build -L association-r5a --output-on-failure
cmake --build build --target benchmark_association_r5a
```
