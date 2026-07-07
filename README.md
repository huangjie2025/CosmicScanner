# CosmicScanner

CosmicScanner 是 Cosmic 基因组分析套件的全基因组关联分析（GWAS）模块。它专为大规模生物样本库设计，支持高性能并行计算和现代基因组数据格式。

详细使用说明请参阅：**[CosmicScanner 中文教程](docs/Tutorial_CN.md)**

## 核心功能

*   **高性能数据读取**：
    *   原生支持 **PLINK 2.0 (`.pgen`)**，高效处理百万级样本。
    *   支持 **BGEN (`.bgen`)** 和传统 **PLINK 1.9 (`.bed`)** 格式。
    *   采用内存映射（Memory Mapping）与多线程解压技术。

*   **先进的 LMM 引擎**：
    *   **LOCO (Leave-One-Chromosome-Out)**：自动剔除本染色体效应，避免近端污染，提升检验效力。
    *   **Fast-LMM 近似**：通过固定方差比（Delta）显著加速全基因组扫描。
    *   **多性状分析 (MvLMM)**：支持多表型联合分析，利用性状间相关性提高统计效力。
    *   **重复测量模型**：支持同一特定个体的多次观测数据分析。

*   **并行计算**：
    *   基于 OpenMP 的多线程扫描与 GRM 构建。
    *   针对现代 CPU 指令集（AVX2/AVX-512）优化的矩阵运算。

## 快速使用示例

```bash
# 1. 计算 GRM
cosmicscanner --bfile my_data --make-grm --out my_grm

# 2. 运行 GWAS (使用 LOCO)
cosmicscanner --bfile my_data --pheno my_pheno.txt --gwas --out my_results

# 3. 运行多性状 GWAS
cosmicscanner --pfile my_ukb_data --pheno my_pheno.txt --mpheno-list 1,2,3 --multivariate --out my_mv_results
```

## 重复力 GWAS

`--repeatability-model` 适用于同一受试对象存在多次观测的长格式表型。当前默认实现采用记录层精确重复力模型：

```bash
cosmicscanner \
  --bfile my_data \
  --pheno repeated_pheno.txt \
  --mpheno 1 \
  --qcovar covar_q.txt \
  --covar covar_c.txt \
  --repeatability-model \
  --gwas \
  --out repeated_assoc
```

输入要求：
- `--pheno` 使用长格式，前两列必须是 `FID IID`，后续为表型列。
- 同一 `FID IID` 可出现多行，每一行代表一条重复记录。
- 当前协变量文件仍按个体对齐，因此稳定支持时不变协变量；时变协变量暂未纳入本版扫描。

输出文件：
- `*.assoc.txt`：重复力 GWAS 结果。
- `*.repeatability.null.txt`：null model 的 `Vg/Vpe/Ve/h2/r2` 摘要。

快速模式：
- 追加 `--fast-lmm` 后，重复力 GWAS 会切换到“个体均值聚合 + 加权近似”的 fast mode。
- 该模式更快，但对不平衡重复记录属于近似方案；默认精确模式仍推荐用于基准分析。
- fast mode 会额外输出 `*.agg.pheno`，并在 `*.repeatability.null.txt` 中将 `scan_mode` 标记为 `fast_repeatability_mean_model`。

工程优化选项：
- `--repeatability-mixed-precision`：在 exact 模式下使用 `Float32` 构建记录层核矩阵，并在最终求解前提升到 `Float64`，以降低中间内存占用。
- `--repeatability-grm-cache <prefix>`：将重复力 exact 模式的 subject-level GRM 缓存为 `Float32` 二进制文件；再次运行时会优先通过 `mmap` 复用缓存，作为后续 `out-of-core` 路线的基础接口。

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
  --rrm-basis legendre \
  --rrm-order 3 \
  --tclass Sex_class \
  --rrm-ai-iter 20 \
  --no-loco \
  --out unbalanced_rrm
```

该模式面向每个个体观测时间点不完全一致或记录数不同的 long-format 表型。内部先按个体拟合 Legendre 随机回归系数，再用 BLUP 系数上的 Kronecker multivariate LMM 扫描 SNP，避免构建大规模 observation-level 协方差矩阵。`--tclass` 指定时间依赖分类协变量，`--tcovar` 指定时间依赖连续协变量；LOCO 模式支持 PLINK/PGEN/BGEN，但每条染色体会重新估计 multivariate null model，速度慢于 `--no-loco`。

RRM 基函数和 LOCO 选项：
- `--rrm-basis legendre`：默认 Legendre 多项式基函数，用于和 GMA/GMAT 对齐。
- `--rrm-basis bspline`：open-uniform B-spline 基函数；`--rrm-order n` 表示输出 `n+1` 个 basis 系数。
- `--rrm-loco-global-null`：快速近似 LOCO，只按染色体分组扫描 SNP，复用 global BLUP null model；默认不开启，默认 LOCO 仍逐染色体重估 null model。

非平衡模式输出：
- `*.rrm_unbal.assoc.txt`：非平衡 RRM 联合检验结果，包含 `JOINT_CHISQ`、`JOINT_P`、`BETA_L0...` 和对应标准误。
- `*.rrm_unbal.null.txt`：BLUP coefficient LMM 的 null 方差组分矩阵，记录遗传协方差 `K_*` 和残差/永久环境合并协方差 `PE_*`。
- `*.rrm_unbal.dynamic_curve.txt`：top SNP 在时间网格上的动态效应曲线。
- `*.rrm_unbal.dynamic_summary.txt`：top SNP 的关键时间、显著窗口和轨迹分类。

GMAT mouse_long 非平衡数据验证结果：1212 个体、11833 SNP、19392 条纵向观测、16 个时间点、`--rrm-order 3`。CosmicScanner BLUP-based RRM 在 `--no-loco` 下运行约 7.9 秒、峰值内存约 195 MB；与 GMAT `longwas --unbalance` 对齐达到 `Pearson r(-log10P)=0.951`、`Spearman rho(P)=0.940`，top SNP 同为 `UNC18848064`。相比早期 observation-level moment 路径约 1682 秒和 17.7 GB，速度提升约 213 倍、内存降低约 93 倍。

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
- `tests/run_gmat_unbalanced_alignment.sh`：使用 GMAT mouse_long 非平衡数据，比较 GMAT `longwas --unbalance` 与 CosmicScanner `--rrm-unbalance`；默认要求 `Pearson_log10P >= 0.94`、`Spearman_P >= 0.93` 且 top SNP 一致。
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
- [x] REGENIE 风格 Biobank 两步法 (`--biobank-mode`)
- [x] GEMMA-style eigen 工具（`--gemma-eigen`，输出 `.eigenD.txt/.eigenU.txt`）
- [x] GEMMA-style summary-stat 工具（`--gs/-gs` 输出 `.S.txt/.size.txt/.snps.txt`，含 `--cat/--mcat/--wcat/--wsnp` 多 component S；`--gq/-gq --beta` 支持 `Z`、`BETA+SE`、`CHISQ` 或 `P`，输出 `.Vq.txt/.q.txt/.size.txt`；`--study/--ref` 和 `--mstudy/--mref --vc/--ci` 输出 summary VC/CI；`--beta/--cor --vc` 支持 GEMMA `.cor.txt` 直接 summary VC）
- [x] GEMMA-style SNP correlation 工具（`--calccor --windowbp/--windowcm`，输出 `.cor.txt`）
- [x] GEMMA-style variance component 工具（`--vc 1|2|3`，输出 `.vc.txt/.vc.summary.txt`）
- [x] 输入适配层 MVP（csv/tsv/txt → prepared input，支持 repeatability/rrm/lmm/multivariate）

### GWAS 模型套件（含 GEMMA/EMmvGWAS 兼容路线）
- [x] GEMMA-style 普通线性模型 GWAS（`--lm --lm-mode 1|2|3|4`，对应 GEMMA `-lm` Wald/LRT/Score/all；输出 `.lm.assoc.txt`）
- [x] 单性状 LMM GWAS（默认 `--gwas`；`--lmm-gemma-mode 1|2|3|4` 对应 GEMMA `-lmm` Wald/LRT/Score/all 输出；`--lmm-gemma-refit` 启用逐 SNP alternative profile refit）
- [x] GEMMA-style prediction（`--predict-effect <assoc/effect.tsv> --predict-score <out>`，基于 SNP effect 表输出个体预测值）
- [x] GEMMA-style multi-locus 输出（`--bslmm 2` ridge/GBLUP；`--bslmm 1/3` deterministic sparse+ridge/probit approximation；`--bslmm 4/5` DAP BF/HYP 与 fine-mapping summary 近似；兼容 `-hmin/-hmax/-rmin/-rmax/-pmin/-pmax/-gmean/-hscale/-rscale/-pscale/-mh/-smin/-smax/-bf/-hyp`；输出 `.param.txt/.bv.txt/.hyp.txt` 或 `.bf.txt/.gamma.txt/.coef.txt`）
- [x] 多性状 GWAS（`--multivariate`，EMMV、EMmvGWAS-compatible `vc_estimation/full_gwas/semi_gwas` alias、GEMMA mode 1/2/3/4、GEMMA-full、GEMMA-exact-style、GEMMA-GxE、hybrid、EB）
- [x] 重复力模型 GWAS（`--repeatability-model`，exact/fast/STCG/refit-validation 四种策略）
- [x] 随机回归 RRM GWAS（`--long-pheno` + `--rrm-order`，Legendre/B-spline 基函数，平衡/非平衡纵向）
- [x] 母体效应模型（`--maternal-model`，ad/am/adm/c/e 多组分 VCE）
- [x] 单步 ssGWAS（`--single-step`，H⁻¹ 矩阵，系谱+基因组联合）
- [x] 阈值/Probit 模型（`--threshold-model --link probit|logit`）
- [x] 计数/二项 PSR-GWAS（`--count-model --family poisson|binomial`，PSR/PQL + weighted LMM）
- [x] 序数 PSR-GWAS（`--ordinal-model`，cumulative-binary + stacked-psr，dense-zkz + grid-reml VC 估计）
- [x] 生存 Cox 模型（`--survival-model`，Cox PH + 鞅残差 LMM）
- [x] GxE 反应范数模型（`--gxe-model`，Vg0+Vg1 双 VC，2-df 联合 Wald）

### 检验与扩展
- [x] 罕见变异检验（SKAT + Burden，beta 权重，混合卡方 p 值）
- [x] SPA 鞍点近似（二元性状，CGF 求解）
- [x] 经验贝叶斯 EB-GWAS（单性状/重复力/多性状稳定辅助路线）
- [x] PSR-GWAS v0.5 序数 stacked-psr 落地（Wald/pseudo-score 校准通过）

### 当前开发焦点
- PSR-GWAS v0.5 多分支整合收尾（大 N 校准、score test 诊断已完成）
- Stage 2C AI-REML 优化与大规模 N 支持规划中
- 重复力/RRM PSR-GWAS 与 exact GLMM 为后续路线图

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
