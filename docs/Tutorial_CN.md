# CosmicScanner 官方教程 (Tutorial)

**作者**: 黄杰 (Jie Huang)  
**版本**: 1.1.0  
**GitHub**: [https://github.com/huangjie2025/CosmicScanner](https://github.com/huangjie2025/CosmicScanner)

---

## 1. 软件简介
**CosmicScanner** 是一款专为全基因组关联分析（GWAS）打造的高性能计算工具，它是 **Cosmic** 全组学分析套件的核心组件之一。它打破了传统工具在计算速度、内存占用与模型支持上的瓶颈，能够轻松应对超大规模的群体基因型数据和复杂的表型数据。

**支持平台**:
- **Linux**: 原生支持各大 Linux 发行版 (x86_64 架构)。
- **Windows**: 推荐通过 WSL2 (Windows Subsystem for Linux) 运行。

---

## 2. 核心功能概览
CosmicScanner 目前支持以下核心功能：
1. **基础单性状 LMM GWAS**: 采用精确方差组分估计与留一染色体 (LOCO) 策略，彻底解决群体分层与近端污染。
2. **多性状联合 GWAS (Multivariate LMM)**: 直接对多个表型进行联合关联检验，提升检测效能 (Power)。
3. **重复力模型 GWAS (Repeatability GWAS)**: 支持同一个体存在多次测量记录的长格式数据，提供记录层精确 GLS 扫描与均值聚合近似扫描。
4. **随机回归模型 GWAS (RRM GWAS)**: 利用 Legendre 多项式对随时间或环境梯度的纵向表型进行联合轨迹关联检验。
5. **罕见变异分析**: 集成了 SKAT、Burden 检验，并支持针对二元性状的 SPA（鞍点近似）检验。
6. **PSR-GWAS 非正态性状分析**: 支持 binary、binomial event/trial 和 Poisson 计数性状的 PSR approximation + weighted LMM 扫描。
7. **超大规模队列分析**: 内置 REGENIE 风格的 Ridge Regression (Biobank Mode)，应对百万级样本。
8. **极速 GRM 构建**: 采用优化的 Bit-Parallel 算法，并支持零拷贝 `mmap` 缓存复用。

---

## 3. 安装与获取

### 3.1 预编译版本 (推荐)
您可以在 GitHub 的 Release 页面直接下载最新的预编译二进制文件，解压即用：
```bash
# 假设下载了最新版压缩包
tar -zxvf CosmicScanner-1.1.0-Linux-x86_64.tar.gz
chmod +x cosmicscanner
./cosmicscanner --help
```

### 3.2 源码编译
确保系统中安装了 `CMake (>= 3.15)` 和支持 C++17 的编译器 (如 GCC 9+)：
```bash
git clone --recurse-submodules https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j4
```
编译成功后，可执行文件位于 `build/cosmicscanner`。

---

## 4. 数据准备

CosmicScanner 需要以下标准格式的数据作为输入：
- **基因型文件**: 支持 PLINK 1.9 (`.bed/.bim/.fam`)、PLINK 2.0 (`.pgen/.pvar/.psam`) 以及 BGEN (`.bgen/.sample`) 格式。
- **表型文件 (`--pheno`)**: 包含三列或以上的文本文件 (FID, IID, Pheno1, Pheno2...)，缺失值用 `NA` 或 `-9` 表示。对于重复力或 RRM 模型，支持包含多行的长格式数据。
- **协变量文件 (`--covar` / `--qcovar`)**:
  - `--covar`: 离散型协变量（如性别、批次）。
  - `--qcovar`: 连续型协变量（如年龄、体重）。

### 4.1 输入适配层 MVP：原始表格转 prepared input

对于重复力 GWAS，CosmicScanner 现在提供官方输入准备器。首版支持 `csv/tsv/txt`，采用显式列声明，先把真实业务表格转换为可审计的中间格式，再交给现有重复力模型。

```bash
./cosmicscanner \
  --bfile data/genotype \
  --input-file data/semen.csv \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration,motility \
  --record-qcovars age,temp \
  --subject-ccovars farm,line \
  --time-parser date \
  --time-index-multiplier auto \
  --prepare-input-only \
  --prepared-input-dir prep/semen
```

准备目录会包含：

```text
sample_manifest.tsv
record_manifest.tsv
pheno_long.tsv
traits/<trait>.pheno.tsv
covar_record_q.tsv / covar_record_c.tsv
covar_subject_q.tsv / covar_subject_c.tsv
model_qcovar.tsv / model_covar.tsv
id_match_report.tsv
input_audit.json
summary.tsv
trait_manifest.tsv
covariate_audit.tsv
subject_covar_check.tsv
unresolved_columns.tsv
```

复用 prepared input 运行单个性状：

```bash
./cosmicscanner \
  --bfile data/genotype \
  --prepared-input-dir prep/semen \
  --input-trait volume \
  --repeatability-model \
  --gwas \
  --out results/volume_repeatability
```

也可以准备后直接运行，未指定 `--input-trait` 时默认使用 `--trait-cols` 的第一个性状：

```bash
./cosmicscanner \
  --bfile data/genotype \
  --input-file data/semen.csv \
  --input-model repeatability \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration,motility \
  --input-trait volume \
  --gwas \
  --out results/volume_repeatability
```

当前 MVP 不直接读取 Excel；请先另存为 `csv/tsv`。`--time-unique-policy within-time-index` 会为同一个体同一时间的多条记录生成唯一 `TIME`，并在 `record_manifest.tsv` 中保留原始时间和组内序号。

常用增强参数：

- `--input-config input_prepare.json`：从配置文件读取输入适配设置，适合正式项目复用；CLI 参数优先级高于配置文件。
- `--write-template-config input_prepare.json`：生成可编辑模板配置后退出，不需要提供基因型文件。
- `--time-index-multiplier auto|<n>`：控制 `TIME = TIME_BASE * multiplier + WITHIN_TIME_INDEX`，默认 `auto` 会按同个体同时间最大重复数自动选择 10、100、1000 等安全倍率。
- `--preview-schema`：只输出 `schema_preview.tsv` 并退出，用于检查列映射、ID 匹配数、重复 ID-time 键和 trait 非缺失数。
- `--infer-columns`：输出 `column_inference.tsv` 和 `inferred_input_config.json`，基于列名、数值比例、日期比例和 genotype ID 匹配率给出 ID/TIME/trait/covariate 候选建议；当前只做建议，不自动替代显式列声明。
- `--wide-to-long`：将 `volume_d1,day2_volume,concentration_1,concentration_T2` 这类宽表列转换为 `FID IID TIME TRAIT VALUE` 长表；此模式下 `--trait-cols` 填宽表列名，支持 `trait_1`、`trait_d1`、`day1_trait`、`trait_day1`、`trait_T1` 等形式。

模型接入范围：

- `--input-model repeatability`：使用 `traits/<trait>.pheno.tsv` 和 record-level covariate。
- `--input-model rrm`：使用 `rrm/<trait>.pheno.tsv`，轨迹时间为 `TIME_BASE`，record q/c covariate 自动接入 `tcovar/tclass`。
- `--input-model lmm`：使用 `subject/<trait>.pheno.tsv`，多记录按个体均值聚合。
- `--input-model multivariate`：使用 `subject/multivariate.pheno.tsv`，多性状按个体均值聚合。

Excel 当前策略：C++ 输入层会安全拒绝 `xlsx/xls`，请先另存为 `csv/tsv`。后续如引入 Excel，会作为可选依赖或官方转换器处理。

生成模板配置：

```bash
./cosmicscanner \
  --input-file data/semen.csv \
  --write-template-config input_prepare.json
```

配置文件示例：

```json
{
  "input_file": "data/semen.csv",
  "input_format": "csv",
  "input_model": "repeatability",
  "prepared_input_dir": "prep/semen",
  "id_col": "pig_id",
  "time_col": "collect_date",
  "trait_cols": ["volume", "concentration", "motility"],
  "record_qcovars": ["age", "temp"],
  "subject_ccovars": ["farm", "line"],
  "time_parser": "date",
  "time_unique_policy": "within-time-index",
  "time_index_multiplier": "auto"
}
```

使用配置文件准备输入：

```bash
./cosmicscanner \
  --bfile data/genotype \
  --input-config input_prepare.json \
  --prepare-input-only
```

只做列推断预览：

```bash
./cosmicscanner \
  --bfile data/genotype \
  --input-file data/semen.csv \
  --infer-columns \
  --preview-schema \
  --prepared-input-dir prep/preview
```

输入适配专项测试：

```bash
bash tests/run_test_input_adaptation.sh
```

该测试覆盖 long CSV、配置文件、列推断、增强审计、同日重复记录、strict 重复键失败、schema preview、wide-to-long；如果真实公猪精液 CSV 路径可访问，还会自动与旧 Python 脚本生成的真实数据 `summary.tsv` 做关键指标回归对齐。

正式接口规范见 `docs/Input_Adaptation_Spec_CN.md`。RRM 后续复用 prepared input 的设计见 `docs/RRM_Prepared_Input_Reuse_Design_CN.md`。

---

## 5. 详细功能与使用示例

### 5.1 标准单性状 LMM GWAS
利用留一染色体 (LOCO) 策略进行标准的 LMM GWAS 扫描：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --mpheno 1   --covar data/my_covar.txt   --out results/standard_gwas
```
**输出文件**: `*.assoc.txt`（包含 CHR, SNP, BP, A1, A2, AF, BETA, SE, P 等列）。

### 5.2 多性状联合 GWAS (Multivariate LMM)
同时分析多个表型，捕捉表型间的遗传协方差：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2,3   --multivariate   --mv-method emmv   --out results/multi_gwas
```
默认 `--mv-method emmv` 采用 EMmvGWAS 风格的半精准多性状 Wald 扫描：null model 通过 GRM 特征分解和 EM 估计性状协方差，扫描阶段复用 null 方差结构，并对每个 SNP 更新残差协方差后计算多自由度 Wald 检验。`--mv-method wald` 保留固定 null 方差结构的基础 Wald 路线，主要用于调试对照。`--mv-method gemma` 提供 GEMMA 风格逐 SNP 方差组分 refit Wald 路线，速度较慢，建议用于小样本精准验证或 top SNP 精修。`--mv-method gemma-full --mv-test wald,lrt,score` 输出 GEMMA-full MVP 的 Wald、LRT、Score 三类检验 p 值。

`--mv-method hybrid` 是当前推荐的性能优化路线：先用 EMmvGWAS 快扫全基因组，可选对 Stage 1 top SNP 做 Score-only prepass，再将 topK、p-threshold 和 Score/Wald discrepancy 候选 SNP 的并集送入 GEMMA-full 精修。不加 `--no-loco` 时，hybrid 默认启用 LOCO-aware screening 和 per-chrom LOCO GEMMA-full refinement。

```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2   --multivariate   --mv-method hybrid   --mv-test wald,lrt,score   --mv-refine-top 100   --mv-refine-p 1   --mv-refine-score-prepass 500   --mv-refine-discrepancy-top 50   --mv-refine-discrepancy 1   --gwas   --out results/mv_hybrid_loco
```

Hybrid 输出 `*.assoc.mv.txt`、`*.refine.snplist`、`*.refined.assoc.mv.txt` 和 `*.hybrid.summary.tsv`。若启用 `--mv-refine-score-prepass`，还会输出 `*.score_prepass.refined.assoc.mv.txt`。

如果需要和 EMmvGWAS R 包进行一致性验证，可使用：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2   --multivariate   --mv-method emmv   --mv-compatible emmvgwas   --no-loco   --out results/emmvgwas_compatible
```
`--mv-compatible emmvgwas` 会自动启用 EMmvGWAS-compatible GRM，即内部等价于 `--make-grm-alg 6`：使用 `Z = dosage - 1` 的 `-1/0/1` 编码，并按 `mean(diag(Z Z'))` 归一化。该模式主要用于和 EMmvGWAS R 包做方法学兼容验证，不建议作为默认 GRM 替代 VanRaden 或 GEMMA 风格 GRM。常规 GWAS 仍建议使用默认 `--make-grm-alg 1`。

`--mv-method gemma-full` 默认使用 GEMMA 风格的快速路径（null 方差组分），与 GEMMA mvLMM 行为一致。对于需要更高精度的场景，可用 `--mv-refit-p` 启用逐位点方差组分 refit：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/multi_pheno.txt   --mpheno-list 1,2,3   --multivariate   --mv-method gemma-full   --mv-test wald,lrt,score   --mv-refit-p 0.05   --no-loco   --out results/mv_gemma_full_refit
```
`--mv-refit-p` 默认 `1.0` 快速路径（null 组分）；设为小于 1 时，先计算 null Score 检验指标，仅对 Score p 值不超过阈值的 SNP 执行完整 Wald/LRT refit，其余 SNP 保留 Score 输出。大规模全基因组分析优先推荐 `--mv-method hybrid`。

### 5.3 重复力模型 GWAS (Repeatability Model)
针对同一受试对象有多条测量记录的数据（长格式表型文件）：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/repeat_pheno.txt   --repeatability-model   --repeatability-grm-cache grm_cache   --out results/repeat_gwas
```
**说明**: `--repeatability-grm-cache` 会将个体层 GRM 缓存至磁盘并通过 `mmap` 复用，大幅加速重复扫描。
生产分析推荐使用 `--repeatability-mode exact`，该模式保留记录层重复测量模型，并使用个体层 Woodbury 求解。`--repeatability-mode fast` 是个体均值聚合近似，只建议用于记录数平衡且不依赖记录层协变量的快速筛查。`--vce-mode auto` 默认使用可扩展的 subject-level Woodbury VCE；`--vce-mode full-record-exact` 仅建议小数据验证使用。

### 5.4 PSR-GWAS 非正态性状模块

PSR-GWAS 面向非正态性状的快速全基因组扫描。当前实现采用 PSR/PQL 线性化，生成 pseudo response 和 working residual variance，并通过 weighted LMM 扫描 SNP。它不是 exact GLMM；显著位点建议再用 candidate refinement 或外部 exact GLMM/GMMAT/POLMM 做复核。

二元性状可以使用 probit/logit threshold 模型：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/binary.pheno \
  --mpheno 1 \
  --gwas \
  --threshold-model \
  --link probit \
  --out results/psr_binary
```

Poisson 计数性状支持 log-offset/exposure。`--offset` 文件为每个样本一行，顺序与 genotype sample 顺序一致，数值应为 log-offset：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/count.pheno \
  --mpheno 1 \
  --gwas \
  --count-model \
  --family poisson \
  --offset data/count.log_offset \
  --out results/psr_poisson
```

Binomial 比例性状请使用 event/trial 两列，不要只输入比例。文件格式为 `FID IID EVENT TRIAL`：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --count-model \
  --family binomial \
  --binomial-pheno data/binomial_event_trial.pheno \
  --out results/psr_binomial
```

当前限制和建议：

- `--family negbin` 仍是实验入口，overdispersion 估计尚未完成系统验证，暂不建议作为正式结论报告。
- 官方 PSR-GLMM/R 站点当前在开发环境中访问超时，暂未接入官方源码和 demo 数据；状态记录见 `tests/PSR_GLMM_R_REFERENCE_STATUS.md`。
- 候选位点 refit MVP 位于 `tests/psr_candidate_refine.R`，默认使用 R `glm`，并提供可选 `firth`、`glmer` 和 `glmmTMB` hook。这些可选方法分别依赖 `logistf`、`lme4` 和 `glmmTMB`；缺包时不会导致验证脚本崩溃，而是逐 SNP 输出 `converged=FALSE` 和 `error_message`。

遇到难收敛或极端稀疏性状时，可以用 `--psr-maxit`、`--psr-tol`、`--psr-mu-eps`、`--psr-min-weight` 和 `--psr-max-weight` 调整 PSR/PQL 迭代与 clipping。默认值分别为 `50`、`1e-6`、`1e-6`、`1e-8` 和 `1e8`。

Convergence warning 是诊断信息，不等同于程序失败。`status_class=CONVERGED_WITH_WARNING` 表示迭代达到收敛标准，但存在边界方差、quasi-likelihood 非单调或 clipping 等诊断现象；`status_class=NOT_CONVERGED` 才表示正式分析中需要谨慎处理。`WARNING_PSR_MAXIT_REACHED` 表示达到迭代上限，`WARNING_PSR_CLIPPED` 表示 `mu` 或 working weight 触及 clipping 边界，`WARNING_PSR_SIGMA_BOUNDARY` 表示最终遗传方差接近 0，`WARNING_PSR_QUASILIK_DECREASED` 表示 quasi-likelihood proxy 出现下降。

PSR-GWAS smoke、R baseline 和 candidate refinement 验证：

```bash
bash tests/run_psr_gwas_validation.sh
```

该脚本会额外生成：

- `*.psr.convergence.tsv`：PSR/PQL 迭代诊断，包括 `mu` 范围、weights 范围、clipping 次数和收敛标记。
- `psr_convergence_summary.tsv`：每个 trait 一行的 convergence 汇总和 warning flags。
- `weighted_lmm_equivalence.tsv`：验证 `V = Kσg² + R` 与加权变换后的 weighted LMM 在小样本固定方差组分下等价。
- `psr_type1_summary.tsv`、`psr_power_summary.tsv`、`psr_qq_points.tsv`：由真实 CosmicScanner PSR-GWAS 运行生成的小规模 simulation smoke 摘要；正式论文验证仍需更大规模仿真和外部 GLMM/PSR-GLMM/R 基准。
- `method_metrics.tsv`：simulation 每次运行的机器可读指标，包括 `engine=cosmicscanner_psr`、seed、trait type、scenario、runtime、convergence flag 和 `status_class`。
- `psr_validation_report.md`：Markdown 验证报告，合并 convergence、weighted equivalence、simulation、candidate refinement 和已知限制。

PSR-GWAS 路线图和验收标准见 `docs/PSR_GWAS_ROADMAP.md`；v0.2 验证记录见 `docs/PSR_GWAS_V0.2_VALIDATION_SUMMARY.md`；stable sprint 总结见 `docs/PSR_GWAS_V0.3_STABLE_SUMMARY.md`。

将 PSR-GWAS 加入核心回归：

```bash
bash tests/run_gwas_core_regression.sh --fast --with-psr
```

### 5.5 随机回归模型 GWAS (RRM GWAS)
用于检测随时间或环境梯度变化的动态遗传效应（纵向数据）：
```bash
./cosmicscanner   --bfile data/my_geno   --long-pheno data/long_pheno.txt   --mpheno 2   --tpheno 1   --rrm-order 2   --gwas   --out results/rrm_gwas
```
**输出文件**:
- `*.rrm_joint.assoc.txt`: 观测层多自由度联合检验 P 值。
- `*.rrm.top_curve.txt`: Top 位点在时间轴上的重构轨迹效应。

### 5.6 罕见变异分析 (SKAT / Burden)
基于基因或自定义集合进行罕见变异关联检验：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --set-file data/gene_sets.txt   --skat   --out results/skat_gwas
```

### 5.7 Biobank 级别群体分析
针对数十万大样本，开启 REGENIE 模式的 Ridge Regression：
```bash
./cosmicscanner   --bfile data/my_geno   --pheno data/my_pheno.txt   --biobank-mode   --block-size 1000   --out results/biobank_gwas
```

---

## 6. 完整参数说明

| 参数类别 | 命令 | 说明 |
| --- | --- | --- |
| **输入/输出** | `--bfile <prefix>` | PLINK 1.9 二进制前缀 |
| | `--pfile <prefix>` | PLINK 2.0 二进制前缀 |
| | `--bgen <file>` | BGEN 格式基因型文件 |
| | `--out <prefix>` | 输出文件前缀，默认 `cosmic_gwas` |
| **表型/协变量** | `--pheno <file>` | 表型文件路径 |
| | `--mpheno <n>` | 单性状分析时指定的表型列号 |
| | `--mpheno-list <list>` | 多性状分析的表型列号 (如 1,2,3) |
| | `--covar <file>` | 离散型协变量文件 |
| | `--qcovar <file>` | 连续型协变量文件 |
| **模型/算法** | `--multivariate` | 开启多性状联合 GWAS |
| | `--mv-method <name>` | 多性状扫描方法：`emmv`、`wald`、`gemma`、`gemma-full` 或 `hybrid`，默认 `emmv` |
| | `--mv-compatible <name>` | 多性状兼容模式；`emmvgwas` 会启用 EMmvGWAS-compatible GRM，主要用于和 EMmvGWAS R 对标 |
| | `--mv-test <list>` | GEMMA-full 检验列表：`wald,lrt,score`，默认全开 |
| | `--mv-refine-top <n>` | hybrid/GEMMA 精修 SNP 数量上限，默认 100 |
| | `--mv-refine-p <p>` | hybrid/GEMMA 精修 p 值阈值，默认 1e-3 |
| | `--mv-refit-p <p>` | GEMMA-full 快速模式的完整 refit 阈值，默认 1 表示全部 SNP refit |
| | `--threshold-model` | 二元性状 PSR/probit 或 PSR/logit GWAS |
| | `--link probit\|logit` | threshold 模型 link，默认 `probit` |
| | `--count-model` | 计数/binomial PSR-GWAS 入口 |
| | `--family poisson\|binomial\|negbin` | count 模型分布；`negbin` 仍为实验入口 |
| | `--offset <file>` | Poisson log-offset/exposure，每个样本一行 |
| | `--binomial-pheno <file>` | Binomial event/trial 表型文件，格式 `FID IID EVENT TRIAL` |
| | `--psr-maxit <n>` | PSR/PQL 最大迭代次数，默认 50 |
| | `--psr-tol <f>` | PSR/PQL 收敛阈值，默认 1e-6 |
| | `--psr-mu-eps <f>` | PSR mean clipping epsilon，默认 1e-6 |
| | `--psr-min-weight <f>` | PSR working weight 下界，默认 1e-8 |
| | `--psr-max-weight <f>` | PSR working weight 上界，默认 1e8 |
| | `--repeatability-model` | 开启重复测量/重复力模型 |
| | `--long-pheno <file>` | 开启随机回归模型，并指定长格式表型文件 |
| | `--rrm-order <n>` | 随机回归模型的 Legendre 多项式阶数，需要配合 `--long-pheno` |
| | `--fast-lmm` | 开启近似 LMM 加速扫描 |
| | `--spa` | 开启鞍点近似检验 (Saddle Point Approximation) |
| | `--skat` / `--burden` | 开启罕见变异聚合检验 |
| | `--biobank-mode` | 开启大规模 Ridge Regression 预处理模式 |
| **工程优化** | `--make-grm` | 仅计算 GRM 矩阵并退出 |
| | `--make-grm-alg 6` | EMmvGWAS-compatible GRM：`Z=dosage-1` 且按 `mean(diag(Z Z'))` 归一化，仅建议用于兼容验证 |
| | `--repeatability-mixed-precision`| 重复力模型中使用 Float32 降低内存占用 |
| | `--repeatability-grm-cache <pre>`| 开启 mmap 磁盘缓存复用机制 |
# PSR-GWAS v0.4 对标说明

CosmicScanner PSR-GWAS 是 PSR/PQL 近似方法，不是 exact GLMM 推断。

- normal 性状：作为 Gaussian identity-link 特例，直接复用现有 LMM GWAS。
- binary 性状：`--threshold-model --link probit|logit`。
- Poisson 性状：`--count-model --family poisson`。
- binomial event/trial 性状：`--count-model --family binomial --binomial-pheno FID_IID_EVENT_TRIAL.txt`。
- ordinal 性状：`--ordinal-model --ordinal-pheno FID_IID_ORDINAL.txt`；当前 v0.4 是 `experimental_cumulative_binary`，不是 full stacked proportional-odds PSR-GLMM parity。

完整 v0.4 验证：

```bash
bash tests/run_psr_full_parity_validation.sh
```

官方 PSR-GLMM/R 只做 availability tracking，不作为 CI 阻塞依赖。替代 benchmark 使用 R GLM baseline，并在安装时记录 optional Firth、glmer、glmmTMB、ordinal 工具结果。
