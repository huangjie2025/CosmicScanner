# CosmicScanner 官方使用说明书 (Tutorial)

**作者**: 黄杰 (Jie Huang)
**开发树文档版本**: 1.1.0
**GitHub**: [https://github.com/huangjie2025/CosmicScanner](https://github.com/huangjie2025/CosmicScanner)
**更新日期**: 2026-08-24

> 本教程描述开发树中的扩展功能。论文纳入范围与软件成熟度是两个维度；
> 发布候选的软件成熟度仅使用 `Stable`、`Experimental validated` 和
> `Experimental approximate`，规范分类见
> [论文模型支持矩阵](PAPER_MODEL_SUPPORT_MATRIX.md)。

---

## 目录

- [1. 软件简介](#1-软件简介)
- [2. 核心功能概览](#2-核心功能概览)
- [3. 安装与获取](#3-安装与获取)
- [4. 数据准备](#4-数据准备)
- [5. 快速开始](#5-快速开始)
- [6. 详细功能与使用示例](#6-详细功能与使用示例)
  - [6.1 单性状 LMM GWAS](#61-单性状-lmm-gwas)
  - [6.2 多性状联合 GWAS (Multivariate LMM)](#62-多性状联合-gwas-multivariate-lmm)
  - [6.3 重复力模型 GWAS (Repeatability Model)](#63-重复力模型-gwas-repeatability-model)
  - [6.4 随机回归模型 GWAS (RRM GWAS)](#64-随机回归模型-gwas-rrm-gwas)
  - [6.5 母体效应模型 GWAS](#65-母体效应模型-gwas)
  - [6.6 单步 ssGWAS](#66-单步-ssgwas)
  - [6.7 阈值/Probit 模型 GWAS](#67-阈值probit-模型-gwas)
  - [6.8 计数/二项 PSR-GWAS](#68-计数二项-psr-gwas)
  - [6.9 序数 PSR-GWAS](#69-序数-psr-gwas)
  - [6.10 生存 Cox 模型 GWAS](#610-生存-cox-模型-gwas)
  - [6.11 GxE 反应范数模型 GWAS](#611-gxe-反应范数模型-gwas)
  - [6.12 罕见变异分析 (SKAT / Burden)](#612-罕见变异分析-skat--burden)
  - [6.13 Biobank 级别群体分析](#613-biobank-级别群体分析)
  - [6.14 经验贝叶斯 EB-GWAS](#614-经验贝叶斯-eb-gwas)
  - [6.15 系谱动物模型 GWAS](#615-系谱动物模型-gwas)
- [7. 输出文件格式说明](#7-输出文件格式说明)
- [8. 完整参数说明](#8-完整参数说明)
- [9. 性能调优指南](#9-性能调优指南)
- [10. 常见问题 (FAQ)](#10-常见问题-faq)
- [11. 相关文档](#11-相关文档)

---

## 1. 软件简介

**CosmicScanner** 是一款专为全基因组关联分析（GWAS）打造的高性能计算工具，它是 **Cosmic** 全组学分析套件的核心组件之一。它打破了传统工具在计算速度、内存占用与模型支持上的瓶颈，能够轻松应对超大规模的群体基因型数据和复杂的表型数据。

**支持平台**:
- **Linux**: 原生支持各大 Linux 发行版 (x86_64 架构)。
- **Windows**: 推荐通过 WSL2 (Windows Subsystem for Linux) 运行。

---

## 2. 核心功能概览

CosmicScanner 当前支持以下核心功能：

### GWAS 模型套件（核心模型与兼容路线）

| 编号 | 模型 | 适用场景 | 关键参数 |
|------|------|---------|---------|
| 1 | 单性状 LMM GWAS | 标准单表型关联分析；GEMMA `-lmm` mode/refit 输出 | `--gwas` / `--lmm-gemma-mode` / `--lmm-gemma-refit` |
| 2 | GEMMA-style LM GWAS | 无混合效应的普通线性模型；GEMMA `-lm` 对标 | `--lm --lm-mode 1\|2\|3\|4` |
| 3 | 多性状联合 GWAS | 多表型联合检验，提升 Power | `--multivariate --mv-method` |
| 4 | 重复力模型 GWAS | 同个体多次测量记录 | `--repeatability-model` |
| 5 | 随机回归 RRM GWAS | 纵向/动态性状轨迹 | `--long-pheno --rrm-order` |
| 6 | 母体效应模型 GWAS | 出生重/断奶重等母体性状 | `--maternal-model` |
| 7 | 单步 ssGWAS | 系谱+基因组联合分析 | `--single-step --ped` |
| 8 | 阈值/Probit 模型 GWAS | 二分类性状 | `--threshold-model --link` |
| 9 | 计数/二项 PSR-GWAS | Poisson/Binomial 性状 | `--count-model --family` |
| 10 | 序数 PSR-GWAS | 3+ 有序分类性状 | `--ordinal-model` |
| 11 | 生存 Cox 模型 GWAS | 生存/时间到事件性状 | `--survival-model` |
| 12 | GxE 反应范数模型 GWAS | 环境梯度交互 | `--gxe-model --env-covar` |
| 13 | GEMMA-style 工具路线 | VC、summary-stat VC/CI、GRM eigen、SNP correlation、effect prediction、multi-locus 输出 | `--vc` / `--ci` / `--gs` / `--gq` / `--beta` / `--cor` / `--cat` / `--mcat` / `--wcat` / `--wsnp` / `--study` / `--ref` / `--mstudy` / `--mref` / `--gemma-eigen` / `--calccor` / `--predict-effect` / `--bslmm` |

### 论文模型的发布候选状态

| 模型/路线 | 软件状态 |
|---|---|
| General LMM；多性状 LMM | Stable |
| 重复力；RRM Method A；动物模型；母体模型；共同环境模型 | Experimental validated |
| H-matrix single-step GWAS（Method 1） | Experimental validated |
| ssGBLUP back-solving GWAS（Method 2） | Experimental validated |
| Exact back-solved SNP P-value | Stable |
| APY-exact back-solved SNP P-value | Experimental validated |
| APY-approx back-solved SNP P-value | Experimental approximate |

Method 1 的固定 SNP 回归系数与 Method 2 的反推标记效应是不同估计量，
论文和结果解释中不得仅用未限定的 `ssGWAS` 混称。APY-approx 原始 P 值
不保证严格 I 类错误控制；CosmicScanner 不对该路线应用 genomic control
或 P-value rescaling。

### 引擎与基础设施

| 编号 | 功能 | 关键参数 |
|------|------|---------|
| 14 | 罕见变异分析 (SKAT/Burden/SPA) | `--skat` / `--burden` / `--spa` |
| 14 | Biobank 两步法 (REGENIE 风格) | `--biobank-mode` |
| 15 | IDUL/AIDUL 加速求解器 | `--lmm-null-solver` |
| 16 | 输入适配层 MVP | `--input-file --input-model` |
| 17 | 极速 GRM 构建 | `--make-grm --make-grm-alg` |
| 18 | 经验贝叶斯 EB-GWAS（稳定辅助路线） | `--model lmm-eb` |

---

## 3. 安装与获取

### 3.1 预编译版本（推荐）

正式发布后可从对应 Release 页面下载预编译二进制文件。当前本地发布候选包名为：

```bash
tar -zxvf CosmicGWAS-v1.1.0-linux-x86_64.tar.gz
cd CosmicGWAS-v1.1.0-linux-x86_64
bash install.sh
bin/cosmicscanner --help
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

### 3.3 验证安装

```bash
./cosmicscanner --help
```

若输出完整的参数列表，说明安装成功。

---

## 4. 数据准备

CosmicScanner 需要以下标准格式的数据作为输入。

### 4.1 基因型文件

支持三种主流格式：

| 格式 | 参数 | 文件组成 |
|------|------|---------|
| PLINK 1.9 | `--bfile <prefix>` | `.bed` / `.bim` / `.fam` |
| PLINK 2.0 | `--pfile <prefix>` 或 `--pgen`/`--pvar`/`--psam` | `.pgen` / `.pvar` / `.psam` |
| BGEN | `--bgen <file> --sample <file>` | `.bgen` / `.sample` (Oxford 格式) |

**位点过滤**：
- `--maf <float>`：最小等位基因频率过滤（默认 0.01）
- `--hwe <float>`：HWE 精确检验 p 值过滤（默认 0.0，即关闭）
- `--extract <file>`：仅包含指定 SNP ID 列表
- `--exclude <file>`：排除指定 SNP ID 列表
- `--chr <chr>`：限制扫描单条染色体
- `--from-bp <pos>` / `--to-bp <pos>`：限制扫描区间

### 4.2 表型文件

**标准格式**（单性状 / 多性状）：

```
FID    IID    Pheno1    Pheno2    Pheno3
IND001 IND001 1.23      2.45      NA
IND002 IND002 0.88      NA        1.67
```

- 缺失值用 `NA` 或 `-9` 表示。
- `--mpheno <n>`：指定表型列号（1-indexed），默认 1。
- `--mpheno-list <list>`：多性状分析指定列号（如 `1,2,3`）。

**长格式**（重复力 / RRM 模型）：

```
FID    IID    TIME    VALUE
IND001 IND001 1       1.23
IND001 IND001 2       1.35
IND002 IND002 1       0.88
```

- `--tpheno <n>`：指定 TIME 列号。

**生存性状格式**（Cox 模型）：

```
FID    IID    TIME    STATUS
IND001 IND001 365     1
IND002 IND002 720     0
```

- `STATUS`：1 = 事件发生，0 = 删失。

**母体效应格式**：

```
FID    IID    DAM_ID    PHENO
IND001 IND001 DAM001    1.23
```

- 第 3 列为母体 ID（可用 `--maternal-id <col>` 指定列号）。

**Binomial event/trial 格式**：

```
FID    IID    EVENT    TRIAL
IND001 IND001 45       100
```

- 通过 `--binomial-pheno <file>` 指定，不要只输入比例值。

### 4.3 协变量文件

| 参数 | 类型 | 格式 |
|------|------|------|
| `--covar <file>` | 离散型（性别、批次） | `FID IID CAT1 CAT2 ...` |
| `--qcovar <file>` | 连续型（年龄、体重） | `FID IID VAL1 VAL2 ...` |

### 4.4 输入适配层 MVP：原始表格转 prepared input

对于重复力 GWAS、RRM 等长格式模型，CosmicScanner 提供官方输入准备器，将 `csv/tsv/txt` 业务表格转换为可审计的中间格式。

**基本用法**：

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

**准备目录输出**：

```
sample_manifest.tsv          # 样本清单
record_manifest.tsv          # 记录清单
pheno_long.tsv               # 长格式表型
traits/<trait>.pheno.tsv     # 按性状拆分的表型
covar_record_q.tsv           # 记录层连续协变量
covar_record_c.tsv           # 记录层离散协变量
covar_subject_q.tsv          # 个体层连续协变量
covar_subject_c.tsv          # 个体层离散协变量
model_qcovar.tsv             # 模型连续协变量
model_covar.tsv              # 模型离散协变量
id_match_report.tsv          # ID 匹配报告
input_audit.json             # 输入审计日志
summary.tsv                  # 汇总
trait_manifest.tsv           # 性状清单
covariate_audit.tsv          # 协变量审计
subject_covar_check.tsv      # 个体协变量检查
unresolved_columns.tsv       # 未解析列
```

**复用 prepared input 运行单个性状**：

```bash
./cosmicscanner \
  --bfile data/genotype \
  --prepared-input-dir prep/semen \
  --input-trait volume \
  --repeatability-model \
  --gwas \
  --out results/volume_repeatability
```

**使用配置文件**（适合正式项目复用）：

```bash
# 生成模板配置
./cosmicscanner --input-file data/semen.csv --write-template-config input_prepare.json

# 使用配置文件准备输入
./cosmicscanner --bfile data/genotype --input-config input_prepare.json --prepare-input-only
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

**模型接入范围**：

| `--input-model` | 使用文件 | 说明 |
|----------------|---------|------|
| `repeatability` | `traits/<trait>.pheno.tsv` + record covariate | 记录层重复力模型 |
| `rrm` | `rrm/<trait>.pheno.tsv` | 轨迹时间为 `TIME_BASE` |
| `lmm` | `subject/<trait>.pheno.tsv` | 多记录按个体均值聚合 |
| `multivariate` | `subject/multivariate.pheno.tsv` | 多性状按个体均值聚合 |

**增强参数**：

| 参数 | 说明 |
|------|------|
| `--preview-schema` | 只输出 `schema_preview.tsv` 并退出 |
| `--infer-columns` | 输出列推断建议（不自动替代显式声明） |
| `--wide-to-long` | 将宽表列（`trait_d1`, `day1_trait`）转换为长表 |
| `--time-unique-policy strict\|within-time-index` | 同个体同时间重复记录处理策略 |
| `--time-index-multiplier auto\|<n>` | TIME 编码倍率，`auto` 自动选择安全倍率 |

> **注意**：当前 MVP 不直接读取 Excel (`xlsx/xls`)，请先另存为 `csv/tsv`。

正式接口规范见 [docs/INPUT_ADAPTATION_SPEC_CN.md](INPUT_ADAPTATION_SPEC_CN.md)。

---

## 5. 快速开始

以下是最基础的单性状 LMM GWAS 示例：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --covar data/my_covar.txt \
  --gwas \
  --out results/standard_gwas
```

输出文件 `results/standard_gwas.assoc.txt` 包含每个 SNP 的关联检验结果。

---

## 6. 详细功能与使用示例

### 6.1 单性状 LMM GWAS

利用留一染色体 (LOCO) 策略进行标准的 LMM GWAS 扫描：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --covar data/my_covar.txt \
  --qcovar data/my_qcovar.txt \
  --gwas \
  --scan-tier semi \
  --out results/standard_gwas
```

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--no-loco` | 禁用 LOCO，使用全局 null 模型 | 启用 LOCO |
| `--scan-tier <m>` | `approx`/`semi`/`exact`，推荐全基因组使用 `semi` | 显式指定 |
| `--fast-lmm` | 单性状 fixed-null 旧兼容入口；新命令使用 `--scan-tier approx` | 关闭 |
| `--lmm-null-solver <s>` | null 模型求解器：`brent`/`idul`/`idul_plus`/`aidul` | `brent` |
| `--maf <f>` | 最小 MAF 过滤 | 0.01 |
| `--hwe <f>` | HWE p 值过滤 | 0.0（关闭） |
| `--threads <n>` | OpenMP 线程数 | 最大可用 |

**输出文件**：`*.assoc.txt`

### 6.2 多性状联合 GWAS (Multivariate LMM)

同时分析多个表型，捕捉表型间的遗传协方差：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2,3 \
  --multivariate \
  --gwas \
  --scan-tier semi \
  --out results/multi_gwas
```

普通多性状 LMM 使用统一档位：`approx` 为 fixed-null Wald，`semi` 为固定协方差比例并逐 SNP 更新残差协方差，`exact` 为逐 SNP full covariance refit。`exact` 当前要求 `--no-loco`。不要在同一命令中混用 `--scan-tier` 和 `--mv-method`。

**高级方法选择**：

| 方法 | 说明 | 适用场景 |
|------|------|---------|
| `emmv`（默认） | EMmvGWAS 风格半精准多性状 Wald | 快速全基因组扫描 |
| `wald` | 固定 null 方差结构基础 Wald | 调试对照 |
| `gemma` | GEMMA 风格逐 SNP refit Wald | 小样本精准验证 |
| `gemma-full` | GEMMA-style Wald/LRT/Score，支持 `--mv-gemma-mode 1/2/3/4` | GEMMA mvLMM 常规对标 |
| `gemma-exact` | 使用 GEMMA EM/NR 控制参数的 full alternative refit 路线 | 更接近 GEMMA exact 流程的验证 |
| `gemma-gxe` | 多性状 SNP×环境交互 Wald | GEMMA 多性状 GxE 快速扫描 |
| `gemma-gxe-full` | 多性状 GxE Wald/LRT/Score | GEMMA 多性状 GxE 多检验 |
| `hybrid`（推荐） | EMmv 快扫 + GEMMA-full 精修 | 大规模全基因组 |
| `eb` | 经验贝叶斯随机 SNP 效应 | 稳定辅助信号增强 |

`gemma-exact` 是 GEMMA exact-style 路线：支持 GEMMA 的 EM/NR 迭代控制、Wald/LRT/Score mode 选择和逐 SNP full alternative refit。加 `--mv-crt` 后启用 GEMMA-style Edgeworth/PCRT 小样本校正，输出 `P_WALD_CRT`、`P_LRT_CRT`、`P_SCORE_CRT`，主 `P` 列按 `--mv-gemma-mode` 使用校正后的 p 值。该实现为 CosmicScanner/Eigen 原生矩阵实现，没有直接链接 GEMMA 的 GSL 代码。

**hybrid 方法（推荐）**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2 \
  --multivariate \
  --mv-method hybrid \
  --mv-test wald,lrt,score \
  --mv-refine-top 100 \
  --mv-refine-p 1e-3 \
  --mv-refine-score-prepass 500 \
  --mv-refine-discrepancy-top 50 \
  --mv-refine-discrepancy 1.0 \
  --gwas \
  --out results/mv_hybrid_loco
```

**hybrid 输出文件**：
- `*.assoc.mv.txt`：Stage 1 全基因组扫描结果
- `*.refine.snplist`：精修 SNP 列表
- `*.refined.assoc.mv.txt`：精修结果
- `*.hybrid.summary.tsv`：hybrid 流程汇总
- `*.score_prepass.refined.assoc.mv.txt`：Score prepass 精修结果（如启用）

**GEMMA-full refit 模式**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2,3 \
  --multivariate \
  --mv-method gemma-full \
  --mv-test wald,lrt,score \
  --mv-refit-p 0.05 \
  --no-loco \
  --gwas \
  --out results/mv_gemma_full_refit
```

`--mv-refit-p` 默认 `1.0`（快速路径，null 组分）；设为小于 1 时，仅对 Score p ≤ 阈值的 SNP 执行完整 refit。

**GEMMA 多性状 GxE 路线**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_gxe_pheno.txt \
  --mpheno-list 2,3 \
  --multivariate \
  --mv-method gemma-gxe-full \
  --env-covar 1 \
  --mv-gemma-mode 4 \
  --no-loco \
  --gwas \
  --out results/mv_gemma_gxe
```

输出为 `*.assoc.mv.gxe.txt`，包含 `P_WALD_GXE`、`P_LRT_GXE`、`P_SCORE_GXE` 和各性状 `BETA_GXE_*`/`SE_GXE_*`。该路线支持 PLINK/PGEN/BGEN 输入；当前要求 `--no-loco`。

**EMmvGWAS R 包兼容验证**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/multi_pheno.txt \
  --mpheno-list 1,2 \
  --emmvgwas \
  --emmvgwas-output \
  --out results/emmvgwas_compatible
```

`--emmvgwas` 是 EMmvGWAS-compatible 路线的一键别名，等价于 `--multivariate --mv-method emmv --mv-compatible emmvgwas --no-loco`。它启用 EMmvGWAS-compatible GRM（等价 `--make-grm-alg 6`：`Z=dosage-1` 编码，按 `mean(diag(Z Z'))` 归一化）。`--emmvgwas-output` 写出 EMmvGWAS 风格 `G.csv`、`R.csv`、`beta.csv` 和 association CSV。

如需对接 EMmvGWAS R 包的分步流程，`--emmvgwas-vc-only` 对应 `vc_estimation`，只输出 G/R/beta 后停止；随后用 `--emmvgwas --emmvgwas-G <G.csv> --emmvgwas-R <R.csv>` 读入外部方差组分，执行 `semi_gwas` 等价扫描。直接 `--emmvgwas --emmvgwas-output` 对应 `full_gwas` 风格的一步扫描和输出。

**VCE 求解器控制**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--mv-vce-solver <s>` | `em` 或 `aireml` | `em` |
| `--mv-em-min-iter <n>` | EM 最小迭代次数 | 10 |
| `--mv-ai-max-iter <n>` | AI-REML 最大迭代次数 | 50 |
| `--mv-ai-tol <f>` | AI-REML 收敛阈值 | 1e-6 |

### 6.3 重复力模型 GWAS (Repeatability Model)

针对同一受试对象有多条测量记录的数据（长格式表型文件）：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat_pheno.txt \
  --tpheno 1 --mpheno 2 \
  --repeatability-model \
  --repeatability-mode exact \
  --gwas \
  --scan-tier semi \
  --repeatability-grm-cache grm_cache \
  --out results/repeat_gwas
```

> **命令规范**：当前多性状和重复力模型标志会兼容性地隐式启动扫描，但所有正式命令仍统一显式写 `--gwas`、表型列和 `--scan-tier`，避免依赖旧默认行为。

**策略选择**（`--repeatability-mode <m>`）：

| 策略 | 说明 | 适用场景 |
|------|------|---------|
| `auto` | 自动选择（记录数 > 80000 时自动切 STCG） | 默认 |
| `exact` | 记录层精确 GLS + 个体层 Woodbury | **生产分析推荐** |
| `fast` | 个体均值聚合近似 | 平衡数据快速筛查 |
| `stcg` | STCG 矩阵自由求解 | Biobank 级纵向数据 |

> **注意**：`refit-validation` 不是 `--repeatability-mode` 的合法取值，逐 SNP VCE refit 请使用独立参数 `--repeatability-refit-validation`。

**VCE 模式**（`--vce-mode <m>`）：

| 模式 | 说明 |
|------|------|
| `auto`（默认） | 可扩展的 subject-level Woodbury VCE |
| `fast` | 快速 VCE，跳过精确 REML 精修 |
| `subject-level-exact` | 强制 subject-level 精确 REML 精修 |
| `full-record-exact` | 完整记录层精确 VCE（仅小数据验证） |

**关键参数**：

| 参数 | 说明 |
|------|------|
| `--repeatability-grm-cache <prefix>` | 将个体层 GRM 缓存至磁盘，通过 `mmap` 复用 |
| `--repeatability-mixed-precision` | 使用 Float32 降低内存占用 |
| `--repeatability-weight <file>` | 记录权重：`FID IID WEIGHT` 或 `FID IID TIME WEIGHT` |
| `--repeatability-refit-validation` | 启用逐 SNP VCE refit |
| `--repeatability-refit-limit <n>` | refit SNP 数量上限（0 = 全部） |

**输出文件**：
- `*.assoc.txt`：主关联结果
- `*.refit.assoc.txt`：refit-validation 结果（如启用）
- `*.stcg.assoc.txt`：STCG 策略结果（如启用）

### 6.4 随机回归模型 GWAS (RRM GWAS)

用于检测随时间或环境梯度变化的动态遗传效应（纵向数据）：

**平衡纵向数据**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --long-pheno data/long_pheno.txt \
  --mpheno 2 \
  --tpheno 1 \
  --rrm-order 2 \
  --gwas \
  --out results/rrm_gwas
```

**非平衡纵向数据**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/unbal_pheno.txt \
  --rrm-unbalance \
  --rrm-gwas-method record \   # 或 coefficient（默认）
  --tpoint time_point \
  --tcovar age,bw \
  --tclass farm,parity \
  --rrm-order 2 \
  --rrm-basis legendre \
  --out results/rrm_unbal_gwas
```

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--rrm-order <n>` | 基函数阶数 (1-4) | — |
| `--rrm-basis <m>` | 基函数：`legendre` 或 `bspline` | `legendre` |
| `--rrm-gwas-method <m>` | `record`：记录级 `x⊗Φ` 联合检验（Type B，`Experimental validated`）；`coefficient`：BLUP 系数空间多性状 Wald（与 GMAT 对齐，默认） | `coefficient` |
| `--rrm-curve-grid <n>` | 动态曲线时间网格点数 | 21 |
| `--rrm-top <n>` | 输出 top SNP 曲线数 | 100 |
| `--rrm-ai-iter <n>` | null AI 精修迭代数 (-1=auto, 0=fast) | -1 |
| `--rrm-sig-threshold <p>` | 动态时间窗口 p 阈值 | 0.05 |
| `--rrm-loco-global-null` | 使用全局 null 的快速近似 LOCO | 关闭 |
| `--tpoint <col>` | 非平衡模式时间点列名 | — |
| `--tcovar <cols>` | 时间依赖连续协变量 | — |
| `--tclass <cols>` | 时间依赖离散协变量 | — |
| `--rrm-compatible <m>` | 兼容模式：`gma-fixreg` 或 `gmat-unbalance` | — |
| `--rrm-compatible-varcom <file>` | GMA `.var` 文件 | — |

**输出文件**：
- `*.rrm_joint.assoc.txt`：观测层多自由度联合检验 P 值
- `*.rrm.top_curve.txt`：Top 位点在时间轴上的重构轨迹效应
- `*.rrm_unbal.assoc.txt`：非平衡模式结果

非平衡 RRM 的两条扫描路径统计语义不同：

- `--rrm-gwas-method record`：记录级联合检验，SNP 固定设计为 `x_{subject(r)} × φ_k(t_r)`（genotype × 基函数），是真正的 time-varying SNP 效应检验（Type B），带 `STATUS/DF` 显式诊断（`OK/RANK_DEFICIENT/ILL_CONDITIONED/NUMERICAL_FAILURE`）。软件状态为 `Experimental validated`，证据包括平衡 joint 路径等价（Pearson ≥ 0.99999）与 A–G 模拟场景。
- `--rrm-gwas-method coefficient`（默认）：先对个体 Legendre 系数做 BLUP，再在系数空间做 Kronecker 多性状 Wald 扫描；`--rrm-compatible gmat-unbalance` 复现 GMAT `longwas --unbalance` 口径（`run_gmat_unbalanced_alignment.sh`，全集 Pearson=0.9968、Spearman=0.9959、top10 全一致），为历史兼容路径。

完整五层审计、对齐数值表与旧版 X 染色体差异复盘见 `docs/RRM_GWAS_METHOD_AUDIT.md`。

### 6.5 母体效应模型 GWAS

建模直接加性、母体加性、直接-母体协方差、母体永久环境等多组分方差结构，适用于出生重/断奶重等母体效应性状：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/birth_weight.txt \
  --maternal-model \
  --maternal-id 3 \
  --covar data/covar.txt \
  --gwas \
  --out results/maternal_gwas
```

**表型文件格式**：

```
FID    IID    DAM_ID    PHENO
IND001 IND001 DAM001    1.35
IND002 IND002 DAM001    1.28
IND003 IND003 DAM002    1.42
```

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--maternal-model` | 开启母体效应模型 | — |
| `--maternal-id <col>` | 母体 ID 列号 (1-indexed) | 3 |
| `--maternal-snp-test <mode>` | 模型专属 SNP 检验：`direct`（子代基因型，默认）/ `maternal`（母亲基因型单列）/ `joint`（direct+maternal 2-df 联合） | `direct` |

**方差组分**：直接加性 (Va)、母体加性 (Vm)、直接-母体协方差 (Cov_am)、母体永久环境 (Vpe)、残差 (Ve)。

**模型专属 SNP 检验（`--maternal-snp-test`，experimental）**：

```bash
# maternal：检验母亲基因型的母体 SNP 效应（H0: β_m=0）
./cosmicscanner --bfile data/my_geno --pheno data/birth_weight.txt \
  --maternal-model --maternal-id 3 --covar data/covar.txt --gwas \
  --maternal-snp-test maternal --out results/maternal_snp
# joint：逐 SNP 做 2-df Wald 检验 H0: β_d=β_m=0
./cosmicscanner --bfile data/my_geno --pheno data/birth_weight.txt \
  --maternal-model --maternal-id 3 --covar data/covar.txt --gwas \
  --maternal-snp-test joint --out results/joint_snp
```

- eligible 样本 = 有表型 **且** dam 为基因型个体（dam 必须是 PLINK 家系文件中的 subject）；非基因型 dam 被显式过滤，并在输出与 `*.maternal.snp.meta.txt` 中以 `N_TOTAL / N_WITH_DAM_GENOTYPE / N_DROPPED_MISSING_DAM` 报告。
- maternal 扫描以 dam 的 ALT dosage（`x_dam`，逐 SNP 中心化）进入设计；null 协方差与 maternal VCE 一致：`σ²ad·K(子代) + σ²am·K(母亲) + σ²e·I`，每染色体分解一次复用；LOCO 与 direct 扫描同策略。
- 默认 `direct` 与旧 `--maternal-model` 输出 byte-identical（回归已验证）。
- 输出：`direct` → `*.assoc.txt`；`maternal` → `*.maternal.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA_MATERNAL SE_MATERNAL P LOG10P` + N 列）；`joint` → `*.joint.assoc.txt`（另含 `BETA_DIRECT SE_DIRECT BETA_MATERNAL SE_MATERNAL COV_BETA_DM CORR_BETA_DM STAT_JOINT DF P_JOINT COND RANK STATUS`；信息矩阵秩亏/病态（如子代即母体）时 `status=COLLINEAR`）。
- BGEN/PGEN 尚未接入 `maternal|joint`（仅 PLINK，报错退出）。

**方法学定义（命名精确性）**：默认（`direct`）模型是 **maternal-adjusted GWAS（母体校正 GWAS）**。随机效应背景为 `Xb + Z_d a_d + Z_m a_m + e`（直接与母体加性分量；当前 `G_m = G_d`，`Cov(a_d,a_m)` 尚未估计）。SNP 扫描使用**个体自身基因型**作单一 direct 固定 SNP 效应（`x_j β_j`），检验 `H0: β_j = 0`，输出与标准 LMM 完全一致的 `BETA/SE/P`。`--maternal-snp-test maternal|joint` 令母亲基因型进入 SNP 设计（母体 SNP 效应 `x_dam β_m` 与 2-df 联合检验）。软件状态为 `Experimental validated`，证据为包含 null 校准和 direct/maternal/joint 真值恢复的 10 场景模拟套件；不声称外部 maternal-SNP reference。详见 [GWAS_MODEL_DEFINITIONS.md](GWAS_MODEL_DEFINITIONS.md#5-maternal-model--maternal-model)。

**输出文件**：`*.assoc.txt`；`maternal`/`joint` 模式另有 `*.maternal.assoc.txt` / `*.joint.assoc.txt` 与 `*.maternal.snp.meta.txt`

### 6.6 单步 ssGWAS

通过标准 A22-corrected H⁻¹ 矩阵联合系谱与基因组信息，支持系谱-基因组联合关联扫描：

**方法学定义（命名精确性）**：提供两条路径。

- **Method 1 — H-matrix single-step GWAS**（默认；`Experimental validated`）。以个体基因型 `x_j` 为设计、`H⁻¹ = A⁻¹ + [0; τG⁻¹ − ωA22⁻¹]` 导出的 H 子块为背景协方差（分析样本），检验固定 SNP 效应 `H0: β_j = 0`，输出 `BETA/SE/P`。
- **Method 2 — ssGBLUP back-solving GWAS**（`--single-step-method backsolve`；`Experimental validated`）。对全系谱（含未基因型动物）求解单步 MME 得 EBV/GEBV，再反解 SNP 效应 `ĝ = Z'G⁻¹ û_g`（VanRaden 2008; Strandén-Garrick 2009; Wang 2012）并输出滑动窗口遗传方差（`%Var_window`，10 SNP 窗口、不跨染色体）。已通过解析门（EBV max diff 5e-7，对封闭解析参考）与模拟门（Top10=5/5 因果恢复；系谱行序洗牌逐位一致）。未实现 iterative weighted ssGBLUP。
- **Method 2 可选 SNP P 值**：`exact` 为 `Stable`；`apy-exact` 为 `Experimental validated` scalable APY uncertainty path（需 core）；`apy-approx` 为 `Experimental approximate`（`STRICT_TYPE_I_CONTROL=NOT_GUARANTEED`，raw P-values only，无 genomic control/rescale，且 scalability 仅指 uncertainty layer），不得描述为 calibrated inference。

详见 [GWAS_MODEL_DEFINITIONS.md](GWAS_MODEL_DEFINITIONS.md#8-ssgwas--single-step)。

back-solving 命令：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --single-step --single-step-method backsolve \
  --ped data/pedigree.txt \
  --genotyped-list data/genotyped_ids.txt \
  --reml-vg 0.5 --reml-ve 1.0 \
  --out results/ssgwas_bs
```

输出：`*.ss_backsolve.ebv.txt`（全部动物 EBV/GEBV）、`*.ss_backsolve.snp.txt`（逐 SNP 效应：`EFFECT_ALLELE = (z'q/σ_z)/S`（等位替换刻度）、`EFFECT_STD = EFFECT_ALLELE·σ_z`（标准化刻度，满足 Z·STD = M·ALLELE = GEBV）、`SNP_VAR = 2p(1−p)·EFFECT_ALLELE²`）、`*.ss_backsolve.window.txt`（WINDOW_VARIANCE = 窗口内 Σ 2p(1−p)·EFFECT_ALLELE²，定义 A）、`*.ss_backsolve.meta.txt`（含 SNP_VARIANCE_DEFINITION/WINDOW_VARIANCE_DEFINITION/PCT_VARIANCE_DENOMINATOR）。

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --single-step \
  --ped data/pedigree.txt \
  --genotyped-list data/genotyped_ids.txt \
  --gwas \
  --out results/ssgwas_out
```

**系谱文件格式**（`--ped`，无表头）：

```
ID       SIRE     DAM
IND001   SIRE001  DAM001
IND002   SIRE001  DAM002
```

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--single-step` | 开启 ssGWAS 模型 | — |
| `--single-step-method lmm\|backsolve` | Method 1（H 核 LMM 关联）或 Method 2（育种值反解） | `lmm` |
| `--make-agh` | 仅构建并输出 A⁻¹、G、H⁻¹，不要求表型 | — |
| `--write-agh` | ssGWAS 同时输出 A⁻¹、G、H⁻¹ | 关闭 |
| `--ped <file>` | 系谱文件：`ID SIRE DAM` | — |
| `--genotyped-list <file>` | 有基因型个体 ID 列表 | — |
| `--h2-tau <f>` | A22 混合的 Tau 缩放（backsolve 时随 `--h2-omega` 固定 τ=h²、ω=(1−h²)τ） | 1.0 |
| `--h2-omega <f>` | A22 混合的 Omega 缩放 | 1.0 |
| `--blend-g <f>` | G 与 I 的混合权重 | 0.05 |
| `--ssgwas-window-snps <n>` | 滑动窗口 SNP 数（backsolve 的 %Var_window） | 10 |
| `--ssgwas-solver ldlt\|pcg` | MME 求解器（backsolve） | 自动 |
| `--ssgwas-pvalue-method off\|exact\|apy-exact\|apy-approx` | 可选 backsolve SNP P-value 路径；`apy-approx` 为 experimental approximate，非 strict type-I controlled | `off` |
| `--ssgwas-apy-core-size <n>` | `apy-exact` / `apy-approx` 所需 APY core size | — |

**矩阵-only 命令**：

```bash
./cosmicscanner --bfile data/my_geno --make-agh \
  --ped data/pedigree.txt --genotyped-list data/genotyped_ids.txt \
  --out results/agh
```

矩阵输出为 `*.Ainv.bin/.id`、`*.G.bin/.id`、`*.Hinv.bin/.id` 和 `*.agh.meta.txt`。ssGWAS 输出为 `*.assoc.txt` 和 `*_ssgwas_null.txt`。当前参与 VCE/扫描的表型个体必须具有基因型；额外系谱个体可作为祖先结构。

### 6.7 阈值/Probit 模型 GWAS

支持 probit/logit 链接函数的二分类性状 PQL 关联分析：

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

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--threshold-model` | 开启二元性状 PSR/probit 模型 | — |
| `--link probit\|logit` | 链接函数 | `probit` |
| `--spa` | 鞍点近似检验（改善二元性状校准） | 关闭 |

> **建议**：二元性状推荐加 `--spa` 改善校准。

### 6.8 计数/二项 PSR-GWAS

支持 Poisson、Binomial (event/trial) 性状的 PSR/PQL 近似 + weighted LMM 扫描。

**Poisson 计数性状**：

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

`--offset` 文件为每个样本一行 log-offset/exposure 值，顺序与 genotype sample 顺序一致。

**Binomial 比例性状**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --count-model \
  --family binomial \
  --binomial-pheno data/binomial_event_trial.pheno \
  --out results/psr_binomial
```

> **注意**：Binomial 请使用 event/trial 两列，不要只输入比例值。

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--count-model` | 开启计数/binomial PSR-GWAS | — |
| `--family poisson\|binomial\|negbin` | 分布族 | `poisson` |
| `--offset <file>` | Poisson log-offset | — |
| `--binomial-pheno <file>` | Binomial event/trial 文件 | — |

> **注意**：`--family negbin` 仍为实验入口，overdispersion 估计未完成系统验证，暂不建议作为正式结论。

### 6.9 序数 PSR-GWAS

支持 3+ 有序分类性状，提供两种引擎：

**推荐：stacked proportional-odds PSR（dense-zkz + grid-reml）**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --ordinal-model \
  --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode stacked-psr \
  --ordinal-stacked-engine dense-zkz \
  --ordinal-vc-method grid-reml \
  --out results/psr_ordinal
```

**基线：cumulative-binary**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --gwas \
  --ordinal-model \
  --ordinal-pheno data/ordinal.pheno \
  --ordinal-mode experimental-cumulative-binary \
  --out results/psr_ordinal_baseline
```

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--ordinal-model` | 开启序数 PSR-GWAS | — |
| `--ordinal-pheno <file>` | 序数表型文件 `FID IID ORDINAL` | — |
| `--ordinal-mode <name>` | `experimental-cumulative-binary` 或 `stacked-psr` | `experimental-cumulative-binary` |
| `--ordinal-stacked-engine <name>` | `fixed-gls` 或 `dense-zkz` | `fixed-gls` |
| `--ordinal-vc-method <name>` | `fixed` 或 `grid-reml` | `fixed` |
| `--ordinal-test <name>` | `wald`/`score`/`robust-score`/`all` | `wald` |
| `--ordinal-allow-fallback` | stacked-psr 失败时回退到 cumulative-binary | 关闭 |
| `--ordinal-sigma2-g <f>` | 固定遗传方差 | — |
| `--ordinal-vc-grid-min <f>` | grid-reml 搜索下界 | — |
| `--ordinal-vc-grid-max <f>` | grid-reml 搜索上界 | — |
| `--ordinal-vc-grid-points <n>` | grid-reml 搜索点数 | — |

**输出文件**：
- `*.gwas.tsv`：主结果（修正后的 Wald / pseudo-score）
- `*.ordinal.score.tsv`：诊断输出（`--ordinal-test all` 时）

### 6.10 生存 Cox 模型 GWAS

基于 Cox 比例风险模型与鞅残差 LMM 的生存性状关联分析：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/survival.pheno \
  --survival-model \
  --covar data/covar.txt \
  --gwas \
  --out results/survival_gwas
```

**表型文件格式**（`FID IID TIME STATUS`）：

```
FID    IID    TIME    STATUS
IND001 IND001 365     1
IND002 IND002 720     0
IND003 IND001 180     1
```

- `TIME`：随访时间或生存时间。
- `STATUS`：1 = 事件发生，0 = 删失。

> **实现说明**：当前生存模型为简化实现——使用 Nelson-Aalen 累积危险估计鞅残差，再对残差执行 fast-LMM 扫描，并非完整 Cox PH 似然求解。适合初筛，显著位点建议用 R `coxme`/`GMMAT` 复核。

**输出文件**：`*.assoc.txt`

### 6.11 GxE 反应范数模型 GWAS

建模随环境梯度变化的 SNP 效应（Vg0+Vg1 双方差组分），提供 2-df 联合 Wald 检验：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/gxe_pheno.txt \
  --gxe-model \
  --env-covar 3 \
  --covar data/covar.txt \
  --gwas \
  --out results/gxe_gwas
```

**表型文件格式**（`FID IID ENV PHENO`）：

```
FID    IID    ENV    PHENO
IND001 IND001 25.3   1.23
IND002 IND002 30.1   1.45
```

- `--env-covar <col>`：环境变量在表型文件中的列号 (1-indexed)。

**输出文件**：`*.gxe.assoc.txt`

**输出列**：`CHR SNP BP A1 A2 FREQ N BETA_G SE_G BETA_GxE SE_GxE JOINT_CHISQ JOINT_P`

### 6.12 罕见变异分析 (SKAT / Burden)

基于基因或自定义集合进行罕见变异关联检验：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --set-file data/gene_sets.txt \
  --skat \
  --gwas \
  --out results/skat_gwas
```

**Burden 检验**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --set-file data/gene_sets.txt \
  --burden \
  --gwas \
  --out results/burden_gwas
```

> **限制**：SKAT/Burden 仅在默认单性状 LMM 路径下可用，不能与 `--multivariate`、`--repeatability-model` 等其他模型标志同时使用。

**集合文件格式**（`--set-file`）：

```
GENE1    rs123    rs456    rs789
GENE2    rs234    rs567
```

**二元性状 SPA（鞍点近似）**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/binary.pheno \
  --gwas \
  --spa \
  --out results/spa_gwas
```

### 6.13 Biobank 级别群体分析

针对数十万大样本，开启 REGENIE 风格的 Ridge Regression 两步法：

**完整两步**：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --block-size 1000 \
  --gwas \
  --out results/biobank_gwas
```

**分步运行**：

```bash
# Step 1: 生成 LOCO 预测
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --regenie-step 1 \
  --write-pred results/loco_pred.txt \
  --gwas \
  --out results/biobank_step1

# Step 2: 使用预测进行关联分析
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --biobank-mode \
  --regenie-step 2 \
  --use-pred results/loco_pred.txt \
  --gwas \
  --out results/biobank_step2
```

> **限制**：Biobank 模式仅在默认单性状 LMM 路径下可用，不能与其他模型标志同时使用。

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--biobank-mode` | 开启 REGENIE 风格两步法 | — |
| `--block-size <n>` | Ridge 分块大小 | 1000 |
| `--regenie-model-snps <file>` | Step 1 建模 SNP 列表 | — |
| `--write-pred <file>` | 输出 Step 1 LOCO 预测 | — |
| `--use-pred <file>` | 使用已有 LOCO 预测 | — |
| `--regenie-step <n>` | 0=两步, 1=仅预测, 2=仅关联 | 0 |

**输出文件**：`*.biobank.assoc.txt`

### 6.14 经验贝叶斯 EB-GWAS

实验性 random-SNP effect 信号增强，支持单性状/重复力/多性状：

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --mpheno 1 \
  --model lmm-eb \
  --test eb-wald \
  --adjust eb-mbonf \
  --gwas \
  --out results/eb_gwas
```

**关键参数**：

| 参数 | 说明 |
|------|------|
| `--model lmm-eb` | 开启 EB-LMM（别名 `--gwas-method eb-lmm`） |
| `--test eb-wald` | 启用 EB Wald 检验 |
| `--adjust eb-mbonf` | 输出 `.eb.summary.txt` 中的实验性阈值 |

> **注意**：EB-GWAS 是稳定辅助路线，但正式分析请同时报告标准 LMM 结果和普通 Bonferroni/FDR。

详细说明见 [docs/EB_GWAS_MODULE_CN.md](EB_GWAS_MODULE_CN.md)。

---

### 6.15 系谱动物模型 GWAS

用系谱亲缘关系矩阵 A 替代基因组 GRM 描述多基因背景遗传效应：

```
y = Xb + x_j * beta_j + Z*a + e,   a ~ N(0, A * sigma2_a)
```

```bash
./cosmicscanner \
  --bfile data/my_geno \
  --pheno data/my_pheno.txt \
  --animal-model \
  --ped data/pedigree.txt \
  --gwas \
  --out results/animal_assoc
```

**系谱文件**（`--ped <f>`，空白分隔，无表头）：

| 列 | 内容 |
|----|------|
| 1 | 个体 ID（必须包含全部分析样本） |
| 2 | 父本 ID（未知写 `NA`/`.`/`0`） |
| 3 | 母本 ID（未知写 `NA`/`.`/`0`） |

ID 匹配顺序：先 `FID:IID`，再裸 `IID`（当 `.fam` 的 FID 为占位符时回退）。系谱可包含额外未测序祖先；仅分析子集用于构建 `A22` 块。

**关键参数**：

| 参数 | 说明 |
|------|------|
| `--animal-model` | 启用系谱动物模型 |
| `--ped <f>` | 系谱文件 `ID SIRE DAM`（必需） |
| `--write-ped-kernel` | 导出分析子集的精确 A22 到 `*.animal_model.k.txt` + `*.animal_model.k.id`（可直接给 GEMMA `-k` 用） |
| `--no-loco` | 关闭 LOCO（多染色体数据需加此参数以消除提示） |

**输出文件**：
- `*.assoc.txt`：`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`
- `*.animal_model.null.txt`：`sigma2_a / sigma2_e / h2` 及系谱解析/`Ainv`/`A22`/VCE/扫描耗时
- `*.animal_model.k.txt` / `*.animal_model.k.id`：精确稠密 A22 核与样本顺序（仅 `--write-ped-kernel` 时输出）

**说明**：
- 动物模型与 ssGWAS 共享同一 `PedigreeRelationship` 系谱解析和 `Ainv`/`A22` 构建（`buildA22` 经由 Ainv 的 Schur complement 求解，不物化稠密 A）。ssGWAS 加 `--h2-tau 0 --h2-omega 0 --no-loco` 后即退化为此模型，逐 SNP 的 BETA/SE/-log10P 完全一致。
- LOCO 不适用（只有一个关系矩阵、无法按染色体分离）；多染色体数据会打印提示，除非加 `--no-loco`。
- GEMMA 外部对标：GEMMA 0.98.5 可用 `-k` 读取稠密 `*.animal_model.k.txt` 并运行 `-lmm 4`；无需 GRM 归一化（刻度被吸收进 `sigma2_g`）。
- 验证：`tests/run_test_animal_model.sh`（模拟套件）和 `tests/run_animal_model_gemma_validation.sh`（外部 GEMMA `-k` 基准；无 GEMMA 时 SKIP）。

---

## 7. 输出文件格式说明

### 7.1 标准 LMM GWAS — `*.assoc.txt`

| 列名 | 说明 |
|------|------|
| CHR | 染色体 |
| SNP | SNP ID |
| BP | 物理位置 (bp) |
| A1 | 效应等位基因 |
| A2 | 非效应等位基因 |
| FREQ | A1 等位基因频率 |
| N | 样本数 |
| BETA | 效应估计值 |
| SE | 标准误 |
| P | p 值 |
| LOG10P | -log10(P) |

**EB-LMM 模式额外列**：`EB_BETA`, `EB_SE`, `EB_WALD`

### 7.2 多性状 GWAS — `*.assoc.mv.txt`

| 列名 | 说明 |
|------|------|
| CHR | 染色体 |
| SNP | SNP ID |
| BP | 物理位置 |
| A1 / A2 | 等位基因 |
| FREQ | A1 频率 |
| WALD_STAT | Wald 统计量 |
| P | p 值 |

GEMMA-full 模式额外输出 LRT、Score 列。

### 7.3 RRM GWAS — `*.rrm_joint.assoc.txt`

| 列名 | 说明 |
|------|------|
| CHR / SNP / BP / A1 / A2 | 基本信息同上 |
| FREQ | A1 频率 |
| N | 样本数 |
| BETA_L0 | 基线效应 |
| SE_L0 | 基线标准误 |
| JOINT_CHISQ | 联合检验卡方统计量 |
| JOINT_P | 联合检验 p 值 |

### 7.4 GxE GWAS — `*.gxe.assoc.txt`

| 列名 | 说明 |
|------|------|
| CHR / SNP / BP / A1 / A2 / FREQ / N | 基本信息同上 |
| BETA_G | 主效应（基因型） |
| SE_G | 主效应标准误 |
| BETA_GxE | 交互效应 |
| SE_GxE | 交互效应标准误 |
| JOINT_CHISQ | 2-df 联合 Wald 卡方 |
| JOINT_P | 联合 p 值 |

### 7.5 序数 PSR-GWAS — `*.gwas.tsv`

默认输出修正后的 Wald / pseudo-score 主结果。`--ordinal-test all` 时额外输出 `*.ordinal.score.tsv` 诊断文件。

### 7.6 Biobank GWAS — `*.biobank.assoc.txt`

列格式同标准 `*.assoc.txt`。

### 7.7 EBI GWAS-SSF 格式

添加 `--out-ssf` 参数可输出 EBI GWAS-SSF 标准格式。

---

## 8. 完整参数说明

### 8.1 输入/输出

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--bfile <prefix>` | PLINK 1.9 二进制前缀 | — |
| `--bgen <file>` | BGEN 基因型文件 | — |
| `--sample <file>` | BGEN sample 文件 (Oxford) | — |
| `--pfile <prefix>` | PLINK 2.0 二进制前缀 | — |
| `--pgen <file>` / `--pvar <file>` / `--psam <file>` | PLINK 2.0 单独文件 | — |
| `--out <prefix>` | 输出文件前缀 | `cosmic_gwas` |
| `--out-ssf` | 输出 EBI GWAS-SSF 格式 | 关闭 |
| `--extract <file>` | 仅包含指定 SNP | — |
| `--exclude <file>` | 排除指定 SNP | — |

### 8.2 表型/协变量

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--pheno <file>` | 表型文件 | — |
| `--mpheno <n>` | 表型列号 (1-indexed) | 1 |
| `--mpheno-list <list>` | 多性状列号 (如 `1,2,3`) | — |
| `--tpheno <n>` | TIME 列号 (重复力/RRM) | — |
| `--covar <file>` | 离散型协变量 | — |
| `--qcovar <file>` | 连续型协变量 | — |
| `--long-pheno <file>` | 长格式表型（开启 RRM） | — |
| `--binomial-pheno <file>` | Binomial event/trial 文件 | — |
| `--offset <file>` | Poisson log-offset | — |

### 8.3 位点过滤

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--maf <f>` | 最小 MAF | 0.01 |
| `--hwe <f>` | HWE p 值过滤 | 0.0（关闭） |
| `--chr <chr>` | 限制单染色体 | — |
| `--from-bp <pos>` | 扫描起始位置 | — |
| `--to-bp <pos>` | 扫描结束位置 | — |

### 8.4 LMM 引擎

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--gwas` | 开启 GWAS 扫描（**所有模型均需显式指定**） | 关闭 |
| `--no-loco` | 禁用 LOCO | 启用 LOCO |
| `--fast-lmm` | 固定 delta 近似加速 | 关闭 |
| `--lmm-null-solver <s>` | null 求解器：`brent`/`idul`/`idul_plus`/`aidul` | `brent` |
| `--no-covar-proj` | 禁用协变量投影 | 启用 |
| `--force-slow` | 强制慢速路径 | 关闭 |
| `--spa` | 鞍点近似 (二元性状) | 关闭 |

> **重要**：`--gwas` 是触发扫描的必要开关。仅指定 `--repeatability-model`、`--multivariate`、`--maternal-model` 等模型标志而不加 `--gwas`，程序会正常退出但不输出关联结果。例外：`--rrm-unbalance` 会自动隐式启用 `--gwas`。

### 8.5 REML / VCE 控制

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--reml-maxit <n>` | VCE 最大迭代次数 | 50 |
| `--reml-tol <f>` | VCE 收敛阈值 | 1e-6 |
| `--reml-vg <f>` | 固定遗传方差（跳过 VCE） | — |
| `--reml-vpe <f>` | 固定永久环境方差 | — |
| `--reml-ve <f>` | 固定残差方差 | — |

### 8.6 多性状 GWAS

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--multivariate` | 开启多性状 GWAS | — |
| `--mv-method <m>` | `emmv`/`wald`/`eb`/`gemma`/`gemma-full`/`hybrid` | `emmv` |
| `--mv-compatible <name>` | 兼容模式（`emmvgwas`） | — |
| `--mv-test <list>` | GEMMA-full 检验：`wald,lrt,score` | 全开 |
| `--mv-refine-top <n>` | 精修 SNP 数上限 | 100 |
| `--mv-refine-p <p>` | 精修 p 阈值 | 1e-3 |
| `--mv-refine-score-prepass <n>` | Score prepass top n | — |
| `--mv-refine-discrepancy-top <n>` | 差异 SNP 精修数 | — |
| `--mv-refine-discrepancy <x>` | 最小 \|log10P 差异\| | 1.0 |
| `--mv-refit-p <p>` | GEMMA-full refit 阈值 | 1.0 |
| `--mv-vce-solver <s>` | `em` 或 `aireml` | `em` |
| `--mv-em-min-iter <n>` | EM 最小迭代 | 10 |
| `--mv-ai-max-iter <n>` | AI-REML 最大迭代 | 50 |
| `--mv-ai-tol <f>` | AI-REML 收敛阈值 | 1e-6 |

### 8.7 重复力模型

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--repeatability-model` | 开启重复力模型 | — |
| `--repeatability-mode <m>` | `auto`/`exact`/`fast`/`stcg` | `auto` |
| `--repeatability-grm-cache <prefix>` | GRM 磁盘缓存 | — |
| `--repeatability-mixed-precision` | Float32 降低内存 | 关闭 |
| `--repeatability-weight <file>` | 记录权重 | — |
| `--repeatability-refit-validation` | 逐 SNP VCE refit | 关闭 |
| `--repeatability-refit-limit <n>` | refit 上限 (0=全部) | 0 |
| `--vce-mode <m>` | `auto`/`fast`/`subject-level-exact`/`full-record-exact` | `auto` |

### 8.8 随机回归模型 (RRM)

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--rrm-order <n>` | 基函数阶数 (1-4) | — |
| `--rrm-basis <m>` | `legendre` 或 `bspline` | `legendre` |
| `--rrm-curve-grid <n>` | 曲线网格点数 | 21 |
| `--rrm-top <n>` | top SNP 曲线数 | 100 |
| `--rrm-ai-iter <n>` | AI 精修迭代 (-1=auto, 0=fast) | -1 |
| `--rrm-sig-threshold <p>` | 动态时间窗口 p 阈值 | 0.05 |
| `--rrm-unbalance` | 非平衡纵向模式 | 关闭 |
| `--rrm-loco-global-null` | 全局 null 快速 LOCO | 关闭 |
| `--tpoint <col>` | 非平衡模式时间列名 | — |
| `--tcovar <cols>` | 时间依赖连续协变量 | — |
| `--tclass <cols>` | 时间依赖离散协变量 | — |
| `--rrm-compatible <m>` | `gma-fixreg` 兼容模式 | — |
| `--rrm-compatible-varcom <file>` | GMA `.var` 文件 | — |

### 8.9 母体效应模型

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--maternal-model` | 开启母体效应模型 | — |
| `--maternal-id <col>` | 母体 ID 列号 (1-indexed) | 3 |

### 8.10 单步 ssGWAS

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--single-step` | 开启 ssGWAS | — |
| `--ped <file>` | 系谱文件 `ID SIRE DAM` | — |
| `--genotyped-list <file>` | 有基因型 ID 列表 | — |
| `--write-ebv` | 预留开关（legacy）；`--single-step-method backsolve` 始终输出 `*.ss_backsolve.ebv.txt` | 关闭 |
| `--h2-tau <f>` | Tau 缩放 | 1.0 |
| `--h2-omega <f>` | Omega 缩放 | 1.0 |
| `--blend-g <f>` | G 与 I 混合权重 | 0.05 |

### 8.11 阈值/Probit 模型

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--threshold-model` | 开启二元 PSR/probit | — |
| `--link probit\|logit` | 链接函数 | `probit` |

### 8.12 计数/二项 PSR-GWAS

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--count-model` | 开启计数/binomial PSR | — |
| `--family poisson\|binomial\|negbin` | 分布族 | `poisson` |
| `--offset <file>` | Poisson log-offset | — |
| `--binomial-pheno <file>` | Binomial event/trial | — |

### 8.13 序数 PSR-GWAS

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--ordinal-model` | 开启序数 PSR | — |
| `--ordinal-pheno <file>` | 序数表型文件 | — |
| `--ordinal-mode <name>` | `experimental-cumulative-binary` / `stacked-psr` | `experimental-cumulative-binary` |
| `--ordinal-stacked-engine <name>` | `fixed-gls` / `dense-zkz` | `fixed-gls` |
| `--ordinal-vc-method <name>` | `fixed` / `grid-reml` | `fixed` |
| `--ordinal-test <name>` | `wald`/`score`/`robust-score`/`all` | `wald` |
| `--ordinal-allow-fallback` | 回退到 cumulative-binary | 关闭 |
| `--ordinal-sigma2-g <f>` | 固定遗传方差 | — |
| `--ordinal-vc-grid-min <f>` | grid 搜索下界 | — |
| `--ordinal-vc-grid-max <f>` | grid 搜索上界 | — |
| `--ordinal-vc-grid-points <n>` | grid 搜索点数 | — |

### 8.14 PSR 迭代控制

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--psr-maxit <n>` | PSR/PQL 最大迭代 | 50 |
| `--psr-tol <f>` | 收敛阈值 | 1e-6 |
| `--psr-mu-eps <f>` | mu clipping epsilon | 1e-6 |
| `--psr-min-weight <f>` | working weight 下界 | 1e-8 |
| `--psr-max-weight <f>` | working weight 上界 | 1e8 |

### 8.15 生存 Cox 模型

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--survival-model` | 开启 Cox 模型 | — |

### 8.16 GxE 反应范数模型

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--gxe-model` | 开启 GxE 模型 | — |
| `--env-covar <col>` | 环境变量列号 (1-indexed) | — |

### 8.17 罕见变异

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--skat` | SKAT 检验 | — |
| `--burden` | Burden 检验 | — |
| `--set-file <file>` | 变异集合定义 | — |

### 8.18 Biobank 模式

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--biobank-mode` | REGENIE 风格两步法 | — |
| `--block-size <n>` | Ridge 分块大小 | 1000 |
| `--regenie-model-snps <file>` | Step 1 建模 SNP | — |
| `--write-pred <file>` | 输出 LOCO 预测 | — |
| `--use-pred <file>` | 使用已有预测 | — |
| `--regenie-step <n>` | 0=两步, 1=预测, 2=关联 | 0 |

### 8.19 GRM 构建

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--make-grm` | 仅计算 GRM 并退出 | — |
| `--make-grm-alg <n>` | 1=VanRaden, 2=Zeng, 3=Yang, 4=Vitezica, 5=TrueVar, 6=EMmvGWAS-compatible | 1 |
| `--make-grm-sparse` | 稀疏 GRM | 关闭 |
| `--sparse-cutoff <f>` | 稀疏截断值 | 0.05 |
| `--snp-weight <file>` | SNP 权重 | — |
| `--pop-class <file>` | 群体分层分类 | — |
| `--force-slow` | 强制 float 慢速路径（验证用） | 关闭 |

> **说明**：`--make-grm-alg 6` 仅在 `--mv-compatible emmvgwas` 时启用（由多性状模型自动设置），单独使用 `--make-grm --make-grm-alg 6` 不会触发 EMmvGWAS 编码。

### 8.20 EB-GWAS

| 参数 | 说明 |
|------|------|
| `--model lmm-eb` | 开启 EB-LMM（别名 `--gwas-method eb-lmm`） |
| `--test eb-wald` | EB Wald 检验 |
| `--adjust eb-mbonf` | 实验性阈值 |

### 8.21 输入适配层

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input-file <file>` | 原始 csv/tsv/txt 表格 | — |
| `--input-config <json>` | 从 JSON 配置读取 | — |
| `--write-template-config <json>` | 生成模板配置并退出 | — |
| `--input-format <fmt>` | `auto`/`csv`/`tsv`/`txt` | `auto` |
| `--input-model <model>` | `repeatability`/`rrm`/`lmm`/`multivariate` | `repeatability` |
| `--prepared-input-dir <dir>` | prepared input 目录 | — |
| `--prepare-input-only` | 仅准备输入并退出 | — |
| `--id-col <name>` | ID 列名 | — |
| `--iid-col <name>` | IID 列名（当 id-col 为 FID） | — |
| `--time-col <name>` | 时间列名 | — |
| `--trait-cols <list>` | 性状列名（逗号分隔） | — |
| `--input-trait <name>` | 准备后运行的性状 | 第一个 |
| `--record-qcovars <list>` | 记录层连续协变量 | — |
| `--record-ccovars <list>` | 记录层离散协变量 | — |
| `--subject-qcovars <list>` | 个体层连续协变量 | — |
| `--subject-ccovars <list>` | 个体层离散协变量 | — |
| `--time-parser <mode>` | `auto`/`date`/`integer`/`batch` | `auto` |
| `--time-unique-policy <p>` | `strict`/`within-time-index` | `within-time-index` |
| `--time-index-multiplier <n\|auto>` | TIME 编码倍率 | `auto` |
| `--preview-schema` | 输出 schema 预览并退出 | — |
| `--infer-columns` | 列推断建议 | — |
| `--wide-to-long` | 宽表转长表 | — |

### 8.22 工程优化

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--threads <n>` | OpenMP 线程数 | 最大可用 |
| `--omp <n>` | `--threads` 别名 | 最大可用 |

---

## 9. 性能调优指南

### 9.1 加速 null 模型求解

```bash
# 使用 IDUL 加速（推荐中等规模数据）
--lmm-null-solver idul

# 使用 IDUL+ 或 AIDUL（推荐大规模数据）
--lmm-null-solver aidul
```

详见 [docs/IDUL_AIDUL_ROADMAP.md](IDUL_AIDUL_ROADMAP.md)。

### 9.2 重复力模型优化

```bash
# 启用 GRM 磁盘缓存复用（大幅加速重复扫描）
--repeatability-grm-cache grm_cache

# 使用 Float32 降低内存（大规模数据）
--repeatability-mixed-precision
```

### 9.3 多线程并行

```bash
# 指定线程数
--threads 16
```

### 9.4 Biobank 大规模数据

```bash
# 使用 Biobank 两步法处理百万级样本
--biobank-mode --block-size 1000
```

### 9.5 快速筛查

```bash
# 使用 fast-lmm 近似加速
--fast-lmm

# 重复力模型使用 fast 模式
--repeatability-mode fast
```

---

## 10. 常见问题 (FAQ)

### Q1: PSR-GWAS 收敛警告是什么意思？

Convergence warning 是诊断信息，**不等同于程序失败**：

| `status_class` | 含义 |
|----------------|------|
| `CONVERGED` | 正常收敛 |
| `CONVERGED_WITH_WARNING` | 收敛但存在边界方差、quasi-likelihood 非单调或 clipping |
| `NOT_CONVERGED` | 未收敛，正式分析需谨慎 |

**Warning flags**：

| Flag | 含义 |
|------|------|
| `WARNING_PSR_MAXIT_REACHED` | 达到迭代上限 |
| `WARNING_PSR_CLIPPED` | mu 或 working weight 触及 clipping 边界 |
| `WARNING_PSR_SIGMA_BOUNDARY` | 遗传方差接近 0 |
| `WARNING_PSR_QUASILIK_DECREASED` | quasi-likelihood proxy 下降 |

### Q2: 遇到难收敛或极端稀疏性状怎么办？

调整 PSR/PQL 迭代与 clipping 参数：

```bash
--psr-maxit 100 \
--psr-tol 1e-8 \
--psr-mu-eps 1e-8 \
--psr-min-weight 1e-10 \
--psr-max-weight 1e10
```

### Q3: PSR-GWAS 显著位点需要复核吗？

**需要**。PSR-GWAS 不是 exact GLMM，显著位点建议用 candidate refinement 或外部 exact GLMM (GMMAT/POLMM) 复核。

### Q4: 如何选择多性状方法？

- **大规模全基因组**：`--mv-method hybrid`（推荐）
- **小样本精准验证**：`--mv-method gemma-full`
- **快速扫描**：`--mv-method emmv`（默认）

### Q5: 如何与 EMmvGWAS R 包对标？

使用 `--mv-compatible emmvgwas --no-loco`，会启用 EMmvGWAS-compatible GRM（`--make-grm-alg 6`）。

### Q6: Binomial 性状为什么不能只输入比例值？

只输入比例值会丢失 trial 大小信息，导致方差估计错误。必须使用 event/trial 两列格式。

### Q7: 如何输出 EBI GWAS-SSF 标准格式？

添加 `--out-ssf` 参数。

### Q8: 输入适配层支持 Excel 吗？

当前 MVP 不支持 `xlsx/xls`，请先另存为 `csv/tsv`。

### Q9: 如何只计算 GRM 不做 GWAS？

```bash
./cosmicscanner --bfile data/my_geno --make-grm --out results/grm
```

### Q10: 为什么指定了模型标志却没有输出关联结果？

**必须显式加 `--gwas` 才会执行扫描**。例如以下命令不会输出结果：

```bash
# 错误：缺少 --gwas
./cosmicscanner --bfile geno --pheno pheno.txt --repeatability-model --out out
```

正确写法：

```bash
./cosmicscanner --bfile geno --pheno pheno.txt --repeatability-model --gwas --out out
```

例外：`--rrm-unbalance` 会自动隐式启用 `--gwas`。

### Q11: 如何限制只扫描一条染色体？

```bash
--chr 1
```

或限制区间：

```bash
--chr 1 --from-bp 1000000 --to-bp 2000000
```

---

## 11. 相关文档

| 文档 | 说明 |
|------|------|
| [Tutorial_EN.md](Tutorial_EN.md) | 英文教程 |
| [PSR_GWAS_ROADMAP.md](PSR_GWAS_ROADMAP.md) | PSR-GWAS 路线图 |
| [PSR_GWAS_VALIDATION_REPORT_V0.5.md](PSR_GWAS_VALIDATION_REPORT_V0.5.md) | PSR-GWAS v0.5 验证报告 |
| [PSR_GWAS_FULL_ORDINAL_DESIGN.md](PSR_GWAS_FULL_ORDINAL_DESIGN.md) | 完整序数 PSR 设计 |
| [PSR_GWAS_BENCHMARK_GUIDE.md](PSR_GWAS_BENCHMARK_GUIDE.md) | PSR-GWAS 基准测试指南 |
| [tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md](tutorials/PSR_GWAS_ORDINAL_TUTORIAL.md) | 序数 PSR 专项教程 |
| [REPEATABILITY_GWAS_RELEASE_CN.md](REPEATABILITY_GWAS_RELEASE_CN.md) | 重复力 GWAS 发布说明 |
| [REPEATABILITY_GWAS_DEVELOPMENT_SUMMARY_20260529_CN.md](REPEATABILITY_GWAS_DEVELOPMENT_SUMMARY_20260529_CN.md) | 重复力 GWAS 全景开发总结 |
| [REPEATABILITY_GWAS_MANUSCRIPT_CN.md](REPEATABILITY_GWAS_MANUSCRIPT_CN.md) | 重复力 GWAS 论文草稿 |
| [RRM_GWAS_VALIDATION_CN.md](RRM_GWAS_VALIDATION_CN.md) | RRM 方法与验证 |
| [RRM_GWAS_INPUT_DESIGN_CN.md](RRM_GWAS_INPUT_DESIGN_CN.md) | RRM 输入复用设计 |
| [IDUL_AIDUL_ROADMAP.md](IDUL_AIDUL_ROADMAP.md) | IDUL/AIDUL 加速路线图 |
| [EB_GWAS_MODULE_CN.md](EB_GWAS_MODULE_CN.md) | EB-GWAS 模块 |
| [INPUT_ADAPTATION_SPEC_CN.md](INPUT_ADAPTATION_SPEC_CN.md) | 输入适配层规范 |
| [ANIMAL_MODEL_GWAS.md](ANIMAL_MODEL_GWAS.md) | 系谱动物模型 GWAS 说明与验证 |
| [SUPPORTED_MODES_CN.md](SUPPORTED_MODES_CN.md) | 支持模式与核矩阵总表 |
| [models/NEW_MODELS_DESIGN.md](models/NEW_MODELS_DESIGN.md) | 新模型设计 |
| [GWAS_RELEASE_CLOSEOUT_CN.md](GWAS_RELEASE_CLOSEOUT_CN.md) | GWAS 发布收口验收记录 |
