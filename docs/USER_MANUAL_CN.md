# CosmicScanner v1.1.3 用户手册

## 概述

CosmicScanner v1.1.3 提供四项核心功能：

- **单性状 GWAS**：针对单个连续性状的线性混合模型关联分析
- **多性状 GWAS**：多个相关连续性状的联合关联检验
- **重复力 GWAS**：重复记录表型的混合模型关联分析
- **A/G/H 矩阵构建**：基于系谱的 A⁻¹、G 和 H⁻¹ 矩阵构建

## 安装

### 系统要求
Linux x86_64、macOS 或 Windows x86_64；源码构建需要 CMake 3.15+、
C++17 编译器、Eigen3、OpenMP、zlib 和 zstd。

### 源码构建
```bash
git clone --branch v1.1.3 --depth 1 https://github.com/huangjie2025/CosmicScanner.git
cd CosmicScanner
bash scripts/build.sh
```

GitHub Release 提供 Linux、macOS 和 Windows 压缩包。解压后运行
`bin/cosmicscanner`，Windows 为 `bin/cosmicscanner.exe`。Windows 包包含
UCRT64 运行时 DLL，macOS 包包含所需的非系统 dylib。所有包均由 GitHub
Actions 在对应平台构建并测试。

## 输入格式

### 基因型
PLINK 1 BED: `--bfile <前缀>`。PLINK 2 PGEN: `--pfile <前缀>`。BGEN: `--bgen <文件> --sample <文件>`。

### 表型
单性状/多性状: `FID IID 性状1 [性状2 ...]`。重复力: `FID IID 时间 表型`。缺失值: NA 或 -9。

### 协变量
`--qcovar`: 定量协变量。`--covar`: 分类协变量（虚拟编码）。

## 主要参数

| 参数 | 说明 | 必需 |
|---|---|---|
| --bfile / --pfile / --bgen | 基因型输入 | 三选一 |
| --pheno <文件> | 表型文件 | 是 |
| --out <前缀> | 输出前缀 | 是 |
| --gwas | 启用 GWAS | 或 --make-agh |
| --scan-tier approx\|semi\|exact | approx 固定 null；semi 固定比例/profile scale；exact 每 SNP 重估方差 | 省略时 approx/默认方法 |
| --no-loco | 使用全局 GRM，不做 LOCO | 多染色体时默认 LOCO |
| --lm 1\|2\|3\|4 | 线性模型 Wald/LRT/Score/all | 否 |
| --lmm-gemma-mode 1\|2\|3\|4 | LMM Wald/LRT/Score/all；4 输出 P_WALD/P_LRT/P_SCORE | 否 |
| --lmm-null-solver brent\|idul\|idul_plus\|aidul | 仅影响 null-model solver；IDUL 系列失败可回退 Brent | brent |
| --multivariate | 多性状 GWAS | 配合 --mpheno-list |
| --mv-method wald\|emmv\|gemma-compat\|eb | 推荐的多性状方法；gemma-compat 与 EB 为 experimental | emmv |
| --repeatability-model | 重复力 GWAS | 配合 --tpheno |
| --repeatability-mode auto\|exact\|fast\|stcg | 重复力模式 | auto |
| --make-agh | A/G/H 构建 | 配合 --ped |
| --maf / --hwe | 频率/HWE 过滤 | 否 |
| --extract / --exclude | SNP 选择 | 否 |
| --chr / --from-bp / --to-bp | 区间过滤 | 否 |
| --threads <n> | 线程数 | 系统默认 |

## 使用示例

### 单性状
```bash
cosmicscanner --bfile data --pheno trait.txt --mpheno 1 --gwas --scan-tier semi --out out
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier approx --out out.approx
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi --out out.semi
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --lmm-gemma-mode 4 --no-loco --out out.exact
```

### 多性状
```bash
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --mv-method wald --no-loco --out mv.wald
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --mv-method emmv --no-loco --out mv.emmv
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --mv-method gemma-compat --no-loco --out mv.gemma_compat
```
`gemma-compat` 目标是复现 GEMMA 默认有限 per-SNP REML EM Wald schedule；`gemma`、`gemma-default`、`gemma-full` 和 `gemma-exact-full` 仅作为 deprecated experimental aliases 保留一个版本。metadata 会记录 `requested_method`、`resolved_method` 和 `deprecated_alias`。`gemma-exact-full` 不应理解为 fully converged exact REML optimum；`reml-exact` 尚未公开。

### 重复力
```bash
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --out rep
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode fast --gwas --out rep.fast
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode stcg --gwas --out rep.stcg
```

### A/G/H
```bash
cosmicscanner --bfile data --make-agh --ped ped.txt --genotyped-list geno.txt --out agh
```

## 输出文件

| 文件 | 内容 |
|---|---|
| .assoc.txt | CHR, SNP, BP, A1, A2, FREQ, N, BETA, SE, P, LOG10P |
| .assoc.mv.txt | 多性状: WALD_STAT, P, BETA_t, SE_t |
| .repeatability.null.txt | Vg, Vpe, Ve, h2, r2, 时间统计 |
| .scan.meta.txt | release_version, git_commit, release_status, analysis_model, model, scan_tier, statistical_test, null_solver, loco, fallback, multivariate_method, compatibility_mode, vc_strategy |
| .Ainv.bin/.id | A⁻¹ 矩阵下三角 + ID |
| .G.bin/.id | G 矩阵下三角 + ID |
| .Hinv.bin/.id | H⁻¹ 矩阵下三角 + ID |
| .agh.meta.txt | 构建参数和验证信息 |

## 扫描档位

| 档位 | 策略 | 推荐场景 |
|---|---|---|
| approx | 固定 null 方差 | 高通量筛选 |
| semi | 固定比例, profile 残差 | 全基因组主扫描 |
| exact | 每 SNP 重估方差 | 小数据/候选位点 |

多性状/重复力 exact 要求 --no-loco；重复力 exact 额外要求 PLINK 输入和无记录权重。
Hybrid 和 GxE 多性状入口在 v1.1.3 focused release 中不支持，会在参数解析阶段报错。
完整 Stable/Experimental/Internal/Unsupported 矩阵见 docs/SUPPORTED_MODES.md。

## 版本

v1.1.3。引用信息见 CITATION.cff。
