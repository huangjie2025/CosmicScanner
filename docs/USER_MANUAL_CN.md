# CosmicScanner v1.1.2 用户手册

## 概述

CosmicScanner v1.1.2 提供四项核心功能：

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
git clone --branch v1.1.2 --depth 1 https://github.com/huangjie2025/CosmicScanner.git
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
| --scan-tier approx\|semi\|exact | 扫描档位 | 推荐 |
| --no-loco | 禁用 LOCO | LOCO 默认 |
| --multivariate | 多性状 GWAS | 配合 --mpheno-list |
| --repeatability-model | 重复力 GWAS | 配合 --tpheno |
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
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --out out.exact
```

### 多性状
```bash
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --scan-tier semi --no-loco --out mv
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 --multivariate --gwas --scan-tier exact --no-loco --out mv.exact
```

### 重复力
```bash
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --out rep
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 --repeatability-model --repeatability-mode fast --gwas --scan-tier approx --out rep.fast
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
| .scan.meta.txt | analysis_model, scan_tier, vc_strategy |
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

## 版本

v1.1.2。引用信息见 CITATION.cff。
