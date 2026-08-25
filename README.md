# CosmicScanner v1.0.0

CosmicScanner 是 PanOmics 基因组分析套件的全基因组关联分析（GWAS）模块。它专为大规模生物样本库设计，支持高性能并行计算和现代基因组数据格式。

详细使用说明请参阅：**[中文教程](docs/Tutorial_CN.md)** | **[English Tutorial](docs/Tutorial_EN.md)**

---

## 版本与发布状态

**当前版本**: v1.0.0 (PanOmics v1.0.0 组件)

**正式支持模型** (v1.0 范围):

| 模型 | CLI 标志 | 状态 | 说明 |
|------|----------|------|------|
| 单性状 LMM GWAS | `--gwas` | **Stable** | 标准高斯 LMM，支持 LOCO、Fast-LMM 近似 |
| 多性状 mvLMM GWAS | `--multivariate --gwas` | **Stable** | 多表型联合分析，支持 EMmv/GEMMA/hybrid 方法 |
| 重复力 GWAS | `--repeatability-model` | **Experimental validated** | 记录层精确模型，支持 exact/fast/STCG 扫描 |
| 随机回归模型 GWAS | `--rrm-order` | **Experimental validated** | Method A (record-level) + Method B (coefficient) |
| 系谱动物模型 GWAS | `--animal-model --ped` | **Experimental validated** | 纯系谱 A 矩阵，LOCO 不适用 |
| 母体效应模型 | `--maternal-model` | **Experimental validated** | direct/maternal/joint 检验 |
| 共同环境模型 | `--common-env-model` | **Experimental validated** | 基因组 GRM + 分组共享环境随机截距 |
| ssGWAS (H-matrix) | `--single-step --single-step-method lmm` | **Experimental validated** | Method 1: 固定 SNP 关联分析 |
| ssGBLUP back-solving | `--single-step-method backsolve` | **Experimental validated** | Method 2: EBV/GEBV → SNP 效应反解 |
| APY 关系矩阵 | `--apy` | **Experimental validated** | 核心/非核心分块求逆 |

**实验性模型** (默认关闭，需 `-DCOSMICSCANNER_ENABLE_EXPERIMENTAL=ON`):

| 模型 | CLI 标志 | 状态 | 说明 |
|------|----------|------|------|
| Ordinal PSR-GWAS | `--ordinal-model` | Experimental | 累积二元/堆叠 PSR |
| Threshold/Probit | `--threshold-model` | Experimental | 二元性状 PQL + probit/logit |
| Count/Binomial PSR | `--count-model` | Experimental | Poisson/NegBin/Binomial |
| Survival/Cox | `--survival-model` | Experimental | Cox PH + martingale residual |
| GxE Reaction Norm | `--gxe-model` | Experimental | SNP×环境交互 |
| BSLMM 兼容 | `--bslmm` | Experimental | GEMMA 风格多位点模型 |
| 稀有变异测试 | `--skat` / `--burden` | Experimental | SKAT/Burden 集合检验 |
| Biobank REGENIE | `--biobank-mode` | Experimental | UKB 规模两步回归 |

> **注意**: 实验性模型源码保留在 `src/experimental_models/`，不参与默认构建，不保证生产级稳定性。

---

## 安装

### 预编译版本（推荐）
```bash
tar -zxvf CosmicScanner-v1.0.0-linux-x86_64.tar.gz
cd CosmicScanner-v1.0.0-linux-x86_64
bash install.sh
bin/cosmicscanner --help
```

### 源码编译
```bash
# 默认构建 (仅正式模型)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# 启用实验性模型
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOSMICSCANNER_ENABLE_EXPERIMENTAL=ON
cmake --build build -j$(nproc)
```

编译成功后，可执行文件位于 `build/cosmicscanner`。

### 依赖
- CMake >= 3.15
- 支持 C++17 的编译器 (GCC 9+ / Clang 10+)
- OpenMP（编译时自动检测）
- Eigen3, ZLIB, zstd (通过 CMake 自动查找)
- 可选：外部 BLAS/LAPACK（默认关闭，以 Eigen 原生求解器保证跨平台数值复现性）

---

## 快速使用示例

```bash
# 1. 计算 GRM
cosmicscanner --bfile my_data --make-grm --out my_grm

# 2. 运行单性状 GWAS（默认 LOCO，推荐 semi）
cosmicscanner --bfile my_data --pheno my_pheno.txt --mpheno 1 --gwas --scan-tier semi --out my_results

# 3. 运行多性状 GWAS
cosmicscanner --pfile my_ukb_data --pheno my_pheno.txt --mpheno-list 1,2,3 --multivariate --gwas --scan-tier semi --out my_mv_results

# 4. 重复力 GWAS
cosmicscanner --bfile my_data --pheno my_pheno.txt --repeatability-model --mpheno 1 --gwas --out my_rep

# 5. ssGWAS
cosmicscanner --bfile my_data --ped pedigree.txt --single-step --gwas --out my_ssgwas
```

---

## Gaussian LMM 扫描档位

单性状、多性状和重复力 Gaussian LMM 支持统一的 `--scan-tier approx|semi|exact`。三档分别对应 fixed-null、fixed-ratio/profile-scale 和 per-SNP full variance-component refit。详细命令及当前边界见 [docs/SCAN_TIERS_CN.md](docs/SCAN_TIERS_CN.md) / [docs/SCAN_TIERS_EN.md](docs/SCAN_TIERS_EN.md)。

---

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

---

## 文档

| 文档 | 说明 |
|------|------|
| `docs/Tutorial_EN.md` / `docs/Tutorial_CN.md` | 完整使用教程 |
| `docs/PAPER_MODEL_SUPPORT_MATRIX.md` | 论文模型支持矩阵 |
| `docs/FINAL_VALIDATION_REPORT.md` | 最终验证报告 |
| `docs/COSMICGWAS_V1.1.0_RELEASE_NOTES.md` | 发布说明 |
| `docs/GWAS_MODEL_DEFINITIONS.md` | 模型数学定义 |
| `docs/SCAN_TIERS_EN.md` / `docs/SCAN_TIERS_CN.md` | 扫描档位说明 |

---

## 示例数据

`examples/mouse_clean/` 包含标准 PLINK 格式示例数据（`.bed/.bim/.fam/.pheno`），可直接用于快速上手测试。

---

## 许可证与贡献

- 许可证：见 `LICENSE`
- 贡献指南：见 `CONTRIBUTING.md`
- 变更记录：见 `CHANGELOG.md`