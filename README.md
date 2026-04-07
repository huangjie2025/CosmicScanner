# CosmicScanner

CosmicScanner 是 Cosmic 基因组分析套件的全基因组关联分析（GWAS）模块。它专为大规模生物样本库设计，支持高性能并行计算和现代基因组数据格式。

详细使用说明请参阅：**[CosmicScanner 使用手册 (中文)](docs/USAGE_cn.md)**

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

## 开发状态
- [x] PLINK 1.9/2.0 & BGEN 读取支持
- [x] 高效 GRM 构建 (Dense & Sparse)
- [x] Exact LMM & Fast-LMM 扫描引擎
- [x] LOCO 策略实现
- [x] 多性状 (Multivariate) LMM
- [x] 重复测量模型 (Repeatability Model)
