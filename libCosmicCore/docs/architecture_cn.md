# Cosmic/Common 模块架构文档

CosmicCommon 是 Cosmic 项目（包含 CosmicSolver 和 CosmicScanner）的共享核心库，提供了矩阵运算、混合模型方程（MME）构建、方程组求解、方差组分估计（VCE）以及多组学数据解析等基础功能。

## 1. 核心数学模块 (Math Core)

### 1.1 共轭梯度法求解器 (PCG Solvers)
- **文件**: `pcg_solver.h`, `matrix_free_pcg.h`
- **功能**: 求解大规模稀疏线性方程组 $Ax = b$。
- **特点**:
  - `pcg_solver.h`: 基于模板的高效预处理 PCG 求解器，支持多种预处理策略（如对角预处理 `diag`、块对角等），适用于舒尔补等显式矩阵。
  - `matrix_free_pcg.h`: 无矩阵（Matrix-Free）PCG 求解器，通过传入闭包函数（std::function）进行矩阵向量乘积，专为隐式 MME 优化，避免显式构建巨大矩阵。

### 1.2 Schur 补 (Schur Complement)
- **文件**: `schur.h`, `schur.cpp`
- **功能**: 计算 Schur 补矩阵 $S = Z'Z - Z'X(X'X)^{-1}X'Z$。
- **用途**: 在方差组分估计中用于降维，将混合模型转换为关于随机效应的简化模型。

### 1.3 内存映射与统一矩阵接口 (Matrix Data Structures)
- **文件**: `mmap_matrix.h`, `matrix_adapter.h`
- **功能**: 提供基于内存映射（Memory Mapping）的文件读取接口，以及统一稠密、稀疏和内存映射矩阵的抽象适配器。
- **用途**: 允许求解器在不关心底层存储格式的情况下，进行统一的高效矩阵乘法和访问操作，避免 OOM (Out-Of-Memory)。

## 2. 混合模型构建 (MME Builder)

### 2.1 隐式 MME (Implicit MME)
- **文件**: `implicit_mme.h`
- **功能**: 实现无矩阵（Matrix-Free）的 MME 乘法操作。
- **原理**: 结合 `matrix_free_pcg.h`，不显式构建巨大的 LHS 矩阵（$C$），而是动态计算 $C v = X'X v + Z'Z v + \lambda A^{-1} v$。
- **优势**: 极大降低内存消耗，适合超大规模生物样本库数据。

### 2.2 显式 MME 构建与设计矩阵 (Explicit MME & Design)
- **文件**: `mme_builder.h`, `mme_builder.cpp`, `design.h`, `design.cpp`
- **功能**: 管理固定效应（$X$）和随机效应（$Z$）结构，并显式构建稀疏 LHS 矩阵。
- **特点**: 支持快速更新 $\lambda$（方差比），无需重新扫描原始数据，适合 AI-REML 等需要频繁构建和求逆的算法。

## 3. 方差组分估计与全基因组扫描 (VCE & LMM Engines)

### 3.1 迭代 VCE 求解器 (Iterative VCE Solvers)
- **文件**: `vce.h`, `vce.cpp`
- **功能**: 针对复杂 LMM 的方差组分估计核心迭代算法。
- **算法**: 包括 AI-REML (Average Information REML)、EM (Expectation Maximization) 以及用于超大规模稀疏矩阵的 MC-AI/MC-EM (基于蒙特卡洛随机迹估计)。

### 3.2 精确 LMM 引擎 (Exact LMM & VCE)
- **文件**: `eigen_vce.h`, `eigen_vce.cpp`, `lmm_engine.h`, `mvlmm_engine.h`
- **功能**: 基于特征值分解（Eigen Decomposition）的精确方差组分估计和全基因组扫描引擎。
- **特点**:
  - `eigen_vce`: 封装了基于 Eigen 库的精确 VCE，服务于 GWAS 的 Null Model 计算以及小规模 GS 的方差计算。
  - `lmm_engine`: 单变量全基因组扫描核心，支持 Brent 优化和 Fast-LMM 模式。
  - `mvlmm_engine`: 多变量（Multi-trait）扫描引擎，支持 LOCO（留一染色体交叉）优化和块矩阵乘法加速。

## 4. 基因组数据解析与预处理 (Genomics & Parsers)

### 4.1 多格式基因型读取器 (Genotype Readers)
- **文件**: `plink_reader.h`, `pgen_reader.h`, `bgen_reader.h`, `sparse_grm_reader.h`
- **功能**: 支持业界主流遗传数据的极速解析。
- **特点**:
  - `plink_reader`: 标准 PLINK 1.9 `.bed` 格式内存映射读取。
  - `pgen_reader`: PLINK 2.0 PGEN 格式解析（支持 zstd 解压）。
  - `bgen_reader`: 英国生物样本库 BGEN 格式支持。
  - `sparse_grm_reader`: 稀疏 GRM 矩阵的高效导入。

### 4.2 数据对齐与转换 (Harmonization & Conversion)
- **文件**: `allele_harmonizer.h`, `genotype_converter.h`, `pheno_reader.h`, `covar_reader.h`
- **功能**: 负责表型、协变量的读取，以及等位基因的跨队列对齐（如链翻转处理）与文本基因型向数值剂量的转换。

### 4.3 GRM 与遗传矩阵构建 (Kinship Matrices)
- **文件**: `grm_builder.h`, `a_matrix.h`, `h_matrix.h`, `apy.h`
- **功能**: 提供 G矩阵 (VanRaden)、A矩阵 (系谱关系) 和 H矩阵 (一步法 ssGBLUP) 的构建。支持 APY (Algorithm for Proven and Young) 的核心/非核心群体分块求逆优化。

## 5. 辅助工具 (Utilities)

- **文件**: `rng.h`, `logger.h`, `qc_utils.h`
- **功能**: 封装随机数生成器、日志格式化输出和基因组数据质控（如 MAF 过滤、缺失率处理）基础函数。
