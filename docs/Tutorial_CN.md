# CosmicScanner 官方教程 (Tutorials)

CosmicScanner 是 Cosmic 生态系统中专为**大规模生物样本库 (Biobank)** 级别设计的高性能全基因组关联分析 (GWAS) 扫描引擎。它不仅在计算精度上对标经典软件 GEMMA，还在底层架构上实现了对现代高压缩基因组格式 (PLINK 2.0 PGEN、BGEN) 的原生支持，能够轻松应对百万级样本的关联分析。

本教程将详细介绍 CosmicScanner 的输入格式、数据过滤、基因组关系矩阵 (GRM) 构建、混合模型关联分析 (Exact LMM/Fast LMM) 以及多表型联合分析。

---

## 1. 基因组数据读入与过滤

CosmicScanner 的数据 I/O 模块经过高度优化，能够极速解析以下三种主流格式，无需任何中间格式转换：

1.  **PLINK 1.9 二进制格式 (`.bed`, `.bim`, `.fam`)**:
    使用 `--bfile <prefix>` 指定文件前缀。
2.  **PLINK 2.0 高压缩格式 (`.pgen`, `.pvar`, `.psam`)**:
    使用 `--pfile <prefix>` 指定文件前缀（推荐用于数十万样本级别的数据）。
3.  **BGEN 格式 (Layout 2 zlib压缩)**:
    使用 `--bgen <file>` 配合 `--sample <file>` 指定。

### 数据过滤 (Quality Control)
在分析之前，您可以直接使用 CosmicScanner 内置的质控参数对位点进行实时过滤：
*   `--maf <float>`：最小等位基因频率过滤（默认：0.01）。
*   `--hwe <float>`：哈迪-温伯格平衡 (Hardy-Weinberg Equilibrium) P 值过滤（默认：0.0，即禁用）。

```bash
# 示例：读取 BGEN 数据并过滤低频位点
./cosmic_scanner --bgen ukb_chr22.bgen \
                 --sample ukb_chr22.sample \
                 --maf 0.05 \
                 --hwe 1e-6
```

---

## 2. 基因组关系矩阵 (GRM) 构建

在进行混合线性模型 (LMM) 分析之前，需要构建基因组关系矩阵来校正群体分层和隐性亲缘关系。

CosmicScanner 提供了独立的 GRM 构建模式：
```bash
./cosmic_scanner --bfile demo \
                 --make-grm \
                 --make-grm-alg 1 \
                 --out demo_grm
```

| 算法 (`--make-grm-alg`) | 描述 |
| :--- | :--- |
| **0 (Yang)** | GCTA 默认的 GRM 构建算法。 |
| **1 (VanRaden)** | 默认加性模型编码，更适用于农本育种与常规关联分析。 |

> **进阶功能**：如果您正在分析混合群体或需要加权 GRM，您还可以提供 `--pop-class <file>` 或 `--snp-weight <file>` 来构建特定结构的亲缘关系矩阵。

---

## 3. 全基因组关联分析 (GWAS Scan)

CosmicScanner 提供了多种模式以适应不同的样本规模和精度需求。

### 标准精确混合模型 (Exact LMM)
精确模型会在每一个测试的 SNP 位点处，重新评估其固定效应，同时校正协变量和 GRM 的多基因背景效应。它与 GEMMA 的精度完全一致。

```bash
./cosmic_scanner --bfile demo \
                 --pheno phenotype.txt \
                 --mpheno 1 \
                 --covar discrete_covar.txt \
                 --qcovar quantitative_covar.txt \
                 --gwas \
                 --out my_exact_results
```
*   `--pheno`：表型文件（FID IID Pheno...）。
*   `--mpheno <n>`：指定表型文件中的列号。
*   `--covar` / `--qcovar`：离散（分类）协变量与定量（连续）协变量文件。

### 极度不平衡表型的鞍点近似校正 (SPA)
在病例对照 (Case-Control) 数据极度不平衡（如患病率仅有 1%）的情况下，传统的 LMM 会产生严重的 P 值假阳性膨胀。添加 `--spa` 参数，CosmicScanner 将自动使用鞍点近似 (Saddle Point Approximation) 对边际显著的位点进行极高精度的 P 值重校准，彻底解决假阳性问题 (对标 SAIGE)。

### 近似加速混合模型 (Fast LMM)
对于数十万样本的 Biobank 数据，Exact LMM 会极其耗时。此时推荐使用 `--fast-lmm`。该模式下，程序会先拟合一次仅包含协变量的零模型（Null Model），固定遗传方差与残差方差的比例 ($\delta = \sigma_e^2 / \sigma_g^2$)，然后在全基因组扫描时复用此方差比例，速度可提升数十倍。

### Biobank 级岭回归降维模式 (REGENIE-style)
当样本量达到数百万时，连计算 GRM 都会面临内存灾难。您可以开启 `--biobank-mode`。
该模式下，CosmicScanner 将完全采用类似 **REGENIE** 的分块交叉验证岭回归算法（Block-wise Ridge Regression），绕过 $N \times N$ 的亲缘关系矩阵，直接计算出每个人的多基因遗传背景并生成留一染色体 (LOCO) 预测值。

```bash
./cosmic_scanner --pfile ukb_data \
                 --pheno ukb_pheno.txt \
                 --gwas \
                 --biobank-mode \
                 --block-size 1000 \
                 --out my_biobank_results
```

### 多表型联合分析 (Multivariate GWAS)
在面临复杂疾病或高度相关的多个生理指标时，CosmicScanner 支持直接进行多表型联合 GWAS 扫描。
```bash
./cosmic_scanner --bfile demo \
                 --pheno multi_pheno.txt \
                 --mpheno-list 1,2,3 \
                 --multivariate \
                 --gwas \
                 --out my_mv_results
```

---

## 4. 罕见变异区域检验 (Rare Variant Association Test)

当变异的次等位基因频率 (MAF) 极低时，单点关联检验效力极差。CosmicScanner 支持基于区域 (例如基因) 的罕见变异聚合检验，全面对标 SAIGE-GENE 和 SKAT 包。

您需要提供一个 `.set` 格式的基因区域文件 (两列：`SetID` 和 `SNP_ID`)，并配合 `--skat` 或 `--burden` 参数使用。

```bash
./cosmic_scanner --bfile demo \
                 --pheno phenotype.txt \
                 --set-file genes.set \
                 --skat \
                 --out my_rare_variants
```
*   `--skat`：序列核关联检验 (Sequence Kernel Association Test)，假设同一基因内不同变异的效应方向有正有负。
*   `--burden`：Burden 负荷检验，假设同一基因内所有罕见变异具有相同的效应方向。

---

## 6. 与 CosmicSolver 的生态联动 (GS/GEBV & SBLUP)

作为 PanOmics 泛组学双引擎，CosmicScanner 输出的 GWAS 汇总统计量文件 (`.assoc.txt` 或 GWAS-SSF 格式) **原生无缝兼容** CosmicSolver 的高级预测与汇总级别分析模块！

*   **SBLUP 汇总级别分析**: 
    您可以直接将 CosmicScanner 输出的结果文件传给 CosmicSolver 的 `--sumstat` 参数，结合参考面板计算联合 SNP 效应：
    ```bash
    ./cosmicsolver --sblup --sumstat my_exact_results.assoc.txt --bfile reference_panel --h2 0.3 --out sblup_res
    ```
*   **GEBV 多基因风险评分**:
    如果您需要基于 CosmicScanner 跑出的 SNP 效应为新队列打分，同样可以直接使用：
    ```bash
    ./cosmicsolver --pred --score my_exact_results.assoc.txt --bfile new_cohort --out new_gebv
    ```
*(注：CosmicSolver 的 `--sumstat` 和 `--score` 引擎已具备智能表头识别能力，能自动解析 CosmicScanner 输出的 `SNP`, `A1`, `A2`, `BETA` 等关键列，无需用户手动提取列或转换格式。)*

默认情况下，CosmicScanner 会输出类似于 GEMMA/GCTA 的制表符分隔文件 (`.assoc.txt`)，包含：
`CHR`, `SNP`, `BP`, `A1`, `A2`, `FREQ`, `N`, `BETA`, `SE`, `P`, `LOG10P`。

**EBI GWAS-SSF 格式支持**：
为了方便您将研究成果直接提交至国际公开数据库 (如 NHGRI-EBI GWAS Catalog)，您可以添加 `--out-ssf` 参数，程序将自动输出符合 GWAS-SSF 官方标准的关联汇总统计文件。
