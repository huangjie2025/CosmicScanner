# 系谱动物模型 GWAS（Animal Model）说明

`--animal-model` 提供基于系谱亲缘关系矩阵 **A** 的加性遗传背景模型的 VCE 与全基因组扫描，是 `ssGWAS` 的"纯系谱"退化形式：不加基因组信息，用 A 描述亲缘结构。

## 统计模型

```
y = Xb + x_j * beta_j + Z*a + e
a ~ N(0, A * sigma2_a)
e ~ N(0, I * sigma2_e)
```

- `X`：固定效应设计（intercept + `--covar/--qcovar`）。
- `x_j`：第 j 个 SNP 的基因型列（剂量 0/1/2）。
- `Z`：样本 → 个体映射（恒等，逐样本表型）。
- `A`：系谱加性亲缘矩阵（自身=1+近交，个体间=2×共祖系数）。
- 方差分量 `(sigma2_a, sigma2_e)` 由 REML（`EigenVce`）估计，`h2 = sigma2_a / (sigma2_a + sigma2_e)`。

## 用法

```bash
cosmicscanner \
  --bfile my_data \
  --pheno my_pheno.txt \
  --animal-model \
  --ped pedigree.txt \
  --gwas \
  --out animal_assoc
```

系谱文件（空白分隔、无表头）：`ID SIRE DAM`，`SIRE`/`DAM` 未知可写 `NA`/`.`/`0`。`ID` 必须覆盖所有分析样本；可额外包含未测序祖先。

| 参数 | 说明 |
|------|------|
| `--animal-model` | 启用系谱动物模型 |
| `--ped <f>` | 系谱文件 `ID SIRE DAM`（必需） |
| `--write-ped-kernel` | 导出分析子集精确 A22 到 `*.animal_model.k.txt` + `*.animal_model.k.id`（可直接给 GEMMA `-k`） |
| `--no-loco` | 关闭 LOCO；多染色体数据请加上以消除 LOCO-not-applicable 提示 |

## 输出

- `*.assoc.txt`：标准 `CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`。
- `*.animal_model.null.txt`：`sigma2_a / sigma2_e / h2`，以及系谱解析 / Ainv / A22 / VCE / scan 的耗时摘要。
- `*.animal_model.k.txt` / `*.animal_model.k.id`（仅 `--write-ped-kernel`）：精确稠密 A22 与样本顺序（`FID:IID`）。

## 与 ssGWAS 的关系

两者共享同一 `PedigreeRelationship` 系谱解析与 `A⁻¹`/`A22` 构建（`buildA22` 对 A⁻¹ 做 Schur complement 求解，不物化稠密 A）：

- 单步模型使用 `H⁻¹ = A⁻¹ + (G⁻¹ − A22⁻¹)` 的 A22 校正。
- Animal Model 直接用 A22（即纯 `A` 核）。
- 数值等价性：ssGWAS 加 `--h2-tau 0 --h2-omega 0 --no-loco` 后退化为 animal model，逐 SNP 的 BETA/SE/-log10P 与直接 `--animal-model` 完全一致（Pearson=1.0）。

## 验证

- `tests/run_test_animal_model.sh`：7 段模拟套件，全部 PASS（causal top10 + LOCO 提示、ssGWAS tau=0 等价、系谱乱序不变性 <1e-12、σ²a=0 时 h2<0.1、高遗传力 h2>0.5、不完整系谱 causal top10、无祖先系谱 600 行有限值且收敛）。
- `tests/run_animal_model_gemma_validation.sh`：把 `--write-ped-kernel` 导出的 A22 交给 GEMMA 0.98.5 `-k -lmm 4`，比较 β/SE/-log10P 与 null VCE；无 GEMMA 时报 SKIP（exit 0）。实测：β Pearson≈0.99996、SE≈0.9992、-log10P(wald)≈0.9984，top5 完全重合、因果 SNP 排名 1/1；GEMMA null VCE `vg=0.500265, ve=0.862614` 与 CosmicScanner 一致（原 A22 无需归一化，刻度被 σ²g 吸收）。
- A⁻¹/A22 数值对照：对同一系谱文件，A⁻¹ 与 Python 暴力 inv(A) 的 maxabs 差=0.0；Schur 法 A22 与暴力 A22 的 maxabs 差≈5e-8。

## 边界

- LOCO 不适用：A 只有一个关系矩阵，无法按染色体分离；多染色体且未指定 `--no-loco` 时打印提示（输出警告，非错误）。
- 分析样本必须出现在系谱中（ID 按 `FID:IID` 再裸 `IID` 匹配）。
- 该模型在开发分支提供，正式报告建议与标准 LMM/GRM 结果并列解释。
