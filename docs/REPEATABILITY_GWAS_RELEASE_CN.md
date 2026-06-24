# CosmicScanner 重复力 GWAS 可运行发布说明

> 面向课题组内部科研使用。当前文档只总结 `--repeatability-model` 重复力 GWAS，不包括 PSR-GWAS、RRM 随机回归 GWAS 或 exact GLMM。

## 1. 当前结论

CosmicScanner 的重复力 GWAS 板块已经具备可投入课题组内部分析的基础条件：

- 支持同一个体多条记录的长格式表型。
- 支持记录层重复力 LMM：加性遗传效应 `g`、永久环境效应 `pe`、记录残差 `e`。
- 支持 `exact`、`fast`、`stcg` 三种路线，其中正式科研分析优先使用 `exact`。
- 支持 PLINK、PGEN、BGEN 基因型输入，PLINK 路径最成熟。
- 支持记录级连续协变量和分类协变量，带时间列时按 `FID IID TIME` 对齐。
- 输出标准 GWAS 结果、null model 方差组分摘要和 GRM cache 摘要。
- 已有 smoke test、edge cases、exact-vs-fast simulation、format consistency、performance benchmark 和 RepeatABEL 对标脚本。

推荐发布定位：

- `--repeatability-mode exact`：正式主分析路线。
- `--repeatability-mode fast`：快速预筛路线，不能替代正式结果。
- `--repeatability-mode stcg`：大规模矩阵自由原型，暂不作为课题组默认生产路线。

## 2. 模型说明

重复力 GWAS 适用于同一个体有多次重复测量的表型，例如多次采样、多个批次、多个时间点但不强调动态轨迹形状的性状。

模型为：

```text
y_ij = x_ij' b + m_i beta + g_i + pe_i + e_ij
g  ~ N(0, Vg G)
pe ~ N(0, Vpe I)
e  ~ N(0, Ve I)
```

其中：

- `y_ij`：个体 `i` 的第 `j` 条记录。
- `x_ij`：记录层固定效应或协变量。
- `m_i`：个体 `i` 的 SNP 剂量，同一个体多条记录共享同一基因型。
- `beta`：待检 SNP 效应。
- `G`：由基因型构建的个体 GRM。
- `Vg`：加性遗传方差。
- `Vpe`：永久环境方差。
- `Ve`：记录层残差方差。

输出中：

```text
h2 = Vg / (Vg + Vpe + Ve)
r2 = (Vg + Vpe) / (Vg + Vpe + Ve)
```

这里 `r2` 是重复力估计。

## 3. 输入格式

### 3.1 基因型

推荐首选 PLINK 1 binary：

```text
data/my_geno.bed
data/my_geno.bim
data/my_geno.fam
```

也支持：

- PLINK 2：`--pfile prefix`
- BGEN：`--bgen file.bgen --sample file.sample`

正式课题组使用建议先从 PLINK 路径跑通，再按需要测试 PGEN/BGEN。

### 3.2 长格式表型

不带显式时间列时：

```text
FID IID PHENO
F001 I001 12.3
F001 I001 13.1
F002 I002 10.8
F002 I002 11.4
```

带显式时间列时：

```text
FID IID TIME PHENO
F001 I001 1 12.3
F001 I001 2 13.1
F002 I002 1 10.8
F002 I002 2 11.4
```

对应参数：

- 不带时间列：`--mpheno 1`
- 带时间列：`--tpheno 1 --mpheno 2`

缺失值使用 `NA` 或 `-9`。基因型中不存在的样本会被自动丢弃。

### 3.3 协变量

如果表型带 `TIME`，记录级协变量也应使用：

```text
FID IID TIME COV1 COV2 ...
```

参数：

- 连续协变量：`--qcovar file`
- 分类协变量：`--covar file`

如果表型不带 `TIME`，协变量按个体层 `FID IID` 对齐。

## 4. 推荐命令

### 4.1 正式主分析

```bash
./build/cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat.pheno \
  --tpheno 1 \
  --mpheno 2 \
  --qcovar data/record.qcovar \
  --repeatability-model \
  --repeatability-mode exact \
  --repeatability-grm-cache results/my_trait_grm_cache \
  --gwas \
  --out results/my_trait_repeatability
```

如果没有时间列，去掉 `--tpheno 1`，并将 `--mpheno` 改为表型列号，例如：

```bash
./build/cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat.pheno \
  --mpheno 1 \
  --repeatability-model \
  --repeatability-mode exact \
  --repeatability-grm-cache results/my_trait_grm_cache \
  --gwas \
  --out results/my_trait_repeatability
```

### 4.2 快速预筛

```bash
./build/cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat.pheno \
  --tpheno 1 \
  --mpheno 2 \
  --repeatability-model \
  --repeatability-mode fast \
  --gwas \
  --out results/my_trait_repeatability_fast
```

`fast` 会把多条记录聚合到个体均值层面，是近似路线。正式报告建议以 `exact` 为准。

### 4.3 小数据 refit 验证

小样本或候选 SNP 验证时可打开 per-SNP refit：

```bash
./build/cosmicscanner \
  --bfile data/my_geno \
  --pheno data/repeat.pheno \
  --tpheno 1 \
  --mpheno 2 \
  --repeatability-model \
  --repeatability-mode exact \
  --repeatability-refit-validation \
  --repeatability-refit-limit 100 \
  --gwas \
  --out results/my_trait_repeatability_refit_check
```

## 5. 输出文件

主结果：

```text
<out>.assoc.txt
```

常规列：

```text
CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P
```

Null model 摘要：

```text
<out>.repeatability.null.txt
```

重点检查字段：

- `n_records`
- `n_subjects`
- `Vg`
- `Vpe`
- `Ve`
- `h2`
- `r2`
- `scan_mode`

GRM cache 摘要：

```text
<out>.repeatability.cache.txt
```

如果启用 refit validation，还会输出：

```text
<out>.refit.assoc.txt
```

## 6. 发布前验证

在 `/disk228/huangjie/PanOmics/CosmicScanner` 中执行：

```bash
make -C build cosmicscanner -j4
bash tests/run_repeatability_release_validation.sh
```

快速验证脚本会运行：

- `tests/run_test_repeatability.sh`：exact、LOCO、cache、fast smoke。
- `tests/run_test_repeatability_edge_cases.sh`：边界输入检查。
- `tests/run_repeatability_exact_vs_fast_simulation.sh`：小规模 exact/fast/mean simulation。

如果本机有 PLINK2，还可以额外跑：

```bash
bash tests/run_repeatability_format_consistency.sh
```

该脚本验证 PLINK/PGEN/BGEN 输出一致性。

## 7. 课题组使用建议

正式交付给同门时建议按以下流程：

1. 先准备一个性状的小数据子集，确认 ID、TIME、协变量和重复记录格式无误。
2. 使用 `--repeatability-mode exact --no-loco` 跑一次快速检查。
3. 检查 `<out>.repeatability.null.txt` 中 `n_records`、`n_subjects`、`Vg`、`Vpe`、`Ve`、`h2`、`r2` 是否合理。
4. 正式全量分析去掉 `--no-loco`，使用默认 LOCO。
5. 每个性状单独设置 `--out` 前缀，避免结果覆盖。
6. 大量性状重复运行时打开 `--repeatability-grm-cache`。
7. top SNP 建议用 `--repeatability-refit-validation` 或外部工具做候选验证。

## 8. 已知限制

- 当前重复力 GWAS 是连续性状 LMM，不适合直接分析二分类、计数、有序分类性状。
- `fast` 是均值聚合近似，不建议作为最终论文主结果。
- `stcg` 是大规模扩展原型，暂不作为默认生产路线。
- RRM 随机回归 GWAS 是另一个模块，不能把 RRM 结果直接当作重复力 GWAS 汇报。
- 极端不平衡重复次数、严重缺失、协变量完全共线时，需要先清洗数据并检查日志警告。
- 当前文档面向内部可运行发布，不等同于软件正式公开 release。

## 9. 相关材料

- 方法学中文稿：`docs/Repeatability_GWAS_Methods_论文稿_CN.md`
- 全景开发总结：`docs/Repeatability_GWAS_全景开发总结_20260529.md`
- Methods/Results 草稿：`docs/Repeatability_GWAS_Manuscript_Methods_Results_CN.md`
- Smoke test：`tests/run_test_repeatability.sh`
- Edge cases：`tests/run_test_repeatability_edge_cases.sh`
- Exact-vs-fast simulation：`tests/run_repeatability_exact_vs_fast_simulation.sh`
- Format consistency：`tests/run_repeatability_format_consistency.sh`
- Performance benchmark：`tests/run_repeatability_performance_benchmark.sh`
