# Gaussian 混合模型 GWAS 三档扫描

`--scan-tier` 统一控制单性状 LMM、多性状 LMM 和重复力 LMM 在逐 SNP 扫描时如何处理方差参数：

| 档位 | 方差处理 | 定位 |
|---|---|---|
| `approx` | 固定 null model 的全部方差参数 | 最高吞吐量 |
| `semi` | 固定方差比例或协方差形状，每个 SNP profile 残差尺度 | 推荐全基因组主扫描 |
| `exact` | 每个 SNP 重估 alternative 方差参数 | 小中型数据或候选位点精细分析 |

这里的档位描述 nuisance variance refit 程度，不代表 fixed-null Score 检验无效，也不代表有限样本 P 值“绝对精确”。重复力的 `--repeatability-mode exact|fast` 与扫描档位在概念上是不同维度；当前 `--scan-tier` 实现仅接入 `--repeatability-mode exact`。

## 单性状

```bash
# fixed null
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier approx --out out.approx

# fixed variance ratio + SNP-specific residual scale
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier semi --out out.semi

# per-SNP variance-ratio and scale refit
cosmicscanner --bfile data --pheno trait.txt --gwas --scan-tier exact --out out.exact
```

`--fast-lmm` 保留为旧 fixed-null 命令兼容入口。新分析应使用 `--scan-tier approx`。

## 多性状

```bash
cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier approx --out mv.approx

cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier semi --out mv.semi

cosmicscanner --bfile data --pheno traits.txt --mpheno-list 1,2 \
  --multivariate --gwas --scan-tier exact --no-loco --out mv.exact
```

映射关系为 `approx -> fixed-null Wald`、`semi -> EMMV`、`exact -> GEMMA exact-style full alternative refit`。严格复现 EMmvGWAS 的半精准公式时增加 `--mv-compatible emmvgwas --no-loco`。

## 重复力

```bash
cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas --scan-tier approx --out rep.approx

cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas --scan-tier semi --out rep.semi

cosmicscanner --bfile data --pheno repeated.txt --tpheno 1 --mpheno 2 \
  --repeatability-model --repeatability-mode exact --gwas --scan-tier exact --no-loco --out rep.exact
```

重复力 `semi` 固定 `Vg/Ve` 和 `Vpe/Ve`，通过 Schur 更新对每个 SNP 闭式 profile `Ve`，不重复进行 Woodbury 分解。当前重复力 `exact` 生产入口限 PLINK、no-LOCO、无记录权重；不支持的组合会明确报错。

当前重复力 LOCO 的 `approx/semi` 路线使用 global null ratios 配合 LOCO GRM，并在 null summary 中标记 `null_vc_scope=global_vc_with_loco_grm`。需要逐染色体重估 null ratios 的严格 LOCO 仍应作为后续增强，不能与当前 global-ratio LOCO 混写。

三类分析都会输出 `<out>.scan.meta.txt`，重复力还会在 `<out>.repeatability.null.txt` 中记录 `scan_tier` 和 `snp_vc_strategy`。不指定 `--scan-tier` 时保留旧命令行为。
