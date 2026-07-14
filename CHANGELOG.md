# CosmicScanner Source Package Changelog

## Version history note

- `v1.1.0` is the first public CosmicScanner focused release.
- `1.0.0` and `1.0.5` below are retrospective internal milestones inferred from
  Git history, not earlier public releases or Git tags.

## v1.1.0-source (2026-07-14)

- 首发源码包只保留单性状 GWAS、多性状 GWAS、重复力 GWAS 和 A/G/H 矩阵构建四项功能。
- 源码目录统一为单一扁平 `src/`，不再保留 `libCosmicCore/`、`src/core/` 或
  `src/models/` 子目录。
- 删除首发源码包无关模型和工具源码，包括 RRM、PSR/count/ordinal、ssGWAS
  关联扫描、maternal、survival、GxE、rare-variant、biobank/REGENIE、
  BSLMM/DAP 和 summary-stat 工具。
- 新增 focused CLI 和 focused CMake，仅构建四项功能所需入口、I/O、GRM、
  LMM/MvLMM、repeatability、A/G/H 和输出代码。
- 包内 release validation 覆盖单性状、多性状、重复力 semi/fast 和 A/G/H 构建。
- 重复力验证包含个体内中心化记录级协变量，覆盖 VCE 的零右端项求解回归。
- 默认 repeatability Woodbury 路径跳过无用的完整 GRM 求逆，并在 null summary
  报告实际分析受试者数。
- 新增文章范围、方法、输入输出、验证、复现、源码索引和引用元数据。

## 1.0.5 - retrospective internal milestone (2026-05-22 to 2026-06-24)

- Dedicated repeatability model with exact and fast analysis routes.
- Expanded multivariate compatibility and reference-validation paths.
- Consolidated single-trait, multivariate, and repeatability regression tests.

## 1.0.0 - retrospective internal milestone (2026-04-07)

- Independent GWAS baseline with single-trait LMM/LOCO and an initial
  multivariate LMM path.
- Early repeated-measure analysis used subject-level aggregation rather than
  the later exact record-level model.
