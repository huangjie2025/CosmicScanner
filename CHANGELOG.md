# CosmicScanner Source Package Changelog

## Version history note

- `v1.1.0` is the first public CosmicScanner focused release.
- `1.0.0` and `1.0.5` below are retrospective internal milestones inferred from
  Git history, not earlier public releases or Git tags.

## v1.1.2 (2026-07-15)

- 增加 Linux、macOS 和 Windows UCRT64 的 GitHub Actions 构建与回归测试。
- 每次 CI 构建生成平台压缩包；版本 tag 自动创建 GitHub Release 并上传资产。
- 修复 Windows 资源统计、只读内存映射、路径和可选 OpenBLAS 符号兼容问题。
- 修复 macOS peak RSS 单位，并显式配置 Apple Clang + Homebrew libomp。
- Windows 压缩包携带运行时 DLL；macOS 压缩包携带并重写非系统 dylib。

## v1.1.1 (2026-07-15)

- 修复 PLINK GRM 并行构建中 `std::vector<bool>` 位压缩写入导致的数据竞争、
  有效 SNP 计数不一致和未初始化矩阵列问题。
- 修复 `compute_homo_hete` 在 SNP 并行循环中更新样本计数器的数据竞争。
- 重复力 GWAS 的 exact/approx/semi/fast/STCG 路径现在统一遵守
  `--chr`、`--from-bp`、`--to-bp`、`--extract` 和 `--exclude`。
- 增加 1/2/4/8 线程 GRM 确定性测试、重复力区域过滤回归测试、CTest 和
  GitHub Actions CI，并加入 PLINK/PGEN/BGEN 格式等价性测试数据。
- 修复 PGEN 1-bit + difflist 解码将 missing code 当作剂量值的问题。
- 增加 `--version`、基因型输入互斥校验，并将 BGEN sample 数量不匹配改为错误。

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
