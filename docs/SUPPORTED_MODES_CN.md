# CosmicScanner 支持模式与核矩阵总表

本文档按"模型入口 × 随机效应核"整理 CosmicScanner 的 GWAS 模式覆盖，方便选择分析路线与外部工具对标。

## 模型入口总表

| 模型 | 入口 | 随机效应核 | 类型 | 方差分量 | LOCO | 输出 |
|------|------|-----------|------|---------|------|------|
| 单性状 LMM | `--gwas` | 基因组 GRM **G**（VanRaden/Zeng/Yang/Vitezica/TrueVar） | 标准 | σ²g, σ²e | 支持 | `*.assoc.txt` + null |
| 多性状 mvLMM | `--multivariate` | **G**（EMmvGWAS-compatible/GEMMA 路线） | 标准 | Σg, Σe | 支持（emmv/hybrid） | `*.assoc.mv.txt` 系列 |
| 重复力 GWAS | `--repeatability-model` | 个体均值 GRM + 记录层 I | 标准 | σ²g, σ²pe, σ²e | 支持 | `*.assoc.txt` + null |
| RRM 动态 GWAS | `--long-pheno --rrm-order` | Legendre/B-spline 核 × GRM | 标准/近似 | 逐阶 σ² | 支持 | `*.rrm_joint.assoc.txt` 等 |
| 母体效应模型 | `--maternal-model` | GRM（母体配对） | 标准 | Va, Vm, Cov_am, Vpe, Ve | — | `*.assoc.txt` |
| 单步 ssGWAS | `--single-step` | **H⁻¹ = A⁻¹ + (G⁻¹ − A22⁻¹)**（A22 校正） | 标准 | σ²a, σ²e | 支持 | `*.assoc.txt` + null |
| 系谱动物模型 | `--animal-model` | 系谱 **A**（A22 子块） | 标准 | σ²a, σ²e | 不适用（提示） | `*.assoc.txt` + null + 可选 A22 核 |
| 公共环境模型 | `--common-env-model` | GRM（遗传）+ I（组共享） | 实验性 | σ²a, σ²c, σ²e | 支持 | `*.assoc.txt` + null |
| 阈值/Probit | `--threshold-model` | GRM（working LMM） | 近似 | σ²a, σ²e | 支持 | `*.assoc.txt` |
| 计数/二项 PSR | `--count-model` | GRM（weighted LMM） | 近似 | σ²a, σ²e | 支持 | `*.assoc.txt` |
| 序数 PSR | `--ordinal-model` | GRM（stacked/binary） | 近似 | σ²a, σ²e | 支持 | `*.gwas.tsv` |
| 生存 Cox | `--survival-model` | GRM（鞅残差 LMM） | 近似 | σ²a, σ²e | 支持 | `*.assoc.txt` |
| GxE 反应范数 | `--gxe-model` | GRM × 环境（Vg0+Vg1） | 标准 | σ²g0, σ²g1, σ²e | 支持 | `*.gxe.assoc.txt` |
| Biobank 两步法 | `--biobank-mode` | 分块 ridge LOCO 预测 | 近似 | 正则化 | 内建 LOCO | `*.biobank.assoc.txt` |

## 核矩阵（G / A / H）与工具

| 核 | 构建 | 输出 | 用途 |
|----|------|------|------|
| A⁻¹ | `--make-agh`（或 ssGWAS `--write-agh`） | `*.Ainv.bin/.id` + `*.agh.meta.txt` | ssGWAS H⁻¹ 组件 |
| G | `--make-grm`（alg 1-6） | `*.grm.bin/.id` | 单性状/多性状/重复力/RRM 等 |
| H⁻¹ | `--make-agh`（或 ssGWAS `--write-agh`） | `*.Hinv.bin/.id` | ssGWAS 使用 |
| A22 | `--animal-model --write-ped-kernel` | `*.animal_model.k.txt` + `*.k.id` | Animal Model 核；可直接给 GEMMA `-k` |

> 外部对标提示：GEMMA `-gk 1` 输出可被本项目 GRM 路线复用；本项目 A22 全矩阵导出可被 GEMMA `-k -lmm 4` 直接读取（刻度被 σ²g 吸收，无需归一化）。

## 类型说明

- **标准**：直接基于对应随机效应模型的 REML 求解，结果可作为正式结论。
- **近似/接口**：通过 working/残差替换或固定方差比实现，建议按对应方法边界解释结果并给出校准说明。
- **实验性**：开发分支实现，暂行限于开发树，不进入 CosmicGWAS v1.1.0 focused 发布范围。
