# CosmicScanner v2.0.0 命令行指南

安装和完整操作说明见中文用户手册。请按所选模型核对输入和输出字段。

## 1. 推荐命令结构

保留 PLINK/GCTA 用户熟悉的长参数：输入 → 表型/协变量 → 任务/模型 → 计算设置 → 输出。

```bash
```

也可以使用等价的可选任务入口：

```bash
cosmicscanner gwas --model lmm --bfile data --pheno trait.txt --threads 4 --out results/run
cosmicscanner grm --bfile data --threads 4 --out results/grm
```

`gwas` 只补充 `--gwas`；`grm` 只补充 `--make-grm`。已有平铺命令继续支持。
推荐新关联分析命令统一使用 `gwas --model`，不要依赖部分旧模型参数隐式启动分析。
旧 `--model lmm/fixed-lmm/lmm-eb/eb-lmm` 在平铺脚本中的历史 EB 开关作用保留；旧脚本不用立即迁移。新任务入口会检查互相冲突的模型选择。

## 2. 安装和帮助

以下示例假设你已获得包含此更新的源码和匹配 Core 子模块；不假设某个公开下载地址已经可用。
需要 C++17、CMake ≥3.15、Eigen、OpenMP、zlib、zstd，以及构建选择的 BLAS/LAPACK 运行库。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target cosmicscanner -j2
BIN=./build/cosmicscanner
"$BIN" --version
"$BIN" --help
"$BIN" --help models
"$BIN" help input
"$BIN" --help repeatability
"$BIN" --help-all
```

默认帮助是简明入口；`--help-all` 或 `--help all` 保留完整旧参数和高级/内部选项。
分组包括 gwas、models、input、grm、repeatability、rrm、ssgwas；all/advanced 是完整参考。
当前在项目 Linux 服务器构建验证；不能据此承诺其他系统或干净发行包安装已经验收。

## 3. 通用参数和列号

| 参数 | 含义 |
|---|---|
| `--bfile PREFIX` | BED/BIM/FAM 前缀 |
| `--pfile PREFIX` | PGEN/PVAR/PSAM 前缀 |
| `--bgen FILE --sample FILE` | BGEN 与样本文件，支持范围依模型而异 |
| `--pheno FILE` | FID IID 后接表型字段 |
| `--mpheno N` | **去掉 FID/IID 后**从 1 开始的表型列号 |
| `--mpheno-list 1,2` | 多性状列号，逗号分隔，使用相同计数规则 |
| `--tpheno N` | 时间列号，同样从 FID/IID 后开始计数 |
| `--covar FILE` / `--qcovar FILE` | 分类 / 连续协变量 |
| `--maf F --hwe P` | SNP 过滤；不是模型选择 |
| `--extract FILE --exclude FILE` | 纳入/排除的 SNP ID 文件 |
| `--threads N --out PREFIX` | 线程数和输出前缀；提前创建输出目录 |
| `--pedigree FILE` | ID SIRE DAM 系谱；限 animal、ssGWAS、AGH 场景 |
| `--ped-genotypes FILE --map FILE` | PED/MAP 基因型输入，与系谱不同 |

`FID IID VALUE` 使用 `--mpheno 1`；`FID IID TIME VALUE` 使用 `--tpheno 1 --mpheno 2`。
这些不是 PLINK2 `--pheno-col-nums` 的绝对文件列号。旧 `--ped` 继续支持，但含义依任务变化，推荐明确别名。
`--time-col` 属于原始表格适配流程，不是 `--tpheno` 的新名称。

## 4. 模型与计算策略

| 统一入口 | 对应旧入口 | 范围说明 |
|---|---|---|
| `--model lmm` | 默认单性状 LMM | 限定参考路径 |
| `--model lm` | `--lm` | 线性回归 |
| `--model multivariate` | `--multivariate` | 必须提供 `--mpheno-list`；参考 BED，其他格式依范围契约 |
| `--model repeatability` | `--repeatability-model` | 重复记录；明确记录策略和扫描档位 |
| `--model rrm` | RRM 选择与原参数 | 实验；平衡宽表与非平衡长表仍需各自参数 |
| `--model animal/maternal/common-env` | 相应 `--*-model` | 实验，保持原模型定义 |
| `--model threshold/count/ordinal/gxe` | 相应 `--*-model` | 实验，不扩大输入或推断承诺 |
| `--model ssgwas` | `--single-step` | 默认 Method 1；Method 2 需明确 backsolve |

`--model` 是模型家族；`--gwas-model additive|additive_dominance|dominance` 是单性状模型内部的效应结构，两者不同。
`--scan-tier approx|semi|exact` 控制逐 SNP 方差处理，不选择模型；exact 不代表有限样本精确 P 值。
原生多性状 semi 尚未实现，显式请求会失败；`--mv-compatible emmvgwas --scan-tier semi` 是独立兼容口径，不能无条件替代原生分析。
多性状 exact 和重复力 exact 扫描档位要求 no-LOCO；重复力该档位还限 BED、无记录权重。
`--repeatability-mode exact` 选择记录级计算，`--scan-tier exact` 选择逐 SNP 全重估，两者不是同一个开关。
不提供扫描档位仍保留原行为；旧原生 emmv 的 fixed-null 计算如实记录为 approx。
Survival/内部校准选项不纳入新的推荐模型入口，完整帮助中的可见选项不等于正式支持。

## 5. 从源码目录运行的完整合成示例

示例数据随当前源码位于 `examples/data/r33_synthetic/dataset`：60 个样本、180 个 SNP、两条染色体。
下面的 `$BIN` 和 `$DATA` 都是 Bash 变量，每次新终端需重新设置。

```bash
BIN=./build/cosmicscanner
DATA=examples/data/r33_synthetic/dataset
mkdir -p demo_results

# 单性状 LMM，默认 LOCO
"$BIN" gwas --model lmm --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" \
  --mpheno 1 --threads 1 --out demo_results/lmm

# 三性状 native fixed-null
"$BIN" gwas --model multivariate --bfile "$DATA/synth" \
  --pheno "$DATA/multivariate.pheno" --mpheno-list 1,2,3 \
  --scan-tier approx --reml-maxit 1200 --reml-tol 1e-10 --no-loco \
  --threads 1 --out demo_results/mv

# 记录级重复力；仅为可复现合成示例固定方差
"$BIN" gwas --model repeatability --bfile "$DATA/synth" \
  --pheno "$DATA/long.pheno" --tpheno 1 --mpheno 2 \
  --repeatability-mode exact --scan-tier approx --no-loco \
  --reml-vg 0.5 --reml-vpe 0.2 --reml-ve 0.3 --threads 1 --out demo_results/rep

# 关系矩阵
"$BIN" grm --bfile "$DATA/synth" --threads 1 --out demo_results/grm
```

固定方差 0.5/0.2/0.3 是演示设定，不是实际研究的通用估计值。真实数据可省略三项由模型估计，并检查拟合诊断。

| 输出 | 核对内容 |
|---|---|
| `lmm.assoc.txt` | 180 个 SNP，chr1/chr2；当前统一表使用 BETA_A/SE_A/P_A 等命名 |
| `lmm.cosmic.sumstats` | 180 个相同 SNP，保留表头和两条染色体 |
| `lmm.scan.meta.txt`、`lmm.run_manifest.json` | 实际分析配置、样本与原始命令 |
| `mv.assoc.mv.txt` | 180 行、联合检验与各性状效应/标准误 |
| `rep.assoc.txt`、`rep.repeatability.null.txt` | 180 行与记录模型诊断；与 LMM 不保证相同 schema |
| `grm.grm.bin`、`grm.grm.id` | 矩阵与匹配样本 ID |

```bash
head -2 demo_results/lmm.assoc.txt
head -7 demo_results/lmm.cosmic.sumstats
cat demo_results/lmm.scan.meta.txt
```

不要仅凭文件存在判断分析成功。检查退出状态、记录数、SNP 数、收敛诊断和非有限数值。

## 6. 常见错误与迁移

- 任务冲突：`grm --model lm`、`gwas --make-grm` 会拒绝；矩阵与关联分别运行。
- 模型冲突：新任务中不要同时写 `--model lmm --multivariate`；使用一个统一选择。
- 缺表型列：`--model multivariate` 必须带 `--mpheno-list`。
- 未收敛：多性状拟合失败会返回非零；检查数据/模型，必要时增加 `--reml-maxit`，不要忽略错误。
- 格式不支持：读取器存在不意味着所有模型支持 VCF/PGEN/BGEN；查看范围契约。
- 相对路径错误：从源码根目录运行上述示例；输出目录先用 `mkdir -p` 创建。
- `--help-all` 的实验/内部参数没有变为正式支持；不要根据算法风格名称推断等价外部软件。

迁移只替换模型选择：`--multivariate` → `--model multivariate`，`--repeatability-model` → `--model repeatability`，`--lm` → `--model lm`。
保留原输入、列号、档位、LOCO、方差、线程和输出设置。新旧入口应给出相同数值结果。
旧平铺命令继续可用，不需要批量重写历史脚本；接口升级不修改历史结果的版本身份。

## 7. 设计参考与验证

保留领域通用输入/输出长参数，参考 [PLINK2](https://www.cog-genomics.org/plink/2.0/general_usage)、[GCTA](https://yanglab.westlake.edu.cn/software/gcta/)、[GEMMA](https://github.com/genetics-statistics/GEMMA)、[REGENIE](https://rgcgithub.github.io/regenie/options/)、[HIBLUP](https://www.hiblup.com/tutorials)。
仅借鉴命令组织习惯，不声称参数或统计方法与这些软件全部等价。

回归入口：`python3 tests/run_cli_regressions.py`，或 CTest `cli_regressions`。
验证帮助、错误组合、新旧命令和 GRM 的数值一致性、系谱别名以及原始命令来源记录。
