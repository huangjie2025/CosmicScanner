# CosmicScanner v2.0.0 用户说明书

基因组关联分析软件

开发者：黄杰

联系邮箱：huangjie2025@zju.edu.cn

主页：待上线

平台：Windows / Linux / macOS

macOS 软件包：Apple 芯片 arm64 / Intel x86_64；macOS 15 或更新

本手册面向首次使用 CosmicScanner 的遗传育种与关联分析研究者，说明软件获取、安装、数据准备、模型选择、分析操作和结果核对。先完成第4章随包示例，再按照第5章选择适合研究设计的模块。

本手册适用于 CosmicScanner v2.0.0。命令示例分别采用 Linux Bash 和 Windows PowerShell；请使用与操作系统匹配的程序包，并保留包内目录结构。

## 目录

第1章 软件概述
第2章 安装与运行
第3章 输入文件格式
第4章 快速入门
第5章 主要模块详解
第6章 官网文档与结果查看
第7章 输出文件说明
第8章 常见问题与结果复现
附录A 完整命令行帮助
附录B 常用命令行参数
附录C 参数迁移与分析记录

## 第1章 软件概述

### 1.1 软件简介

CosmicScanner 是用于基因组关联分析和关系矩阵计算的命令行软件。典型流程是读取基因型、表型及协变量，完成样本匹配和位点筛选，拟合选定模型，输出位点效应、标准误、检验统计量与结果记录。软件采用本地计算方式，研究数据位于用户自己的计算机或服务器。

软件同时提供单性状、多性状和重复记录分析入口。用户应先根据表型类型、重复测量情况、群体结构与亲缘关系选择模型，再选择计算策略。计算更快、模型更复杂或名称带有 exact，均不能单独证明统计结果更适合研究问题。

### 1.2 主要功能

| 功能类别 | 主要用途 | 使用提示 |
|---|---|---|
| 线性模型与混合模型 | 连续性状单标记关联分析 | 首次使用优先完成 BED 示例 |
| 加性与显性分析 | 区分加性、显性及联合效应 | 编码、输入格式与随机效应结构需同时确认 |
| 多性状分析 | 同时分析多个相关性状 | 需要明确性状列与扫描策略 |
| 重复力与随机回归 | 分析多次测量及时间变化 | 记录、个体、时间三个层次不能混淆 |
| 关系矩阵与方差组分 | 构建 GRM 和拟合方差结构 | 矩阵必须与样本 ID 一一对应 |
| 系谱与单步分析 | 使用系谱和基因型信息 | 各 Method 的估计对象不同 |
| 二分类及其他分布模型 | 分析二分类、计数、等级性状 | 按性状分布选择链接函数并检查拟合 |
| 区域及后关联分析 | 集合检验、条件分析、COJO、Meta 和精细定位 | 核对等位基因、LD 与统计假设 |
| 输入适配及转换 | 整理表格和多种基因型输入 | 输入支持不等于每种模型均支持 |

### 1.3 如何选择分析模块

根据表型和研究设计选择模型。连续性状且每个个体只有一条记录时，先使用 LM 或 LMM；多个相关连续性状使用多性状模型；同一个体多次测量使用重复力模型，需要描述随时间变化的遗传效应时再考虑随机回归。

| 数据与研究问题 | 建议模块 | 分析前确认 |
|---|---|---|
| 单个连续性状 | LM 或 LMM | 是否需要控制亲缘关系与群体结构 |
| 多个连续性状 | 多性状模型 | 性状列、缺失分布与协方差拟合 |
| 同一个体多次测量 | 重复力模型 | 个体 ID、时间和记录层协变量 |
| 时间变化的遗传效应 | RRM | 时间单位、测量覆盖与基函数阶数 |
| 0/1 二分类性状 | Logistic GLM | 病例编码、病例数与位点 MAC |
| 系谱与基因型联合利用 | 单步分析 | 系谱匹配、Method 和方差处理 |
| 汇总统计后续分析 | COJO Meta 或精细定位 | 等位基因、统计字段与 LD 来源 |

低秩、矩阵自由、两阶段预测和扩展分布模型应在完成基础流程后按第5章的适用条件选择。使用近似算法时，先用可计算的样本子集比较基础模型结果，并检查拟合、数值稳定性和回退记录。

### 1.4 平台与计算资源

| 平台 | 当前状态 | 系统要求与说明 |
|---|---|---|
| Linux x86-64 | 提供对应平台程序包 | glibc 2.34 或更新；保留包内 lib 和启动脚本 |
| Windows x86-64 | 提供对应平台程序包 | 中文路径建议 Windows 10 1903 或更新；完整保留 DLL |
| macOS Apple 芯片 | arm64 程序包 | macOS 15 或更新；保留包内 bin 和 lib |
| macOS Intel 芯片 | x86_64 程序包 | macOS 15 或更新；保留包内 bin 和 lib |

稠密关系矩阵的存储随样本数平方增长，多个矩阵、分解和多性状结构会进一步增加内存。内存需求取决于样本数、性状数和所选模型。先用小规模数据试运行，观察内存和耗时，再确定正式分析规模。--memory-gib 是诊断设置，不是保证进程不超额占用内存的硬限制。

## 第2章 安装与运行

### 2.1 获取与核验程序包

获取与本手册匹配的 v2.0.0 平台程序包及 SHA-256 校验文件。保存 PACKAGE_ID.json，便于核对程序构建信息和复现分析。

Linux 程序包名称为 CosmicScanner-v2.0.0-validation-linux-x86_64.tar.gz；Windows 包名称为 CosmicScanner-v2.0.0-windows-x86_64.zip。下列命令中的文件名应与实际下载文件一致。

```bash
sha256sum -c CosmicScanner-v2.0.0-validation-linux-x86_64.tar.gz.sha256
```

```powershell
Get-FileHash -Algorithm SHA256 -LiteralPath CosmicScanner-v2.0.0-windows-x86_64.zip
```

Windows 将输出值与官方同名校验文件比较。校验不一致时重新下载，不继续分析。不要只复制单个 exe 或 Linux 主程序；启动器、依赖库和示例应作为完整目录保存。

### 2.2 Linux 安装与首次运行

解压后进入包根目录，调用 bin/cosmicscanner。直接运行程序包不需要安装 Python、Conda 或编译器。

```bash
tar -xzf CosmicScanner-v2.0.0-validation-linux-x86_64.tar.gz
# 进入解压生成的目录后运行
./bin/cosmicscanner --version
./bin/cosmicscanner --help
```

包内启动器负责加载随包依赖。若出现 glibc 版本不足，应更换满足要求的系统或采用经核验的源码构建方案，不要直接替换系统核心库。

### 2.3 Windows 安装与首次运行

将整个 ZIP 解压到有足够空间的目录，例如 E 盘的软件目录，在该目录打开 PowerShell。直接从压缩包预览界面运行，可能导致依赖库缺失。

```powershell
& '.\bin\cosmicscanner.exe' --version
& '.\bin\cosmicscanner.exe' --help
```

PowerShell 运行路径变量或带引号的可执行文件路径时需要调用运算符 &。含空格的输入和输出路径必须加引号；不需要为了运行软件把数据复制到 C 盘。

### 2.4 从源码构建

源码构建适合开发者和需要特定平台支持的用户。应获取与主项目匹配的 libCosmicCore、匹配的第三方依赖及完整构建说明，再按该版本的 NATIVE_INSTALL_CN.md 操作。孤立下载主目录或混用不同 Core 版本会造成编译或行为差异。

源码构建得到的常见程序位置为 ./build/cosmicscanner；源码示例位置为 examples/data/r33_synthetic/dataset。程序包使用 ./bin/cosmicscanner 和 examples/dataset。两组路径不要混用。构建后先运行第4章示例，确认程序能够读取数据并输出预期结果。

### 2.5 命令组织与帮助

推荐命令顺序为任务、模型、输入、模型专用参数、资源和输出。任务词 gwas 与 grm 用于提高可读性，原有平铺命令仍可使用。

```text
cosmicscanner gwas --model lmm --bfile PREFIX --pheno FILE --out PREFIX
cosmicscanner grm --bfile PREFIX --out PREFIX
cosmicscanner help models
cosmicscanner help input
cosmicscanner --help-all
```

--model 选择模型家族；--gwas-model 选择单性状内部效应结构；--scan-tier 选择逐位点方差处理策略。三个参数分别控制模型、效应结构和计算策略。常用参数见附录A；输入与组合要求见第5章。

## 第3章 输入文件格式

### 3.1 基因型文件

首次分析建议使用 PLINK BED。PREFIX.bed、PREFIX.bim、PREFIX.fam 必须来自同一数据集，样本顺序和位点顺序必须匹配。--bfile 只填写前缀，不加 .bed 后缀。

| 输入 | 主要参数 | 文件及限制 |
|---|---|---|
| PLINK BED | --bfile PREFIX | BED 加 BIM 加 FAM，建议作为入门格式 |
| PLINK PGEN | --pfile PREFIX | PGEN 加 PVAR 加 PSAM；各模型支持需分别确认 |
| BGEN | --bgen FILE --sample FILE | 支持实现范围内的 Layout 2 zlib；概率输入路径与硬调用不同 |
| VCF | --vcf FILE | 压缩索引读取需匹配索引；参见 --vcf-tbi |
| HapMap | --hmp FILE | 约定的11个元信息列；需检查等位基因编码 |
| 文本矩阵 | --geno-matrix FILE | 另配 --geno-var 和 --geno-sample，明确矩阵方向 |
| NPY | --npy FILE | float32 小端，另配变异与样本信息 |
| PED MAP | --ped-genotypes FILE --map FILE | 该 PED 是基因型，不能当作系谱 |

BGEN 剂量、硬调用、非加性编码的处理不是同一个数据路径。某格式能用于加性模型，不意味着显性、多性状、记录模型或精确扫描也支持该格式。

### 3.2 单性状表型

推荐使用空格或制表符分隔的纯文本。前两列是 FID 和 IID，后续为表型。下例使用 --pheno quantitative.pheno --mpheno 1。

```text
S001 S001 0.246761
S002 S002 1.519086
S003 S003 0.831204
```

--mpheno 从 FID/IID 之后开始以1计数，不是整张文件的绝对列号。FID IID HEIGHT WEIGHT 中，HEIGHT 对应1，WEIGHT 对应2。不要照搬其他软件的表型列号。

缺失值需按对应读取路径使用 NA 或默认 -9 等约定；分析前检查实际保留样本数。0 在二分类和计数中有真实含义，不能作为通用缺失占位符。不得仅靠行号假定表型与基因型匹配。

### 3.3 多性状表型

在两个 ID 后排列多个连续性状。例如 FID IID HEIGHT WEIGHT LENGTH，使用 --mpheno-list 1,2,3。各性状单位、变换和缺失分布应在正式分析前核对。高度共线、重复或零方差性状可能导致方差估计或矩阵分解失败。

多性状联合检验 P 值检验联合假设；单个性状的 BETA 和 SE 解释该性状对应效应。不能把联合 P 值当成每个性状都显著的证据。

### 3.4 重复测量与时间列

长表允许同一 FID/IID 多行，时间字段和数值字段分开。以下输入使用 --tpheno 1 --mpheno 2。

```text
FID IID TIME VALUE
S001 S001 1 1.947087
S001 S001 2 1.463607
S002 S002 1 0.723190
```

同一个体多条记录并不增加基因型个体数。核对记录数、独立个体数、各个体时间范围与缺失分布；时间标签、实际年龄和批次不应在不说明的情况下混用。RRM 对时间尺度和设计有更强要求，参见5.11。

### 3.5 协变量

--qcovar 用于连续协变量，--covar 用于分类协变量。协变量文件需使用与基因型和表型一致的 FID/IID。年龄、主成分等通常作为连续协变量；批次和分组可根据设计作为分类协变量。

```text
FID IID AGE PC1
S001 S001 120 0.013
S002 S002 125 -0.008
```

不要同时引入完全重复的列、由其他列线性组合得到的列或与截距完全重合的常数列。分类变量的编码和参考水平影响系数解释。不同记录模型对个体层和记录层协变量的支持不同，不要把个体一行的文件直接当成记录一行的文件。

### 3.6 关系矩阵与系谱

GRM 的二进制矩阵必须配套其 ID 文件，不能在矩阵生成后单独重新排列 ID。FAM 的样本顺序、矩阵 ID 与模型实际匹配样本应一致。GRM 不同构建方法具有不同缩放和方差解释。

系谱参数推荐 --pedigree FILE，内容按 ID SIRE DAM 的要求准备。父母缺失编码及无基因型个体的处理应遵循对应模块要求。应检查重复个体、自身作父母、循环系谱和无法匹配的个体。--ped 的含义随任务变化，推荐明确区分 --pedigree 与 --ped-genotypes。

### 3.7 集合与汇总统计输入

集合检验需提供位点分组文件；COJO、Meta 和精细定位需要符合各模块要求的字段格式 的汇总统计及必要的 LD 信息。位点 ID、染色体、物理位置、效应等位基因和另一等位基因须可核对。

不要只将外部结果列改名就认定已经完成 harmonization。需要处理等位基因交换、链方向、回文位点、基因组版本差异、重复位点和频率不一致。具体字段以当前帮助和对应模块要求为准，提交正式分析前先以少量位点检查读取日志。

### 3.8 原始表格适配

原始 CSV、TSV 或 TXT 可通过 --input-file 及配置流程整理为模型输入。--preview-schema 和 --infer-columns 用于预览与辅助识别；--prepare-input-only 可只准备文件，不直接分析。当前流程不是 Excel 工作簿读取器，xlsx 应先导出为纯文本。

--id-col、--iid-col、--time-col、--trait-cols 属于原始表格字段选择；--tpheno 属于已整理表型文件的时间列，两者不是同义参数。保存适配配置、重复记录处理策略及输出目录，人工检查日期解析、宽转长和个体去重后再运行模型。

## 第4章 快速入门

### 4.1 示例数据与运行位置

本章命令从完整程序包根目录运行。examples/dataset 中包含60个个体、180个 SNP、两条染色体的合成数据，以及单性状、多性状和重复测量表型。示例用于验证程序与流程，不用于证明真实研究效应或性能。

Linux 和 Windows 各有自己的终端语法。不要把 Bash 的反斜杠续行和变量赋值直接复制到 PowerShell。每次开启新终端，应重新设置 BIN 和 DATA。

### 4.2 Linux 完整操作



下面的 `$BIN` 和 `$DATA` 都是 Bash 变量，每次新终端需重新设置。

```bash
BIN=./bin/cosmicscanner
DATA=examples/dataset
mkdir -p demo_results

# 单性状 LMM，默认 LOCO
"$BIN" --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" \
  --mpheno 1 --gwas --model lmm --threads 1 --out demo_results/lmm

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
| `rep.assoc.txt`、`rep.repeatability.null.txt` | 180 行与记录模型诊断；与 LMM 不保证相同字段格式 |
| `grm.grm.bin`、`grm.grm.id` | 矩阵与匹配样本 ID |

```bash
head -2 demo_results/lmm.assoc.txt
head -7 demo_results/lmm.cosmic.sumstats
cat demo_results/lmm.scan.meta.txt
```

不要仅凭文件存在判断分析成功。检查退出状态、记录数、SNP 数、收敛诊断和非有限数值。



### 4.3 Windows 完整操作
~~~powershell
$BIN = '.\bin\cosmicscanner.exe'
$DATA = 'examples/dataset'
New-Item -ItemType Directory -Force demo_results | Out-Null
& $BIN --version
& $BIN gwas --model lmm --bfile "$DATA/synth" --pheno "$DATA/quantitative.pheno" --mpheno 1 --threads 1 --out demo_results/lmm
& $BIN gwas --model multivariate --bfile "$DATA/synth" --pheno "$DATA/multivariate.pheno" --mpheno-list 1,2,3 --scan-tier approx --reml-maxit 1200 --reml-tol 1e-10 --no-loco --threads 1 --out demo_results/mv
& $BIN gwas --model repeatability --bfile "$DATA/synth" --pheno "$DATA/long.pheno" --tpheno 1 --mpheno 2 --repeatability-mode exact --scan-tier approx --no-loco --reml-vg 0.5 --reml-vpe 0.2 --reml-ve 0.3 --threads 1 --out demo_results/rep
& $BIN grm --bfile "$DATA/synth" --threads 1 --out demo_results/grm
(Import-Csv -LiteralPath demo_results/lmm.assoc.txt -Delimiter ([char]9)).Count
(Import-Csv -LiteralPath demo_results/mv.assoc.mv.txt -Delimiter ([char]9)).Count
(Import-Csv -LiteralPath demo_results/rep.assoc.txt -Delimiter ([char]9)).Count
~~~

最后三个计数应均为180；GRM的ID文件应包含60个样本。重复力固定方差仅是合成示例设定。


### 4.4 成功判据与下一步

首先检查程序退出状态；Linux 可在命令执行后立即查看 echo $?，PowerShell 可立即查看 $LASTEXITCODE。0表示程序正常结束，还需检查结果文件、有效样本量、位点数、收敛和有限数值。不要用下一条命令覆盖退出状态后再查询。

示例的三个关联表均应有180条数据记录，GRM ID 应有60行。单性状统一关联表使用 P_A，多性状和重复力表使用各自字段，不能按同一个 P 列强行解析。180与60仅是随包合成数据的核对值，自己的数据受缺失和过滤影响。

完成示例后，更换为自己的基因型前缀、表型路径和新输出前缀，先运行小规模分析。真实重复力数据不应直接沿用演示固定方差0.5、0.2、0.3；应估计或使用具有依据的方差并检查诊断。

## 第5章 主要模块详解

本章说明各模块的用途、输入、主要参数与结果。第4章可直接使用随包数据运行；本章命令中的 data/cohort、data/trait.txt 等为用户数据位置，请替换为实际路径。运行 Bash 示例前设置 BIN 为程序路径，并创建 results 目录。

### 5.1 基因型质量控制

用途是排除不符合研究质控要求的样本与位点。--maf 控制最小次要等位基因频率，--hwe 控制 Hardy-Weinberg 筛选，BED 路径可使用 --geno 和 --mind 控制位点及个体缺失。--keep/--remove 控制样本，--extract/--exclude 控制位点。

```bash
"$BIN" gwas --model lmm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 \
  --threads 4 --out results/qc_lmm
```

上述阈值仅用于说明参数，并非适合所有物种、群体与设计的统一标准。育种群体、近交群体和病例对照研究的 HWE 处理应有研究依据。核对每一步保留样本数和位点数，保存筛选条件。

### 5.2 关系矩阵构建

用途是估计个体间基因组关系，供混合模型或后续矩阵分析使用。推荐入口为 grm，也可使用 --make-grm。默认和可选算法的缩放不同，报告结果时写明构建方法及位点筛选。

```bash
"$BIN" grm --bfile data/cohort --threads 4 --out results/cohort
```

主要输出为 .grm.bin、.grm.id 和 .grm.N.bin。--make-grm-sparse 可生成实现支持的稀疏关系结果；稀疏截断会改变关系信息，不能默认视为完整矩阵的无损替代。--make-grm-alg 指定算法：1为 VanRaden，2为 Zeng，3为 Yang，4为 Vitezica，5为 TrueVar，6为 emmvgwas 兼容算法。参数取值见附录A。

### 5.3 线性模型 LM

LM 适用于连续性状的固定效应关联分析。它不通过随机遗传效应自动控制亲缘关系，用户需根据群体结构与设计选择协变量和模型。入口为 --model lm，--lm-mode 控制支持的检验模式。

```bash
"$BIN" gwas --model lm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --qcovar data/qcov.txt --threads 4 --out results/lm
```

常见输出为 .lm.assoc.txt，具体列按文件表头核对。检查表型方差、协变量秩和效应方向；LM 与 LMM 的 P 值差异可能反映亲缘结构，不应只按显著位点数量选择结果。

### 5.4 线性混合模型 LMM 与 LOCO

LMM 用随机遗传效应刻画样本相关性，并逐位点检验固定效应。入口为 --model lmm。多染色体标准路径默认采用 LOCO，即扫描某条染色体时，关系结构排除该染色体信息；--no-loco 明确关闭。

```bash
"$BIN" gwas --model lmm --bfile data/cohort --pheno data/trait.txt \
  --mpheno 1 --qcovar data/qcov.txt --threads 4 --out results/lmm
```

当前统一输出 .assoc.txt 的加性字段为 BETA_A、SE_A、P_A 和 LOG10P_A。LOCO 改变零模型与计算量，比较两个分析时必须记录是否启用。仅一条染色体或受限路径应核对日志的实际行为。

--reml-maxit 与 --reml-tol 控制相应拟合迭代；提高迭代上限只能解决部分未收敛情况。先检查输入、协变量、矩阵和方差边界。--lmm-gemma-mode 的 Wald、LRT 和 score 模式具有不同检验假设，兼容命名不等于所有实现细节完全一致。

### 5.5 加性 显性与位点编码

--gwas-model additive、additive_dominance、dominance 选择单性状内部效应结构；加性显性联合模型和显性模型不能仅用一个表型字段解释。加性显性联合模型使用 BED/PGEN 硬调用；显性模型还可使用对应 BGEN 概率路径。

```bash
"$BIN" gwas --model lmm --gwas-model additive_dominance \
  --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --threads 4 --out results/additive_dominance
```

--snp-model additive、dominant、recessive、heterozygote 指定测试位点的编码方式，区别于上述随机效应结构。输出中分别核对加性、显性与联合字段，不能把联合检验当作单一加性 P 值。

BGEN 显性路径还需核对缺失、INFO_D 与概率熵诊断。某些统一输出压缩或附加格式在加性显性路径存在限制，应遵循程序的显式检查，不要绕过检查强行拼接。

### 5.6 低秩与矩阵自由 LMM

--fastlmm、--lmm-low-rank 和 --fastlmm-matrix-free 用于低秩或矩阵自由计算，适合需要控制矩阵存储与计算量的分析。近似误差与秩的选择有关，应与可计算的完整模型结果进行比较。--fastlmm-trace-threshold 控制相关低秩近似目标，最小最大秩、随机种子和误差阈值也会影响结果。

--fast-lmm 与 --fastlmm 不是同一个开关：前者为相关固定零模型扫描别名，后者涉及低秩计算。不要凭名称相近替换命令。使用前在可计算的小数据上比较标准 LMM 的效应、标准误、P 值和拟合诊断，并记录是否发生回退。

### 5.7 FarmCPU 与 BLINK

FarmCPU 和 BLINK 以迭代或位点选择方式处理多位点关联，入口分别包括 --farmcpu 与 --blink。FarmCPU 的 --farmcpu-max-iterations 和 --farmcpu-pseudo-qtn 控制迭代与候选伪 QTN 数量。使用时应检查迭代停止原因和所选位点。

使用时保存迭代日志、入选位点和停止原因。改变迭代数、候选数或数据筛选会改变最终模型，不应将某个模型得到更多显著结果作为准确性证明。建议同时保留同一数据的 LMM 结果，用于比较效应、显著性与群体结构控制情况。


FarmCPU 命令示例。请替换为实际数据路径。

```bash
"$BIN" --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --farmcpu --farmcpu-max-iterations 5 --out results/farmcpu
```

### 5.8 二分类 Logistic 与 Firth

0/1 logistic GLM 适用于独立样本二分类结局，--logistic score 或 --logistic wald 选择支持的检验。确认病例编码为1、对照为0，以及缺失与协变量处理正确。

```bash
"$BIN" --bfile data/cohort --pheno data/binary.txt --mpheno 1 \
  --logistic score --threads 4 --out results/logistic
```

完全分离、极少病例及低 MAC 可能使普通估计不稳定。--firth-logistic 可选择 Firth 校正路径，需检查估计与收敛。--logistic-glmm 采用 PQL 近似；--saige-style 选择相应混合模型 score 与 SPA 计算，并要求 --biobank-mode。三者的拟合与检验不同，应按具体分析方法记录。对病例极少或低 MAC 位点，需额外检查估计和检验的稳定性。

### 5.9 多性状混合模型

用途是联合分析相关连续性状。推荐入口为 --model multivariate，并通过 --mpheno-list 选择两个 ID 后的性状。第4章给出三性状、BED、fixed-null 的已核对示例。

--mv-method 选择具体实现，--mv-test 选择支持的检验，--scan-tier 控制逐位点方差处理。原生 approx 表示相应固定零模型近似；原生多性状计算使用 approx 或满足条件的 exact；semi 只适用于 --mv-compatible emmvgwas，使用前应确认需要该兼容模型及其方差处理方法。

exact 扫描档位要求 --no-loco，并按路径满足输入条件。exact 表示逐位点重新估计方差组分，P 值仍按对应检验的参考分布计算。输出 .assoc.mv.txt 中 WALD_STAT 与 P 为联合检验，各性状效应分别位于 BETA_1/SE_1 等列。检查残差和遗传协方差拟合、数值秩与收敛。

### 5.10 重复力模型

用途是分析同一个体多次测量，区分遗传效应、永久环境与残差。推荐 --model repeatability，配合 --tpheno 与 --mpheno 读取长表。第4章采用记录级计算和固定零模型扫描，演示方差为人为设定。

--repeatability-mode auto、exact、fast、stcg 控制记录模型计算路线；--scan-tier approx、semi、exact 控制逐位点方差更新。--repeatability-mode exact 不等于 --scan-tier exact。逐位点 exact 当前限制 BED、无记录权重、关闭 LOCO。

--repeatability-weight 支持规定的个体或记录权重文件；权重应有测量精度或设计依据。缓存应与基因型、筛选和配置一致，不要跨数据强行复用。输出 .repeatability.null.txt 应与 .assoc.txt 一起保存，检查个体数、记录数、方差和拟合状态。

### 5.11 随机回归模型 RRM

RRM 用于建模遗传效应随时间变化。入口 --model rrm；--rrm-order 设置实现支持的1至4阶，--rrm-basis 选择 legendre 或 bspline。平衡宽表和非平衡长表需要对应数据与参数，不能只替换模型名称。

--rrm-gwas-method record 或 coefficient 选择记录或系数层扫描，--rrm-curve-grid 与 --rrm-top 控制曲线输出相关设置，--rrm-unbalance 用于相应非平衡流程。兼容参数 gma-fixreg 和 gmat-unbalance 对应特定数据与方差处理约定，使用时需匹配输入。

分析前确认时间单位、覆盖范围、每个个体的记录数与基函数阶数；过高阶数和稀疏时间点可能导致不可识别或不稳定。曲线与位点效应必须结合基函数和时间尺度解释。

### 5.12 系谱 Animal 母体与共同环境模型

--model animal 使用系谱关系结构，需 --pedigree；其 A 核模型与基因组 LOCO 概念不同。--model maternal 增加母体相关结构，要求按模块提供母体 ID；--maternal-id 的计数规则不能套用 --mpheno 的相对表型列号。

--model common-env 使用共享环境分组结构，需要正确的环境分组文件，并区分个体层和记录层。母体、栏位、批次等因素若与其他结构完全重合，方差可能难以分开估计。

应先检查系谱和分组是否匹配，分析拟合诊断，并用具有已知结构的数据或独立实现核对。不要把未测量的共同环境解释成遗传贡献。


Animal 命令示例。请替换为实际数据路径。

```bash
"$BIN" gwas --model animal --pedigree data/pedigree.txt \
  --bfile data/cohort --pheno data/trait.txt --mpheno 1 --out results/animal
```

### 5.13 单步分析 Method 1

统一入口为 --model ssgwas，默认对应 Method 1。单步模型结合系谱与基因型关系；--single-step、--make-agh、--write-agh 可用于对应分析和矩阵输出。tau、omega 和 blend 改变矩阵缩放或混合方式，应有明确依据。

Method 1 在有基因型个体的 H 子集上拟合模型并扫描固定 SNP 效应，不输出 Method 1 EBV。需要个体 EBV 和位点回代时使用 Method 2。准备系谱、基因型个体清单和表型时，要分别核对有基因型与无基因型个体。

### 5.14 单步回代 Method 2

--single-step-method backsolve 明确选择 Method 2。本节示例采用 BED、VanRaden 关系结构、LDLT 与固定 Va/Ve；可输出个体 EBV、位点效应和窗口汇总。Method 1 与 Method 2 的估计对象不同，结果不能不加说明地混用。

Method 2 默认不输出位点显著性 P 值。显式启用 --ssgwas-pvalue-method exact 的支持路径才会增加相应不确定性与 P 值输出。该 exact 选项通过 LDLT 下的选定逆计算位点不确定性，P 值按对应参考分布计算。PCG 与 APY 使用不同求解或近似策略，应单独检查收敛、核心个体选择和不确定性计算。

主要文件 .ss_backsolve.ebv.txt、.ss_backsolve.snp.txt、.ss_backsolve.window.txt、.ss_backsolve.meta.txt；支持的 P 值流程另有 .pvalue.txt。报告效应时保留等位基因定义，窗口解释需说明划分方式，EBV 应按个体层输出解释。


Method 2 固定方差示例。这些数值只演示命令形态，固定方差应替换为有依据的研究值。

```bash
"$BIN" gwas --model ssgwas --single-step-method backsolve \
  --bfile data/cohort --pedigree data/pedigree.txt \
  --pheno data/trait.txt --mpheno 1 --reml-vg 0.5 --reml-ve 0.5 \
  --ssgwas-solver ldlt --out results/backsolve
```

### 5.15 阈值 计数与等级模型

--model threshold 用于二分类阈值模型，可使用 probit 或 logit 链接。连续性状只有在研究设计有明确依据时才显式二分；--threshold-dichotomize 不是普遍推荐的预处理。

--model count 提供 poisson、negbin 或 binomial 的相应路径。--offset 为规定的逐样本数值对数暴露量输入，不能误作 FID/IID 表；二项输入 --binomial-pheno 按 FID IID EVENTS TRIALS 准备。核对计数非负、事件数不超过试验数及暴露量定义。

--model ordinal 用于等级性状。experimental-cumulative-binary 与 stacked-psr 是不同实现；fixed-gls 和 dense-zkz、固定或网格方差策略也应明确记录。--ordinal-allow-fallback 允许的回退应查日志，不能把回退结果按原模型报告。

阈值模型的汇总收敛标志不能单独作为拟合成功依据，应同时检查迭代、方差边界和数值结果。等级间距不能默认当作等距连续值。


计数模型命令示例。请替换为实际数据路径。

```bash
"$BIN" gwas --model count --bfile data/cohort \
  --pheno data/count.txt --mpheno 1 --family poisson --out results/count
```

### 5.16 基因与环境互作 GxE

入口 --model gxe，--env-covar 指定该模型约定的环境协变量。相关路径包含基因主效应和基因环境互作，联合检验为相应2自由度假设。

主要字段包括 BETA_G、SE_G、BETA_GxE、SE_GxE、JOINT_CHISQ 与 JOINT_P。联合 P 显著不表示互作项单独显著；环境变量的单位、中心化和编码会影响主效应解释。需检查环境覆盖与群体混杂。

### 5.17 方差组分与 REML

--reml 或 --vce 提供对应方差估计入口；--reml-method 的可用求解器包括 ai、em、brent、idul、idul_plus 和 aidul 等，实际允许组合取决于模型。扫描内部零模型的求解器与独立方差分析不能仅凭参数同名认定相同。

--reml-maxit 和 --reml-tol 控制相关迭代，--reml-vg、--reml-vpe、--reml-ve 在对应路径可指定方差。固定方差有研究依据时才使用，并在报告中说明，不能只为消除未收敛而任意填值。

核对方差非负边界、似然或目标值、梯度或迭代状态、矩阵分解与样本规模。不同关系矩阵缩放、模型固定项和样本子集会改变方差解释，遗传力不应脱离这些条件横向比较。

### 5.18 BSLMM BVSR 与 EB 扩展

--bslmm、--bvsr 和 --bslmm-mcmc 对应不同的近似或采样路径。默认近似输出用于对应近似模型，后验采样需显式选择 --bslmm-mcmc 或 --bvsr；--bslmm 2 的岭回归路径也不能与采样模型混为一谈。

采样应记录 --mcmc-burnin、--mcmc-samples、--mcmc-thin 和 --seed，检查混合、有效样本量与不同种子的稳定性。某些输出的位点排序或近似证据不等于已校准的后验包含概率。

EB LMM 与 EB-Wald 使用经验贝叶斯相关估计和检验。选择 EB 调整时，保留标准 LMM 作为对照，分别记录原始检验和调整指标；两类指标应分列报告。EB 调整的显著性解释需根据所用方法评估，不能直接沿用普通 Wald 阈值。

### 5.19 Burden SKAT 与集合检验

集合检验将预先定义的位点组作为分析对象。--burden 和 --skat 配合 --set-file 使用；--skat-o 使用当前实现的固定网格，--acat-v 使用集合内 Cauchy 合并。二者的权重、相关性处理和 P 值校准应根据研究设计核对。集合划分、位点权重和等位基因编码应在分析前确定。

--rare-max-maf 的当前默认值为0.5，不会自动把输入限定为通常意义的稀有变异。若研究设计要求 MAF 不高于0.01，应显式设定相应阈值；--rare-min-mac 也需结合样本量与设计设置。

Burden 对方向一致的聚合效应有特定解释；SKAT 检验不同集合假设，当前相关尾概率为近似实现。保存集合中实际纳入的位点、过滤条件和统计方法，不把不同集合大小的结果只按原始统计量比较。


Burden 命令示例。请替换为实际数据路径。

```bash
"$BIN" --bfile data/cohort --pheno data/trait.txt --mpheno 1 \
  --burden --set-file data/sets.txt --rare-max-maf 0.01 --out results/burden
```

### 5.20 LD 剪枝 聚类与条件分析

--ld-prune 提供实现支持的 LD 剪枝，--ld-window 为物理距离窗口，--ld-step 与 --ld-r2 控制步长和相关阈值。--ld-clump 使用含 SNP、P、CHR、BP 的统计表，--clump-p 等参数控制筛选。单位与默认值见附录A。

条件分析通过 --condition-snp 指定条件位点，--conditional-joint 控制对应联合处理。条件位点缺失、ID 冲突和高度相关可能导致不稳定，应检查实际进入模型的位点。

这些扩展不应与 PLINK 等工具的完整行为默认等同。LD 使用参考群体时，要核对祖源、位点匹配与基因组坐标版本；未核对的外部 LD 会影响后续解释。

### 5.21 COJO

--cojo、--cojo-slct 和 --cojo-joint 对应支持的条件联合分析流程。--cojo-summary 提供汇总统计，--cojo-ld-ref 提供所需参考 LD；--cojo-p、--cojo-r2 和 --cojo-wind 控制相关筛选、共线性和范围条件。

LD 参考中的 signed-r 与位点顺序、等位基因方向必须一致；不可把 r平方矩阵当作 signed-r 矩阵。显式条件位点通过 --cojo-cond 指定。输出应解释为对应参考 LD 和选择条件下的条件效应，而不是已证实的独立因果位点。


COJO 命令示例。请替换为实际数据路径。

```bash
"$BIN" --cojo --cojo-slct --cojo-summary data/summary.txt \
  --cojo-ld-ref data/signed_ld.txt --out results/cojo
```

### 5.22 Meta 分析

--meta-analysis 结合 --meta-studies 的文件列表读取多个研究结果。--meta-method 支持 fixed、random-dl 和 random-reml，分别对应相应固定或随机效应策略。

```bash
"$BIN" --meta-analysis --meta-studies data/studies.list \
  --meta-method fixed --meta-min-studies 2 --out results/meta
```

列表中每个文件必须满足该模块要求的字段格式。先做等位基因、单位、频率和基因组版本核对，再解释合并效应与异质性。--meta-gc 的校正会改变结果，是否启用及 lambda 设置必须记录。研究样本重叠需要额外考虑，不能仅靠软件自动消除。

### 5.23 精细定位

--fine-map 通过 --finemap-summary 和 --finemap-region 指定统计输入与区域。--credible-level、--finemap-prior-var 和 --finemap-max-signals 控制支持的可信集与信号假设。

该路径要求实际基因型 LD；不能只凭一张汇总表宣称已完成充分的多信号定位。可信集依赖先验、LD、候选范围与模型假设，不代表集合内所有位点具有相同因果可能性，也不替代功能实验。


精细定位命令示例。请替换为实际数据路径。

```bash
"$BIN" --fine-map --bfile data/cohort --finemap-summary data/summary.txt \
  --finemap-region 1:1000000-2000000 --credible-level 0.95 --out results/finemap
```

### 5.24 两阶段岭回归预测与关联分析

--biobank-mode、--regenie-model-snps、--write-pred、--use-pred 和 --regenie-step 控制相应预测与扫描流程。步骤0为相应合并流程，1和2分别对应预测与扫描阶段，具体输入要求见实际帮助。

第一阶段用分块岭回归生成 LOCO 预测，第二阶段读取相应预测进行关联扫描。保存训练设置、所选位点、预测文件和样本对应关系。样本或基因型变更后，应重新生成匹配的预测文件。使用该近似流程时，先在可计算的数据子集与标准 LMM 比较结果。

### 5.25 PCA 与格式转换

--pca 与 --pca-k 提供 BED/GRM 主成分分析。核对主成分数、样本顺序和矩阵构建方法。主成分可帮助描述群体结构，但包含多少主成分应由研究设计决定。

--convert 配合 --to 将支持的输入转换为 text、npy、ped 或 vcf。输入格式仍需对应的基因型参数。转换前后核对样本数、位点数、ID、坐标和等位基因；概率、剂量与硬调用之间可能丢失信息。不要把转换成功当作后续每种模型均支持该格式的证明。

## 第6章 官网文档与结果查看

### 6.1 官网入口与文档导航

官网提供软件介绍、平台下载信息和在线文档。首次使用请阅读安装说明，随后按第4章运行随包示例并核对结果。

CosmicScanner 在本地计算机或 Linux 服务器运行。在终端执行分析命令，使用 --out 指定结果位置；结果文件保存在该计算环境中。

### 6.2 读取结果表

可使用纯文本查看器、R、Python 或支持制表符的表格软件打开关联表。读取时按表头选列，不按固定列位置。Excel 等工具可能自动改变位点 ID、日期样式或科学计数法表示，建议保留原始输出并用副本查看。

Windows 第4章的 Import-Csv 示例按制表符读取。对超大文件，优先使用流式读取或分析工具，避免将整张表强行打开在表格界面导致内存耗尽。

### 6.3 Manhattan 图与 QQ 图

Manhattan 图展示位置与负对数 P 值，QQ 图比较观测和参考分布。单性状 LMM 当前应选择 P_A，多性状选择联合 P，重复力选择其 P 字段。先过滤或单独记录非有限值和缺失值；数值下溢得到0的 P 值需结合 LOG10P 字段处理，不能盲目计算负对数。

显著性阈值应依据研究设计和多重检验策略预先确定。QQ 图偏离可能涉及真实多基因效应、群体结构、模型失配或数值问题，不能仅凭图形判断是哪一种原因。绘图文件不是统计结果本身，应同时保存原始表和绘图参数。

### 6.4 下载与文档版本核对

选择与操作系统匹配的下载包，核对版本、SHA-256 和 PACKAGE_ID.json。升级程序时同时更新手册与示例数据。macOS 用户请先确认下载页是否提供对应架构的软件包。

## 第7章 输出文件说明

### 7.1 输出前缀与目录

--out 是输出前缀，不只是文件夹名称。--out results/study 会生成 results 目录下的 study 开头文件。运行前创建父目录，并给不同模型、筛选条件和重复运行使用独立前缀，避免覆盖。

### 7.2 常见文件

| 文件后缀 | 对应内容 | 主要核对点 |
|---|---|---|
| .assoc.txt | 统一单性状或某些记录模型关联表 | 必须结合模型读取表头 |
| .lm.assoc.txt | LM 关联结果 | 以实际检验字段为准 |
| .cosmic.sumstats | 支持路径的汇总统计 | 位点完整性和等位基因定义 |
| .assoc.mv.txt | 多性状联合结果 | 联合统计量及各性状效应 |
| .repeatability.null.txt | 重复力零模型诊断 | 个体数、记录数、方差与状态 |
| .grm.bin | 二进制关系矩阵 | 与 ID 和计数文件配套 |
| .grm.id | 关系矩阵样本 ID | 顺序与样本匹配 |
| .grm.N.bin | 关系矩阵位点计数信息 | 与同前缀矩阵配套 |
| .scan.meta.txt | 支持扫描路径的元信息 | 实际模型、策略与样本规模 |
| .run_manifest.json | 支持路径的运行记录 | 命令、配置与运行身份 |
| .ss_backsolve.ebv.txt | Method 2 个体 EBV | 个体 ID 与估计对象 |
| .ss_backsolve.snp.txt | Method 2 位点回代 | 效应及等位基因定义 |
| .ss_backsolve.window.txt | Method 2 窗口汇总 | 窗口范围与汇总方式 |
| .ss_backsolve.meta.txt | Method 2 配置记录 | 方差、矩阵与求解路线 |

不同路径不保证产生所有文件。不应编写只检查统一后缀的脚本，然后把缺失文件当成所有模块失败；应按选定模块的约定核对。

### 7.3 单性状字段

| 字段 | 含义 | 解释要求 |
|---|---|---|
| CHR SNP BP | 染色体、位点 ID 与位置 | 确认基因组版本与重复位点 |
| A1 A2 | 输出等位基因 | 按该模型约定确认效应方向 |
| BETA_A | 加性效应估计 | 单位随表型和编码变化 |
| SE_A | 加性效应标准误 | 非有限或异常值需检查 |
| P_A | 加性检验 P 值 | 需结合检验和多重比较 |
| LOG10P_A | 对应负对数概率字段 | 极端值结合数值记录核对 |

显性和联合模型会增加对应字段。不能将第一列 BETA 或 P 的位置写死，也不要把 A1 一概当作所有输出、所有软件中的同一效应等位基因约定。

### 7.4 多性状与重复力字段

多性状表的示例字段为 CHR SNP BP A1 A2 FREQ WALD_STAT P，以及 BETA_1 SE_1、BETA_2 SE_2、BETA_3 SE_3。P 属于联合检验，各性状的效应应按 --mpheno-list 的顺序对应。

重复力参考表的字段为 CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P。N 的具体统计层次按该路径元信息核对，不应直接与独立个体数混用。记录模型中个体数和测量记录数都需要报告。

### 7.5 日志与异常结果

检查输出 P 值在0至1范围内且具有有效数值。NA、无穷值、负标准误、异常频率或大量边界方差需要排查。全体位点结果完全相同、有效样本数为0或未收敛时，不可仅凭输出文件存在认定成功。

保存标准输出和错误输出可以帮助定位问题，但日志不替代原始结果。对分染色体 LOCO 流程，检查汇总文件包含所有预期染色体，避免把最后一个分块当作完整结果。

## 第8章 常见问题与结果复现

### 8.1 安装与路径问题

| 现象 | 常见原因 | 建议操作 |
|---|---|---|
| 找不到程序或 DLL | 未完整解压或只复制 exe | 重新完整解压并从包根运行 |
| Linux 提示 glibc 不兼容 | 系统低于包要求 | 使用满足要求的环境或匹配源码方案 |
| 找不到 BED 或表型 | 当前目录错误或前缀带错后缀 | 检查完整路径，--bfile 不加 .bed |
| 含空格路径失败 | 路径未加引号 | 为输入、输出和程序路径加引号 |
| 无法写出结果 | 父目录不存在或权限不足 | 先创建目录并使用可写位置 |
| Mac 无法运行现有包 | 下载包与系统不匹配 | 使用下载页列出的对应平台程序 |

### 8.2 输入与模型问题

| 现象 | 检查内容 |
|---|---|
| 有效样本比预期少 | FID/IID 匹配、缺失、样本筛选与协变量完整性 |
| 选错性状或时间 | --mpheno 与 --tpheno 从两个 ID 后计数 |
| 协变量矩阵奇异 | 常数、重复、共线列与分类编码 |
| 多性状未收敛 | 性状方差、相关性、缺失和协变量，随后才调整迭代 |
| 原生 semi 被拒绝 | 原生多性状选择 approx 或 exact；semi 需匹配 emmvgwas 兼容模型 |
| exact 扫描被拒绝 | 检查 no-LOCO、BED、权重及具体模块限制 |
| 找不到 P 列 | 按模型选择 P_A、P 或其他明确字段 |
| Method 2 没有 P 值 | 默认关闭，需要支持路径的显式参数 |

### 8.3 性能与资源问题

先在可控规模上运行，确认算法路线和内存需求。线程数增大不保证线性加速，多个分析同时占满线程可能降低吞吐。关系矩阵与多性状模型可能需要大量内存，避免只根据基因型文件压缩体积估计资源。

非零退出状态意味着需要检查错误信息。不要通过删除输入检查、强行指定任意方差或忽略收敛失败继续输出正式结论。采用近似、低秩或矩阵自由策略时，应说明与参考计算的差异。

### 8.4 分析复现记录

每次分析至少保存软件版本与包身份、输入文件校验值、完整命令、平台、线程数、模型与效应结构、LOCO、扫描档位、样本和位点筛选、方差处理、输出表、诊断与绘图设置。使用采样、迭代或近似计算时，还需保存种子、迭代、回退和近似参数。

完整数据可能受隐私或项目协议限制。反馈问题时提供可共享的小型复现数据、完整命令、错误日志和规模说明，不把私有原始数据、密码或内部服务器记录放入公开页面。

## 附录A 完整命令行帮助

以下收录 CosmicScanner v2.0.0 的简要帮助与完整命令行帮助输出。按参数名称检索具体选项；中文分类说明见附录B，模型输入与组合条件见第5章。

### A.1 简要帮助

```text
CosmicScanner: genetic association analysis
Usage: ./cosmicscanner [gwas|grm] [options]

Common options (PLINK/GCTA-style long flags):
  --bfile PREFIX         PLINK BED/BIM/FAM input
  --pfile PREFIX         PLINK PGEN/PVAR/PSAM input
  --bgen FILE --sample FILE  BGEN input (model-specific support)
  --pheno FILE           FID IID followed by phenotype columns
  --mpheno N             1-based column AFTER FID/IID; default 1
  --mpheno-list 1,2      Phenotype columns for multivariate analysis
  --covar FILE           Categorical covariates (FID IID ...)
  --qcovar FILE          Quantitative covariates (FID IID ...)
  --gwas                 Association task; optional after 'gwas' command
  --model NAME           lmm, lm, multivariate, repeatability, ...
  --scan-tier TIER       approx|semi|exact; see help models for limits
  --no-loco              Global scan; otherwise LOCO where supported
  --maf F --hwe P        Variant filters
  --extract FILE --exclude FILE  Variant-ID inclusion/exclusion
  --threads N            Number of threads
  --out PREFIX           Output file prefix

Example:
  cosmicscanner --bfile data --pheno trait.txt --gwas --model lmm --out run
  cosmicscanner gwas --model multivariate --bfile data --pheno traits.txt --mpheno-list 1,2 --scan-tier approx --reml-maxit 1000 --out mv


Help: --help [topic], help [topic], --help-all (complete legacy/advanced reference)
Manual: docs/CLI_GUIDE_EN.md / docs/CLI_GUIDE_CN.md
Existing flat commands, aliases and statistical defaults are retained.
```

### A.2 完整参数帮助

```text
CosmicScanner v2.0.0
  --version             Print software version and exit
Public scope: A=scoped official contract, B=experimental, C=internal/legacy.
See COSMICSCANNER_V2_SCOPE_FREEZE_PROPOSAL.md; final acceptance is pending.
Formats, engines and output schemas are model-specific; no universal equivalence.
Usage: ./cosmicscanner --bfile <prefix> [options]

=== INPUT/OUTPUT ===
  --pca --pca-k <n>     Experimental BED GRM PCA; default k=10
  --geno <f> --mind <f> Maximum variant/sample missing fraction [0,1]; BED QC
  --keep <f> --remove <f> BED sample lists: FID IID; remove takes precedence
  --ld-prune            Experimental BED window LD pruning
  --ld-clump <f>        Experimental clumping; summary header SNP P CHR BP
  --ld-window <bp> --ld-step <n> --ld-r2 <f> Defaults: 250000, 1, 0.2
  --clump-p <f>         Lead P threshold (default 5e-8)
  --skat-o              Fixed-grid dependence-adjusted SKAT-O
  --acat-v              Experimental set-wise Cauchy combination (--set-file)
  --rare-max-maf <f>    Maximum set-test MAF in (0,0.5] (default 0.5)
  --rare-min-mac <n>    Minimum set-test minor allele count (default 1)
  --snp-model <m>       additive (default), dominant, recessive, heterozygote
                        Fixed-effect marker coding for standard LM/LMM;
                        non-additive --snp-model requires hard-call genotypes.
                        --gwas-model dominance/additive_dominance requires
                        additive --snp-model.
  --gwas-model <m>      additive (default) | additive_dominance | dominance
                        additive:          a + A
                        additive_dominance: a + d + A + D LMM scan (PLINK BED or PGEN hardcall)
                        dominance:         d + D LMM scan (BED/PGEN hardcall or BGEN probability)
  --bgen-dom-max-missing <f>  BGEN dominance variant missingness filter (default 0.05)
  --bgen-dom-min-info <f>     Minimum dominance INFO_D (default 0.01)
  --bgen-dom-entropy-warning <f>  Normalized entropy warning threshold (default 0.50)
  --farmcpu             Experimental native FarmCPU multi-locus GWAS
  --farmcpu-max-iterations <n>  Iteration cap (default 5)
  --farmcpu-pseudo-qtn <n>      Pseudo-QTN cap (default 10)
  --blink               Experimental native BLINK multi-locus GWAS
  --condition-snp <f>   Experimental conditional GWAS; one lead SNP ID per line
  --conditional-joint  Experimental LD/BIC stepwise joint selection for --condition-snp
  --firth-logistic     Experimental Firth penalized binary GWAS
  --logistic <test>    Ordinary logistic GWAS: score or wald
  --logistic-glmm      B: PQL mixed null with full projected working-score variance
                        Not an exact mixed-likelihood fit.
  --saige-style        Experimental binary mixed score/SPA scan; requires --biobank-mode
  --cojo               COJO-style conditional analysis from summary stats + LD reference
  --cojo-slct          COJO stepwise conditional analysis (forward selection)
  --cojo-joint         COJO joint model fit of selected SNPs
  --cojo-summary <f>   GWAS summary statistics file
  --cojo-ld-ref <f>    Dense signed-r LD: SNP CHR BP A1 A2 <ordered SNP IDs>
  --cojo-p <f>         P-value threshold for COJO selection (default: 5e-8)
  --cojo-r2 <f>        LD r^2 threshold for collinearity (default: 0.9)
  --cojo-wind <bp>     Window size in bp for LD (default: 10000000)
  --cojo-cond <f>      Condition on these SNP IDs (one per line)
  --fine-map           Native summary-statistic fine-mapping with real genotype LD
  --finemap-summary <f> Canonical .cosmic.sumstats input
  --finemap-region <r> Locus in chr:start-end form
  --credible-level <f> Credible-set posterior level (default: 0.95)
  --finemap-prior-var <f> Prior effect variance W (default: 0.04)
  --finemap-max-signals <n> Maximum conditional signals (default: 5)
  --finemap-signal-p <f> Conditional signal threshold (default: 5e-8)
  --meta-analysis      Summary-statistics meta-analysis (fixed/random effects)
  --meta-studies <f>   File listing study summary stat files
  --meta-method <m>    Meta-analysis method: fixed, random-dl, random-reml (default: fixed)
  --meta-het-p <f>     Heterogeneity threshold p-value (default: 0.05)
  --meta-min-studies <n> Min studies required (default: 2)
  --meta-gc            Apply genomic control correction
  --meta-gc-lambda <f> Genomic control lambda (default: 1.0)
  --meta-out-format <f> Output format: metal, simple (default: metal)
  --bfile <prefix>       PLINK 1.x binary prefix (.bed/.bim/.fam)
  --bgen <file>          BGEN file (Layout 2, zlib compressed)
  --sample <file>        Sample file for BGEN (Oxford format)
  --pfile <prefix>       PLINK 2.0 prefix (.pgen/.pvar/.psam)
  --pgen <file>          PLINK 2.0 genotype file
  --pvar <file>          PLINK 2.0 variant file
  --psam <file>          PLINK 2.0 sample file
  --vcf <file>           VCF file (bgzip compressed, tabix indexed)
  --vcf-tbi <file>       Tabix index file for VCF (optional if .tbi alongside)
  --vcf-sample <file>    Sample mapping file for VCF (optional)
  --hmp <file>           Standard HapMap file (11 metadata columns + samples)
  --geno-matrix <file>   Text genotype matrix (default: variant rows, sample columns)
  --geno-var <file>      Text variant metadata: ID CHROM POS REF ALT
  --geno-sample <file>   Optional text sample metadata: IID or FID IID
  --geno-orientation <o> variants-by-samples (default) or samples-by-variants
  --npy <file>           NPY genotype matrix (float32, little-endian)
  --npy-var <file>       NPY variant metadata: ID CHROM POS REF ALT
  --npy-sample <file>    NPY sample metadata: FID IID
  --npy-orientation <o>  variants-by-samples (default) or samples-by-variants
  --ped <file>           PLINK PED file
  --map <file>           PLINK MAP file (must be used with --ped)
  --convert <format>     Convert genotype format (input: --bfile/--pfile/--bgen/--vcf/--hmp/--geno-matrix/--npy/--ped+--map)
  --to <format>          Target format for conversion: text, npy, ped, vcf (requires --convert)
  --out <prefix>         Output prefix for conversion
  --batch-models <list>  C: developer batch runner: lm,lmm,fastlmm,saige,farmcpu,blink
  --batch-out <prefix>   Output prefix for batch run (each model writes its own output)
  --out/-o <prefix>      Output prefix (default: cosmic_gwas)
  --out-ssf              EBI GWAS-SSF output on supporting paths only
  --out-gz               Gzip on supporting paths; AD/dominance guards reject it
                          No all-model output schema; sumstats INFO=1 is not measured INFO.
  --extract <file>       File listing SNP IDs to include in analysis
  --exclude <file>       File listing SNP IDs to exclude from analysis

=== PHENOTYPE / COVARIATES ===
  --pheno/-p <file>      Phenotype file: FID IID [TIME] VALUE
  --mpheno/-n <n>        Phenotype column number, 1-indexed (default: 1)
  --tpheno <n>           Time column number for repeatability/RRM
  --qcovar <file>        Quantitative covariates: FID IID VAL1 VAL2 ...
  --covar <file>         Categorical covariates: FID IID CAT1 CAT2 ...

=== INPUT ADAPTATION MVP ===
  --input-file <file>    Raw csv/tsv/txt table to prepare before modeling
  --input-config <json>  Load input adaptation settings from JSON config
  --write-template-config <json> Write an editable input config and exit
  --input-format <fmt>   auto|csv|tsv|txt (default: auto; xlsx not in MVP)
  --input-model <model>  repeatability|rrm|lmm|multivariate (default: repeatability)
  --prepared-input-dir <dir> Directory for prepared manifests and model files
  --prepare-input-only   Write prepared inputs and exit before GWAS
  --id-col <name>        Raw ID column; matched to genotype IID or FID:IID
  --iid-col <name>       Optional IID column when --id-col is FID
  --time-col <name>      Record time/date/batch column
  --trait-cols <list>    Comma-separated trait columns
  --input-trait <name>   Trait to run after preparation (default: first trait)
  --record-qcovars <l>   Record-level quantitative covariates
  --record-ccovars <l>   Record-level categorical covariates
  --subject-qcovars <l>  Subject-level quantitative covariates
  --subject-ccovars <l>  Subject-level categorical covariates
  --time-parser <mode>   auto|date|integer|batch (default: auto)
  --time-unique-policy <p> strict|within-time-index (default)
  --time-index-multiplier <n|auto> Multiplier for duplicate-time keys
  --preview-schema       Write schema_preview.tsv and exit before preparation
  --infer-columns        Suggest ID/time/trait/covariate columns in preview
  --wide-to-long         Treat --trait-cols as wide trait-time columns
  --maf <float>          Minor allele frequency filter (default: 0.01)
  --hwe <float>          HWE exact test p-value filter (default: 0.0, disabled)
  --gwas                 Run standard GWAS; some model flags enable scans implicitly
  --model <name>         Model family: lmm, lm, multivariate, repeatability, rrm,
                          animal, maternal, ssgwas, common-env, threshold, count, ordinal, gxe
                          Legacy --model lmm-eb/eb-lmm remains an EB modifier.
                          Documented analysis commands should include it explicitly.
  --no-loco              Disable leave-one-chromosome-out (run global)
  --scan-tier <tier>     Per-SNP VC handling: approx, semi, exact
                        Multivariate semi requires --mv-compatible emmvgwas; native uses approx
                          approx=fixed null; semi=fixed ratio/profile scale;
                          exact=full per-SNP variance-component refit
                          Not finite-sample exact P. mv/repeatability exact tier needs no-LOCO;
                          repeatability exact tier also requires BED and no record weights.
  --fast-lmm             Legacy alias for fixed-null single-trait scan
  --fastlmm              B: BED Gaussian trace-ranked randomized low-rank LMM
  --fastlmm-trace-threshold <f> Minimum retained GRM trace fraction (default: 0.95)
  --fastlmm-min-rank <n> Minimum automatic rank (default: 1)
  --fastlmm-max-rank <n> Maximum automatic rank (default: 256)
  --fastlmm-matrix-free  B: BED matrix-free low-rank path; rank adequacy is data-dependent
  --fastlmm-no-fallback  Don't fall back to full rank if trace threshold not met
  --lmm-low-rank <k>     Opt-in randomized rank-k LMM (0 keeps exact Eigen default)
  --lmm-low-rank-oversample <n> Randomized eigensolver oversampling (default: 10)
  --lmm-low-rank-seed <n> Deterministic randomized eigensolver seed (default: 42)
  --lmm-low-rank-error-tol <f> Relative Frobenius error limit; exact fallback above it
  --lmm-gemma-mode/-lmm <n> GEMMA LMM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --lmm-gemma-refit      Per-SNP alternative REML refit for GEMMA LMM Wald/LRT
  --lm                   Ordinary linear model flag; use --lm-mode <n>
  -lm <n>                GEMMA-style short LM alias
  --lm-mode <n>          GEMMA LM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --model lmm-eb         Empirical Bayes random-SNP LMM GWAS (experimental; report standard LMM too)
  --gwas-method eb-lmm   Alias for --model lmm-eb
  --test eb-wald         Alias enabling EB-LMM Wald test
  --adjust eb-mbonf      Experimental threshold in .eb.summary.txt; do not report alone
  --spa                  Use Saddle Point Approximation for binary traits
  --threads <n>          Number of OpenMP threads (default: max)
  --memory-gib <value>   SolverPlanner diagnostic memory budget in GiB
  --chr <chr>            Restrict scan to one chromosome
  --from-bp <pos>        Restrict scan start position
  --to-bp <pos>          Restrict scan end position

=== VCE OPTIONS (model-dependent) ===
  --lmm-null-solver <s>  Null model solver: brent, idul, idul_plus, aidul
                          (default: brent)
  --reml-maxit <n>       Max VCE iterations (default: 50)
  --reml-tol <f>         VCE convergence tolerance (default: 1e-6)
  --reml-vg <f>          Fix genetic variance (skip VCE)
  --reml-vpe <f>         Fix permanent environment variance (skip VCE)
  --reml-ve <f>          Fix residual variance (skip VCE)

=== MODEL SELECTION (choose exactly one) ===

--- Standard LMM (default, no flag needed) ---
  Exact single-trait mixed model GWAS with Eigen decomposition.
  Standard for quantitative traits with one record per individual.

--- Repeatability Model ---
  --repeatability-model  Multi-record per individual (long format).
  --pheno format: FID IID TIME VALUE  (--tpheno <col>)
  VCE: additive genetic + permanent environment + residual.
  --repeatability-mode <m> Strategy: auto, exact, fast, stcg
                          A: BED record-level approx/semi; B: fast/STCG/full refit.
                          Auto may switch strategy; exact mode differs from exact tier.
  --repeatability-grm-cache <p> Save/load subject GRM
  --repeatability-weight <file> Record weights: FID IID WEIGHT or FID IID TIME WEIGHT
  --vce-mode <m>         Repeatability VCE refinement: auto, fast,
                          subject-level-exact, full-record-exact
                          auto uses scalable subject-level Woodbury VCE;
                          full-record-exact is validation-only.
  --repeatability-refit-validation Per-SNP VCE refit (small data)
  --repeatability-refit-limit <n> Refit limit (0 = all)

--- Common-Environment Model (Experimental) ---
  --common-env-model     y = Xb + x_j*beta_j + Z_a*a + Z_c*c + e
                         a ~ N(0, K*sigma2_a) (genomic GRM), c ~ N(0, I*sigma2_c)
  --pheno format: FID IID [TIME] VALUE  (long format, --tpheno <col>)
  --common-env-file <f>  FID IID GROUP per subject, or per record with
  --common-env-record-level
  --reml-vg/-vc/-ve <v>  Fix variance components (skip VCE)

--- Random Regression Model (RRM, B experimental) ---
  --rrm-order <n>        Basis order/df-1 (Legendre default, 1-4)
  --rrm-basis <m>        Basis: legendre (default) or bspline
  --rrm-curve-grid <n>   Time grid points for dynamic curves (default: 21)
  --rrm-top <n>          Number of top SNPs for RRM curve outputs (default: 100)
  --rrm-ai-iter <n>      RRM null AI refinement iterations (-1=auto, 0=fast)
  --rrm-sig-threshold <p> P threshold for dynamic time windows (default: 0.05)
  --rrm-unbalance         Unbalanced longitudinal GWAS (long-format --pheno)
  --rrm-unbalance        Unbalanced longitudinal RRM scan (default: coefficient method)
  --rrm-gwas-method <m>  Unbalanced RRM-GWAS: record (record-level x*Phi Method A) | coefficient (two-stage Method B, default)
  --rrm-loco-global-null  Fast approximate RRM LOCO using global null model
  --tpoint <time_col>     Column name for time points (unbalanced mode)
  --tclass <cols>         Time-dependent class covariates (comma-sep column names)
  --tcovar <cols>         Time-dependent q-covariates (comma-sep column names)
  --rrm-compatible <m>   Compatibility mode: gma-fixreg (balanced PLINK) | gmat-unbalance
  --rrm-compatible-varcom <file> GMA .var file for gma-fixreg mode
  --long-pheno <file>    Long-format pheno: FID IID TIME VALUE
  Fits subject-level Legendre coefficients with GRM-based
  random regression. Outputs joint multi-df chi-square test.

--- Multivariate Model ---
  --multivariate         Requires --mpheno-list <c1,c2,...>
                          A: scoped BED reference; PGEN/BGEN extensions B; VCF unsupported.
                          Hybrid requires BED; exact scan tier requires no-LOCO.
  --mpheno-list <list>   Comma-separated phenotype column numbers
  --mv-method <m>        emmv, wald, eb, gemma, gemma-full, gemma-exact, gemma-gxe, gemma-gxe-full, hybrid
  --mv-test <list>       GEMMA-full tests: wald,lrt,score
  --mv-gemma-mode <n>    GEMMA mvLMM mode: 1 Wald, 2 LRT, 3 Score, 4 all
  --mv-pnr <p>           GEMMA mvLMM refinement threshold (default 0.001)
  --mv-emi <n>           GEMMA exact EM iterations (default 10000)
  --mv-nri <n>           GEMMA exact refinement iterations (default 100)
  --mv-emp <f>           GEMMA exact EM precision (default 0.0001)
  --mv-nrp <f>           GEMMA exact refinement precision (default 0.0001)
  --mv-crt              GEMMA-style Edgeworth/PCRT p-value correction after refinement
  --mv-refine-discrepancy-top <n> Hybrid: add top Score/Wald discrepancy SNPs
  --mv-refine-discrepancy <x>     Hybrid: min |log10P score-wald diff| (1.0)
  --mv-refine-score-prepass <n>   Hybrid: Score-only audit top n Stage1 SNPs
  --mv-vce-solver <s>    em (default) or aireml
  --mv-em-min-iter <n>   Min EM iterations before AI-REML (10)
  --mv-ai-max-iter <n>   Max AI-REML iterations (50)
  --mv-output-vbeta      Add upper-triangle beta covariance columns
  --mv-emmvgwas-output   Also write EMmvGWAS-style G/R/beta/association CSV files
  --mv-emmvgwas-G <csv>  EMmvGWAS semi_gwas mode: external G.csv from vc_estimation
  --mv-emmvgwas-R <csv>  EMmvGWAS semi_gwas mode: external R.csv from vc_estimation
  --emmvgwas             Alias: --multivariate --mv-method emmv --mv-compatible emmvgwas --no-loco
  --emmvgwas-vc-only     EMmvGWAS vc_estimation route: write G/R/beta and stop before scan
  --emmvgwas-output      Alias for --mv-emmvgwas-output
  --emmvgwas-G/R <csv>   Aliases for EMmvGWAS semi_gwas external G/R inputs
  --env-covar <col>      Environment column for --mv-method gemma-gxe or --gxe-model

--- Maternal Effects Model (B experimental) ---
  --maternal-model       Direct + maternal + PE multi-component model
  --pheno format: FID IID DAM_ID PHENO  (4+ columns)
  --maternal-id <col>    Maternal ID column (default: 3, 1-indexed)
  VCE: direct additive (ad), maternal additive (am),
       direct-maternal covariance (adm), PE (c), residual (e)
  Scan uses effective VCs: Vg=ad+am+(c) with LOCO GRM.

--- Single-Step GWAS (ssGWAS) Model ---
  --single-step          B Method 1 by default: genotyped H-subset fixed-SNP LMM
                          Not full-pedigree MME association; no Method-1 EBV writer.
  --make-agh             Build A-inverse, G, and H-inverse matrices and exit
  --write-agh            Write A-inverse, G, and H-inverse during ssGWAS
  --ped <file>           Pedigree: ID SIRE DAM (no header)
  --genotyped-list <f>   File listing genotyped individual IDs
  --write-ebv            C: reserved Method-1 output switch, writer not implemented
  --h2-tau <f>          Genomic inverse scaling (default: 1.0)
  --h2-omega <f>        A22 inverse correction scaling (default: 1.0)
  --blend-g <f>         G blend weight with I (default: 0.05)
  Method 1 estimates VCs on the genotyped H subset, then scans fixed SNPs.

--- Method 2: full-pedigree ssGBLUP back-solving ---
  --single-step-method lmm|backsolve
                        lmm       = Method 1 fixed-SNP association (default, backward compatible)
                        backsolve = Method 2: ssGBLUP genetic evaluation -> GEBV -> SNP
                                    back-solving -> SNP variance sums/window summaries
                        A reference: BED, fixed positive Va/Ve, VanRaden, LDLT; P default off.
                        PGEN/BGEN backsolve rejected; VC/PCG/operator extensions are B.
  --reml-vg <f> --reml-ve <f>
                        Fixed variance components for Method 2 (canonical/deterministic path)
  --ssgwas-vce <mode>   VCE estimation mode for Method 2 (ai|emai|em|hi|exact|...).
                        If omitted, --reml-vg/--reml-ve are required (fixed-VC path)
  --ssgwas-window-snps <n>
                        SNP-count windows (default: 10); denominator total SNP variance,
                        not model Va or LD-adjusted window variance.
  --ssgwas-tune-g       Align G to A22 (G = b*G + a, preGSf90-style) before inversion
  --ssgwas-blend-a22    Blend G* with A22 instead of I (HIBLUP/preGSf90 reference; backsolve only)
  --ssgwas-solver ldlt|pcg
                        MME solver (default: ldlt reference; pcg = matrix-free PCG)
  --ssgwas-pcg-tol <x>  Matrix-free MME PCG tolerance (default: 1e-8)
  --ssgwas-pcg-maxit <n> Matrix-free MME PCG maximum iterations (default: 2000)
  --ssgwas-pvalue-method off|exact|apy-exact|apy-approx
                        Frequentist P-values for back-solved SNP effects (Phases C1/C2/D).
                        Default off. exact = C1 reference (selected MME inverse via LDLT,
                        O(n_genotyped^2) memory). apy-exact = APY route (Leite et al.
                        2024 GSE 56:59): APY G^-1 in the MME + exact core-slice selected
                        inverse (requires --ssgwas-solver ldlt; O(n_core^2) uncertainty
                        block; choose an explicit estimator below).
  --ssgwas-apy-exact-estimator <core-v2|legacy-hybrid-v1>
                        core-v2: matched projected/core SNP effects and exact selected-core
                        covariance (Bermann 2022 Eq. 21-22), kappa=1/S. Initial scope:
                        fixed VC, tau=omega=1, no tuning, identity blend 0.05, dense Gcc
                        solver and fixed core selection. B: limited validation, no universal calibration.
                        Default legacy-hybrid-v1 preserves full-point/core-variance outputs;
                        C compatibility: unmatched pairing, not supported calibrated inference.
                        apy-approx = Phase D Approx_GinvAPY (Leite Eq. 8/9):
                        APY G^-1 + approximate core PEV C_core = (D_c + lambda*Gcc^-1 +
                        Q*S*Q')^-1 (Schur/Woodbury form of Eq. 8; D = Eq. 9 diagonal
                        effective records); NO selected MME inverse and NO dense n_g^2
                        uncertainty block; works with --ssgwas-solver pcg; status
                        EXPERIMENTAL_APPROXIMATE: raw P-values, no genomic control,
                        strict type-I control NOT guaranteed. Use for research/exploration,
                        not as a calibrated inferential P-value path.
  --ssgwas-apy-core-size <n>
                        Core animals for the APY G^-1 (apy-* modes; >= 1, <= n_genotyped;
                        n = n_genotyped reproduces the C1 exact path)
                        'auto' = C developer-only sizing; excluded from v2 supported inference
                        due to the trace defect; no C98/C99 coverage guarantee.
  --ssgwas-apy-core-file <file>
                        Core-animal ID file (one ID per line, pedigree IDs); overrides
                        the seeded random core selection
  --ssgwas-apy-seed <s>  Seed of the deterministic random core selection (default 12345)
  --ssgwas-apy-core-target-ratio <x>
                        Harness-assigned label core_size / K98 (recorded in meta only; 0 = unset)
  --ssgwas-apy-core-variance <x>
                        Target proportion of genomic variance for automatic core sizing
                        (C98=0.98, C99=0.99; used when --ssgwas-apy-core-size auto)
  --ssgwas-apy-core-eigensolver <method>
                        Eigensolver for automatic core sizing: 'auto' (dense for small,
                        Lanczos for large), 'exact' (dense always), 'lanczos' (matrix-free)
  --ssgwas-grm-mode <mode>
                        VanRaden G construction: 'auto' (dense for small n_g, block-streamed
                        for large), 'dense' (full n_g^2 G), 'streamed' (block-streamed,
                        never materializes the full G; APY methods only)
  --ssgwas-apy-inverse-mode <dense|operator>
                        APY effect inverse application backend (default dense); operator
                        streams Gcn blocks and does not materialize Q in the operator
  --ssgwas-mme-apy-mode <dense|operator>
                        APY precision backend used by the MME (default dense); operator
                        requires --ssgwas-solver pcg and fixed variance components
  --ssgwas-pev-apy-mode <dense|operator>
                        APY Eq. 8 core PEV backend (default dense); operator streams
                        Gcn and avoids Q, explicit Gcc^-1, Kc^-1, and core Vhat
  --ssgwas-hinv-mode <dense|operator>
                        Pedigree A22^-1 backend (default dense); operator uses the
                        exact sparse Schur action and requires PCG + fixed VC
  --ssgwas-pev-scan-mode <scalar|batched>
                        SNP uncertainty scan backend (default batched; scalar is reference)
                        Exact/core-v2 use batched only when explicitly requested; omission
                        preserves historical scalar execution. Legacy hybrid stays scalar.
  --ssgwas-pev-scan-block-size <n>
                        Marker block width for batched uncertainty scan (default 256)
  --allow-large-exact-ssgwas-pvalue
                        Bypass the 8 GiB selected-inverse memory guard of
                        --ssgwas-pvalue-method exact (and apy-exact core slice)
  Method 2 outputs: <out>.ss_backsolve.ebv.txt / .snp.txt / .window.txt / .meta.txt
  With --ssgwas-pvalue-method exact additionally:
    <out>.ss_backsolve.pvalue.txt (EFFECT_ALLELE / SE_ALLELE / Z / P / LOG10P)

--- Pedigree Animal Model (B experimental) ---
  --animal-model        Pedigree-based Simple Animal Model GWAS: y = Xb + x*beta + Za + e,
                        a ~ N(0, A*sigma2_a) with full-pedigree A22 for the analysis sample set.
                        Requires --ped. A is built from the pedigree only (no genomic GRM);
                        LOCO does not apply (pedigree kernel is the same for all chromosomes).
                        Outputs *.assoc.txt (standard format) and *.animal_model.null.txt.
  --write-ped-kernel    Write the analysis A22 kernel (*.animal_model.k.txt + .k.id) for
                        external benchmarking (e.g. GEMMA -k).

--- Threshold/Probit Model (B experimental) ---
  --threshold-model      Binary trait GWAS (PQL + probit/logit)
  --threshold-dichotomize  Allow --threshold-model to dichotomise a non-binary
                         response at 0.5 (off by default; binary input is the contract)
  --link probit|logit    Link function (default: probit)
  PQL with fixed Ve=1 (probit) or Ve=pi^2/3 (logit).
  SPA is not a universal calibration guarantee; null-summary convergence flag is unreliable.

--- Count/Binomial PSR-GWAS Model (B experimental) ---
  --count-model          Count/binomial trait GWAS (PSR/PQL + weighted LMM)
  --family poisson|negbin|binomial Distribution (default: poisson)
  --offset <file>        One numeric log-offset/exposure value per sample for count traits
  --binomial-pheno <file> Binomial FID IID EVENT TRIAL file for event/trial traits

--- Ordinal PSR-GWAS Model (B experimental) ---
  --ordinal-model        Ordinal trait GWAS (cumulative-binary or stacked PSR)
  --ordinal-pheno <file> Ordinal FID IID ORDINAL phenotype file
  --ordinal-mode <name>  experimental-cumulative-binary (default) or stacked-psr
  --ordinal-stacked-engine <name> fixed-gls (default) or dense-zkz
  --ordinal-vc-method <name> fixed (default) or grid-reml
  --ordinal-test <name>  wald (default), score, robust-score, or all
  --ordinal-allow-fallback Fall back from stacked-psr to cumulative-binary
  --psr-null-only        Output null model (Vg, Ve, fixed effects, BLUP, quasi-likelihood)
                        then exit without SNP scan (available for threshold/count/ordinal)
  --psr-maxit <n>       Max PSR/PQL iterations (default: 50)
  --psr-tol <f>         PSR/PQL convergence tolerance (default: 1e-6)
  --psr-mu-eps <f>      Mu clipping epsilon for PSR (default: 1e-6)
  --psr-min-weight <f>  Minimum PSR working weight (default: 1e-8)
  --psr-max-weight <f>  Maximum PSR working weight (default: 1e8)

--- Survival prototype (C developer-only, not supported inference) ---
  --survival-model       Incomplete residual-like LMM prototype; no fitted Cox null
  --pheno format: FID IID TIME STATUS  (2 columns after FID/IID)
  runCoxNull only prints a message; do not use for screening or inference.

--- GxE Reaction Norm Model (B experimental) ---
  --gxe-model            SNP-by-environment interaction GWAS
  --env-covar <col>      Environment variable column in --pheno
  VCE: Vg0 (intercept) + Vg1 (slope on env) + residual.
  Scan: 2-df Wald test for beta_g + beta_gxe jointly.
  Output: BETA_G, SE_G, BETA_GxE, SE_GxE, JOINT_CHISQ, JOINT_P

=== INFRASTRUCTURE ===
  --make-grm             Calculate GRM and exit
  -gk [n]                GEMMA-compatible GRM alias; optional algorithm
  --make-grm-sparse      PLINK-only sparse VanRaden GRM and exit
  --sparse-cutoff <f>    Sparse GRM absolute-value cutoff (default: 0.05)
  --gs/-gs <n>           GEMMA-style S matrix for summary-stat VC (mode 1)
  --gq/-gq <n>           GEMMA-style q vector from --beta summary stats (mode 1)
  --gemma-eigen/-eigen   Calculate GRM eigen decomposition and write eigenD/eigenU
  --calccor/-calccor/-gw GEMMA-style SNP correlation windows; writes <out>.cor.txt
  --windowbp <n>         BP window for --calccor (default: 1000000; 0=no BP limit)
  --windowcm <f>         cM window for --calccor (default: 0; PLINK BIM only)
  --windowns <n>         SNP-count window for --calccor (default: 0=no SNP-count limit)
  --vc/-vc <n>           GEMMA-style VC estimate: 1 HE, 2 REML, 3 AI/IDUL REML
  --reml, --vce          Standalone production REML/VCE estimation
  --reml-method <m>      ai, em, brent, idul, idul_plus, or aidul
  --ci/-ci <n>           GEMMA-style summary-stat PVE confidence interval
  --beta/-beta <file>    Summary statistics for --gq or direct --vc with --cor
  --cor/-cor <file>      GEMMA .cor.txt for direct summary-stat --vc with --beta
  --cat/-cat <file>      SNP category file for summary-stat routes
  --mcat/-mcat <file>    List of SNP category files for summary-stat routes
  --wcat/-wcat <file>    SNP-by-category weights for summary-stat routes
  --wsnp/-wsnp <file>    SNP weights for summary-stat routes
  --study/-study <p>     Study prefix for summary-stat --vc (<p>.Vq/q/size.txt)
  --ref/-ref <p>         Reference prefix for summary-stat --vc (<p>.S/size.txt)
  --mstudy/-mstudy <f>   File listing study prefixes for summary-stat --vc/--ci
  --mref/-mref <f>       File listing reference prefixes for summary-stat --vc/--ci
  --make-grm-alg <n>     Method: 1=VanRaden, 2=Zeng, 3=Yang, 4=Vitezica,
                          5=TrueVar, 6=EMmvGWAS-compatible (default: 1)
                          Note: alg 6 is auto-set by --mv-compatible emmvgwas
  --mass-regress <file>  Multi-phenotype mass regressor (testing)
  --biobank-mode         B: native REGENIE-style two-step research path, not reference equivalence
  --regenie-model-snps <f> Model SNP list for Step 1
  --write-pred <file>    Write step 1 LOCO predictions
  --use-pred <file>      Use existing LOCO predictions
  --regenie-step <n>     0=both, 1=predictions only, 2=association only
  --predict-effect/-epm <f> GEMMA-style prediction from SNP effect table
  --predict-score <f>    Output prediction file (default: <out>.predict.txt)
  --predict-snp-col <c>  SNP ID column in effect table (default: SNP/rs/id)
  --predict-effect-col <c> Effect column (default: BETA/beta/alpha/effect)
  --bslmm/-bslmm <n>     B: native/compatibility modes: sparse ridge, ridge, probit, DAP-style
                          Native sampling is separate; ridge parity does not validate posterior parity.
  --bslmm-top <n>        Sparse top-SNP count for --bslmm 1/3 approximation (default: 100)
  --bslmm-h2 <f>         Ridge/GBLUP shrinkage h2 in (0,1), default 0.5
  --bslmm-mcmc           Native posterior-sampling BSLMM
  --bvsr                 Native posterior-sampling BVSR
  --mcmc-burnin <n>      MCMC burn-in iterations
  --mcmc-samples <n>     MCMC post-burn-in iterations
  --mcmc-thin <n>        MCMC thinning interval
  --seed <n>             Deterministic BSLMM/BVSR seed
  --bf/-bf <file>        DAP BF file for --bslmm 5
  --hyp/-hyp <file>      DAP hyperparameter file for --bslmm 5
  --bslmm-hgrid <n>      DAP h grid size for --bslmm 4 (default: 5)
  --bslmm-rho-grid <n>   DAP rho grid size for --bslmm 4 (default: 5)
  -hmin/-hmax <f>        GEMMA BSLMM/DAP h prior range
  -rmin/-rmax <f>        GEMMA BSLMM/DAP rho prior range
  -pmin/-pmax <f>        GEMMA BSLMM/DAP log10(pi) prior range
  -gmean <f>             GEMMA BSLMM gamma-count prior mean hint
  -hscale/-rscale/-pscale <f> GEMMA BSLMM proposal scale hints
  -mh <n>                GEMMA BSLMM MH step-count hint
  -smin/-smax <n>        GEMMA sparse gamma size bounds (-smax aliases --bslmm-top)
  --block-size <n>       Ridge block size (default: 1000)
  --skat                 Scoped set test with tail approximation, not exact-tail equivalence
  --burden               Run Burden rare-variant test (--set-file)
  --set-file <file>      Variant set definition for SKAT/Burden

Model notes:
  Models are mutually exclusive. Use exactly one --*-model flag.
  Default (no flag): standard univariate LMM (single-trait).
  Report bugs at: https://github.com/anomalyco/PanOmics/issues
```

## 附录B 常用命令行参数

本附录按任务列出常用参数。表中的 PREFIX 为文件前缀，FILE 为完整文件路径，N 为整数，F 为数值。输入格式和参数组合仍需满足对应模块的要求。更多参数可在终端运行 help models、help input 或 --help-all 查看。

### B.1 输入 输出与资源

| 参数 | 作用与用法 |
|---|---|
| --bfile PREFIX | BED BIM FAM 前缀，不含后缀 |
| --pfile PREFIX | PGEN PVAR PSAM 前缀 |
| --bgen FILE --sample FILE | BGEN 和配套样本文件 |
| --pheno FILE | 表型文件，前两列 FID IID |
| --mpheno N | 两个 ID 后从1计数的表型列 |
| --mpheno-list 1,2,3 | 多性状列，按给定顺序输出效应 |
| --tpheno N | 两个 ID 后从1计数的时间列 |
| --qcovar FILE | 连续协变量文件 |
| --covar FILE | 分类协变量文件 |
| --pedigree FILE | 系谱 ID SIRE DAM |
| --ped-genotypes FILE --map FILE | PED MAP 基因型输入 |
| --threads N | 计算线程数；未指定时默认使用可用最大线程 |
| --memory-gib F | 内存预算诊断值，不限制进程实际占用 |
| --out PREFIX | 输出前缀，先创建父目录 |

### B.2 筛选与关系矩阵

| 参数 | 作用与默认值 |
|---|---|
| --maf F | 最小 MAF，默认0.01 |
| --hwe F | HWE P 值筛选，默认0为关闭 |
| --geno F --mind F | BED 位点与个体最大缺失比例，范围0至1 |
| --keep FILE --remove FILE | BED 样本 ID 列表；remove 优先 |
| --extract FILE --exclude FILE | 纳入或排除位点 ID 列表 |
| --chr CHR | 指定扫描染色体 |
| --from-bp N --to-bp N | 指定扫描物理位置范围 |
| --make-grm-alg N | 默认1；1 VanRaden，2 Zeng，3 Yang，4 Vitezica，5 TrueVar，6 EMmvGWAS 兼容算法 |
| --make-grm-sparse | BED 稀疏 VanRaden GRM |
| --sparse-cutoff F | 稀疏矩阵截断，默认0.05 |

### B.3 模型 效应与方差估计

| 参数 | 作用与选择 |
|---|---|
| --model NAME | lmm、lm、multivariate、repeatability、rrm、animal、maternal、ssgwas、common-env、threshold、count、ordinal、gxe |
| --gwas-model NAME | additive、additive_dominance、dominance；默认 additive |
| --snp-model NAME | additive、dominant、recessive、heterozygote；非加性编码需硬调用；显性结构需 additive 位点编码 |
| --no-loco | 关闭 LOCO，使用全基因组模型 |
| --scan-tier NAME | approx 固定零模型；semi 固定比例并处理尺度；exact 逐位点重估；按模块满足使用条件 |
| --lmm-null-solver NAME | brent、idul、idul_plus、aidul；默认 brent |
| --reml-maxit N | 相应 VCE 最大迭代次数，默认50 |
| --reml-tol F | 相应 VCE 收敛容差，默认1e-6 |
| --reml-vg F --reml-vpe F --reml-ve F | 对应模型的固定遗传、永久环境和残差方差 |
| --lm-mode N | 1 Wald，2 LRT，3 score，4全部 |
| --lmm-gemma-mode N | 对应 LMM 检验选择：1 Wald，2 LRT，3 score，4全部 |
| --lmm-gemma-refit | 对应检验逐位点重估 |

### B.4 多性状 重复测量与随机回归

| 参数 | 作用与选择 |
|---|---|
| --mv-method NAME | 选择多性状计算方法，按第5.9节和程序帮助设置 |
| --mv-test LIST | 对应多性状方法的 wald、lrt、score |
| --mv-compatible emmvgwas | 选择 EMmvGWAS 兼容模型；semi 需要该设置 |
| --mv-vce-solver NAME | em 或 aireml，默认 em |
| --mv-output-vbeta | 增加效应协方差输出 |
| --repeatability-mode NAME | auto、exact、fast、stcg；默认 auto；与 scan-tier 分开设置 |
| --repeatability-grm-cache PREFIX | 保存或读取匹配数据和设置的个体 GRM 缓存 |
| --repeatability-weight FILE | FID IID WEIGHT 或 FID IID TIME WEIGHT |
| --rrm-order N | 基函数阶数参数，范围1至4，定义为对应基函数自由度减1 |
| --rrm-basis NAME | legendre 或 bspline，默认 legendre |
| --rrm-unbalance | 非平衡长表流程 |
| --rrm-gwas-method NAME | record 或 coefficient，默认 coefficient |
| --rrm-curve-grid N | 曲线时间网格点数，默认21 |
| --rrm-top N | 曲线输出位点数量，默认100 |

### B.5 单步与其他表型模型

| 参数 | 作用与选择 |
|---|---|
| --single-step-method NAME | lmm 为 Method 1；backsolve 为 Method 2 |
| --genotyped-list FILE | 系谱中的基因型个体 ID 列表 |
| --ssgwas-solver NAME | ldlt 或 pcg，默认 ldlt |
| --ssgwas-pvalue-method off | Method 2 默认不计算位点 P 值 |
| --ssgwas-pvalue-method exact | 在固定方差等对应条件下，通过 LDLT 选定逆计算位点不确定性 |
| --ssgwas-window-snps N | Method 2 每个窗口位点数，默认10 |
| --logistic NAME | score 或 wald；表型按0/1编码 |
| --firth-logistic | Firth 校正路径，检查估计与收敛 |
| --family NAME | count 模型的 poisson、negbin 或 binomial；默认 poisson |
| --offset FILE | count 模型逐样本数值对数暴露量 |
| --binomial-pheno FILE | FID IID EVENT TRIAL 二项表型 |
| --ordinal-pheno FILE | FID IID ORDINAL 等级表型 |
| --link NAME | threshold 的 probit 或 logit；默认 probit |
| --env-covar COL | 对应 GxE 模型的表型环境列 |

### B.6 集合与后关联分析

| 参数 | 作用与默认值 |
|---|---|
| --burden --set-file FILE | Burden 集合检验 |
| --skat --set-file FILE | SKAT 集合检验，尾概率使用相应近似 |
| --rare-max-maf F | 集合位点最大 MAF，默认0.5；稀有变异研究应显式设置 |
| --rare-min-mac N | 集合位点最小 MAC，默认1 |
| --ld-prune --ld-window N --ld-step N --ld-r2 F | BED LD 剪枝；窗口默认250000 bp，步长1，相关阈值0.2 |
| --ld-clump FILE --clump-p F | 使用 SNP P CHR BP 统计表聚类；P 默认5e-8 |
| --condition-snp FILE | 条件位点，每行一个 SNP ID |
| --cojo-summary FILE --cojo-ld-ref FILE | COJO 汇总统计与 signed-r LD |
| --cojo-p F --cojo-r2 F --cojo-wind N | 默认5e-8、0.9、10000000 bp |
| --meta-studies FILE --meta-method NAME | 研究文件列表；fixed、random-dl 或 random-reml，默认 fixed |
| --meta-min-studies N | 最少研究数，默认2 |
| --finemap-summary FILE --finemap-region REGION | 汇总统计及 chr:start-end 区域；需要实际基因型 LD |
| --credible-level F | 可信集概率水平，默认0.95 |
| --finemap-prior-var F | 精细定位效应先验方差，默认0.04 |
| --finemap-max-signals N | 最多条件信号数，默认5 |

### B.7 输入整理与运行帮助

| 参数或命令 | 作用 |
|---|---|
| --input-file FILE | 原始 CSV TSV TXT 表格 |
| --input-config FILE | 输入整理 JSON 配置 |
| --write-template-config FILE | 写出可编辑配置 |
| --preview-schema --infer-columns | 预览字段并辅助识别 |
| --id-col NAME --iid-col NAME | 原始表格个体字段 |
| --time-col NAME --trait-cols LIST | 原始表格时间与性状字段 |
| --prepare-input-only | 只整理输入，先检查结果再分析 |
| --prepared-input-dir DIR | 整理文件输出目录 |
| --version | 查看程序版本 |
| --help | 查看常用帮助 |
| help models | 查看模型帮助 |
| help input | 查看输入帮助 |
| --help-all | 查看当前程序完整参数 |

## 附录C 参数迁移与分析记录

### C.1 常用入口迁移

| 原有入口 | 推荐写法 | 保留事项 |
|---|---|---|
| --lm | gwas --model lm | 保留输入与检验模式 |
| --multivariate | gwas --model multivariate | 保留性状列、方法与扫描档位 |
| --repeatability-model | gwas --model repeatability | 保留时间列、记录策略和方差 |
| --single-step | gwas --model ssgwas | 显式核对 Method 1 或 backsolve |
| --make-grm | grm | 保留矩阵算法与筛选 |
| --ped 用于系谱 | --pedigree FILE | 不与基因型 PED 混用 |
| --ped 用于基因型 | --ped-genotypes FILE --map FILE | 不与系谱混用 |

原有平铺命令仍可使用。迁移入口时保留模型定义、表型列、LOCO 与检验策略；按当前输出表头更新结果读取脚本，并分别保存升级前后的程序身份与结果。

### C.2 推荐分析记录模板

```text
研究与分析编号
程序版本及 PACKAGE_ID
下载包 SHA-256
操作系统与硬件资源
基因型与表型文件及 SHA-256
样本匹配前后数量
位点筛选前后数量
模型与效应结构
LOCO 与 scan-tier
协变量与时间编码
方差估计或固定值及依据
完整运行命令
运行退出状态与收敛诊断
主要输出文件
多重检验及绘图设置
近似参数 种子与回退记录
```
