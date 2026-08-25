# RRM 复用 Prepared Input 设计

## 背景

Repeatability GWAS 和 RRM 都需要纵向记录结构：

```text
FID IID TIME TRAIT VALUE
```

当前输入适配层已经为 repeatability 生成 `pheno_long.tsv`、`traits/<trait>.pheno.tsv`、`model_qcovar.tsv`、`model_covar.tsv`。RRM 后续应复用同一 prepared input，而不是重新实现一套输入逻辑。

## 目标

下一阶段支持：

```bash
cosmicscanner \
  --bfile data/genotype \
  --prepared-input-dir prep/semen \
  --input-model rrm \
  --input-trait volume \
  --rrm-order 3 \
  --gwas \
  --out result/volume_rrm
```

或一步到位：

```bash
cosmicscanner \
  --bfile data/genotype \
  --input-file data/semen.csv \
  --input-model rrm \
  --id-col pig_id \
  --time-col collect_date \
  --trait-cols volume,concentration \
  --input-trait volume \
  --rrm-order 3 \
  --gwas \
  --out result/volume_rrm
```

## 复用文件

RRM 可直接复用：

```text
traits/<trait>.pheno.tsv
record_manifest.tsv
pheno_long.tsv
covar_record_q.tsv
covar_record_c.tsv
input_audit.json
```

首版 RRM prepared 接入建议使用：

```text
traits/<trait>.pheno.tsv
```

格式：

```text
FID IID TIME VALUE
```

这样和当前 repeatability 文件一致，也能最小化 RRM reader 改动。

## TIME 语义

RRM 比 repeatability 更依赖 TIME 的连续尺度。需要区分：

```text
TIME
TIME_BASE
WITHIN_TIME_INDEX
RAW_TIME
```

当前 repeatability 为保证记录键唯一，默认：

```text
TIME = TIME_BASE * multiplier + WITHIN_TIME_INDEX
```

RRM 中如果把 `TIME` 当连续轨迹变量，这种唯一化会轻微改变时间尺度。因此 RRM 接入必须明确时间策略。

建议：

- repeatability 默认使用 `TIME`
- RRM 默认使用 `TIME_BASE`
- 如果 RRM 需要唯一记录键，内部用 `RECORD_ID` 或 `(FID,IID,TIME,WITHIN_TIME_INDEX)`，不要把 `WITHIN_TIME_INDEX` 混入轨迹时间

这意味着 RRM prepared 接入不能只盲读 `traits/<trait>.pheno.tsv` 中的 `TIME`。更稳的方式是通过 `record_manifest.tsv` 把 `TIME` 映射回 `TIME_BASE`。

## 建议实现路径

### 阶段 1：显式 RRM prepared reader

新增一个轻量转换步骤：

```text
prepared input
-> traits/<trait>.pheno.tsv
-> record_manifest.tsv
-> rrm_trait.pheno.tsv using TIME_BASE
-> current RRM path
```

输出临时/派生文件：

```text
rrm/<trait>.pheno.tsv
```

格式：

```text
FID IID TIME_BASE VALUE
```

再把它作为当前 RRM 的 long phenotype 输入。

### 阶段 2：record covariate 接入

当前 RRM 支持：

```text
--tclass <cols>
--tcovar <cols>
--tpoint <time_col>
```

输入适配层应把：

```text
covar_record_q.tsv
covar_record_c.tsv
```

转换为 RRM 可消费的 record-level covariate 表，并保留列名。

### 阶段 3：统一 PreparedInputRunner

最终应抽象：

```text
PreparedInputRunner::applyToRepeatability(cfg, prepared_dir, trait)
PreparedInputRunner::applyToRrm(cfg, prepared_dir, trait)
```

避免在 `gwas_app.cpp` 里堆模型特定逻辑。

## CLI 行为建议

`--input-model repeatability`：

- 生成 prepared input
- 默认选 `traits/<input_trait>.pheno.tsv`
- `TIME` 使用唯一化后的模型时间

`--input-model rrm`：

- 生成同样 prepared input
- 派生 `rrm/<input_trait>.pheno.tsv`
- RRM 轨迹时间默认使用 `TIME_BASE`
- `record_manifest.tsv` 作为时间追溯源

## 风险

主要风险：

- RRM TIME 语义被唯一化污染
- 同一天多条记录如何进入轨迹模型需要统计解释
- record-level categorical covariate 与 RRM 现有 `tclass` 接口要对齐
- 宽表时间列解析可能与 RRM time grid 不一致

## 验收标准

RRM 接入完成时应满足：

- `--prepared-input-dir prep --input-model rrm --input-trait trait` 可运行
- RRM 使用 `TIME_BASE` 而不是唯一化 `TIME`
- `record_manifest.tsv` 中的 `WITHIN_TIME_INDEX` 不丢失
- small smoke 输出 RRM assoc 文件
- 与手工准备 long pheno 的 RRM 结果一致或可解释一致
