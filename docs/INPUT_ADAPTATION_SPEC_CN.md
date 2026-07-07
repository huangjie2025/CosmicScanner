# CosmicScanner 输入适配层规范

## 目标

输入适配层负责把真实业务表格转换为 CosmicScanner 内部稳定消费的 prepared input。模型层不直接读取原始 CSV/TSV/TXT，而是读取输入适配层生成的标准文件。

当前规范覆盖 repeatability GWAS、RRM、standard LMM 和 multivariate LMM 的 prepared input 首版接入。

## 当前支持范围

支持：

- 输入格式：`csv`、`tsv`、`txt`、`auto`
- 输入模型：`repeatability`、`rrm`、`lmm`、`multivariate`
- 显式列映射
- JSON 配置文件
- schema preview
- column inference preview
- long table
- 简单 wide-to-long
- genotype FID/IID 匹配
- date/integer/batch 时间解析
- 同 ID 同 TIME 多记录唯一化
- record-level 和 subject-level covariate 分流
- prepared input 复用
- 真实公猪精液数据回归测试

暂不支持：

- Excel 直接读取
- YAML 配置
- 自动确认列映射并直接运行
- BLUPF90/ASReml/GCTA/GEMMA profile
- threshold/count/survival 等非 LMM 模型公共接入

## CLI 参数

核心参数：

```text
--input-file <file>
--input-config <json>
--write-template-config <json>
--input-format auto|csv|tsv|txt
--input-model repeatability|rrm|lmm|multivariate
--prepared-input-dir <dir>
--prepare-input-only
```

列映射：

```text
--id-col <name>
--iid-col <name>
--time-col <name>
--trait-cols <comma-list>
--input-trait <name>
--record-qcovars <comma-list>
--record-ccovars <comma-list>
--subject-qcovars <comma-list>
--subject-ccovars <comma-list>
```

时间与结构：

```text
--time-parser auto|date|integer|batch
--time-unique-policy strict|within-time-index
--time-index-multiplier auto|<n>
--wide-to-long
```

检查与推断：

```text
--preview-schema
--infer-columns
```

## JSON 配置格式

配置文件使用轻量 JSON，不依赖第三方 YAML 库。CLI 参数优先级高于配置文件。

示例：

```json
{
  "input_file": "data/semen.csv",
  "input_format": "csv",
  "input_model": "repeatability",
  "prepared_input_dir": "prep/semen",
  "id_col": "pig_id",
  "iid_col": "",
  "time_col": "collect_date",
  "input_trait": "volume",
  "trait_cols": ["volume", "concentration", "motility"],
  "record_qcovars": ["age", "temp"],
  "record_ccovars": [],
  "subject_qcovars": [],
  "subject_ccovars": ["farm", "line"],
  "time_parser": "date",
  "time_unique_policy": "within-time-index",
  "time_index_multiplier": "auto",
  "wide_to_long": false,
  "infer_columns": false
}
```

模板生成：

```bash
cosmicscanner --input-file data/semen.csv --write-template-config input_prepare.json
```

## Prepared Input 文件

标准输出目录包含：

```text
sample_manifest.tsv
record_manifest.tsv
pheno_long.tsv
traits/<trait>.pheno.tsv
rrm/<trait>.pheno.tsv
subject/<trait>.pheno.tsv
subject/multivariate.pheno.tsv
covar_record_q.tsv
covar_record_c.tsv
covar_subject_q.tsv
covar_subject_c.tsv
model_qcovar.tsv
model_covar.tsv
model_subject_qcovar.tsv
model_subject_covar.tsv
id_match_report.tsv
input_audit.json
summary.tsv
trait_manifest.tsv
covariate_audit.tsv
subject_covar_check.tsv
unresolved_columns.tsv
schema_preview.tsv
column_inference.tsv
inferred_input_config.json
validation_report.tsv
```

其中 `schema_preview.tsv` 只在 `--preview-schema` 下生成，`column_inference.tsv` 只在 `--infer-columns` 下生成。

### `pheno_long.tsv`

```text
FID IID TIME TRAIT VALUE
```

这是输入层的标准中间表型格式。

### `traits/<trait>.pheno.tsv`

```text
FID IID TIME VALUE
```

这是当前 repeatability reader 直接消费的性状拆分文件。

### `rrm/<trait>.pheno.tsv`

```text
FID IID TIME VALUE
```

这是当前 RRM reader 直接消费的性状拆分文件。这里的 `TIME` 使用 `record_manifest.tsv` 中的 `TIME_BASE`，不是 repeatability 用于唯一记录键的 `TIME`。

### `record_manifest.tsv`

```text
FID IID TIME RECORD_ID RAW_TIME TIME_BASE WITHIN_TIME_INDEX RAW_ROW_INDEX
```

用于追溯原始行、原始时间、标准化时间和同时间内序号。

### `model_qcovar.tsv` / `model_covar.tsv`

```text
FID IID TIME covariate...
```

这是当前 repeatability 模型直接消费的协变量文件。record-level 和 subject-level covariate 会在这里统一展开到记录级。

### `subject/<trait>.pheno.tsv`

```text
FID IID VALUE
```

这是 standard LMM 直接消费的 subject-level 性状文件。如果原始输入包含多记录，当前首版按 `FID/IID/trait` 取非缺失均值。

### `subject/multivariate.pheno.tsv`

```text
FID IID trait1 trait2 ...
```

这是 multivariate LMM 直接消费的 subject-level 多性状文件。每个性状按个体取均值。

### `model_subject_qcovar.tsv` / `model_subject_covar.tsv`

```text
FID IID covariate...
```

这是 standard LMM 和 multivariate LMM 直接消费的无表头 subject-level 协变量文件。

## 时间规则

`--time-parser`：

- `date`：将日期转为 `days_since_min_date`
- `integer`：直接读取数值时间
- `batch`：按出现顺序编号
- `auto`：全数值时按 integer，否则优先按 date

`--time-unique-policy`：

- `strict`：同一 `FID/IID/TIME_BASE` 出现多条记录时报 `DUPLICATE_RECORD_KEY`
- `within-time-index`：使用 `TIME = TIME_BASE * multiplier + WITHIN_TIME_INDEX`

`--time-index-multiplier auto` 会根据最大同 ID 同时间重复数自动选择 10、100、1000 等安全倍率。

## Wide-To-Long 规则

`--wide-to-long` 当前支持：

```text
trait_1
trait_d1
```

例子：

```text
volume_d1 volume_d2 concentration_1 concentration_2
```

会生成：

```text
TRAIT=volume, TIME=1/2
TRAIT=concentration, TIME=1/2
```

无法识别的列会进入 `unresolved_columns.tsv`。

增强 wide-to-long 还支持：

```text
day1_trait
d1_trait
trait_day1
trait_T1
trait_20240101
```

## Column Inference

`--infer-columns` 只输出建议，不自动替代显式列声明。

输出文件：

```text
column_inference.tsv
inferred_input_config.json
```

列包括：

```text
column suggested_role confidence reason nonmissing missing_rate unique_values numeric_fraction id_match_rate
```

推断依据：

- 列名关键词
- genotype ID 匹配率
- 日期可解析比例
- 数值比例
- 唯一值数量

## Audit 与报告

`input_audit.json` 记录：

- 输入文件
- 列映射
- 时间策略
- `time_index_multiplier`
- `wide_to_long`
- `infer_columns`
- 匹配记录数
- 未匹配记录数
- duplicate ID-time 键数量
- mixed subject covariate warning 数量
- 每个 trait 非缺失数
- 报告文件列表

`summary.tsv` 记录：

- matched IDs
- matched rows
- repeat 分布
- duplicate time keys
- time multiplier
- trait 非缺失数

`covariate_audit.tsv` 记录：

- covariate 名称
- record/subject 层级
- quantitative/categorical 类型
- 非缺失数
- 缺失率
- 唯一值数量

`subject_covar_check.tsv` 记录 subject-level covariate 在同一 IID 内是否出现多个值。

`validation_report.tsv` 记录：

- `ID_MATCH_RATE_LOW`
- `TRAIT_ALL_MISSING`
- `COVARIATE_HIGH_MISSING`
- `CATEGORICAL_TOO_SPARSE`
- `MIXED_SUBJECT_COVARIATE`

## 错误码

输入适配层错误信息必须包含稳定错误码。当前规范错误码：

```text
UNSUPPORTED_INPUT_FORMAT
UNSUPPORTED_INPUT_MODEL
INPUT_FILE_NOT_READABLE
INPUT_FILE_EMPTY
DUPLICATE_INPUT_COLUMN
MISSING_REQUIRED_COLUMN
TIME_NOT_PARSEABLE
DUPLICATE_RECORD_KEY
INVALID_TIME_MULTIPLIER
WIDE_COLUMN_PARSE_FAILED
PREPARED_PHENOTYPE_NOT_FOUND
INPUT_CONFIG_NOT_READABLE
TEMPLATE_CONFIG_NOT_WRITABLE
OUTPUT_DIR_NOT_WRITABLE
TRAIT_NOT_PREPARED
```

warning 类报告码：

```text
ID_NOT_MATCHED
MIXED_SUBJECT_COVARIATE
UNMAPPED_COLUMN
```

## Repeatability 对接

`--input-file` 流程：

```text
raw csv/tsv/txt
-> InputPreparer
-> prepared input directory
-> cfg.pheno_file = traits/<selected>.pheno.tsv
-> cfg.qcovar_file = model_qcovar.tsv
-> cfg.covar_file = model_covar.tsv
-> RepeatabilityModel
```

`--prepared-input-dir` 流程：

```text
prepared input directory
-> traits/<input_trait>.pheno.tsv
-> model_qcovar.tsv / model_covar.tsv if present
-> RepeatabilityModel
```

## RRM 对接

`--input-model rrm` 流程：

```text
raw csv/tsv/txt 或 prepared input directory
-> rrm/<selected>.pheno.tsv
-> cfg.long_pheno_file
-> RrmModel
```

RRM 接入规则：

- 使用 `TIME_BASE` 作为轨迹时间
- 不使用 repeatability 唯一化后的 `TIME`
- 当前首版不自动传入 record-level covariate
- `record_qcovars` 自动作为 `tcovar`
- `record_ccovars` 自动作为 `tclass`
- 后续仍需增强更复杂的 RRM covariate 编码和诊断

## Standard LMM 对接

`--input-model lmm` 流程：

```text
raw csv/tsv/txt 或 prepared input directory
-> subject/<selected>.pheno.tsv
-> model_subject_qcovar.tsv / model_subject_covar.tsv
-> UnivariateModel
```

## Multivariate LMM 对接

`--input-model multivariate` 流程：

```text
raw csv/tsv/txt 或 prepared input directory
-> subject/multivariate.pheno.tsv
-> model_subject_qcovar.tsv / model_subject_covar.tsv
-> MultivariateModel
```

## Excel 策略

当前 C++ 输入层仍安全拒绝 `xlsx/xls`，错误码为：

```text
UNSUPPORTED_INPUT_FORMAT
```

原因是直接引入 Excel 解析库会影响 release 构建和依赖可控性。当前正式建议是先另存为 `csv/tsv` 后使用输入适配层。后续 Excel 支持应做成可选依赖或官方转换器。

## 验证

专项测试：

```bash
bash tests/run_test_input_adaptation.sh
```

覆盖：

- long CSV
- JSON config
- template config
- column inference
- schema preview
- audit reports
- strict duplicate failure
- wide-to-long
- 公猪精液真实数据 prepared summary 回归
- 公猪精液 prepared-to-GWAS 小规模 smoke
- RRM prepared input 小规模 smoke
- RRM record-level covariate smoke
- Standard LMM prepared input smoke
- Multivariate LMM prepared input smoke
