# CosmicScanner 公共环境模型 GWAS 开发总结

> 实验性模型，开发分支实现；不在 CosmicGWAS v1.1.0 focused 发布范围内。

## 一、统计模型与原理

同一受试对象存在多条记录、且记录来自不同环境（组）时，公共环境模型在重复力模型基础上把"受试对象永久环境效应"替换为"组共享环境效应"：

```
y_r = X_r b + x_jr·beta_j + (Z_a a)_r + (Z_c c)_r + e_r

a ~ N(0, K·sigma2_a)      # 跨受试对象的基因组 GRM 遗传效应
c ~ N(0, I·sigma2_c)      # 按组共享的随机环境截距（组间独立）
e ~ N(0, I·sigma2_e)
```

- `Z_a`：记录→受试对象映射（subject rank）；`K` 为 subjects×subjects GRM（默认 VanRaden）。
- `Z_c`：记录→组映射；同组记录共享 `sigma2_c` 相关。组可以由受试对象文件定义（subject 级，该受试对象全部记录进同一组），也可以按记录文件逐行定义（record 级）。
- 记录层协方差 `V = vg·Z_a·K·Z_a' + vc·Z_c·Z_c' + ve·I`（n_records×n_records），由通用 multi-random-effect LMM 层做 LDLT + jitter 因子化。
- 方差分量 `(sigma2_a, sigma2_c, sigma2_e)` 默认 REML 估计；`--reml-vg/-vc/-ve` 可固定后跳过 VCE。

## 二、实现流程路线

1. **解析与分组**：长格式表型 `FID IID [TIME] VALUE`（`--tpheno <col>` 指定值列）按 `FID IID` 与基因型 `GenRecord` 对齐；subject 级分组按 `FID IID GROUP` 读入、以 GROUP LABEL 建 level（同一 label 的所有受试对象共享一个 level）；record 级分组按行顺序 `FID IID [TIME] GROUP` 对齐到记录，错行记警告并剔除。
2. **GRM**：`GrmBuilder::compute` 构建 subjects GRM（PLINK/BGEN/PGEN 三条路径），扫描前 `subsetByRank` 按保留 subject 排序。
3. **VCE**：记录层 REML（3 分量），`runVCE` 调用通用 LMM 层；失败/异常时打印诊断。
4. **扫描**：块式（block_size=256）统一管线——PLINK 走 `getSnpGenotypes`、BGEN/PGEN 走 `readVariant` 顺序读取 → 按记录映射到 dosage → 均值中心化 → GLS score 检验（`scores = G'Py`、`var_score = gVg - xvX·inv(XtVinvX)·xvX'`）→ Wald p 值（`erfc(|z|/√2)`）→ `appendAssocResult` 追加输出。
5. **LOCO**：`computeSplit` 得到全局与分染色体标准化基因型交叉积；`G_loco = G*(S_global/S_loco) - raw_chr/S_loco`（`raw_chr = Z_chr·Z_chr'` 未归一化，S 为 snp 数或 sum_2pq），逐染色体重估 VCE 并扫描。
6. **输出**：`*.assoc.txt`（`CHR SNP BP A1 A2 FREQ N BETA SE P LOG10P`）与 `*.common_env.null.txt`。

## 三、代码文件清单

| 文件 | 说明 |
|------|------|
| `src/models/common_environment_model.h/.cpp` | 模型入口：分组解析、prepare、VCE、扫描、输出（新增） |
| `src/models/multi_random_effect_lmm.h/.cpp` | 通用 multi-random-effect LMM 层（新增） |
| `src/gwas_app.cpp` | 模型互斥检查 + dispatch（修改） |
| `src/main.cpp` | CLI 参数与 help（修改） |
| `src/config.h` | `common_env_model/common_env_file/common_env_record_level/reml_vc` 字段（修改） |
| `tests/run_test_common_environment.sh` | 模拟验证脚本（新增） |

## 四、CLI 参数完整清单

```
--common-env-model            启用公共环境模型（隐式 --gwas）
--common-env-file <f>         分组文件（默认 subject 级：FID IID GROUP）
--common-env-record-level     分组文件按记录解析（FID IID [TIME] GROUP）
--reml-vg/-vc/-ve <v>         固定方差分量，跳过 VCE
--no-loco                     关闭逐染色体 LOCO
--min-maf / --min-hwe         SNP 过滤（复用单性状路径语义）
--bfile / --bgen / --pgen     基因型输入
--covar / --qcovar            记录级协变量
```

## 五、测试数据集与验证方法

`tests/run_test_common_environment.sh` 自包含模拟：N=120、M=160（chr1=100, chr2=60）、12 组×10 个体、360 条记录（2+i%3/受试对象）、真值 `sigma2a/c/e = 0.45/0.25/0.30`、因果 SNP rs8（β=0.7）、seed 20240815、GRM 加 1e-6 ridge 保证正定。

四组运行：

| 运行 | 参数 | 目的 |
|------|------|------|
| ce_global | `--no-loco` | VCE + 全局扫描基线 |
| ce_loco | 默认 | 逐染色体 LOCO |
| ce_record_level | `--common-env-record-level --no-loco` | 记录级分组路径 |
| ce_fixedvc | `--reml-vg 0.45 --reml-vc 0.25 --reml-ve 0.30 --no-loco` | 固定 VC 路径 |

## 六、验证结果

| 检查项 | 结果 |
|--------|------|
| VCE 回收（global） | 0.4877 / 0.3475 / 0.3227（真值 0.45/0.25/0.30，误差 0.0377/0.0975/0.0227） |
| 因果 SNP（global） | top hit，beta=-0.642, p=6.5e-07 |
| 因果 SNP（LOCO） | p=5.7e-10 |
| 因果 SNP（record-level） | p=6.5e-07（与 subject 级一致） |
| 关联行数 | global/LOCO/record-level/fixed-VC 均为 160/160 |
| 固定 VC | `--reml-vg/-vc/-ve` 精确生效 |
| 分组语义 | 120 受试对象 → 12 组 level（按 label 而非 subject 建组） |
| 回归 | `run_test_repeatability.sh`、`run_rrm_dynamic_simulation.sh` 通过 |

## 七、开发中的关键修正

1. **score 检验矩阵方向错误**（导致 global 缺 9 个 SNP、LOCO 全部被丢弃）：`term2` 应为 `xvX · (X'V^-1X)^-1 · xvX'` 的逐 SNP 对角线，写成 `(X'V^-1X)^-1 · xvX` 属于维数不匹配的未定义行为。修正为 `xvX * XtVinvXInv` + `rowwise().sum()`，与 `repeatability_scan_exact.cpp` 的 `xvg'·XtVinvX_inv·xvg` 一致。
2. **LOCO 公式错误**：`G_chr` 是未归一化交叉积 `Z_chr·Z_chr'`，不能乘 `w2 = S_chr/S_loco`；正确为 `G_loco = G_used*(S_global/S_loco) - G_chr/S_loco`，与 `gwas_scanner.cpp` PLINK-style LOCO 一致。修正后 LOCO 关联行数与因果 p 恢复正常。
3. **subject 级分组按 subject 建 level**（120 组而非 12 组）：改为按 GROUP LABEL 建 level。
4. **编译期 API 对齐**：`RandomComponent` 需 `#include "mme_builder.h"`；`factorize_with_jitter` 只接受 `Eigen::LDLT<MatrixXd>&`（`XtVinvX_llt_` 由 LLT 改 LDLT）；GRM 适配器改 `unique_ptr<AbstractMatrix>`。

## 八、输出文件格式

`*.assoc.txt`（与单性状 LMM 同构，`appendAssocResult`）：

```
CHR  SNP  BP  A1  A2  FREQ  N  BETA  SE  P  LOG10P
```

`*.common_env.null.txt`：

```
Item  Value
n_records  ...
n_subjects  ...
n_groups  ...
sigma2_a  ...
sigma2_common  ...
sigma2_e  ...
h2  ...
common_environment_fraction  ...
converged  ...
vce_best_logL  ...
solver  record_level_dense_ldlt
scan_mode  record_level_gls | record_level_gls_loco
```

## 九、当前限制与风险

- **实验性模型**：未经过与外部软件（如 h2-GREML/GEMMA 对照）的正式对标，建议以模拟验证口径为准。
- 记录层密集 LDLT：内存随记录数 O(n²)，大样本需 Woodbury/fast 化路线。
- LOCO 每染色体重估 VCE 耗时随染色体数线性增长；BGEN/PGEN 的 LOCO 依赖按染色体逐 SNP 顺序读取。
- 分组文件无 GROUP 标签自动推断（缺失行/多余行剔除并告警）；record 级模式依赖表型记录顺序一致。
- 只支持连续表型（无 GLMM/序数扩展）。

## 十、验证脚本与结果索引

- 脚本：`tests/run_test_common_environment.sh`（`--build` 或复用现有 `build/cosmicscanner`）
- 数据：`tests/data/common_env_simulation/`（PLINK bed/bim/fam、长格式表型、subject/record 分组文件、truth.npz）
- 结果：`tests/results/common_env_simulation/`（ce_global/ce_loco/ce_record_level/ce_fixedvc 的 assoc/null/log）
- 回归：`bash CosmicScanner/tests/run_test_repeatability.sh`、`bash CosmicScanner/tests/run_rrm_dynamic_simulation.sh`
