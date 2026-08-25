# [OPEN] repeatability-woodbury

## 问题概述

- 症状：公猪精液真实数据 `vol` 性状在 `repeatability-model + exact + --no-loco` 路线下失败。
- 已知现象：
  - 输入读取成功
  - GRM 构建完成
  - subject-level VCE 三组初始化均返回 `logL=-inf`
  - exact scan 阶段报错 `failed to factorize Woodbury subject matrix`
  - `assoc` 文件只有表头，无有效结果

## 当前结论状态

- 状态：根因已定位，修复已完成首轮验证
- 业务逻辑修改：已实施最小修复
- 首要目标：继续做大样本真实数据验证

## 可证伪假设

1. 记录级协变量与表型记录未完全对齐，导致设计矩阵中存在大量全 NA/异常行，进而使 VCE 初始化退化。
2. 分类协变量哑变量展开后出现严重共线或近共线，导致 `X'V^-1X` 或 Woodbury 主矩阵不可逆/不可稳定分解。
3. 某些个体的重复记录数、权重或残差结构在真实数据上触发极端不平衡，导致 subject-level VCE 初始化直接落入 `logL=-inf`。
4. `exact` 路线在大规模真实数据下对 `vg/vpe/ve` 的初始化或边界约束不稳，导致方差分量全部塌缩到 0。
5. Woodbury 求解中实际失败的不是整体模型，而是某个中间矩阵出现非有限值、负对角或极端条件数。

## 证据计划

- 读取 `repeatability_model.cpp` 中 VCE 初始化与 Woodbury 分解路径
- 对关键矩阵维度、有限值比例、对角范围、最小/最大特征信息做插桩
- 复现 `vol` smoke test，先定位失败点，再决定最小修复

## 记录

- 2026-06-12：建立调试会话，准备进入插桩阶段
- 2026-06-12：运行时证据确认两层问题：
  - 预处理把 `TIME` 简化为 `days_since_min_date`，导致同一公猪同一天多次记录共享同一个 `FID/IID/TIME` 键；真实数据中存在最多 6 次同日重复。
  - 修正 `TIME` 唯一化后，`repeatability_model.cpp` 中固定效应设计矩阵仍然秩亏；插桩显示 `x_cols=34, x_rank=31`，`XtVinvX` 的 LDLT 对角出现极小负值并触发 `logL=-inf`。
- 2026-06-12：已完成两项修复：
  - 在预处理脚本中将 `TIME` 改为 `days_since_min_date * time_multiplier + within_day_index`，保证记录键唯一。
  - 在 `RepeatabilityModel::loadCovariates()` 中自动删除别名固定效应列，避免秩亏设计矩阵直接击穿 `exact` 路线。
- 2026-06-12：验证结果：
  - 大样本真实数据已不再“秒崩”，能稳定进入 `VCE` 计算。
  - 真实子集验证成功：`300` 头公猪、`28589` 条记录、`5000` SNP 的 `exact + --no-loco` 运行完成，得到有限方差分量和 `4252` 条 SNP 结果。
