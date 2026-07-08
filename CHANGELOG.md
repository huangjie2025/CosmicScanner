# 变更记录

## v1.1.0
- 发布阻断项清理：移除 repeatability 路径中的 `curl`/`std::system` 调试上报，改为本地诊断输出。
- 发布阻断项清理：程序启动版本改为使用 CMake `PROJECT_VERSION` 注入的 `COSMICSCANNER_VERSION`。
- 测试门禁：将 `test_sparse_reml`、`test_multi_vc`、`test_spa_accuracy`、`test_regenie`、`test_input_adaptation_units` 纳入 CTest。
- Biobank/PGEN：为 PGEN biobank matrix-free scan 增加 keep-index 子集适配，避免过滤样本后的全样本矩阵乘法错配。
- 发布收口：补充 CosmicGWAS 发布状态、稳定/条件/实验功能边界和最小验收命令。
- 文档：更新 `README.md` 与 `docs/GWAS_RELEASE_CLOSEOUT_CN.md`，记录阻断项清理结果、CTest 覆盖范围和发布边界。
- 重构：统一字符串处理，移除 `trim_copy2`，规范缺失值判断。
- 读取器：`readPhenGeneric` 取消 ID 列隐式回退，未找到将报错。
- 路径：统一使用 `std::filesystem::path` 进行拼接，提升跨平台性。
- 求解：PCG 在分母过小时增加一次回退到对角预条件并重启。
- 文档：完善 `docs/USAGE.md`（安装指南、输入输出、FAQ），新增顶层 `README.md`。
- 构建：`CMakeLists.txt` 增加安装与测试（ctest 占位）。

## v1.0.0
- 初始版本：支持 PBLUP/GBLUP/ssGBLUP，输出 `*.beta/*.rand/*.log`。
