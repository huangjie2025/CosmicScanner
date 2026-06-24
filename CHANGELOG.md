# 变更记录

## v1.1.0
- 重构：统一字符串处理，移除 `trim_copy2`，规范缺失值判断。
- 读取器：`readPhenGeneric` 取消 ID 列隐式回退，未找到将报错。
- 路径：统一使用 `std::filesystem::path` 进行拼接，提升跨平台性。
- 求解：PCG 在分母过小时增加一次回退到对角预条件并重启。
- 文档：完善 `docs/USAGE.md`（安装指南、输入输出、FAQ），新增顶层 `README.md`。
- 构建：`CMakeLists.txt` 增加安装与测试（ctest 占位）。

## v1.0.0
- 初始版本：支持 PBLUP/GBLUP/ssGBLUP，输出 `*.beta/*.rand/*.log`。