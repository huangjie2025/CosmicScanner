# 贡献指南

- 提交前运行本地构建与测试：`cmake -S . -B build && cmake --build build && ctest --test-dir build`
- 代码风格：保持现有风格与库选择；避免在发布路径中加入测试辅助逻辑
- 提交内容：说明变更、影响范围与验证方法；如涉及性能优化请附基准结果