# 贡献指南

感谢有兴趣为 LiteQuantumComputingSimulator 做贡献！

## 开发环境

- C++20 编译器（GCC 12+ / Clang 15+ / MSVC v143+）
- CMake ≥ 3.20
- 可选：Python 3.9+（构建 `pylqcs` 绑定）、clang-format 18（代码格式）

## 构建与测试

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

提交前请确保本地通过以下检查（与 CI 一致）：

| 检查 | 命令 |
|------|------|
| Release 构建 + 测试 | `cmake -B build && cmake --build build -j && ctest --test-dir build` |
| Sanitizers（ASan+UBSan） | `cmake -B b -DLQCS_ENABLE_SANITIZERS=ON -DLQCS_ENABLE_OPENMP=OFF && cmake --build b -j && ctest --test-dir b` |
| 标量回退（无 AVX2） | `cmake -B b -DLQCS_ENABLE_AVX2=OFF && cmake --build b -j && ctest --test-dir b` |
| Python 绑定 | `pip install . && pytest tests/python -q` |
| 代码格式 | `clang-format -i <改动的文件>`（配置见 `.clang-format`） |

## 代码约定

- **比特序**：little-endian（qubit 0 为最低有效位），与 Qiskit 一致。
- **正确性优先**：任何新内核/算法都应有与已知解析解或既有后端的交叉验证测试。
  状态向量、密度矩阵、轨迹法、stabilizer 各后端在 Clifford/理想电路上应给出一致分布。
- **错误尽早暴露**：非法输入在电路构建期抛异常，而非运行期。
- **性能内核**：保留标量回退路径；SIMD/并行仅作为快路径，由 `#if` 或运行时阈值守护。
- 命名空间 `lqcs`；公共头放 `include/lqcs/`，实现放 `src/`，二者镜像。

## 提交 PR

1. 从最新主分支切出特性分支。
2. 保持提交聚焦，提交信息说明「为什么」。
3. 新功能配套测试；确保上表所有检查通过。
4. CI 全绿后请求评审。

## 架构

整体设计见 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)。
