# LiteQuantumComputingSimulator

[![CI](https://github.com/JumboGod/LiteQuantumComputingSimulator/actions/workflows/ci.yml/badge.svg)](https://github.com/JumboGod/LiteQuantumComputingSimulator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

轻量级、高性能的量子计算模拟器，使用现代 C++（C++20）编写，API 设计参考 Qiskit。

详细架构见 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)，贡献指南见 [CONTRIBUTING.md](CONTRIBUTING.md)。

## 当前进度

- ✅ **M1 核心骨架**：基础设施、单比特门全集 + CX/CZ/CP/SWAP、电路 IR、
  状态向量模拟器、末态测量采样
- ✅ **M2 完整门集**：统一多控制门模型（任意基门 + n 控制位）、置换门内核、
  iSWAP/RXX/RYY/RZZ、自定义幺正、中间测量/Reset、compose/inverse/power、
  ASCII 电路图
- ✅ **M3 算法库**：Shor 质因数分解（QFT⁻¹ + 受控模幂 + 连分数后处理）、
  Deutsch-Jozsa、Bernstein-Vazirani、Simon、量子相位估计 QPE、Grover、
  QFT/oracle 构建器、经典数论工具
- ✅ **M4 Python 绑定与可视化**：pybind11 绑定（pylqcs 包）、numpy 零拷贝
  状态向量、matplotlib 可视化（直方图/态振幅图/电路图/Shor 分布图）
- ✅ **M5 性能优化**：OpenMP 并行内核（带规模阈值）、单比特门融合 pass
  （实测 5.3~5.7x 加速）、冗余门消除 pass、基准测试套件
- ✅ **M6 进阶能力**：Pauli 期望值、VQE（Rotosolve 优化器，H₂ 基态收敛到
  1e-12）、OpenQASM 2.0 导入/导出、密度矩阵噪声后端（去极化/振幅阻尼/
  相位阻尼等 Kraus 通道，含噪 Bell 保真度与理论解析值一致）
- ✅ **M7 轨迹法噪声**（对标 qsim/Aer）：StatevectorSimulator 直接支持
  噪声模型——逐 shot 按 ‖Kᵢ|ψ⟩‖² 采样 Kraus 分支，内存 O(2ⁿ)，
  含噪容量从密度矩阵的 13 比特提升到 26+ 比特；与密度矩阵精确解
  交叉验证统计一致
- ✅ **M8 参数化电路 + 梯度**（对标 Qulacs）：多比特 Pauli 旋转门
  exp(-iθ/2·P) 专用 O(2ⁿ) 内核（UCCSD ansatz 主力门）、ParametricCircuit
  参数槽原位绑定（bind 不重建结构）、parameter-shift 解析梯度、
  梯度下降 VQE
- ✅ **M9 多比特门融合**（对标 qsim/Aer）：贪心地把相邻门（含 CX/CCX）
  聚合成 ≤k 比特稠密块（默认 4），把多次全状态向量扫描合并为一次；
  实测相比无融合 ~30x、相比单比特融合再快 ~5x
- ✅ **M10 Stabilizer 模拟器**（对标 Stim/Aer）：Aaronson-Gottesman
  tableau 算法，Clifford 电路 O(n²) 内存/时间，实测 5000 比特 GHZ
  秒级完成（状态向量法需 2⁵⁰⁰⁰ 振幅）；与状态向量分布交叉验证一致，
  量子纠错码/图态模拟的关键能力
- ✅ **M11 显式 SIMD 内核**（对标 Qulacs）：AVX2+FMA 复数内核，单比特
  门按「块」结构遍历（两半各自连续，利于向量化），实测无融合路径
  ~2.3x、QFT(24) ~1.2x；状态向量模拟受内存带宽限制，增益主要体现在
  计算密集的对角/高位目标门。`LQCS_ENABLE_AVX2=OFF` 可回退纯标量

实测性能（4 核容器，`benchmarks/lqcs_bench`）：26 比特单门 135~275 ms，
QFT(24) 3.4 s，门融合对单比特门密集电路加速 5.3~5.7 倍。

## 快速开始

```cpp
#include <lqcs/lqcs.hpp>
using namespace lqcs;

int main() {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();      // Bell 态

    StatevectorSimulator sim;
    Result result = sim.run(qc, 1024);
    result.print_counts();               // {"00": ~512, "11": ~512}
}
```

## Python 接口（pylqcs）

```python
import pylqcs as lq
from pylqcs.visualization import plot_histogram, plot_shor

# Shor 分解 15
res = lq.algorithms.shor(15, seed=42)
print(f"{res.N} = {res.factor1} x {res.factor2}")
plot_shor(res)                      # 测量分布 + 理论峰位标注

# 电路构建与采样（与 C++ API 一致的链式调用）
qc = lq.QuantumCircuit(2, 2)
qc.h(0).cx(0, 1).measure_all()
result = lq.StatevectorSimulator(seed=7).run(qc, shots=1024)
plot_histogram(result.counts())

# 状态向量零拷贝映射为 numpy 数组
sv = lq.StatevectorSimulator().run_statevector(qc)
amps = sv.to_numpy()                # complex128, 无拷贝
```

构建 Python 绑定：

```bash
cmake -B build -DLQCS_BUILD_PYTHON=ON && cmake --build build -j
PYTHONPATH=python python3 python/pylqcs/examples/shor_demo.py
```

## 构建与测试

要求：CMake ≥ 3.20，支持 C++20 的编译器（GCC 12+ / Clang 15+）。

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# 运行单元测试（GoogleTest 由 FetchContent 自动获取）
ctest --test-dir build --output-on-failure

# 运行示例
./build/examples/bell_state
./build/examples/ghz_state
```

### CMake 选项

| 选项 | 默认 | 说明 |
|------|:----:|------|
| `LQCS_BUILD_TESTS` | ON | 构建 GoogleTest 单元测试 |
| `LQCS_BUILD_EXAMPLES` | ON | 构建示例程序 |
| `LQCS_BUILD_BENCHMARKS` | ON | 构建性能基准 |
| `LQCS_BUILD_PYTHON` | OFF | 构建 pybind11 Python 绑定 |
| `LQCS_ENABLE_OPENMP` | ON | OpenMP 并行内核 |
| `LQCS_ENABLE_AVX2` | ON | AVX2+FMA SIMD 内核（x86-64；自动探测，不支持则跳过） |
| `LQCS_ENABLE_SANITIZERS` | OFF | AddressSanitizer + UBSan（调试/CI） |

持续集成（[CI](.github/workflows/ci.yml)）在 Linux/macOS/Windows × Debug/Release
上构建并测试，另含 ASan+UBSan、无 AVX2 标量回退、Python 绑定与 clang-format 检查。

## Windows / Visual Studio 2022

项目已适配 MSVC（v143 工具集，C++20）。两种方式任选：

**方式一（推荐）：VS2022 直接打开 CMake 工程**

1. 安装 VS2022 时勾选「使用 C++ 的桌面开发」工作负载
   （含 MSVC v143、Windows SDK、CMake 工具）。
2. 启动 VS2022 → 「打开本地文件夹」→ 选择仓库根目录，
   VS 会自动识别 CMakeLists.txt 与 CMakePresets.json 并开始配置
   （首次配置需联网下载 GoogleTest）。
3. 顶部工具栏选择配置预设（vs2022 / Debug）→ 「生成」→「全部生成」。
4. 「选择启动项」下拉框选一个可执行文件（如 shor.exe、bell_state.exe）
   → F5 调试运行。
5. 「测试」→「测试资源管理器」可直接浏览并运行全部 GoogleTest 用例。

**方式二：生成 .sln 解决方案**

```bat
cmake --preset vs2022
:: 打开 build\vs2022\LiteQuantumComputingSimulator.sln
:: 在解决方案资源管理器中右键示例工程 → 设为启动项目 → F5
```

提示：
- 调试内核代码时建议关闭 OpenMP 以便单步：配置时加
  `-DLQCS_ENABLE_OPENMP=OFF`。
- Python 绑定需本机安装 Python 3.9+：`-DLQCS_BUILD_PYTHON=ON`，
  构建后 `set PYTHONPATH=python` 即可 `import pylqcs`。

## 设计约定

- **比特序与 Qiskit 一致（little-endian）**：qubit 0 为最低有效位，
  测量比特串最高位在左。
- 电路是纯数据（IR），执行由后端完成：`电路构建 → backend.run → Result`。
- 固定 `StatevectorSimulator::Options::seed` 可完全复现采样结果。
