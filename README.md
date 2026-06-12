# LiteQuantumComputingSimulator

轻量级、高性能的量子计算模拟器，使用现代 C++（C++20）编写，API 设计参考 Qiskit。

详细架构见 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)。

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

## 设计约定

- **比特序与 Qiskit 一致（little-endian）**：qubit 0 为最低有效位，
  测量比特串最高位在左。
- 电路是纯数据（IR），执行由后端完成：`电路构建 → backend.run → Result`。
- 固定 `StatevectorSimulator::Options::seed` 可完全复现采样结果。
