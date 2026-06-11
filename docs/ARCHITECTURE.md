# LiteQuantumComputingSimulator 架构设计文档

> 一个轻量级、高性能的量子计算模拟器，使用现代 C++（C++20）编写，设计理念参考 Qiskit，
> 但聚焦于核心模拟能力，保持代码精简、易于理解和扩展。

---

## 1. 设计目标与原则

### 1.1 目标

| 目标 | 说明 |
|------|------|
| 正确性 | 严格遵循量子力学公设，所有门操作保持幺正性，测量遵循 Born 规则 |
| 性能 | 状态向量模拟支持 ~28 个量子比特（单机 16GB 内存），核心循环 SIMD 友好 + OpenMP 并行 |
| 易用性 | 提供类似 Qiskit 的流式（fluent）电路构建 API |
| 可扩展性 | 门库、后端、噪声模型均通过抽象接口解耦，便于后续添加密度矩阵模拟、张量网络等 |
| 轻量级 | 核心零第三方依赖（仅标准库），测试框架与构建工具除外 |

### 1.2 设计原则

- **分层架构**：上层 API 与底层数值内核严格分离，电路（做什么）与后端（怎么算）解耦。
- **电路即数据（Circuit as IR）**：`QuantumCircuit` 是一个纯数据结构（中间表示），
  不持有任何量子态；执行由后端完成。这与 Qiskit 的 `circuit → transpile → backend.run` 流程一致。
- **值语义优先**：小对象（门、指令）使用值语义；大对象（状态向量）使用移动语义，避免隐式深拷贝。
- **编译期安全**：尽量用 `enum class`、`std::span`、`constexpr` 等手段在编译期捕获错误。

---

## 2. 总体架构

```
┌─────────────────────────────────────────────────────────────┐
│                        用户层 (User API)                      │
│   QuantumCircuit 流式构建 · 预置算法 (QFT/Grover/...) · 示例    │
├─────────────────────────────────────────────────────────────┤
│                      电路中间表示层 (IR)                       │
│   QuantumCircuit · Instruction · QuantumRegister              │
│   ClassicalRegister · 参数化门 (Parameter)                     │
├─────────────────────────────────────────────────────────────┤
│                     编译/优化层 (Transpiler)                   │
│   Pass 管线：门分解 → 相邻门融合 → 冗余门消除 → 门重排            │
├─────────────────────────────────────────────────────────────┤
│                      执行后端层 (Backends)                     │
│   Backend 抽象接口                                            │
│   ├─ StatevectorSimulator   (v1, 核心)                        │
│   ├─ DensityMatrixSimulator (v2, 噪声模拟)                     │
│   └─ StabilizerSimulator    (v3, Clifford 电路)                │
├─────────────────────────────────────────────────────────────┤
│                     数值内核层 (Kernels)                       │
│   状态向量门作用内核 (单比特/双比特/受控/对角门特化)               │
│   SIMD + OpenMP 并行 · 测量采样 · 期望值计算                    │
├─────────────────────────────────────────────────────────────┤
│                      基础设施层 (Core)                         │
│   complex 类型 · 对齐内存分配器 · 随机数引擎 · 位运算工具          │
└─────────────────────────────────────────────────────────────┘
```

### 2.1 数据流

```
用户代码                 IR              编译优化            执行              结果
QuantumCircuit  ──►  Instruction  ──►  Transpiler  ──►  Backend::run  ──►  Result
  qc.h(0);            列表 (门序列)      (可选 pass)       状态向量演化        counts /
  qc.cx(0,1);                                            + 测量采样          statevector /
  qc.measure_all();                                                        expectation
```

---

## 3. 目录结构

```
LiteQuantumComputingSimulator/
├── CMakeLists.txt                  # 顶层构建脚本
├── README.md
├── docs/
│   └── ARCHITECTURE.md             # 本文档
├── include/lqcs/                   # 公共头文件（对外 API）
│   ├── lqcs.hpp                    # 总入口头文件
│   ├── core/
│   │   ├── types.hpp               # complex_t、qubit_t 等基础类型
│   │   ├── aligned_allocator.hpp   # 64 字节对齐分配器（SIMD/缓存行）
│   │   ├── bit_utils.hpp           # 位操作工具（索引展开等）
│   │   └── random.hpp              # 可设种子的随机数引擎封装
│   ├── circuit/
│   │   ├── gate.hpp                # 门的种类、矩阵定义
│   │   ├── instruction.hpp         # 指令 = 门 + 作用比特 + 参数
│   │   ├── quantum_circuit.hpp     # 电路 IR 与流式 API
│   │   └── parameter.hpp           # 参数化门支持（变分算法用）
│   ├── transpiler/
│   │   ├── pass.hpp                # Pass 抽象基类
│   │   ├── gate_fusion.hpp         # 相邻单比特门融合
│   │   └── transpiler.hpp          # Pass 管线调度
│   ├── backend/
│   │   ├── backend.hpp             # Backend 抽象接口
│   │   ├── statevector.hpp         # 状态向量容器
│   │   └── statevector_simulator.hpp
│   ├── result/
│   │   └── result.hpp              # 运行结果（counts、态、期望值）
│   └── visualization/
│       └── circuit_drawer.hpp      # ASCII 电路图绘制
├── src/                            # 实现文件（与 include 镜像）
│   ├── circuit/...
│   ├── transpiler/...
│   ├── backend/
│   │   ├── statevector_simulator.cpp
│   │   └── kernels/                # 性能关键的门作用内核
│   │       ├── single_qubit_kernel.cpp
│   │       ├── two_qubit_kernel.cpp
│   │       └── measure_kernel.cpp
│   └── visualization/...
├── examples/
│   ├── bell_state.cpp              # Bell 态制备与测量
│   ├── ghz_state.cpp               # GHZ 态
│   ├── qft.cpp                     # 量子傅里叶变换
│   ├── grover.cpp                  # Grover 搜索
│   └── vqe_h2.cpp                  # 变分算法示例（v2）
├── tests/                          # 单元测试（GoogleTest，FetchContent 引入）
│   ├── test_gates.cpp              # 门矩阵幺正性、已知真值表
│   ├── test_circuit.cpp            # IR 构建、寄存器边界检查
│   ├── test_statevector.cpp        # 态演化与解析解对比
│   ├── test_measurement.cpp        # 采样统计检验（卡方）
│   └── test_transpiler.cpp         # 优化前后等价性验证
└── benchmarks/
    └── bench_gates.cpp             # 门作用吞吐量基准
```

命名空间统一为 `lqcs`（**L**ite **Q**uantum **C**omputing **S**imulator）。

---

## 4. 各模块详细设计

### 4.1 基础设施层 (`core/`)

```cpp
namespace lqcs {
    using complex_t = std::complex<double>;   // 振幅类型
    using qubit_t   = std::uint32_t;          // 量子比特索引
    using clbit_t   = std::uint32_t;          // 经典比特索引

    // 64 字节对齐的向量，保证 AVX-512 / 缓存行对齐
    template <typename T>
    using aligned_vector = std::vector<T, AlignedAllocator<T, 64>>;
}
```

**关键约定 —— 比特序（与 Qiskit 一致，little-endian）**：
qubit 0 是最低有效位。状态 `|q_{n-1} … q_1 q_0⟩` 对应状态向量下标
`index = q_0·2⁰ + q_1·2¹ + … + q_{n-1}·2ⁿ⁻¹`。
该约定必须在文档、内核与可视化中全程一致，是最容易出错的地方，将用专门的单元测试锁定。

**随机数**：封装 `std::mt19937_64`，全局可设种子，保证测量采样可复现（测试必需）。

### 4.2 门库 (`circuit/gate.hpp`)

```cpp
enum class GateType {
    // 单比特无参门
    I, X, Y, Z, H, S, Sdg, T, Tdg, SX,
    // 单比特参数门
    RX, RY, RZ, P /*相位门*/, U /*通用 U(θ,φ,λ)*/,
    // 双比特门
    CX, CY, CZ, CH, CP, CRX, CRY, CRZ, SWAP, iSWAP, RXX, RYY, RZZ,
    // 三比特门
    CCX /*Toffoli*/, CSWAP /*Fredkin*/,
    // 非幺正操作
    Measure, Reset, Barrier,
    // 自定义幺正矩阵
    Unitary,
};

struct Gate {
    GateType            type;
    std::vector<double> params;            // 旋转角等
    // 返回 2^k × 2^k 门矩阵（行主序），k 为作用比特数
    std::vector<complex_t> matrix() const;
    Gate inverse() const;                  // 用于电路求逆 / 验证
    std::string name() const;
};
```

设计要点：

- **门矩阵按需生成**而非常驻存储——指令对象保持轻量，矩阵只在内核需要时构造。
- 无参门矩阵以 `constexpr` 表格给出；参数门由解析公式生成。
- `Unitary` 类型允许用户注入任意自定义幺正矩阵（构造时校验 `U·U† = I`，容差 1e-10）。
- 受控门**不单独存矩阵**：内核对"控制位为 1 的子空间"直接应用基门矩阵，节省内存与计算。

### 4.3 电路 IR (`circuit/`)

```cpp
struct Instruction {
    Gate                  gate;
    std::vector<qubit_t>  qubits;    // 作用的量子比特（受控门：控制位在前）
    std::vector<clbit_t>  clbits;    // 仅 Measure 使用
};

class QuantumCircuit {
public:
    explicit QuantumCircuit(std::size_t num_qubits, std::size_t num_clbits = 0);

    // —— 流式 API（返回 *this，支持链式调用）——
    QuantumCircuit& h(qubit_t q);
    QuantumCircuit& x(qubit_t q);
    QuantumCircuit& rx(double theta, qubit_t q);
    QuantumCircuit& cx(qubit_t control, qubit_t target);
    QuantumCircuit& ccx(qubit_t c1, qubit_t c2, qubit_t target);
    QuantumCircuit& unitary(std::span<const complex_t> matrix,
                            std::span<const qubit_t> qubits);
    QuantumCircuit& measure(qubit_t q, clbit_t c);
    QuantumCircuit& measure_all();
    QuantumCircuit& barrier();
    // ...其余门同理

    // —— 组合与变换 ——
    QuantumCircuit& compose(const QuantumCircuit& other,
                            std::span<const qubit_t> qubit_map = {});
    QuantumCircuit  inverse() const;       // 全电路求逆（U† 序列反转）
    QuantumCircuit  power(int n) const;

    // —— 访问器 ——
    std::size_t num_qubits() const;
    std::size_t depth() const;             // 电路深度（关键路径长度）
    std::span<const Instruction> instructions() const;

    std::string draw() const;              // ASCII 电路图
private:
    std::size_t              num_qubits_, num_clbits_;
    std::vector<Instruction> instructions_;
};
```

设计要点：

- 所有添加门的方法都做**边界检查**（比特索引越界、控制位与目标位重复立即抛 `std::out_of_range` / `std::invalid_argument`），错误在构建期暴露而不是运行期。
- `compose` 支持把子电路（如 QFT 模块）映射到大电路的任意比特子集上，这是搭建算法库的基础。
- 使用示例：

```cpp
#include <lqcs/lqcs.hpp>
using namespace lqcs;

int main() {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();        // Bell 态

    StatevectorSimulator sim;
    Result result = sim.run(qc, /*shots=*/1024);
    result.print_counts();                 // {"00": ~512, "11": ~512}
}
```

### 4.4 编译/优化层 (`transpiler/`)

```cpp
class Pass {
public:
    virtual ~Pass() = default;
    virtual QuantumCircuit run(const QuantumCircuit& circuit) const = 0;
};

class Transpiler {
public:
    Transpiler& add_pass(std::unique_ptr<Pass> pass);
    QuantumCircuit run(const QuantumCircuit& circuit) const;
};
```

v1 实现两个最有价值的 pass：

| Pass | 作用 | 收益 |
|------|------|------|
| `SingleQubitGateFusion` | 同一比特上连续的单比特门合并为一个 2×2 `Unitary` | 内存带宽是状态向量模拟瓶颈，每融合一个门就省一次全状态向量扫描 |
| `RedundantGateElimination` | 消除互逆相邻门对（H·H、X·X、CX·CX 等） | 直接减少门数 |

正确性保障：每个 pass 的单元测试用小规模电路（≤ 6 比特）验证优化前后状态向量逐元素相等（容差 1e-12）。

### 4.5 执行后端层 (`backend/`)

```cpp
class Backend {
public:
    virtual ~Backend() = default;
    virtual Result run(const QuantumCircuit& circuit, std::size_t shots = 1024) = 0;
    virtual std::string name() const = 0;
};

class Statevector {
public:
    explicit Statevector(std::size_t num_qubits);    // 初始化为 |0…0⟩
    complex_t& operator[](std::size_t i);
    std::size_t num_qubits() const;
    double norm() const;
    std::vector<double> probabilities() const;
    // 期望值 ⟨ψ|P|ψ⟩，P 为 Pauli 串（如 "ZZI"）
    double expectation_value(std::string_view pauli) const;
private:
    aligned_vector<complex_t> data_;       // 长度 2^n
};

class StatevectorSimulator final : public Backend {
public:
    struct Options {
        std::optional<std::uint64_t> seed;
        bool        fuse_gates  = true;    // 自动启用门融合
        int         num_threads = 0;       // 0 = OpenMP 默认
    };
    Result run(const QuantumCircuit&, std::size_t shots) override;
    Statevector run_statevector(const QuantumCircuit&);   // 返回末态（去掉测量）
};
```

**测量的两种路径**（重要的性能设计）：

1. **电路末尾才测量**（最常见）：演化一次到末态，对概率分布做 `shots` 次**别名采样
   （alias sampling）或二分查找采样**，O(2ⁿ + shots·n)。无需重复演化电路。
2. **电路中间有测量/重置**（测量后还有门）：每个 shot 独立演化，测量时按 Born 规则
   坍缩并归一化。自动检测电路属于哪种情况。

### 4.6 数值内核层 (`backend/kernels/`) —— 性能核心

**单比特门内核**（最热的代码路径）：

对作用在 qubit `k` 上的门 `[[a,b],[c,d]]`，状态向量中下标仅在第 `k` 位不同的振幅对
`(α₀, α₁)` 做 2×2 变换。通过位运算直接枚举所有振幅对，无需构造 2ⁿ×2ⁿ 矩阵：

```cpp
void apply_single_qubit_gate(complex_t* state, std::size_t n_amps,
                             qubit_t k, const complex_t m[4]) {
    const std::size_t stride = 1uz << k;
    #pragma omp parallel for schedule(static)
    for (std::size_t g = 0; g < n_amps / 2; ++g) {
        // 把 g 展开成「第 k 位为 0」的下标 i0
        std::size_t i0 = ((g >> k) << (k + 1)) | (g & (stride - 1));
        std::size_t i1 = i0 | stride;
        complex_t a0 = state[i0], a1 = state[i1];
        state[i0] = m[0] * a0 + m[1] * a1;
        state[i1] = m[2] * a0 + m[3] * a1;
    }
}
```

复杂度与特化策略：

| 内核 | 复杂度 | 特化优化 |
|------|--------|----------|
| 通用单比特门 | O(2ⁿ) | OpenMP 并行；k 较小时内层连续访存，编译器自动向量化 |
| 对角门 (Z/S/T/RZ/P) | O(2ⁿ) | 只乘相位、不交换振幅，访存减半 |
| X 门 | O(2ⁿ) | 纯 swap，无复数乘法 |
| 受控门 (CX/CZ/CP…) | O(2ⁿ⁻ᶜ) | 仅遍历控制位全为 1 的子空间，c 个控制位省 2ᶜ 倍工作量 |
| 通用双比特门 | O(2ⁿ) | 4×4 矩阵作用于振幅四元组 |

**内存估算**（`complex<double>` = 16 字节）：

| 比特数 n | 状态向量大小 |
|---------|-------------|
| 20 | 16 MB |
| 24 | 256 MB |
| 28 | 4 GB |
| 30 | 16 GB |

构造 `Statevector` 时检查 n 上限（默认 30），超限抛出带清晰内存估算信息的异常。

### 4.7 结果层 (`result/`)

```cpp
class Result {
public:
    const std::map<std::string, std::size_t>& counts() const;  // {"01": 498, ...}
    std::map<std::string, double> probabilities() const;
    const Statevector* statevector() const;    // 仅在保存末态时非空
    void print_counts(std::ostream& os = std::cout) const;
};
```

比特串格式与 Qiskit 一致：最高位在左（`"q_{n-1} … q_1 q_0"`）。

### 4.8 可视化 (`visualization/`)

ASCII 电路图（`qc.draw()`），v1 实现：

```
q_0: ─[H]──●──[M]─
           │   ║
q_1: ─────[X]─[M]─
               ║
c: 2/══════════╩══
```

---

## 5. 构建、测试与质量保障

### 5.1 构建系统

- **CMake ≥ 3.20**，C++20，目标分为 `lqcs`（静态库）、`lqcs_tests`、`lqcs_examples`、`lqcs_bench`。
- 编译选项：`-O3 -march=native -fopenmp`（Release）；`-fsanitize=address,undefined`（Debug/CI）。
- GoogleTest 通过 `FetchContent` 引入，核心库本身零依赖。

### 5.2 测试策略（正确性是模拟器的生命线）

| 层级 | 方法 |
|------|------|
| 门矩阵 | 幺正性检查 `U·U† = I`；与解析真值对比（如 H|0⟩ = |+⟩） |
| 态演化 | 小电路（≤ 6 比特）与手算/解析解逐振幅对比，容差 1e-12 |
| 等价性 | 随机电路：优化前 vs 优化后、内核特化版 vs 通用矩阵版，结果必须一致 |
| 测量统计 | 固定种子 + 卡方检验，验证采样分布符合理论概率 |
| 物理性质 | 任意门序列后范数恒为 1；Bell 态关联、GHZ 奇偶性等已知物理结论 |
| 边界 | 0 门电路、单比特电路、索引越界、非幺正自定义矩阵拒绝 |

### 5.3 性能基准

`benchmarks/` 跟踪关键指标：每门作用时间随 n 的伸缩曲线、QFT(24) 端到端耗时、门融合开关对比。

---

## 6. 分阶段实施路线图

| 阶段 | 内容 | 验收标准 |
|------|------|----------|
| **M1 核心骨架** | core 基础设施、门库（单比特+CX）、QuantumCircuit IR、Statevector、通用内核、末态测量采样、Result | Bell/GHZ 示例输出正确分布；单元测试全绿 |
| **M2 完整门集** | 全部门类型、受控门子空间内核、对角门特化、中间测量/Reset、ASCII 绘图 | QFT、Grover 示例正确；随机电路等价性测试通过 |
| **M3 性能优化** | OpenMP 并行、门融合 pass、冗余消除 pass、基准测试 | 26 比特电路可在数秒级单门作用；融合带来可测量加速 |
| **M4 进阶能力** | Pauli 期望值、参数化电路（VQE 支撑）、OpenQASM 2.0 导入/导出 | vqe_h2 示例收敛到基态能量 |
| **M5 扩展后端** | DensityMatrixSimulator + 噪声通道（去极化/振幅阻尼）、StabilizerSimulator | 含噪 Bell 态保真度符合理论值 |

---

## 7. 关键设计决策记录 (ADR)

| 决策 | 选择 | 理由 |
|------|------|------|
| 态的表示 | 状态向量优先（而非密度矩阵） | 内存 O(2ⁿ) vs O(4ⁿ)，覆盖最常见的理想电路模拟场景 |
| 比特序 | little-endian（同 Qiskit） | 用户从 Qiskit 迁移零认知成本；位运算实现自然 |
| 门作用方式 | 振幅对原位变换（而非矩阵乘法） | O(2ⁿ) vs O(4ⁿ)，且零额外内存 |
| 受控门实现 | 子空间遍历（不展开成大矩阵） | 计算量降 2ᶜ 倍，且支持任意多控制位 |
| 并行模型 | OpenMP（而非手写线程池） | 数据并行场景下最简洁可靠，单行 pragma 即可 |
| 精度 | double（不提供 float 切换，v1） | 长电路误差累积；简化 API。后续可模板化 |
| 错误处理 | 异常（构建期校验为主） | 电路构建是用户代码边界，fail-fast 最友好 |

---

## 8. 后续展望（不在 v1 范围）

- GPU 后端（CUDA/SYCL）：内核已按"对状态向量的数据并行扫描"组织，天然可移植。
- 张量网络后端：用于低纠缠深电路。
- Python 绑定（pybind11）：复刻 Qiskit 的使用体验。
- MPI 分布式状态向量：突破单机 30+ 比特限制。
