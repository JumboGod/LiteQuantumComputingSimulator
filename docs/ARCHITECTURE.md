# LiteQuantumComputingSimulator 架构设计文档

> 一个轻量级、高性能的量子计算模拟器，使用现代 C++（C++20）编写，设计理念参考 Qiskit，
> 但聚焦于核心模拟能力，保持代码精简、易于理解和扩展。
>
> **v2 修订**：将「完整模拟 Shor 算法」与「Python 可视化」确立为一级需求，新增算法库层、
> Python 绑定层、置换门（Permutation Gate）专用内核与经典后处理模块。

---

## 1. 设计目标与原则

### 1.1 目标

| 目标 | 说明 |
|------|------|
| 正确性 | 严格遵循量子力学公设，所有门操作保持幺正性，测量遵循 Born 规则 |
| **Shor 算法（基线验收）** | 端到端分解 N=15、N=21：QFT⁻¹ + 受控模幂 + 测量 + 连分数后处理，全流程内置 |
| **Python 可视化** | pybind11 绑定，Python 侧用 matplotlib 绘制测量直方图、态振幅图、电路图 |
| 性能 | 状态向量模拟支持 ~28 个量子比特（单机 16GB 内存），核心循环 SIMD 友好 + OpenMP 并行 |
| 易用性 | C++ 与 Python 双 API，均提供类似 Qiskit 的流式电路构建体验 |
| 可扩展性 | 门库、后端、算法、可视化均通过抽象接口解耦 |
| 轻量级 | C++ 核心零第三方依赖；Python 绑定（pybind11）为可选构建项 |

**容量与 Shor 的关系**：教科书式 Shor 电路约需 `4n+2` 个量子比特（n 为待分解数 N 的位数）。
28 比特容量对应 n=6，即可分解 N < 64 的所有合数（15、21、33、35、39、51、55、57…），
满足教学与验证需求。

### 1.2 设计原则

- **分层架构**：上层 API 与底层数值内核严格分离，电路（做什么）与后端（怎么算）解耦。
- **电路即数据（Circuit as IR）**：`QuantumCircuit` 是一个纯数据结构（中间表示），
  不持有任何量子态；执行由后端完成。与 Qiskit 的 `circuit → transpile → backend.run` 流程一致。
- **C++ 做计算，Python 做交互与可视化**：所有数值计算留在 C++ 核心；Python 层只是薄绑定 +
  matplotlib 绘图，状态向量通过 numpy 零拷贝暴露。
- **经典-量子混合算法是一等公民**：Shor、VQE 等算法需要经典预/后处理（gcd、连分数、优化器），
  在 C++ 算法库中与电路构建器放在一起。
- **值语义优先**：小对象（门、指令）使用值语义；大对象（状态向量）使用移动语义。
- **编译期安全**：尽量用 `enum class`、`std::span`、`constexpr` 在编译期捕获错误。

---

## 2. 总体架构

```
┌─────────────────────────────────────────────────────────────┐
│                  Python 层 (pylqcs, 可选构建)                  │
│   pybind11 绑定 · matplotlib 可视化 (直方图/振幅图/电路图)       │
│   numpy 零拷贝状态向量 · Jupyter 友好                          │
├─────────────────────────────────────────────────────────────┤
│                     算法库层 (Algorithms)                      │
│   电路构建器: QFT/QFT⁻¹ · 受控模幂 · Grover 算子                │
│   完整算法: Shor (含连分数后处理) · Grover · VQE(v2)            │
│   经典工具: gcd · 模幂 · 连分数展开 · 周期验证                   │
├─────────────────────────────────────────────────────────────┤
│                        用户层 (User API)                      │
│   QuantumCircuit 流式构建 · 示例程序                           │
├─────────────────────────────────────────────────────────────┤
│                      电路中间表示层 (IR)                       │
│   QuantumCircuit · Instruction · Parameter                    │
│   门库: 标准门 · 多控制门 (MCX/MCP) · 置换门 · 自定义幺正         │
├─────────────────────────────────────────────────────────────┤
│                     编译/优化层 (Transpiler)                   │
│   Pass 管线：门分解 → 相邻门融合 → 冗余门消除                    │
├─────────────────────────────────────────────────────────────┤
│                      执行后端层 (Backends)                     │
│   Backend 抽象接口                                            │
│   ├─ StatevectorSimulator   (v1, 核心)                        │
│   ├─ DensityMatrixSimulator (v2, 噪声模拟)                     │
│   └─ StabilizerSimulator    (v3, Clifford 电路)                │
├─────────────────────────────────────────────────────────────┤
│                     数值内核层 (Kernels)                       │
│   单比特/双比特/多控制/对角/置换门特化内核                        │
│   SIMD + OpenMP 并行 · 测量采样 · 期望值计算                    │
├─────────────────────────────────────────────────────────────┤
│                      基础设施层 (Core)                         │
│   complex 类型 · 对齐内存分配器 · 随机数引擎 · 位运算工具          │
└─────────────────────────────────────────────────────────────┘
```

### 2.1 数据流（以 Shor 为例）

```
ShorResult shor = algorithms::shor(15);          // 用户一行调用
   │
   ├─ 经典预处理: 随机选 a, gcd(a,N)==1 ?
   ├─ 电路构建:   counting 寄存器 H^⊗2n
   │              + 受控模幂 (置换门序列)          ──► QuantumCircuit (IR)
   │              + QFT⁻¹ + 测量                        │
   ├─ 执行:       StatevectorSimulator::run  ◄──────────┘
   ├─ 经典后处理: 测量值 → 连分数 → 候选周期 r → 验证 a^r ≡ 1 (mod N)
   └─ 输出:       因子 {3, 5} + 中间数据 (供 Python 可视化)
```

---

## 3. 目录结构

```
LiteQuantumComputingSimulator/
├── CMakeLists.txt                  # 顶层构建脚本 (LQCS_BUILD_PYTHON 开关)
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
│   │   ├── gate.hpp                # 门的种类、矩阵定义（含置换门/多控制门）
│   │   ├── instruction.hpp         # 指令 = 门 + 作用比特 + 参数
│   │   ├── quantum_circuit.hpp     # 电路 IR 与流式 API
│   │   └── parameter.hpp           # 参数化门支持（变分算法用）
│   ├── algorithms/
│   │   ├── qft.hpp                 # QFT / QFT⁻¹ 电路构建器
│   │   ├── arithmetic.hpp          # 受控模幂等可逆算术电路构建器
│   │   ├── shor.hpp                # Shor 算法（电路 + 经典后处理）
│   │   ├── grover.hpp              # Grover 搜索
│   │   └── number_theory.hpp       # gcd、模幂、连分数（纯经典工具）
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
│       └── circuit_drawer.hpp      # ASCII 电路图绘制（C++ 侧兜底）
├── src/                            # 实现文件（与 include 镜像）
│   ├── circuit/...
│   ├── algorithms/...
│   ├── transpiler/...
│   ├── backend/
│   │   ├── statevector_simulator.cpp
│   │   └── kernels/                # 性能关键的门作用内核
│   │       ├── single_qubit_kernel.cpp
│   │       ├── two_qubit_kernel.cpp
│   │       ├── controlled_kernel.cpp   # 多控制门子空间内核
│   │       ├── permutation_kernel.cpp  # 置换门内核（模幂的关键）
│   │       └── measure_kernel.cpp
│   └── visualization/...
├── python/                         # Python 绑定与可视化包
│   ├── CMakeLists.txt
│   ├── pyproject.toml              # pip install -e . （scikit-build-core）
│   ├── src/bindings.cpp            # pybind11 绑定（单翻译单元起步）
│   └── pylqcs/
│       ├── __init__.py             # 重导出 C++ 类型 + Python 工具
│       ├── visualization.py        # plot_histogram / plot_state / 电路图
│       └── examples/
│           ├── bell.py
│           └── shor_demo.py        # Shor 分解 15/21 + 全套可视化
├── examples/                       # C++ 示例
│   ├── bell_state.cpp
│   ├── ghz_state.cpp
│   ├── qft.cpp
│   ├── grover.cpp
│   └── shor.cpp                    # Shor 分解 N=15 端到端
├── tests/                          # 单元测试（GoogleTest，FetchContent 引入）
│   ├── test_gates.cpp
│   ├── test_circuit.cpp
│   ├── test_statevector.cpp
│   ├── test_measurement.cpp
│   ├── test_transpiler.cpp
│   ├── test_qft.cpp                # QFT(n)|x⟩ 与解析 DFT 对比
│   ├── test_arithmetic.cpp         # 模幂置换门真值表验证
│   ├── test_shor.cpp               # 固定种子下分解 15/21 成功
│   └── python/test_bindings.py     # pytest：绑定层冒烟测试
└── benchmarks/
    └── bench_gates.cpp
```

命名空间统一为 `lqcs`；Python 包名 `pylqcs`。

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
该约定必须在文档、内核、算法库与可视化中全程一致，将用专门的单元测试锁定。

**随机数**：封装 `std::mt19937_64`，全局可设种子，保证测量采样与 Shor 的随机选 a 均可复现。

### 4.2 门库 (`circuit/gate.hpp`)

```cpp
enum class GateType {
    // 单比特无参门
    I, X, Y, Z, H, S, Sdg, T, Tdg, SX,
    // 单比特参数门
    RX, RY, RZ, P /*相位门*/, U /*通用 U(θ,φ,λ)*/,
    // 双比特门
    CX, CY, CZ, CH, CP, CRX, CRY, CRZ, SWAP, iSWAP, RXX, RYY, RZZ,
    // 三比特及多控制门
    CCX /*Toffoli*/, CSWAP /*Fredkin*/,
    MCX, MCP,        // 任意多控制 X / 相位门（QFT、Grover、算术电路必需）
    // 非幺正操作
    Measure, Reset, Barrier,
    // 自定义算符
    Unitary,         // 任意稠密幺正矩阵（构造时校验幺正性）
    Permutation,     // 经典可逆函数 f: i ↦ σ(i)（模幂的高效表示，见 4.4）
};

struct Gate {
    GateType                 type;
    std::vector<double>      params;       // 旋转角等
    std::vector<std::size_t> perm;         // 仅 Permutation 使用：作用子空间上的置换表
    std::vector<complex_t>   matrix() const;   // Permutation 也可展开成矩阵（仅测试用）
    Gate inverse() const;
    std::string name() const;
};
```

设计要点：

- **门矩阵按需生成**而非常驻存储；无参门矩阵以 `constexpr` 表格给出。
- **受控门不单独存矩阵**：内核对「控制位全为 1 的子空间」直接应用基门矩阵。
  `MCX`/`MCP` 支持任意数量控制位，受控模幂天然就是「计数比特控制 + 工作寄存器置换」。
- **`Permutation` 是 Shor 的关键抽象**（详见 4.4 与 4.7）：模幂 `|y⟩ ↦ |a·y mod N⟩`
  是计算基的置换，若用稠密 `Unitary` 表示需要 4ᵏ 复数存储和 O(4ᵏ) 计算；
  置换门只需一张大小 2ᵏ 的索引表，作用复杂度 O(2ⁿ)、零浮点乘法。
- 构造时校验：`Unitary` 检查 `U·U† = I`（容差 1e-10）；`Permutation` 检查是双射。

### 4.3 电路 IR (`circuit/`)

```cpp
struct Instruction {
    Gate                  gate;
    std::vector<qubit_t>  qubits;    // 受控门：控制位在前，目标位在后
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
    QuantumCircuit& cp(double theta, qubit_t control, qubit_t target);  // QFT 主力门
    QuantumCircuit& ccx(qubit_t c1, qubit_t c2, qubit_t target);
    QuantumCircuit& mcx(std::span<const qubit_t> controls, qubit_t target);
    QuantumCircuit& mcp(double theta, std::span<const qubit_t> controls, qubit_t target);
    QuantumCircuit& swap(qubit_t a, qubit_t b);
    QuantumCircuit& unitary(std::span<const complex_t> matrix,
                            std::span<const qubit_t> qubits);
    QuantumCircuit& permutation(std::span<const std::size_t> table,
                                std::span<const qubit_t> qubits);
    QuantumCircuit& measure(qubit_t q, clbit_t c);
    QuantumCircuit& measure_all();
    QuantumCircuit& barrier();
    // ...其余门同理

    // —— 组合与变换 ——
    QuantumCircuit& compose(const QuantumCircuit& other,
                            std::span<const qubit_t> qubit_map = {});
    // 整个子电路加上控制位（Shor 的 controlled-U^(2^j) 直接由此得到）
    QuantumCircuit  controlled(std::size_t num_controls = 1) const;
    QuantumCircuit  inverse() const;
    QuantumCircuit  power(int n) const;

    // —— 访问器 ——
    std::size_t num_qubits() const;
    std::size_t depth() const;
    std::span<const Instruction> instructions() const;

    std::string draw() const;              // ASCII 电路图
};
```

设计要点：

- 所有添加门的方法都做**边界检查**，错误在构建期暴露而不是运行期。
- `compose` + `controlled` + `inverse` 三件套是搭建 Shor 的基础：
  `模幂子电路 → controlled() → compose 到主电路`。

### 4.4 算法库 (`algorithms/`) —— 新增，Shor 的家

```cpp
namespace lqcs::algorithms {

// —— 电路构建器（返回可 compose 的子电路）——
QuantumCircuit qft(std::size_t n, bool do_swaps = true);
QuantumCircuit qft_inverse(std::size_t n, bool do_swaps = true);

// 受控模幂电路：|c⟩|y⟩ ↦ |c⟩|a^c · y mod N⟩
// 实现方式：对每个计数比特 j，构造置换 σ_j(y) = a^(2^j)·y mod N 的
// 单控制置换门。a^(2^j) mod N 由经典快速模幂预先算好——这是
// 状态向量模拟器的合理捷径：物理电路需要可逆加法器分解，但模拟器
// 可以直接应用置换，正确性完全等价，速度快几个数量级。
QuantumCircuit modular_exponentiation(std::uint64_t a, std::uint64_t N,
                                      std::size_t num_counting_qubits);

// —— 经典数论工具（纯经典，独立可测）——
namespace number_theory {
    std::uint64_t gcd(std::uint64_t a, std::uint64_t b);
    std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exp, std::uint64_t mod);
    // 连分数展开：phase ≈ s/r，返回所有收敛分数的分母（候选周期）
    std::vector<std::uint64_t> continued_fraction_denominators(
        double phase, std::uint64_t max_denominator);
}

// —— 完整 Shor 算法 ——
struct ShorOptions {
    std::size_t shots        = 64;     // 每轮采样数
    std::size_t max_attempts = 10;     // 最多随机尝试多少个 a
    std::optional<std::uint64_t> seed; // 可复现
};

struct ShorResult {
    bool          success;
    std::uint64_t N, factor1, factor2;
    std::uint64_t a, period;                       // 成功那一轮的底数与周期
    std::map<std::string, std::size_t> counts;     // 计数寄存器测量分布（供可视化）
    QuantumCircuit circuit;                        // 成功那一轮的电路（供绘图）
};

ShorResult shor(std::uint64_t N, const ShorOptions& opts = {});

QuantumCircuit grover(/* oracle 构建参数 */);
}
```

**Shor 主流程**（在 C++ 内闭环，Python 只负责调用与画图）：

1. 经典检查：N 偶数 / 素数幂直接返回；随机选 a，若 `gcd(a,N) > 1` 直接得因子。
2. 构建电路：2n 个计数比特 `H^⊗2n` → 受控模幂 → 计数寄存器 `QFT⁻¹` → 测量计数寄存器。
3. 运行 `shots` 次采样，对每个测量值 m：相位 `φ = m / 2^(2n)` → 连分数 → 候选 r。
4. 验证 `a^r ≡ 1 (mod N)` 且 r 为偶数且 `a^(r/2) ≢ -1`，则 `gcd(a^(r/2) ± 1, N)` 给出因子。
5. 失败则换 a 重试，直到成功或达到 `max_attempts`。

**资源估算**：N=15（n=4）→ 12 比特，状态向量 64 KB，毫秒级；
N=57（n=6）→ 18 比特，4 MB，亚秒级。完全在容量内。

### 4.5 编译/优化层 (`transpiler/`)

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
| `SingleQubitGateFusion` | 同一比特上连续的单比特门合并为一个 2×2 `Unitary` | 每融合一个门省一次全状态向量扫描（内存带宽是瓶颈） |
| `RedundantGateElimination` | 消除互逆相邻门对（H·H、X·X、CX·CX 等） | 直接减少门数 |

正确性保障：每个 pass 用小规模电路（≤ 6 比特）验证优化前后状态向量逐元素相等（容差 1e-12）。

### 4.6 执行后端层 (`backend/`)

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
    double expectation_value(std::string_view pauli) const;
    // Python 绑定用：暴露连续内存指针，零拷贝映射为 numpy 数组
    complex_t*       data();
    const complex_t* data() const;
};

class StatevectorSimulator final : public Backend {
public:
    struct Options {
        std::optional<std::uint64_t> seed;
        bool        fuse_gates  = true;
        int         num_threads = 0;       // 0 = OpenMP 默认
    };
    Result run(const QuantumCircuit&, std::size_t shots) override;
    Statevector run_statevector(const QuantumCircuit&);   // 返回末态（去掉测量）
};
```

**测量的两种路径**：

1. **电路末尾才测量**（Shor 即此类）：演化一次到末态，对概率分布做 `shots` 次
   别名采样，O(2ⁿ + shots·n)，无需重复演化。
2. **电路中间有测量/重置**：每个 shot 独立演化，按 Born 规则坍缩并归一化。自动检测。

### 4.7 数值内核层 (`backend/kernels/`) —— 性能核心

**单比特门内核**（最热路径）：

```cpp
void apply_single_qubit_gate(complex_t* state, std::size_t n_amps,
                             qubit_t k, const complex_t m[4]) {
    const std::size_t stride = 1uz << k;
    #pragma omp parallel for schedule(static)
    for (std::size_t g = 0; g < n_amps / 2; ++g) {
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
| 多控制门 (CX/CP/MCX/MCP…) | O(2ⁿ⁻ᶜ) | 仅遍历控制位全为 1 的子空间，c 个控制位省 2ᶜ 倍 |
| 通用双比特门 | O(2ⁿ) | 4×4 矩阵作用于振幅四元组 |
| **置换门（k 比特）** | **O(2ⁿ)** | 按置换的环分解原位搬移振幅，零浮点乘法；带控制位时同样走子空间遍历 |

**置换门内核**是 Shor 高效模拟的支柱：受控模幂若展开为基本门，N=15 就需要数千个
Toffoli；用置换门表示后整个模幂只是 2n 个受控置换，每个 O(2ⁿ) 一遍扫完。

**内存估算**（`complex<double>` = 16 字节）：

| 比特数 n | 状态向量大小 | 备注 |
|---------|-------------|------|
| 12 | 64 KB | Shor N=15 |
| 18 | 4 MB | Shor N<64 |
| 24 | 256 MB | |
| 28 | 4 GB | 默认上限附近 |

构造 `Statevector` 时检查 n 上限（默认 30），超限抛出带内存估算信息的异常。

### 4.8 结果层 (`result/`)

```cpp
class Result {
public:
    const std::map<std::string, std::size_t>& counts() const;  // {"01": 498, ...}
    std::map<std::string, double> probabilities() const;
    const Statevector* statevector() const;
    void print_counts(std::ostream& os = std::cout) const;
};
```

比特串格式与 Qiskit 一致：最高位在左。`counts`/`probabilities` 在 Python 侧
直接映射为 `dict`，喂给 matplotlib 即可画图。

### 4.9 Python 绑定与可视化 (`python/`) —— 新增

**绑定层（pybind11，薄封装）**，暴露的类型一一对应 C++ API：

```python
import pylqcs as lq

qc = lq.QuantumCircuit(2, 2)
qc.h(0).cx(0, 1).measure_all()

sim = lq.StatevectorSimulator()
result = sim.run(qc, shots=1024)
print(result.counts())                  # {'00': 498, '11': 526}

sv = sim.run_statevector(qc)
amps = sv.to_numpy()                    # 零拷贝 numpy complex128 数组
```

关键绑定设计：

| 内容 | 方案 |
|------|------|
| 状态向量 ↔ numpy | `py::array_t<complex_t>` 包装 `Statevector::data()`，零拷贝 + keep-alive 引用保证生命周期 |
| counts ↔ dict | `std::map<std::string,size_t>` 自动转换 |
| 流式 API | 绑定返回 `*this` 的方法，Python 侧同样支持链式调用 |
| GIL | `run()` 等长计算释放 GIL，不阻塞 Jupyter |
| 构建 | `LQCS_BUILD_PYTHON=ON` 时由 CMake 构建；`pip install -e python/` 走 scikit-build-core |

**可视化模块（纯 Python，`pylqcs.visualization`）**：

```python
from pylqcs.visualization import plot_histogram, plot_state, plot_circuit

plot_histogram(result.counts())          # 测量结果柱状图（类 Qiskit）
plot_state(sv)                           # 振幅模长 + 相位（双联图）
plot_circuit(qc)                         # 电路图（matplotlib 绘制）
```

| 函数 | 内容 |
|------|------|
| `plot_histogram(counts)` | 比特串-频次柱状图，自动按概率排序，支持多组对比 |
| `plot_state(sv)` | 上图：各基态振幅模方柱状图；下图：相位色环。Shor 演示中可直观看到周期峰 |
| `plot_circuit(qc)` | 读取电路 IR 指令列表，matplotlib 画线与门框（C++ 的 ASCII `draw()` 作为无 matplotlib 时的兜底） |
| `plot_shor(shor_result)` | Shor 专用：计数寄存器测量分布 + 标注理论峰位置 s·2²ⁿ/r |

**Shor 端到端 Python 演示**（`python/pylqcs/examples/shor_demo.py`）：

```python
import pylqcs as lq
from pylqcs.visualization import plot_shor

res = lq.algorithms.shor(15, seed=42)
print(f"15 = {res.factor1} × {res.factor2}  (a={res.a}, r={res.period})")
plot_shor(res)        # 测量分布图，峰落在 0, 256, 512, 768 (r=4)
```

依赖策略：`pylqcs` 核心绑定仅依赖 numpy；`visualization` 模块惰性导入 matplotlib，
未安装时给出清晰提示，不影响计算功能。

---

## 5. 构建、测试与质量保障

### 5.1 构建系统

- **CMake ≥ 3.20**，C++20。目标：`lqcs`（静态库）、`lqcs_tests`、`lqcs_examples`、
  `lqcs_bench`、`_pylqcs`（pybind11 模块，`LQCS_BUILD_PYTHON=ON` 时）。
- 编译选项：`-O3 -march=native -fopenmp`（Release）；`-fsanitize=address,undefined`（Debug/CI）。
- GoogleTest 与 pybind11 通过 `FetchContent` 引入，C++ 核心库本身零依赖。

### 5.2 测试策略

| 层级 | 方法 |
|------|------|
| 门矩阵 | 幺正性检查 `U·U† = I`；与解析真值对比 |
| 态演化 | 小电路（≤ 6 比特）与解析解逐振幅对比，容差 1e-12 |
| QFT | `QFT|x⟩` 与解析 DFT 公式对比；`QFT·QFT⁻¹ = I` |
| 算术电路 | 模幂置换门对所有输入做经典真值表验证 |
| **Shor 端到端** | 固定种子下分解 15、21、35 必须成功；统计分布峰位符合 s·2²ⁿ/r |
| 等价性 | 随机电路：优化前 vs 优化后、特化内核 vs 通用矩阵版结果一致 |
| 测量统计 | 固定种子 + 卡方检验 |
| 物理性质 | 范数恒为 1；Bell 关联、GHZ 奇偶性 |
| Python 绑定 | pytest 冒烟测试：链式 API、counts dict、numpy 零拷贝（修改 numpy 数组反映到 C++ 侧） |

### 5.3 性能基准

`benchmarks/` 跟踪：每门作用时间随 n 的伸缩曲线、QFT(24) 端到端耗时、
Shor(N=57) 端到端耗时、门融合开关对比。

---

## 6. 分阶段实施路线图

| 阶段 | 内容 | 验收标准 |
|------|------|----------|
| **M1 核心骨架** | core 基础设施、基本门（单比特+CX）、QuantumCircuit IR、Statevector、通用内核、末态测量采样、Result | Bell/GHZ 示例输出正确分布；单元测试全绿 |
| **M2 完整门集** | 全部门类型、多控制门子空间内核、**置换门内核**、对角门特化、中间测量、ASCII 绘图 | QFT、Grover 示例正确；随机电路等价性测试通过 |
| **M3 Shor 算法** ⭐ | algorithms 模块：QFT 构建器、受控模幂、数论工具、Shor 主流程 | **固定种子分解 15/21/35 成功**；测量分布峰位正确 |
| **M4 Python 绑定与可视化** ⭐ | pybind11 绑定、numpy 零拷贝、plot_histogram/plot_state/plot_circuit/plot_shor、shor_demo.py | Jupyter 中跑通 Shor 演示并出图；pytest 通过 |
| **M5 性能优化** | OpenMP 并行、门融合 pass、冗余消除 pass、基准测试 | 26 比特电路单门秒级以内；融合带来可测加速 |
| **M6 进阶能力** | Pauli 期望值、参数化电路（VQE）、OpenQASM 2.0 导入/导出、噪声后端 | vqe_h2 收敛；含噪 Bell 态保真度符合理论 |

> M3、M4 为本项目的基线验收（⭐）：Shor 端到端 + Python 可视化跑通即达成最基本要求。

---

## 7. 关键设计决策记录 (ADR)

| 决策 | 选择 | 理由 |
|------|------|------|
| 态的表示 | 状态向量优先 | 内存 O(2ⁿ) vs O(4ⁿ)，覆盖理想电路模拟场景 |
| 比特序 | little-endian（同 Qiskit） | 迁移零认知成本；位运算实现自然 |
| 门作用方式 | 振幅对原位变换（非矩阵乘法） | O(2ⁿ) vs O(4ⁿ)，零额外内存 |
| 受控门实现 | 子空间遍历（不展开大矩阵） | 计算量降 2ᶜ 倍，支持任意多控制位 |
| **模幂的表示** | **置换门（不分解为 Toffoli 网络）** | 模拟器无需物理可实现性；置换表示快几个数量级，且真值表可经典验证。物理级分解留作 M6 可选项 |
| **Shor 后处理位置** | C++ 算法库内（非 Python 脚本） | 算法自洽、C++/Python 两侧都能一行调用；Python 只做可视化 |
| **Python 集成方式** | pybind11 薄绑定 + 纯 Python 可视化模块 | 计算逻辑不重复实现；matplotlib 生态留在 Python 侧 |
| 状态向量传递 | numpy 零拷贝视图 | 28 比特 4GB 的态拷贝一次代价过高 |
| 并行模型 | OpenMP | 数据并行场景最简洁可靠 |
| 精度 | double（v1 不提供 float） | 长电路误差累积；简化 API |
| 错误处理 | 异常（构建期校验为主） | 电路构建是用户代码边界，fail-fast 最友好 |

---

## 8. 后续展望（不在当前范围）

- GPU 后端（CUDA/SYCL）：内核按数据并行扫描组织，天然可移植。
- 张量网络后端：用于低纠缠深电路。
- MPI 分布式状态向量：突破单机 30+ 比特限制。
- 模幂的物理级门分解（可逆加法器/Beauregard 电路）：用于研究真实资源开销。
