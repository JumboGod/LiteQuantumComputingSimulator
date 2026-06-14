#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "../core/types.hpp"

namespace lqcs {

// 基门类型。受控门不单列：任何幺正基门都可以通过 Gate::n_controls
// 附加任意数量的控制位（CX = X + 1 控制位，CCX = X + 2，MCP = P + k，
// 受控模幂 = Permutation + 控制位）。
enum class GateType {
    // 单比特无参门
    I,
    X,
    Y,
    Z,
    H,
    S,
    Sdg,
    T,
    Tdg,
    SX,
    // 单比特参数门
    RX,
    RY,
    RZ,
    P,
    U,
    // 双比特门
    SWAP,
    iSWAP,
    RXX,
    RYY,
    RZZ,
    // 多比特 Pauli 旋转 exp(-iθ/2·P)，P 为 Pauli 串（VQE ansatz 主力门）
    PauliRotation,
    // 自定义算符
    Unitary,      // 任意幺正矩阵（构造时校验 U·U† = I）
    Permutation,  // 计算基置换 |l> → |perm[l]>（模幂的高效表示）
    // 非幺正操作
    Measure,
    Reset,
    Barrier,
};

struct Gate {
    GateType type;
    std::vector<double> params;     // 旋转角等
    std::size_t n_controls;         // 控制位数量（0 = 不受控）
    std::vector<complex_t> mat;     // 仅 Unitary：基门矩阵（行主序）
    std::vector<std::size_t> perm;  // 仅 Permutation：基门置换表
    std::string pauli;              // 仅 PauliRotation：Pauli 串
                        // （最右字符 ↔ 指令的首个目标位）

    Gate(GateType t, std::vector<double> p = {}, std::size_t nc = 0,
         std::vector<complex_t> m = {}, std::vector<std::size_t> pm = {},
         std::string pl = {})
        : type(t),
          params(std::move(p)),
          n_controls(nc),
          mat(std::move(m)),
          perm(std::move(pm)),
          pauli(std::move(pl)) {}

    // PauliRotation 的局部位掩码（局部位 j 对应 pauli[k-1-j]）
    struct PauliMasks {
        std::size_t x = 0;   // X/Y 位（翻转）
        std::size_t zy = 0;  // Z/Y 位（符号）
        unsigned n_y = 0;    // Y 的个数
    };
    PauliMasks pauli_masks() const;

    std::size_t base_qubits() const;  // 基门作用的比特数（不含控制位）
    std::size_t num_qubits() const { return n_controls + base_qubits(); }
    bool is_unitary() const;
    bool is_diagonal() const;  // 基门是否对角（内核走相位乘法快路径）

    // 基门矩阵 2^k × 2^k（不含控制位）
    std::vector<complex_t> base_matrix() const;

    // 含控制位的完整矩阵。局部比特序为 little-endian：局部位 j 对应
    // 指令 qubits[j]，控制位排在前（低位）。仅用于测试与参考实现，
    // 模拟内核走子空间路径，不经过该矩阵。
    std::vector<complex_t> matrix() const;

    Gate inverse() const;
    std::string name() const;  // 组合命名：x / cx / ccx / mcx ...
};

}  // namespace lqcs
