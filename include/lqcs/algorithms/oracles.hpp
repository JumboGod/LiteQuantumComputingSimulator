#pragma once

#include <cstdint>
#include <functional>
#include <span>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

// XOR oracle：U_f |x>|y> = |x>|y ⊕ f(x)>。
// 输入寄存器 [0, n_in)，输出寄存器 [n_in, n_in + n_out)。
// 任意经典函数 f 都是该形式下的计算基置换，用单个置换门实现。
QuantumCircuit xor_oracle(std::size_t n_in, std::size_t n_out,
                          const std::function<std::uint64_t(std::uint64_t)>& f);

// 相位 oracle：对每个 marked 基态翻转相位（|m> → -|m>）
QuantumCircuit phase_oracle(std::size_t n,
                            std::span<const std::uint64_t> marked);

}  // namespace lqcs::algorithms
