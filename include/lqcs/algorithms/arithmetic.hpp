#pragma once

#include <cstddef>
#include <cstdint>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

// 受控模幂电路：|c>|y> → |c>|a^c · y mod N>（要求 gcd(a, N) = 1）。
// 计数寄存器 [0, n_counting)，工作寄存器 [n_counting, n_counting + n_work)。
//
// 实现：对每个计数比特 j 施加单控制置换 σ_j(y) = a^(2^j)·y mod N
// （y >= N 的子空间保持恒等）。a^(2^j) 由经典快速模幂预先算好——
// 这是状态向量模拟器的合理捷径：物理电路需要可逆加法器网络，
// 而模拟器直接应用置换，正确性等价、快几个数量级。
QuantumCircuit modular_exponentiation(std::uint64_t a, std::uint64_t N,
                                      std::size_t n_counting, std::size_t n_work);

}  // namespace lqcs::algorithms
