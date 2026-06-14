#pragma once

#include <cstdint>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct QPEResult {
    double phase;            // 估计的相位 φ ∈ [0, 1)
    std::uint64_t measured;  // 计数寄存器最可能的读数 m（φ ≈ m / 2^t）
    std::size_t n_counting;
    QuantumCircuit circuit;
};

// 量子相位估计：U|ψ> = e^{2πiφ}|ψ>，估计 φ（t 个计数比特给出 t 位精度）。
// u 为不带控制位的幺正门（U^(2^j) 由矩阵反复平方经典预计算）；
// eigenstate_prep 在工作寄存器上制备本征态 |ψ>。
// 寄存器布局：计数 [0, t)，工作 [t, t + k)。
QPEResult phase_estimation(const Gate& u, const QuantumCircuit& eigenstate_prep,
                           std::size_t n_counting);

}  // namespace lqcs::algorithms
