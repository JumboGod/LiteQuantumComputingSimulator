#pragma once

#include <cstddef>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

// 标准 QFT 电路（little-endian）：|x> → (1/√N) Σ_k e^{2πi·x·k/N} |k>
// do_swaps = false 时省略末尾的比特序翻转（输出为比特反序）
QuantumCircuit qft(std::size_t n, bool do_swaps = true);

// QFT 的逆电路
QuantumCircuit qft_inverse(std::size_t n, bool do_swaps = true);

}  // namespace lqcs::algorithms
