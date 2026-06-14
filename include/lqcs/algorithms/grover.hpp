#pragma once

#include <cstdint>
#include <optional>
#include <span>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct GroverResult {
    std::uint64_t best_state;    // 末态中概率最大的基态
    double success_probability;  // 测得任一 marked 态的概率
    std::size_t iterations;      // 实际迭代次数
    QuantumCircuit circuit;      // 含 measure_all，可直接 run
};

// Grover 搜索：在 2^n 个状态中放大 marked 态的振幅。
// iterations 缺省时取最优值 floor(π / (4·asin(√(M/2^n))))。
GroverResult grover(std::size_t n, std::span<const std::uint64_t> marked,
                    std::optional<std::size_t> iterations = std::nullopt);

}  // namespace lqcs::algorithms
