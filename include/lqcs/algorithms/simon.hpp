#pragma once

#include <cstdint>
#include <functional>
#include <optional>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct SimonResult {
    std::uint64_t secret;  // 恢复出的周期串 s（f(x) = f(x ⊕ s)）
    std::size_t quantum_samples;  // 实际使用的量子采样数
    QuantumCircuit circuit;
};

// Simon 算法：f: {0,1}^n → {0,1}^n 满足 f(x) = f(y) ⇔ y = x ⊕ s（s ≠ 0）。
// 量子部分采样满足 y·s = 0 的随机 y；经典部分解 GF(2) 线性方程组得 s。
// 指数级加速：经典需要 Ω(2^(n/2)) 次查询，量子只需 O(n) 次。
SimonResult simon(std::size_t n,
                  const std::function<std::uint64_t(std::uint64_t)>& f,
                  std::optional<std::uint64_t> seed = std::nullopt);

}  // namespace lqcs::algorithms
