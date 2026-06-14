#pragma once

#include <cstdint>
#include <functional>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct BernsteinVaziraniResult {
    std::uint64_t secret;  // 恢复出的隐藏串 s
    QuantumCircuit circuit;
};

// Bernstein-Vazirani：f(x) = s·x mod 2（承诺为该形式），
// 单次 oracle 查询恢复隐藏串 s（经典需要 n 次查询）。
BernsteinVaziraniResult bernstein_vazirani(
    std::size_t n, const std::function<bool(std::uint64_t)>& f);

}  // namespace lqcs::algorithms
