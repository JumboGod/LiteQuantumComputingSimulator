#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <string>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct ShorOptions {
    std::size_t shots = 64;             // 每轮量子采样数
    std::size_t max_attempts = 10;      // 最多随机尝试多少个底数 a
    std::optional<std::uint64_t> seed;  // 固定后全流程可复现
};

struct ShorResult {
    bool success = false;
    std::uint64_t N = 0;
    std::uint64_t factor1 = 0, factor2 = 0;
    std::uint64_t a = 0;       // 成功（或经典捷径）使用的底数
    std::uint64_t period = 0;  // 找到的周期 r（经典捷径时为 0）
    std::size_t attempts = 0;  // 实际尝试轮数
    std::map<std::string, std::size_t> counts;  // 计数寄存器测量分布
    std::optional<QuantumCircuit> circuit;      // 成功那一轮的电路
};

// Shor 质因数分解。
// 流程：经典预处理（偶数/素数/素数幂/幸运 gcd）→ 量子求周期
// （H^⊗2n + 受控模幂 + QFT⁻¹ + 测量）→ 连分数后处理 → 因子验证。
// N 为素数时返回 success = false；N < 4 抛出异常。
ShorResult shor(std::uint64_t N, const ShorOptions& options = {});

}  // namespace lqcs::algorithms
