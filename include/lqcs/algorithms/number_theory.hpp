#pragma once

#include <cstdint>
#include <vector>

// 纯经典数论工具：Shor 的预处理与后处理
namespace lqcs::algorithms::number_theory {

std::uint64_t gcd(std::uint64_t a, std::uint64_t b);

// (a * b) mod m，128 位中间结果防溢出
std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b, std::uint64_t mod);

// (base ^ exp) mod m，快速模幂
std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exp, std::uint64_t mod);

bool is_prime(std::uint64_t n);  // 试除法（n 较小时足够）

// 若 n = b^k（k >= 2）返回 b，否则返回 0
std::uint64_t perfect_power_base(std::uint64_t n);

// num/den 的连分数收敛分数分母序列（<= max_denominator）。
// Shor 后处理：相位 m/2^t ≈ s/r，候选周期 r 即这些分母
std::vector<std::uint64_t> continued_fraction_denominators(
    std::uint64_t num, std::uint64_t den, std::uint64_t max_denominator);

}  // namespace lqcs::algorithms::number_theory
