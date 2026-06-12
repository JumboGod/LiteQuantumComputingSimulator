#include "lqcs/algorithms/number_theory.hpp"

#include <cmath>

namespace lqcs::algorithms::number_theory {

std::uint64_t gcd(std::uint64_t a, std::uint64_t b) {
    while (b != 0) {
        const std::uint64_t r = a % b;
        a = b;
        b = r;
    }
    return a;
}

std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b, std::uint64_t mod) {
#ifdef __SIZEOF_INT128__
    return static_cast<std::uint64_t>(
        static_cast<unsigned __int128>(a) * b % mod);
#else
    // MSVC 无 128 位整型：双倍-累加（俄国农民乘法），O(64) 次迭代
    std::uint64_t result = 0;
    a %= mod;
    while (b > 0) {
        if (b & 1) result = (result + a) % mod;
        a = (a * 2) % mod;
        b >>= 1;
    }
    return result;
#endif
}

std::uint64_t pow_mod(std::uint64_t base, std::uint64_t exp, std::uint64_t mod) {
    std::uint64_t result = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mul_mod(result, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

bool is_prime(std::uint64_t n) {
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;
    for (std::uint64_t d = 3; d * d <= n; d += 2) {
        if (n % d == 0) return false;
    }
    return true;
}

std::uint64_t perfect_power_base(std::uint64_t n) {
    if (n < 4) return 0;
    for (unsigned k = 2; (1ULL << k) <= n; ++k) {
        // 整数 k 次根（浮点估计 + 邻域校正）
        auto b = static_cast<std::uint64_t>(
            std::round(std::pow(static_cast<double>(n), 1.0 / k)));
        for (std::uint64_t cand = (b > 1 ? b - 1 : 2); cand <= b + 1; ++cand) {
            std::uint64_t p = 1;
            bool overflow = false;
            for (unsigned i = 0; i < k; ++i) {
                if (p > n / cand) { overflow = true; break; }
                p *= cand;
            }
            if (!overflow && p == n) return cand;
        }
    }
    return 0;
}

std::vector<std::uint64_t> continued_fraction_denominators(
    std::uint64_t num, std::uint64_t den, std::uint64_t max_denominator) {
    // 收敛分数分母递推：k_i = a_i·k_{i-1} + k_{i-2}，k_{-1} = 0，k_{-2} = 1
    std::vector<std::uint64_t> out;
    std::uint64_t km1 = 0, km2 = 1;
    while (den != 0) {
        const std::uint64_t a = num / den;
        const std::uint64_t rem = num % den;
        num = den;
        den = rem;
        const std::uint64_t k = a * km1 + km2;
        if (k > max_denominator) break;
        out.push_back(k);
        km2 = km1;
        km1 = k;
    }
    return out;
}

}  // namespace lqcs::algorithms::number_theory
