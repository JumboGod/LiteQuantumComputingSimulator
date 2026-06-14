#include "lqcs/algorithms/shor.hpp"

#include <algorithm>
#include <bit>
#include <random>
#include <stdexcept>
#include <vector>

#include "lqcs/algorithms/arithmetic.hpp"
#include "lqcs/algorithms/number_theory.hpp"
#include "lqcs/algorithms/qft.hpp"
#include "lqcs/backend/statevector_simulator.hpp"
#include "lqcs/core/random.hpp"

namespace lqcs::algorithms {

namespace nt = number_theory;

namespace {

// 由候选周期 r 验证并提取因子；成功则填充 result 并返回 true
bool try_period(std::uint64_t N, std::uint64_t a, std::uint64_t r,
                ShorResult& result) {
    if (r == 0 || r >= N || nt::pow_mod(a, r, N) != 1) return false;
    if (r % 2 != 0) return false;
    const std::uint64_t x = nt::pow_mod(a, r / 2, N);
    if (x == N - 1) return false;  // a^(r/2) ≡ -1：无效
    const std::uint64_t f1 = nt::gcd(x - 1, N);
    const std::uint64_t f2 = nt::gcd(x + 1, N);
    for (std::uint64_t f : {f1, f2}) {
        if (f != 1 && f != N) {
            result.factor1 = std::min(f, N / f);
            result.factor2 = std::max(f, N / f);
            result.period = r;
            return true;
        }
    }
    return false;
}

}  // namespace

ShorResult shor(std::uint64_t N, const ShorOptions& options) {
    ShorResult result;
    result.N = N;
    if (N < 4) {
        throw std::invalid_argument("shor: N must be >= 4");
    }

    // —— 经典预处理 ——
    if (N % 2 == 0) {
        result.success = true;
        result.factor1 = 2;
        result.factor2 = N / 2;
        return result;
    }
    if (nt::is_prime(N)) {
        return result;  // 素数无非平凡因子：success = false
    }
    if (const std::uint64_t b = nt::perfect_power_base(N)) {
        result.success = true;
        result.factor1 = std::min(b, N / b);
        result.factor2 = std::max(b, N / b);
        return result;
    }

    const auto n_work =
        static_cast<std::size_t>(std::bit_width(N));  // 2^n_work > N
    const std::size_t n_count = 2 * n_work;
    const std::uint64_t count_dim = std::uint64_t{1} << n_count;

    Rng rng(options.seed);
    std::uniform_int_distribution<std::uint64_t> pick_a(2, N - 2);

    for (std::size_t attempt = 0; attempt < options.max_attempts; ++attempt) {
        result.attempts = attempt + 1;
        const std::uint64_t a = pick_a(rng.engine());
        result.a = a;

        // 幸运捷径：a 与 N 不互素时 gcd 直接给出因子
        const std::uint64_t g = nt::gcd(a, N);
        if (g != 1) {
            result.success = true;
            result.factor1 = std::min(g, N / g);
            result.factor2 = std::max(g, N / g);
            result.period = 0;
            return result;
        }

        // —— 量子求周期 ——
        QuantumCircuit qc(n_count + n_work, n_count);
        qc.x(static_cast<qubit_t>(n_count));  // 工作寄存器 |1>
        for (std::size_t j = 0; j < n_count; ++j) {
            qc.h(static_cast<qubit_t>(j));
        }
        qc.compose(modular_exponentiation(a, N, n_count, n_work));
        qc.compose(qft_inverse(n_count));  // 作用于计数寄存器（前 n_count 位）
        for (std::size_t j = 0; j < n_count; ++j) {
            qc.measure(static_cast<qubit_t>(j), static_cast<clbit_t>(j));
        }

        StatevectorSimulator::Options sv_opt;
        sv_opt.seed = rng.engine()();
        StatevectorSimulator sim(sv_opt);
        const Result run = sim.run(qc, options.shots);

        // —— 连分数后处理：按出现频次从高到低处理测量值 ——
        std::vector<std::pair<std::size_t, std::uint64_t>> outcomes;
        for (const auto& [key, count] : run.counts()) {
            outcomes.emplace_back(count, std::stoull(key, nullptr, 2));
        }
        std::sort(outcomes.rbegin(), outcomes.rend());

        for (const auto& [count, m] : outcomes) {
            if (m == 0) continue;  // 相位 0 无信息
            for (const std::uint64_t r_base :
                 nt::continued_fraction_denominators(m, count_dim, N)) {
                // 收敛分数可能给出 r 的约数：尝试小倍数
                for (std::uint64_t k = 1; k <= 4; ++k) {
                    if (try_period(N, a, k * r_base, result)) {
                        result.success = true;
                        result.counts = run.counts();
                        result.circuit = std::move(qc);
                        return result;
                    }
                }
            }
        }
    }
    return result;  // 所有尝试失败
}

}  // namespace lqcs::algorithms
