#include "lqcs/algorithms/simon.hpp"

#include <bit>
#include <stdexcept>
#include <vector>

#include "lqcs/algorithms/oracles.hpp"
#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

namespace {

unsigned highest_bit(std::uint64_t v) {
    return static_cast<unsigned>(std::bit_width(v)) - 1;
}

// GF(2) 高斯消元：把 y 化简进 echelon 基（主元 = 最高置位），满秩前返回基
void reduce_into_basis(std::vector<std::uint64_t>& basis, std::uint64_t y) {
    for (const std::uint64_t row : basis) {
        if (y == 0) return;
        if ((y >> highest_bit(row)) & 1u) y ^= row;
    }
    if (y != 0) basis.push_back(y);
}

// 由 rank = n-1 的方程组 y·s = 0 解出唯一非零 s
std::uint64_t solve_nullspace(std::vector<std::uint64_t> basis, std::size_t n) {
    std::uint64_t pivot_mask = 0;
    for (const std::uint64_t row : basis) {
        pivot_mask |= std::uint64_t{1} << highest_bit(row);
    }
    // 自由变量：取最低的非主元位，置 1
    std::uint64_t s = 0;
    for (std::size_t b = 0; b < n; ++b) {
        if (!((pivot_mask >> b) & 1u)) {
            s = std::uint64_t{1} << b;
            break;
        }
    }
    // 按主元从低到高回代：s_p = parity(row 的其余位 & s)
    std::sort(basis.begin(), basis.end());
    for (const std::uint64_t row : basis) {
        const unsigned p = highest_bit(row);
        const std::uint64_t rest = row ^ (std::uint64_t{1} << p);
        if (std::popcount(rest & s) % 2 == 1) s |= std::uint64_t{1} << p;
    }
    return s;
}

}  // namespace

SimonResult simon(std::size_t n,
                  const std::function<std::uint64_t(std::uint64_t)>& f,
                  std::optional<std::uint64_t> seed) {
    if (n < 2) {
        throw std::invalid_argument("simon: needs n >= 2");
    }
    // 输入寄存器 [0,n)，输出寄存器 [n,2n)
    QuantumCircuit qc(2 * n, n);
    for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));
    qc.compose(xor_oracle(n, n, f));
    for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));
    for (std::size_t q = 0; q < n; ++q) {
        qc.measure(static_cast<qubit_t>(q), static_cast<clbit_t>(q));
    }

    // 采样 y（必满足 y·s = 0），凑满 rank n-1 后解出 s
    StatevectorSimulator::Options sv_opt;
    sv_opt.seed = seed;
    StatevectorSimulator sim(sv_opt);
    std::vector<std::uint64_t> basis;
    std::size_t samples = 0;
    constexpr std::size_t kMaxRounds = 8;
    for (std::size_t round = 0; round < kMaxRounds; ++round) {
        const std::size_t shots = 8 * n;
        samples += shots;
        const Result run = sim.run(qc, shots);
        for (const auto& [key, count] : run.counts()) {
            const std::uint64_t y = std::stoull(key, nullptr, 2);
            if (y != 0) reduce_into_basis(basis, y);
        }
        if (basis.size() >= n - 1) break;
    }
    if (basis.size() < n - 1) {
        throw std::runtime_error(
            "simon: failed to collect n-1 independent equations");
    }
    return {solve_nullspace(std::move(basis), n), samples, std::move(qc)};
}

}  // namespace lqcs::algorithms
