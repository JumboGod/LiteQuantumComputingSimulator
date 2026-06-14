#include "lqcs/algorithms/oracles.hpp"

#include <numbers>
#include <stdexcept>
#include <vector>

namespace lqcs::algorithms {

namespace {

// 对全 1 态 |11...1> 翻转相位（n=1 退化为 Z）
void flip_phase_all_ones(QuantumCircuit& qc, std::size_t n) {
    if (n == 1) {
        qc.z(0);
        return;
    }
    std::vector<qubit_t> controls(n - 1);
    for (std::size_t j = 0; j + 1 < n; ++j)
        controls[j] = static_cast<qubit_t>(j);
    qc.mcp(std::numbers::pi, controls, static_cast<qubit_t>(n - 1));
}

}  // namespace

QuantumCircuit xor_oracle(
    std::size_t n_in, std::size_t n_out,
    const std::function<std::uint64_t(std::uint64_t)>& f) {
    const std::size_t in_dim = std::size_t{1} << n_in;
    const std::size_t out_mask = (std::size_t{1} << n_out) - 1;
    const std::size_t dim = std::size_t{1} << (n_in + n_out);

    // 置换表：局部低 n_in 位为 x，高 n_out 位为 y
    std::vector<std::size_t> table(dim);
    for (std::size_t x = 0; x < in_dim; ++x) {
        const std::size_t fx = static_cast<std::size_t>(f(x));
        if (fx > out_mask) {
            throw std::invalid_argument(
                "xor_oracle: f(x) exceeds output register");
        }
        for (std::size_t y = 0; y <= out_mask; ++y) {
            table[(y << n_in) | x] = ((y ^ fx) << n_in) | x;
        }
    }

    QuantumCircuit qc(n_in + n_out);
    std::vector<qubit_t> qs(n_in + n_out);
    for (std::size_t i = 0; i < qs.size(); ++i) qs[i] = static_cast<qubit_t>(i);
    qc.permutation(table, qs);
    return qc;
}

QuantumCircuit phase_oracle(std::size_t n,
                            std::span<const std::uint64_t> marked) {
    if (marked.empty()) {
        throw std::invalid_argument("phase_oracle: no marked states");
    }
    QuantumCircuit qc(n);
    for (std::uint64_t m : marked) {
        if (m >= (std::uint64_t{1} << n)) {
            throw std::invalid_argument(
                "phase_oracle: marked state out of range");
        }
        // X 共轭把 |m> 映到 |11...1>，翻转相位后映回
        for (std::size_t b = 0; b < n; ++b) {
            if (!((m >> b) & 1u)) qc.x(static_cast<qubit_t>(b));
        }
        flip_phase_all_ones(qc, n);
        for (std::size_t b = 0; b < n; ++b) {
            if (!((m >> b) & 1u)) qc.x(static_cast<qubit_t>(b));
        }
    }
    return qc;
}

}  // namespace lqcs::algorithms
