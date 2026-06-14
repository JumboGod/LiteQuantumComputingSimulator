#include "lqcs/algorithms/grover.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <unordered_set>

#include "lqcs/algorithms/oracles.hpp"
#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

GroverResult grover(std::size_t n, std::span<const std::uint64_t> marked,
                    std::optional<std::size_t> iterations) {
    const std::size_t dim = std::size_t{1} << n;
    if (marked.empty() || marked.size() >= dim) {
        throw std::invalid_argument(
            "grover: marked set must be non-empty and a strict subset");
    }
    std::unordered_set<std::uint64_t> marked_set(marked.begin(), marked.end());
    if (marked_set.size() != marked.size()) {
        throw std::invalid_argument("grover: duplicate marked states");
    }

    // 最优迭代次数：floor(π / (4·asin(√(M/N))))
    const double theta = std::asin(std::sqrt(
        static_cast<double>(marked.size()) / static_cast<double>(dim)));
    const std::size_t iters = iterations.value_or(std::max<std::size_t>(
        1, static_cast<std::size_t>(std::numbers::pi / (4.0 * theta))));

    const QuantumCircuit oracle = phase_oracle(n, marked);
    const std::vector<std::uint64_t> all_zero = {0};
    const QuantumCircuit flip_zero = phase_oracle(n, all_zero);

    QuantumCircuit qc(n, n);
    for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));
    for (std::size_t it = 0; it < iters; ++it) {
        qc.compose(oracle);
        // 扩散算子 = H^⊗n · (关于 |0...0> 的反射) · H^⊗n
        for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));
        qc.compose(flip_zero);
        for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));
    }

    // 末态统计（测量前）
    auto sv = StatevectorSimulator().run_statevector(qc);
    double p_success = 0.0;
    std::uint64_t best = 0;
    double best_p = -1.0;
    for (std::size_t i = 0; i < dim; ++i) {
        const double p = std::norm(sv[i]);
        if (marked_set.count(i)) p_success += p;
        if (p > best_p) {
            best_p = p;
            best = i;
        }
    }

    qc.measure_all();
    return {best, p_success, iters, std::move(qc)};
}

}  // namespace lqcs::algorithms
