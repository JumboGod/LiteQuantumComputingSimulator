#include "lqcs/algorithms/bernstein_vazirani.hpp"

#include "lqcs/algorithms/oracles.hpp"
#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

BernsteinVaziraniResult bernstein_vazirani(
    std::size_t n, const std::function<bool(std::uint64_t)>& f) {
    // 与 Deutsch-Jozsa 相同的相位回踢结构：末态输入寄存器恰为 |s>
    QuantumCircuit qc(n + 1);
    const auto ancilla = static_cast<qubit_t>(n);
    qc.x(ancilla);
    for (std::size_t q = 0; q <= n; ++q) qc.h(static_cast<qubit_t>(q));
    qc.compose(xor_oracle(n, 1, [&](std::uint64_t x) {
        return static_cast<std::uint64_t>(f(x) ? 1 : 0);
    }));
    for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));

    // 输入寄存器的边缘分布峰值即 s（理论上确定性）
    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t in_dim = std::size_t{1} << n;
    std::vector<double> marginal(in_dim, 0.0);
    for (std::size_t i = 0; i < sv.size(); ++i) {
        marginal[i & (in_dim - 1)] += std::norm(sv[i]);
    }
    std::uint64_t best = 0;
    for (std::size_t v = 1; v < in_dim; ++v) {
        if (marginal[v] > marginal[best]) best = v;
    }
    return {best, std::move(qc)};
}

}  // namespace lqcs::algorithms
