#include "lqcs/algorithms/deutsch_jozsa.hpp"

#include "lqcs/algorithms/oracles.hpp"
#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

DeutschJozsaResult deutsch_jozsa(std::size_t n,
                                 const std::function<bool(std::uint64_t)>& f) {
    // 标准结构：辅助比特置 |->，oracle 相位回踢，再对输入寄存器做 H
    QuantumCircuit qc(n + 1);
    const auto ancilla = static_cast<qubit_t>(n);
    qc.x(ancilla);
    for (std::size_t q = 0; q <= n; ++q) qc.h(static_cast<qubit_t>(q));
    qc.compose(xor_oracle(n, 1, [&](std::uint64_t x) {
        return static_cast<std::uint64_t>(f(x) ? 1 : 0);
    }));
    for (std::size_t q = 0; q < n; ++q) qc.h(static_cast<qubit_t>(q));

    // 输入寄存器全 0 的概率：常数函数 = 1，平衡函数 = 0（确定性判据）
    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t in_mask = (std::size_t{1} << n) - 1;
    double p_zero = 0.0;
    for (std::size_t i = 0; i < sv.size(); ++i) {
        if ((i & in_mask) == 0) p_zero += std::norm(sv[i]);
    }
    return {p_zero > 0.5, p_zero, std::move(qc)};
}

}  // namespace lqcs::algorithms
