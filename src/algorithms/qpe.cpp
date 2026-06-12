#include "lqcs/algorithms/qpe.hpp"

#include <stdexcept>
#include <vector>

#include "lqcs/algorithms/qft.hpp"
#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

namespace {

std::vector<complex_t> matmul(const std::vector<complex_t>& a,
                              const std::vector<complex_t>& b, std::size_t dim) {
    std::vector<complex_t> c(dim * dim, complex_t{0, 0});
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t l = 0; l < dim; ++l)
            for (std::size_t j = 0; j < dim; ++j)
                c[i * dim + j] += a[i * dim + l] * b[l * dim + j];
    return c;
}

}  // namespace

QPEResult phase_estimation(const Gate& u, const QuantumCircuit& eigenstate_prep,
                           std::size_t n_counting) {
    if (u.n_controls != 0 || !u.is_unitary()) {
        throw std::invalid_argument(
            "phase_estimation: u must be an uncontrolled unitary gate");
    }
    const std::size_t k = u.base_qubits();
    if (eigenstate_prep.num_qubits() != k) {
        throw std::invalid_argument(
            "phase_estimation: eigenstate_prep must act on the same qubits as u");
    }
    const std::size_t t = n_counting;
    const std::size_t dim = std::size_t{1} << k;

    QuantumCircuit qc(t + k, t);

    // 工作寄存器制备本征态
    std::vector<qubit_t> work_map(k);
    for (std::size_t w = 0; w < k; ++w) work_map[w] = static_cast<qubit_t>(t + w);
    qc.compose(eigenstate_prep, work_map);

    // 计数寄存器叠加 + 受控 U^(2^j)（矩阵经典反复平方）
    for (std::size_t j = 0; j < t; ++j) qc.h(static_cast<qubit_t>(j));
    std::vector<complex_t> u_pow = u.base_matrix();
    std::vector<qubit_t> qs(1 + k);
    for (std::size_t w = 0; w < k; ++w) qs[1 + w] = work_map[w];
    for (std::size_t j = 0; j < t; ++j) {
        qs[0] = static_cast<qubit_t>(j);
        qc.append({Gate{GateType::Unitary, {}, 1, u_pow}, qs, {}});
        u_pow = matmul(u_pow, u_pow, dim);
    }

    qc.compose(qft_inverse(t));
    for (std::size_t j = 0; j < t; ++j) {
        qc.measure(static_cast<qubit_t>(j), static_cast<clbit_t>(j));
    }

    // 读出：计数寄存器边缘分布的峰值 m，φ ≈ m / 2^t
    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t count_dim = std::size_t{1} << t;
    std::vector<double> marginal(count_dim, 0.0);
    for (std::size_t i = 0; i < sv.size(); ++i) {
        marginal[i & (count_dim - 1)] += std::norm(sv[i]);
    }
    std::uint64_t best = 0;
    for (std::size_t m = 1; m < count_dim; ++m) {
        if (marginal[m] > marginal[best]) best = m;
    }
    const double phase =
        static_cast<double>(best) / static_cast<double>(count_dim);
    return {phase, best, t, std::move(qc)};
}

}  // namespace lqcs::algorithms
