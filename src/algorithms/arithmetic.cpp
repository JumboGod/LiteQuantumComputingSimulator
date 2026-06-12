#include "lqcs/algorithms/arithmetic.hpp"

#include <stdexcept>
#include <vector>

#include "lqcs/algorithms/number_theory.hpp"

namespace lqcs::algorithms {

QuantumCircuit modular_exponentiation(std::uint64_t a, std::uint64_t N,
                                      std::size_t n_counting, std::size_t n_work) {
    namespace nt = number_theory;
    if (N < 2 || (std::uint64_t{1} << n_work) < N) {
        throw std::invalid_argument(
            "modular_exponentiation: work register too small for N");
    }
    if (nt::gcd(a % N, N) != 1) {
        throw std::invalid_argument(
            "modular_exponentiation: gcd(a, N) must be 1 (else not a bijection)");
    }

    const std::size_t dim = std::size_t{1} << n_work;
    QuantumCircuit qc(n_counting + n_work);

    std::vector<qubit_t> qs(1 + n_work);  // {控制位, 工作寄存器...}
    for (std::size_t w = 0; w < n_work; ++w) {
        qs[1 + w] = static_cast<qubit_t>(n_counting + w);
    }

    for (std::size_t j = 0; j < n_counting; ++j) {
        const std::uint64_t factor = nt::pow_mod(a, std::uint64_t{1} << j, N);
        std::vector<std::size_t> table(dim);
        for (std::size_t y = 0; y < dim; ++y) {
            table[y] = (y < N) ? static_cast<std::size_t>(
                                     nt::mul_mod(factor, y, N))
                               : y;  // y >= N 的子空间恒等
        }
        qs[0] = static_cast<qubit_t>(j);
        qc.append({Gate{GateType::Permutation, {}, 1, {}, std::move(table)},
                   qs, {}});
    }
    return qc;
}

}  // namespace lqcs::algorithms
