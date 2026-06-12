#include "lqcs/algorithms/qft.hpp"

#include <numbers>

namespace lqcs::algorithms {

QuantumCircuit qft(std::size_t n, bool do_swaps) {
    QuantumCircuit qc(n);
    for (std::size_t i = n; i-- > 0;) {
        qc.h(static_cast<qubit_t>(i));
        for (std::size_t j = i; j-- > 0;) {
            const double angle =
                std::numbers::pi / static_cast<double>(1ULL << (i - j));
            qc.cp(angle, static_cast<qubit_t>(j), static_cast<qubit_t>(i));
        }
    }
    if (do_swaps) {
        for (std::size_t i = 0; i < n / 2; ++i) {
            qc.swap(static_cast<qubit_t>(i), static_cast<qubit_t>(n - 1 - i));
        }
    }
    return qc;
}

QuantumCircuit qft_inverse(std::size_t n, bool do_swaps) {
    return qft(n, do_swaps).inverse();
}

}  // namespace lqcs::algorithms
