#include "lqcs/backend/statevector.hpp"

#include <bit>
#include <stdexcept>
#include <string>

#include "lqcs/core/bit_utils.hpp"

namespace lqcs {

Statevector::Statevector(std::size_t num_qubits) : num_qubits_(num_qubits) {
    if (num_qubits == 0) {
        throw std::invalid_argument("Statevector: num_qubits must be >= 1");
    }
    if (num_qubits > kMaxQubits) {
        const double gib = static_cast<double>(sizeof(complex_t)) *
                           (1ULL << (num_qubits - 30));
        throw std::invalid_argument(
            "Statevector: " + std::to_string(num_qubits) +
            " qubits would need " + std::to_string(gib) + " GiB (max " +
            std::to_string(kMaxQubits) + " qubits)");
    }
    data_.assign(std::size_t{1} << num_qubits, complex_t{0.0, 0.0});
    data_[0] = complex_t{1.0, 0.0};
}

double Statevector::norm() const {
    double sum = 0.0;
    for (const auto& amp : data_) sum += std::norm(amp);
    return std::sqrt(sum);
}

std::vector<double> Statevector::probabilities() const {
    std::vector<double> probs(data_.size());
    for (std::size_t i = 0; i < data_.size(); ++i)
        probs[i] = std::norm(data_[i]);
    return probs;
}

double Statevector::expectation_value(std::string_view pauli) const {
    if (pauli.size() != num_qubits_) {
        throw std::invalid_argument(
            "expectation_value: pauli string length must equal num_qubits");
    }
    // P|i> = i^{nY} · (-1)^{popcount(i & zy_mask)} · |i ⊕ x_mask>
    std::size_t x_mask = 0, zy_mask = 0;
    int n_y = 0;
    for (std::size_t j = 0; j < pauli.size(); ++j) {
        const std::size_t bit = std::size_t{1} << (num_qubits_ - 1 - j);
        switch (pauli[j]) {
            case 'I': break;
            case 'X': x_mask |= bit; break;
            case 'Z': zy_mask |= bit; break;
            case 'Y':
                x_mask |= bit;
                zy_mask |= bit;
                ++n_y;
                break;
            default:
                throw std::invalid_argument(
                    "expectation_value: pauli must contain only I/X/Y/Z");
        }
    }
    static constexpr complex_t kIPow[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    const complex_t y_phase = kIPow[n_y % 4];

    const std::size_t n = data_.size();
    double sum_re = 0.0, sum_im = 0.0;
#pragma omp parallel for schedule(static) \
    reduction(+ : sum_re, sum_im) if (n >= (std::size_t{1} << 14))
    for (std::size_t i = 0; i < n; ++i) {
        const double sign = (std::popcount(i & zy_mask) % 2 == 0) ? 1.0 : -1.0;
        const complex_t term = sign * data_[i] * std::conj(data_[i ^ x_mask]);
        sum_re += term.real();
        sum_im += term.imag();
    }
    return (y_phase * complex_t{sum_re, sum_im}).real();
}

}  // namespace lqcs
