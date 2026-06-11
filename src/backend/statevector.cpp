#include "lqcs/backend/statevector.hpp"

#include <stdexcept>
#include <string>

namespace lqcs {

Statevector::Statevector(std::size_t num_qubits) : num_qubits_(num_qubits) {
    if (num_qubits == 0) {
        throw std::invalid_argument("Statevector: num_qubits must be >= 1");
    }
    if (num_qubits > kMaxQubits) {
        const double gib =
            static_cast<double>(sizeof(complex_t)) * (1ULL << (num_qubits - 30));
        throw std::invalid_argument(
            "Statevector: " + std::to_string(num_qubits) + " qubits would need " +
            std::to_string(gib) + " GiB (max " + std::to_string(kMaxQubits) +
            " qubits)");
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
    for (std::size_t i = 0; i < data_.size(); ++i) probs[i] = std::norm(data_[i]);
    return probs;
}

}  // namespace lqcs
