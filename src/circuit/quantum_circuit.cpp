#include "lqcs/circuit/quantum_circuit.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace lqcs {

QuantumCircuit::QuantumCircuit(std::size_t num_qubits, std::size_t num_clbits)
    : num_qubits_(num_qubits), num_clbits_(num_clbits) {
    if (num_qubits == 0) {
        throw std::invalid_argument("QuantumCircuit: num_qubits must be >= 1");
    }
}

void QuantumCircuit::check_qubit(qubit_t q) const {
    if (q >= num_qubits_) {
        throw std::out_of_range("qubit index " + std::to_string(q) +
                                " out of range (num_qubits = " +
                                std::to_string(num_qubits_) + ")");
    }
}

void QuantumCircuit::check_clbit(clbit_t c) const {
    if (c >= num_clbits_) {
        throw std::out_of_range("clbit index " + std::to_string(c) +
                                " out of range (num_clbits = " +
                                std::to_string(num_clbits_) + ")");
    }
}

QuantumCircuit& QuantumCircuit::add_1q(GateType type, qubit_t q,
                                       std::vector<double> params) {
    check_qubit(q);
    instructions_.push_back({Gate{type, std::move(params)}, {q}, {}});
    return *this;
}

QuantumCircuit& QuantumCircuit::add_2q(GateType type, qubit_t a, qubit_t b,
                                       std::vector<double> params) {
    check_qubit(a);
    check_qubit(b);
    if (a == b) {
        throw std::invalid_argument("two-qubit gate requires distinct qubits");
    }
    instructions_.push_back({Gate{type, std::move(params)}, {a, b}, {}});
    return *this;
}

QuantumCircuit& QuantumCircuit::i(qubit_t q)   { return add_1q(GateType::I, q); }
QuantumCircuit& QuantumCircuit::x(qubit_t q)   { return add_1q(GateType::X, q); }
QuantumCircuit& QuantumCircuit::y(qubit_t q)   { return add_1q(GateType::Y, q); }
QuantumCircuit& QuantumCircuit::z(qubit_t q)   { return add_1q(GateType::Z, q); }
QuantumCircuit& QuantumCircuit::h(qubit_t q)   { return add_1q(GateType::H, q); }
QuantumCircuit& QuantumCircuit::s(qubit_t q)   { return add_1q(GateType::S, q); }
QuantumCircuit& QuantumCircuit::sdg(qubit_t q) { return add_1q(GateType::Sdg, q); }
QuantumCircuit& QuantumCircuit::t(qubit_t q)   { return add_1q(GateType::T, q); }
QuantumCircuit& QuantumCircuit::tdg(qubit_t q) { return add_1q(GateType::Tdg, q); }
QuantumCircuit& QuantumCircuit::sx(qubit_t q)  { return add_1q(GateType::SX, q); }

QuantumCircuit& QuantumCircuit::rx(double theta, qubit_t q) {
    return add_1q(GateType::RX, q, {theta});
}
QuantumCircuit& QuantumCircuit::ry(double theta, qubit_t q) {
    return add_1q(GateType::RY, q, {theta});
}
QuantumCircuit& QuantumCircuit::rz(double theta, qubit_t q) {
    return add_1q(GateType::RZ, q, {theta});
}
QuantumCircuit& QuantumCircuit::p(double lambda, qubit_t q) {
    return add_1q(GateType::P, q, {lambda});
}
QuantumCircuit& QuantumCircuit::u(double theta, double phi, double lambda, qubit_t q) {
    return add_1q(GateType::U, q, {theta, phi, lambda});
}

QuantumCircuit& QuantumCircuit::cx(qubit_t control, qubit_t target) {
    return add_2q(GateType::CX, control, target);
}
QuantumCircuit& QuantumCircuit::cz(qubit_t control, qubit_t target) {
    return add_2q(GateType::CZ, control, target);
}
QuantumCircuit& QuantumCircuit::cp(double lambda, qubit_t control, qubit_t target) {
    return add_2q(GateType::CP, control, target, {lambda});
}
QuantumCircuit& QuantumCircuit::swap(qubit_t a, qubit_t b) {
    return add_2q(GateType::SWAP, a, b);
}

QuantumCircuit& QuantumCircuit::measure(qubit_t q, clbit_t c) {
    check_qubit(q);
    check_clbit(c);
    instructions_.push_back({Gate{GateType::Measure, {}}, {q}, {c}});
    return *this;
}

QuantumCircuit& QuantumCircuit::measure_all() {
    num_clbits_ = std::max(num_clbits_, num_qubits_);
    for (qubit_t q = 0; q < num_qubits_; ++q) {
        measure(q, static_cast<clbit_t>(q));
    }
    return *this;
}

QuantumCircuit& QuantumCircuit::barrier() {
    instructions_.push_back({Gate{GateType::Barrier, {}}, {}, {}});
    return *this;
}

std::size_t QuantumCircuit::depth() const {
    std::vector<std::size_t> level(num_qubits_, 0);
    for (const auto& inst : instructions_) {
        if (inst.gate.type == GateType::Barrier) continue;
        std::size_t d = 0;
        for (qubit_t q : inst.qubits) d = std::max(d, level[q]);
        ++d;
        for (qubit_t q : inst.qubits) level[q] = d;
    }
    return *std::max_element(level.begin(), level.end());
}

}  // namespace lqcs
