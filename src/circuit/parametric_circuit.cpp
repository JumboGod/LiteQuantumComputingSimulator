#include "lqcs/circuit/parametric_circuit.hpp"

#include <stdexcept>

namespace lqcs {

ParametricCircuit::ParametricCircuit(std::size_t num_qubits, std::size_t num_clbits)
    : template_(num_qubits, num_clbits) {}

ParametricCircuit& ParametricCircuit::x(qubit_t q) { template_.x(q); return *this; }
ParametricCircuit& ParametricCircuit::y(qubit_t q) { template_.y(q); return *this; }
ParametricCircuit& ParametricCircuit::z(qubit_t q) { template_.z(q); return *this; }
ParametricCircuit& ParametricCircuit::h(qubit_t q) { template_.h(q); return *this; }
ParametricCircuit& ParametricCircuit::s(qubit_t q) { template_.s(q); return *this; }
ParametricCircuit& ParametricCircuit::sdg(qubit_t q) { template_.sdg(q); return *this; }
ParametricCircuit& ParametricCircuit::t(qubit_t q) { template_.t(q); return *this; }
ParametricCircuit& ParametricCircuit::sx(qubit_t q) { template_.sx(q); return *this; }
ParametricCircuit& ParametricCircuit::cx(qubit_t c, qubit_t t) {
    template_.cx(c, t); return *this;
}
ParametricCircuit& ParametricCircuit::cy(qubit_t c, qubit_t t) {
    template_.cy(c, t); return *this;
}
ParametricCircuit& ParametricCircuit::cz(qubit_t c, qubit_t t) {
    template_.cz(c, t); return *this;
}
ParametricCircuit& ParametricCircuit::swap(qubit_t a, qubit_t b) {
    template_.swap(a, b); return *this;
}
ParametricCircuit& ParametricCircuit::ccx(qubit_t c1, qubit_t c2, qubit_t t) {
    template_.ccx(c1, c2, t); return *this;
}
ParametricCircuit& ParametricCircuit::barrier() { template_.barrier(); return *this; }

std::size_t ParametricCircuit::add_param_gate(Gate gate,
                                              std::vector<qubit_t> qubits) {
    const std::size_t inst_idx = template_.num_instructions();
    template_.append({std::move(gate), std::move(qubits), {}});
    param_inst_.push_back(inst_idx);
    return param_inst_.size() - 1;
}

std::size_t ParametricCircuit::rx(qubit_t q) {
    return add_param_gate({GateType::RX, {0.0}}, {q});
}
std::size_t ParametricCircuit::ry(qubit_t q) {
    return add_param_gate({GateType::RY, {0.0}}, {q});
}
std::size_t ParametricCircuit::rz(qubit_t q) {
    return add_param_gate({GateType::RZ, {0.0}}, {q});
}

std::size_t ParametricCircuit::rp(std::string_view pauli,
                                  std::span<const qubit_t> qubits) {
    if (pauli.empty() || pauli.size() != qubits.size()) {
        throw std::invalid_argument(
            "ParametricCircuit::rp: pauli length must equal qubit count");
    }
    for (char c : pauli) {
        if (c != 'I' && c != 'X' && c != 'Y' && c != 'Z') {
            throw std::invalid_argument("rp: pauli must contain only I/X/Y/Z");
        }
    }
    Gate g{GateType::PauliRotation, {0.0}, 0, {}, {}, std::string(pauli)};
    return add_param_gate(std::move(g),
                          std::vector<qubit_t>(qubits.begin(), qubits.end()));
}
std::size_t ParametricCircuit::rp(std::string_view pauli,
                                  std::initializer_list<qubit_t> qubits) {
    return rp(pauli, std::span<const qubit_t>(qubits.begin(), qubits.size()));
}

QuantumCircuit ParametricCircuit::bind(std::span<const double> values) const {
    if (values.size() != num_parameters()) {
        throw std::invalid_argument(
            "bind: number of values must equal num_parameters()");
    }
    // 参数指令下标 → 参数 id（反查表）
    std::vector<std::size_t> param_of(template_.num_instructions(),
                                      static_cast<std::size_t>(-1));
    for (std::size_t p = 0; p < param_inst_.size(); ++p) {
        param_of[param_inst_[p]] = p;
    }

    QuantumCircuit out(template_.num_qubits(), template_.num_clbits());
    std::size_t i = 0;
    for (const auto& inst : template_.instructions()) {
        Instruction copy = inst;
        if (param_of[i] != static_cast<std::size_t>(-1)) {
            copy.gate.params[0] = values[param_of[i]];
        }
        out.append(std::move(copy));
        ++i;
    }
    return out;
}

}  // namespace lqcs
