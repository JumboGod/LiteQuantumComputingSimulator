#include "lqcs/circuit/quantum_circuit.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace lqcs {

namespace {

// 幺正性校验：U·U† = I，容差 1e-10
void validate_unitary(std::span<const complex_t> m, std::size_t dim) {
    constexpr double kTol = 1e-10;
    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) {
            complex_t sum{0, 0};
            for (std::size_t l = 0; l < dim; ++l) {
                sum += m[i * dim + l] * std::conj(m[j * dim + l]);
            }
            const complex_t expected = (i == j) ? complex_t{1, 0} : complex_t{0, 0};
            if (std::abs(sum - expected) > kTol) {
                throw std::invalid_argument(
                    "unitary: matrix is not unitary (U*Udag != I)");
            }
        }
    }
}

bool is_power_of_two(std::size_t v) { return v != 0 && (v & (v - 1)) == 0; }

}  // namespace

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

QuantumCircuit& QuantumCircuit::add_gate(Gate gate, std::vector<qubit_t> qubits) {
    if (qubits.size() != gate.num_qubits()) {
        throw std::invalid_argument(gate.name() + ": expects " +
                                    std::to_string(gate.num_qubits()) +
                                    " qubits, got " + std::to_string(qubits.size()));
    }
    std::unordered_set<qubit_t> seen;
    for (qubit_t q : qubits) {
        check_qubit(q);
        if (!seen.insert(q).second) {
            throw std::invalid_argument(gate.name() + ": duplicate qubit " +
                                        std::to_string(q));
        }
    }
    instructions_.push_back({std::move(gate), std::move(qubits), {}});
    return *this;
}

QuantumCircuit& QuantumCircuit::append(Instruction inst) {
    if (inst.gate.type == GateType::Measure) {
        if (inst.qubits.size() != 1 || inst.clbits.size() != 1) {
            throw std::invalid_argument("append: malformed measure instruction");
        }
        check_clbit(inst.clbits[0]);
        check_qubit(inst.qubits[0]);
        instructions_.push_back(std::move(inst));
        return *this;
    }
    return add_gate(std::move(inst.gate), std::move(inst.qubits));
}

QuantumCircuit& QuantumCircuit::i(qubit_t q)   { return add_gate({GateType::I, {}}, {q}); }
QuantumCircuit& QuantumCircuit::x(qubit_t q)   { return add_gate({GateType::X, {}}, {q}); }
QuantumCircuit& QuantumCircuit::y(qubit_t q)   { return add_gate({GateType::Y, {}}, {q}); }
QuantumCircuit& QuantumCircuit::z(qubit_t q)   { return add_gate({GateType::Z, {}}, {q}); }
QuantumCircuit& QuantumCircuit::h(qubit_t q)   { return add_gate({GateType::H, {}}, {q}); }
QuantumCircuit& QuantumCircuit::s(qubit_t q)   { return add_gate({GateType::S, {}}, {q}); }
QuantumCircuit& QuantumCircuit::sdg(qubit_t q) { return add_gate({GateType::Sdg, {}}, {q}); }
QuantumCircuit& QuantumCircuit::t(qubit_t q)   { return add_gate({GateType::T, {}}, {q}); }
QuantumCircuit& QuantumCircuit::tdg(qubit_t q) { return add_gate({GateType::Tdg, {}}, {q}); }
QuantumCircuit& QuantumCircuit::sx(qubit_t q)  { return add_gate({GateType::SX, {}}, {q}); }

QuantumCircuit& QuantumCircuit::rx(double theta, qubit_t q) {
    return add_gate({GateType::RX, {theta}}, {q});
}
QuantumCircuit& QuantumCircuit::ry(double theta, qubit_t q) {
    return add_gate({GateType::RY, {theta}}, {q});
}
QuantumCircuit& QuantumCircuit::rz(double theta, qubit_t q) {
    return add_gate({GateType::RZ, {theta}}, {q});
}
QuantumCircuit& QuantumCircuit::p(double lambda, qubit_t q) {
    return add_gate({GateType::P, {lambda}}, {q});
}
QuantumCircuit& QuantumCircuit::u(double theta, double phi, double lambda, qubit_t q) {
    return add_gate({GateType::U, {theta, phi, lambda}}, {q});
}

QuantumCircuit& QuantumCircuit::cx(qubit_t c, qubit_t t) {
    return add_gate({GateType::X, {}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::cy(qubit_t c, qubit_t t) {
    return add_gate({GateType::Y, {}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::cz(qubit_t c, qubit_t t) {
    return add_gate({GateType::Z, {}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::ch(qubit_t c, qubit_t t) {
    return add_gate({GateType::H, {}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::cp(double lambda, qubit_t c, qubit_t t) {
    return add_gate({GateType::P, {lambda}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::crx(double theta, qubit_t c, qubit_t t) {
    return add_gate({GateType::RX, {theta}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::cry(double theta, qubit_t c, qubit_t t) {
    return add_gate({GateType::RY, {theta}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::crz(double theta, qubit_t c, qubit_t t) {
    return add_gate({GateType::RZ, {theta}, 1}, {c, t});
}
QuantumCircuit& QuantumCircuit::ccx(qubit_t c1, qubit_t c2, qubit_t t) {
    return add_gate({GateType::X, {}, 2}, {c1, c2, t});
}
QuantumCircuit& QuantumCircuit::cswap(qubit_t c, qubit_t a, qubit_t b) {
    return add_gate({GateType::SWAP, {}, 1}, {c, a, b});
}

QuantumCircuit& QuantumCircuit::mcx(std::span<const qubit_t> controls, qubit_t t) {
    if (controls.empty()) {
        throw std::invalid_argument("mcx: needs at least one control qubit");
    }
    std::vector<qubit_t> qs(controls.begin(), controls.end());
    qs.push_back(t);
    return add_gate({GateType::X, {}, controls.size()}, std::move(qs));
}
QuantumCircuit& QuantumCircuit::mcx(std::initializer_list<qubit_t> controls,
                                    qubit_t t) {
    return mcx(std::span<const qubit_t>(controls.begin(), controls.size()), t);
}
QuantumCircuit& QuantumCircuit::mcp(double lambda,
                                    std::span<const qubit_t> controls, qubit_t t) {
    if (controls.empty()) {
        throw std::invalid_argument("mcp: needs at least one control qubit");
    }
    std::vector<qubit_t> qs(controls.begin(), controls.end());
    qs.push_back(t);
    return add_gate({GateType::P, {lambda}, controls.size()}, std::move(qs));
}
QuantumCircuit& QuantumCircuit::mcp(double lambda,
                                    std::initializer_list<qubit_t> controls,
                                    qubit_t t) {
    return mcp(lambda, std::span<const qubit_t>(controls.begin(), controls.size()), t);
}

QuantumCircuit& QuantumCircuit::swap(qubit_t a, qubit_t b) {
    return add_gate({GateType::SWAP, {}}, {a, b});
}
QuantumCircuit& QuantumCircuit::iswap(qubit_t a, qubit_t b) {
    return add_gate({GateType::iSWAP, {}}, {a, b});
}
QuantumCircuit& QuantumCircuit::rxx(double theta, qubit_t a, qubit_t b) {
    return add_gate({GateType::RXX, {theta}}, {a, b});
}
QuantumCircuit& QuantumCircuit::ryy(double theta, qubit_t a, qubit_t b) {
    return add_gate({GateType::RYY, {theta}}, {a, b});
}
QuantumCircuit& QuantumCircuit::rzz(double theta, qubit_t a, qubit_t b) {
    return add_gate({GateType::RZZ, {theta}}, {a, b});
}

QuantumCircuit& QuantumCircuit::rp(double theta, std::string_view pauli,
                                   std::span<const qubit_t> qubits) {
    if (pauli.empty() || pauli.size() != qubits.size()) {
        throw std::invalid_argument(
            "rp: pauli length must equal qubit count and be non-empty");
    }
    for (char c : pauli) {
        if (c != 'I' && c != 'X' && c != 'Y' && c != 'Z') {
            throw std::invalid_argument("rp: pauli must contain only I/X/Y/Z");
        }
    }
    Gate g{GateType::PauliRotation, {theta}, 0, {}, {}, std::string(pauli)};
    return add_gate(std::move(g), std::vector<qubit_t>(qubits.begin(), qubits.end()));
}
QuantumCircuit& QuantumCircuit::rp(double theta, std::string_view pauli,
                                   std::initializer_list<qubit_t> qubits) {
    return rp(theta, pauli, std::span<const qubit_t>(qubits.begin(), qubits.size()));
}

QuantumCircuit& QuantumCircuit::unitary(std::span<const complex_t> matrix,
                                        std::span<const qubit_t> qubits) {
    const std::size_t dim = std::size_t{1} << qubits.size();
    if (matrix.size() != dim * dim) {
        throw std::invalid_argument(
            "unitary: matrix size must be 2^k x 2^k for k qubits");
    }
    validate_unitary(matrix, dim);
    Gate g{GateType::Unitary, {}, 0,
           std::vector<complex_t>(matrix.begin(), matrix.end()), {}};
    return add_gate(std::move(g), std::vector<qubit_t>(qubits.begin(), qubits.end()));
}
QuantumCircuit& QuantumCircuit::unitary(std::span<const complex_t> matrix,
                                        std::initializer_list<qubit_t> qubits) {
    return unitary(matrix, std::span<const qubit_t>(qubits.begin(), qubits.size()));
}

QuantumCircuit& QuantumCircuit::permutation(std::span<const std::size_t> table,
                                            std::span<const qubit_t> qubits) {
    const std::size_t dim = std::size_t{1} << qubits.size();
    if (table.size() != dim || !is_power_of_two(table.size())) {
        throw std::invalid_argument("permutation: table size must be 2^k for k qubits");
    }
    std::vector<bool> hit(dim, false);
    for (std::size_t v : table) {
        if (v >= dim || hit[v]) {
            throw std::invalid_argument("permutation: table is not a bijection");
        }
        hit[v] = true;
    }
    Gate g{GateType::Permutation, {}, 0, {},
           std::vector<std::size_t>(table.begin(), table.end())};
    return add_gate(std::move(g), std::vector<qubit_t>(qubits.begin(), qubits.end()));
}
QuantumCircuit& QuantumCircuit::permutation(std::span<const std::size_t> table,
                                            std::initializer_list<qubit_t> qubits) {
    return permutation(table, std::span<const qubit_t>(qubits.begin(), qubits.size()));
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

QuantumCircuit& QuantumCircuit::reset(qubit_t q) {
    check_qubit(q);
    instructions_.push_back({Gate{GateType::Reset, {}}, {q}, {}});
    return *this;
}

QuantumCircuit& QuantumCircuit::barrier() {
    instructions_.push_back({Gate{GateType::Barrier, {}}, {}, {}});
    return *this;
}

QuantumCircuit& QuantumCircuit::compose(const QuantumCircuit& other,
                                        std::span<const qubit_t> qubit_map) {
    std::vector<qubit_t> map(qubit_map.begin(), qubit_map.end());
    if (map.empty()) {
        if (other.num_qubits() > num_qubits_) {
            throw std::invalid_argument("compose: sub-circuit has more qubits");
        }
        map.resize(other.num_qubits());
        for (std::size_t q = 0; q < map.size(); ++q) map[q] = static_cast<qubit_t>(q);
    }
    if (map.size() != other.num_qubits()) {
        throw std::invalid_argument(
            "compose: qubit_map size must equal sub-circuit qubit count");
    }
    std::unordered_set<qubit_t> seen;
    for (qubit_t q : map) {
        check_qubit(q);
        if (!seen.insert(q).second) {
            throw std::invalid_argument("compose: duplicate qubit in map");
        }
    }
    if (other.num_clbits() > num_clbits_) {
        throw std::invalid_argument("compose: sub-circuit has more clbits");
    }
    for (const auto& inst : other.instructions()) {
        Instruction mapped = inst;
        for (auto& q : mapped.qubits) q = map[q];
        instructions_.push_back(std::move(mapped));
    }
    return *this;
}

QuantumCircuit QuantumCircuit::inverse() const {
    QuantumCircuit inv(num_qubits_, num_clbits_);
    for (auto it = instructions_.rbegin(); it != instructions_.rend(); ++it) {
        if (it->gate.type == GateType::Barrier) {
            inv.barrier();
            continue;
        }
        if (!it->gate.is_unitary()) {
            throw std::invalid_argument(
                "inverse: circuit contains non-unitary operation '" +
                it->gate.name() + "'");
        }
        inv.instructions_.push_back({it->gate.inverse(), it->qubits, {}});
    }
    return inv;
}

QuantumCircuit QuantumCircuit::power(int n) const {
    QuantumCircuit result(num_qubits_, num_clbits_);
    const QuantumCircuit base = (n < 0) ? inverse() : *this;
    const int reps = (n < 0) ? -n : n;
    for (int r = 0; r < reps; ++r) result.compose(base);
    return result;
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
