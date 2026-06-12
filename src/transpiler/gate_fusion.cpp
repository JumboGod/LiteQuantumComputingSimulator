#include "lqcs/transpiler/gate_fusion.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <vector>

namespace lqcs {

namespace {

using Mat2 = std::array<complex_t, 4>;

Mat2 matmul2(const Mat2& a, const Mat2& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

bool is_identity2(const Mat2& m) {
    constexpr double kTol = 1e-12;
    return std::abs(m[0] - complex_t{1, 0}) < kTol && std::abs(m[1]) < kTol &&
           std::abs(m[2]) < kTol && std::abs(m[3] - complex_t{1, 0}) < kTol;
}

bool is_fusable(const Instruction& inst) {
    return inst.gate.is_unitary() && inst.gate.n_controls == 0 &&
           inst.gate.base_qubits() == 1;
}

}  // namespace

QuantumCircuit SingleQubitGateFusion::run(const QuantumCircuit& circuit) const {
    QuantumCircuit out(circuit.num_qubits(), circuit.num_clbits());
    // 每个比特上待合并的 2×2 矩阵累积
    std::vector<std::optional<Mat2>> pending(circuit.num_qubits());

    auto flush = [&](qubit_t q) {
        if (!pending[q]) return;
        if (!is_identity2(*pending[q])) {
            out.append({Gate{GateType::Unitary, {}, 0,
                             {pending[q]->begin(), pending[q]->end()}},
                        {q},
                        {}});
        }
        pending[q].reset();
    };

    for (const auto& inst : circuit.instructions()) {
        if (is_fusable(inst)) {
            const qubit_t q = inst.qubits[0];
            const auto m = inst.gate.base_matrix();
            const Mat2 g = {m[0], m[1], m[2], m[3]};
            // 后作用的门左乘
            pending[q] = pending[q] ? matmul2(g, *pending[q]) : g;
            continue;
        }
        // Barrier 阻断所有比特的融合；其余指令只阻断其触及的比特
        if (inst.gate.type == GateType::Barrier) {
            for (qubit_t q = 0; q < circuit.num_qubits(); ++q) flush(q);
        } else {
            for (qubit_t q : inst.qubits) flush(q);
        }
        out.append(inst);
    }
    for (qubit_t q = 0; q < circuit.num_qubits(); ++q) flush(q);
    return out;
}

}  // namespace lqcs
