#include "lqcs/transpiler/redundant_elimination.hpp"

#include <cmath>
#include <vector>

namespace lqcs {

namespace {

bool params_negated(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == 1 && b.size() == 1 && a[0] == -b[0];
}

// b 是否为 a 的逆（结构化判断，避免大矩阵乘法）
bool is_inverse_pair(const Gate& a, const Gate& b) {
    if (a.n_controls != b.n_controls) return false;
    switch (a.type) {
        // 自逆门
        case GateType::I:
        case GateType::X:
        case GateType::Y:
        case GateType::Z:
        case GateType::H:
        case GateType::SWAP:
            return b.type == a.type;
        // 参数互为相反数
        case GateType::RX:
        case GateType::RY:
        case GateType::RZ:
        case GateType::P:
        case GateType::RXX:
        case GateType::RYY:
        case GateType::RZZ:
            return b.type == a.type && params_negated(a.params, b.params);
        case GateType::S:   return b.type == GateType::Sdg;
        case GateType::Sdg: return b.type == GateType::S;
        case GateType::T:   return b.type == GateType::Tdg;
        case GateType::Tdg: return b.type == GateType::T;
        case GateType::U:
            return b.type == GateType::U && a.params.size() == 3 &&
                   b.params.size() == 3 && b.params[0] == -a.params[0] &&
                   b.params[1] == -a.params[2] && b.params[2] == -a.params[1];
        case GateType::Permutation: {
            if (b.type != GateType::Permutation ||
                a.perm.size() != b.perm.size()) {
                return false;
            }
            for (std::size_t l = 0; l < a.perm.size(); ++l) {
                if (b.perm[a.perm[l]] != l) return false;
            }
            return true;
        }
        default:
            return false;  // SX/iSWAP/Unitary：留给门融合处理
    }
}

}  // namespace

QuantumCircuit RedundantGateElimination::run(const QuantumCircuit& circuit) const {
    std::vector<Instruction> insts(circuit.instructions().begin(),
                                   circuit.instructions().end());

    // 迭代到不动点，处理级联消除（X H H X → X X → 空）
    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<bool> removed(insts.size(), false);
        // 每个比特上最近一条未删除指令的下标（-1 = 无）
        std::vector<std::ptrdiff_t> last(circuit.num_qubits(), -1);

        for (std::size_t i = 0; i < insts.size(); ++i) {
            const Instruction& cur = insts[i];
            if (cur.gate.type == GateType::Barrier) {
                for (auto& l : last) l = -1;  // 屏障阻断跨越优化
                continue;
            }
            // 候选：所有触及比特的「上一条指令」是同一条，且比特列表一致
            bool cancelled = false;
            if (cur.gate.is_unitary() && !cur.qubits.empty()) {
                const std::ptrdiff_t j = last[cur.qubits[0]];
                if (j >= 0 && !removed[j] && insts[j].qubits == cur.qubits) {
                    bool all_match = true;
                    for (qubit_t q : cur.qubits) {
                        if (last[q] != j) { all_match = false; break; }
                    }
                    if (all_match && is_inverse_pair(insts[j].gate, cur.gate)) {
                        removed[i] = removed[j] = true;
                        for (qubit_t q : cur.qubits) last[q] = -1;
                        cancelled = true;
                        changed = true;
                    }
                }
            }
            if (!cancelled) {
                for (qubit_t q : cur.qubits) last[q] = static_cast<std::ptrdiff_t>(i);
            }
        }

        if (changed) {
            std::vector<Instruction> kept;
            kept.reserve(insts.size());
            for (std::size_t i = 0; i < insts.size(); ++i) {
                if (!removed[i]) kept.push_back(std::move(insts[i]));
            }
            insts = std::move(kept);
        }
    }

    QuantumCircuit out(circuit.num_qubits(), circuit.num_clbits());
    for (auto& inst : insts) out.append(std::move(inst));
    return out;
}

}  // namespace lqcs
