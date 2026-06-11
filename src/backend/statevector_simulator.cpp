#include "lqcs/backend/statevector_simulator.hpp"

#include <algorithm>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

#include "lqcs/core/bit_utils.hpp"
#include "lqcs/core/random.hpp"
#include "kernels/kernels.hpp"

namespace lqcs {

namespace {

void apply_instruction(Statevector& sv, const Instruction& inst) {
    const Gate& g = inst.gate;
    switch (g.type) {
        case GateType::CX:
        case GateType::CZ:
        case GateType::CP: {
            // 受控门走子空间路径：对 control=1 的子空间应用基门矩阵
            const GateType base = (g.type == GateType::CX)   ? GateType::X
                                  : (g.type == GateType::CZ) ? GateType::Z
                                                             : GateType::P;
            const auto m = Gate{base, g.params}.matrix();
            kernels::apply_controlled_single_qubit(sv.data(), sv.size(),
                                                   inst.qubits[0], inst.qubits[1],
                                                   m.data());
            break;
        }
        case GateType::SWAP:
            kernels::apply_swap(sv.data(), sv.size(), inst.qubits[0], inst.qubits[1]);
            break;
        default: {
            const auto m = g.matrix();
            kernels::apply_single_qubit(sv.data(), sv.size(), inst.qubits[0],
                                        m.data());
            break;
        }
    }
}

}  // namespace

Statevector StatevectorSimulator::run_statevector(const QuantumCircuit& circuit) const {
    Statevector sv(circuit.num_qubits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.is_unitary()) apply_instruction(sv, inst);
    }
    return sv;
}

Result StatevectorSimulator::run(const QuantumCircuit& circuit, std::size_t shots) {
    if (shots == 0) {
        throw std::invalid_argument("run: shots must be >= 1");
    }

    // M1 限制：只支持末尾测量。校验测量之后不再出现幺正门。
    bool seen_measure = false;
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            seen_measure = true;
        } else if (seen_measure && inst.gate.is_unitary()) {
            throw std::runtime_error(
                "mid-circuit measurement is not supported yet (planned for M2): "
                "gates found after a measure");
        }
    }
    if (!seen_measure) {
        throw std::invalid_argument(
            "run: circuit has no measurements; use run_statevector() to inspect "
            "the final state");
    }

    // 演化一次到末态
    Statevector sv = run_statevector(circuit);

    // 测量映射：clbit <- qubit（同一 clbit 被多次写入时，以最后一次为准）
    std::vector<std::optional<qubit_t>> clbit_source(circuit.num_clbits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            clbit_source[inst.clbits[0]] = inst.qubits[0];
        }
    }

    // 对概率分布累积求和后做二分采样
    std::vector<double> cumulative = sv.probabilities();
    std::partial_sum(cumulative.begin(), cumulative.end(), cumulative.begin());
    const double total = cumulative.back();

    Rng rng(options_.seed);
    std::map<std::string, std::size_t> counts;
    const std::size_t n_clbits = circuit.num_clbits();
    std::string key(n_clbits, '0');
    for (std::size_t s = 0; s < shots; ++s) {
        const double r = rng.uniform() * total;
        const std::size_t idx = static_cast<std::size_t>(
            std::upper_bound(cumulative.begin(), cumulative.end(), r) -
            cumulative.begin());
        // 比特串与 Qiskit 一致：最高位经典比特在左
        for (std::size_t c = 0; c < n_clbits; ++c) {
            const bool bit =
                clbit_source[c] && bits::test_bit(idx, *clbit_source[c]);
            key[n_clbits - 1 - c] = bit ? '1' : '0';
        }
        ++counts[key];
    }
    return Result(std::move(counts), shots);
}

}  // namespace lqcs
