#include "lqcs/backend/stabilizer_simulator.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "lqcs/core/random.hpp"

namespace lqcs {

StabilizerTableau StabilizerSimulator::run_tableau(
    const QuantumCircuit& circuit) const {
    Rng rng(options_.seed);
    StabilizerTableau tab(circuit.num_qubits());
    for (const auto& inst : circuit.instructions()) {
        switch (inst.gate.type) {
            case GateType::Barrier:
                break;
            case GateType::Measure:
                tab.measure(inst.qubits[0], rng.uniform());
                break;
            case GateType::Reset:
                tab.reset(inst.qubits[0], rng.uniform());
                break;
            default:
                tab.apply(inst);
                break;
        }
    }
    return tab;
}

Result StabilizerSimulator::run(const QuantumCircuit& circuit, std::size_t shots) {
    if (shots == 0) {
        throw std::invalid_argument("run: shots must be >= 1");
    }
    // 预检查门集（在采样前给出清晰报错）
    for (const auto& inst : circuit.instructions()) {
        const auto t = inst.gate.type;
        if (t != GateType::Measure && t != GateType::Reset &&
            t != GateType::Barrier && !StabilizerTableau::is_clifford(inst.gate)) {
            throw std::invalid_argument(
                "StabilizerSimulator: non-Clifford gate '" + inst.gate.name() +
                "'");
        }
    }

    bool has_measure = false;
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) { has_measure = true; break; }
    }
    if (!has_measure) {
        throw std::invalid_argument(
            "run: circuit has no measurements; use run_tableau() to inspect "
            "the stabilizer state");
    }

    Rng rng(options_.seed);
    const std::size_t n_clbits = circuit.num_clbits();
    std::map<std::string, std::size_t> counts;
    std::string key(n_clbits, '0');

    for (std::size_t shot = 0; shot < shots; ++shot) {
        StabilizerTableau tab(circuit.num_qubits());
        std::vector<bool> creg(n_clbits, false);
        for (const auto& inst : circuit.instructions()) {
            switch (inst.gate.type) {
                case GateType::Barrier:
                    break;
                case GateType::Measure:
                    creg[inst.clbits[0]] =
                        tab.measure(inst.qubits[0], rng.uniform());
                    break;
                case GateType::Reset:
                    tab.reset(inst.qubits[0], rng.uniform());
                    break;
                default:
                    tab.apply(inst);
                    break;
            }
        }
        for (std::size_t c = 0; c < n_clbits; ++c) {
            key[n_clbits - 1 - c] = creg[c] ? '1' : '0';  // 最高位 clbit 在左
        }
        ++counts[key];
    }
    return Result(std::move(counts), shots);
}

}  // namespace lqcs
