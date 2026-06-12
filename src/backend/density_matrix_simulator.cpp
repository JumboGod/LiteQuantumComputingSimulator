#include "lqcs/backend/density_matrix_simulator.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "lqcs/core/bit_utils.hpp"
#include "lqcs/core/random.hpp"

namespace lqcs {

DensityMatrix DensityMatrixSimulator::run_density_matrix(
    const QuantumCircuit& circuit) const {
    DensityMatrix rho(circuit.num_qubits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Reset) {
            throw std::invalid_argument(
                "DensityMatrixSimulator: Reset is not supported yet");
        }
        if (!inst.gate.is_unitary()) continue;  // Measure/Barrier
        rho.apply_instruction(inst);
        for (const auto& channel : options_.noise.channels()) {
            for (qubit_t q : inst.qubits) rho.apply_channel(channel, q);
        }
    }
    return rho;
}

Result DensityMatrixSimulator::run(const QuantumCircuit& circuit,
                                   std::size_t shots) {
    if (shots == 0) {
        throw std::invalid_argument("run: shots must be >= 1");
    }
    // 仅支持末尾测量
    bool seen_measure = false;
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            seen_measure = true;
        } else if (seen_measure && inst.gate.is_unitary()) {
            throw std::invalid_argument(
                "DensityMatrixSimulator: mid-circuit measurement is not "
                "supported");
        }
    }
    if (!seen_measure) {
        throw std::invalid_argument(
            "run: circuit has no measurements; use run_density_matrix()");
    }

    const DensityMatrix rho = run_density_matrix(circuit);

    std::vector<std::optional<qubit_t>> clbit_source(circuit.num_clbits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            clbit_source[inst.clbits[0]] = inst.qubits[0];
        }
    }

    // 对角元即测量分布，与状态向量路径相同的累积采样
    std::vector<double> cumulative = rho.probabilities();
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
