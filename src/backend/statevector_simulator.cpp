#include "lqcs/backend/statevector_simulator.hpp"

#include <algorithm>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "lqcs/core/bit_utils.hpp"
#include "lqcs/core/random.hpp"
#include "lqcs/transpiler/gate_fusion.hpp"
#include "kernels/kernels.hpp"

namespace lqcs {

namespace {

void apply_instruction(Statevector& sv, const Instruction& inst) {
    const Gate& g = inst.gate;
    const std::size_t nc = g.n_controls;
    const std::size_t k = g.base_qubits();
    const qubit_t* controls = inst.qubits.data();
    const qubit_t* targets = inst.qubits.data() + nc;

    if (g.type == GateType::Permutation) {
        kernels::apply_controlled_permutation(sv.data(), sv.size(), controls, nc,
                                              targets, k, g.perm.data());
        return;
    }

    // 无控制位的快路径
    if (nc == 0) {
        if (k == 1) {
            const auto m = g.base_matrix();
            if (g.is_diagonal()) {
                kernels::apply_diagonal_single_qubit(sv.data(), sv.size(),
                                                     targets[0], m[0], m[3]);
            } else {
                kernels::apply_single_qubit(sv.data(), sv.size(), targets[0],
                                            m.data());
            }
            return;
        }
        if (g.type == GateType::SWAP) {
            kernels::apply_swap(sv.data(), sv.size(), targets[0], targets[1]);
            return;
        }
        if (k == 2) {
            const auto m = g.base_matrix();
            kernels::apply_two_qubit(sv.data(), sv.size(), targets[0], targets[1],
                                     m.data());
            return;
        }
    }

    // 单控制单比特快路径
    if (nc == 1 && k == 1) {
        const auto m = g.base_matrix();
        kernels::apply_controlled_single_qubit(sv.data(), sv.size(), controls[0],
                                               targets[0], m.data());
        return;
    }

    // 通用路径：任意控制位数 + 任意基门
    const auto m = g.base_matrix();
    kernels::apply_controlled_unitary(sv.data(), sv.size(), controls, nc, targets,
                                      k, m.data());
}

// 末尾测量电路的快速路径：演化一次，对概率分布采样 shots 次
Result run_sampling(const QuantumCircuit& circuit, std::size_t shots,
                    const Statevector& sv, Rng& rng) {
    // 测量映射：clbit <- qubit（同一 clbit 被多次写入时，以最后一次为准）
    std::vector<std::optional<qubit_t>> clbit_source(circuit.num_clbits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            clbit_source[inst.clbits[0]] = inst.qubits[0];
        }
    }

    std::vector<double> cumulative = sv.probabilities();
    std::partial_sum(cumulative.begin(), cumulative.end(), cumulative.begin());
    const double total = cumulative.back();

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

// 中间测量/Reset 电路：每个 shot 独立演化并按 Born 规则坍缩
Result run_per_shot(const QuantumCircuit& circuit, std::size_t shots, Rng& rng) {
    std::map<std::string, std::size_t> counts;
    const std::size_t n_clbits = circuit.num_clbits();
    std::string key(n_clbits, '0');
    for (std::size_t s = 0; s < shots; ++s) {
        Statevector sv(circuit.num_qubits());
        std::vector<bool> creg(n_clbits, false);
        for (const auto& inst : circuit.instructions()) {
            switch (inst.gate.type) {
                case GateType::Measure:
                    creg[inst.clbits[0]] = kernels::collapse_qubit(
                        sv.data(), sv.size(), inst.qubits[0], rng.uniform());
                    break;
                case GateType::Reset:
                    kernels::reset_qubit(sv.data(), sv.size(), inst.qubits[0],
                                         rng.uniform());
                    break;
                case GateType::Barrier:
                    break;
                default:
                    apply_instruction(sv, inst);
                    break;
            }
        }
        for (std::size_t c = 0; c < n_clbits; ++c) {
            key[n_clbits - 1 - c] = creg[c] ? '1' : '0';
        }
        ++counts[key];
    }
    return Result(std::move(counts), shots);
}

}  // namespace

namespace {

// 无融合、无线程设置的纯演化（融合与线程配置由公共入口统一处理）
Statevector evolve(const QuantumCircuit& circuit) {
    Statevector sv(circuit.num_qubits());
    for (const auto& inst : circuit.instructions()) {
        if (inst.gate.type == GateType::Reset) {
            throw std::invalid_argument(
                "run_statevector: Reset is non-deterministic; use run() instead");
        }
        if (inst.gate.is_unitary()) apply_instruction(sv, inst);
    }
    return sv;
}

}  // namespace

void StatevectorSimulator::configure_threads() const {
#ifdef _OPENMP
    if (options_.num_threads > 0) omp_set_num_threads(options_.num_threads);
#endif
}

QuantumCircuit StatevectorSimulator::maybe_fuse(const QuantumCircuit& circuit) const {
    if (!options_.fuse_gates) return circuit;
    return SingleQubitGateFusion().run(circuit);
}

Statevector StatevectorSimulator::run_statevector(const QuantumCircuit& circuit) const {
    configure_threads();
    return evolve(maybe_fuse(circuit));
}

Result StatevectorSimulator::run(const QuantumCircuit& circuit, std::size_t shots) {
    if (shots == 0) {
        throw std::invalid_argument("run: shots must be >= 1");
    }
    configure_threads();
    const QuantumCircuit opt = maybe_fuse(circuit);

    // 自动选择执行路径：含 Reset 或「测量后还有门」时走逐 shot 演化
    bool seen_measure = false;
    bool per_shot = false;
    for (const auto& inst : opt.instructions()) {
        if (inst.gate.type == GateType::Measure) {
            seen_measure = true;
        } else if (inst.gate.type == GateType::Reset) {
            per_shot = true;
        } else if (seen_measure && inst.gate.is_unitary()) {
            per_shot = true;
        }
    }
    if (!seen_measure) {
        throw std::invalid_argument(
            "run: circuit has no measurements; use run_statevector() to inspect "
            "the final state");
    }

    Rng rng(options_.seed);
    if (per_shot) {
        return run_per_shot(opt, shots, rng);
    }
    const Statevector sv = evolve(opt);
    return run_sampling(opt, shots, sv, rng);
}

}  // namespace lqcs
