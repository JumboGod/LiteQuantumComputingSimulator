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

    // 多比特 Pauli 旋转（无控制位走专用 O(2^n) 内核）
    if (g.type == GateType::PauliRotation && nc == 0) {
        const auto pm = g.pauli_masks();
        std::size_t x_global = 0, zy_global = 0;
        for (std::size_t j = 0; j < k; ++j) {
            if (pm.x & (std::size_t{1} << j)) x_global |= std::size_t{1} << targets[j];
            if (pm.zy & (std::size_t{1} << j)) zy_global |= std::size_t{1} << targets[j];
        }
        kernels::apply_pauli_rotation(sv.data(), sv.size(), x_global, zy_global,
                                      pm.n_y, g.params.at(0));
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

// 轨迹法噪声：对一个比特按 ||K_i|ψ>||² 采样并应用一个 Kraus 分支
void sample_kraus_branch(Statevector& sv, const KrausChannel& channel,
                         qubit_t q, Rng& rng) {
    const double r = rng.uniform();
    double cum = 0.0;
    // 兜底：浮点残差导致累计概率略小于 1 时，落到最后一个非零分支
    const std::array<complex_t, 4>* chosen = nullptr;
    double chosen_p = 0.0;
    for (const auto& kraus : channel.ops) {
        const double p =
            kernels::kraus_probability(sv.data(), sv.size(), q, kraus.data());
        if (p <= 0.0) continue;
        chosen = &kraus;
        chosen_p = p;
        cum += p;
        if (r < cum) break;
    }
    if (!chosen) return;  // 理论上不可达（完备性保证 Σp = 1）
    // 归一化 1/√p 融进矩阵，复用单比特内核
    const double inv = 1.0 / std::sqrt(chosen_p);
    const complex_t m[4] = {(*chosen)[0] * inv, (*chosen)[1] * inv,
                            (*chosen)[2] * inv, (*chosen)[3] * inv};
    kernels::apply_single_qubit(sv.data(), sv.size(), q, m);
}

void apply_trajectory_noise(Statevector& sv, const std::vector<qubit_t>& qubits,
                            const NoiseModel& noise, Rng& rng) {
    for (const auto& channel : noise.channels()) {
        for (qubit_t q : qubits) {
            sample_kraus_branch(sv, channel, q, rng);
        }
    }
}

// 中间测量/Reset/含噪电路：每个 shot 独立演化（轨迹）
Result run_per_shot(const QuantumCircuit& circuit, std::size_t shots, Rng& rng,
                    const NoiseModel& noise) {
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
                    if (!noise.empty()) {
                        apply_trajectory_noise(sv, inst.qubits, noise, rng);
                    }
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
    if (!options_.noise.empty()) {
        throw std::invalid_argument(
            "run_statevector: noisy simulation is stochastic (one trajectory is "
            "not the quantum state); use run() or DensityMatrixSimulator");
    }
    configure_threads();
    return evolve(maybe_fuse(circuit));
}

Result StatevectorSimulator::run(const QuantumCircuit& circuit, std::size_t shots) {
    if (shots == 0) {
        throw std::invalid_argument("run: shots must be >= 1");
    }
    configure_threads();
    const bool noisy = !options_.noise.empty();
    // 噪声绑定在门上：融合会合并/消去门，改变噪声作用点，故含噪时禁用
    const QuantumCircuit opt = noisy ? circuit : maybe_fuse(circuit);

    // 自动选择执行路径：含噪声、Reset 或「测量后还有门」时走逐 shot 演化
    bool seen_measure = false;
    bool per_shot = noisy;
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
        return run_per_shot(opt, shots, rng, options_.noise);
    }
    const Statevector sv = evolve(opt);
    return run_sampling(opt, shots, sv, rng);
}

}  // namespace lqcs
