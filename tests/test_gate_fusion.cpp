// 多比特门融合（GateFusion）：等价性、块约束、与单比特融合的关系
#include <cmath>
#include <random>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>
#include <lqcs/transpiler/gate_fusion.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-11;

void expect_same_state(const QuantumCircuit& a, const QuantumCircuit& b) {
    StatevectorSimulator sim({.fuse_gates = false});
    const auto sa = sim.run_statevector(a);
    const auto sb = sim.run_statevector(b);
    ASSERT_EQ(sa.size(), sb.size());
    for (std::size_t i = 0; i < sa.size(); ++i) {
        EXPECT_NEAR(std::abs(sa[i] - sb[i]), 0.0, kTol) << "amplitude " << i;
    }
}

QuantumCircuit random_circuit(std::size_t n, std::size_t n_gates,
                              std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle(0, 6.28);
    QuantumCircuit qc(n);
    for (std::size_t i = 0; i < n_gates; ++i) {
        const auto q = static_cast<qubit_t>(rng() % n);
        switch (rng() % 8) {
            case 0: qc.h(q); break;
            case 1: qc.t(q); break;
            case 2: qc.rx(angle(rng), q); break;
            case 3: qc.u(angle(rng), angle(rng), angle(rng), q); break;
            case 4: qc.sx(q); break;
            case 5: {
                auto b = static_cast<qubit_t>((q + 1) % n);
                qc.cx(q, b);
            } break;
            case 6: {
                auto b = static_cast<qubit_t>((q + 1) % n);
                qc.cp(angle(rng), q, b);
            } break;
            default: {
                auto b = static_cast<qubit_t>((q + 1) % n);
                auto c = static_cast<qubit_t>((q + 2) % n);
                if (q != b && b != c && q != c)
                    qc.ccx(q, b, c);
                else
                    qc.h(q);
                break;
            }
        }
    }
    return qc;
}

}  // namespace

TEST(GateFusion, PreservesStateAcrossBlockSizes) {
    for (std::size_t max_block : {1u, 2u, 3u, 4u}) {
        for (std::uint64_t seed : {1u, 2u, 3u}) {
            const auto qc = random_circuit(5, 80, seed);
            const auto fused = GateFusion(max_block).run(qc);
            SCOPED_TRACE("max_block " + std::to_string(max_block) + " seed " +
                         std::to_string(seed));
            expect_same_state(qc, fused);
        }
    }
}

TEST(GateFusion, ReducesInstructionCount) {
    // 单比特门密集电路：融合应显著减少指令数
    QuantumCircuit qc(2);
    qc.h(0).t(0).sx(0).h(0).cx(0, 1).rx(0.3, 1).t(1).h(1);
    const auto fused = GateFusion(2).run(qc);
    EXPECT_LT(fused.num_instructions(), qc.num_instructions());
    expect_same_state(qc, fused);
}

TEST(GateFusion, FusesControlledGatesIntoBlocks) {
    // 全部落在 2 比特内：应融合成单个 2 比特 Unitary
    QuantumCircuit qc(2);
    qc.h(0).cx(0, 1).rz(0.5, 1).cx(0, 1).h(0);
    const auto fused = GateFusion(2).run(qc);
    EXPECT_EQ(fused.num_instructions(), 1u);
    EXPECT_EQ(fused.instructions()[0].gate.type, GateType::Unitary);
    EXPECT_EQ(fused.instructions()[0].qubits.size(), 2u);
    expect_same_state(qc, fused);
}

TEST(GateFusion, RespectsMaxBlockSize) {
    // 3 比特电路，max_block=2 时不能把全部融成一块
    QuantumCircuit qc(3);
    qc.h(0).cx(0, 1).cx(1, 2).h(2);
    const auto fused = GateFusion(2).run(qc);
    for (const auto& inst : fused.instructions()) {
        EXPECT_LE(inst.qubits.size(), 2u);
    }
    expect_same_state(qc, fused);
}

TEST(GateFusion, DoesNotCrossBarrierOrMeasure) {
    QuantumCircuit qc(1, 1);
    qc.h(0).t(0).barrier().h(0).t(0).measure(0, 0);
    const auto fused = GateFusion(4).run(qc);
    // barrier 与 measure 保留；两段各自融合
    bool has_barrier = false, has_measure = false;
    for (const auto& inst : fused.instructions()) {
        if (inst.gate.type == GateType::Barrier) has_barrier = true;
        if (inst.gate.type == GateType::Measure) has_measure = true;
    }
    EXPECT_TRUE(has_barrier);
    EXPECT_TRUE(has_measure);
}

TEST(GateFusion, LeavesLargePermutationUntouched) {
    // 模幂置换门（大门）不应被融合成稠密矩阵
    QuantumCircuit qc(3);
    const std::vector<std::size_t> cycle = {1, 2, 3, 4, 5, 6, 7, 0};
    qc.h(0).permutation(cycle, {0, 1, 2}).h(0);
    const auto fused = GateFusion(4).run(qc);
    bool has_perm = false;
    for (const auto& inst : fused.instructions()) {
        if (inst.gate.type == GateType::Permutation) has_perm = true;
    }
    EXPECT_TRUE(has_perm);
    expect_same_state(qc, fused);
}

TEST(GateFusion, BlockSizeOneMatchesSingleQubitFusion) {
    const auto qc = random_circuit(4, 50, 42);
    expect_same_state(GateFusion(1).run(qc), SingleQubitGateFusion().run(qc));
}

TEST(GateFusion, SimulatorResultsIndependentOfBlockSize) {
    QuantumCircuit qc(4, 4);
    qc.h(0).cx(0, 1).t(1).rx(0.4, 2).cx(1, 2).sx(3).cx(2, 3).h(0).measure_all();
    Result ref;
    for (std::size_t mb : {1u, 2u, 3u, 4u}) {
        StatevectorSimulator sim(
            {.seed = 5, .fuse_gates = true, .fusion_max_qubits = mb});
        auto r = sim.run(qc, 1024);
        if (mb == 1)
            ref = r;
        else
            EXPECT_EQ(r.counts(), ref.counts()) << "block size " << mb;
    }
}
