// 编译优化 pass：优化前后状态向量逐元素相等（容差 1e-12）
#include <cmath>
#include <memory>
#include <random>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>
#include <lqcs/transpiler/gate_fusion.hpp>
#include <lqcs/transpiler/redundant_elimination.hpp>
#include <lqcs/transpiler/transpiler.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

void expect_same_state(const QuantumCircuit& a, const QuantumCircuit& b) {
    StatevectorSimulator sim({.fuse_gates = false});
    const auto sa = sim.run_statevector(a);
    const auto sb = sim.run_statevector(b);
    ASSERT_EQ(sa.size(), sb.size());
    for (std::size_t i = 0; i < sa.size(); ++i) {
        EXPECT_NEAR(std::abs(sa[i] - sb[i]), 0.0, kTol) << "amplitude " << i;
    }
}

QuantumCircuit random_1q_heavy_circuit(std::size_t n, std::size_t n_gates,
                                       std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle(0, 6.28);
    QuantumCircuit qc(n);
    for (std::size_t i = 0; i < n_gates; ++i) {
        const auto q = static_cast<qubit_t>(rng() % n);
        switch (rng() % 6) {
            case 0: qc.h(q); break;
            case 1: qc.t(q); break;
            case 2: qc.rx(angle(rng), q); break;
            case 3: qc.u(angle(rng), angle(rng), angle(rng), q); break;
            case 4: qc.sx(q); break;
            default: {
                const auto q2 = static_cast<qubit_t>((q + 1) % n);
                qc.cx(q, q2);
                break;
            }
        }
    }
    return qc;
}

}  // namespace

TEST(GateFusion, PreservesStateOnRandomCircuits) {
    for (std::uint64_t seed : {1u, 2u, 3u}) {
        const auto qc = random_1q_heavy_circuit(4, 60, seed);
        const auto fused = SingleQubitGateFusion().run(qc);
        SCOPED_TRACE("seed " + std::to_string(seed));
        EXPECT_LT(fused.num_instructions(), qc.num_instructions());
        expect_same_state(qc, fused);
    }
}

TEST(GateFusion, MergesRunsIntoSingleUnitary) {
    QuantumCircuit qc(1);
    qc.h(0).t(0).sx(0).rz(0.3, 0);
    const auto fused = SingleQubitGateFusion().run(qc);
    ASSERT_EQ(fused.num_instructions(), 1u);
    EXPECT_EQ(fused.instructions()[0].gate.type, GateType::Unitary);
    expect_same_state(qc, fused);
}

TEST(GateFusion, IdentityProductIsDropped) {
    QuantumCircuit qc(1);
    qc.h(0).h(0);  // H·H = I
    EXPECT_EQ(SingleQubitGateFusion().run(qc).num_instructions(), 0u);
}

TEST(GateFusion, DoesNotCrossMeasureOrBarrier) {
    QuantumCircuit qc(1, 1);
    qc.h(0).measure(0, 0).h(0);
    const auto fused = SingleQubitGateFusion().run(qc);
    // 两个 H 不能合并：中间隔着测量
    EXPECT_EQ(fused.num_instructions(), 3u);

    QuantumCircuit qb(1);
    qb.h(0).barrier().h(0);
    // barrier 同样阻断融合：两个 H 各自成为独立 Unitary
    EXPECT_EQ(SingleQubitGateFusion().run(qb).num_instructions(), 3u);
}

TEST(GateFusion, ControlledGatesAreNotFused) {
    QuantumCircuit qc(2);
    qc.cx(0, 1).cx(0, 1);
    EXPECT_EQ(SingleQubitGateFusion().run(qc).num_instructions(), 2u);
}

TEST(RedundantElimination, RemovesInversePairs) {
    QuantumCircuit qc(2);
    qc.h(0).h(0).x(1).x(1).cx(0, 1).cx(0, 1).s(0).sdg(0).rz(0.7, 1).rz(-0.7, 1);
    EXPECT_EQ(RedundantGateElimination().run(qc).num_instructions(), 0u);
}

TEST(RedundantElimination, CascadesToFixpoint) {
    QuantumCircuit qc(1);
    qc.x(0).h(0).h(0).x(0);  // 消掉 HH 后 XX 相邻 → 全消
    EXPECT_EQ(RedundantGateElimination().run(qc).num_instructions(), 0u);
}

TEST(RedundantElimination, BlockedByInterveningGate) {
    QuantumCircuit qc(2);
    qc.h(0).cx(0, 1).h(0);  // CX 触及 q0：两个 H 不相邻
    EXPECT_EQ(RedundantGateElimination().run(qc).num_instructions(), 3u);

    QuantumCircuit qc2(2);
    qc2.cx(0, 1).cx(1, 0);  // 方向不同：不是逆对
    EXPECT_EQ(RedundantGateElimination().run(qc2).num_instructions(), 2u);
}

TEST(RedundantElimination, PreservesStateOnRandomCircuits) {
    const auto qc = random_1q_heavy_circuit(4, 50, 9);
    expect_same_state(qc, RedundantGateElimination().run(qc));
}

TEST(Transpiler, PipelineAppliesAllPasses) {
    Transpiler tp;
    tp.add_pass(std::make_unique<RedundantGateElimination>())
        .add_pass(std::make_unique<SingleQubitGateFusion>());
    EXPECT_EQ(tp.num_passes(), 2u);

    const auto qc = random_1q_heavy_circuit(5, 80, 21);
    const auto opt = tp.run(qc);
    EXPECT_LT(opt.num_instructions(), qc.num_instructions());
    expect_same_state(qc, opt);
}

TEST(Simulator, FusionOptionPreservesResults) {
    // 末尾测量路径：同种子下融合开/关 counts 完全一致
    QuantumCircuit qc(3, 3);
    qc.h(0).t(0).cx(0, 1).rx(0.4, 2).sx(2).cx(1, 2).measure_all();
    auto on =
        StatevectorSimulator({.seed = 5, .fuse_gates = true}).run(qc, 512);
    auto off =
        StatevectorSimulator({.seed = 5, .fuse_gates = false}).run(qc, 512);
    EXPECT_EQ(on.counts(), off.counts());

    // 逐 shot 路径（中间测量）：rng 消耗序列不变，结果同样一致
    QuantumCircuit qm(2, 2);
    qm.h(0).t(0).measure(0, 0).h(1).cx(0, 1).measure(1, 1);
    auto pon =
        StatevectorSimulator({.seed = 9, .fuse_gates = true}).run(qm, 256);
    auto poff =
        StatevectorSimulator({.seed = 9, .fuse_gates = false}).run(qm, 256);
    EXPECT_EQ(pon.counts(), poff.counts());
}
