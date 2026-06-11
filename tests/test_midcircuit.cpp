// 中间测量与 Reset：逐 shot 演化路径
#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

StatevectorSimulator seeded(std::uint64_t seed) {
    return StatevectorSimulator({.seed = seed});
}

}  // namespace

TEST(MidCircuit, MeasureCollapsesState) {
    // 测量后 CX：两个比特必然相同 —— 坍缩正确性的直接验证
    QuantumCircuit qc(2, 2);
    qc.h(0).measure(0, 0).cx(0, 1).measure(1, 1);
    auto result = seeded(42).run(qc, 2048);
    for (const auto& [key, count] : result.counts()) {
        EXPECT_TRUE(key == "00" || key == "11") << "unexpected outcome: " << key;
    }
}

TEST(MidCircuit, RemeasureAfterHadamardIsUniform) {
    // 第一次测量坍缩后再做 H，第二次测量重新 50/50：四种组合各约 1/4
    QuantumCircuit qc(1, 2);
    qc.h(0).measure(0, 0).h(0).measure(0, 1);
    const std::size_t shots = 8192;
    auto result = seeded(7).run(qc, shots);
    ASSERT_EQ(result.counts().size(), 4u);
    for (const auto& [key, count] : result.counts()) {
        EXPECT_NEAR(static_cast<double>(count), shots / 4.0, 280.0)  // ~6σ
            << "outcome " << key;
    }
}

TEST(MidCircuit, ResetForcesZero) {
    {   // |1> 经 reset 后测得 0
        QuantumCircuit qc(1, 1);
        qc.x(0).reset(0).measure(0, 0);
        auto result = seeded(3).run(qc, 100);
        EXPECT_EQ(result.counts().at("0"), 100u);
    }
    {   // 叠加态经 reset 后确定为 0，再做 X 确定为 1
        QuantumCircuit qc(1, 1);
        qc.h(0).reset(0).x(0).measure(0, 0);
        auto result = seeded(4).run(qc, 100);
        EXPECT_EQ(result.counts().at("1"), 100u);
    }
}

TEST(MidCircuit, ResetOnEntangledQubitLeavesPartnerDistributed) {
    // Bell 态中 reset 一个比特：另一个仍 50/50
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).reset(0).measure_all();
    const std::size_t shots = 4096;
    auto result = seeded(11).run(qc, shots);
    std::size_t q1_ones = 0;
    for (const auto& [key, count] : result.counts()) {
        EXPECT_EQ(key[1], '0') << "q0 should always be 0 after reset";
        if (key[0] == '1') q1_ones += count;  // key[0] 是 clbit1（最高位在左）
    }
    EXPECT_NEAR(static_cast<double>(q1_ones), shots / 2.0, 200.0);
}

TEST(MidCircuit, SameSeedReproduces) {
    QuantumCircuit qc(2, 2);
    qc.h(0).measure(0, 0).h(1).cx(1, 0).reset(1).measure_all();
    auto r1 = seeded(123).run(qc, 500);
    auto r2 = seeded(123).run(qc, 500);
    EXPECT_EQ(r1.counts(), r2.counts());
}

TEST(MidCircuit, RunStatevectorRejectsReset) {
    QuantumCircuit qc(1);
    qc.h(0).reset(0);
    EXPECT_THROW(StatevectorSimulator().run_statevector(qc), std::invalid_argument);
}
