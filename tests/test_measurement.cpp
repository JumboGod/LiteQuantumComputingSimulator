// 测量采样：确定性结果、统计分布、可复现性、比特串格式
#include <stdexcept>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

StatevectorSimulator seeded(std::uint64_t seed) {
    return StatevectorSimulator({.seed = seed});
}

}  // namespace

TEST(Measurement, DeterministicOutcome) {
    QuantumCircuit qc(1, 1);
    qc.x(0).measure(0, 0);
    auto result = seeded(7).run(qc, 100);
    ASSERT_EQ(result.counts().size(), 1u);
    EXPECT_EQ(result.counts().at("1"), 100u);
}

TEST(Measurement, BellDistribution) {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();
    const std::size_t shots = 4096;
    auto result = seeded(42).run(qc, shots);

    std::size_t total = 0;
    for (const auto& [key, count] : result.counts()) {
        EXPECT_TRUE(key == "00" || key == "11") << "unexpected outcome: " << key;
        total += count;
    }
    EXPECT_EQ(total, shots);
    // 双侧 5σ 余量内（σ = √(N·p·(1-p)) = 32）
    EXPECT_NEAR(static_cast<double>(result.counts().at("00")), shots / 2.0, 160.0);
    EXPECT_NEAR(static_cast<double>(result.counts().at("11")), shots / 2.0, 160.0);
}

TEST(Measurement, SameSeedReproduces) {
    QuantumCircuit qc(2, 2);
    qc.h(0).h(1).measure_all();
    auto r1 = seeded(123).run(qc, 1000);
    auto r2 = seeded(123).run(qc, 1000);
    EXPECT_EQ(r1.counts(), r2.counts());
}

TEST(Measurement, ClbitMappingAndStringFormat) {
    // 把 qubit0 的结果写到 clbit1：x(0) 后应得 "10"（最高位 clbit 在左）
    QuantumCircuit qc(2, 2);
    qc.x(0).measure(0, 1);
    auto result = seeded(5).run(qc, 10);
    ASSERT_EQ(result.counts().size(), 1u);
    EXPECT_EQ(result.counts().begin()->first, "10");
}

TEST(Measurement, ProbabilitiesSumToOne) {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();
    auto result = seeded(9).run(qc, 2048);
    double sum = 0;
    for (const auto& [key, p] : result.probabilities()) sum += p;
    EXPECT_NEAR(sum, 1.0, 1e-12);
}

TEST(Measurement, RejectsUnsupportedCircuits) {
    {   // 无测量
        QuantumCircuit qc(1, 1);
        qc.h(0);
        EXPECT_THROW(StatevectorSimulator().run(qc, 10), std::invalid_argument);
    }
    {   // shots = 0
        QuantumCircuit qc(1, 1);
        qc.measure(0, 0);
        EXPECT_THROW(StatevectorSimulator().run(qc, 0), std::invalid_argument);
    }
}

TEST(Measurement, GroverFindsMarkedState) {
    // 3 比特 Grover 搜索 |101>，2 次迭代后理论命中率 ~94.5%
    const double pi = 3.14159265358979323846;
    QuantumCircuit qc(3, 3);
    for (qubit_t q = 0; q < 3; ++q) qc.h(q);
    for (int iter = 0; iter < 2; ++iter) {
        qc.x(1);
        qc.mcp(pi, {0, 1}, 2);  // oracle: 翻转 |101> 相位
        qc.x(1);
        for (qubit_t q = 0; q < 3; ++q) qc.h(q);
        for (qubit_t q = 0; q < 3; ++q) qc.x(q);
        qc.mcp(pi, {0, 1}, 2);  // 扩散算子
        for (qubit_t q = 0; q < 3; ++q) qc.x(q);
        for (qubit_t q = 0; q < 3; ++q) qc.h(q);
    }
    qc.measure_all();

    const std::size_t shots = 4096;
    auto result = seeded(2026).run(qc, shots);
    const double p101 =
        static_cast<double>(result.counts().at("101")) / static_cast<double>(shots);
    EXPECT_GT(p101, 0.90);
    EXPECT_LT(p101, 0.99);
}
