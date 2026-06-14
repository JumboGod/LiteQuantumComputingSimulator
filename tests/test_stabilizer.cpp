// Stabilizer 模拟器：与状态向量交叉验证、已知 stabilizer 群、大规模容量
#include <random>
#include <string>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

StabilizerSimulator seeded(std::uint64_t seed) {
    return StabilizerSimulator({.seed = seed});
}

// 随机 Clifford 电路（仅 Clifford 门，无测量）
QuantumCircuit random_clifford(std::size_t n, std::size_t n_gates,
                               std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    QuantumCircuit qc(n);
    for (std::size_t i = 0; i < n_gates; ++i) {
        const auto q = static_cast<qubit_t>(rng() % n);
        switch (rng() % 8) {
            case 0: qc.h(q); break;
            case 1: qc.s(q); break;
            case 2: qc.sdg(q); break;
            case 3: qc.x(q); break;
            case 4: qc.z(q); break;
            case 5: qc.sx(q); break;
            default: {
                const auto b =
                    static_cast<qubit_t>((q + 1 + rng() % (n - 1)) % n);
                if (b != q) {
                    if (rng() & 1)
                        qc.cx(q, b);
                    else
                        qc.cz(q, b);
                }
                break;
            }
        }
    }
    return qc;
}

}  // namespace

TEST(Stabilizer, MatchesStatevectorDistribution) {
    // 金标准：随机 Clifford 电路的测量分布 vs 状态向量精确概率
    for (std::uint64_t seed : {1u, 2u, 3u}) {
        const std::size_t n = 5;
        QuantumCircuit qc = random_clifford(n, 40, seed);

        // 状态向量精确概率
        const auto sv = StatevectorSimulator().run_statevector(qc);
        const auto probs = sv.probabilities();

        // stabilizer 采样
        QuantumCircuit measured = qc;
        measured = QuantumCircuit(n, n);
        measured.compose(qc);
        measured.measure_all();
        const std::size_t shots = 20000;
        auto counts = seeded(seed).run(measured, shots);

        SCOPED_TRACE("seed " + std::to_string(seed));
        for (std::size_t idx = 0; idx < (std::size_t{1} << n); ++idx) {
            std::string key(n, '0');
            for (std::size_t b = 0; b < n; ++b) {
                if ((idx >> b) & 1) key[n - 1 - b] = '1';
            }
            const auto it = counts.counts().find(key);
            const double freq = it == counts.counts().end()
                                    ? 0.0
                                    : static_cast<double>(it->second) / shots;
            EXPECT_NEAR(freq, probs[idx], 0.02) << "state " << key;
        }
    }
}

TEST(Stabilizer, BellStateStabilizers) {
    QuantumCircuit qc(2);
    qc.h(0).cx(0, 1);
    auto tab = StabilizerSimulator().run_tableau(qc);
    const auto gens = tab.stabilizers();
    // |Φ+> 由 +XX 与 +ZZ 稳定
    ASSERT_EQ(gens.size(), 2u);
    EXPECT_EQ(gens[0], "+XX");
    EXPECT_EQ(gens[1], "+ZZ");
}

TEST(Stabilizer, GHZMeasurementCorrelated) {
    const std::size_t n = 6;
    QuantumCircuit qc(n, n);
    qc.h(0);
    for (qubit_t q = 0; q + 1 < n; ++q) qc.cx(q, q + 1);
    qc.measure_all();
    auto counts = seeded(42).run(qc, 4096);
    // GHZ 只出现全 0 / 全 1
    for (const auto& [key, count] : counts.counts()) {
        EXPECT_TRUE(key == std::string(n, '0') || key == std::string(n, '1'))
            << key;
    }
}

TEST(Stabilizer, DeterministicMeasurement) {
    // |0> 测量必得 0；X|0>=|1> 必得 1
    QuantumCircuit q0(1, 1);
    q0.measure(0, 0);
    EXPECT_EQ(seeded(1).run(q0, 100).counts().at("0"), 100u);

    QuantumCircuit q1(1, 1);
    q1.x(0).measure(0, 0);
    EXPECT_EQ(seeded(1).run(q1, 100).counts().at("1"), 100u);
}

TEST(Stabilizer, PlusStateIsUniform) {
    QuantumCircuit qc(1, 1);
    qc.h(0).measure(0, 0);
    const std::size_t shots = 8000;
    auto counts = seeded(7).run(qc, shots);
    EXPECT_NEAR(static_cast<double>(counts.counts().at("0")), shots / 2.0, 250);
    EXPECT_NEAR(static_cast<double>(counts.counts().at("1")), shots / 2.0, 250);
}

TEST(Stabilizer, ResetForcesZero) {
    QuantumCircuit qc(1, 1);
    qc.h(0).reset(0).measure(0, 0);
    EXPECT_EQ(seeded(3).run(qc, 200).counts().at("0"), 200u);
}

TEST(Stabilizer, ScalesToThousandQubits) {
    // 1000 比特 GHZ：状态向量法需要 2^1000 振幅，绝无可能；
    // stabilizer 法 O(n²) 内存秒级完成
    const std::size_t n = 1000;
    QuantumCircuit qc(n, n);
    qc.h(0);
    for (qubit_t q = 0; q + 1 < n; ++q) qc.cx(q, q + 1);
    qc.measure_all();
    auto counts = seeded(11).run(qc, 4);
    EXPECT_EQ(counts.shots(), 4u);
    for (const auto& [key, count] : counts.counts()) {
        EXPECT_TRUE(key == std::string(n, '0') || key == std::string(n, '1'));
    }
}

TEST(Stabilizer, RejectsNonClifford) {
    QuantumCircuit qc(1, 1);
    qc.t(0).measure(0, 0);
    EXPECT_THROW(StabilizerSimulator().run(qc, 10), std::invalid_argument);

    QuantumCircuit qc2(3, 3);
    qc2.ccx(0, 1, 2).measure_all();
    EXPECT_THROW(StabilizerSimulator().run(qc2, 10), std::invalid_argument);
}

TEST(Stabilizer, ReproducibleWithSeed) {
    QuantumCircuit qc = random_clifford(6, 50, 99);
    QuantumCircuit m(6, 6);
    m.compose(qc);
    m.measure_all();
    EXPECT_EQ(seeded(123).run(m, 500).counts(),
              seeded(123).run(m, 500).counts());
}

TEST(Stabilizer, GateSetMatchesStatevector) {
    // 每种 Clifford 门逐一与状态向量末态分布对比
    QuantumCircuit qc(3);
    qc.h(0).s(1).sdg(2).sx(0).x(1).y(2).z(0).cx(0, 1).cy(1, 2).cz(0, 2).swap(0,
                                                                             2);
    const auto probs =
        StatevectorSimulator().run_statevector(qc).probabilities();

    QuantumCircuit m(3, 3);
    m.compose(qc);
    m.measure_all();
    const std::size_t shots = 40000;
    auto counts = seeded(5).run(m, shots);
    for (std::size_t idx = 0; idx < 8; ++idx) {
        std::string key(3, '0');
        for (std::size_t b = 0; b < 3; ++b)
            if ((idx >> b) & 1) key[2 - b] = '1';
        const auto it = counts.counts().find(key);
        const double freq =
            it == counts.counts().end() ? 0.0 : double(it->second) / shots;
        EXPECT_NEAR(freq, probs[idx], 0.02) << "state " << key;
    }
}
