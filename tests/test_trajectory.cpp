// 轨迹法噪声模拟：与密度矩阵精确解交叉验证、解析案例、容量
#include <cmath>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

StatevectorSimulator noisy_sim(NoiseModel noise, std::uint64_t seed) {
    return StatevectorSimulator({.seed = seed, .noise = std::move(noise)});
}

NoiseModel depolarizing(double p) {
    NoiseModel n;
    n.add_all_qubit_channel(KrausChannel::depolarizing(p));
    return n;
}

}  // namespace

TEST(Trajectory, MatchesDensityMatrixExactly) {
    // 金标准：同一含噪电路，轨迹法频率 vs 密度矩阵精确概率
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).t(1).cx(1, 0).measure_all();

    const double p = 0.15;
    const std::size_t shots = 40000;

    auto traj = noisy_sim(depolarizing(p), 99).run(qc, shots);
    auto dm = DensityMatrixSimulator({.seed = 1, .noise = depolarizing(p)})
                  .run_density_matrix(qc);
    const auto exact = dm.probabilities();

    // 统计容差：5σ，σ = √(p(1-p)/N) ≤ √(0.25/40000) = 0.0025
    for (std::size_t i = 0; i < exact.size(); ++i) {
        std::string key = {i & 2 ? '1' : '0', i & 1 ? '1' : '0'};
        const auto it = traj.counts().find(key);
        const double freq =
            it == traj.counts().end()
                ? 0.0
                : static_cast<double>(it->second) / static_cast<double>(shots);
        EXPECT_NEAR(freq, exact[i], 0.0125) << "state " << key;
    }
}

TEST(Trajectory, AmplitudeDampingGammaOneIsDeterministic) {
    NoiseModel n;
    n.add_all_qubit_channel(KrausChannel::amplitude_damping(1.0));
    QuantumCircuit qc(1, 1);
    qc.x(0).measure(0, 0);  // X 后必被完全弛豫回 |0>
    auto result = noisy_sim(std::move(n), 5).run(qc, 200);
    EXPECT_EQ(result.counts().at("0"), 200u);
}

TEST(Trajectory, BitFlipRateMatchesProbability) {
    const double p = 0.3;
    NoiseModel n;
    n.add_all_qubit_channel(KrausChannel::bit_flip(p));
    QuantumCircuit qc(1, 1);
    qc.i(0).measure(0, 0);  // 恒等门只为触发一次噪声
    const std::size_t shots = 20000;
    auto result = noisy_sim(std::move(n), 11).run(qc, shots);
    const double freq = static_cast<double>(result.counts().at("1")) /
                        static_cast<double>(shots);
    EXPECT_NEAR(freq, p, 0.017);  // ~5σ
}

TEST(Trajectory, PhaseDampingAnalytic) {
    // H [PD] H [PD]：P(0) = 1/2 + (1/2)·√(1-γ)（第二次 PD 不改变布居）
    const double gamma = 0.36;
    NoiseModel n;
    n.add_all_qubit_channel(KrausChannel::phase_damping(gamma));
    QuantumCircuit qc(1, 1);
    qc.h(0).h(0).measure(0, 0);
    const std::size_t shots = 40000;
    auto result = noisy_sim(std::move(n), 17).run(qc, shots);
    const double freq = static_cast<double>(result.counts().at("0")) /
                        static_cast<double>(shots);
    EXPECT_NEAR(freq, 0.5 + 0.5 * std::sqrt(1.0 - gamma), 0.0125);
}

TEST(Trajectory, SeedReproduces) {
    QuantumCircuit qc(3, 3);
    qc.h(0).cx(0, 1).cx(1, 2).measure_all();
    auto r1 = noisy_sim(depolarizing(0.1), 123).run(qc, 500);
    auto r2 = noisy_sim(depolarizing(0.1), 123).run(qc, 500);
    EXPECT_EQ(r1.counts(), r2.counts());
}

TEST(Trajectory, ScalesBeyondDensityMatrixLimit) {
    // 18 比特含噪 GHZ：密度矩阵后端（上限 13）不可能做到
    const std::size_t n = 18;
    QuantumCircuit qc(n, n);
    qc.h(0);
    for (qubit_t q = 0; q + 1 < n; ++q) qc.cx(q, q + 1);
    qc.measure_all();
    auto result = noisy_sim(depolarizing(0.01), 7).run(qc, 16);
    EXPECT_EQ(result.shots(), 16u);
    std::size_t total = 0;
    for (const auto& [key, count] : result.counts()) total += count;
    EXPECT_EQ(total, 16u);
}

TEST(Trajectory, RunStatevectorRejectsNoise) {
    QuantumCircuit qc(1);
    qc.h(0);
    auto sim = noisy_sim(depolarizing(0.1), 1);
    EXPECT_THROW(sim.run_statevector(qc), std::invalid_argument);
}

TEST(Trajectory, ZeroNoiseChannelIsIdentity) {
    // p=0 的通道不应改变结果：与无噪声同种子完全一致
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();
    auto noisy = noisy_sim(depolarizing(0.0), 42).run(qc, 1000);
    for (const auto& [key, count] : noisy.counts()) {
        EXPECT_TRUE(key == "00" || key == "11") << key;
    }
}
