// 密度矩阵后端：纯态一致性、噪声通道理论值、含噪 Bell 保真度（M6 验收）
#include <cmath>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

}  // namespace

TEST(DensityMatrix, PureEvolutionMatchesStatevector) {
    QuantumCircuit qc(3);
    qc.h(0).cx(0, 1).t(1).cswap(0, 1, 2).rx(0.7, 2).ccx(0, 2, 1);

    const auto sv = StatevectorSimulator({.fuse_gates = false}).run_statevector(qc);
    const auto rho = DensityMatrixSimulator().run_density_matrix(qc);

    // ρ = |ψ><ψ|：逐矩阵元对比
    for (std::size_t r = 0; r < sv.size(); ++r) {
        for (std::size_t c = 0; c < sv.size(); ++c) {
            EXPECT_NEAR(std::abs(rho.element(r, c) - sv[r] * std::conj(sv[c])),
                        0.0, kTol)
                << "(" << r << "," << c << ")";
        }
    }
    EXPECT_NEAR(rho.trace(), 1.0, kTol);
    EXPECT_NEAR(rho.purity(), 1.0, kTol);
    EXPECT_NEAR(rho.fidelity(sv), 1.0, kTol);
}

TEST(DensityMatrix, AmplitudeDampingDecaysExcitedState) {
    // |1> 经振幅阻尼 γ：P(1) = 1-γ，且为对角混态
    const double gamma = 0.3;
    QuantumCircuit qc(1);
    qc.x(0);
    auto rho = DensityMatrixSimulator().run_density_matrix(qc);
    rho.apply_channel(KrausChannel::amplitude_damping(gamma), 0);

    const auto probs = rho.probabilities();
    EXPECT_NEAR(probs[1], 1.0 - gamma, kTol);
    EXPECT_NEAR(probs[0], gamma, kTol);
    EXPECT_NEAR(rho.trace(), 1.0, kTol);
}

TEST(DensityMatrix, PhaseDampingKillsCoherence) {
    // |+> 经相位阻尼 γ：对角不变，相干项 × √(1-γ)
    const double gamma = 0.5;
    QuantumCircuit qc(1);
    qc.h(0);
    auto rho = DensityMatrixSimulator().run_density_matrix(qc);
    rho.apply_channel(KrausChannel::phase_damping(gamma), 0);

    EXPECT_NEAR(rho.element(0, 0).real(), 0.5, kTol);
    EXPECT_NEAR(rho.element(1, 1).real(), 0.5, kTol);
    EXPECT_NEAR(rho.element(0, 1).real(), 0.5 * std::sqrt(1 - gamma), kTol);
}

TEST(DensityMatrix, DepolarizedBellFidelityMatchesTheory) {
    // M6 验收：完美 Bell 态两比特各过一次去极化通道，
    // F = (1-p)² + p(2-p)/4（边缘分布最大混合，逐次应用可解析求出）
    for (const double p : {0.05, 0.1, 0.3}) {
        QuantumCircuit bell(2);
        bell.h(0).cx(0, 1);
        const auto ideal = StatevectorSimulator().run_statevector(bell);

        auto rho = DensityMatrixSimulator().run_density_matrix(bell);
        rho.apply_channel(KrausChannel::depolarizing(p), 0);
        rho.apply_channel(KrausChannel::depolarizing(p), 1);

        const double expected = (1 - p) * (1 - p) + p * (2 - p) / 4.0;
        EXPECT_NEAR(rho.fidelity(ideal), expected, 1e-12) << "p=" << p;
        EXPECT_NEAR(rho.trace(), 1.0, kTol);
        EXPECT_LT(rho.purity(), 1.0);
    }
}

TEST(DensityMatrix, BitFlipChannel) {
    // |0> 经比特翻转 p：P(1) = p
    const double p = 0.2;
    DensityMatrix rho(1);
    rho.apply_channel(KrausChannel::bit_flip(p), 0);
    EXPECT_NEAR(rho.probabilities()[1], p, kTol);
}

TEST(DensityMatrixSimulator, NoisyCountsLeakIntoErrorStates) {
    QuantumCircuit bell(2, 2);
    bell.h(0).cx(0, 1).measure_all();

    NoiseModel noise;
    noise.add_all_qubit_channel(KrausChannel::depolarizing(0.1));
    auto noisy =
        DensityMatrixSimulator({.seed = 7, .noise = noise}).run(bell, 4096);

    // 理想情况下 01/10 不出现；噪声使其以可观概率出现
    const auto& counts = noisy.counts();
    const std::size_t errors = (counts.count("01") ? counts.at("01") : 0) +
                               (counts.count("10") ? counts.at("10") : 0);
    EXPECT_GT(errors, 100u);
    EXPECT_LT(errors, 2000u);

    // 无噪声时与状态向量后端一致地只出现 00/11
    auto clean = DensityMatrixSimulator({.seed = 7}).run(bell, 1024);
    for (const auto& [key, count] : clean.counts()) {
        EXPECT_TRUE(key == "00" || key == "11");
    }
}

TEST(DensityMatrix, ValidationAndLimits) {
    EXPECT_THROW(DensityMatrix(14), std::invalid_argument);
    // 非完备 Kraus 集合被拒绝
    EXPECT_THROW(KrausChannel({{1, 0, 0, 1}, {0, 1, 1, 0}}),
                 std::invalid_argument);
    EXPECT_THROW(KrausChannel::depolarizing(1.5), std::invalid_argument);
}
