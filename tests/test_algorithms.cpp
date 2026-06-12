// 经典量子算法：Deutsch-Jozsa、Bernstein-Vazirani、Simon、QPE、Grover
#include <cmath>
#include <numbers>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

// —— Deutsch-Jozsa ——

TEST(DeutschJozsa, DetectsConstantFunctions) {
    auto r0 = deutsch_jozsa(4, [](std::uint64_t) { return false; });
    EXPECT_TRUE(r0.is_constant);
    EXPECT_NEAR(r0.p_all_zero, 1.0, 1e-12);
    auto r1 = deutsch_jozsa(4, [](std::uint64_t) { return true; });
    EXPECT_TRUE(r1.is_constant);
    EXPECT_NEAR(r1.p_all_zero, 1.0, 1e-12);
}

TEST(DeutschJozsa, DetectsBalancedFunctions) {
    // 奇偶函数与「最低位」函数都是平衡的
    auto parity = deutsch_jozsa(
        4, [](std::uint64_t x) { return std::popcount(x) % 2 == 1; });
    EXPECT_FALSE(parity.is_constant);
    EXPECT_NEAR(parity.p_all_zero, 0.0, 1e-12);
    auto lowbit = deutsch_jozsa(3, [](std::uint64_t x) { return x & 1; });
    EXPECT_FALSE(lowbit.is_constant);
}

// —— Bernstein-Vazirani ——

TEST(BernsteinVazirani, RecoversSecretInOneQuery) {
    for (std::uint64_t secret : {0b1011ULL, 0b0001ULL, 0b1111ULL, 0b0000ULL}) {
        auto res = bernstein_vazirani(4, [secret](std::uint64_t x) {
            return std::popcount(x & secret) % 2 == 1;
        });
        EXPECT_EQ(res.secret, secret);
    }
}

TEST(BernsteinVazirani, WorksOnWiderRegisters) {
    const std::uint64_t secret = 0b10110101;
    auto res = bernstein_vazirani(8, [secret](std::uint64_t x) {
        return std::popcount(x & secret) % 2 == 1;
    });
    EXPECT_EQ(res.secret, secret);
}

// —— Simon ——

namespace {

// 标准 Simon oracle：f(x) = min(x, x^s)，保证 f(x) = f(x⊕s) 且别无碰撞
std::function<std::uint64_t(std::uint64_t)> simon_oracle(std::uint64_t s) {
    return [s](std::uint64_t x) { return std::min(x, x ^ s); };
}

}  // namespace

TEST(Simon, RecoversPeriodString) {
    for (std::uint64_t secret : {0b011ULL, 0b101ULL, 0b110ULL}) {
        auto res = simon(3, simon_oracle(secret), 42);
        EXPECT_EQ(res.secret, secret) << "secret=" << secret;
    }
    auto res4 = simon(4, simon_oracle(0b1101), 7);
    EXPECT_EQ(res4.secret, 0b1101u);
    auto res5 = simon(5, simon_oracle(0b10011), 7);
    EXPECT_EQ(res5.secret, 0b10011u);
}

// —— 量子相位估计 ——

TEST(QPE, ExactPhaseForTGate) {
    // T|1> = e^{iπ/4}|1>，φ = 1/8，3 个计数比特精确表示
    QuantumCircuit prep(1);
    prep.x(0);
    auto res = phase_estimation(Gate{GateType::T, {}}, prep, 3);
    EXPECT_DOUBLE_EQ(res.phase, 0.125);
    EXPECT_EQ(res.measured, 1u);
}

TEST(QPE, ExactPhaseForSAndZ) {
    QuantumCircuit prep(1);
    prep.x(0);
    EXPECT_DOUBLE_EQ(phase_estimation(Gate{GateType::S, {}}, prep, 4).phase, 0.25);
    EXPECT_DOUBLE_EQ(phase_estimation(Gate{GateType::Z, {}}, prep, 2).phase, 0.5);
}

TEST(QPE, ApproximatesIrrationalPhase) {
    // P(2π·0.3)|1> = e^{2πi·0.3}|1>：8 比特精度内误差 < 2^-8
    const double target = 0.3;
    QuantumCircuit prep(1);
    prep.x(0);
    auto res = phase_estimation(
        Gate{GateType::P, {2 * std::numbers::pi * target}}, prep, 8);
    EXPECT_NEAR(res.phase, target, 1.0 / 256.0);
}

TEST(QPE, TwoQubitUnitary) {
    // CZ 的本征态 |11>，本征值 -1 → φ = 1/2
    QuantumCircuit prep(2);
    prep.x(0).x(1);
    auto res = phase_estimation(Gate{GateType::RZZ, {std::numbers::pi}}, prep, 3);
    // RZZ(π)|11> = e^{-iπ/2}|11> → φ = 3/4
    EXPECT_DOUBLE_EQ(res.phase, 0.75);
}

// —— Grover ——

TEST(Grover, FindsSingleMarkedState) {
    const std::vector<std::uint64_t> marked = {5};
    auto res = grover(3, marked);
    EXPECT_EQ(res.best_state, 5u);
    EXPECT_GT(res.success_probability, 0.90);
    EXPECT_EQ(res.iterations, 2u);
}

TEST(Grover, FindsMultipleMarkedStates) {
    const std::vector<std::uint64_t> marked = {3, 12};
    auto res = grover(4, marked);
    EXPECT_TRUE(res.best_state == 3 || res.best_state == 12);
    EXPECT_GT(res.success_probability, 0.90);
}

TEST(Grover, LargerSearchSpace) {
    const std::vector<std::uint64_t> marked = {123};
    auto res = grover(8, marked);  // 256 个状态中找 1 个
    EXPECT_EQ(res.best_state, 123u);
    EXPECT_GT(res.success_probability, 0.99);
}

TEST(Grover, RejectsInvalidInput) {
    const std::vector<std::uint64_t> empty;
    EXPECT_THROW(grover(3, empty), std::invalid_argument);
    const std::vector<std::uint64_t> out_of_range = {8};
    EXPECT_THROW(grover(3, out_of_range), std::invalid_argument);
}

// —— QFT 构建器（库版本）——

TEST(AlgorithmsQFT, LibraryQftMatchesAnalytic) {
    const std::size_t n = 4, x = 9;
    QuantumCircuit qc(n);
    qc.x(0).x(3);
    qc.compose(qft(n));
    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t dim = 1u << n;
    for (std::size_t k = 0; k < dim; ++k) {
        const double ph = 2.0 * std::numbers::pi * static_cast<double>(x * k) / dim;
        const complex_t expected =
            complex_t{std::cos(ph), std::sin(ph)} / std::sqrt(static_cast<double>(dim));
        EXPECT_NEAR(std::abs(sv[k] - expected), 0.0, 1e-12) << "k=" << k;
    }
}
