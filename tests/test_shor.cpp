// Shor 算法端到端：固定种子下分解 15/21/35 必须成功
#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;
namespace nt = lqcs::algorithms::number_theory;

namespace {

void expect_factors(std::uint64_t N, std::uint64_t f1, std::uint64_t f2,
                    std::uint64_t seed) {
    const auto res = shor(N, {.seed = seed});
    ASSERT_TRUE(res.success) << "N=" << N << " seed=" << seed;
    EXPECT_EQ(res.factor1, f1);
    EXPECT_EQ(res.factor2, f2);
    EXPECT_EQ(res.factor1 * res.factor2, N);
}

}  // namespace

TEST(Shor, Factors15) { expect_factors(15, 3, 5, 42); }
TEST(Shor, Factors21) { expect_factors(21, 3, 7, 42); }
TEST(Shor, Factors35) { expect_factors(35, 5, 7, 7); }

TEST(Shor, PeriodIsValidWhenQuantum) {
    const auto res = shor(15, {.seed = 1234});
    ASSERT_TRUE(res.success);
    if (res.period != 0) {  // 量子路径（而非幸运 gcd 捷径）
        EXPECT_EQ(res.period % 2, 0u);
        EXPECT_EQ(nt::pow_mod(res.a, res.period, 15), 1u);
        EXPECT_FALSE(res.counts.empty());
        ASSERT_TRUE(res.circuit.has_value());
        EXPECT_EQ(res.circuit->num_qubits(), 12u);  // 2*4 计数 + 4 工作
    }
}

TEST(Shor, QuantumPathFactors21) {
    // 扫种子直到底数 a 与 21 互素（绕开幸运 gcd 捷径），强制走量子求周期
    bool quantum_success = false;
    for (std::uint64_t seed = 1; seed <= 20; ++seed) {
        const auto res = shor(21, {.seed = seed});
        ASSERT_TRUE(res.success) << "seed=" << seed;
        EXPECT_EQ(res.factor1 * res.factor2, 21u);
        if (res.period != 0) {
            EXPECT_EQ(res.period % 2, 0u);
            EXPECT_EQ(nt::pow_mod(res.a, res.period, 21), 1u);
            quantum_success = true;
            break;
        }
    }
    EXPECT_TRUE(quantum_success) << "no seed exercised the quantum path";
}

TEST(Shor, ClassicalShortcuts) {
    // 偶数
    auto even = shor(12);
    EXPECT_TRUE(even.success);
    EXPECT_EQ(even.factor1, 2u);
    // 素数幂 27 = 3^3
    auto pp = shor(27, {.seed = 5});
    EXPECT_TRUE(pp.success);
    EXPECT_EQ(pp.factor1, 3u);
    EXPECT_EQ(pp.factor2, 9u);
    // 素数无非平凡因子
    EXPECT_FALSE(shor(17, {.seed = 5}).success);
    // 过小的 N
    EXPECT_THROW(shor(3), std::invalid_argument);
}

TEST(Shor, ReproducibleWithSeed) {
    const auto r1 = shor(15, {.seed = 99});
    const auto r2 = shor(15, {.seed = 99});
    ASSERT_TRUE(r1.success);
    EXPECT_EQ(r1.a, r2.a);
    EXPECT_EQ(r1.period, r2.period);
    EXPECT_EQ(r1.attempts, r2.attempts);
}
