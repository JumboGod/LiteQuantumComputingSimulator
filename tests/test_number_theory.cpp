// 经典数论工具：gcd、模幂、素性、完全幂、连分数
#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs::algorithms::number_theory;

TEST(NumberTheory, GcdAndPowMod) {
    EXPECT_EQ(gcd(48, 18), 6u);
    EXPECT_EQ(gcd(17, 5), 1u);
    EXPECT_EQ(gcd(0, 7), 7u);
    EXPECT_EQ(pow_mod(7, 4, 15), 1u);    // 7 的阶是 4 (mod 15)
    EXPECT_EQ(pow_mod(2, 10, 1000), 24u);
    EXPECT_EQ(pow_mod(5, 0, 13), 1u);
    // 大数不溢出（mul_mod 用 128 位中间结果或双倍-累加）
    // 0xFFFFFFFF^2 = 0xFFFFFFFE00000001；对 0xFFFFFFFFFFFF 取模 = 0xFFFE00010000
    EXPECT_EQ(mul_mod(0xFFFFFFFFULL, 0xFFFFFFFFULL, 0xFFFFFFFFFFFFULL),
              0xFFFE00010000ULL);
}

TEST(NumberTheory, Primality) {
    EXPECT_TRUE(is_prime(2));
    EXPECT_TRUE(is_prime(17));
    EXPECT_TRUE(is_prime(97));
    EXPECT_FALSE(is_prime(1));
    EXPECT_FALSE(is_prime(15));
    EXPECT_FALSE(is_prime(91));  // 7 × 13
}

TEST(NumberTheory, PerfectPower) {
    EXPECT_EQ(perfect_power_base(8), 2u);    // 2^3
    EXPECT_EQ(perfect_power_base(27), 3u);   // 3^3
    EXPECT_EQ(perfect_power_base(49), 7u);   // 7^2
    EXPECT_EQ(perfect_power_base(15), 0u);
    EXPECT_EQ(perfect_power_base(21), 0u);
}

TEST(NumberTheory, ContinuedFractionRecoversPeriod) {
    // Shor 典型场景：N=15, a=7, r=4。测得 m=192（s=3, 2^8=256）：
    // 192/256 = 3/4 → 分母里应有 4
    const auto dens = continued_fraction_denominators(192, 256, 15);
    EXPECT_NE(std::find(dens.begin(), dens.end(), 4u), dens.end());
    // m=64：64/256 = 1/4 → 同样恢复 4
    const auto dens2 = continued_fraction_denominators(64, 256, 15);
    EXPECT_NE(std::find(dens2.begin(), dens2.end(), 4u), dens2.end());
    // 分母超过上限的收敛分数被截断
    for (const auto d : continued_fraction_denominators(123, 4096, 20)) {
        EXPECT_LE(d, 20u);
    }
}
