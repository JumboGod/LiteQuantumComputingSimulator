// 门矩阵正确性：幺正性、解析真值、特殊参数等价
#include <cmath>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

// C = A · B（k×k 行主序）
std::vector<complex_t> matmul(const std::vector<complex_t>& a,
                              const std::vector<complex_t>& b, std::size_t k) {
    std::vector<complex_t> c(k * k, complex_t{0, 0});
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t l = 0; l < k; ++l)
            for (std::size_t j = 0; j < k; ++j)
                c[i * k + j] += a[i * k + l] * b[l * k + j];
    return c;
}

std::vector<complex_t> dagger(const std::vector<complex_t>& a, std::size_t k) {
    std::vector<complex_t> d(k * k);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) d[j * k + i] = std::conj(a[i * k + j]);
    return d;
}

void expect_identity(const std::vector<complex_t>& a, std::size_t k) {
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) {
            const complex_t expected = (i == j) ? complex_t{1, 0} : complex_t{0, 0};
            EXPECT_NEAR(std::abs(a[i * k + j] - expected), 0.0, kTol)
                << "entry (" << i << "," << j << ")";
        }
}

void expect_matrix_eq(const std::vector<complex_t>& a,
                      const std::vector<complex_t>& b) {
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_NEAR(std::abs(a[i] - b[i]), 0.0, kTol) << "entry " << i;
    }
}

}  // namespace

TEST(Gates, AllUnitaryGatesSatisfyUnitarity) {
    const double th = 0.7321;  // 任意非特殊角度
    const std::vector<Gate> gates = {
        {GateType::I, {}},   {GateType::X, {}},    {GateType::Y, {}},
        {GateType::Z, {}},   {GateType::H, {}},    {GateType::S, {}},
        {GateType::Sdg, {}}, {GateType::T, {}},    {GateType::Tdg, {}},
        {GateType::SX, {}},  {GateType::RX, {th}}, {GateType::RY, {th}},
        {GateType::RZ, {th}}, {GateType::P, {th}}, {GateType::U, {th, 0.3, 1.1}},
        {GateType::CX, {}},  {GateType::CZ, {}},   {GateType::CP, {th}},
        {GateType::SWAP, {}},
    };
    for (const auto& g : gates) {
        const std::size_t k = std::size_t{1} << g.num_qubits();
        const auto m = g.matrix();
        SCOPED_TRACE("gate: " + g.name());
        expect_identity(matmul(m, dagger(m, k), k), k);
    }
}

TEST(Gates, KnownTruthValues) {
    // H|0> = (|0> + |1>)/√2：取 H 矩阵第一列
    const auto h = Gate{GateType::H, {}}.matrix();
    const double s = std::numbers::sqrt2 / 2.0;
    EXPECT_NEAR(std::abs(h[0] - complex_t{s, 0}), 0.0, kTol);
    EXPECT_NEAR(std::abs(h[2] - complex_t{s, 0}), 0.0, kTol);

    // S·S = Z, T·T = S
    const auto z = Gate{GateType::Z, {}}.matrix();
    const auto ss = matmul(Gate{GateType::S, {}}.matrix(),
                           Gate{GateType::S, {}}.matrix(), 2);
    expect_matrix_eq(ss, z);
    const auto tt = matmul(Gate{GateType::T, {}}.matrix(),
                           Gate{GateType::T, {}}.matrix(), 2);
    expect_matrix_eq(tt, Gate{GateType::S, {}}.matrix());

    // SX·SX = X
    const auto sxsx = matmul(Gate{GateType::SX, {}}.matrix(),
                             Gate{GateType::SX, {}}.matrix(), 2);
    expect_matrix_eq(sxsx, Gate{GateType::X, {}}.matrix());
}

TEST(Gates, UGateSpecialCases) {
    const double pi = std::numbers::pi;
    // U(π/2, 0, π) = H（精确，无全局相位差）
    expect_matrix_eq(Gate{GateType::U, {pi / 2, 0, pi}}.matrix(),
                     Gate{GateType::H, {}}.matrix());
    // U(0, 0, λ) = P(λ)
    expect_matrix_eq(Gate{GateType::U, {0, 0, 0.9}}.matrix(),
                     Gate{GateType::P, {0.9}}.matrix());
}

TEST(Gates, RZIsDiagonalPhase) {
    const double th = 1.234;
    const auto m = Gate{GateType::RZ, {th}}.matrix();
    EXPECT_NEAR(std::abs(m[1]), 0.0, kTol);
    EXPECT_NEAR(std::abs(m[2]), 0.0, kTol);
    EXPECT_NEAR(std::abs(m[0] - std::exp(complex_t{0, -th / 2})), 0.0, kTol);
    EXPECT_NEAR(std::abs(m[3] - std::exp(complex_t{0, th / 2})), 0.0, kTol);
}

TEST(Gates, MeasureHasNoMatrix) {
    EXPECT_FALSE((Gate{GateType::Measure, {}}.is_unitary()));
    EXPECT_THROW((Gate{GateType::Measure, {}}.matrix()), std::logic_error);
}
