// 门矩阵正确性：幺正性、解析真值、逆门、特殊参数等价
#include <cmath>
#include <numbers>
#include <vector>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

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
        for (std::size_t j = 0; j < k; ++j)
            d[j * k + i] = std::conj(a[i * k + j]);
    return d;
}

void expect_identity(const std::vector<complex_t>& a, std::size_t k) {
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) {
            const complex_t expected =
                (i == j) ? complex_t{1, 0} : complex_t{0, 0};
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

// 测试用门集合：覆盖所有类型与控制位组合
std::vector<Gate> all_test_gates() {
    const double th = 0.7321;
    const std::vector<complex_t> rot_mat = {
        // RY(0.6) 的矩阵：任意自定义幺正
        std::cos(0.3), -std::sin(0.3), std::sin(0.3), std::cos(0.3)};
    const std::vector<std::size_t> cycle4 = {1, 2, 3, 0};
    return {
        {GateType::I, {}},
        {GateType::X, {}},
        {GateType::Y, {}},
        {GateType::Z, {}},
        {GateType::H, {}},
        {GateType::S, {}},
        {GateType::Sdg, {}},
        {GateType::T, {}},
        {GateType::Tdg, {}},
        {GateType::SX, {}},
        {GateType::RX, {th}},
        {GateType::RY, {th}},
        {GateType::RZ, {th}},
        {GateType::P, {th}},
        {GateType::U, {th, 0.3, 1.1}},
        {GateType::SWAP, {}},
        {GateType::iSWAP, {}},
        {GateType::RXX, {th}},
        {GateType::RYY, {th}},
        {GateType::RZZ, {th}},
        {GateType::X, {}, 1},     // CX
        {GateType::Z, {}, 1},     // CZ
        {GateType::Y, {}, 1},     // CY
        {GateType::H, {}, 1},     // CH
        {GateType::P, {th}, 1},   // CP
        {GateType::RX, {th}, 1},  // CRX
        {GateType::X, {}, 2},     // CCX
        {GateType::SWAP, {}, 1},  // CSWAP
        {GateType::X, {}, 3},     // MCX
        {GateType::P, {th}, 3},   // MCP
        {GateType::Unitary, {}, 0, rot_mat},
        {GateType::Unitary, {}, 2, rot_mat},  // 受控自定义幺正
        {GateType::Permutation, {}, 0, {}, cycle4},
        {GateType::Permutation, {}, 1, {}, cycle4},  // 受控置换
    };
}

}  // namespace

TEST(Gates, AllUnitaryGatesSatisfyUnitarity) {
    for (const auto& g : all_test_gates()) {
        const std::size_t k = std::size_t{1} << g.num_qubits();
        const auto m = g.matrix();
        SCOPED_TRACE("gate: " + g.name());
        ASSERT_EQ(m.size(), k * k);
        expect_identity(matmul(m, dagger(m, k), k), k);
    }
}

TEST(Gates, InverseGivesIdentity) {
    for (const auto& g : all_test_gates()) {
        const std::size_t k = std::size_t{1} << g.num_qubits();
        SCOPED_TRACE("gate: " + g.name());
        expect_identity(matmul(g.matrix(), g.inverse().matrix(), k), k);
    }
}

TEST(Gates, KnownTruthValues) {
    const double s = std::numbers::sqrt2 / 2.0;
    // H|0> = (|0> + |1>)/√2：取 H 矩阵第一列
    const auto h = Gate{GateType::H, {}}.matrix();
    EXPECT_NEAR(std::abs(h[0] - complex_t{s, 0}), 0.0, kTol);
    EXPECT_NEAR(std::abs(h[2] - complex_t{s, 0}), 0.0, kTol);

    // S·S = Z, T·T = S, SX·SX = X
    expect_matrix_eq(matmul(Gate{GateType::S, {}}.matrix(),
                            Gate{GateType::S, {}}.matrix(), 2),
                     Gate{GateType::Z, {}}.matrix());
    expect_matrix_eq(matmul(Gate{GateType::T, {}}.matrix(),
                            Gate{GateType::T, {}}.matrix(), 2),
                     Gate{GateType::S, {}}.matrix());
    expect_matrix_eq(matmul(Gate{GateType::SX, {}}.matrix(),
                            Gate{GateType::SX, {}}.matrix(), 2),
                     Gate{GateType::X, {}}.matrix());

    // CCX·CCX = I（8×8）
    const Gate ccx{GateType::X, {}, 2};
    expect_identity(matmul(ccx.matrix(), ccx.matrix(), 8), 8);
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

TEST(Gates, ControlledMatrixActsOnlyOnControlOneSubspace) {
    // CX 完整矩阵（控制位为局部低位）：控制位为 0 的列保持恒等
    const auto m = Gate{GateType::X, {}, 1}.matrix();
    // 列 0 (c=0,t=0)、列 2 (c=0,t=1) 为恒等
    EXPECT_NEAR(std::abs(m[0 * 4 + 0] - complex_t{1, 0}), 0.0, kTol);
    EXPECT_NEAR(std::abs(m[2 * 4 + 2] - complex_t{1, 0}), 0.0, kTol);
    // 列 1 (c=1,t=0) → 行 3 (c=1,t=1)
    EXPECT_NEAR(std::abs(m[3 * 4 + 1] - complex_t{1, 0}), 0.0, kTol);
    EXPECT_NEAR(std::abs(m[1 * 4 + 3] - complex_t{1, 0}), 0.0, kTol);
}

TEST(Gates, CompositeNames) {
    EXPECT_EQ((Gate{GateType::X, {}, 0}.name()), "x");
    EXPECT_EQ((Gate{GateType::X, {}, 1}.name()), "cx");
    EXPECT_EQ((Gate{GateType::X, {}, 2}.name()), "ccx");
    EXPECT_EQ((Gate{GateType::X, {}, 3}.name()), "mcx");
    EXPECT_EQ((Gate{GateType::SWAP, {}, 1}.name()), "cswap");
    EXPECT_EQ((Gate{GateType::Permutation, {}, 1, {}, {0, 1}}.name()), "cperm");
}

TEST(Gates, NonUnitaryOpsHaveNoMatrixOrInverse) {
    for (auto type : {GateType::Measure, GateType::Reset, GateType::Barrier}) {
        EXPECT_FALSE((Gate{type, {}}.is_unitary()));
        EXPECT_THROW((Gate{type, {}}.matrix()), std::logic_error);
        EXPECT_THROW((Gate{type, {}}.inverse()), std::logic_error);
    }
}

TEST(Gates, CircuitValidatesCustomOperators) {
    QuantumCircuit qc(2);
    // 非幺正矩阵被拒绝
    const std::vector<complex_t> bad = {1, 1, 0, 1};
    EXPECT_THROW(qc.unitary(bad, {0}), std::invalid_argument);
    // 尺寸不匹配被拒绝
    const std::vector<complex_t> id2 = {1, 0, 0, 1};
    EXPECT_THROW(qc.unitary(id2, {0, 1}), std::invalid_argument);
    // 非双射置换被拒绝
    const std::vector<std::size_t> not_bijection = {0, 0, 1, 2};
    EXPECT_THROW(qc.permutation(not_bijection, {0, 1}), std::invalid_argument);
    // 合法的通过
    const std::vector<std::size_t> cycle = {1, 2, 3, 0};
    EXPECT_NO_THROW(qc.permutation(cycle, {0, 1}));
    EXPECT_NO_THROW(qc.unitary(id2, {0}));
}
