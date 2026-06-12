// 态演化：与解析解逐振幅对比（容差 1e-12），范数守恒
#include <cmath>
#include <numbers>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;
const double kS = std::numbers::sqrt2 / 2.0;

void expect_amp(const Statevector& sv, std::size_t i, complex_t expected) {
    EXPECT_NEAR(std::abs(sv[i] - expected), 0.0, kTol) << "amplitude " << i;
}

}  // namespace

TEST(Statevector, InitialStateIsAllZeros) {
    Statevector sv(3);
    EXPECT_EQ(sv.size(), 8u);
    expect_amp(sv, 0, {1, 0});
    for (std::size_t i = 1; i < 8; ++i) expect_amp(sv, i, {0, 0});
    EXPECT_NEAR(sv.norm(), 1.0, kTol);
}

TEST(Statevector, RejectsOversizedRegister) {
    EXPECT_THROW(Statevector(31), std::invalid_argument);
    EXPECT_THROW(Statevector(0), std::invalid_argument);
}

TEST(Statevector, HadamardCreatesSuperposition) {
    QuantumCircuit qc(1);
    qc.h(0);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 0, {kS, 0});
    expect_amp(sv, 1, {kS, 0});
}

TEST(Statevector, XFlips) {
    QuantumCircuit qc(2);
    qc.x(0);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 1, {1, 0});  // |01> = idx 1（qubit 0 为最低位）
}

TEST(Statevector, BellState) {
    QuantumCircuit qc(2);
    qc.h(0).cx(0, 1);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 0, {kS, 0});
    expect_amp(sv, 1, {0, 0});
    expect_amp(sv, 2, {0, 0});
    expect_amp(sv, 3, {kS, 0});
}

TEST(Statevector, GHZState) {
    QuantumCircuit qc(3);
    qc.h(0).cx(0, 1).cx(1, 2);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 0, {kS, 0});
    expect_amp(sv, 7, {kS, 0});
    for (std::size_t i = 1; i < 7; ++i) expect_amp(sv, i, {0, 0});
}

TEST(Statevector, CXControlDirectionMatters) {
    {   // control=0：|01> → |11>
        QuantumCircuit qc(2);
        qc.x(0).cx(0, 1);
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 3, {1, 0});
    }
    {   // control=1 但 qubit1 为 0：CX 不动作，仍是 |01>
        QuantumCircuit qc(2);
        qc.x(0).cx(1, 0);
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 1, {1, 0});
    }
    {   // control=1 为 1：|10> → |11>
        QuantumCircuit qc(2);
        qc.x(1).cx(1, 0);
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 3, {1, 0});
    }
}

TEST(Statevector, SwapMovesAmplitude) {
    QuantumCircuit qc(2);
    qc.x(0).swap(0, 1);  // |01> → |10>
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 2, {1, 0});
}

TEST(Statevector, CZAndCPApplyPhase) {
    {
        QuantumCircuit qc(2);
        qc.x(0).x(1).cz(0, 1);  // |11> → -|11>
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 3, {-1, 0});
    }
    {
        const double th = 0.81;
        QuantumCircuit qc(2);
        qc.x(0).x(1).cp(th, 0, 1);  // |11> → e^{iθ}|11>
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 3, std::exp(complex_t{0, th}));
    }
}

TEST(Statevector, HHIsIdentity) {
    QuantumCircuit qc(1);
    qc.h(0).h(0);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 0, {1, 0});
}

TEST(Statevector, NormIsPreservedByRandomishCircuit) {
    QuantumCircuit qc(4);
    qc.h(0).rx(0.3, 1).ry(1.1, 2).rz(2.2, 3)
      .cx(0, 2).cz(1, 3).cp(0.5, 2, 3).swap(0, 1)
      .t(0).sdg(2).sx(3).u(0.4, 1.5, 2.6, 1);
    auto sv = StatevectorSimulator().run_statevector(qc);
    EXPECT_NEAR(sv.norm(), 1.0, 1e-12);
}

TEST(Statevector, ToffoliTruthTable) {
    {   // 两个控制位都为 1 才翻转：|011> → |111>
        QuantumCircuit qc(3);
        qc.x(0).x(1).ccx(0, 1, 2);
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 7, {1, 0});
    }
    {   // 只有一个控制位为 1：不动作
        QuantumCircuit qc(3);
        qc.x(0).ccx(0, 1, 2);
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 1, {1, 0});
    }
}

TEST(Statevector, FredkinSwapsWhenControlSet) {
    QuantumCircuit qc(3);
    qc.x(0).x(1).cswap(0, 1, 2);  // 控制位 q0=1：交换 q1,q2，|011> → |101>
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 5, {1, 0});
}

TEST(Statevector, PermutationGateMovesBasisStates) {
    // 2 比特循环移位 |l> → |l+1 mod 4>
    const std::vector<std::size_t> cycle = {1, 2, 3, 0};
    {
        QuantumCircuit qc(2);
        qc.x(0).permutation(cycle, {0, 1});  // |01>=idx1 → idx2
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 2, {1, 0});
    }
    {   // 受控置换：控制位为 0 时不动作（Shor 受控模幂的雏形）
        QuantumCircuit qc(3);
        qc.x(0);  // 工作位 q0=1，控制位 q2=0
        qc.append({Gate{GateType::Permutation, {}, 1, {}, cycle}, {2, 0, 1}, {}});
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 1, {1, 0});  // 不变
    }
    {   // 控制位为 1 时动作
        QuantumCircuit qc(3);
        qc.x(0).x(2);
        qc.append({Gate{GateType::Permutation, {}, 1, {}, cycle}, {2, 0, 1}, {}});
        auto sv = StatevectorSimulator().run_statevector(qc);
        expect_amp(sv, 6, {1, 0});  // 工作寄存器 |01>→|10>，加上 q2=1 → idx 6
    }
}

TEST(Statevector, CustomUnitaryMatchesBuiltinGate) {
    // 自定义 H 矩阵的行为与内置 H 完全一致
    const auto h_mat = Gate{GateType::H, {}}.matrix();
    QuantumCircuit qc1(2), qc2(2);
    qc1.unitary(h_mat, {1});
    qc2.h(1);
    auto sv1 = StatevectorSimulator().run_statevector(qc1);
    auto sv2 = StatevectorSimulator().run_statevector(qc2);
    for (std::size_t i = 0; i < sv1.size(); ++i) {
        EXPECT_NEAR(std::abs(sv1[i] - sv2[i]), 0.0, kTol);
    }
}

TEST(Statevector, TargetOnHighQubit) {
    // 在高位比特上作用门，验证 stride 展开正确
    QuantumCircuit qc(5);
    qc.x(4);
    auto sv = StatevectorSimulator().run_statevector(qc);
    expect_amp(sv, 16, {1, 0});
}
