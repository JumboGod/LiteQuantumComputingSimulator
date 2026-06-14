// Pauli 期望值：解析真值与厄米性
#include <cmath>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

Statevector prepare(const QuantumCircuit& qc) {
    return StatevectorSimulator().run_statevector(qc);
}

}  // namespace

TEST(Expectation, ComputationalBasisStates) {
    QuantumCircuit qc(2);
    qc.x(0);  // |01>
    auto sv = prepare(qc);
    EXPECT_NEAR(sv.expectation_value("IZ"), -1.0, kTol);  // Z on q0
    EXPECT_NEAR(sv.expectation_value("ZI"), +1.0, kTol);  // Z on q1
    EXPECT_NEAR(sv.expectation_value("ZZ"), -1.0, kTol);
    EXPECT_NEAR(sv.expectation_value("II"), +1.0, kTol);
    EXPECT_NEAR(sv.expectation_value("IX"), 0.0, kTol);
}

TEST(Expectation, SuperpositionAndPhase) {
    {  // |+>：<X> = 1，<Z> = 0
        QuantumCircuit qc(1);
        qc.h(0);
        auto sv = prepare(qc);
        EXPECT_NEAR(sv.expectation_value("X"), 1.0, kTol);
        EXPECT_NEAR(sv.expectation_value("Z"), 0.0, kTol);
        EXPECT_NEAR(sv.expectation_value("Y"), 0.0, kTol);
    }
    {  // |+i> = S·H|0>：<Y> = 1
        QuantumCircuit qc(1);
        qc.h(0).s(0);
        auto sv = prepare(qc);
        EXPECT_NEAR(sv.expectation_value("Y"), 1.0, kTol);
        EXPECT_NEAR(sv.expectation_value("X"), 0.0, kTol);
    }
}

TEST(Expectation, BellCorrelations) {
    QuantumCircuit qc(2);
    qc.h(0).cx(0, 1);
    auto sv = prepare(qc);
    // |Φ+>：<XX> = <ZZ> = 1，<YY> = -1，单比特期望全 0
    EXPECT_NEAR(sv.expectation_value("XX"), +1.0, kTol);
    EXPECT_NEAR(sv.expectation_value("ZZ"), +1.0, kTol);
    EXPECT_NEAR(sv.expectation_value("YY"), -1.0, kTol);
    EXPECT_NEAR(sv.expectation_value("ZI"), 0.0, kTol);
    EXPECT_NEAR(sv.expectation_value("IX"), 0.0, kTol);
}

TEST(Expectation, RotationAngleSweep) {
    // RY(θ)|0>：<Z> = cos θ，<X> = sin θ
    for (const double theta : {0.3, 1.1, 2.5}) {
        QuantumCircuit qc(1);
        qc.ry(theta, 0);
        auto sv = prepare(qc);
        EXPECT_NEAR(sv.expectation_value("Z"), std::cos(theta), kTol);
        EXPECT_NEAR(sv.expectation_value("X"), std::sin(theta), kTol);
    }
}

TEST(Expectation, HamiltonianAggregation) {
    using namespace lqcs::algorithms;
    QuantumCircuit qc(2);
    qc.x(0);
    const Hamiltonian h = {{0.5, "IZ"}, {2.0, "ZI"}, {-1.0, "II"}};
    // 0.5·(-1) + 2.0·(+1) - 1.0 = 0.5
    EXPECT_NEAR(expectation(prepare(qc), h), 0.5, kTol);
}

TEST(Expectation, RejectsInvalidInput) {
    Statevector sv(2);
    EXPECT_THROW(sv.expectation_value("Z"), std::invalid_argument);   // 长度
    EXPECT_THROW(sv.expectation_value("ZA"), std::invalid_argument);  // 字符
}
