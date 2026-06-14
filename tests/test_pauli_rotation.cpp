// 多比特 Pauli 旋转门、参数化电路、parameter-shift 梯度
#include <cmath>
#include <numbers>
#include <random>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

namespace {

constexpr double kTol = 1e-12;

void expect_same_state(const QuantumCircuit& a, const QuantumCircuit& b) {
    StatevectorSimulator sim({.fuse_gates = false});
    const auto sa = sim.run_statevector(a);
    const auto sb = sim.run_statevector(b);
    ASSERT_EQ(sa.size(), sb.size());
    for (std::size_t i = 0; i < sa.size(); ++i) {
        EXPECT_NEAR(std::abs(sa[i] - sb[i]), 0.0, kTol) << "amplitude " << i;
    }
}

}  // namespace

TEST(PauliRotation, MatchesNamedTwoQubitRotations) {
    // rp("XX") ≡ RXX, rp("YY") ≡ RYY, rp("ZZ") ≡ RZZ
    for (const double th : {0.3, 1.7, -2.1}) {
        {
            QuantumCircuit a(2), b(2);
            a.h(0).rp(th, "XX", {0, 1});
            b.h(0).rxx(th, 0, 1);
            expect_same_state(a, b);
        }
        {
            QuantumCircuit a(2), b(2);
            a.h(0).h(1).rp(th, "YY", {0, 1});
            b.h(0).h(1).ryy(th, 0, 1);
            expect_same_state(a, b);
        }
        {
            QuantumCircuit a(2), b(2);
            a.h(0).rp(th, "ZZ", {0, 1});
            b.h(0).rzz(th, 0, 1);
            expect_same_state(a, b);
        }
    }
}

TEST(PauliRotation, MatchesSingleQubitRotations) {
    for (const double th : {0.5, -1.3}) {
        {
            QuantumCircuit a(1), b(1);
            a.rp(th, "X", {0});
            b.rx(th, 0);
            expect_same_state(a, b);
        }
        {
            QuantumCircuit a(1), b(1);
            a.rp(th, "Y", {0});
            b.ry(th, 0);
            expect_same_state(a, b);
        }
        {
            QuantumCircuit a(1), b(1);
            a.rp(th, "Z", {0});
            b.rz(th, 0);
            expect_same_state(a, b);
        }
    }
}

TEST(PauliRotation, KernelMatchesDenseMatrixForLongStrings) {
    // 专用内核 vs base_matrix（稠密参考）在多比特 Pauli 串上一致
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> angle(-3.0, 3.0);
    for (const std::string pauli : {"XYZ", "ZXY", "IYX", "YYYY", "ZIZI"}) {
        const std::size_t k = pauli.size();
        const double th = angle(rng);
        std::vector<qubit_t> qs(k);
        for (std::size_t i = 0; i < k; ++i) qs[i] = static_cast<qubit_t>(i);

        // 内核路径
        QuantumCircuit qc(k);
        for (std::size_t i = 0; i < k; ++i) qc.h(static_cast<qubit_t>(i));
        qc.rp(th, pauli, qs);
        const auto sv_kernel =
            StatevectorSimulator({.fuse_gates = false}).run_statevector(qc);

        // 稠密参考路径（同样初态 + Unitary 门）
        QuantumCircuit ref(k);
        for (std::size_t i = 0; i < k; ++i) ref.h(static_cast<qubit_t>(i));
        const Gate rp_gate{GateType::PauliRotation, {th}, 0, {}, {}, pauli};
        ref.unitary(rp_gate.base_matrix(), qs);
        const auto sv_ref =
            StatevectorSimulator({.fuse_gates = false}).run_statevector(ref);

        SCOPED_TRACE("pauli " + pauli);
        for (std::size_t i = 0; i < sv_kernel.size(); ++i) {
            EXPECT_NEAR(std::abs(sv_kernel[i] - sv_ref[i]), 0.0, 1e-11)
                << "amp " << i;
        }
    }
}

TEST(PauliRotation, GateIsUnitaryAndInvertible) {
    const Gate g{GateType::PauliRotation, {0.9}, 0, {}, {}, "XYZ"};
    const std::size_t dim = 8;
    const auto m = g.matrix();
    const auto inv = g.inverse().matrix();
    // m · inv = I
    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) {
            complex_t sum{0, 0};
            for (std::size_t l = 0; l < dim; ++l)
                sum += m[i * dim + l] * inv[l * dim + j];
            const complex_t expected =
                (i == j) ? complex_t{1, 0} : complex_t{0, 0};
            EXPECT_NEAR(std::abs(sum - expected), 0.0, kTol);
        }
    }
}

TEST(PauliRotation, RejectsInvalidPauli) {
    QuantumCircuit qc(2);
    EXPECT_THROW(qc.rp(0.5, "XAB", {0, 1}), std::invalid_argument);  // 字符
    EXPECT_THROW(qc.rp(0.5, "X", {0, 1}), std::invalid_argument);  // 长度不匹配
    EXPECT_THROW(qc.rp(0.5, "", {}), std::invalid_argument);       // 空
}

// —— ParametricCircuit ——

TEST(ParametricCircuit, BindMatchesDirectConstruction) {
    // 固定门返回引用可链式；参数门返回参数索引，需分开调用
    ParametricCircuit pc(2);
    EXPECT_EQ(pc.ry(0), 0u);  // 参数 0
    EXPECT_EQ(pc.rx(1), 1u);  // 参数 1
    pc.cx(0, 1);
    EXPECT_EQ(pc.rz(0), 2u);  // 参数 2
    EXPECT_EQ(pc.num_parameters(), 3u);

    const std::vector<double> theta = {0.4, 1.2, -0.7};
    const QuantumCircuit bound = pc.bind(theta);

    QuantumCircuit expected(2);
    expected.ry(0.4, 0).rx(1.2, 1).cx(0, 1).rz(-0.7, 0);
    expect_same_state(bound, expected);
}

TEST(ParametricCircuit, RebindChangesAngles) {
    ParametricCircuit pc(1);
    pc.ry(0);
    auto sv0 = StatevectorSimulator().run_statevector(pc.bind({0.0}));
    EXPECT_NEAR(std::abs(sv0[0] - complex_t{1, 0}), 0.0,
                kTol);  // RY(0)|0> = |0>
    auto sv1 =
        StatevectorSimulator().run_statevector(pc.bind({std::numbers::pi}));
    EXPECT_NEAR(std::abs(sv1[1] - complex_t{1, 0}), 0.0,
                kTol);  // RY(π)|0> = |1>
}

TEST(ParametricCircuit, ParametricPauliRotation) {
    ParametricCircuit pc(2);
    pc.h(0).rp("ZZ", {0, 1});
    QuantumCircuit expected(2);
    expected.h(0).rzz(0.6, 0, 1);
    expect_same_state(pc.bind({0.6}), expected);
}

TEST(ParametricCircuit, RejectsWrongParamCount) {
    ParametricCircuit pc(1);
    pc.rx(0);
    pc.rz(0);
    EXPECT_THROW(pc.bind({0.1}), std::invalid_argument);
    EXPECT_THROW(pc.bind({0.1, 0.2, 0.3}), std::invalid_argument);
}

// —— parameter-shift 梯度 ——

TEST(Gradient, MatchesFiniteDifference) {
    const Hamiltonian h = {{1.0, "ZZ"}, {0.5, "XI"}, {-0.3, "IZ"}};
    ParametricCircuit pc(2);
    pc.ry(0);
    pc.ry(1);
    pc.cx(0, 1);
    pc.rz(0);
    pc.rp("XX", {0, 1});

    const std::vector<double> theta = {0.3, -0.8, 1.1, 0.5};
    const auto grad = parameter_shift_gradient(h, pc, theta);
    ASSERT_EQ(grad.size(), 4u);

    StatevectorSimulator sim;
    const auto energy = [&](std::vector<double> p) {
        return expectation(sim.run_statevector(pc.bind(p)), h);
    };
    const double eps = 1e-6;
    for (std::size_t k = 0; k < 4; ++k) {
        auto plus = theta, minus = theta;
        plus[k] += eps;
        minus[k] -= eps;
        const double fd = (energy(plus) - energy(minus)) / (2 * eps);
        EXPECT_NEAR(grad[k], fd, 1e-5) << "param " << k;
    }
}

TEST(Gradient, VanishesAtOptimum) {
    // H = Z，ansatz RY(θ)|0>，最优 θ=π（基态 |1>）梯度应为 0
    const Hamiltonian h = {{1.0, "Z"}};
    ParametricCircuit pc(1);
    pc.ry(0);
    const std::vector<double> at_opt = {std::numbers::pi};
    const auto grad = parameter_shift_gradient(h, pc, at_opt);
    EXPECT_NEAR(grad[0], 0.0, 1e-9);
}

TEST(GradientVQE, ConvergesToH2GroundState) {
    const Hamiltonian h2 = {
        {-1.052373245772859, "II"},   {+0.39793742484318045, "IZ"},
        {-0.39793742484318045, "ZI"}, {-0.01128010425623538, "ZZ"},
        {+0.18093119978423156, "XX"},
    };
    // X(0)·RY(θ,1)·CX(1,0) 单参数 ansatz（覆盖基态所在子空间）
    ParametricCircuit pc(2);
    pc.x(0);
    pc.ry(1);
    pc.cx(1, 0);

    auto res = vqe_gradient_descent(
        h2, pc, {.max_iterations = 500, .learning_rate = 0.3, .tol = 1e-7});
    EXPECT_NEAR(res.energy, -1.8572750302023824, 1e-5);
}
