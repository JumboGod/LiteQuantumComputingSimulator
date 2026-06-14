// VQE：H2 基态收敛（M6 验收标准）
#include <span>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

namespace {

const Hamiltonian kH2 = {
    {-1.052373245772859, "II"},   {+0.39793742484318045, "IZ"},
    {-0.39793742484318045, "ZI"}, {-0.01128010425623538, "ZZ"},
    {+0.18093119978423156, "XX"},
};
constexpr double kH2Exact = -1.8572750302023824;

// 教科书标准 H2 ansatz：X(0)·RY(θ,1)·CX(1,0)
// 制备 cos(θ/2)|01> + sin(θ/2)|10>，精确覆盖基态所在的奇宇称子空间
QuantumCircuit h2_ansatz(std::span<const double> theta) {
    QuantumCircuit qc(2);
    qc.x(0);
    qc.ry(theta[0], 1);
    qc.cx(1, 0);
    return qc;
}

}  // namespace

TEST(VQE, ConvergesToH2GroundState) {
    const auto res = vqe(kH2, h2_ansatz, 1, {.max_iterations = 50});
    EXPECT_NEAR(res.energy, kH2Exact, 1e-9);
    // 变分原理：能量不低于精确基态
    EXPECT_GE(res.energy, kH2Exact - 1e-9);
    EXPECT_EQ(res.parameters.size(), 1u);
    EXPECT_FALSE(res.history.empty());
}

TEST(VQE, HistoryIsMonotonicallyDecreasing) {
    const auto res = vqe(kH2, h2_ansatz, 1, {.max_iterations = 20});
    for (std::size_t i = 1; i < res.history.size(); ++i) {
        EXPECT_LE(res.history[i], res.history[i - 1] + 1e-12);
    }
}

TEST(VQE, SingleQubitAnalyticCase) {
    // H = Z：基态 |1>，能量 -1；ansatz RY(θ)|0> 在 θ=π 处达到
    const Hamiltonian h = {{1.0, "Z"}};
    const auto res = vqe(h,
                         [](std::span<const double> t) {
                             QuantumCircuit qc(1);
                             qc.ry(t[0], 0);
                             return qc;
                         },
                         1, {.max_iterations = 5});
    EXPECT_NEAR(res.energy, -1.0, 1e-10);
}

TEST(VQE, RespectsInitialParameters) {
    const auto res = vqe(kH2, h2_ansatz, 1,
                         {.max_iterations = 50, .initial_parameters = {0.5}});
    EXPECT_NEAR(res.energy, kH2Exact, 1e-9);
    EXPECT_THROW(vqe(kH2, h2_ansatz, 1, {.initial_parameters = {1.0, 2.0}}),
                 std::invalid_argument);
}
