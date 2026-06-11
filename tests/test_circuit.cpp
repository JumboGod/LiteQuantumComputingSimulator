// 电路 IR：构建期校验、流式 API、深度计算
#include <stdexcept>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

TEST(Circuit, FluentChainBuildsInstructions) {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();
    EXPECT_EQ(qc.num_instructions(), 4u);  // h + cx + 2 次 measure
    EXPECT_EQ(qc.instructions()[0].gate.type, GateType::H);
    EXPECT_EQ(qc.instructions()[1].gate.type, GateType::X);  // CX = X + 1 控制位
    EXPECT_EQ(qc.instructions()[1].gate.n_controls, 1u);
    EXPECT_EQ(qc.instructions()[1].qubits, (std::vector<qubit_t>{0, 1}));
}

TEST(Circuit, MultiControlAndCustomGates) {
    QuantumCircuit qc(4);
    qc.ccx(0, 1, 2).cswap(0, 1, 2).mcx({0, 1, 2}, 3).mcp(0.5, {0, 1}, 3);
    EXPECT_EQ(qc.instructions()[2].gate.n_controls, 3u);
    EXPECT_EQ(qc.instructions()[2].qubits, (std::vector<qubit_t>{0, 1, 2, 3}));
    EXPECT_THROW(qc.mcx({0, 1}, 1), std::invalid_argument);  // 控制/目标重叠
    EXPECT_THROW(qc.mcx({}, 0), std::invalid_argument);      // 无控制位
}

TEST(Circuit, ComposeMapsQubits) {
    QuantumCircuit sub(2);
    sub.h(0).cx(0, 1);
    QuantumCircuit qc(4);
    const std::vector<qubit_t> map = {2, 3};
    qc.compose(sub, map);
    EXPECT_EQ(qc.instructions()[0].qubits, (std::vector<qubit_t>{2}));
    EXPECT_EQ(qc.instructions()[1].qubits, (std::vector<qubit_t>{2, 3}));
    // 非法映射
    const std::vector<qubit_t> bad = {2};
    EXPECT_THROW(qc.compose(sub, bad), std::invalid_argument);
}

TEST(Circuit, InverseReversesAndInverts) {
    QuantumCircuit qc(2);
    qc.s(0).cx(0, 1);
    auto inv = qc.inverse();
    EXPECT_EQ(inv.instructions()[0].gate.type, GateType::X);   // cx 在前
    EXPECT_EQ(inv.instructions()[1].gate.type, GateType::Sdg); // s → sdg
    // 含测量的电路不可逆
    QuantumCircuit qm(1, 1);
    qm.h(0).measure(0, 0);
    EXPECT_THROW(qm.inverse(), std::invalid_argument);
}

TEST(Circuit, PowerRepeatsCircuit) {
    QuantumCircuit qc(1);
    qc.t(0);
    EXPECT_EQ(qc.power(3).num_instructions(), 3u);
    EXPECT_EQ(qc.power(0).num_instructions(), 0u);
    EXPECT_EQ(qc.power(-2).instructions()[0].gate.type, GateType::Tdg);
}

TEST(Circuit, DrawProducesReadableOutput) {
    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();
    const std::string art = qc.draw();
    EXPECT_NE(art.find("q_0:"), std::string::npos);
    EXPECT_NE(art.find("q_1:"), std::string::npos);
    EXPECT_NE(art.find("[h]"), std::string::npos);
    EXPECT_NE(art.find("[x]"), std::string::npos);  // cx 的目标位
    EXPECT_NE(art.find('*'), std::string::npos);    // cx 的控制位
    EXPECT_NE(art.find('M'), std::string::npos);
    EXPECT_NE(art.find("c: 2/"), std::string::npos);
}

TEST(Circuit, BoundsAreCheckedAtBuildTime) {
    QuantumCircuit qc(2, 1);
    EXPECT_THROW(qc.h(2), std::out_of_range);
    EXPECT_THROW(qc.cx(0, 5), std::out_of_range);
    EXPECT_THROW(qc.measure(0, 3), std::out_of_range);
    EXPECT_THROW(qc.cx(1, 1), std::invalid_argument);
    EXPECT_THROW(QuantumCircuit(0), std::invalid_argument);
}

TEST(Circuit, MeasureAllExtendsClassicalRegister) {
    QuantumCircuit qc(3, 0);
    qc.measure_all();
    EXPECT_EQ(qc.num_clbits(), 3u);
    EXPECT_EQ(qc.num_instructions(), 3u);
}

TEST(Circuit, DepthIsCriticalPath) {
    QuantumCircuit qc(3);
    EXPECT_EQ(qc.depth(), 0u);
    qc.h(0);          // q0 层 1
    qc.h(1);          // q1 层 1（与上面并行）
    EXPECT_EQ(qc.depth(), 1u);
    qc.cx(0, 1);      // q0/q1 层 2
    qc.x(2);          // q2 层 1
    EXPECT_EQ(qc.depth(), 2u);
    qc.barrier();     // 不计入深度
    EXPECT_EQ(qc.depth(), 2u);
    qc.cx(1, 2);      // 层 3
    EXPECT_EQ(qc.depth(), 3u);
}
