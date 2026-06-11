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
    EXPECT_EQ(qc.instructions()[1].gate.type, GateType::CX);
    EXPECT_EQ(qc.instructions()[1].qubits, (std::vector<qubit_t>{0, 1}));
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
