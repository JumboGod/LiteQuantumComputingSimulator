#pragma once

#include "pass.hpp"

namespace lqcs {

// 单比特门融合：同一比特上连续的单比特门（无控制位）合并为一个
// 2×2 Unitary。状态向量模拟的瓶颈是内存带宽——每融合一个门就省
// 一次全状态向量扫描。乘积为恒等时直接消去（HH、XX 等顺带优化）。
// 不跨越 Measure/Reset/Barrier 及任何触及该比特的多比特门。
class SingleQubitGateFusion final : public Pass {
public:
    QuantumCircuit run(const QuantumCircuit& circuit) const override;
    std::string name() const override { return "single_qubit_gate_fusion"; }
};

}  // namespace lqcs
