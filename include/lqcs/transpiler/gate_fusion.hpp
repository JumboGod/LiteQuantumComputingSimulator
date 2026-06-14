#pragma once

#include <cstddef>

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

// 多比特门融合（qsim/Aer 风格）：贪心地把相邻门（含 CX/CCX 等受控门）
// 聚合成不超过 max_block_size 比特的稠密 Unitary 块，把多次全状态
// 向量扫描合并为一次。块内门数 ≥ 2 才融合（单门保留以走快速内核）。
// 不跨越 Measure/Reset/Barrier；过大的门（如模幂置换门）原样保留。
class GateFusion final : public Pass {
   public:
    explicit GateFusion(std::size_t max_block_size = 4)
        : max_block_size_(max_block_size) {}
    QuantumCircuit run(const QuantumCircuit& circuit) const override;
    std::string name() const override { return "gate_fusion"; }

   private:
    std::size_t max_block_size_;
};

}  // namespace lqcs
