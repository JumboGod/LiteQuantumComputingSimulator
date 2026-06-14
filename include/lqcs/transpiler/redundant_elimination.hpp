#pragma once

#include "pass.hpp"

namespace lqcs {

// 冗余门消除：删除相邻的互逆门对（H·H、X·X、CX·CX、S·Sdg、
// RZ(θ)·RZ(-θ)…）。"相邻"指两门作用于完全相同的比特列表，且其间
// 没有任何触及这些比特的指令。迭代到不动点（处理级联，如 X H H X）。
class RedundantGateElimination final : public Pass {
   public:
    QuantumCircuit run(const QuantumCircuit& circuit) const override;
    std::string name() const override { return "redundant_gate_elimination"; }
};

}  // namespace lqcs
