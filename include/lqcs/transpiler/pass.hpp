#pragma once

#include <string>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs {

// 编译优化 pass：输入电路 → 等价的优化电路
class Pass {
   public:
    virtual ~Pass() = default;
    virtual QuantumCircuit run(const QuantumCircuit& circuit) const = 0;
    virtual std::string name() const = 0;
};

}  // namespace lqcs
