#pragma once

#include <vector>

#include "gate.hpp"

namespace lqcs {

struct Instruction {
    Gate gate;
    std::vector<qubit_t> qubits;  // 受控门约定：控制位在前，目标位在后
    std::vector<clbit_t> clbits;  // 仅 Measure 使用
};

}  // namespace lqcs
