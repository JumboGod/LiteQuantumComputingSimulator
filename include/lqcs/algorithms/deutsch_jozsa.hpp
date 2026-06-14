#pragma once

#include <cstdint>
#include <functional>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

struct DeutschJozsaResult {
    bool is_constant;  // f 是常数函数还是平衡函数
    double p_all_zero;  // 输入寄存器测得全 0 的概率（常数 = 1，平衡 = 0）
    QuantumCircuit circuit;
};

// Deutsch-Jozsa：单次 oracle 查询判定 f: {0,1}^n → {0,1}
// 是常数函数还是平衡函数（承诺二者之一）。
DeutschJozsaResult deutsch_jozsa(std::size_t n,
                                 const std::function<bool(std::uint64_t)>& f);

}  // namespace lqcs::algorithms
