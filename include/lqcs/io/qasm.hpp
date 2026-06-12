#pragma once

#include <string>
#include <string_view>

#include "../circuit/quantum_circuit.hpp"

namespace lqcs::io {

// 导出为 OpenQASM 2.0（qelib1.inc 门集）。
// 支持：全部单比特门（P→u1、U→u3）、cx/cy/cz/ch/cu1/crx/cry/crz、
// swap/ccx/cswap、measure/reset/barrier。
// iSWAP/RXX/RYY/RZZ/Unitary/Permutation/3+ 控制门不在 OpenQASM 2.0
// 标准门集内，遇到时抛出异常。
std::string to_qasm(const QuantumCircuit& circuit);

// 解析 OpenQASM 2.0 子集（to_qasm 的输出 + Qiskit 常见导出）。
// 要求恰好一个 qreg、至多一个 creg；支持 pi 四则运算的参数表达式。
QuantumCircuit from_qasm(std::string_view source);

}  // namespace lqcs::io
