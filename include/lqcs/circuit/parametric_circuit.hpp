#pragma once

#include <cstddef>
#include <initializer_list>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "quantum_circuit.hpp"

namespace lqcs {

// 参数化电路（Qulacs 风格）：固定门与参数化旋转门混合，参数槽
// 原位绑定（bind 不重建电路结构）。每个参数化门引入一个自由参数，
// 按添加顺序编号；bind(θ) 把第 k 个参数门的旋转角设为 θ[k]。
//
// 参数化门限定为生成元本征值 ±1 的旋转（RX/RY/RZ/PauliRotation），
// 因此 parameter-shift 梯度精确成立。
class ParametricCircuit {
   public:
    explicit ParametricCircuit(std::size_t num_qubits,
                               std::size_t num_clbits = 0);

    // —— 固定（非参数）门 ——
    ParametricCircuit& x(qubit_t q);
    ParametricCircuit& y(qubit_t q);
    ParametricCircuit& z(qubit_t q);
    ParametricCircuit& h(qubit_t q);
    ParametricCircuit& s(qubit_t q);
    ParametricCircuit& sdg(qubit_t q);
    ParametricCircuit& t(qubit_t q);
    ParametricCircuit& sx(qubit_t q);
    ParametricCircuit& cx(qubit_t control, qubit_t target);
    ParametricCircuit& cy(qubit_t control, qubit_t target);
    ParametricCircuit& cz(qubit_t control, qubit_t target);
    ParametricCircuit& swap(qubit_t a, qubit_t b);
    ParametricCircuit& ccx(qubit_t c1, qubit_t c2, qubit_t target);
    ParametricCircuit& barrier();

    // —— 参数化旋转门：返回新引入的参数索引 ——
    std::size_t rx(qubit_t q);
    std::size_t ry(qubit_t q);
    std::size_t rz(qubit_t q);
    // 多比特 Pauli 旋转 exp(-iθ/2·P)，pauli 最右字符对应 qubits[0]
    std::size_t rp(std::string_view pauli, std::span<const qubit_t> qubits);
    std::size_t rp(std::string_view pauli,
                   std::initializer_list<qubit_t> qubits);

    // 绑定参数 → 具体电路（结构不变，仅填入旋转角）
    QuantumCircuit bind(std::span<const double> values) const;
    QuantumCircuit bind(std::initializer_list<double> values) const {
        return bind(std::span<const double>(values.begin(), values.size()));
    }

    std::size_t num_parameters() const { return param_inst_.size(); }
    std::size_t num_qubits() const { return template_.num_qubits(); }
    std::size_t num_clbits() const { return template_.num_clbits(); }

   private:
    std::size_t add_param_gate(Gate gate, std::vector<qubit_t> qubits);

    QuantumCircuit template_;              // 含占位角（0）的指令模板
    std::vector<std::size_t> param_inst_;  // 参数 id → 模板指令下标
};

}  // namespace lqcs
