#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "instruction.hpp"

namespace lqcs {

// 电路中间表示（IR）：纯数据结构，只记录门序列，不持有量子态。
// 比特序约定（与 Qiskit 一致，little-endian）：qubit 0 为最低有效位。
class QuantumCircuit {
public:
    explicit QuantumCircuit(std::size_t num_qubits, std::size_t num_clbits = 0);

    // —— 单比特无参门 ——
    QuantumCircuit& i(qubit_t q);
    QuantumCircuit& x(qubit_t q);
    QuantumCircuit& y(qubit_t q);
    QuantumCircuit& z(qubit_t q);
    QuantumCircuit& h(qubit_t q);
    QuantumCircuit& s(qubit_t q);
    QuantumCircuit& sdg(qubit_t q);
    QuantumCircuit& t(qubit_t q);
    QuantumCircuit& tdg(qubit_t q);
    QuantumCircuit& sx(qubit_t q);

    // —— 单比特参数门 ——
    QuantumCircuit& rx(double theta, qubit_t q);
    QuantumCircuit& ry(double theta, qubit_t q);
    QuantumCircuit& rz(double theta, qubit_t q);
    QuantumCircuit& p(double lambda, qubit_t q);
    QuantumCircuit& u(double theta, double phi, double lambda, qubit_t q);

    // —— 双比特门 ——
    QuantumCircuit& cx(qubit_t control, qubit_t target);
    QuantumCircuit& cz(qubit_t control, qubit_t target);
    QuantumCircuit& cp(double lambda, qubit_t control, qubit_t target);
    QuantumCircuit& swap(qubit_t a, qubit_t b);

    // —— 非幺正操作 ——
    QuantumCircuit& measure(qubit_t q, clbit_t c);
    // 依次测量所有量子比特到同号经典比特；经典寄存器不足时自动扩展
    QuantumCircuit& measure_all();
    QuantumCircuit& barrier();

    // —— 访问器 ——
    std::size_t num_qubits() const { return num_qubits_; }
    std::size_t num_clbits() const { return num_clbits_; }
    std::size_t num_instructions() const { return instructions_.size(); }
    std::span<const Instruction> instructions() const { return instructions_; }
    std::size_t depth() const;  // 关键路径长度（Barrier 不计入）

private:
    QuantumCircuit& add_1q(GateType type, qubit_t q, std::vector<double> params = {});
    QuantumCircuit& add_2q(GateType type, qubit_t a, qubit_t b,
                           std::vector<double> params = {});
    void check_qubit(qubit_t q) const;
    void check_clbit(clbit_t c) const;

    std::size_t              num_qubits_;
    std::size_t              num_clbits_;
    std::vector<Instruction> instructions_;
};

}  // namespace lqcs
