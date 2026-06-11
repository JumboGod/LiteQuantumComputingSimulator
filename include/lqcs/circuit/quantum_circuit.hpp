#pragma once

#include <cstddef>
#include <initializer_list>
#include <span>
#include <string>
#include <vector>

#include "instruction.hpp"

namespace lqcs {

// 电路中间表示（IR）：纯数据结构，只记录门序列，不持有量子态。
// 比特序约定（与 Qiskit 一致，little-endian）：qubit 0 为最低有效位。
// 受控门约定：指令的 qubits 中控制位在前、目标位在后。
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

    // —— 受控门（控制位在前）——
    QuantumCircuit& cx(qubit_t control, qubit_t target);
    QuantumCircuit& cy(qubit_t control, qubit_t target);
    QuantumCircuit& cz(qubit_t control, qubit_t target);
    QuantumCircuit& ch(qubit_t control, qubit_t target);
    QuantumCircuit& cp(double lambda, qubit_t control, qubit_t target);
    QuantumCircuit& crx(double theta, qubit_t control, qubit_t target);
    QuantumCircuit& cry(double theta, qubit_t control, qubit_t target);
    QuantumCircuit& crz(double theta, qubit_t control, qubit_t target);
    QuantumCircuit& ccx(qubit_t c1, qubit_t c2, qubit_t target);  // Toffoli
    QuantumCircuit& cswap(qubit_t control, qubit_t a, qubit_t b); // Fredkin
    QuantumCircuit& mcx(std::span<const qubit_t> controls, qubit_t target);
    QuantumCircuit& mcx(std::initializer_list<qubit_t> controls, qubit_t target);
    QuantumCircuit& mcp(double lambda, std::span<const qubit_t> controls,
                        qubit_t target);
    QuantumCircuit& mcp(double lambda, std::initializer_list<qubit_t> controls,
                        qubit_t target);

    // —— 双比特门 ——
    QuantumCircuit& swap(qubit_t a, qubit_t b);
    QuantumCircuit& iswap(qubit_t a, qubit_t b);
    QuantumCircuit& rxx(double theta, qubit_t a, qubit_t b);
    QuantumCircuit& ryy(double theta, qubit_t a, qubit_t b);
    QuantumCircuit& rzz(double theta, qubit_t a, qubit_t b);

    // —— 自定义算符 ——
    // 任意幺正矩阵（行主序 2^k × 2^k，构造时校验幺正性，容差 1e-10）。
    // 局部比特序 little-endian：矩阵的局部位 j 对应 qubits[j]。
    QuantumCircuit& unitary(std::span<const complex_t> matrix,
                            std::span<const qubit_t> qubits);
    QuantumCircuit& unitary(std::span<const complex_t> matrix,
                            std::initializer_list<qubit_t> qubits);
    // 计算基置换 |l> → |table[l]>（构造时校验双射）
    QuantumCircuit& permutation(std::span<const std::size_t> table,
                                std::span<const qubit_t> qubits);
    QuantumCircuit& permutation(std::span<const std::size_t> table,
                                std::initializer_list<qubit_t> qubits);

    // —— 非幺正操作 ——
    QuantumCircuit& measure(qubit_t q, clbit_t c);
    // 依次测量所有量子比特到同号经典比特；经典寄存器不足时自动扩展
    QuantumCircuit& measure_all();
    QuantumCircuit& reset(qubit_t q);
    QuantumCircuit& barrier();

    // 直接追加一条完整指令（算法库构建受控置换门等用）
    QuantumCircuit& append(Instruction inst);

    // —— 组合与变换 ——
    // 把子电路按 qubit_map 映射后追加到本电路；map 为空时按同号映射
    QuantumCircuit& compose(const QuantumCircuit& other,
                            std::span<const qubit_t> qubit_map = {});
    QuantumCircuit inverse() const;  // 全电路求逆（含 Measure/Reset 时抛异常）
    QuantumCircuit power(int n) const;

    // —— 访问器 ——
    std::size_t num_qubits() const { return num_qubits_; }
    std::size_t num_clbits() const { return num_clbits_; }
    std::size_t num_instructions() const { return instructions_.size(); }
    std::span<const Instruction> instructions() const { return instructions_; }
    std::size_t depth() const;  // 关键路径长度（Barrier 不计入）

    std::string draw() const;   // ASCII 电路图

private:
    QuantumCircuit& add_gate(Gate gate, std::vector<qubit_t> qubits);
    void check_qubit(qubit_t q) const;
    void check_clbit(clbit_t c) const;

    std::size_t              num_qubits_;
    std::size_t              num_clbits_;
    std::vector<Instruction> instructions_;
};

}  // namespace lqcs
