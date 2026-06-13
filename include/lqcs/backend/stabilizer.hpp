#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../circuit/instruction.hpp"
#include "../core/types.hpp"

namespace lqcs {

// 稳定子表（Aaronson-Gottesman / CHP 算法）：用 2n×(2n+1) 二元表表示
// n 比特 stabilizer 态。Clifford 门 O(n) 更新，测量 O(n²)，内存 O(n²)，
// 可模拟数千比特的 Clifford 电路（状态向量法的指数内存在此完全规避）。
//
// 行 0..n-1 为 destabilizer 生成元，n..2n-1 为 stabilizer 生成元，
// 第 2n 行为 scratch。每行 (x_0..x_{n-1}, z_0..z_{n-1}, r) 表示一个
// Pauli 算符及其符号。
class StabilizerTableau {
public:
    explicit StabilizerTableau(std::size_t num_qubits);  // 初始 |0...0>

    std::size_t num_qubits() const { return n_; }

    // —— Clifford 门 ——
    void h(qubit_t q);
    void s(qubit_t q);
    void sdg(qubit_t q);
    void sx(qubit_t q);
    void x(qubit_t q);
    void y(qubit_t q);
    void z(qubit_t q);
    void cx(qubit_t control, qubit_t target);
    void cy(qubit_t control, qubit_t target);
    void cz(qubit_t control, qubit_t target);
    void swap(qubit_t a, qubit_t b);

    // Z 基测量：random01 ∈ [0,1) 决定随机分支，返回测量比特并坍缩
    bool measure(qubit_t q, double random01);
    // 重置到 |0>：测量后若为 1 则翻转
    void reset(qubit_t q, double random01);

    // 应用电路指令（非 Clifford 抛异常）
    void apply(const Instruction& inst);
    static bool is_clifford(const Gate& gate);

    // n 个 stabilizer 生成元的 Pauli 串（带符号，如 "+XX"），用于检验
    std::vector<std::string> stabilizers() const;

private:
    std::uint8_t& xat(std::size_t i, std::size_t j) { return xs_[i * n_ + j]; }
    std::uint8_t& zat(std::size_t i, std::size_t j) { return zs_[i * n_ + j]; }
    // row h ← row i · row h（Pauli 乘积，含符号）
    void rowsum(std::size_t h, std::size_t i);

    std::size_t               n_;
    std::vector<std::uint8_t> xs_;  // (2n+1) × n
    std::vector<std::uint8_t> zs_;  // (2n+1) × n
    std::vector<std::uint8_t> rs_;  // (2n+1) 相位位
};

}  // namespace lqcs
