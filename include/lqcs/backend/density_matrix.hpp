#pragma once

#include <cstddef>
#include <vector>

#include "../circuit/instruction.hpp"
#include "../core/types.hpp"
#include "noise.hpp"
#include "statevector.hpp"

namespace lqcs {

// 密度矩阵 ρ（2^n × 2^n），初始为 |0...0><0...0|。
// 内部按向量化表示存储：|ρ>> 视作 2n 比特状态向量
// （低 n 位 = ket 侧，高 n 位 = bra 侧），门作用 UρU† 即
// 对 ket 侧施加 U、bra 侧施加 U*，完全复用状态向量内核。
class DensityMatrix {
public:
    static constexpr std::size_t kMaxQubits = 13;  // 4^13 × 16B = 1 GiB

    explicit DensityMatrix(std::size_t num_qubits);

    std::size_t num_qubits() const { return num_qubits_; }
    std::size_t dim() const { return std::size_t{1} << num_qubits_; }

    // ρ 的矩阵元 <row|ρ|col>
    complex_t element(std::size_t row, std::size_t col) const;

    // ρ → U ρ U†（支持任意受控门/置换门指令）
    void apply_instruction(const Instruction& inst);

    // ρ → Σ K_i ρ K_i†，通道作用于比特 q
    void apply_channel(const KrausChannel& channel, qubit_t q);

    double trace() const;                  // 物理态恒为 1
    double purity() const;                 // Tr(ρ²)：纯态 1，混态 < 1
    std::vector<double> probabilities() const;  // 对角元（计算基测量分布）
    // 与纯态的保真度 <ψ|ρ|ψ>
    double fidelity(const Statevector& psi) const;

private:
    std::size_t num_qubits_;
    aligned_vector<complex_t> data_;  // 长度 4^n 的向量化 |ρ>>
};

}  // namespace lqcs
