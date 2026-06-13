#pragma once

#include <cstddef>

#include "lqcs/core/types.hpp"

// 状态向量门作用内核：对长度 n_amps = 2^n 的振幅数组做原位变换。
// 全部基于位运算枚举振幅组，不构造 2^n × 2^n 矩阵。
// 局部比特序约定：基门矩阵/置换表的局部位 j 对应 targets[j]。
namespace lqcs::kernels {

// 通用单比特门：m 为行主序 2×2 矩阵
void apply_single_qubit(complex_t* state, std::size_t n_amps, qubit_t target,
                        const complex_t m[4]);

// 对角单比特门快路径：只乘相位，不配对访存
void apply_diagonal_single_qubit(complex_t* state, std::size_t n_amps,
                                 qubit_t target, complex_t d0, complex_t d1);

// 单控制单比特门：仅作用于 control 位为 1 的子空间
void apply_controlled_single_qubit(complex_t* state, std::size_t n_amps,
                                   qubit_t control, qubit_t target,
                                   const complex_t m[4]);

// SWAP：纯振幅交换，无复数乘法
void apply_swap(complex_t* state, std::size_t n_amps, qubit_t a, qubit_t b);

// 通用双比特门：m 为行主序 4×4 矩阵（局部位 0 ↔ q0，位 1 ↔ q1）
void apply_two_qubit(complex_t* state, std::size_t n_amps, qubit_t q0, qubit_t q1,
                     const complex_t m[16]);

// 通用受控幺正：在 controls 全为 1 的子空间上，对 k 个 targets
// 应用 2^k × 2^k 基门矩阵。复杂度 O(2^(n-c) · 4^k)。
void apply_controlled_unitary(complex_t* state, std::size_t n_amps,
                              const qubit_t* controls, std::size_t n_controls,
                              const qubit_t* targets, std::size_t n_targets,
                              const complex_t* m);

// 通用受控置换：|l> → |perm[l]>，零浮点乘法。Shor 受控模幂的核心。
void apply_controlled_permutation(complex_t* state, std::size_t n_amps,
                                  const qubit_t* controls, std::size_t n_controls,
                                  const qubit_t* targets, std::size_t n_targets,
                                  const std::size_t* perm);

// 多比特 Pauli 旋转 exp(-iθ/2·P)：O(2^n) 振幅对变换，不构造稠密矩阵。
// x_global/zy_global 为已映射到全局比特的 X/Y 位与 Z/Y 位掩码，n_y 为 Y 个数。
void apply_pauli_rotation(complex_t* state, std::size_t n_amps,
                          std::size_t x_global, std::size_t zy_global,
                          unsigned n_y, double theta);

// 按 Born 规则坍缩一个比特：random01 ∈ [0,1) 决定分支。
// 原位投影并归一化，返回测量结果。
bool collapse_qubit(complex_t* state, std::size_t n_amps, qubit_t q,
                    double random01);

// 重置比特到 |0>：先按 Born 规则坍缩，若得 1 再翻转回 0
void reset_qubit(complex_t* state, std::size_t n_amps, qubit_t q, double random01);

// 轨迹法噪声：单比特 Kraus 算符的分支概率 ||K|ψ>||²（只读，不修改态）。
// 分支应用复用 apply_single_qubit（把 1/√p 归一化融进矩阵即可）。
double kraus_probability(const complex_t* state, std::size_t n_amps,
                         qubit_t target, const complex_t k[4]);

}  // namespace lqcs::kernels
