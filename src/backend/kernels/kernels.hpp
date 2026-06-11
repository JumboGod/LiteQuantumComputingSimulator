#pragma once

#include <cstddef>

#include "lqcs/core/types.hpp"

// 状态向量门作用内核：对长度 n_amps = 2^n 的振幅数组做原位变换。
// 全部基于位运算枚举振幅组，不构造 2^n × 2^n 矩阵。
namespace lqcs::kernels {

// 通用单比特门：m 为行主序 2×2 矩阵
void apply_single_qubit(complex_t* state, std::size_t n_amps, qubit_t target,
                        const complex_t m[4]);

// 单控制单比特门：仅作用于 control 位为 1 的子空间，O(2^(n-1))
void apply_controlled_single_qubit(complex_t* state, std::size_t n_amps,
                                   qubit_t control, qubit_t target,
                                   const complex_t m[4]);

// SWAP：交换两比特，纯振幅交换，无复数乘法
void apply_swap(complex_t* state, std::size_t n_amps, qubit_t a, qubit_t b);

}  // namespace lqcs::kernels
