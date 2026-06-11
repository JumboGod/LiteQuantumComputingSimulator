#include "kernels.hpp"

#include <algorithm>
#include <utility>

#include "lqcs/core/bit_utils.hpp"

namespace lqcs::kernels {

void apply_single_qubit(complex_t* state, std::size_t n_amps, qubit_t target,
                        const complex_t m[4]) {
    for (std::size_t g = 0; g < n_amps / 2; ++g) {
        const std::size_t i0 = bits::insert_zero_bit(g, target);
        const std::size_t i1 = i0 | (std::size_t{1} << target);
        const complex_t a0 = state[i0];
        const complex_t a1 = state[i1];
        state[i0] = m[0] * a0 + m[1] * a1;
        state[i1] = m[2] * a0 + m[3] * a1;
    }
}

void apply_controlled_single_qubit(complex_t* state, std::size_t n_amps,
                                   qubit_t control, qubit_t target,
                                   const complex_t m[4]) {
    const unsigned lo = std::min(control, target);
    const unsigned hi = std::max(control, target);
    const std::size_t c_mask = std::size_t{1} << control;
    const std::size_t t_mask = std::size_t{1} << target;
    // 枚举其余 n-2 个自由比特，展开成 control=1、target=0 的下标
    for (std::size_t g = 0; g < n_amps / 4; ++g) {
        const std::size_t base =
            bits::insert_zero_bit(bits::insert_zero_bit(g, lo), hi);
        const std::size_t i0 = base | c_mask;
        const std::size_t i1 = i0 | t_mask;
        const complex_t a0 = state[i0];
        const complex_t a1 = state[i1];
        state[i0] = m[0] * a0 + m[1] * a1;
        state[i1] = m[2] * a0 + m[3] * a1;
    }
}

void apply_swap(complex_t* state, std::size_t n_amps, qubit_t a, qubit_t b) {
    const unsigned lo = std::min(a, b);
    const unsigned hi = std::max(a, b);
    const std::size_t a_mask = std::size_t{1} << a;
    const std::size_t b_mask = std::size_t{1} << b;
    // 仅 |..a=1,b=0..> 与 |..a=0,b=1..> 交换
    for (std::size_t g = 0; g < n_amps / 4; ++g) {
        const std::size_t base =
            bits::insert_zero_bit(bits::insert_zero_bit(g, lo), hi);
        std::swap(state[base | a_mask], state[base | b_mask]);
    }
}

}  // namespace lqcs::kernels
