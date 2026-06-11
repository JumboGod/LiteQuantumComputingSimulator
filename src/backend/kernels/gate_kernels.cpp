#include "kernels.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "lqcs/core/bit_utils.hpp"

namespace lqcs::kernels {

namespace {

// 把固定比特位置（升序）之外的紧凑编号 g 展开成这些位置为 0 的完整下标
std::size_t expand(std::size_t g, const std::vector<unsigned>& sorted_positions) {
    for (unsigned pos : sorted_positions) g = bits::insert_zero_bit(g, pos);
    return g;
}

}  // namespace

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

void apply_diagonal_single_qubit(complex_t* state, std::size_t n_amps,
                                 qubit_t target, complex_t d0, complex_t d1) {
    for (std::size_t i = 0; i < n_amps; ++i) {
        state[i] *= bits::test_bit(i, target) ? d1 : d0;
    }
}

void apply_controlled_single_qubit(complex_t* state, std::size_t n_amps,
                                   qubit_t control, qubit_t target,
                                   const complex_t m[4]) {
    const unsigned lo = std::min(control, target);
    const unsigned hi = std::max(control, target);
    const std::size_t c_mask = std::size_t{1} << control;
    const std::size_t t_mask = std::size_t{1} << target;
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

void apply_two_qubit(complex_t* state, std::size_t n_amps, qubit_t q0, qubit_t q1,
                     const complex_t m[16]) {
    const unsigned lo = std::min(q0, q1);
    const unsigned hi = std::max(q0, q1);
    const std::size_t m0 = std::size_t{1} << q0;
    const std::size_t m1 = std::size_t{1} << q1;
    for (std::size_t g = 0; g < n_amps / 4; ++g) {
        const std::size_t base =
            bits::insert_zero_bit(bits::insert_zero_bit(g, lo), hi);
        const std::size_t idx[4] = {base, base | m0, base | m1, base | m0 | m1};
        const complex_t a[4] = {state[idx[0]], state[idx[1]], state[idx[2]],
                                state[idx[3]]};
        for (std::size_t r = 0; r < 4; ++r) {
            state[idx[r]] =
                m[r * 4] * a[0] + m[r * 4 + 1] * a[1] + m[r * 4 + 2] * a[2] +
                m[r * 4 + 3] * a[3];
        }
    }
}

namespace {

// 受控多比特门的公共枚举骨架：对 controls 全 1 的每组振幅调用 body(idx 表)
template <typename Body>
void for_each_controlled_group(std::size_t n_amps, const qubit_t* controls,
                               std::size_t n_controls, const qubit_t* targets,
                               std::size_t n_targets, Body&& body) {
    std::vector<unsigned> fixed(controls, controls + n_controls);
    fixed.insert(fixed.end(), targets, targets + n_targets);
    std::sort(fixed.begin(), fixed.end());

    std::size_t c_mask = 0;
    for (std::size_t c = 0; c < n_controls; ++c) c_mask |= std::size_t{1} << controls[c];

    const std::size_t dim = std::size_t{1} << n_targets;
    // 局部目标编号 l → 完整下标偏移
    std::vector<std::size_t> offset(dim, 0);
    for (std::size_t l = 0; l < dim; ++l) {
        for (std::size_t j = 0; j < n_targets; ++j) {
            if (bits::test_bit(l, static_cast<unsigned>(j))) {
                offset[l] |= std::size_t{1} << targets[j];
            }
        }
    }

    std::vector<std::size_t> idx(dim);
    const std::size_t n_groups = n_amps >> fixed.size();
    for (std::size_t g = 0; g < n_groups; ++g) {
        const std::size_t base = expand(g, fixed) | c_mask;
        for (std::size_t l = 0; l < dim; ++l) idx[l] = base | offset[l];
        body(idx);
    }
}

}  // namespace

void apply_controlled_unitary(complex_t* state, std::size_t n_amps,
                              const qubit_t* controls, std::size_t n_controls,
                              const qubit_t* targets, std::size_t n_targets,
                              const complex_t* m) {
    const std::size_t dim = std::size_t{1} << n_targets;
    std::vector<complex_t> a(dim);
    for_each_controlled_group(
        n_amps, controls, n_controls, targets, n_targets,
        [&](const std::vector<std::size_t>& idx) {
            for (std::size_t l = 0; l < dim; ++l) a[l] = state[idx[l]];
            for (std::size_t r = 0; r < dim; ++r) {
                complex_t sum{0, 0};
                for (std::size_t c = 0; c < dim; ++c) sum += m[r * dim + c] * a[c];
                state[idx[r]] = sum;
            }
        });
}

void apply_controlled_permutation(complex_t* state, std::size_t n_amps,
                                  const qubit_t* controls, std::size_t n_controls,
                                  const qubit_t* targets, std::size_t n_targets,
                                  const std::size_t* perm) {
    const std::size_t dim = std::size_t{1} << n_targets;
    std::vector<complex_t> a(dim);
    for_each_controlled_group(
        n_amps, controls, n_controls, targets, n_targets,
        [&](const std::vector<std::size_t>& idx) {
            for (std::size_t l = 0; l < dim; ++l) a[l] = state[idx[l]];
            for (std::size_t l = 0; l < dim; ++l) state[idx[perm[l]]] = a[l];
        });
}

bool collapse_qubit(complex_t* state, std::size_t n_amps, qubit_t q,
                    double random01) {
    double p1 = 0.0;
    for (std::size_t i = 0; i < n_amps; ++i) {
        if (bits::test_bit(i, q)) p1 += std::norm(state[i]);
    }
    const bool outcome = random01 < p1;
    const double keep_p = outcome ? p1 : 1.0 - p1;
    const double scale = 1.0 / std::sqrt(keep_p);
    for (std::size_t i = 0; i < n_amps; ++i) {
        if (bits::test_bit(i, q) == outcome) {
            state[i] *= scale;
        } else {
            state[i] = complex_t{0, 0};
        }
    }
    return outcome;
}

void reset_qubit(complex_t* state, std::size_t n_amps, qubit_t q, double random01) {
    if (!collapse_qubit(state, n_amps, q, random01)) return;
    // 坍缩到 |1> 分支：把振幅搬回 |0> 分支
    const std::size_t mask = std::size_t{1} << q;
    for (std::size_t g = 0; g < n_amps / 2; ++g) {
        const std::size_t i0 = bits::insert_zero_bit(g, q);
        state[i0] = state[i0 | mask];
        state[i0 | mask] = complex_t{0, 0};
    }
}

}  // namespace lqcs::kernels
