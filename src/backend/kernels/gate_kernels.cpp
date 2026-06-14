#include "kernels.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <utility>
#include <vector>

#include "lqcs/core/bit_utils.hpp"
#include "simd.hpp"

namespace lqcs::kernels {

namespace {

// 小于该规模不并行：线程启动开销超过收益
constexpr std::size_t kParallelThreshold = std::size_t{1} << 14;

// 把固定比特位置（升序）之外的紧凑编号 g 展开成这些位置为 0 的完整下标
std::size_t expand(std::size_t g, const std::vector<unsigned>& sorted_positions) {
    for (unsigned pos : sorted_positions) g = bits::insert_zero_bit(g, pos);
    return g;
}

}  // namespace

// 单比特门按「块」结构遍历：目标位 k 把态切成 n_amps/(2·stride) 个块，
// 每块前 stride 个振幅为 |0> 半、后 stride 个为 |1> 半，块内两半各自连续，
// 利于 SIMD 与缓存。stride = 2^k。
void apply_single_qubit(complex_t* state, std::size_t n_amps, qubit_t target,
                        const complex_t m[4]) {
    const std::size_t stride = std::size_t{1} << target;
    const std::size_t block = 2 * stride;
    const std::size_t n_blocks = n_amps / block;
    const complex_t m0 = m[0], m1 = m[1], m2 = m[2], m3 = m[3];
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t b = 0; b < n_blocks; ++b) {
        complex_t* p0 = state + b * block;
        complex_t* p1 = p0 + stride;
        std::size_t off = 0;
#if defined(__AVX2__)
        // 两半各自连续：每次处理 2 个复振幅对
        const double m0r = m0.real(), m0i = m0.imag();
        const double m1r = m1.real(), m1i = m1.imag();
        const double m2r = m2.real(), m2i = m2.imag();
        const double m3r = m3.real(), m3i = m3.imag();
        for (; off + 2 <= stride; off += 2) {
            const __m256d a0 = _mm256_loadu_pd(reinterpret_cast<double*>(p0 + off));
            const __m256d a1 = _mm256_loadu_pd(reinterpret_cast<double*>(p1 + off));
            const __m256d r0 = simd::cmul_add(simd::cmul(a0, m0r, m0i), a1, m1r, m1i);
            const __m256d r1 = simd::cmul_add(simd::cmul(a0, m2r, m2i), a1, m3r, m3i);
            _mm256_storeu_pd(reinterpret_cast<double*>(p0 + off), r0);
            _mm256_storeu_pd(reinterpret_cast<double*>(p1 + off), r1);
        }
#endif
        for (; off < stride; ++off) {
            const complex_t a0 = p0[off];
            const complex_t a1 = p1[off];
            p0[off] = m0 * a0 + m1 * a1;
            p1[off] = m2 * a0 + m3 * a1;
        }
    }
}

void apply_diagonal_single_qubit(complex_t* state, std::size_t n_amps,
                                 qubit_t target, complex_t d0, complex_t d1) {
    const std::size_t stride = std::size_t{1} << target;
    const std::size_t block = 2 * stride;
    const std::size_t n_blocks = n_amps / block;
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t b = 0; b < n_blocks; ++b) {
        complex_t* p0 = state + b * block;  // |0> 半 × d0
        complex_t* p1 = p0 + stride;        // |1> 半 × d1
        std::size_t off = 0;
#if defined(__AVX2__)
        const double d0r = d0.real(), d0i = d0.imag();
        const double d1r = d1.real(), d1i = d1.imag();
        for (; off + 2 <= stride; off += 2) {
            _mm256_storeu_pd(
                reinterpret_cast<double*>(p0 + off),
                simd::cmul(_mm256_loadu_pd(reinterpret_cast<double*>(p0 + off)),
                           d0r, d0i));
            _mm256_storeu_pd(
                reinterpret_cast<double*>(p1 + off),
                simd::cmul(_mm256_loadu_pd(reinterpret_cast<double*>(p1 + off)),
                           d1r, d1i));
        }
#endif
        for (; off < stride; ++off) {
            p0[off] *= d0;
            p1[off] *= d1;
        }
    }
}

void apply_controlled_single_qubit(complex_t* state, std::size_t n_amps,
                                   qubit_t control, qubit_t target,
                                   const complex_t m[4]) {
    const unsigned lo = std::min(control, target);
    const unsigned hi = std::max(control, target);
    const std::size_t c_mask = std::size_t{1} << control;
    const std::size_t t_mask = std::size_t{1} << target;
    const std::size_t n_groups = n_amps / 4;
    const complex_t m0 = m[0], m1 = m[1], m2 = m[2], m3 = m[3];
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t g = 0; g < n_groups; ++g) {
        const std::size_t base =
            bits::insert_zero_bit(bits::insert_zero_bit(g, lo), hi);
        const std::size_t i0 = base | c_mask;
        const std::size_t i1 = i0 | t_mask;
        const complex_t a0 = state[i0];
        const complex_t a1 = state[i1];
        state[i0] = m0 * a0 + m1 * a1;
        state[i1] = m2 * a0 + m3 * a1;
    }
}

void apply_swap(complex_t* state, std::size_t n_amps, qubit_t a, qubit_t b) {
    const unsigned lo = std::min(a, b);
    const unsigned hi = std::max(a, b);
    const std::size_t a_mask = std::size_t{1} << a;
    const std::size_t b_mask = std::size_t{1} << b;
    const std::size_t n_groups = n_amps / 4;
    // 仅 |..a=1,b=0..> 与 |..a=0,b=1..> 交换
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t g = 0; g < n_groups; ++g) {
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
    const std::size_t n_groups = n_amps / 4;
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t g = 0; g < n_groups; ++g) {
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

// 受控多比特门的公共准备工作：固定位、控制掩码、局部偏移表
struct ControlledLayout {
    std::vector<unsigned> fixed;            // 所有固定比特位置（升序）
    std::size_t c_mask = 0;                 // 控制位掩码
    std::vector<std::size_t> offset;        // 局部目标编号 l → 完整下标偏移
    std::size_t n_groups = 0;
};

ControlledLayout make_layout(std::size_t n_amps, const qubit_t* controls,
                             std::size_t n_controls, const qubit_t* targets,
                             std::size_t n_targets) {
    ControlledLayout lay;
    lay.fixed.assign(controls, controls + n_controls);
    lay.fixed.insert(lay.fixed.end(), targets, targets + n_targets);
    std::sort(lay.fixed.begin(), lay.fixed.end());
    for (std::size_t c = 0; c < n_controls; ++c) {
        lay.c_mask |= std::size_t{1} << controls[c];
    }
    const std::size_t dim = std::size_t{1} << n_targets;
    lay.offset.assign(dim, 0);
    for (std::size_t l = 0; l < dim; ++l) {
        for (std::size_t j = 0; j < n_targets; ++j) {
            if (bits::test_bit(l, static_cast<unsigned>(j))) {
                lay.offset[l] |= std::size_t{1} << targets[j];
            }
        }
    }
    lay.n_groups = n_amps >> lay.fixed.size();
    return lay;
}

}  // namespace

void apply_controlled_unitary(complex_t* state, std::size_t n_amps,
                              const qubit_t* controls, std::size_t n_controls,
                              const qubit_t* targets, std::size_t n_targets,
                              const complex_t* m) {
    const auto lay = make_layout(n_amps, controls, n_controls, targets, n_targets);
    const std::size_t dim = std::size_t{1} << n_targets;
#pragma omp parallel if (n_amps >= kParallelThreshold)
    {
        std::vector<complex_t> a(dim);   // 线程私有 scratch
        std::vector<std::size_t> idx(dim);
#pragma omp for schedule(static)
        for (std::size_t g = 0; g < lay.n_groups; ++g) {
            const std::size_t base = expand(g, lay.fixed) | lay.c_mask;
            for (std::size_t l = 0; l < dim; ++l) {
                idx[l] = base | lay.offset[l];
                a[l] = state[idx[l]];
            }
            for (std::size_t r = 0; r < dim; ++r) {
                complex_t sum{0, 0};
                for (std::size_t c = 0; c < dim; ++c) sum += m[r * dim + c] * a[c];
                state[idx[r]] = sum;
            }
        }
    }
}

void apply_controlled_permutation(complex_t* state, std::size_t n_amps,
                                  const qubit_t* controls, std::size_t n_controls,
                                  const qubit_t* targets, std::size_t n_targets,
                                  const std::size_t* perm) {
    const auto lay = make_layout(n_amps, controls, n_controls, targets, n_targets);
    const std::size_t dim = std::size_t{1} << n_targets;
#pragma omp parallel if (n_amps >= kParallelThreshold)
    {
        std::vector<complex_t> a(dim);
#pragma omp for schedule(static)
        for (std::size_t g = 0; g < lay.n_groups; ++g) {
            const std::size_t base = expand(g, lay.fixed) | lay.c_mask;
            for (std::size_t l = 0; l < dim; ++l) {
                a[l] = state[base | lay.offset[l]];
            }
            for (std::size_t l = 0; l < dim; ++l) {
                state[base | lay.offset[perm[l]]] = a[l];
            }
        }
    }
}

void apply_pauli_rotation(complex_t* state, std::size_t n_amps,
                          std::size_t x_global, std::size_t zy_global,
                          unsigned n_y, double theta) {
    const double h = theta / 2;
    const complex_t c = std::cos(h);
    const complex_t mis{0.0, -std::sin(h)};  // -i·sin(θ/2)
    static constexpr complex_t kIPow[4] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    const complex_t yphase = mis * kIPow[n_y % 4];

    if (x_global == 0) {
        // 对角情形（纯 I/Z）：|i> → (c + yphase·sign(i))·|i>
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
        for (std::size_t i = 0; i < n_amps; ++i) {
            const double sign =
                (std::popcount(i & zy_global) % 2 == 0) ? 1.0 : -1.0;
            state[i] *= c + yphase * sign;
        }
        return;
    }

    // 非对角：枚举振幅对 (i, i⊕x)，只在 i 的最低 X 位为 0 时处理避免重复
    const std::size_t lsb = x_global & (~x_global + 1);
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t i = 0; i < n_amps; ++i) {
        if (i & lsb) continue;
        const std::size_t j = i ^ x_global;
        const double si = (std::popcount(i & zy_global) % 2 == 0) ? 1.0 : -1.0;
        const double sj = (std::popcount(j & zy_global) % 2 == 0) ? 1.0 : -1.0;
        const complex_t ai = state[i];
        const complex_t aj = state[j];
        // a'_i = c·a_i + yphase·sign(j)·a_j;  a'_j = yphase·sign(i)·a_i + c·a_j
        state[i] = c * ai + yphase * sj * aj;
        state[j] = yphase * si * ai + c * aj;
    }
}

bool collapse_qubit(complex_t* state, std::size_t n_amps, qubit_t q,
                    double random01) {
    double p1 = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : p1) \
    if (n_amps >= kParallelThreshold)
    for (std::size_t i = 0; i < n_amps; ++i) {
        if (bits::test_bit(i, q)) p1 += std::norm(state[i]);
    }
    const bool outcome = random01 < p1;
    const double keep_p = outcome ? p1 : 1.0 - p1;
    const double scale = 1.0 / std::sqrt(keep_p);
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
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
    const std::size_t n_groups = n_amps / 2;
#pragma omp parallel for schedule(static) if (n_amps >= kParallelThreshold)
    for (std::size_t g = 0; g < n_groups; ++g) {
        const std::size_t i0 = bits::insert_zero_bit(g, q);
        state[i0] = state[i0 | mask];
        state[i0 | mask] = complex_t{0, 0};
    }
}

double kraus_probability(const complex_t* state, std::size_t n_amps,
                         qubit_t target, const complex_t k[4]) {
    const std::size_t n_groups = n_amps / 2;
    const complex_t k0 = k[0], k1 = k[1], k2 = k[2], k3 = k[3];
    double p = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : p) \
    if (n_amps >= kParallelThreshold)
    for (std::size_t g = 0; g < n_groups; ++g) {
        const std::size_t i0 = bits::insert_zero_bit(g, target);
        const std::size_t i1 = i0 | (std::size_t{1} << target);
        const complex_t a0 = state[i0];
        const complex_t a1 = state[i1];
        p += std::norm(k0 * a0 + k1 * a1) + std::norm(k2 * a0 + k3 * a1);
    }
    return p;
}

}  // namespace lqcs::kernels
