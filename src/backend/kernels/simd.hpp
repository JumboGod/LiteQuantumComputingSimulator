#pragma once

#include "lqcs/core/types.hpp"

// 复数 SIMD 辅助：AVX2 下一个 __m256d 打包 2 个 complex<double>
// （内存布局 [re0, im0, re1, im1]）。无 AVX2 时本头不提供任何符号，
// 内核走标量回退路径。
#if defined(__AVX2__)
#include <immintrin.h>

namespace lqcs::kernels::simd {

// 打包向量 × 标量复数 (mr + i·mi)：每个 lane 做复数乘法
//   (re + i·im)(mr + i·mi) = (re·mr − im·mi) + i·(im·mr + re·mi)
inline __m256d cmul(__m256d a, double mr, double mi) {
    const __m256d vmr = _mm256_set1_pd(mr);
    const __m256d vmi = _mm256_set1_pd(mi);
    const __m256d a_swap = _mm256_permute_pd(a, 0x5);  // [im,re,im,re]
    // addsub(a·mr, a_swap·mi) = [re·mr − im·mi, im·mr + re·mi, ...]
    return _mm256_addsub_pd(_mm256_mul_pd(a, vmr), _mm256_mul_pd(a_swap, vmi));
}

// fused 复数乘加：acc + a·(mr + i·mi)
inline __m256d cmul_add(__m256d acc, __m256d a, double mr, double mi) {
    return _mm256_add_pd(acc, cmul(a, mr, mi));
}

}  // namespace lqcs::kernels::simd
#endif
