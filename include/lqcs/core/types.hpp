#pragma once

#include <complex>
#include <cstdint>
#include <vector>

#include "aligned_allocator.hpp"

namespace lqcs {

using complex_t = std::complex<double>;
using qubit_t = std::uint32_t;
using clbit_t = std::uint32_t;

template <typename T>
using aligned_vector = std::vector<T, AlignedAllocator<T, 64>>;

}  // namespace lqcs
