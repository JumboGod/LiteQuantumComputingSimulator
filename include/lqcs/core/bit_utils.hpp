#pragma once

#include <cstddef>

namespace lqcs::bits {

// 在 pos 位插入一个 0 比特：低 pos 位保持不动，pos 及以上整体左移一位。
// 用于把「自由比特的紧凑编号」展开成「指定比特为 0 的完整下标」。
constexpr std::size_t insert_zero_bit(std::size_t value, unsigned pos) {
    const std::size_t low_mask = (std::size_t{1} << pos) - 1;
    return ((value & ~low_mask) << 1) | (value & low_mask);
}

constexpr bool test_bit(std::size_t value, unsigned pos) {
    return (value >> pos) & std::size_t{1};
}

}  // namespace lqcs::bits
