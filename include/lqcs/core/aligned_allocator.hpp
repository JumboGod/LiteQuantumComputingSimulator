#pragma once

#include <cstddef>
#include <new>

namespace lqcs {

// 按固定字节数对齐的分配器：保证状态向量按缓存行/SIMD 宽度对齐
template <typename T, std::size_t Alignment>
struct AlignedAllocator {
    static_assert((Alignment & (Alignment - 1)) == 0,
                  "Alignment must be a power of two");
    static_assert(Alignment >= alignof(T),
                  "Alignment must not be weaker than alignof(T)");

    using value_type = T;

    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    AlignedAllocator() noexcept = default;

    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    T* allocate(std::size_t n) {
        return static_cast<T*>(
            ::operator new(n * sizeof(T), std::align_val_t{Alignment}));
    }

    void deallocate(T* p, std::size_t) noexcept {
        ::operator delete(p, std::align_val_t{Alignment});
    }

    friend bool operator==(const AlignedAllocator&,
                           const AlignedAllocator&) noexcept {
        return true;
    }
};

}  // namespace lqcs
