#pragma once

#include <cstdint>
#include <optional>
#include <random>

namespace lqcs {

// 可设种子的随机数引擎：固定种子时测量采样完全可复现
class Rng {
public:
    explicit Rng(std::optional<std::uint64_t> seed = std::nullopt)
        : engine_(seed ? *seed : std::random_device{}()) {}

    // [0, 1) 上均匀分布
    double uniform() { return dist_(engine_); }

    std::mt19937_64& engine() { return engine_; }

private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> dist_{0.0, 1.0};
};

}  // namespace lqcs
