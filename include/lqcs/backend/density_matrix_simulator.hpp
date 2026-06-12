#pragma once

#include <cstdint>
#include <optional>

#include "backend.hpp"
#include "density_matrix.hpp"
#include "noise.hpp"

namespace lqcs {

// 密度矩阵模拟器：支持噪声模型（每个门后对触及比特施加 Kraus 通道）。
// 内存 O(4^n)，适合 ≤ 13 比特的含噪模拟。仅支持末尾测量。
class DensityMatrixSimulator final : public Backend {
public:
    struct Options {
        std::optional<std::uint64_t> seed;
        NoiseModel noise;  // 空 = 理想模拟
    };

    DensityMatrixSimulator() = default;
    explicit DensityMatrixSimulator(Options options)
        : options_(std::move(options)) {}

    Result run(const QuantumCircuit& circuit, std::size_t shots = 1024) override;
    std::string name() const override { return "density_matrix_simulator"; }

    // 返回末态 ρ（跳过 Measure/Barrier；噪声照常施加）
    DensityMatrix run_density_matrix(const QuantumCircuit& circuit) const;

private:
    Options options_;
};

}  // namespace lqcs
