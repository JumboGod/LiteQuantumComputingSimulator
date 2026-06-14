#pragma once

#include <cstdint>
#include <optional>

#include "backend.hpp"
#include "noise.hpp"
#include "statevector.hpp"

namespace lqcs {

// 状态向量模拟器：
//  - 末尾测量：演化一次到末态，对概率分布采样 shots 次
//  - 中间测量/Reset：逐 shot 演化坍缩
//  - 含噪声模型时走量子轨迹法（per-shot 采样 Kraus 分支）：
//    内存 O(2^n)，统计上与密度矩阵精确解一致
class StatevectorSimulator final : public Backend {
   public:
    struct Options {
        std::optional<std::uint64_t> seed;  // 固定种子使采样可复现
        bool fuse_gates = true;  // 自动门融合（省全向量扫描）
        std::size_t fusion_max_qubits = 4;  // 融合块最大比特数（qsim 风格）
        int num_threads = 0;                // OpenMP 线程数（0 = 默认）
        // 噪声模型（轨迹法）。非空时：每个幺正门后对触及比特采样通道分支；
        // 自动禁用门融合（噪声绑定在门上，融合会改变物理语义）
        NoiseModel noise;
    };

    StatevectorSimulator() = default;
    explicit StatevectorSimulator(Options options) : options_(options) {}

    Result run(const QuantumCircuit& circuit,
               std::size_t shots = 1024) override;
    std::string name() const override { return "statevector_simulator"; }

    // 返回末态（跳过 Measure/Barrier），用于直接检查振幅
    Statevector run_statevector(const QuantumCircuit& circuit) const;

   private:
    void configure_threads() const;
    QuantumCircuit maybe_fuse(const QuantumCircuit& circuit) const;

    Options options_;
};

}  // namespace lqcs
