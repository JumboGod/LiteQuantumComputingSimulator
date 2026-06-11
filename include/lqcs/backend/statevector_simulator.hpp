#pragma once

#include <cstdint>
#include <optional>

#include "backend.hpp"
#include "statevector.hpp"

namespace lqcs {

// 状态向量模拟器（M1）：
//  - 支持电路末尾测量：演化一次到末态，再对概率分布采样 shots 次
//  - 中间测量（测量后还有门）留待 M2
class StatevectorSimulator final : public Backend {
public:
    struct Options {
        std::optional<std::uint64_t> seed;  // 固定种子使采样可复现
    };

    StatevectorSimulator() = default;
    explicit StatevectorSimulator(Options options) : options_(options) {}

    Result run(const QuantumCircuit& circuit, std::size_t shots = 1024) override;
    std::string name() const override { return "statevector_simulator"; }

    // 返回末态（跳过 Measure/Barrier），用于直接检查振幅
    Statevector run_statevector(const QuantumCircuit& circuit) const;

private:
    Options options_;
};

}  // namespace lqcs
