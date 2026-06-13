#pragma once

#include <cstdint>
#include <optional>

#include "backend.hpp"
#include "stabilizer.hpp"

namespace lqcs {

// Stabilizer 模拟器：仅支持 Clifford 电路（H/S/CX/... + 测量/重置），
// 内存与时间 O(n²)，可模拟数千比特——状态向量法的指数内存在此规避。
// 适用于量子纠错码、GHZ/图态、随机 Clifford 基准等。
// 遇到非 Clifford 门（T、任意角旋转、Toffoli）抛异常。
class StabilizerSimulator final : public Backend {
public:
    struct Options {
        std::optional<std::uint64_t> seed;
    };

    StabilizerSimulator() = default;
    explicit StabilizerSimulator(Options options) : options_(options) {}

    Result run(const QuantumCircuit& circuit, std::size_t shots = 1024) override;
    std::string name() const override { return "stabilizer_simulator"; }

    // 演化一次并返回末态 tableau（含中间测量的随机坍缩）
    StabilizerTableau run_tableau(const QuantumCircuit& circuit) const;

private:
    Options options_;
};

}  // namespace lqcs
