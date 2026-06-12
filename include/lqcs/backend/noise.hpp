#pragma once

#include <array>
#include <vector>

#include "../core/types.hpp"

namespace lqcs {

// 单比特 Kraus 通道：ρ → Σ_i K_i ρ K_i†（构造时校验 Σ K_i†K_i = I）
struct KrausChannel {
    std::vector<std::array<complex_t, 4>> ops;

    explicit KrausChannel(std::vector<std::array<complex_t, 4>> kraus_ops);

    // —— 常用噪声通道 ——
    // 去极化：ρ → (1-p)ρ + p·I/2
    static KrausChannel depolarizing(double p);
    // 振幅阻尼（能量弛豫 T1）：|1> 以概率 γ 衰减到 |0>
    static KrausChannel amplitude_damping(double gamma);
    // 相位阻尼（退相干 T2）：相干项以 √(1-γ) 衰减
    static KrausChannel phase_damping(double gamma);
    // 比特翻转：以概率 p 施加 X
    static KrausChannel bit_flip(double p);
    // 相位翻转：以概率 p 施加 Z
    static KrausChannel phase_flip(double p);
};

// 噪声模型（M6 范围）：每个门作用后，对其触及的每个比特施加同一通道
class NoiseModel {
public:
    NoiseModel() = default;

    NoiseModel& add_all_qubit_channel(KrausChannel channel) {
        channels_.push_back(std::move(channel));
        return *this;
    }

    const std::vector<KrausChannel>& channels() const { return channels_; }
    bool empty() const { return channels_.empty(); }

private:
    std::vector<KrausChannel> channels_;
};

}  // namespace lqcs
