#include "lqcs/backend/noise.hpp"

#include <cmath>
#include <stdexcept>

namespace lqcs {

namespace {

void check_probability(double p, const char* name) {
    if (p < 0.0 || p > 1.0) {
        throw std::invalid_argument(std::string(name) +
                                    ": probability must be in [0, 1]");
    }
}

}  // namespace

KrausChannel::KrausChannel(std::vector<std::array<complex_t, 4>> kraus_ops)
    : ops(std::move(kraus_ops)) {
    if (ops.empty()) {
        throw std::invalid_argument(
            "KrausChannel: needs at least one operator");
    }
    // 完备性校验：Σ K_i† K_i = I
    std::array<complex_t, 4> sum = {0, 0, 0, 0};
    for (const auto& k : ops) {
        // K†K
        sum[0] += std::conj(k[0]) * k[0] + std::conj(k[2]) * k[2];
        sum[1] += std::conj(k[0]) * k[1] + std::conj(k[2]) * k[3];
        sum[2] += std::conj(k[1]) * k[0] + std::conj(k[3]) * k[2];
        sum[3] += std::conj(k[1]) * k[1] + std::conj(k[3]) * k[3];
    }
    constexpr double kTol = 1e-10;
    if (std::abs(sum[0] - complex_t{1, 0}) > kTol || std::abs(sum[1]) > kTol ||
        std::abs(sum[2]) > kTol || std::abs(sum[3] - complex_t{1, 0}) > kTol) {
        throw std::invalid_argument(
            "KrausChannel: completeness violated (sum K_i^dag K_i != I)");
    }
}

KrausChannel KrausChannel::depolarizing(double p) {
    check_probability(p, "depolarizing");
    const double k0 = std::sqrt(1.0 - 3.0 * p / 4.0);
    const double kp = std::sqrt(p / 4.0);
    const complex_t i{0, 1};
    return KrausChannel({
        {k0, 0, 0, k0},           // √(1-3p/4)·I
        {0, kp, kp, 0},           // √(p/4)·X
        {0, -kp * i, kp * i, 0},  // √(p/4)·Y
        {kp, 0, 0, -kp},          // √(p/4)·Z
    });
}

KrausChannel KrausChannel::amplitude_damping(double gamma) {
    check_probability(gamma, "amplitude_damping");
    return KrausChannel({
        {1, 0, 0, std::sqrt(1.0 - gamma)},  // K0
        {0, std::sqrt(gamma), 0, 0},        // K1: |0><1|·√γ
    });
}

KrausChannel KrausChannel::phase_damping(double gamma) {
    check_probability(gamma, "phase_damping");
    return KrausChannel({
        {1, 0, 0, std::sqrt(1.0 - gamma)},
        {0, 0, 0, std::sqrt(gamma)},
    });
}

KrausChannel KrausChannel::bit_flip(double p) {
    check_probability(p, "bit_flip");
    const double a = std::sqrt(1.0 - p), b = std::sqrt(p);
    return KrausChannel({{a, 0, 0, a}, {0, b, b, 0}});
}

KrausChannel KrausChannel::phase_flip(double p) {
    check_probability(p, "phase_flip");
    const double a = std::sqrt(1.0 - p), b = std::sqrt(p);
    return KrausChannel({{a, 0, 0, a}, {b, 0, 0, -b}});
}

}  // namespace lqcs
