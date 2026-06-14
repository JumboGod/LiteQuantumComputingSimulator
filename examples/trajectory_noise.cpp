// 轨迹法噪声模拟：
// 1) 小规模下与密度矩阵精确解对照（验证统计一致性）
// 2) 16 比特含噪 GHZ（密度矩阵后端无法企及的规模）
#include <cstdio>
#include <string>

#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

NoiseModel depolarizing(double p) {
    NoiseModel n;
    n.add_all_qubit_channel(KrausChannel::depolarizing(p));
    return n;
}

}  // namespace

int main() {
    // —— 1. 交叉验证：轨迹法 vs 密度矩阵 ——
    const double p = 0.1;
    QuantumCircuit bell(2, 2);
    bell.h(0).cx(0, 1).measure_all();

    const std::size_t shots = 50000;
    auto traj = StatevectorSimulator({.seed = 7, .noise = depolarizing(p)})
                    .run(bell, shots);
    auto rho = DensityMatrixSimulator({.noise = depolarizing(p)})
                   .run_density_matrix(bell);
    const auto exact = rho.probabilities();

    std::printf(
        "Bell + depolarizing(p=%.2f): trajectory (%zu shots) vs exact\n", p,
        shots);
    std::printf("  state   trajectory   density-matrix\n");
    const char* labels[4] = {"00", "01", "10", "11"};
    for (int i = 0; i < 4; ++i) {
        const auto it = traj.counts().find(labels[i]);
        const double freq =
            it == traj.counts().end()
                ? 0.0
                : static_cast<double>(it->second) / static_cast<double>(shots);
        std::printf("   %s     %.4f       %.4f\n", labels[i], freq, exact[i]);
    }

    // —— 2. 大规模含噪：16 比特 GHZ 噪声扫描 ——
    const std::size_t n = 16;
    std::printf(
        "\nGHZ(%zu) under depolarizing noise (O(2^n) trajectory method;\n"
        "density matrix would need O(4^n) = 64 TiB):\n",
        n);
    std::printf("  p        P(all-0) + P(all-1)\n");
    const std::string zeros(n, '0'), ones(n, '1');
    for (const double pn : {0.0, 0.002, 0.01, 0.03}) {
        QuantumCircuit ghz(n, n);
        ghz.h(0);
        for (qubit_t q = 0; q + 1 < n; ++q) ghz.cx(q, q + 1);
        ghz.measure_all();
        auto r = StatevectorSimulator({.seed = 11, .noise = depolarizing(pn)})
                     .run(ghz, 512);
        std::size_t good = 0;
        if (auto it = r.counts().find(zeros); it != r.counts().end())
            good += it->second;
        if (auto it = r.counts().find(ones); it != r.counts().end())
            good += it->second;
        std::printf("  %-7.3g  %.3f\n", pn, good / 512.0);
    }
    return 0;
}
