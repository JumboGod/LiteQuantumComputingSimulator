// 含噪 Bell 态：去极化噪声下保真度随噪声强度的衰减
#include <cstdio>

#include <lqcs/lqcs.hpp>

using namespace lqcs;

int main() {
    // 理想 Bell 态（保真度参考）
    QuantumCircuit bell(2, 2);
    bell.h(0).cx(0, 1);
    const Statevector ideal = StatevectorSimulator().run_statevector(bell);

    std::printf("Bell state under depolarizing noise (channel after each gate)\n");
    std::printf("  p        fidelity   purity\n");
    for (const double p : {0.0, 0.01, 0.05, 0.1, 0.2}) {
        NoiseModel noise;
        if (p > 0) noise.add_all_qubit_channel(KrausChannel::depolarizing(p));
        DensityMatrixSimulator sim({.seed = 7, .noise = noise});
        const DensityMatrix rho = sim.run_density_matrix(bell);
        std::printf("  %-7.3g  %.6f   %.6f\n", p, rho.fidelity(ideal),
                    rho.purity());
    }

    // 含噪采样
    bell.measure_all();
    NoiseModel noise;
    noise.add_all_qubit_channel(KrausChannel::depolarizing(0.1));
    Result r = DensityMatrixSimulator({.seed = 7, .noise = noise}).run(bell, 2048);
    std::printf("\ncounts at p = 0.1 (errors leak into 01/10):\n");
    r.print_counts();
    return 0;
}
