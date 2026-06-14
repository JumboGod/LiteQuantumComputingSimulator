// VQE 求解 H2 分子基态能量（键长 0.735 Å，STO-3G，约化为 2 比特）。
// 精确基态能量 ≈ -1.857275 Hartree。
#include <cstdio>
#include <span>

#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

int main() {
    // Jordan-Wigner 变换后的 2 比特 H2 哈密顿量（Qiskit 教科书标准值）
    const Hamiltonian h2 = {
        {-1.052373245772859, "II"},   {+0.39793742484318045, "IZ"},
        {-0.39793742484318045, "ZI"}, {-0.01128010425623538, "ZZ"},
        {+0.18093119978423156, "XX"},
    };

    // 教科书标准 H2 ansatz（单参数，精确覆盖基态所在子空间）：
    //   X(0)        → Hartree-Fock 态 |01>
    //   RY(θ, 1)    → cos(θ/2)|01> + sin(θ/2)|11>
    //   CX(1, 0)    → cos(θ/2)|01> + sin(θ/2)|10>
    // 基态 ≈ 0.994|01> − 0.111|10> 就在该曲线上
    const Ansatz ansatz = [](std::span<const double> theta) {
        QuantumCircuit qc(2);
        qc.x(0);
        qc.ry(theta[0], 1);
        qc.cx(1, 0);
        return qc;
    };

    VQEResult res = vqe(h2, ansatz, 1, {.max_iterations = 50, .tol = 1e-12});

    std::printf("VQE ground state of H2 (d = 0.735 A)\n");
    std::printf("------------------------------------\n");
    for (std::size_t i = 0; i < res.history.size(); ++i) {
        std::printf("  sweep %2zu: E = %+.10f\n", i + 1, res.history[i]);
    }
    std::printf("converged: E = %+.10f Hartree (%zu sweeps)\n", res.energy,
                res.iterations);
    std::printf("exact:     E = -1.8572750302 Hartree\n");
    std::printf("error:     %.2e\n", res.energy + 1.8572750302);
    return 0;
}
