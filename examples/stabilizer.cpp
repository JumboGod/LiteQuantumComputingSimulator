// Stabilizer 模拟器：Clifford 电路 O(n²) 模拟，规模远超状态向量法。
#include <cstdio>
#include <string>

#include <lqcs/lqcs.hpp>

using namespace lqcs;

int main() {
    // 1) Bell 态的 stabilizer 群
    {
        QuantumCircuit qc(2);
        qc.h(0).cx(0, 1);
        auto tab = StabilizerSimulator().run_tableau(qc);
        std::printf("Bell state stabilizers:\n");
        for (const auto& g : tab.stabilizers())
            std::printf("  %s\n", g.c_str());
    }

    // 2) 大规模 GHZ：状态向量法不可能，stabilizer 法 O(n²) 完成
    //    （测量两端比特即可看出全局关联，避免 O(n) 次全比特测量的开销）
    std::printf("\nGHZ state correlation (stabilizer simulator):\n");
    for (const std::size_t n : {100u, 1000u, 5000u}) {
        QuantumCircuit qc(static_cast<qubit_t>(n), 2);
        qc.h(0);
        for (qubit_t q = 0; q + 1 < n; ++q) qc.cx(q, q + 1);
        qc.measure(0, 0);
        qc.measure(static_cast<qubit_t>(n - 1), 1);  // 首尾两比特
        auto counts = StabilizerSimulator({.seed = 42}).run(qc, 8);
        std::printf("  n=%4zu: outcomes(q0,q_last) = ", n);
        for (const auto& [key, c] : counts.counts()) {
            std::printf("%s x%zu  ", key.c_str(),
                        c);  // 应只有 00/11（完全关联）
        }
        std::printf("(state vector would need 2^%zu amplitudes)\n", n);
    }

    // 3) 5 比特 [[5,1,3]] 码的稳定子（Clifford 编码电路片段示意）
    std::printf("\nrandom Clifford on 50 qubits, 1000 gates: ");
    QuantumCircuit qc(50, 50);
    for (qubit_t q = 0; q < 50; ++q) qc.h(q);
    for (qubit_t q = 0; q + 1 < 50; ++q) qc.cx(q, q + 1);
    for (qubit_t q = 0; q < 50; ++q) qc.s(q);
    qc.measure_all();
    auto r = StabilizerSimulator({.seed = 7}).run(qc, 4);
    std::printf("%zu shots sampled ok\n", r.shots());
    return 0;
}
