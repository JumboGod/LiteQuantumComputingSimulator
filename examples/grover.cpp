// Grover 搜索：3 比特空间中搜索被标记的 |101>。
// 2 次迭代后测得 "101" 的概率约 94.5%。
#include <iostream>
#include <numbers>

#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

const double kPi = std::numbers::pi;

// Oracle：对 |101> 翻转相位（先把目标态映射到 |111>，再做 CCZ）
void oracle_101(QuantumCircuit& qc) {
    qc.x(1);
    qc.mcp(kPi, {0, 1}, 2);  // CCZ
    qc.x(1);
}

// 扩散算子：关于均匀叠加态的反射
void diffusion(QuantumCircuit& qc) {
    for (qubit_t q = 0; q < 3; ++q) qc.h(q);
    for (qubit_t q = 0; q < 3; ++q) qc.x(q);
    qc.mcp(kPi, {0, 1}, 2);
    for (qubit_t q = 0; q < 3; ++q) qc.x(q);
    for (qubit_t q = 0; q < 3; ++q) qc.h(q);
}

}  // namespace

int main() {
    QuantumCircuit qc(3, 3);
    for (qubit_t q = 0; q < 3; ++q) qc.h(q);
    for (int iter = 0; iter < 2; ++iter) {  // round(pi/4 * sqrt(8)) = 2
        oracle_101(qc);
        diffusion(qc);
    }
    qc.measure_all();

    StatevectorSimulator sim;
    Result result = sim.run(qc, 1024);

    std::cout << "Grover search for |101> in 3 qubits, 2 iterations, 1024 shots:\n";
    result.print_counts();
    std::cout << "(theory: P(101) ~ 94.5%)\n";
    return 0;
}
