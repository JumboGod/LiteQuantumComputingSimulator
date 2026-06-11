// GHZ 态制备与测量：(|000> + |111>)/√2
// 预期输出：约一半 "000"、一半 "111"
#include <iostream>

#include <lqcs/lqcs.hpp>

int main() {
    using namespace lqcs;

    QuantumCircuit qc(3, 3);
    qc.h(0).cx(0, 1).cx(1, 2).measure_all();

    StatevectorSimulator sim;
    Result result = sim.run(qc, 1024);

    std::cout << "GHZ state (|000> + |111>)/sqrt(2), 1024 shots:\n";
    result.print_counts();
    return 0;
}
