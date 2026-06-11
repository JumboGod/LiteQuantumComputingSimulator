// Bell 态制备与测量：|Φ+> = (|00> + |11>)/√2
// 预期输出：约一半 "00"、一半 "11"，不出现 "01"/"10"
#include <iostream>

#include <lqcs/lqcs.hpp>

int main() {
    using namespace lqcs;

    QuantumCircuit qc(2, 2);
    qc.h(0).cx(0, 1).measure_all();

    StatevectorSimulator sim;
    Result result = sim.run(qc, 1024);

    std::cout << "Bell state |Phi+> = (|00> + |11>)/sqrt(2), 1024 shots:\n";
    result.print_counts();
    return 0;
}
