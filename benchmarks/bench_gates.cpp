// 性能基准：单门作用吞吐、QFT 端到端、门融合加速比
#include <chrono>
#include <cstdio>
#include <numbers>
#include <random>

#include <lqcs/lqcs.hpp>

using namespace lqcs;
using Clock = std::chrono::steady_clock;

namespace {

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

// 单门作用平均耗时（毫秒）
template <typename ApplyFn>
double bench_gate(std::size_t n, ApplyFn&& apply, int reps = 5) {
    QuantumCircuit qc(n);
    for (int r = 0; r < reps; ++r) apply(qc);
    StatevectorSimulator sim({.fuse_gates = false});
    const auto t0 = Clock::now();
    auto sv = sim.run_statevector(qc);
    (void)sv[0];
    return ms_since(t0) / reps;
}

void bench_single_gates(std::size_t n) {
    const auto lo = static_cast<qubit_t>(0);
    const auto hi = static_cast<qubit_t>(n - 1);
    const double h_lo = bench_gate(n, [&](auto& qc) { qc.h(lo); });
    const double h_hi = bench_gate(n, [&](auto& qc) { qc.h(hi); });
    const double t_diag = bench_gate(n, [&](auto& qc) { qc.t(hi); });
    const double cx = bench_gate(n, [&](auto& qc) { qc.cx(lo, hi); });
    const double ccx = bench_gate(n, [&](auto& qc) { qc.ccx(lo, 1, hi); });
    std::printf(
        "n=%2zu  | H(q0) %8.2f | H(q%zu) %8.2f | T %8.2f | CX %8.2f | CCX %8.2f\n",
        n, h_lo, n - 1, h_hi, t_diag, cx, ccx);
}

void bench_qft(std::size_t n) {
    QuantumCircuit qc(n);
    qc.compose(algorithms::qft(n));
    StatevectorSimulator sim;
    const auto t0 = Clock::now();
    auto sv = sim.run_statevector(qc);
    (void)sv[0];
    std::printf("QFT(%zu): %.1f ms  (%zu gates)\n", n, ms_since(t0),
                qc.num_instructions());
}

void bench_fusion(std::size_t n, std::size_t gates_per_qubit) {
    // 单比特门密集的随机电路：融合收益的典型场景
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> angle(0, 6.28);
    QuantumCircuit qc(n);
    for (std::size_t g = 0; g < gates_per_qubit; ++g) {
        for (qubit_t q = 0; q < n; ++q) {
            switch (rng() % 4) {
                case 0: qc.h(q); break;
                case 1: qc.rx(angle(rng), q); break;
                case 2: qc.t(q); break;
                default: qc.u(angle(rng), angle(rng), angle(rng), q); break;
            }
        }
        if (g % 8 == 7) {
            for (qubit_t q = 0; q + 1 < n; q += 2) qc.cx(q, q + 1);
        }
    }

    auto bench = [&](StatevectorSimulator sim) {
        const auto t0 = Clock::now();
        auto sv = sim.run_statevector(qc);
        (void)sv[0];
        return ms_since(t0);
    };
    const double t_plain = bench(StatevectorSimulator({.fuse_gates = false}));
    const double t_1q =
        bench(StatevectorSimulator({.fuse_gates = true, .fusion_max_qubits = 1}));
    const double t_k =
        bench(StatevectorSimulator({.fuse_gates = true, .fusion_max_qubits = 4}));

    std::printf(
        "fusion n=%zu (%zu gates): plain %.1f ms | 1-qubit %.1f ms (%.2fx) | "
        "block-4 %.1f ms (%.2fx)\n",
        n, qc.num_instructions(), t_plain, t_1q, t_plain / t_1q, t_k,
        t_plain / t_k);
}

}  // namespace

int main() {
    std::printf("== per-gate latency (ms, avg of 5) ==\n");
    for (std::size_t n : {20, 22, 24, 26}) bench_single_gates(n);

    std::printf("\n== QFT end-to-end ==\n");
    for (std::size_t n : {16, 20, 24}) bench_qft(n);

    std::printf("\n== single-qubit gate fusion ==\n");
    bench_fusion(20, 32);
    bench_fusion(22, 32);
    return 0;
}
