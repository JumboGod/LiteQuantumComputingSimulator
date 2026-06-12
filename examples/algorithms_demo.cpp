// 经典量子算法演示：Deutsch-Jozsa / Bernstein-Vazirani / Simon / QPE / Grover
#include <bit>
#include <iostream>
#include <numbers>
#include <vector>

#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

int main() {
    std::cout << "1) Deutsch-Jozsa: is f constant or balanced? (1 query)\n";
    auto dj_c = deutsch_jozsa(4, [](std::uint64_t) { return true; });
    auto dj_b = deutsch_jozsa(4, [](std::uint64_t x) {
        return std::popcount(x) % 2 == 1;
    });
    std::cout << "   f(x) = 1       -> " << (dj_c.is_constant ? "constant" : "balanced")
              << "\n   f(x) = parity  -> " << (dj_b.is_constant ? "constant" : "balanced")
              << "\n\n";

    std::cout << "2) Bernstein-Vazirani: recover hidden s from f(x) = s.x (1 query)\n";
    const std::uint64_t s_bv = 0b101101;
    auto bv = bernstein_vazirani(6, [](std::uint64_t x) {
        return std::popcount(x & 0b101101ULL) % 2 == 1;
    });
    std::cout << "   hidden s = " << s_bv << ", recovered = " << bv.secret << "\n\n";

    std::cout << "3) Simon: find s with f(x) = f(x^s) (exponential speedup)\n";
    const std::uint64_t s_simon = 0b1101;
    auto sm = simon(4, [s_simon](std::uint64_t x) {
        return std::min(x, x ^ s_simon);
    }, 42);
    std::cout << "   hidden s = " << s_simon << ", recovered = " << sm.secret
              << " (" << sm.quantum_samples << " quantum samples)\n\n";

    std::cout << "4) Quantum Phase Estimation: T|1> = e^{i*pi/4}|1>\n";
    QuantumCircuit prep(1);
    prep.x(0);
    auto qpe = phase_estimation(Gate{GateType::T, {}}, prep, 5);
    std::cout << "   estimated phase = " << qpe.phase << " (exact: 0.125)\n\n";

    std::cout << "5) Grover: search marked |42> among 2^7 = 128 states\n";
    const std::vector<std::uint64_t> marked = {42};
    auto gr = grover(7, marked);
    std::cout << "   best state = " << gr.best_state << ", P(success) = "
              << gr.success_probability << " after " << gr.iterations
              << " iterations\n";
    return 0;
}
