// 量子傅里叶变换：对基态 |x> 作用 QFT，
// 末态振幅应为 e^{2πi·x·k/2^n} / √(2^n)（所有模长相等，相位线性递增）
#include <cmath>
#include <cstdio>
#include <iostream>
#include <numbers>

#include <lqcs/lqcs.hpp>

using namespace lqcs;

// 标准 QFT 电路（little-endian，含末尾翻转 SWAP）
QuantumCircuit build_qft(std::size_t n) {
    QuantumCircuit qc(n);
    for (std::size_t i = n; i-- > 0;) {
        qc.h(static_cast<qubit_t>(i));
        for (std::size_t j = i; j-- > 0;) {
            const double angle =
                std::numbers::pi / static_cast<double>(1ULL << (i - j));
            qc.cp(angle, static_cast<qubit_t>(j), static_cast<qubit_t>(i));
        }
    }
    for (std::size_t i = 0; i < n / 2; ++i) {
        qc.swap(static_cast<qubit_t>(i), static_cast<qubit_t>(n - 1 - i));
    }
    return qc;
}

int main() {
    const std::size_t n = 3;
    const std::size_t x = 5;  // 输入基态 |101>

    QuantumCircuit qc(n);
    qc.x(0).x(2);  // 制备 |101> (qubit0 与 qubit2 置 1)
    qc.compose(build_qft(n));

    std::cout << "QFT(" << n << ") applied to |" << x << ">:\n\n"
              << qc.draw() << "\n";

    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t dim = sv.size();
    std::cout << "amplitudes (expected |amp| = " << 1.0 / std::sqrt(dim)
              << ", phase = 2*pi*" << x << "*k/" << dim << "):\n";
    for (std::size_t k = 0; k < dim; ++k) {
        std::printf("  |%zu>: %+.4f%+.4fi   |amp| = %.4f  phase = %+.4f\n", k,
                    sv[k].real(), sv[k].imag(), std::abs(sv[k]),
                    std::arg(sv[k]));
    }
    return 0;
}
