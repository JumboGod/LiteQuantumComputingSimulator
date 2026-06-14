// QFT 与解析 DFT 公式对比：QFT|x> = (1/√N) Σ_k e^{2πi·x·k/N} |k>
#include <cmath>
#include <numbers>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-12;

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

void check_qft_on_basis_state(std::size_t n, std::size_t x) {
    QuantumCircuit qc(n);
    for (std::size_t b = 0; b < n; ++b) {
        if ((x >> b) & 1u) qc.x(static_cast<qubit_t>(b));
    }
    qc.compose(build_qft(n));

    auto sv = StatevectorSimulator().run_statevector(qc);
    const std::size_t dim = std::size_t{1} << n;
    const double amp = 1.0 / std::sqrt(static_cast<double>(dim));
    for (std::size_t k = 0; k < dim; ++k) {
        const double phase = 2.0 * std::numbers::pi *
                             static_cast<double>(x * k) /
                             static_cast<double>(dim);
        const complex_t expected =
            amp * complex_t{std::cos(phase), std::sin(phase)};
        EXPECT_NEAR(std::abs(sv[k] - expected), 0.0, kTol)
            << "n=" << n << " x=" << x << " k=" << k;
    }
}

}  // namespace

TEST(QFT, MatchesAnalyticDFT) {
    check_qft_on_basis_state(3, 0);
    check_qft_on_basis_state(3, 1);
    check_qft_on_basis_state(3, 5);
    check_qft_on_basis_state(4, 7);
    check_qft_on_basis_state(4, 13);
    check_qft_on_basis_state(5, 21);
}

TEST(QFT, InverseUndoesQFT) {
    const std::size_t n = 4;
    QuantumCircuit qc(n);
    qc.x(1).x(3);  // |1010> = |10>
    const auto qft = build_qft(n);
    qc.compose(qft).compose(qft.inverse());

    auto sv = StatevectorSimulator().run_statevector(qc);
    EXPECT_NEAR(std::abs(sv[10] - complex_t{1, 0}), 0.0, kTol);
}
