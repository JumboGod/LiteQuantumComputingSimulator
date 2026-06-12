// 随机电路等价性：特化内核 vs 朴素稠密矩阵参考实现，结果必须一致
#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include <lqcs/lqcs.hpp>

using namespace lqcs;

namespace {

constexpr double kTol = 1e-10;

// 参考实现：把指令的完整矩阵（含控制位，局部位 j ↔ qubits[j]）
// 朴素地嵌入全空间作用——慢但显然正确
void naive_apply(std::vector<complex_t>& state, const Instruction& inst) {
    const auto m = inst.gate.matrix();
    const std::size_t k = inst.qubits.size();
    const std::size_t dim = std::size_t{1} << k;

    std::vector<complex_t> out(state.size(), complex_t{0, 0});
    for (std::size_t idx = 0; idx < state.size(); ++idx) {
        // 提取局部列编号并清掉相关位
        std::size_t col = 0;
        std::size_t base = idx;
        for (std::size_t j = 0; j < k; ++j) {
            if ((idx >> inst.qubits[j]) & 1u) {
                col |= std::size_t{1} << j;
                base &= ~(std::size_t{1} << inst.qubits[j]);
            }
        }
        for (std::size_t row = 0; row < dim; ++row) {
            const complex_t coeff = m[row * dim + col];
            if (coeff == complex_t{0, 0}) continue;
            std::size_t target = base;
            for (std::size_t j = 0; j < k; ++j) {
                if ((row >> j) & 1u) target |= std::size_t{1} << inst.qubits[j];
            }
            out[target] += coeff * state[idx];
        }
    }
    state = std::move(out);
}

// 生成覆盖所有门种类的随机电路
QuantumCircuit random_circuit(std::size_t n, std::size_t n_gates,
                              std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> angle(0.0, 6.28);
    auto pick_qubits = [&](std::size_t count) {
        std::vector<qubit_t> qs(n);
        for (std::size_t i = 0; i < n; ++i) qs[i] = static_cast<qubit_t>(i);
        std::shuffle(qs.begin(), qs.end(), rng);
        qs.resize(count);
        return qs;
    };

    QuantumCircuit qc(n);
    for (std::size_t i = 0; i < n_gates; ++i) {
        const int choice = static_cast<int>(rng() % 16);
        switch (choice) {
            case 0:  qc.h(pick_qubits(1)[0]); break;
            case 1:  qc.x(pick_qubits(1)[0]); break;
            case 2:  qc.t(pick_qubits(1)[0]); break;
            case 3:  qc.sx(pick_qubits(1)[0]); break;
            case 4:  qc.rx(angle(rng), pick_qubits(1)[0]); break;
            case 5:  qc.u(angle(rng), angle(rng), angle(rng), pick_qubits(1)[0]); break;
            case 6:  { auto q = pick_qubits(2); qc.cx(q[0], q[1]); } break;
            case 7:  { auto q = pick_qubits(2); qc.cp(angle(rng), q[0], q[1]); } break;
            case 8:  { auto q = pick_qubits(2); qc.swap(q[0], q[1]); } break;
            case 9:  { auto q = pick_qubits(2); qc.iswap(q[0], q[1]); } break;
            case 10: { auto q = pick_qubits(2); qc.rxx(angle(rng), q[0], q[1]); } break;
            case 11: { auto q = pick_qubits(2); qc.rzz(angle(rng), q[0], q[1]); } break;
            case 12: { auto q = pick_qubits(3); qc.ccx(q[0], q[1], q[2]); } break;
            case 13: { auto q = pick_qubits(3); qc.cswap(q[0], q[1], q[2]); } break;
            case 14: {
                auto q = pick_qubits(3);
                qc.mcp(angle(rng), {q[0], q[1]}, q[2]);
                break;
            }
            default: {
                // 2 比特随机置换
                auto q = pick_qubits(2);
                std::vector<std::size_t> table = {0, 1, 2, 3};
                std::shuffle(table.begin(), table.end(), rng);
                qc.permutation(table, {q[0], q[1]});
                break;
            }
        }
    }
    return qc;
}

}  // namespace

TEST(Equivalence, KernelsMatchNaiveReferenceOnRandomCircuits) {
    for (std::uint64_t seed : {11u, 22u, 33u}) {
        const std::size_t n = 5;
        QuantumCircuit qc = random_circuit(n, 40, seed);

        auto sv = StatevectorSimulator().run_statevector(qc);

        std::vector<complex_t> ref(std::size_t{1} << n, complex_t{0, 0});
        ref[0] = complex_t{1, 0};
        for (const auto& inst : qc.instructions()) naive_apply(ref, inst);

        SCOPED_TRACE("seed " + std::to_string(seed));
        for (std::size_t i = 0; i < ref.size(); ++i) {
            EXPECT_NEAR(std::abs(sv[i] - ref[i]), 0.0, kTol) << "amplitude " << i;
        }
    }
}

TEST(Equivalence, CircuitTimesItsInverseIsIdentity) {
    QuantumCircuit qc = random_circuit(4, 30, 77);
    qc.compose(qc.inverse());
    auto sv = StatevectorSimulator().run_statevector(qc);
    EXPECT_NEAR(std::abs(sv[0] - complex_t{1, 0}), 0.0, kTol);
    for (std::size_t i = 1; i < sv.size(); ++i) {
        EXPECT_NEAR(std::abs(sv[i]), 0.0, kTol);
    }
}

TEST(Equivalence, NormPreservedOnRandomCircuit) {
    auto sv = StatevectorSimulator().run_statevector(random_circuit(5, 60, 99));
    EXPECT_NEAR(sv.norm(), 1.0, 1e-11);
}
