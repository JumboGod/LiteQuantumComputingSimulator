#pragma once

#include <functional>
#include <span>
#include <string>
#include <vector>

#include "../backend/statevector.hpp"
#include "../circuit/quantum_circuit.hpp"

namespace lqcs::algorithms {

// 哈密顿量 = Pauli 串的实系数线性组合
struct PauliTerm {
    double      coeff;
    std::string pauli;  // 与 Qiskit 一致：最左字符 = 最高位比特
};
using Hamiltonian = std::vector<PauliTerm>;

// ⟨ψ|H|ψ⟩
double expectation(const Statevector& sv, const Hamiltonian& hamiltonian);

// ansatz：参数向量 → 制备 |ψ(θ)> 的电路
using Ansatz = std::function<QuantumCircuit(std::span<const double>)>;

struct VQEOptions {
    std::size_t max_iterations = 100;   // 最多多少轮全参数扫描
    double      tol = 1e-9;             // 一轮能量改善小于 tol 即收敛
    std::vector<double> initial_parameters;  // 缺省为全 0
};

struct VQEResult {
    double              energy;       // 收敛能量（变分上界）
    std::vector<double> parameters;   // 最优参数
    std::size_t         iterations;   // 实际轮数
    std::vector<double> history;      // 每轮结束时的能量
};

// 变分量子本征求解器，Rotosolve 优化器：
// 假定每个参数都是旋转角（RX/RY/RZ/CRX… 形式，能量对单参数为正弦曲线），
// 每个坐标用 3 次能量评估解析求出最小值，无需梯度与学习率。
VQEResult vqe(const Hamiltonian& hamiltonian, const Ansatz& ansatz,
              std::size_t n_params, const VQEOptions& options = {});

}  // namespace lqcs::algorithms
