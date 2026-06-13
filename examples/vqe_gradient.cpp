// 参数化电路 + parameter-shift 梯度下降求解 H2 基态。
// 对比 Qulacs 风格工作流：ParametricCircuit 构建一次，梯度驱动优化。
#include <cstdio>

#include <lqcs/lqcs.hpp>

using namespace lqcs;
using namespace lqcs::algorithms;

int main() {
    const Hamiltonian h2 = {
        {-1.052373245772859, "II"},   {+0.39793742484318045, "IZ"},
        {-0.39793742484318045, "ZI"}, {-0.01128010425623538, "ZZ"},
        {+0.18093119978423156, "XX"},
    };

    // 参数化 ansatz：X(0)·RY(θ,1)·CX(1,0)（覆盖基态所在奇宇称子空间）
    ParametricCircuit ansatz(2);
    ansatz.x(0);
    ansatz.ry(1);
    ansatz.cx(1, 0);
    std::printf("ansatz: %zu parameters\n\n", ansatz.num_parameters());

    // parameter-shift 梯度（解析）
    const std::vector<double> theta0 = {0.1};
    const auto grad = parameter_shift_gradient(h2, ansatz, theta0);
    std::printf("parameter-shift gradient at theta=0.1: dE/dtheta = %+.6f\n\n",
                grad[0]);

    // 梯度下降优化
    auto res = vqe_gradient_descent(
        h2, ansatz, {.max_iterations = 200, .learning_rate = 0.3, .tol = 1e-8});

    std::printf("gradient-descent VQE:\n");
    for (std::size_t i = 0; i < res.history.size(); i += 20) {
        std::printf("  iter %3zu: E = %+.10f\n", i, res.history[i]);
    }
    std::printf("converged: E = %+.10f Hartree (%zu iters, theta = %.6f)\n",
                res.energy, res.iterations, res.parameters[0]);
    std::printf("exact:     E = -1.8572750302 Hartree   error = %.2e\n",
                res.energy + 1.8572750302);
    return 0;
}
