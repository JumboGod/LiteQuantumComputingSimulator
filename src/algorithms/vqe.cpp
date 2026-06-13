#include "lqcs/algorithms/vqe.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "lqcs/backend/statevector_simulator.hpp"

namespace lqcs::algorithms {

double expectation(const Statevector& sv, const Hamiltonian& hamiltonian) {
    double e = 0.0;
    for (const auto& term : hamiltonian) {
        e += term.coeff * sv.expectation_value(term.pauli);
    }
    return e;
}

VQEResult vqe(const Hamiltonian& hamiltonian, const Ansatz& ansatz,
              std::size_t n_params, const VQEOptions& options) {
    if (n_params == 0) {
        throw std::invalid_argument("vqe: n_params must be >= 1");
    }
    std::vector<double> theta = options.initial_parameters;
    if (theta.empty()) {
        theta.assign(n_params, 0.0);
    } else if (theta.size() != n_params) {
        throw std::invalid_argument(
            "vqe: initial_parameters size must equal n_params");
    }

    StatevectorSimulator sim;
    const auto energy = [&](std::span<const double> params) {
        return expectation(sim.run_statevector(ansatz(params)), hamiltonian);
    };

    constexpr double kHalfPi = std::numbers::pi / 2;
    VQEResult result;
    double current = energy(theta);

    for (std::size_t iter = 0; iter < options.max_iterations; ++iter) {
        const double before = current;
        for (std::size_t k = 0; k < n_params; ++k) {
            // Rotosolve：E(θ_k) = a + b·cos(θ_k − φ)，三点解析定位最小值
            const double e0 = current;
            const double saved = theta[k];
            theta[k] = saved + kHalfPi;
            const double e_plus = energy(theta);
            theta[k] = saved - kHalfPi;
            const double e_minus = energy(theta);

            theta[k] = saved - kHalfPi -
                       std::atan2(2.0 * e0 - e_plus - e_minus, e_plus - e_minus);
            // 折回 (-π, π]，避免参数无界漂移
            theta[k] = std::remainder(theta[k], 2 * std::numbers::pi);
            current = energy(theta);
        }
        result.history.push_back(current);
        result.iterations = iter + 1;
        if (before - current < options.tol) break;
    }

    result.energy = current;
    result.parameters = std::move(theta);
    return result;
}

std::vector<double> parameter_shift_gradient(const Hamiltonian& hamiltonian,
                                             const ParametricCircuit& ansatz,
                                             std::span<const double> params) {
    const std::size_t n = ansatz.num_parameters();
    if (params.size() != n) {
        throw std::invalid_argument(
            "parameter_shift_gradient: params size must equal num_parameters");
    }
    StatevectorSimulator sim;
    const auto energy = [&](const std::vector<double>& p) {
        return expectation(sim.run_statevector(ansatz.bind(p)), hamiltonian);
    };

    constexpr double kHalfPi = std::numbers::pi / 2;
    std::vector<double> theta(params.begin(), params.end());
    std::vector<double> grad(n);
    for (std::size_t k = 0; k < n; ++k) {
        const double saved = theta[k];
        theta[k] = saved + kHalfPi;
        const double e_plus = energy(theta);
        theta[k] = saved - kHalfPi;
        const double e_minus = energy(theta);
        theta[k] = saved;
        grad[k] = 0.5 * (e_plus - e_minus);
    }
    return grad;
}

VQEResult vqe_gradient_descent(const Hamiltonian& hamiltonian,
                               const ParametricCircuit& ansatz,
                               const GradientVQEOptions& options) {
    const std::size_t n = ansatz.num_parameters();
    if (n == 0) {
        throw std::invalid_argument("vqe_gradient_descent: ansatz has no parameters");
    }
    std::vector<double> theta = options.initial_parameters;
    if (theta.empty()) {
        theta.assign(n, 0.0);
    } else if (theta.size() != n) {
        throw std::invalid_argument(
            "vqe_gradient_descent: initial_parameters size must equal "
            "num_parameters");
    }

    StatevectorSimulator sim;
    const auto energy = [&](const std::vector<double>& p) {
        return expectation(sim.run_statevector(ansatz.bind(p)), hamiltonian);
    };

    VQEResult result;
    for (std::size_t iter = 0; iter < options.max_iterations; ++iter) {
        const auto grad = parameter_shift_gradient(hamiltonian, ansatz, theta);
        double grad_norm = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            theta[k] -= options.learning_rate * grad[k];
            grad_norm += grad[k] * grad[k];
        }
        result.history.push_back(energy(theta));
        result.iterations = iter + 1;
        if (std::sqrt(grad_norm) < options.tol) break;
    }

    result.energy = result.history.empty() ? energy(theta) : result.history.back();
    result.parameters = std::move(theta);
    return result;
}

}  // namespace lqcs::algorithms
