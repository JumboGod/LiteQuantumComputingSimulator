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

}  // namespace lqcs::algorithms
