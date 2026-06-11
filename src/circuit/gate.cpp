#include "lqcs/circuit/gate.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace lqcs {

namespace {

constexpr complex_t kI{0.0, 1.0};
const double kSqrt1_2 = std::numbers::sqrt2 / 2.0;

complex_t expi(double theta) { return {std::cos(theta), std::sin(theta)}; }

// 单比特门的 2×2 矩阵（行主序）
std::vector<complex_t> matrix_1q(GateType type, const std::vector<double>& params) {
    switch (type) {
        case GateType::I:   return {1, 0, 0, 1};
        case GateType::X:   return {0, 1, 1, 0};
        case GateType::Y:   return {0, -kI, kI, 0};
        case GateType::Z:   return {1, 0, 0, -1};
        case GateType::H:   return {kSqrt1_2, kSqrt1_2, kSqrt1_2, -kSqrt1_2};
        case GateType::S:   return {1, 0, 0, kI};
        case GateType::Sdg: return {1, 0, 0, -kI};
        case GateType::T:   return {1, 0, 0, expi(std::numbers::pi / 4)};
        case GateType::Tdg: return {1, 0, 0, expi(-std::numbers::pi / 4)};
        case GateType::SX:
            return {complex_t{0.5, 0.5}, complex_t{0.5, -0.5},
                    complex_t{0.5, -0.5}, complex_t{0.5, 0.5}};
        case GateType::RX: {
            const double h = params.at(0) / 2;
            return {std::cos(h), -kI * std::sin(h), -kI * std::sin(h), std::cos(h)};
        }
        case GateType::RY: {
            const double h = params.at(0) / 2;
            return {std::cos(h), -std::sin(h), std::sin(h), std::cos(h)};
        }
        case GateType::RZ: {
            const double h = params.at(0) / 2;
            return {expi(-h), 0, 0, expi(h)};
        }
        case GateType::P:
            return {1, 0, 0, expi(params.at(0))};
        case GateType::U: {
            // U(θ,φ,λ) = [[cos(θ/2), -e^{iλ}·sin(θ/2)],
            //             [e^{iφ}·sin(θ/2), e^{i(φ+λ)}·cos(θ/2)]]
            const double h = params.at(0) / 2;
            const double phi = params.at(1), lam = params.at(2);
            return {std::cos(h), -expi(lam) * std::sin(h),
                    expi(phi) * std::sin(h), expi(phi + lam) * std::cos(h)};
        }
        default:
            throw std::logic_error("matrix_1q: not a single-qubit gate");
    }
}

}  // namespace

std::size_t Gate::num_qubits() const {
    switch (type) {
        case GateType::CX:
        case GateType::CZ:
        case GateType::CP:
        case GateType::SWAP:
            return 2;
        case GateType::Measure:
            return 1;
        case GateType::Barrier:
            return 0;
        default:
            return 1;
    }
}

bool Gate::is_unitary() const {
    return type != GateType::Measure && type != GateType::Barrier;
}

std::vector<complex_t> Gate::matrix() const {
    if (!is_unitary()) {
        throw std::logic_error("Gate::matrix: " + name() + " has no matrix");
    }
    switch (type) {
        case GateType::CX:
            return {1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 0, 1,
                    0, 0, 1, 0};
        case GateType::CZ:
            return {1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, -1};
        case GateType::CP:
            return {1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, expi(params.at(0))};
        case GateType::SWAP:
            return {1, 0, 0, 0,
                    0, 0, 1, 0,
                    0, 1, 0, 0,
                    0, 0, 0, 1};
        default:
            return matrix_1q(type, params);
    }
}

std::string Gate::name() const {
    switch (type) {
        case GateType::I:       return "id";
        case GateType::X:       return "x";
        case GateType::Y:       return "y";
        case GateType::Z:       return "z";
        case GateType::H:       return "h";
        case GateType::S:       return "s";
        case GateType::Sdg:     return "sdg";
        case GateType::T:       return "t";
        case GateType::Tdg:     return "tdg";
        case GateType::SX:      return "sx";
        case GateType::RX:      return "rx";
        case GateType::RY:      return "ry";
        case GateType::RZ:      return "rz";
        case GateType::P:       return "p";
        case GateType::U:       return "u";
        case GateType::CX:      return "cx";
        case GateType::CZ:      return "cz";
        case GateType::CP:      return "cp";
        case GateType::SWAP:    return "swap";
        case GateType::Measure: return "measure";
        case GateType::Barrier: return "barrier";
    }
    return "unknown";
}

}  // namespace lqcs
