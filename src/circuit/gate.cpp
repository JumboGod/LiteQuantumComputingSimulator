#include "lqcs/circuit/gate.hpp"

#include <bit>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace lqcs {

namespace {

constexpr complex_t kI{0.0, 1.0};
const double kSqrt1_2 = std::numbers::sqrt2 / 2.0;

complex_t expi(double theta) { return {std::cos(theta), std::sin(theta)}; }

std::vector<complex_t> conj_transpose(const std::vector<complex_t>& m,
                                      std::size_t dim) {
    std::vector<complex_t> d(dim * dim);
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j)
            d[j * dim + i] = std::conj(m[i * dim + j]);
    return d;
}

}  // namespace

std::size_t Gate::base_qubits() const {
    switch (type) {
        case GateType::SWAP:
        case GateType::iSWAP:
        case GateType::RXX:
        case GateType::RYY:
        case GateType::RZZ: return 2;
        case GateType::Unitary: {
            // mat 为 dim×dim，dim = 2^k
            const std::size_t dim =
                static_cast<std::size_t>(std::sqrt(mat.size()));
            return static_cast<std::size_t>(std::countr_zero(dim));
        }
        case GateType::Permutation:
            return static_cast<std::size_t>(std::countr_zero(perm.size()));
        case GateType::PauliRotation: return pauli.size();
        case GateType::Barrier: return 0;
        default: return 1;  // 单比特门、Measure、Reset
    }
}

Gate::PauliMasks Gate::pauli_masks() const {
    PauliMasks m;
    const std::size_t k = pauli.size();
    for (std::size_t j = 0; j < k; ++j) {
        const char c = pauli[k - 1 - j];  // 局部位 j ↔ pauli 最右起第 j 个
        const std::size_t bit = std::size_t{1} << j;
        switch (c) {
            case 'I': break;
            case 'X': m.x |= bit; break;
            case 'Y':
                m.x |= bit;
                m.zy |= bit;
                ++m.n_y;
                break;
            case 'Z': m.zy |= bit; break;
            default:
                throw std::invalid_argument(
                    "PauliRotation: pauli must contain only I/X/Y/Z");
        }
    }
    return m;
}

bool Gate::is_unitary() const {
    return type != GateType::Measure && type != GateType::Reset &&
           type != GateType::Barrier;
}

bool Gate::is_diagonal() const {
    switch (type) {
        case GateType::I:
        case GateType::Z:
        case GateType::S:
        case GateType::Sdg:
        case GateType::T:
        case GateType::Tdg:
        case GateType::RZ:
        case GateType::P:
        case GateType::RZZ: return true;
        case GateType::PauliRotation:
            return pauli_masks().x == 0;  // 仅含 I/Z 时对角
        default: return false;
    }
}

std::vector<complex_t> Gate::base_matrix() const {
    if (!is_unitary()) {
        throw std::logic_error("Gate::base_matrix: " + name() +
                               " has no matrix");
    }
    switch (type) {
        case GateType::I: return {1, 0, 0, 1};
        case GateType::X: return {0, 1, 1, 0};
        case GateType::Y: return {0, -kI, kI, 0};
        case GateType::Z: return {1, 0, 0, -1};
        case GateType::H: return {kSqrt1_2, kSqrt1_2, kSqrt1_2, -kSqrt1_2};
        case GateType::S: return {1, 0, 0, kI};
        case GateType::Sdg: return {1, 0, 0, -kI};
        case GateType::T: return {1, 0, 0, expi(std::numbers::pi / 4)};
        case GateType::Tdg: return {1, 0, 0, expi(-std::numbers::pi / 4)};
        case GateType::SX:
            return {complex_t{0.5, 0.5}, complex_t{0.5, -0.5},
                    complex_t{0.5, -0.5}, complex_t{0.5, 0.5}};
        case GateType::RX: {
            const double h = params.at(0) / 2;
            return {std::cos(h), -kI * std::sin(h), -kI * std::sin(h),
                    std::cos(h)};
        }
        case GateType::RY: {
            const double h = params.at(0) / 2;
            return {std::cos(h), -std::sin(h), std::sin(h), std::cos(h)};
        }
        case GateType::RZ: {
            const double h = params.at(0) / 2;
            return {expi(-h), 0, 0, expi(h)};
        }
        case GateType::P: return {1, 0, 0, expi(params.at(0))};
        case GateType::U: {
            // U(θ,φ,λ) = [[cos(θ/2), -e^{iλ}·sin(θ/2)],
            //             [e^{iφ}·sin(θ/2), e^{i(φ+λ)}·cos(θ/2)]]
            const double h = params.at(0) / 2;
            const double phi = params.at(1), lam = params.at(2);
            return {std::cos(h), -expi(lam) * std::sin(h),
                    expi(phi) * std::sin(h), expi(phi + lam) * std::cos(h)};
        }
        case GateType::SWAP:
            return {1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1};
        case GateType::iSWAP:
            return {1, 0, 0, 0, 0, 0, kI, 0, 0, kI, 0, 0, 0, 0, 0, 1};
        case GateType::RXX: {
            const double h = params.at(0) / 2;
            const complex_t c = std::cos(h), is = kI * std::sin(h);
            return {c, 0, 0, -is, 0, c, -is, 0, 0, -is, c, 0, -is, 0, 0, c};
        }
        case GateType::RYY: {
            const double h = params.at(0) / 2;
            const complex_t c = std::cos(h), is = kI * std::sin(h);
            return {c, 0, 0, is, 0, c, -is, 0, 0, -is, c, 0, is, 0, 0, c};
        }
        case GateType::RZZ: {
            const double h = params.at(0) / 2;
            return {expi(-h), 0, 0,       0, 0, expi(h), 0, 0,
                    0,        0, expi(h), 0, 0, 0,       0, expi(-h)};
        }
        case GateType::PauliRotation: {
            // exp(-iθ/2·P) = cos(θ/2)·I − i·sin(θ/2)·P
            // P|col> = i^{nY}·(-1)^{popcount(col & zy)}·|col ⊕ x>
            const std::size_t k = pauli.size();
            const std::size_t dim = std::size_t{1} << k;
            const auto m = pauli_masks();
            const double h = params.at(0) / 2;
            const complex_t c = std::cos(h);
            const complex_t mis{0.0, -std::sin(h)};  // -i·sin
            static constexpr complex_t kIPow[4] = {
                {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
            const complex_t yphase = kIPow[m.n_y % 4];
            std::vector<complex_t> out(dim * dim, complex_t{0, 0});
            for (std::size_t col = 0; col < dim; ++col) {
                out[col * dim + col] += c;
                const std::size_t row = col ^ m.x;
                const double sign =
                    (static_cast<unsigned>(std::popcount(col & m.zy)) % 2 == 0)
                        ? 1.0
                        : -1.0;
                out[row * dim + col] += mis * yphase * sign;
            }
            return out;
        }
        case GateType::Unitary: return mat;
        case GateType::Permutation: {
            // |l> → |perm[l]>：列 l 在行 perm[l] 处为 1
            const std::size_t dim = perm.size();
            std::vector<complex_t> m(dim * dim, complex_t{0, 0});
            for (std::size_t l = 0; l < dim; ++l) m[perm[l] * dim + l] = 1;
            return m;
        }
        default:
            throw std::logic_error("Gate::base_matrix: unhandled gate type");
    }
}

std::vector<complex_t> Gate::matrix() const {
    const auto base = base_matrix();
    if (n_controls == 0) return base;

    const std::size_t k = base_qubits();
    const std::size_t base_dim = std::size_t{1} << k;
    const std::size_t dim = std::size_t{1} << (n_controls + k);
    const std::size_t c_mask = (std::size_t{1} << n_controls) - 1;

    // 控制位占局部低位：控制位非全 1 的列为恒等，全 1 的列作用基门
    std::vector<complex_t> m(dim * dim, complex_t{0, 0});
    for (std::size_t col = 0; col < dim; ++col) {
        if ((col & c_mask) != c_mask) {
            m[col * dim + col] = 1;
            continue;
        }
        const std::size_t t_col = col >> n_controls;
        for (std::size_t t_row = 0; t_row < base_dim; ++t_row) {
            const std::size_t row = (t_row << n_controls) | c_mask;
            m[row * dim + col] = base[t_row * base_dim + t_col];
        }
    }
    return m;
}

Gate Gate::inverse() const {
    if (!is_unitary()) {
        throw std::logic_error("Gate::inverse: " + name() +
                               " is not invertible");
    }
    Gate inv = *this;
    switch (type) {
        // 自逆门
        case GateType::I:
        case GateType::X:
        case GateType::Y:
        case GateType::Z:
        case GateType::H:
        case GateType::SWAP: return inv;
        // 参数取反
        case GateType::RX:
        case GateType::RY:
        case GateType::RZ:
        case GateType::P:
        case GateType::RXX:
        case GateType::RYY:
        case GateType::RZZ:
        case GateType::PauliRotation:
            inv.params[0] = -inv.params[0];
            return inv;
        case GateType::S: inv.type = GateType::Sdg; return inv;
        case GateType::Sdg: inv.type = GateType::S; return inv;
        case GateType::T: inv.type = GateType::Tdg; return inv;
        case GateType::Tdg: inv.type = GateType::T; return inv;
        case GateType::U:
            // U(θ,φ,λ)† = U(-θ,-λ,-φ)
            inv.params = {-params.at(0), -params.at(2), -params.at(1)};
            return inv;
        case GateType::SX:
        case GateType::iSWAP: {
            // 无专属枚举的逆：退化为自定义幺正（共轭转置）
            const std::size_t dim = std::size_t{1} << base_qubits();
            inv.type = GateType::Unitary;
            inv.mat = conj_transpose(base_matrix(), dim);
            return inv;
        }
        case GateType::Unitary: {
            const std::size_t dim = std::size_t{1} << base_qubits();
            inv.mat = conj_transpose(mat, dim);
            return inv;
        }
        case GateType::Permutation: {
            for (std::size_t l = 0; l < perm.size(); ++l) inv.perm[perm[l]] = l;
            return inv;
        }
        default: throw std::logic_error("Gate::inverse: unhandled gate type");
    }
}

std::string Gate::name() const {
    std::string base;
    switch (type) {
        case GateType::I: base = "id"; break;
        case GateType::X: base = "x"; break;
        case GateType::Y: base = "y"; break;
        case GateType::Z: base = "z"; break;
        case GateType::H: base = "h"; break;
        case GateType::S: base = "s"; break;
        case GateType::Sdg: base = "sdg"; break;
        case GateType::T: base = "t"; break;
        case GateType::Tdg: base = "tdg"; break;
        case GateType::SX: base = "sx"; break;
        case GateType::RX: base = "rx"; break;
        case GateType::RY: base = "ry"; break;
        case GateType::RZ: base = "rz"; break;
        case GateType::P: base = "p"; break;
        case GateType::U: base = "u"; break;
        case GateType::SWAP: base = "swap"; break;
        case GateType::iSWAP: base = "iswap"; break;
        case GateType::RXX: base = "rxx"; break;
        case GateType::RYY: base = "ryy"; break;
        case GateType::RZZ: base = "rzz"; break;
        case GateType::PauliRotation: base = "rp(" + pauli + ")"; break;
        case GateType::Unitary: base = "unitary"; break;
        case GateType::Permutation: base = "perm"; break;
        case GateType::Measure: base = "measure"; break;
        case GateType::Reset: base = "reset"; break;
        case GateType::Barrier: base = "barrier"; break;
    }
    if (n_controls == 0) return base;
    if (n_controls == 1) return "c" + base;
    if (n_controls == 2) return "cc" + base;
    return "mc" + base;
}

}  // namespace lqcs
