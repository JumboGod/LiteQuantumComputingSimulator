#include "lqcs/backend/stabilizer.hpp"

#include <stdexcept>
#include <utility>

namespace lqcs {

StabilizerTableau::StabilizerTableau(std::size_t num_qubits) : n_(num_qubits) {
    if (num_qubits == 0) {
        throw std::invalid_argument("StabilizerTableau: num_qubits must be >= 1");
    }
    const std::size_t rows = 2 * n_ + 1;  // 含 scratch 行
    xs_.assign(rows * n_, 0);
    zs_.assign(rows * n_, 0);
    rs_.assign(rows, 0);
    // destabilizers = X_i (行 i)，stabilizers = Z_i (行 n+i)
    for (std::size_t i = 0; i < n_; ++i) {
        xs_[i * n_ + i] = 1;
        zs_[(n_ + i) * n_ + i] = 1;
    }
}

void StabilizerTableau::h(qubit_t q) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) {
        std::uint8_t& xi = xat(i, q);
        std::uint8_t& zi = zat(i, q);
        rs_[i] ^= xi & zi;
        std::swap(xi, zi);
    }
}

void StabilizerTableau::s(qubit_t q) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) {
        std::uint8_t& xi = xat(i, q);
        std::uint8_t& zi = zat(i, q);
        rs_[i] ^= xi & zi;
        zi ^= xi;
    }
}

void StabilizerTableau::sdg(qubit_t q) {
    s(q);  // Sdg = S^3
    s(q);
    s(q);
}

void StabilizerTableau::sx(qubit_t q) {
    h(q);  // SX = H·S·H（差全局相位，stabilizer 无关）
    s(q);
    h(q);
}

void StabilizerTableau::x(qubit_t q) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) rs_[i] ^= zat(i, q);  // X·Z·X = -Z
}

void StabilizerTableau::z(qubit_t q) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) rs_[i] ^= xat(i, q);  // Z·X·Z = -X
}

void StabilizerTableau::y(qubit_t q) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) rs_[i] ^= xat(i, q) ^ zat(i, q);
}

void StabilizerTableau::cx(qubit_t a, qubit_t b) {
    const std::size_t rows = 2 * n_;
    for (std::size_t i = 0; i < rows; ++i) {
        std::uint8_t& xa = xat(i, a);
        std::uint8_t& xb = xat(i, b);
        std::uint8_t& za = zat(i, a);
        std::uint8_t& zb = zat(i, b);
        rs_[i] ^= xa & zb & (xb ^ za ^ 1);
        xb ^= xa;
        za ^= zb;
    }
}

void StabilizerTableau::cz(qubit_t a, qubit_t b) {
    h(b);
    cx(a, b);
    h(b);
}

void StabilizerTableau::cy(qubit_t a, qubit_t b) {
    sdg(b);  // CY = (I⊗S)·CX·(I⊗Sdg)
    cx(a, b);
    s(b);
}

void StabilizerTableau::swap(qubit_t a, qubit_t b) {
    cx(a, b);
    cx(b, a);
    cx(a, b);
}

void StabilizerTableau::rowsum(std::size_t h, std::size_t i) {
    // g(x1,z1,x2,z2)：两 Pauli 相乘的相位指数贡献
    auto g = [](std::uint8_t x1, std::uint8_t z1, std::uint8_t x2,
                std::uint8_t z2) -> int {
        if (x1 == 0 && z1 == 0) return 0;
        if (x1 == 1 && z1 == 1) return static_cast<int>(z2) - static_cast<int>(x2);
        if (x1 == 1 && z1 == 0) return z2 * (2 * x2 - 1);
        return x2 * (1 - 2 * z2);  // x1==0, z1==1
    };
    int sum = 2 * rs_[h] + 2 * rs_[i];
    for (std::size_t j = 0; j < n_; ++j) {
        sum += g(xat(i, j), zat(i, j), xat(h, j), zat(h, j));
    }
    sum %= 4;
    if (sum < 0) sum += 4;
    rs_[h] = (sum == 2) ? 1 : 0;  // sum ∈ {0,2}
    for (std::size_t j = 0; j < n_; ++j) {
        xat(h, j) ^= xat(i, j);
        zat(h, j) ^= zat(i, j);
    }
}

bool StabilizerTableau::measure(qubit_t a, double random01) {
    // 找 stabilizer 行（n..2n-1）中 x_{p,a}==1
    std::size_t p = 2 * n_;
    bool found = false;
    for (std::size_t i = n_; i < 2 * n_; ++i) {
        if (xat(i, a)) { p = i; found = true; break; }
    }

    if (found) {
        // 随机结果
        for (std::size_t i = 0; i < 2 * n_; ++i) {
            if (i != p && xat(i, a)) rowsum(i, p);
        }
        // destabilizer 行 (p-n) ← 旧 stabilizer 行 p
        const std::size_t d = p - n_;
        for (std::size_t j = 0; j < n_; ++j) {
            xat(d, j) = xat(p, j);
            zat(d, j) = zat(p, j);
        }
        rs_[d] = rs_[p];
        // 行 p ← Z_a，符号为随机测量结果
        for (std::size_t j = 0; j < n_; ++j) { xat(p, j) = 0; zat(p, j) = 0; }
        zat(p, a) = 1;
        const std::uint8_t outcome = (random01 < 0.5) ? 0 : 1;
        rs_[p] = outcome;
        return outcome != 0;
    }

    // 确定性结果：用 scratch 行 2n
    const std::size_t scratch = 2 * n_;
    for (std::size_t j = 0; j < n_; ++j) { xat(scratch, j) = 0; zat(scratch, j) = 0; }
    rs_[scratch] = 0;
    for (std::size_t i = 0; i < n_; ++i) {
        if (xat(i, a)) rowsum(scratch, i + n_);
    }
    return rs_[scratch] != 0;
}

void StabilizerTableau::reset(qubit_t q, double random01) {
    if (measure(q, random01)) x(q);
}

bool StabilizerTableau::is_clifford(const Gate& g) {
    if (!g.is_unitary()) return false;
    if (g.n_controls == 0) {
        switch (g.type) {
            case GateType::I:
            case GateType::X:
            case GateType::Y:
            case GateType::Z:
            case GateType::H:
            case GateType::S:
            case GateType::Sdg:
            case GateType::SX:
            case GateType::SWAP:
                return true;
            default:
                return false;
        }
    }
    if (g.n_controls == 1) {
        return g.type == GateType::X || g.type == GateType::Y ||
               g.type == GateType::Z;
    }
    return false;
}

void StabilizerTableau::apply(const Instruction& inst) {
    const Gate& g = inst.gate;
    const auto& q = inst.qubits;
    if (g.n_controls == 1) {
        switch (g.type) {
            case GateType::X: cx(q[0], q[1]); return;
            case GateType::Y: cy(q[0], q[1]); return;
            case GateType::Z: cz(q[0], q[1]); return;
            default: break;
        }
    } else if (g.n_controls == 0) {
        switch (g.type) {
            case GateType::I: return;
            case GateType::X: x(q[0]); return;
            case GateType::Y: y(q[0]); return;
            case GateType::Z: z(q[0]); return;
            case GateType::H: h(q[0]); return;
            case GateType::S: s(q[0]); return;
            case GateType::Sdg: sdg(q[0]); return;
            case GateType::SX: sx(q[0]); return;
            case GateType::SWAP: swap(q[0], q[1]); return;
            default: break;
        }
    }
    throw std::invalid_argument("StabilizerSimulator: gate '" + g.name() +
                                "' is not a Clifford gate");
}

std::vector<std::string> StabilizerTableau::stabilizers() const {
    std::vector<std::string> out;
    out.reserve(n_);
    for (std::size_t i = n_; i < 2 * n_; ++i) {
        std::string s = rs_[i] ? "-" : "+";
        for (std::size_t j = 0; j < n_; ++j) {
            const std::uint8_t xx = xs_[i * n_ + j], zz = zs_[i * n_ + j];
            s += (xx && zz) ? 'Y' : xx ? 'X' : zz ? 'Z' : 'I';
        }
        out.push_back(std::move(s));
    }
    return out;
}

}  // namespace lqcs
