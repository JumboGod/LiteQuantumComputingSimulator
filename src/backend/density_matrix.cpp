#include "lqcs/backend/density_matrix.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "kernels/kernels.hpp"

namespace lqcs {

namespace {

std::array<complex_t, 4> conj2(const std::array<complex_t, 4>& m) {
    return {std::conj(m[0]), std::conj(m[1]), std::conj(m[2]), std::conj(m[3])};
}

}  // namespace

DensityMatrix::DensityMatrix(std::size_t num_qubits) : num_qubits_(num_qubits) {
    if (num_qubits == 0 || num_qubits > kMaxQubits) {
        throw std::invalid_argument(
            "DensityMatrix: num_qubits must be in [1, " +
            std::to_string(kMaxQubits) + "] (memory is O(4^n))");
    }
    data_.assign(std::size_t{1} << (2 * num_qubits), complex_t{0, 0});
    data_[0] = complex_t{1, 0};  // |0...0><0...0|
}

complex_t DensityMatrix::element(std::size_t row, std::size_t col) const {
    // |ρ>> 编号：ket 在低 n 位，bra 在高 n 位；ρ_{rc} = <r|ρ|c>
    return data_[(col << num_qubits_) | row];
}

void DensityMatrix::apply_instruction(const Instruction& inst) {
    const Gate& g = inst.gate;
    if (!g.is_unitary()) {
        throw std::invalid_argument("DensityMatrix: '" + g.name() +
                                    "' is not unitary");
    }
    const std::size_t nc = g.n_controls;
    const std::size_t k = g.base_qubits();
    const std::size_t n_amps = data_.size();

    // ket 侧比特即原比特；bra 侧比特 = 原比特 + n
    std::vector<qubit_t> ket(inst.qubits);
    std::vector<qubit_t> bra(inst.qubits);
    for (auto& q : bra) q += static_cast<qubit_t>(num_qubits_);

    if (g.type == GateType::Permutation) {
        kernels::apply_controlled_permutation(data_.data(), n_amps, ket.data(),
                                              nc, ket.data() + nc, k,
                                              g.perm.data());
        kernels::apply_controlled_permutation(data_.data(), n_amps, bra.data(),
                                              nc, bra.data() + nc, k,
                                              g.perm.data());
        return;
    }
    const auto m = g.base_matrix();
    std::vector<complex_t> mc(m.size());
    for (std::size_t i = 0; i < m.size(); ++i) mc[i] = std::conj(m[i]);
    // ρ → UρU†：ket 侧 U，bra 侧 U*
    kernels::apply_controlled_unitary(data_.data(), n_amps, ket.data(), nc,
                                      ket.data() + nc, k, m.data());
    kernels::apply_controlled_unitary(data_.data(), n_amps, bra.data(), nc,
                                      bra.data() + nc, k, mc.data());
}

void DensityMatrix::apply_channel(const KrausChannel& channel, qubit_t q) {
    if (q >= num_qubits_) {
        throw std::out_of_range("apply_channel: qubit index out of range");
    }
    const auto bra_q = static_cast<qubit_t>(q + num_qubits_);
    // ρ → Σ_i K_i ρ K_i†：逐项作用后累加
    aligned_vector<complex_t> acc(data_.size(), complex_t{0, 0});
    aligned_vector<complex_t> work(data_.size());
    for (const auto& kraus : channel.ops) {
        work = data_;
        const auto kc = conj2(kraus);
        kernels::apply_single_qubit(work.data(), work.size(), q, kraus.data());
        kernels::apply_single_qubit(work.data(), work.size(), bra_q, kc.data());
        for (std::size_t i = 0; i < acc.size(); ++i) acc[i] += work[i];
    }
    data_ = std::move(acc);
}

double DensityMatrix::trace() const {
    const std::size_t d = dim();
    double tr = 0.0;
    for (std::size_t i = 0; i < d; ++i) tr += element(i, i).real();
    return tr;
}

double DensityMatrix::purity() const {
    // Tr(ρ²) = Σ_{rc} |ρ_{rc}|²（ρ 厄米）
    double sum = 0.0;
    for (const auto& v : data_) sum += std::norm(v);
    return sum;
}

std::vector<double> DensityMatrix::probabilities() const {
    const std::size_t d = dim();
    std::vector<double> probs(d);
    for (std::size_t i = 0; i < d; ++i) probs[i] = element(i, i).real();
    return probs;
}

double DensityMatrix::fidelity(const Statevector& psi) const {
    if (psi.num_qubits() != num_qubits_) {
        throw std::invalid_argument("fidelity: qubit count mismatch");
    }
    // <ψ|ρ|ψ> = Σ_{rc} conj(ψ_r)·ρ_{rc}·ψ_c
    const std::size_t d = dim();
    complex_t f{0, 0};
    for (std::size_t c = 0; c < d; ++c) {
        complex_t col{0, 0};
        for (std::size_t r = 0; r < d; ++r) {
            col += std::conj(psi[r]) * element(r, c);
        }
        f += col * psi[c];
    }
    return f.real();
}

}  // namespace lqcs
