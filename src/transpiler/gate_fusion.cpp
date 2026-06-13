#include "lqcs/transpiler/gate_fusion.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

namespace lqcs {

namespace {

using Mat2 = std::array<complex_t, 4>;

Mat2 matmul2(const Mat2& a, const Mat2& b) {
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

bool is_identity2(const Mat2& m) {
    constexpr double kTol = 1e-12;
    return std::abs(m[0] - complex_t{1, 0}) < kTol && std::abs(m[1]) < kTol &&
           std::abs(m[2]) < kTol && std::abs(m[3] - complex_t{1, 0}) < kTol;
}

bool is_fusable(const Instruction& inst) {
    return inst.gate.is_unitary() && inst.gate.n_controls == 0 &&
           inst.gate.base_qubits() == 1;
}

}  // namespace

QuantumCircuit SingleQubitGateFusion::run(const QuantumCircuit& circuit) const {
    QuantumCircuit out(circuit.num_qubits(), circuit.num_clbits());
    // 每个比特上待合并的 2×2 矩阵累积
    std::vector<std::optional<Mat2>> pending(circuit.num_qubits());

    auto flush = [&](qubit_t q) {
        if (!pending[q]) return;
        if (!is_identity2(*pending[q])) {
            out.append({Gate{GateType::Unitary, {}, 0,
                             {pending[q]->begin(), pending[q]->end()}},
                        {q},
                        {}});
        }
        pending[q].reset();
    };

    for (const auto& inst : circuit.instructions()) {
        if (is_fusable(inst)) {
            const qubit_t q = inst.qubits[0];
            const auto m = inst.gate.base_matrix();
            const Mat2 g = {m[0], m[1], m[2], m[3]};
            // 后作用的门左乘
            pending[q] = pending[q] ? matmul2(g, *pending[q]) : g;
            continue;
        }
        // Barrier 阻断所有比特的融合；其余指令只阻断其触及的比特
        if (inst.gate.type == GateType::Barrier) {
            for (qubit_t q = 0; q < circuit.num_qubits(); ++q) flush(q);
        } else {
            for (qubit_t q : inst.qubits) flush(q);
        }
        out.append(inst);
    }
    for (qubit_t q = 0; q < circuit.num_qubits(); ++q) flush(q);
    return out;
}

// ====================== 多比特门融合 ======================

namespace {

// C = A · B（dim×dim 行主序）
std::vector<complex_t> matmul(const std::vector<complex_t>& a,
                              const std::vector<complex_t>& b, std::size_t dim) {
    std::vector<complex_t> c(dim * dim, complex_t{0, 0});
    for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t l = 0; l < dim; ++l) {
            const complex_t ail = a[i * dim + l];
            if (ail == complex_t{0, 0}) continue;
            for (std::size_t j = 0; j < dim; ++j) c[i * dim + j] += ail * b[l * dim + j];
        }
    return c;
}

// 把作用在 src（局部位 j↔src[j]）上的矩阵嵌入到 dst（src⊆dst，均升序）。
// 非 src 比特上恒等。
std::vector<complex_t> expand(const std::vector<complex_t>& m,
                              const std::vector<qubit_t>& src,
                              const std::vector<qubit_t>& dst) {
    const std::size_t S = src.size(), D = dst.size();
    const std::size_t dim_s = std::size_t{1} << S, dim_d = std::size_t{1} << D;
    // dst 局部位 p → src 局部位编号 a（非 src 为 -1）
    std::vector<int> dst_to_src(D, -1);
    for (std::size_t a = 0; a < S; ++a) {
        for (std::size_t p = 0; p < D; ++p) {
            if (dst[p] == src[a]) { dst_to_src[p] = static_cast<int>(a); break; }
        }
    }
    std::vector<complex_t> out(dim_d * dim_d, complex_t{0, 0});
    for (std::size_t row = 0; row < dim_d; ++row) {
        for (std::size_t col = 0; col < dim_d; ++col) {
            std::size_t r_s = 0, c_s = 0;
            bool match = true;
            for (std::size_t p = 0; p < D; ++p) {
                const std::size_t br = (row >> p) & 1, bc = (col >> p) & 1;
                if (dst_to_src[p] < 0) {
                    if (br != bc) { match = false; break; }
                } else {
                    const auto a = static_cast<std::size_t>(dst_to_src[p]);
                    r_s |= br << a;
                    c_s |= bc << a;
                }
            }
            if (match) out[row * dim_d + col] = m[r_s * dim_s + c_s];
        }
    }
    return out;
}

struct Block {
    std::vector<qubit_t>     qubits;  // 升序
    std::vector<complex_t>   mat;     // 2^k × 2^k
    std::vector<Instruction> gates;   // 原始门（flush 决策与单门还原用）
};

// 一个幺正门是否可参与融合
bool fusable(const Gate& g, std::size_t max_block) {
    return g.is_unitary() && g.type != GateType::Permutation &&
           g.num_qubits() <= max_block;
}

std::vector<qubit_t> sorted_union(const std::vector<qubit_t>& a,
                                  const std::vector<qubit_t>& b) {
    std::vector<qubit_t> u(a);
    for (qubit_t q : b) {
        if (std::find(u.begin(), u.end(), q) == u.end()) u.push_back(q);
    }
    std::sort(u.begin(), u.end());
    return u;
}

}  // namespace

QuantumCircuit GateFusion::run(const QuantumCircuit& circuit) const {
    if (max_block_size_ <= 1) {
        // 退化为单比特融合（保留其恒等消去优化）
        return SingleQubitGateFusion().run(circuit);
    }

    QuantumCircuit out(circuit.num_qubits(), circuit.num_clbits());
    std::vector<Block> blocks;                     // 开放块
    std::vector<int> owner(circuit.num_qubits(), -1);  // qubit → 块索引

    auto flush_block = [&](int bi) {
        Block& b = blocks[bi];
        if (b.gates.size() == 1) {
            out.append(b.gates[0]);  // 单门保留，走快速内核
        } else {
            out.append({Gate{GateType::Unitary, {}, 0, b.mat}, b.qubits, {}});
        }
        for (qubit_t q : b.qubits) owner[q] = -1;
        b.gates.clear();  // 标记失效
    };
    auto flush_overlapping = [&](const std::vector<qubit_t>& qs) {
        std::vector<int> hit;
        for (qubit_t q : qs) {
            if (owner[q] >= 0 &&
                std::find(hit.begin(), hit.end(), owner[q]) == hit.end()) {
                hit.push_back(owner[q]);
            }
        }
        for (int bi : hit) flush_block(bi);
    };
    auto flush_all = [&] {
        for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
            if (!blocks[bi].gates.empty()) flush_block(static_cast<int>(bi));
        }
    };

    for (const auto& inst : circuit.instructions()) {
        const Gate& g = inst.gate;
        if (g.type == GateType::Barrier) {
            flush_all();
            out.append(inst);
            continue;
        }
        if (!fusable(g, max_block_size_)) {
            // Measure/Reset：flush 全部开放块，避免把延迟输出的融合块排到
            // 测量之后（那会破坏「门在前、测量在后」结构，触发逐 shot 路径）。
            // 大幺正门/置换门只与重叠块冲突（彼此可对易），flush 重叠即可。
            if (g.type == GateType::Measure || g.type == GateType::Reset) {
                flush_all();
            } else {
                flush_overlapping(inst.qubits);
            }
            out.append(inst);
            continue;
        }

        // 收集与新门重叠的开放块
        std::vector<int> involved;
        for (qubit_t q : inst.qubits) {
            if (owner[q] >= 0 &&
                std::find(involved.begin(), involved.end(), owner[q]) ==
                    involved.end()) {
                involved.push_back(owner[q]);
            }
        }
        std::vector<qubit_t> new_qubits = inst.qubits;
        std::sort(new_qubits.begin(), new_qubits.end());
        for (int bi : involved) new_qubits = sorted_union(new_qubits, blocks[bi].qubits);

        if (new_qubits.size() > max_block_size_) {
            // 超出块上限：先 flush 重叠块，新门自成一块
            for (int bi : involved) flush_block(bi);
            involved.clear();
            new_qubits = inst.qubits;
            std::sort(new_qubits.begin(), new_qubits.end());
        }

        // 合并 involved 块 + 新门 → 新块（块间不相交可对易，门在块之后左乘）
        const std::size_t dim = std::size_t{1} << new_qubits.size();
        std::vector<complex_t> acc(dim * dim, complex_t{0, 0});
        for (std::size_t i = 0; i < dim; ++i) acc[i * dim + i] = 1;  // 单位阵
        Block merged;
        merged.qubits = new_qubits;
        for (int bi : involved) {
            acc = matmul(expand(blocks[bi].mat, blocks[bi].qubits, new_qubits),
                         acc, dim);
            for (auto& gate : blocks[bi].gates) merged.gates.push_back(std::move(gate));
            blocks[bi].gates.clear();
            for (qubit_t q : blocks[bi].qubits) owner[q] = -1;
        }
        acc = matmul(expand(g.matrix(), inst.qubits, new_qubits), acc, dim);
        merged.gates.push_back(inst);
        merged.mat = std::move(acc);

        // 放入一个空槽或追加，并更新 owner
        int slot = -1;
        for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
            if (blocks[bi].gates.empty()) { slot = static_cast<int>(bi); break; }
        }
        if (slot < 0) {
            slot = static_cast<int>(blocks.size());
            blocks.push_back({});
        }
        blocks[slot] = std::move(merged);
        for (qubit_t q : blocks[slot].qubits) owner[q] = slot;
    }
    flush_all();
    return out;
}

}  // namespace lqcs
