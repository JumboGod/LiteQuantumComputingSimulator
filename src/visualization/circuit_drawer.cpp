// QuantumCircuit::draw() 的实现：纯 ASCII 电路图。
// 布局：每条指令占一列；偶数行为量子比特线（'-'），其间为连接线行，
// 最后一行为经典寄存器线（'='）。控制位 '*'，SWAP 目标 'x'，测量 'M'。
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "lqcs/circuit/quantum_circuit.hpp"

namespace lqcs {

namespace {

std::string format_params(const std::vector<double>& params) {
    if (params.empty()) return "";
    std::string s = "(";
    char buf[32];
    for (std::size_t i = 0; i < params.size(); ++i) {
        std::snprintf(buf, sizeof(buf), "%.4g", params[i]);
        if (i > 0) s += ",";
        s += buf;
    }
    return s + ")";
}

std::string center(const std::string& text, std::size_t width, char fill) {
    const std::size_t pad = width - text.size();
    const std::size_t left = pad / 2;
    return std::string(left, fill) + text + std::string(pad - left, fill);
}

}  // namespace

std::string QuantumCircuit::draw() const {
    const std::size_t nq = num_qubits();
    const bool has_creg = num_clbits() > 0;
    // 行布局：量子比特 q 在第 2q 行，其后是连接线行；经典线在最后
    const std::size_t n_rows = 2 * nq + (has_creg ? 1 : 0);
    const std::size_t creg_row = 2 * nq;

    auto fill_char = [&](std::size_t row) -> char {
        if (has_creg && row == creg_row) return '=';
        return (row % 2 == 0 && row < 2 * nq) ? '-' : ' ';
    };

    // 行前缀
    std::vector<std::string> prefix(n_rows);
    for (std::size_t q = 0; q < nq; ++q) {
        prefix[2 * q] = "q_" + std::to_string(q) + ": ";
    }
    if (has_creg) {
        prefix[creg_row] = "c: " + std::to_string(num_clbits()) + "/";
    }
    std::size_t pw = 0;
    for (const auto& p : prefix) pw = std::max(pw, p.size());
    for (std::size_t r = 0; r < n_rows; ++r) {
        prefix[r] = std::string(pw - prefix[r].size(), ' ') + prefix[r];
    }

    std::vector<std::string> lines(n_rows);
    for (std::size_t r = 0; r < n_rows; ++r) {
        lines[r] = prefix[r] + fill_char(r);
    }

    for (const auto& inst : instructions()) {
        const Gate& g = inst.gate;
        std::map<std::size_t, std::string> cells;  // 行 → 单元格文本
        std::size_t span_lo = n_rows, span_hi = 0;

        if (g.type == GateType::Barrier) {
            for (std::size_t r = 0; r < 2 * nq; ++r) cells[r] = ":";
        } else if (g.type == GateType::Measure) {
            const std::size_t qrow = 2 * inst.qubits[0];
            cells[qrow] = "M";
            cells[creg_row] = std::to_string(inst.clbits[0]);
            span_lo = qrow;
            span_hi = creg_row;
        } else {
            // 控制位在前、目标位在后
            for (std::size_t j = 0; j < g.n_controls; ++j) {
                cells[2 * inst.qubits[j]] = "*";
            }
            const bool is_swap = (g.type == GateType::SWAP);
            const std::string label =
                is_swap ? "x"
                        : "[" + Gate{g.type, g.params, 0, g.mat, g.perm}.name() +
                              format_params(g.params) + "]";
            for (std::size_t j = g.n_controls; j < inst.qubits.size(); ++j) {
                cells[2 * inst.qubits[j]] = label;
            }
            for (qubit_t q : inst.qubits) {
                span_lo = std::min<std::size_t>(span_lo, 2 * q);
                span_hi = std::max<std::size_t>(span_hi, 2 * q);
            }
        }

        std::size_t width = 1;
        for (const auto& [row, text] : cells) {
            width = std::max(width, text.size());
        }

        for (std::size_t r = 0; r < n_rows; ++r) {
            const char fc = fill_char(r);
            std::string cell;
            if (auto it = cells.find(r); it != cells.end()) {
                cell = center(it->second, width, fc);
            } else if (r > span_lo && r < span_hi) {
                cell = center("|", width, fc);  // 跨行连接线
            } else {
                cell = std::string(width, fc);
            }
            lines[r] += cell;
            lines[r] += fc;  // 列间隔
        }
    }

    std::string out;
    for (std::size_t r = 0; r < n_rows; ++r) {
        // 连接线行全空白时跳过，让图更紧凑
        if (r % 2 == 1 && lines[r].find_first_not_of(' ') == std::string::npos) {
            continue;
        }
        out += lines[r];
        out += '\n';
    }
    return out;
}

}  // namespace lqcs
