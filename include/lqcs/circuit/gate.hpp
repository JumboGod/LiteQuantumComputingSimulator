#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "../core/types.hpp"

namespace lqcs {

// M1 门集：全部单比特门 + CX/CZ/CP/SWAP + 测量/屏障。
// 多控制门、置换门、自定义幺正等在 M2 引入。
enum class GateType {
    // 单比特无参门
    I, X, Y, Z, H, S, Sdg, T, Tdg, SX,
    // 单比特参数门
    RX, RY, RZ, P, U,
    // 双比特门
    CX, CZ, CP, SWAP,
    // 非幺正操作
    Measure, Barrier,
};

struct Gate {
    GateType            type;
    std::vector<double> params;  // RX/RY/RZ/P/CP: {theta}; U: {theta, phi, lambda}

    std::size_t num_qubits() const;
    bool is_unitary() const;  // Measure/Barrier 为 false

    // 完整门矩阵（行主序 2^k × 2^k）。受控门返回 4×4 矩阵，
    // 仅用于测试与验证；模拟内核走子空间路径，不经过该矩阵。
    std::vector<complex_t> matrix() const;

    std::string name() const;
};

}  // namespace lqcs
