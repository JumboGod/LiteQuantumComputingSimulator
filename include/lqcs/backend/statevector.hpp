#pragma once

#include <cstddef>
#include <vector>

#include "../core/types.hpp"

namespace lqcs {

// 状态向量容器：长度 2^n 的复振幅数组，初始为 |0...0>
class Statevector {
public:
    static constexpr std::size_t kMaxQubits = 30;  // 30 比特 ≈ 16 GB

    explicit Statevector(std::size_t num_qubits);

    std::size_t num_qubits() const { return num_qubits_; }
    std::size_t size() const { return data_.size(); }

    complex_t&       operator[](std::size_t i) { return data_[i]; }
    const complex_t& operator[](std::size_t i) const { return data_[i]; }

    complex_t*       data() { return data_.data(); }
    const complex_t* data() const { return data_.data(); }

    double norm() const;
    std::vector<double> probabilities() const;

private:
    std::size_t               num_qubits_;
    aligned_vector<complex_t> data_;
};

}  // namespace lqcs
