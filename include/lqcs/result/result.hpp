#pragma once

#include <cstddef>
#include <iostream>
#include <map>
#include <string>

namespace lqcs {

// 一次 run 的测量结果。比特串格式与 Qiskit 一致：最高位经典比特在左。
class Result {
   public:
    Result() = default;
    Result(std::map<std::string, std::size_t> counts, std::size_t shots)
        : counts_(std::move(counts)), shots_(shots) {}

    const std::map<std::string, std::size_t>& counts() const { return counts_; }
    std::size_t shots() const { return shots_; }

    std::map<std::string, double> probabilities() const;
    void print_counts(std::ostream& os = std::cout) const;

   private:
    std::map<std::string, std::size_t> counts_;
    std::size_t shots_ = 0;
};

}  // namespace lqcs
