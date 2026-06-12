// Shor 算法端到端：分解 N = 15 与 N = 21
#include <algorithm>
#include <iostream>
#include <vector>

#include <lqcs/lqcs.hpp>

using namespace lqcs::algorithms;

namespace {

void factor(std::uint64_t N, std::uint64_t seed) {
    const auto res = shor(N, {.seed = seed});
    if (!res.success) {
        std::cout << "N = " << N << ": failed after " << res.attempts
                  << " attempts\n";
        return;
    }
    std::cout << "N = " << N << " = " << res.factor1 << " x " << res.factor2
              << "   (a = " << res.a << ", r = " << res.period
              << ", attempts = " << res.attempts << ")\n";
    if (!res.counts.empty()) {
        // 打印计数寄存器分布的前几个峰：应聚在 s·2^t/r 附近
        std::vector<std::pair<std::size_t, std::string>> top;
        for (const auto& [key, count] : res.counts) top.emplace_back(count, key);
        std::sort(top.rbegin(), top.rend());
        std::cout << "  top measurement peaks:";
        for (std::size_t i = 0; i < std::min<std::size_t>(4, top.size()); ++i) {
            std::cout << "  " << std::stoull(top[i].second, nullptr, 2) << " ("
                      << top[i].first << ")";
        }
        std::cout << "\n";
    }
}

}  // namespace

int main() {
    std::cout << "Shor's algorithm\n================\n";
    factor(15, 42);
    factor(21, 42);
    factor(35, 7);
    return 0;
}
