#include "lqcs/result/result.hpp"

namespace lqcs {

std::map<std::string, double> Result::probabilities() const {
    std::map<std::string, double> probs;
    if (shots_ == 0) return probs;
    for (const auto& [key, count] : counts_) {
        probs[key] = static_cast<double>(count) / static_cast<double>(shots_);
    }
    return probs;
}

void Result::print_counts(std::ostream& os) const {
    os << "{";
    bool first = true;
    for (const auto& [key, count] : counts_) {
        if (!first) os << ", ";
        os << "\"" << key << "\": " << count;
        first = false;
    }
    os << "}\n";
}

}  // namespace lqcs
