#pragma once

#include <memory>
#include <vector>

#include "pass.hpp"

namespace lqcs {

// Pass 管线：按注册顺序依次应用
class Transpiler {
   public:
    Transpiler& add_pass(std::unique_ptr<Pass> pass) {
        passes_.push_back(std::move(pass));
        return *this;
    }

    QuantumCircuit run(const QuantumCircuit& circuit) const {
        QuantumCircuit result = circuit;
        for (const auto& pass : passes_) {
            result = pass->run(result);
        }
        return result;
    }

    std::size_t num_passes() const { return passes_.size(); }

   private:
    std::vector<std::unique_ptr<Pass>> passes_;
};

}  // namespace lqcs
