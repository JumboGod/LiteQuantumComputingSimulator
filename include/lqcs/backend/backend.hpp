#pragma once

#include <cstddef>
#include <string>

#include "../circuit/quantum_circuit.hpp"
#include "../result/result.hpp"

namespace lqcs {

class Backend {
   public:
    virtual ~Backend() = default;

    virtual Result run(const QuantumCircuit& circuit,
                       std::size_t shots = 1024) = 0;
    virtual std::string name() const = 0;
};

}  // namespace lqcs
