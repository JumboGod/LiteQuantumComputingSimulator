cmake_minimum_required(VERSION 3.20)
project(LiteQuantumComputingSimulator VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(LQCS_BUILD_TESTS "Build unit tests" ON)
option(LQCS_BUILD_EXAMPLES "Build examples" ON)

add_library(lqcs STATIC
    src/circuit/gate.cpp
    src/circuit/quantum_circuit.cpp
    src/backend/statevector.cpp
    src/backend/statevector_simulator.cpp
    src/backend/kernels/gate_kernels.cpp
    src/result/result.cpp
    src/visualization/circuit_drawer.cpp
)
target_include_directories(lqcs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_include_directories(lqcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(lqcs PRIVATE -Wall -Wextra)

if(LQCS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(LQCS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
