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
option(LQCS_BUILD_PYTHON "Build Python bindings (pylqcs)" OFF)
option(LQCS_BUILD_BENCHMARKS "Build benchmarks" ON)
option(LQCS_ENABLE_OPENMP "Enable OpenMP parallel kernels" ON)

add_library(lqcs STATIC
    src/algorithms/arithmetic.cpp
    src/algorithms/bernstein_vazirani.cpp
    src/algorithms/deutsch_jozsa.cpp
    src/algorithms/grover.cpp
    src/algorithms/number_theory.cpp
    src/algorithms/oracles.cpp
    src/algorithms/qft.cpp
    src/algorithms/qpe.cpp
    src/algorithms/shor.cpp
    src/algorithms/simon.cpp
    src/algorithms/vqe.cpp
    src/circuit/gate.cpp
    src/circuit/quantum_circuit.cpp
    src/circuit/parametric_circuit.cpp
    src/transpiler/gate_fusion.cpp
    src/transpiler/redundant_elimination.cpp
    src/backend/density_matrix.cpp
    src/backend/density_matrix_simulator.cpp
    src/backend/noise.cpp
    src/backend/statevector.cpp
    src/backend/statevector_simulator.cpp
    src/io/qasm.cpp
    src/backend/kernels/gate_kernels.cpp
    src/result/result.cpp
    src/visualization/circuit_drawer.cpp
)
target_include_directories(lqcs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_include_directories(lqcs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
if(MSVC)
    target_compile_options(lqcs PRIVATE /W4 /permissive-)
else()
    target_compile_options(lqcs PRIVATE -Wall -Wextra)
endif()
# 静态库可能被链接进 Python 扩展（共享库）
set_target_properties(lqcs PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(LQCS_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(lqcs PUBLIC OpenMP::OpenMP_CXX)
        if(MSVC)
            # MSVC 默认 OpenMP 2.0 不接受无符号循环变量；LLVM 模式支持
            target_compile_options(lqcs PRIVATE /openmp:llvm)
        endif()
    else()
        message(STATUS "OpenMP not found - building serial kernels")
    endif()
endif()

if(LQCS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(LQCS_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(LQCS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# python 子目录的 add_test 依赖 enable_testing() 已生效，必须放在其后
if(LQCS_BUILD_PYTHON)
    add_subdirectory(python)
endif()
