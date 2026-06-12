include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# Windows 下保持运行时库一致（无害的标准设置）
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(lqcs_tests
    test_gates.cpp
    test_circuit.cpp
    test_statevector.cpp
    test_measurement.cpp
    test_equivalence.cpp
    test_qft.cpp
    test_midcircuit.cpp
    test_transpiler.cpp
    test_number_theory.cpp
    test_algorithms.cpp
    test_shor.cpp
    test_expectation.cpp
    test_vqe.cpp
    test_qasm.cpp
    test_density.cpp
    test_trajectory.cpp
)
target_link_libraries(lqcs_tests PRIVATE lqcs GTest::gtest_main)
target_include_directories(lqcs_tests PRIVATE ${CMAKE_SOURCE_DIR}/src)

include(GoogleTest)
gtest_discover_tests(lqcs_tests)
