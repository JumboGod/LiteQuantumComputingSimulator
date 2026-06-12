# Python 绑定模块（pybind11 由 FetchContent 引入）
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
set(PYBIND11_FINDPYTHON ON)

include(FetchContent)
FetchContent_Declare(
    pybind11
    URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.6.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_pylqcs src/bindings.cpp)
target_link_libraries(_pylqcs PRIVATE lqcs)

# 开发模式：模块直接输出到 python/pylqcs/，PYTHONPATH=python 即可 import
set_target_properties(_pylqcs PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pylqcs)

# pip / scikit-build-core 安装路径
if(SKBUILD)
    install(TARGETS _pylqcs LIBRARY DESTINATION pylqcs)
endif()

# pytest 集成进 ctest
if(LQCS_BUILD_TESTS)
    add_test(NAME python_bindings
             COMMAND ${Python3_EXECUTABLE} -m pytest
                     ${CMAKE_SOURCE_DIR}/tests/python -q)
    set_tests_properties(python_bindings PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/python")
endif()
