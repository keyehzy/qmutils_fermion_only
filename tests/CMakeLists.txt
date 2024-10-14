cmake_minimum_required(VERSION 3.14)

add_executable(
    qmutils-test
    operator-test.cpp
    term-test.cpp
    expression-test.cpp
    normal_order-test.cpp
    basis-test.cpp
    sparse_matrix-test.cpp
    matrix_elements-test.cpp
)

target_compile_options(
    qmutils-test
    PRIVATE
    "${QMUTILS_CXX_COMPILER_OPTIONS}"
)

target_include_directories(
    qmutils-test
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
    qmutils-test
    PRIVATE
    gmock
    gmock_main
    gtest
    qmutils
)

add_test(
    NAME qmutils-test
    COMMAND qmutils-test
)
