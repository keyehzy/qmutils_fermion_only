cmake_minimum_required(VERSION 3.14)

# find_package(OpenMP REQUIRED)

add_executable(
  main
  main.cpp
)

target_link_libraries(
  main
  qmutils
)

target_compile_options(
  main
  PRIVATE
  "${QMUTILS_CXX_COMPILER_OPTIONS}"
)

add_library(
  qmutils
  operator.cpp
  term.cpp
  expression.cpp
  normal_order.cpp
)

target_include_directories(
  qmutils
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# target_link_libraries(
#   qmutils
#   PUBLIC
#   OpenMP::OpenMP_CXX
# )
