cmake_minimum_required(VERSION 3.14)

project(qmutils VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(QMUTILS_CXX_COMPILER_OPTIONS "" CACHE STRING "")
mark_as_advanced(QMUTILS_CXX_COMPILER_OPTIONS)

add_subdirectory(src)
add_subdirectory(vendor)

if (BUILD_TESTING)
  include(CTest)
  add_subdirectory(tests)
endif()

if (BUILD_BENCHMARK)
  add_subdirectory(benchmark)
endif()

if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
