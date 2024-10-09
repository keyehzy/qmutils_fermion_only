cmake_minimum_required(VERSION 3.14)

set(BUILD_GMOCK TRUE CACHE INTERNAL "")
set(INSTALL_GTEST FALSE CACHE INTERNAL "")
set(gmock_build_tests FALSE CACHE INTERNAL "")
set(gtest_build_samples FALSE CACHE INTERNAL "")
set(gtest_build_tests FALSE CACHE INTERNAL "")
set(gtest_force_shared_crt TRUE CACHE INTERNAL "")

add_subdirectory(googletest)
