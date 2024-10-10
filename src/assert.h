#pragma once

#include <cstdlib>
#include <iostream>
#include <source_location>

namespace qmutils {

[[noreturn]] inline void assert_fail(
    const char* expression,
    const std::source_location& location = std::source_location::current()) {
  std::cerr << location.file_name() << ":" << location.line() << ": "
            << location.function_name() << ": Assertion `" << expression
            << "' failed." << std::endl;
  __builtin_trap();
}

#ifdef NDEBUG
#define QMUTILS_ASSERT(condition) ((void)0)
#else
#define QMUTILS_ASSERT(condition) \
  ((condition)                    \
       ? ((void)0)                \
       : ::qmutils::assert_fail(#condition, std::source_location::current()))
#endif

}  // namespace qmutils
