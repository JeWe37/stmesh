#ifndef SAMPLE_LIBRARY_HPP
#define SAMPLE_LIBRARY_HPP

#include <stmesh/stmesh_export.hpp>

[[nodiscard]] STMESH_EXPORT int factorial(int) noexcept;

[[nodiscard]] constexpr int factorial_constexpr(int input) noexcept {
  if (input == 0) {
    return 1;
  }

  return input * factorial_constexpr(input - 1);
}

#endif
