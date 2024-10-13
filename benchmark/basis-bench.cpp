#include <benchmark/benchmark.h>

#include <algorithm>

#include "qmutils/basis.h"

namespace qmutils {
namespace {

static void BM_BasisConstruction(benchmark::State& state) {
  const size_t orbitals = state.range(0);
  const size_t particles =
      std::min(static_cast<size_t>(state.range(1)), 2 * orbitals);

  for (auto _ : state) {
    Basis basis(orbitals, particles);
    benchmark::DoNotOptimize(basis);
  }

  state.SetComplexityN(orbitals * particles);
}

static void BM_BasisContains(benchmark::State& state) {
  const size_t orbitals = state.range(0);
  const size_t particles =
      std::min(static_cast<size_t>(state.range(1)), 2 * orbitals);

  Basis basis(orbitals, particles);

  // Create a state that's guaranteed to be in the basis
  Basis::operators_type test_state;
  for (size_t i = 0; i < particles; ++i) {
    test_state.push_back(Operator::creation(Operator::Spin::Up, i % orbitals));
  }

  for (auto _ : state) {
    bool result = basis.contains(test_state);
    benchmark::DoNotOptimize(result);
  }

  state.SetComplexityN(orbitals * particles);
}

static void BM_BasisIndex(benchmark::State& state) {
  const size_t orbitals = state.range(0);
  const size_t particles =
      std::min(static_cast<size_t>(state.range(1)), 2 * orbitals);

  Basis basis(orbitals, particles);

  // Create a state that's guaranteed to be in the basis
  Basis::operators_type test_state;
  for (size_t i = 0; i < particles; ++i) {
    Operator::Spin spin =
        (i % 2 == 0) ? Operator::Spin::Up : Operator::Spin::Down;
    test_state.push_back(Operator::creation(spin, i / 2));
  }

  for (auto _ : state) {
    size_t result = basis.index(test_state);
    benchmark::DoNotOptimize(result);
  }

  state.SetComplexityN(orbitals * particles);
}

BENCHMARK(BM_BasisConstruction)
    ->RangeMultiplier(2)
    ->Ranges({{1, 8}, {1, 16}})
    ->Complexity();

BENCHMARK(BM_BasisContains)
    ->RangeMultiplier(2)
    ->Ranges({{1, 8}, {1, 16}})
    ->Complexity();

BENCHMARK(BM_BasisIndex)
    ->RangeMultiplier(2)
    ->Ranges({{1, 8}, {1, 16}})
    ->Complexity();

}  // namespace
}  // namespace qmutils

BENCHMARK_MAIN();
