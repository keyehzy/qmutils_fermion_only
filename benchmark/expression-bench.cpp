#include <benchmark/benchmark.h>

#include <random>

#include "qmutils/expression.h"
#include "qmutils/utils.h"

namespace qmutils {
namespace {

static void BM_LargeExpressionConstruction(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  for (auto _ : state) {
    Expression expr = generate_large_expression(numTerms, opsPerTerm);
    benchmark::DoNotOptimize(expr);
  }

  state.SetComplexityN(numTerms * opsPerTerm);
}

static void BM_LargeExpressionAddition(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  Expression expr1 = generate_large_expression(numTerms, opsPerTerm);
  Expression expr2 = generate_large_expression(numTerms, opsPerTerm);

  for (auto _ : state) {
    Expression result = expr1 + expr2;
    benchmark::DoNotOptimize(result);
  }

  state.SetComplexityN(numTerms * opsPerTerm);
}

static void BM_LargeExpressionMultiplication(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  Expression expr1 = generate_large_expression(numTerms, opsPerTerm);
  Expression expr2 = generate_large_expression(numTerms, opsPerTerm);

  for (auto _ : state) {
    Expression result = expr1 * expr2;
    benchmark::DoNotOptimize(result);
  }

  state.SetComplexityN(numTerms * opsPerTerm);
}

BENCHMARK(BM_LargeExpressionConstruction)
    ->RangeMultiplier(2)
    ->Ranges({{1 << 8, 1 << 12}, {2, 8}})
    ->Complexity();

BENCHMARK(BM_LargeExpressionAddition)
    ->RangeMultiplier(2)
    ->Ranges({{1 << 8, 1 << 12}, {2, 8}})
    ->Complexity();

BENCHMARK(BM_LargeExpressionMultiplication)
    ->RangeMultiplier(2)
    ->Ranges({{1 << 6, 1 << 10}, {2, 8}})
    ->Complexity();

}  // namespace
}  // namespace qmutils
