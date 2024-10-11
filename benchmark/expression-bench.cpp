#include <benchmark/benchmark.h>

#include <random>

#include "qmutils/expression.h"

namespace qmutils {
namespace {

std::complex<float> randomComplex() {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<float> dis(-1.0, 1.0);
  return {dis(gen), dis(gen)};
}

Operator randomOperator() {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_int_distribution<> type_dis(0, 1);
  static std::uniform_int_distribution<> spin_dis(0, 1);
  static std::uniform_int_distribution<> orbital_dis(0, 63);

  return Operator(static_cast<Operator::Type>(type_dis(gen)),
                  static_cast<Operator::Spin>(spin_dis(gen)), orbital_dis(gen));
}

Expression generateLargeExpression(size_t numTerms, size_t opsPerTerm) {
  Expression expr;
  for (size_t i = 0; i < numTerms; ++i) {
    Term::container_type ops;
    for (size_t j = 0; j < opsPerTerm; ++j) {
      ops.push_back(randomOperator());
    }
    expr += Term(randomComplex(), std::move(ops));
  }
  return expr;
}

static void BM_LargeExpressionConstruction(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  for (auto _ : state) {
    Expression expr = generateLargeExpression(numTerms, opsPerTerm);
    benchmark::DoNotOptimize(expr);
  }

  state.SetComplexityN(numTerms * opsPerTerm);
}

static void BM_LargeExpressionAddition(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  Expression expr1 = generateLargeExpression(numTerms, opsPerTerm);
  Expression expr2 = generateLargeExpression(numTerms, opsPerTerm);

  for (auto _ : state) {
    Expression result = expr1 + expr2;
    benchmark::DoNotOptimize(result);
  }

  state.SetComplexityN(numTerms * opsPerTerm);
}

static void BM_LargeExpressionMultiplication(benchmark::State& state) {
  const size_t numTerms = state.range(0);
  const size_t opsPerTerm = state.range(1);

  Expression expr1 = generateLargeExpression(numTerms, opsPerTerm);
  Expression expr2 = generateLargeExpression(numTerms, opsPerTerm);

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
