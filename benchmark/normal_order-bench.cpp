#include <benchmark/benchmark.h>

#include <algorithm>
#include <random>

#include "qmutils/normal_order.h"

namespace qmutils {
namespace {

// Helper function to create a random operator
Operator random_operator(std::mt19937& gen) {
  std::uniform_int_distribution<> type_dist(0, 1);
  std::uniform_int_distribution<> spin_dist(0, 1);
  std::uniform_int_distribution<> orbital_dist(0,
                                               63);  // Assuming 6-bit orbital

  Operator::Type type = static_cast<Operator::Type>(type_dist(gen));
  Operator::Spin spin = static_cast<Operator::Spin>(spin_dist(gen));
  uint8_t orbital = static_cast<uint8_t>(orbital_dist(gen));

  return Operator(type, spin, orbital);
}

// Helper function to create a term with random operators
Term create_random_term(size_t n_operators) {
  std::random_device rd;
  std::mt19937 gen(rd());

  Term::container_type operators;
  operators.reserve(n_operators);

  for (size_t i = 0; i < n_operators; ++i) {
    operators.push_back(random_operator(gen));
  }

  return Term(Term::coefficient_type(1.0f, 0.0f), operators);
}

void PrintCacheStats(benchmark::State& state, const NormalOrderer& orderer) {
  (void)state;
  (void)orderer;
  return;
  // state.PauseTiming();
  // orderer.print_cache_stats();
  // state.ResumeTiming();
}

Term create_noncommuting_term(size_t n_operators) {
  Term::container_type operators;
  operators.reserve(n_operators);

  for (size_t i = 0; i < n_operators / 2; ++i) {
    operators.push_back(Operator::annihilation(Operator::Spin::Up, 0));
  }

  for (size_t i = 0; i < n_operators / 2; ++i) {
    operators.push_back(Operator::creation(Operator::Spin::Up, 0));
  }

  return Term(Term::coefficient_type(1.0f, 0.0f), std::move(operators));
}

Term create_commuting_term(size_t n_operators) {
  Term::container_type operators;
  operators.reserve(n_operators);

  for (size_t i = 0; i < n_operators / 2; ++i) {
    operators.push_back(Operator::creation(Operator::Spin::Up, i % 64));
  }

  for (size_t i = 0; i < n_operators / 2; ++i) {
    operators.push_back(Operator::annihilation(Operator::Spin::Up, i % 64));
  }

  std::reverse(operators.begin(), operators.end());

  return Term(Term::coefficient_type(1.0f, 0.0f), std::move(operators));
}

// Benchmark for normal ordering a single random term
static void BM_NormalOrderRandomTerm(benchmark::State& state) {
  NormalOrderer orderer;
  Term term = create_random_term(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(orderer.normal_order(term));
  }

  PrintCacheStats(state, orderer);
  state.SetComplexityN(state.range(0));
}

// Benchmark for normal ordering a term with many non-commuting operators (worst
// case)
static void BM_NormalOrderNonCommutingTerm(benchmark::State& state) {
  NormalOrderer orderer;
  Term term = create_noncommuting_term(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(orderer.normal_order(term));
  }

  PrintCacheStats(state, orderer);
  state.SetComplexityN(state.range(0));
}

static void BM_NormalOrderCommutingTerm(benchmark::State& state) {
  NormalOrderer orderer;
  Term term = create_commuting_term(state.range(0));

  for (auto _ : state) {
    benchmark::DoNotOptimize(orderer.normal_order(term));
  }

  PrintCacheStats(state, orderer);
  state.SetComplexityN(state.range(0));
}

// Benchmark for normal ordering an expression with multiple terms
static void BM_NormalOrderExpression(benchmark::State& state) {
  NormalOrderer orderer;
  Expression expr;

  // Create an expression with 10 random terms
  for (int i = 0; i < 10; ++i) {
    expr += create_random_term(state.range(0));
  }

  for (auto _ : state) {
    benchmark::DoNotOptimize(orderer.normal_order(expr));
  }

  PrintCacheStats(state, orderer);
  state.SetComplexityN(state.range(0));
}

// Register the benchmarks
BENCHMARK(BM_NormalOrderRandomTerm)->Range(1, 1 << 7)->Complexity();
BENCHMARK(BM_NormalOrderNonCommutingTerm)->Range(1, 1 << 8)->Complexity();
BENCHMARK(BM_NormalOrderCommutingTerm)->Range(1, 1 << 8)->Complexity();
BENCHMARK(BM_NormalOrderExpression)->Range(1, 1 << 7)->Complexity();

}  // namespace
}  // namespace qmutils
