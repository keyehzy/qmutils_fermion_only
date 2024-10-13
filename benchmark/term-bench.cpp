#include <benchmark/benchmark.h>

#include <random>

#include "qmutils/term.h"
#include "qmutils/utils.h"

namespace qmutils {
namespace {

static void BM_ConstructLargeTerm(benchmark::State& state) {
  const size_t largeSize = state.range(0);
  auto largeOps = generate_operator_sequence(largeSize);

  for (auto _ : state) {
    Term largeTerm(random_complex(), largeOps);
    benchmark::DoNotOptimize(largeTerm);
  }
  state.SetItemsProcessed(state.iterations() * largeSize);
}
BENCHMARK(BM_ConstructLargeTerm)->Range(1 << 10, 1 << 20);

static void BM_MultiplyLargeTerms(benchmark::State& state) {
  const size_t size = state.range(0);
  Term term1(random_complex(), generate_operator_sequence(size));
  Term term2(random_complex(), generate_operator_sequence(size));

  for (auto _ : state) {
    Term result = term1 * term2;
    benchmark::DoNotOptimize(result);
  }
  state.SetItemsProcessed(state.iterations() * size * 2);
}
BENCHMARK(BM_MultiplyLargeTerms)->Range(1 << 10, 1 << 15);

static void BM_HashLargeTerm(benchmark::State& state) {
  const size_t largeSize = state.range(0);
  Term largeTerm(random_complex(), generate_operator_sequence(largeSize));
  std::hash<Term> hasher;

  for (auto _ : state) {
    size_t hash = hasher(largeTerm);
    benchmark::DoNotOptimize(hash);
  }
  state.SetItemsProcessed(state.iterations() * largeSize);
}
BENCHMARK(BM_HashLargeTerm)->Range(1 << 10, 1 << 17);

}  // namespace
}  // namespace qmutils
