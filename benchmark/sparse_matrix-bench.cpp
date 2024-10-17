
#include <benchmark/benchmark.h>

#include <random>

#include "qmutils/sparse_matrix.h"
#include "qmutils/utils.h"

namespace qmutils {
namespace {

template <typename T>
void fill_sparse_matrix(SparseMatrix<T>& matrix, double density) {
  static std::mt19937 gen = get_random_generator();
  std::uniform_real_distribution<> dis(0.0, 1.0);
  std::uniform_int_distribution<> row_dis(0, matrix.rows() - 1);
  std::uniform_int_distribution<> col_dis(0, matrix.cols() - 1);

  size_t num_elements =
      static_cast<size_t>(matrix.rows() * matrix.cols() * density);
  for (size_t i = 0; i < num_elements; ++i) {
    size_t row = row_dis(gen);
    size_t col = col_dis(gen);
    matrix(row, col) = static_cast<T>(dis(gen));
  }
}

// Benchmark for element insertion
static void BM_SparseMatrixInsertion(benchmark::State& state) {
  const size_t size = state.range(0);
  const double density = 0.1;  // 1% non-zero elements

  for (auto _ : state) {
    state.PauseTiming();
    qmutils::SparseMatrix<double> matrix(size, size);
    state.ResumeTiming();

    fill_sparse_matrix(matrix, density);
  }

  state.SetComplexityN(size * size * density);
}

// Benchmark for element access (both existing and non-existing elements)
static void BM_SparseMatrixAccess(benchmark::State& state) {
  const size_t size = state.range(0);
  const double density = 0.1;  // 1% non-zero elements

  qmutils::SparseMatrix<double> matrix(size, size);
  fill_sparse_matrix(matrix, density);

  static std::mt19937 gen = get_random_generator();
  std::uniform_int_distribution<> row_dis(0, size - 1);
  std::uniform_int_distribution<> col_dis(0, size - 1);

  for (auto _ : state) {
    size_t row = row_dis(gen);
    size_t col = col_dis(gen);
    benchmark::DoNotOptimize(matrix(row, col));
  }

  state.SetComplexityN(size * size);
}

// Benchmark for iterating over non-zero elements
static void BM_SparseMatrixIteration(benchmark::State& state) {
  const size_t size = state.range(0);
  const double density = 0.1;  // 1% non-zero elements

  qmutils::SparseMatrix<double> matrix(size, size);
  fill_sparse_matrix(matrix, density);

  for (auto _ : state) {
    double sum = 0.0;
    for (const auto& elem : matrix) {
      benchmark::DoNotOptimize(sum += elem.second);
    }
  }

  state.SetComplexityN(size * size * density);
}

BENCHMARK(BM_SparseMatrixInsertion)
    ->RangeMultiplier(2)
    ->Range(1 << 7, 1 << 11)
    ->Complexity();

BENCHMARK(BM_SparseMatrixAccess)
    ->RangeMultiplier(2)
    ->Range(1 << 7, 1 << 11)
    ->Complexity();

BENCHMARK(BM_SparseMatrixIteration)
    ->RangeMultiplier(2)
    ->Range(1 << 7, 1 << 11)
    ->Complexity();

}  // namespace
}  // namespace qmutils
