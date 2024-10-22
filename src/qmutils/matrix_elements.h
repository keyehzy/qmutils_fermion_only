#pragma once

#include <omp.h>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/normal_order.h"

namespace qmutils {

template <typename MatrixType>
MatrixType compute_matrix_elements_serial(const Basis& basis,
                                          const Expression& A) {
  MatrixType matrix_elements(basis.size(), basis.size());
  NormalOrderer orderer;
  for (size_t i = 0; i < basis.size(); ++i) {
    for (size_t j = 0; j < basis.size(); ++j) {
      Expression left(Term(basis.at(i)));
      Expression right(Term(basis.at(j)));
      Expression result = orderer.normal_order(left.adjoint() * A * right);
      matrix_elements(i, j) = result[{}];
    }
  }
  return matrix_elements;
}

template <typename MatrixType>
MatrixType compute_matrix_elements(const Basis& basis, const Expression& A) {
  MatrixType matrix_elements(basis.size(), basis.size());

#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for collapse(2) schedule(dynamic)
    for (size_t i = 0; i < basis.size(); ++i) {
      for (size_t j = 0; j < basis.size(); ++j) {
        Expression left(Term(basis.at(i)));
        Expression right(Term(basis.at(j)));
        Expression result = orderer.normal_order(left.adjoint() * A * right);
        matrix_elements(i, j) = result[{}];
      }
    }
  }

  return matrix_elements;
}

template <typename MatrixType>
void compute_matrix_elements_parallel_subset(
    MatrixType& matrix_elements, const Basis& basis, const Expression& A,
    const std::vector<std::pair<size_t, size_t>>& indices) {
#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for schedule(dynamic)
    for (size_t idx = 0; idx < indices.size(); ++idx) {
      const auto& [i, j] = indices[idx];
      Expression left(Term(basis.at(i)));
      Expression right(Term(basis.at(j)));
      Expression result = orderer.normal_order(left.adjoint() * A * right);
      matrix_elements(i, j) = result[{}];
    }
  }
}

inline size_t get_optimal_chunk_size(size_t total_size,
                                     size_t min_chunk_size = 16) {
  int num_threads = omp_get_max_threads();
  size_t chunk_size = total_size / (4 * static_cast<size_t>(num_threads));
  return std::max(chunk_size, min_chunk_size);
}

template <typename MatrixType>
MatrixType compute_matrix_elements_chunked(const Basis& basis,
                                           const Expression& A,
                                           size_t chunk_size = 0) {
  MatrixType matrix_elements(basis.size(), basis.size());

  if (chunk_size == 0) {
    chunk_size = get_optimal_chunk_size(basis.size() * basis.size());
  }

#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for collapse(2) schedule(dynamic, chunk_size)
    for (size_t i = 0; i < basis.size(); ++i) {
      for (size_t j = 0; j < basis.size(); ++j) {
        Expression left(Term(basis.at(i)));
        Expression right(Term(basis.at(j)));
        Expression result = orderer.normal_order(left.adjoint() * A * right);
        matrix_elements(i, j) = result[{}];
      }
    }
  }

  return matrix_elements;
}

}  // namespace qmutils
