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

  for (size_t j = 0; j < basis.size(); ++j) {
    Expression right(Term(basis.at(j)));
    Expression product = orderer.normal_order(A * right);

    std::erase_if(product.terms(), [&](const auto& item) {
      return !basis.contains(item.first);
    });

    for (const auto& term : product.terms()) {
      size_t i = static_cast<size_t>(basis.index_of(term.first));
      matrix_elements(i, j) = term.second;
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
#pragma omp for schedule(dynamic)
    for (size_t j = 0; j < basis.size(); ++j) {
      Expression right(Term(basis.at(j)));
      Expression product = orderer.normal_order(A * right);

      std::erase_if(product.terms(), [&](const auto& item) {
        return !basis.contains(item.first);
      });

      for (const auto& term : product.terms()) {
        size_t i = static_cast<size_t>(basis.index_of(term.first));
        matrix_elements(i, j) = term.second;
      }
    }
  }

  return matrix_elements;
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
    chunk_size = get_optimal_chunk_size(basis.size());
  }

#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for schedule(dynamic, chunk_size)
    for (size_t j = 0; j < basis.size(); ++j) {
      Expression right(Term(basis.at(j)));
      Expression product = orderer.normal_order(A * right);

      std::erase_if(product.terms(), [&](const auto& item) {
        return !basis.contains(item.first);
      });

      for (const auto& term : product.terms()) {
        size_t i = static_cast<size_t>(basis.index_of(term.first));
        matrix_elements(i, j) = term.second;
      }
    }
  }

  return matrix_elements;
}

}  // namespace qmutils
