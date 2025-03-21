#pragma once

#include <omp.h>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/normal_order.h"

namespace qmutils {

template <typename VectorType, typename Basis>
VectorType compute_vector_elements_serial(const Basis& basis,
                                          const Expression& A) {
  VectorType vector_elements(basis.size());
  NormalOrderer orderer;
  for (size_t i = 0; i < basis.size(); ++i) {
    Expression left(basis.at(i));
    Expression product = orderer.normal_order(left.adjoint() * A);
    vector_elements(i) = product[{}];
  }
  return vector_elements;
}

template <typename VectorType, typename Basis>
VectorType compute_vector_elements(const Basis& basis, const Expression& A) {
  VectorType vector_elements(basis.size());

#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < basis.size(); ++i) {
      Expression left(basis.at(i));
      Expression product = orderer.normal_order(left.adjoint() * A);
      vector_elements(i) = product[{}];
    }
  }

  return vector_elements;
}

template <typename MatrixType, typename Basis>
MatrixType compute_matrix_elements_serial(const Basis& basis,
                                          const Expression& A) {
  MatrixType matrix_elements(basis.size(), basis.size());
  NormalOrderer orderer;

  for (size_t j = 0; j < basis.size(); ++j) {
    Expression right(basis.at(j));
    Expression product = orderer.normal_order(A * right);

    std::erase_if(product.terms(), [&](const auto& item) {
      return !basis.contains(item.first);
    });

    for (const auto& term : product.terms()) {
      size_t i = static_cast<size_t>(basis.index_of(term.first));
      matrix_elements(i, j) = term.second / basis.at(i).coefficient();
    }
  }

  return matrix_elements;
}

template <typename MatrixType, typename Basis>
MatrixType compute_matrix_elements(const Basis& basis, const Expression& A) {
  MatrixType matrix_elements(basis.size(), basis.size());

#pragma omp parallel
  {
    NormalOrderer orderer;
#pragma omp for schedule(dynamic)
    for (size_t j = 0; j < basis.size(); ++j) {
      Expression right(basis.at(j));
      Expression product = orderer.normal_order(A * right);

      std::erase_if(product.terms(), [&](const auto& item) {
        return !basis.contains(item.first);
      });

      std::vector<std::pair<size_t, Expression::coefficient_type>>
          local_updates;
      local_updates.reserve(product.terms().size());

      for (const auto& term : product.terms()) {
        size_t i = static_cast<size_t>(basis.index_of(term.first));
        local_updates.emplace_back(i, term.second / basis.at(i).coefficient());
      }

#pragma omp critical
      {
        for (const auto& [i, val] : local_updates) {
          matrix_elements(i, j) = val;
        }
      }
    }
  }

  return matrix_elements;
}
}  // namespace qmutils
