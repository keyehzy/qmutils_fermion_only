// In a new header file, e.g., qmutils/matrix_elements.h

#pragma once

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/normal_order.h"

namespace qmutils {

template <typename MatrixType>
MatrixType compute_matrix_elements(const Basis& basis, const Expression& A) {
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

}  // namespace qmutils
