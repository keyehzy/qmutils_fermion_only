#pragma once

#include <vector>

#include "qmutils/assert.h"

namespace qmutils {

template <typename VecType, typename MatType>
class EigenSystem {
 public:
  EigenSystem(VecType eigenvalues, MatType eigenvectors)
      : eigenvalues_(eigenvalues), eigenvectors_(eigenvectors) {}

  const VecType& eigenvalues() const { return eigenvalues_; }

  const MatType& eigenvectors() const { return eigenvectors_; }

  auto eigenvalue(size_t index) const {
    QMUTILS_ASSERT(index < size());
    return eigenvalues_[index];
  }

  auto eigenvector(size_t index) const {
    QMUTILS_ASSERT(index < size());
    return eigenvectors_.col(index);
  }

  size_t size() const { return eigenvalues_.size(); }

 private:
  VecType eigenvalues_;
  MatType eigenvectors_;
};

}  // namespace qmutils
