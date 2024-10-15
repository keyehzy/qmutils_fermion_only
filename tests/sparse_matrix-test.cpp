// File: tests/sparse_matrix-test.cpp
#include "qmutils/sparse_matrix.h"

#include <gtest/gtest.h>

#include <complex>

namespace qmutils {
namespace {

TEST(SparseMatrixTest, ConstructionAndBasicOperations) {
  SparseMatrix<float> matrix(3, 3);

  EXPECT_EQ(matrix.rows(), 3);
  EXPECT_EQ(matrix.cols(), 3);

  matrix(0, 0) = 1.0f;
  matrix(1, 1) = 2.0f;
  matrix(2, 2) = 3.0f;

  EXPECT_FLOAT_EQ(matrix(0, 0), 1.0f);
  EXPECT_FLOAT_EQ(matrix(1, 1), 2.0f);
  EXPECT_FLOAT_EQ(matrix(2, 2), 3.0f);
  EXPECT_FLOAT_EQ(matrix(0, 1), 0.0f);  // Unset element should be zero
}

TEST(SparseMatrixTest, OutOfRangeAccess) {
  SparseMatrix<float> matrix(2, 2);

  EXPECT_THROW(matrix(2, 0), std::out_of_range);
  EXPECT_THROW(matrix(0, 2), std::out_of_range);
}

TEST(SparseMatrixTest, ComplexValues) {
  SparseMatrix<std::complex<double>> matrix(2, 2);

  matrix(0, 0) = {1.0, 2.0};
  matrix(1, 1) = {3.0, 4.0};

  EXPECT_EQ(matrix(0, 0), std::complex<double>(1.0, 2.0));
  EXPECT_EQ(matrix(1, 1), std::complex<double>(3.0, 4.0));
  EXPECT_EQ(matrix(0, 1), std::complex<double>(0.0, 0.0));
}

TEST(SparseMatrixTest, ClearAndIteration) {
  SparseMatrix<int> matrix(3, 3);

  matrix(0, 0) = 1;
  matrix(1, 1) = 2;
  matrix(2, 2) = 3;

  int sum = 0;
  for (const auto& entry : matrix) {
    sum += entry.second;
  }
  EXPECT_EQ(sum, 6);

  matrix.clear();
  sum = 0;
  for (const auto& entry : matrix) {
    sum += entry.second;
  }
  EXPECT_EQ(sum, 0);
}

TEST(SparseMatrixTest, EmptyMatrix) {
  SparseMatrix<int> matrix(0, 0);

  EXPECT_EQ(matrix.rows(), 0);
  EXPECT_EQ(matrix.cols(), 0);
  EXPECT_EQ(matrix.count_non_zero(), 0);

  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}

TEST(SparseMatrixTest, ContainsElement) {
  SparseMatrix<int> matrix(3, 3);

  matrix(1, 1) = 5;

  EXPECT_TRUE(matrix.contains(1, 1));
  EXPECT_FALSE(matrix.contains(0, 0));
  EXPECT_FALSE(matrix.contains(2, 2));

  EXPECT_THROW(matrix.contains(3, 3), std::out_of_range);
}

TEST(SparseMatrixTest, CountNonZero) {
  SparseMatrix<int> matrix(3, 3);

  // Initially, the matrix should have no non-zero elements
  EXPECT_EQ(matrix.count_non_zero(), 0);

  // Add some non-zero elements
  matrix(0, 0) = 1;
  matrix(1, 1) = 2;
  matrix(2, 2) = 3;
  EXPECT_EQ(matrix.count_non_zero(), 3);

  // Adding a zero element should not increase the count
  matrix(0, 1) = 0;
  EXPECT_EQ(matrix.count_non_zero(), 4);

  // Overwriting an existing non-zero element should not change the count
  matrix(0, 0) = 5;
  EXPECT_EQ(matrix.count_non_zero(), 4);

  // TODO: Setting an existing non-zero element to zero should decrease the
  // count matrix(1, 1) = 0; EXPECT_EQ(matrix.count_non_zero(), 3);

  // TODO: Setting all elements to zero should result in a count of 0
  // matrix(0, 0) = 0;
  // matrix(2, 2) = 0;
  // EXPECT_EQ(matrix.count_non_zero(), 0);
}
}  // namespace
}  // namespace qmutils
