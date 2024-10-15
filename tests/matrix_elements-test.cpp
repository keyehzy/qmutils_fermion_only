#include "qmutils/matrix_elements.h"

#include <gtest/gtest.h>

#include <complex>
#include <vector>

#include "qmutils/sparse_matrix.h"

namespace qmutils {
namespace {

class MatrixElementsTest : public ::testing::Test {};

TEST_F(MatrixElementsTest, SparseMatrixComputationOffDiagonal1) {
  Basis basis(2, 1);  // 2 orbitals, 1 particle

  Expression H =
      Expression(Term(1.0f, {Operator::creation(Operator::Spin::Up, 0),
                             Operator::annihilation(Operator::Spin::Up, 1)}));

  SparseMatrix<Expression::coefficient_type> matrix(basis.size(), basis.size());
  compute_matrix_elements(matrix, basis, H.adjoint());

  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}

TEST_F(MatrixElementsTest, SparseMatrixComputationOffDiagonal2) {
  Basis basis(2, 1);  // 2 orbitals, 1 particle

  Expression H =
      Expression(Term(1.0f, {Operator::creation(Operator::Spin::Up, 1),
                             Operator::annihilation(Operator::Spin::Up, 0)}));

  SparseMatrix<Expression::coefficient_type> matrix(basis.size(), basis.size());
  compute_matrix_elements(matrix, basis, H.adjoint());

  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}

TEST_F(MatrixElementsTest, SparseMatrixComputationOffDiagonal3) {
  Basis basis(2, 1);  // 2 orbitals, 1 particle

  Expression H;
  H += Expression(Term(1.0f, {Operator::creation(Operator::Spin::Up, 1),
                              Operator::annihilation(Operator::Spin::Up, 0)}));
  H += Expression(Term(1.0f, {Operator::creation(Operator::Spin::Up, 0),
                              Operator::annihilation(Operator::Spin::Up, 1)}));

  SparseMatrix<Expression::coefficient_type> matrix(basis.size(), basis.size());
  compute_matrix_elements(matrix, basis, H.adjoint());

  // FIXME: this doesn't work at the moment
  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}
}  // namespace
}  // namespace qmutils
