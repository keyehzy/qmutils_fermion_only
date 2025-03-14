#include "qmutils/matrix_elements.h"

#include <gtest/gtest.h>

#include <complex>
#include <vector>

#include "qmutils/sparse_matrix.h"

namespace qmutils {
namespace {

class OperatorMatrixElementsTest : public ::testing::Test {
 protected:
  static Operator c(Operator::Spin spin, uint8_t orbital) {
    return Operator::creation(spin, orbital);
  }

  static Operator a(Operator::Spin spin, uint8_t orbital) {
    return Operator::annihilation(spin, orbital);
  }
};

TEST_F(OperatorMatrixElementsTest, SparseMatrixComputationOffUpperDiagonal) {
  Basis basis(2, 1);

  Expression H = Expression(
      Term(1.0f, {c(Operator::Spin::Up, 0), a(Operator::Spin::Up, 1)}));

  auto matrix = compute_matrix_elements_serial<
      SparseMatrix<Expression::coefficient_type>>(basis, H);

  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}

TEST_F(OperatorMatrixElementsTest, SparseMatrixComputationOffLowerDiagonal) {
  Basis basis(2, 1);

  Expression H = Expression(
      Term(1.0f, {c(Operator::Spin::Up, 1), a(Operator::Spin::Up, 0)}));

  auto matrix = compute_matrix_elements_serial<
      SparseMatrix<Expression::coefficient_type>>(basis, H);

  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}

TEST_F(OperatorMatrixElementsTest, SparseMatrixComputationOffDiagonal) {
  Basis basis(2, 1);

  Expression H;
  H += Expression(
      Term(1.0f, {c(Operator::Spin::Up, 1), a(Operator::Spin::Up, 0)}));
  H += Expression(
      Term(1.0f, {c(Operator::Spin::Up, 0), a(Operator::Spin::Up, 1)}));

  auto matrix = compute_matrix_elements_serial<
      SparseMatrix<Expression::coefficient_type>>(basis, H);

  EXPECT_EQ(matrix.rows(), 4);
  EXPECT_EQ(matrix.cols(), 4);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(1.0f, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
}

TEST_F(OperatorMatrixElementsTest,
       SparseMatrixComputationOffUpperDiagonalTwoBody) {
  Basis basis(2, 2);

  Expression H = Expression(
      Term(2.0f, {c(Operator::Spin::Up, 0), c(Operator::Spin::Down, 0),
                  a(Operator::Spin::Up, 0), a(Operator::Spin::Down, 1)}));

  auto matrix = compute_matrix_elements_serial<
      SparseMatrix<Expression::coefficient_type>>(basis, H);

  EXPECT_EQ(matrix.rows(), 6);
  EXPECT_EQ(matrix.cols(), 6);
  EXPECT_EQ(matrix(0, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 1), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(0, 2), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 1), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(1, 2), Expression::coefficient_type(-2.0f, 0));
  EXPECT_EQ(matrix(2, 0), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(2, 1), Expression::coefficient_type(0, 0));
  EXPECT_EQ(matrix(2, 2), Expression::coefficient_type(0, 0));
}

}  // namespace
}  // namespace qmutils
