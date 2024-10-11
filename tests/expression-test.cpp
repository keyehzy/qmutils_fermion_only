#include "qmutils/expression.h"

#include <gtest/gtest.h>

namespace qmutils {
namespace {

class ExpressionTest : public ::testing::Test {
 protected:
  Operator op1 = Operator::creation(Operator::Spin::Up, 0);
  Operator op2 = Operator::annihilation(Operator::Spin::Down, 1);
  Operator op3 = Operator::creation(Operator::Spin::Up, 2);
  std::complex<float> coeff1{1.0f, 0.0f};
  std::complex<float> coeff2{0.0f, 2.0f};
};

TEST_F(ExpressionTest, DefaultConstructor) {
  Expression expr;
  EXPECT_EQ(expr.size(), 0);
}

TEST_F(ExpressionTest, ConstructorWithTerm) {
  Term term(coeff1, {op1, op2});
  Expression expr(term);
  EXPECT_EQ(expr.size(), 1);
  EXPECT_EQ(expr.terms().begin()->second, coeff1);
}

TEST_F(ExpressionTest, CopyConstructor) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(expr1);
  EXPECT_EQ(expr2.size(), 1);
  EXPECT_EQ(expr2.terms(), expr1.terms());
}

TEST_F(ExpressionTest, MoveConstructor) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(std::move(expr1));
  EXPECT_EQ(expr2.size(), 1);
  EXPECT_EQ(expr1.size(), 0);  // expr1 should be empty after move
}

TEST_F(ExpressionTest, AdditionOperator) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(Term(coeff2, {op2, op3}));
  Expression result = expr1 + expr2;
  EXPECT_EQ(result.size(), 2);
}

TEST_F(ExpressionTest, SubtractionOperator) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(Term(coeff1, {op1, op2}));
  Expression result = expr1 - expr2;
  EXPECT_EQ(result.size(), 0);  // Should cancel out
}

TEST_F(ExpressionTest, MultiplicationOperator) {
  Expression expr1(Term(coeff1, {op1}));
  Expression expr2(Term(coeff2, {op2}));
  Expression result = expr1 * expr2;
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second, coeff1 * coeff2);
}

TEST_F(ExpressionTest, ScalarMultiplication) {
  Expression expr(Term(coeff1, {op1, op2}));
  std::complex<float> scalar(2.0f, 1.0f);
  Expression result = expr * scalar;
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second, coeff1 * scalar);
}

TEST_F(ExpressionTest, Normalization) {
  Expression expr;
  expr += Term(coeff1, {op1, op1});  // Should be removed (a^2 = 0)
  expr += Term(std::complex<float>(0.0f, 0.0f),
               {op2, op3});  // Should be removed (zero coefficient)
  expr += Term(coeff2, {op1, op2});
  expr.normalize();
  EXPECT_EQ(expr.size(), 1);
  EXPECT_EQ(expr.terms().begin()->second, coeff2);
}

TEST_F(ExpressionTest, ToString) {
  Expression expr;
  expr += Term(coeff1, {op1, op2});
  expr += Term(coeff2, {op2, op3});
  // Not necessarily ordered
  std::string expected = "(0,2) c(↓,1)c+(↑,2) + (1,0) c+(↑,0)c(↓,1)";
  EXPECT_EQ(expr.to_string(), expected);
}

TEST_F(ExpressionTest, EqualityOperator) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(Term(coeff1, {op1, op2}));
  Expression expr3(Term(coeff2, {op1, op2}));
  EXPECT_EQ(expr1, expr2);
  EXPECT_NE(expr1, expr3);
}

TEST_F(ExpressionTest, ComplexExpression) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(Term(coeff2, {op2, op3}));
  Expression result = expr1 * expr2 + expr1 - expr2;
  EXPECT_EQ(result.size(), 2);
}

TEST_F(ExpressionTest, ZeroExpression) {
  Expression expr1(Term(coeff1, {op1, op2}));
  Expression expr2(Term(-coeff1, {op1, op2}));
  Expression result = expr1 + expr2;
  EXPECT_EQ(result.size(), 0);
}

}  // namespace
}  // namespace qmutils
