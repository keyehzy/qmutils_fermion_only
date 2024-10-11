#include "qmutils/expression.h"

#include <gtest/gtest.h>

namespace qmutils {
namespace {

class ExpressionTest : public ::testing::Test {
 protected:
  Operator op1 = Operator::creation(Operator::Spin::Up, 0);
  Operator op2 = Operator::annihilation(Operator::Spin::Down, 1);
  std::complex<float> coeff{0.5f, -0.5f};
};

TEST_F(ExpressionTest, DefaultConstructor) {
  Expression expr;
  EXPECT_EQ(expr.size(), 0);
}

}  // namespace
}  // namespace qmutils
