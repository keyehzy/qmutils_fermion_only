// File: tests/normal_order_test.cpp

#include "qmutils/normal_order.h"

#include <gtest/gtest.h>

namespace qmutils {
namespace {

class NormalOrderTest : public ::testing::Test {
 protected:
  NormalOrderer orderer;

  // Helper function to create operators
  static Operator c(Operator::Spin spin, uint8_t orbital) {
    return Operator::creation(spin, orbital);
  }
  static Operator a(Operator::Spin spin, uint8_t orbital) {
    return Operator::annihilation(spin, orbital);
  }
};

TEST_F(NormalOrderTest, SingleOperator) {
  Term term(c(Operator::Spin::Up, 0));
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second, Term::coefficient_type(1.0f, 0.0f));
  EXPECT_EQ(result.terms().begin()->first, term.operators());
}

TEST_F(NormalOrderTest, TwoCommutingOperatorsAlreadyOrdered) {
  Term term({c(Operator::Spin::Up, 0), c(Operator::Spin::Down, 1)});
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second, Term::coefficient_type(1.0f, 0.0f));
  EXPECT_EQ(result.terms().begin()->first, term.operators());
}

TEST_F(NormalOrderTest, TwoCommutingOperatorsOutOfOrder) {
  Term term({c(Operator::Spin::Down, 1), c(Operator::Spin::Up, 0)});
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second,
            Term::coefficient_type(-1.0f, 0.0f));
  EXPECT_EQ(result.terms().begin()->first,
            Term::container_type(
                {c(Operator::Spin::Up, 0), c(Operator::Spin::Down, 1)}));
}

TEST_F(NormalOrderTest, TwoNonCommutingOperators) {
  Term term({a(Operator::Spin::Up, 0), c(Operator::Spin::Up, 0)});
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 2);

  auto it = result.terms().begin();
  EXPECT_EQ(it->second, Term::coefficient_type(-1.0f, 0.0f));
  EXPECT_EQ(it->first, Term::container_type({c(Operator::Spin::Up, 0),
                                             a(Operator::Spin::Up, 0)}));
  ++it;
  EXPECT_EQ(it->second, Term::coefficient_type(1.0f, 0.0f));
  EXPECT_EQ(it->first, Term::container_type());
}

TEST_F(NormalOrderTest, MultipleCommutingOperators) {
  Term term({a(Operator::Spin::Up, 1), c(Operator::Spin::Down, 0),
             c(Operator::Spin::Up, 0)});
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 1);
  EXPECT_EQ(result.terms().begin()->second,
            Term::coefficient_type(-1.0f, 0.0f));
  EXPECT_EQ(result.terms().begin()->first,
            Term::container_type({c(Operator::Spin::Up, 0),
                                  c(Operator::Spin::Down, 0),
                                  a(Operator::Spin::Up, 1)}));
}

TEST_F(NormalOrderTest, MultipleNonCommutingOperators) {
  Term term(42.0f, {a(Operator::Spin::Up, 1), a(Operator::Spin::Up, 0),
                    c(Operator::Spin::Up, 0), c(Operator::Spin::Up, 1)});
  Expression result = orderer.normal_order(term);
  Expression expected;
  expected += Term(42.0f, {});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), a(Operator::Spin::Up, 0)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 1), a(Operator::Spin::Up, 1)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), c(Operator::Spin::Up, 1),
                           a(Operator::Spin::Up, 0), a(Operator::Spin::Up, 1)});

  EXPECT_EQ(result, expected);
}

TEST_F(NormalOrderTest, MultipleNonCommutingPairs) {
  Term term(42.0f, {a(Operator::Spin::Up, 0), c(Operator::Spin::Up, 0),
                    a(Operator::Spin::Up, 1), c(Operator::Spin::Up, 1)});
  Expression result = orderer.normal_order(term);
  Expression expected;
  expected += Term(42.0f, {});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), a(Operator::Spin::Up, 0)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 1), a(Operator::Spin::Up, 1)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), c(Operator::Spin::Up, 1),
                           a(Operator::Spin::Up, 0), a(Operator::Spin::Up, 1)});

  EXPECT_EQ(result, expected);
}

TEST_F(NormalOrderTest, MultipleNonCommutingPairs2) {
  Term term(42.0f, {a(Operator::Spin::Up, 1), c(Operator::Spin::Up, 1),
                    a(Operator::Spin::Up, 0), c(Operator::Spin::Up, 0)});
  Expression result = orderer.normal_order(term);
  Expression expected;
  expected += Term(42.0f, {});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), a(Operator::Spin::Up, 0)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 1), a(Operator::Spin::Up, 1)});
  expected -= Term(42.0f, {c(Operator::Spin::Up, 0), c(Operator::Spin::Up, 1),
                           a(Operator::Spin::Up, 0), a(Operator::Spin::Up, 1)});

  EXPECT_EQ(result, expected);
}

TEST_F(NormalOrderTest, EmptyTerm) {
  Term term;
  Expression result = orderer.normal_order(term);
  EXPECT_EQ(result.size(), 0);
}

}  // namespace
}  // namespace qmutils
