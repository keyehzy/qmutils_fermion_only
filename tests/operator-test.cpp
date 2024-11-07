#include "qmutils/operator.h"

#include <gtest/gtest.h>

namespace qmutils {
namespace {

TEST(OperatorTest, ConstructionAndProperties) {
  Operator op1 = Operator::creation(Operator::Spin::Up, 5);
  EXPECT_EQ(op1.type(), Operator::Type::Creation);
  EXPECT_EQ(op1.spin(), Operator::Spin::Up);
  EXPECT_EQ(op1.orbital(), 5);

  Operator op2 = Operator::annihilation(Operator::Spin::Down, 3);
  EXPECT_EQ(op2.type(), Operator::Type::Annihilation);
  EXPECT_EQ(op2.spin(), Operator::Spin::Down);
  EXPECT_EQ(op2.orbital(), 3);
}

TEST(OperatorTest, ComparisonOperators) {
  Operator op1 = Operator::creation(Operator::Spin::Up, 5);
  Operator op2 = Operator::creation(Operator::Spin::Up, 5);
  Operator op3 = Operator::creation(Operator::Spin::Down, 5);
  Operator op4 = Operator::annihilation(Operator::Spin::Up, 5);

  EXPECT_EQ(op1, op2);
  EXPECT_NE(op1, op3);
  EXPECT_NE(op1, op4);
  EXPECT_LT(op1, op3);  // Test operator<
}

TEST(OperatorTest, CommutationRelations) {
  Operator c_up_5 = Operator::creation(Operator::Spin::Up, 5);
  Operator c_down_5 = Operator::creation(Operator::Spin::Down, 5);
  Operator c_up_6 = Operator::creation(Operator::Spin::Up, 6);
  Operator a_up_5 = Operator::annihilation(Operator::Spin::Up, 5);

  EXPECT_TRUE(c_up_5.commutes_with(c_down_5));
  EXPECT_TRUE(c_up_5.commutes_with(c_up_6));
  EXPECT_FALSE(c_up_5.commutes_with(a_up_5));
}

TEST(OperatorTest, COmmutationRelationBetweenOperatorsWithDifferentStatistics) {
  Operator fermion_up_5 = Operator::Fermion::creation(Operator::Spin::Up, 5);
  Operator boson_up_5 = Operator::Boson::annihilation(Operator::Spin::Up, 5);

  EXPECT_TRUE(fermion_up_5.commutes_with(boson_up_5));
}

TEST(OperatorTest, AdjointOperation) {
  Operator c_up_5 = Operator::creation(Operator::Spin::Up, 5);
  Operator a_up_5 = Operator::annihilation(Operator::Spin::Up, 5);

  EXPECT_EQ(c_up_5.adjoint(), a_up_5);
  EXPECT_EQ(a_up_5.adjoint(), c_up_5);
}

TEST(OperatorTest, StringRepresentation) {
  Operator c_up_5 = Operator::creation(Operator::Spin::Up, 5);
  Operator a_down_3 = Operator::annihilation(Operator::Spin::Down, 3);

  EXPECT_EQ(c_up_5.to_string(), "c+(↑,5)");
  EXPECT_EQ(a_down_3.to_string(), "c(↓,3)");
}

TEST(OperatorTest, DataRepresentationAndHashing) {
  Operator op1 = Operator::creation(Operator::Spin::Up, 5);
  Operator op2 = Operator::creation(Operator::Spin::Up, 5);
  Operator op3 = Operator::annihilation(Operator::Spin::Down, 3);

  EXPECT_EQ(op1.data(), op2.data());
  EXPECT_NE(op1.data(), op3.data());

  std::hash<Operator> hasher;
  EXPECT_EQ(hasher(op1), hasher(op2));
  EXPECT_NE(hasher(op1), hasher(op3));
}

TEST(OperatorTest, OperatorFlipSpin) {
  Operator op_up = Operator::creation(Operator::Spin::Up, 0);
  Operator op_down = Operator::annihilation(Operator::Spin::Down, 1);

  EXPECT_EQ(op_up.flip_spin(), Operator::creation(Operator::Spin::Down, 0));
  EXPECT_EQ(op_down.flip_spin(), Operator::annihilation(Operator::Spin::Up, 1));
}

}  // namespace
}  // namespace qmutils
