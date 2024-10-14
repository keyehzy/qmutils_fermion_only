#include "qmutils/basis.h"

#include <gtest/gtest.h>

namespace qmutils {
namespace {

class BasisTest : public ::testing::Test {
 protected:
  static Basis::operators_type create_operators(
      std::initializer_list<std::pair<Operator::Spin, uint8_t>> ops) {
    Basis::operators_type result;
    for (const auto& [spin, orbital] : ops) {
      result.push_back(Operator::creation(spin, orbital));
    }
    return result;
  }
};

TEST_F(BasisTest, ConstructorAndBasicProperties) {
  Basis basis(3, 2);
  EXPECT_EQ(basis.orbitals(), 3);
  EXPECT_EQ(basis.particles(), 2);
  EXPECT_EQ(basis.size(), qmutils_choose(6, 2));
}

TEST_F(BasisTest, EmptyBasis) {
  Basis basis(3, 0);
  EXPECT_EQ(basis.size(), 1);
  EXPECT_TRUE(
      basis.contains(Basis::operators_type{}));  // Only the vacuum state
}

TEST_F(BasisTest, FullBasis) {
  Basis basis(2, 4);
  EXPECT_EQ(basis.size(), 1);  // Only one fully occupied state
}

TEST_F(BasisTest, ContainsAndIndex) {
  Basis basis(2, 2);

  auto state1 =
      create_operators({{Operator::Spin::Up, 0}, {Operator::Spin::Down, 0}});
  auto state2 =
      create_operators({{Operator::Spin::Up, 0}, {Operator::Spin::Up, 1}});
  auto state3 =
      create_operators({{Operator::Spin::Up, 0}, {Operator::Spin::Down, 1}});

  EXPECT_TRUE(basis.contains(state1));
  EXPECT_TRUE(basis.contains(state2));
  EXPECT_TRUE(basis.contains(state3));
}

TEST_F(BasisTest, DoesNotContain) {
  Basis basis(2, 2);

  auto invalid_state1 =
      create_operators({{Operator::Spin::Up, 0},
                        {Operator::Spin::Up, 0}});  // Same orbital, same spin
  auto invalid_state2 =
      create_operators({{Operator::Spin::Up, 0},
                        {Operator::Spin::Up, 1},
                        {Operator::Spin::Down, 1}});  // Too many particles

  EXPECT_FALSE(basis.contains(invalid_state1));
  EXPECT_FALSE(basis.contains(invalid_state2));
}

TEST_F(BasisTest, OrderDependence) {
  Basis basis(2, 2);

  auto state1 =
      create_operators({{Operator::Spin::Up, 0}, {Operator::Spin::Down, 1}});
  auto state2 =
      create_operators({{Operator::Spin::Down, 1}, {Operator::Spin::Up, 0}});

  EXPECT_TRUE(basis.contains(state1));
  EXPECT_FALSE(basis.contains(state2));
}

TEST_F(BasisTest, CorrectSizeForVariousConfigurations) {
  const size_t max_orbital_count = 8;
  for (size_t orbitals = 0; orbitals < max_orbital_count; ++orbitals) {
    for (size_t particles = 0; particles < 2 * orbitals; ++particles) {
      Basis basis(orbitals, particles);
      EXPECT_EQ(basis.size(), qmutils_choose(2 * orbitals, particles));
    }
  }
}

TEST_F(BasisTest, EqualityOperator) {
  Basis basis1(2, 2);
  Basis basis2(2, 2);
  Basis basis3(3, 2);

  EXPECT_EQ(basis1, basis2);
  EXPECT_NE(basis1, basis3);
}

}  // namespace
}  // namespace qmutils
