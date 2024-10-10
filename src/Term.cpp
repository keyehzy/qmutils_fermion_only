#include "Term.h"

#include <numeric>
#include <sstream>

#include "Assert.h"

namespace qmutils {

Term Term::adjoint() const noexcept {
  container_type adjoint_operators;
  adjoint_operators.reserve(m_operators.size());

  for (auto it = m_operators.rbegin(); it != m_operators.rend(); ++it) {
    adjoint_operators.push_back(it->adjoint());
  }

  return Term(std::conj(m_coefficient), std::move(adjoint_operators));
}

std::string Term::to_string() const {
  std::ostringstream oss;
  oss << m_coefficient.real();
  if (m_coefficient.imag() != 0.0f) {
    oss << (m_coefficient.imag() > 0 ? "+" : "") << m_coefficient.imag() << "i";
  }
  oss << " * ";
  for (const auto& op : m_operators) {
    oss << op.to_string() << " ";
  }
  return oss.str();
}

Term Term::creation(Operator::Spin spin, uint8_t orbital) {
  return Term(1.0f, {Operator::creation(spin, orbital)});
}

Term Term::annihilation(Operator::Spin spin, uint8_t orbital) {
  return Term(1.0f, {Operator::annihilation(spin, orbital)});
}

Term Term::number(Operator::Spin spin, uint8_t orbital) {
  return Term(1.0f, {Operator::creation(spin, orbital),
                     Operator::annihilation(spin, orbital)});
}

Term Term::spin_flip(uint8_t orbital) {
  return Term(1.0f, {Operator::creation(Operator::Spin::Up, orbital),
                     Operator::annihilation(Operator::Spin::Down, orbital)});
}

Term Term::hopping(uint8_t from_orbital, uint8_t to_orbital,
                   Operator::Spin spin) {
  QMUTILS_ASSERT(from_orbital != to_orbital);
  return Term(1.0f, {Operator::creation(spin, to_orbital),
                     Operator::annihilation(spin, from_orbital)});
}

}  // namespace qmutils
