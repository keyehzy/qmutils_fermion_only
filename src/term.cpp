#include "qmutils/term.h"

#include <numeric>
#include <sstream>

namespace qmutils {

Term Term::adjoint() const {
  container_type adjoint_operators;
  adjoint_operators.reserve(m_operators.size());

  for (auto it = m_operators.rbegin(); it != m_operators.rend(); ++it) {
    adjoint_operators.push_back(it->adjoint());
  }

  return Term(std::conj(m_coefficient), std::move(adjoint_operators));
}

Term Term::flip_spin() const {
  container_type flipped_operators;
  flipped_operators.reserve(m_operators.size());
  for (const auto& op : m_operators) {
    flipped_operators.push_back(op.flip_spin());
  }
  return Term(m_coefficient, std::move(flipped_operators));
}

std::string Term::to_string() const {
  std::ostringstream oss;
  oss << m_coefficient << " ";
  for (size_t i = 0; i < m_operators.size(); ++i) {
    oss << m_operators[i].to_string();
  }
  return oss.str();
}

Term Term::creation(Operator::Spin spin, size_t orbital) {
  return Term(1.0f, {Operator::creation(spin, orbital)});
}

Term Term::annihilation(Operator::Spin spin, size_t orbital) {
  return Term(1.0f, {Operator::annihilation(spin, orbital)});
}

Term Term::one_body(Operator::Spin spin1, size_t orbital1, Operator::Spin spin2,
                    size_t orbital2) {
  return Term(1.0f, {Operator::creation(spin1, orbital1),
                     Operator::annihilation(spin2, orbital2)});
}

Term Term::density(Operator::Spin spin, size_t orbital) {
  return Term::one_body(spin, orbital, spin, orbital);
}

Term density_density(Operator::Spin spin1, size_t i, Operator::Spin spin2,
                     size_t j) {
  return Term::density(spin1, i) * Term::density(spin2, j);
}

}  // namespace qmutils
