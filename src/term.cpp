#include "term.h"

#include <numeric>
#include <sstream>

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

}  // namespace qmutils
