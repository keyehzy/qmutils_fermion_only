#include "term.h"

#include <numeric>
#include <sstream>

namespace qmutils {

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
