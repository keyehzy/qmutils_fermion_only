#include "qmutils/expression.h"

namespace qmutils {

void Expression::normalize() {
  // 1) Remove terms with coefficients that are effectively zero
  // 2) Remove terms with adjacent identical operators
  for (auto it = m_terms.begin(); it != m_terms.end();) {
    if ((std::abs(it->second) <
         std::numeric_limits<coefficient_type::value_type>::epsilon()) ||
        (std::adjacent_find(it->first.begin(), it->first.end()) !=
         it->first.end())) {
      it = m_terms.erase(it);
    } else {
      ++it;
    }
  }
}

std::string Expression::to_string() const {
  std::ostringstream oss;
  bool first = true;
  for (const auto& [ops, coeff] : m_terms) {
    if (!first) {
      oss << " + ";
    }
    first = false;
    oss << coeff << " ";
    for (size_t i = 0; i < ops.size(); ++i) {
      oss << ops[i].to_string();
    }
  }
  return oss.str();
}

}  // namespace qmutils
