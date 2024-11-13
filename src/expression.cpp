#include "qmutils/expression.h"

#include "qmutils/assert.h"

namespace qmutils {

static bool constexpr both_are_fermions_and_equal(const Operator& a,
                                                  const Operator& b) {
  return Operator::is_fermion(a) && Operator::is_fermion(b) && a == b;
}

void Expression::normalize() {
  // 1) Remove terms with coefficients that are effectively zero
  // 2) In the case of Fermions, remove terms with adjacent identical operators
  for (auto it = m_terms.begin(); it != m_terms.end();) {
    if ((std::abs(it->second) <
         std::numeric_limits<coefficient_type::value_type>::epsilon()) ||
        (std::adjacent_find(it->first.begin(), it->first.end(),
                            both_are_fermions_and_equal) != it->first.end())) {
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

Expression Expression::adjoint() const {
  Expression result;
  for (const auto& [ops, coeff] : terms()) {
    Term::container_type adjoint_operators;
    adjoint_operators.reserve(ops.size());
    for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
      adjoint_operators.push_back(it->adjoint());
    }
    result.m_terms[adjoint_operators] += std::conj(coeff);
  }
  return result;
}

Expression Expression::flip_spin() const {
  Expression result;
  for (const auto& [ops, coeff] : m_terms) {
    Term flipped_term(coeff, ops);
    result += flipped_term.flip_spin();
  }
  return result;
}

Expression Expression::hopping(size_t from_orbital, size_t to_orbital,
                               Operator::Spin spin) {
  QMUTILS_ASSERT(from_orbital != to_orbital);
  Expression result;
  result += Term::one_body(spin, from_orbital, spin, to_orbital);
  result += Term::one_body(spin, to_orbital, spin, from_orbital);
  return result;
}

Expression Expression::hopping(coefficient_type coeff, size_t from_orbital,
                               size_t to_orbital, Operator::Spin spin) {
  QMUTILS_ASSERT(from_orbital != to_orbital);
  Expression result;
  result += coeff * Term::one_body(spin, from_orbital, spin, to_orbital);
  result +=
      std::conj(coeff) * Term::one_body(spin, to_orbital, spin, from_orbital);
  return result;
}

}  // namespace qmutils
