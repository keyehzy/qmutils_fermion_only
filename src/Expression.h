#pragma once

#include <complex>
#include <unordered_map>

#include "Term.h"

namespace qmutils {

class Expression {
 public:
  using coefficient_type = Term::coefficient_type;
  using operators_type = Term::container_type;

  Expression() = default;

  explicit Expression(const Term& term) {
    m_terms[term.operators()] = term.coefficient();
  }

  explicit Expression(Operator op) {
    m_terms[{op}] = coefficient_type(1.0f, 0.0f);
  }

  Expression(const coefficient_type& scalar) { m_terms[{}] = scalar; }

  Expression(const Expression&) = default;

  Expression(Expression&&) noexcept = default;

  Expression& operator=(const Expression&) = default;

  Expression& operator=(Expression&&) noexcept = default;

  size_t size() const { return m_terms.size(); }

  const std::unordered_map<operators_type, coefficient_type>& terms() const {
    return m_terms;
  }

  std::unordered_map<operators_type, coefficient_type>& terms() {
    return m_terms;
  }

  bool operator==(const Expression& other) const {
    return m_terms == other.m_terms;
  }

  bool operator!=(const Expression& other) const { return !(*this == other); }

 private:
  std::unordered_map<operators_type, coefficient_type> m_terms;
};

}  // namespace qmutils
