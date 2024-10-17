#pragma once

#include <algorithm>
#include <complex>
#include <unordered_map>

#include "qmutils/term.h"

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

  // Arithmetic operations
  Expression& operator+=(const Expression& other) {
    for (const auto& [ops, coeff] : other.m_terms) {
      m_terms[ops] += coeff;
    }
    normalize();
    return *this;
  }

  Expression& operator+=(const Term& other) {
    m_terms[other.operators()] += other.coefficient();
    normalize();
    return *this;
  }

  Expression& operator+=(Operator other) {
    m_terms[{other}] += coefficient_type(1.0);
    normalize();
    return *this;
  }

  Expression& operator+=(coefficient_type other) {
    m_terms[{}] += coefficient_type(other);
    normalize();
    return *this;
  }

  Expression& operator-=(const Expression& other) {
    for (const auto& [ops, coeff] : other.m_terms) {
      m_terms[ops] -= coeff;
    }
    normalize();
    return *this;
  }

  Expression& operator-=(const Term& other) {
    m_terms[other.operators()] -= other.coefficient();
    normalize();
    return *this;
  }

  Expression& operator*=(const Expression& other) {
    std::unordered_map<operators_type, coefficient_type> new_terms;
    for (const auto& [ops1, coeff1] : m_terms) {
      for (const auto& [ops2, coeff2] : other.m_terms) {
        operators_type new_ops;
        new_ops.reserve(ops1.size() + ops2.size());
        new_ops.insert(new_ops.end(), ops1.begin(), ops1.end());
        new_ops.insert(new_ops.end(), ops2.begin(), ops2.end());
        new_terms[std::move(new_ops)] += coeff1 * coeff2;
      }
    }
    m_terms = std::move(new_terms);
    normalize();
    return *this;
  }

  Expression& operator*=(const Term& term) {
    std::unordered_map<operators_type, coefficient_type> new_terms;
    for (const auto& [ops, coeff] : m_terms) {
      operators_type new_ops;
      new_ops.reserve(ops.size() + term.operators().size());
      new_ops.insert(new_ops.end(), ops.begin(), ops.end());
      new_ops.insert(new_ops.end(), term.operators().begin(),
                     term.operators().end());
      new_terms[new_ops] += coeff * term.coefficient();
    }

    m_terms = std::move(new_terms);
    normalize();
    return *this;
  }

  Expression& operator*=(const coefficient_type& scalar) {
    for (auto& [ops, coeff] : m_terms) {
      coeff *= scalar;
    }
    normalize();
    return *this;
  }

  friend Expression operator+(Expression lhs, const Expression& rhs) {
    lhs += rhs;
    return lhs;
  }

  friend Expression operator-(Expression lhs, const Expression& rhs) {
    lhs -= rhs;
    return lhs;
  }

  friend Expression operator*(Expression lhs, const Expression& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend Expression operator*(Expression lhs, const coefficient_type& rhs) {
    lhs *= rhs;
    return lhs;
  }

  friend Expression operator*(const coefficient_type& lhs, Expression rhs) {
    rhs *= lhs;
    return rhs;
  }

  const coefficient_type& operator[](const operators_type& ops) {
    return m_terms[ops];
  }

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

  std::string to_string() const;

  void normalize();

  Expression adjoint() const;

  Expression flip_spin() const;

  static Expression hopping(uint8_t from_orbital, uint8_t to_orbital,
                            Operator::Spin spin);

 private:
  std::unordered_map<operators_type, coefficient_type> m_terms;
};

}  // namespace qmutils
