#pragma once

#include <functional>
#include <type_traits>

#include "qmutils/expression.h"
#include "qmutils/operator.h"
#include "qmutils/term.h"

namespace qmutils {

template <typename F, typename... Args>
using is_operator_callable = std::is_invocable<F, const Operator&, Args...>;

template <typename F, typename... Args>
auto transform_operator(F&& f, const Operator& op, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value,
                        decltype(f(op, std::forward<Args>(args)...))> {
  return f(op, std::forward<Args>(args)...);
}

template <typename F, typename... Args>
auto transform_term(F&& f, const Term& term, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result(term.coefficient());
  for (const auto& op : term.operators()) {
    result *=
        transform_operator(std::forward<F>(f), op, std::forward<Args>(args)...);
  }
  return result;
}

template <typename F, typename... Args>
auto transform_expression(F&& f, const Expression& expr, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result;
  for (const auto& [ops, coeff] : expr.terms()) {
    result += transform_term(std::forward<F>(f), Term(coeff, ops),
                             std::forward<Args>(args)...);
  }
  return result;
}

}  // namespace qmutils
