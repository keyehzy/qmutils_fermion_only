#include "qmutils/normal_order.h"

#include <iostream>

namespace qmutils {

Expression NormalOrderer::normal_order(const Term& term) {
  return term.coefficient() * normal_order_recursive(term.operators());
}

Expression NormalOrderer::normal_order(const Expression& expr) {
  Expression result;
  for (const auto& [ops, coeff] : expr.terms()) {
    result += coeff * normal_order_recursive(ops);
  }
  return result;
}

static float phase_factor(const Operator& a, const Operator& b) {
  return Operator::is_fermion(a) && Operator::is_fermion(b) ? -1.0f : 1.0f;
}

Expression NormalOrderer::normal_order_iterative(const operators_type& ops) {
  Expression result;
  m_queue.emplace(ops, 1.0);

  while (!m_queue.empty()) {
    auto [current, phase] = m_queue.top();
    m_queue.pop();

    if (current.size() < 2) {
      result += Term(phase, current);
      continue;
    }

    auto ops_key = std::hash<operators_type>{}(current);
    if (auto cached_result = m_cache.get(ops_key)) {
      m_cache_hits++;
      result += phase * cached_result.value();
      continue;
    }
    m_cache_misses++;

    bool is_sorted = true;

    for (size_t i = 1; i < current.size(); ++i) {
      size_t j = i;
      while (j > 0 && current[j] < current[j - 1]) {
        float sign = phase_factor(current[j], current[j - 1]);
        if (current[j].commutes_with(current[j - 1])) {
          std::swap(current[j], current[j - 1]);
          phase *= sign;
          --j;
        } else {
          is_sorted = false;
          operators_type contracted(current);
          contracted.erase(contracted.begin() + j - 1,
                           contracted.begin() + j + 1);
          m_queue.emplace(std::move(contracted), phase);

          operators_type swapped(current);
          std::swap(swapped[j - 1], swapped[j]);
          m_queue.emplace(std::move(swapped), sign * phase);
          break;
        }
      }
      if (!is_sorted) break;
    }
    if (is_sorted) {
      Expression partial(Term(phase, current));
      m_cache.put(ops_key, partial);
      result += partial;
    }
  }

  return result;
}

Expression NormalOrderer::normal_order_recursive(const operators_type& ops) {
  if (ops.size() < 2) {
    return Expression(Term(ops));
  }

  std::size_t ops_key = std::hash<operators_type>{}(ops);
  return normal_order_recursive(ops, ops_key);
}

Expression NormalOrderer::normal_order_recursive(const operators_type& ops,
                                                 std::size_t ops_hash) {
  if (auto cached_result = m_cache.get(ops_hash)) {
    m_cache_hits++;
    return cached_result.value();
  }
  m_cache_misses++;

  operators_type local_ops = ops;
  coefficient_type phase = 1.0f;

  for (size_t i = 1; i < local_ops.size(); ++i) {
    size_t j = i;
    while (j > 0 && local_ops[j] < local_ops[j - 1]) {
      float sign = phase_factor(local_ops[j], local_ops[j - 1]);
      if (local_ops[j].commutes_with(local_ops[j - 1])) {
        std::swap(local_ops[j], local_ops[j - 1]);
        phase *= sign;
        --j;
      } else {
        Expression result =
            phase * handle_non_commuting(local_ops, j - 1, sign);
        m_cache.put(ops_hash, result);
        return result;
      }
    }
  }

  Expression result(Term(phase, local_ops));
  m_cache.put(ops_hash, result);
  return result;
}

Expression NormalOrderer::handle_non_commuting(const operators_type& ops,
                                               size_t index, float sign) {
  operators_type contracted(ops);
  contracted.erase(contracted.begin() + index, contracted.begin() + index + 2);

  operators_type swapped(ops);
  std::swap(swapped[index], swapped[index + 1]);

  return normal_order_recursive(contracted) +
         sign * normal_order_recursive(swapped);
}

void NormalOrderer::print_cache_stats() const {
  std::cout << "Total entries: " << m_cache.size() << std::endl;
  std::cout << "Cache hits: " << m_cache_hits << std::endl;
  std::cout << "Cache misses: " << m_cache_misses << std::endl;
  double hit_ratio =
      static_cast<double>(m_cache_hits) / (m_cache_hits + m_cache_misses);
  std::cout << "Hit ratio: " << hit_ratio << std::endl;
}

}  // namespace qmutils
