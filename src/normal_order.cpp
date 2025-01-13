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

static constexpr float phase_factor(const Operator& a, const Operator& b) {
  return Operator::is_fermion(a) && Operator::is_fermion(b) ? -1.0f : 1.0f;
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

  if (std::is_sorted(ops.begin(), ops.end())) {
    Expression result(1.0f, ops);
    m_cache.put(ops_hash, result);
    return result;
  }

  std::unique_ptr<operators_type> local_ops_ptr = m_pool.acquire();
  operators_type& local_ops = *local_ops_ptr;
  local_ops.assign(ops.begin(), ops.end());
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
        m_pool.release(std::move(local_ops_ptr));
        return result;
      }
    }
  }

  Expression result(phase, local_ops);
  m_cache.put(ops_hash, result);
  m_pool.release(std::move(local_ops_ptr));
  return result;
}

Expression NormalOrderer::handle_non_commuting(const operators_type& ops,
                                               size_t index, float sign) {
  std::unique_ptr<operators_type> contracted_ptr = m_pool.acquire();
  operators_type& contracted = *contracted_ptr;
  contracted.assign(ops.begin(), ops.begin() + index);
  contracted.insert(contracted.end(), ops.begin() + index + 2, ops.end());

  std::unique_ptr<operators_type> swapped_ptr = m_pool.acquire();
  operators_type& swapped = *swapped_ptr;
  swapped.assign(ops.begin(), ops.end());
  std::swap(swapped[index], swapped[index + 1]);

  Expression result = normal_order_recursive(contracted) +
                      sign * normal_order_recursive(swapped);
  m_pool.release(std::move(contracted_ptr));
  m_pool.release(std::move(swapped_ptr));
  return result;
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
