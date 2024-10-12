#pragma once

#include <queue>
#include <unordered_map>

#include "qmutils/expression.h"
#include "qmutils/term.h"

namespace qmutils {

class NormalOrderer {
 public:
  using coefficient_type = Term::coefficient_type;
  using operators_type = Term::container_type;

  NormalOrderer() = default;

  Expression normal_order(const Term& term);
  Expression normal_order(const Expression& expr);

  void print_cache_stats() const;

 private:
  Expression normal_order_iterative(const operators_type& ops);
  Expression normal_order_recursive(operators_type ops);
  Expression handle_non_commuting(const operators_type& ops, size_t index);

  std::unordered_map<operators_type, Expression> m_cache;

  struct QueueElement {
    operators_type ops;
    coefficient_type phase;

    QueueElement(const operators_type& ops, coefficient_type phase)
        : ops(ops), phase(phase) {}

    bool operator<(const QueueElement& other) const {
      return other.ops.size() < this->ops.size();
    }
  };

  std::priority_queue<QueueElement> m_queue;
  size_t m_cache_hits{0};
  size_t m_cache_misses{0};
};

}  // namespace qmutils
