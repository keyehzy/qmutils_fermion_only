// operator_pool.h

#pragma once
#include <memory>
#include <queue>
#include <vector>

#include "qmutils/term.h"

namespace qmutils {

class OperatorPool {
 public:
  using operators_type = Term::container_type;

  explicit OperatorPool(std::size_t initial_pool_size = 64) {
    resize(initial_pool_size);
  }

  ~OperatorPool() {
    for (auto& ptr : m_pool) {
      ptr->clear();
    }
  }

  OperatorPool(const OperatorPool&) = delete;
  OperatorPool& operator=(const OperatorPool&) = delete;

  OperatorPool(OperatorPool&&) noexcept = default;
  OperatorPool& operator=(OperatorPool&&) noexcept = default;

  void resize(std::size_t new_size) {
    m_pool.reserve(new_size);
    for (std::size_t i = m_pool.size(); i < new_size; ++i) {
      m_pool.push_back(std::make_unique<operators_type>());
    }
  }

  std::unique_ptr<operators_type> acquire() {
    if (!m_pool.empty()) {
      auto ptr = std::move(m_pool.back());
      m_pool.pop_back();
      ptr->clear();
      return ptr;
    }
    return std::make_unique<operators_type>();
  }

  void release(std::unique_ptr<operators_type> ops) {
    ops->clear();
    m_pool.push_back(std::move(ops));
  }

 private:
  std::vector<std::unique_ptr<operators_type>> m_pool;
};

}  // namespace qmutils
