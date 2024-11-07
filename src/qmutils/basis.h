#pragma once

#include <algorithm>
#include <unordered_set>

#include "qmutils/assert.h"
#include "qmutils/operator.h"
#include "qmutils/term.h"

namespace qmutils {

static constexpr uint64_t qmutils_choose(uint64_t n, uint64_t m) {
  if (m > n) return 0;
  if (m == 0 || m == n) return 1;

  if (m > n - m) m = n - m;

  uint64_t result = 1;
  for (uint64_t i = 0; i < m; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }

  return result;
}

class Basis {
 public:
  using operators_type = Term::container_type;

  Basis(size_t orbitals, size_t particles)
      : m_orbitals(orbitals), m_particles(particles) {
    QMUTILS_ASSERT(orbitals <= Operator::max_orbital_size());
    QMUTILS_ASSERT(particles <= 2 * orbitals);
    m_index_map.reserve(qmutils_choose(2 * orbitals, particles));
    generate_basis();
    QMUTILS_ASSERT(m_index_map.size() ==
                   qmutils_choose(2 * orbitals, particles));
    std::sort(m_index_map.begin(), m_index_map.end(), term_sorter);
  }

  Basis(const Basis&) = default;
  Basis& operator=(const Basis&) = default;
  Basis(Basis&&) noexcept = default;
  Basis& operator=(Basis&&) noexcept = default;

  ~Basis() = default;

  size_t orbitals() const noexcept { return m_orbitals; }
  size_t particles() const noexcept { return m_particles; }
  size_t size() const noexcept { return m_index_map.size(); }

  bool operator==(const Basis& other) const {
    return m_orbitals == other.m_orbitals && m_particles == other.m_particles &&
           m_index_map == other.m_index_map;
  }

  bool operator!=(const Basis& other) const { return !(*this == other); }

  bool contains(const operators_type& value) const {
    auto it = std::lower_bound(m_index_map.begin(), m_index_map.end(),
                               Term(value), term_sorter);
    return it != m_index_map.end() && it->operators() == value;
  }

  // TODO: needs testing
  ptrdiff_t index_of(const operators_type& value) const {
    QMUTILS_ASSERT(contains(value));
    auto it = std::lower_bound(m_index_map.begin(), m_index_map.end(),
                               Term(value), term_sorter);
    return std::distance(m_index_map.begin(), it);
  }

  void insert(const operators_type& value) {
    QMUTILS_ASSERT(!contains(value));
    Term term(value);
    auto it = std::lower_bound(m_index_map.begin(), m_index_map.end(), term,
                               term_sorter);
    m_index_map.insert(it, term);
  }

  auto begin() const noexcept { return m_index_map.begin(); }
  auto end() const noexcept { return m_index_map.end(); }

  const Term& at(size_t i) const noexcept {
    QMUTILS_ASSERT(i < m_index_map.size());
    return m_index_map[i];
  }

 private:
  void generate_basis();

  void generate_combinations(operators_type& current, size_t first_orbital,
                             size_t depth);

  static bool term_sorter(const Term& a, const Term& b) noexcept {
    return a.operators() < b.operators();
  }

  size_t m_orbitals;
  size_t m_particles;
  std::vector<Term> m_index_map;
};

}  // namespace qmutils
