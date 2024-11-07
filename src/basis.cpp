#include "qmutils/basis.h"

namespace qmutils {

void BasisBase::generate_basis() {
  operators_type current;
  current.reserve(m_particles);
  generate_combinations(current, 0, 0);
}

void FermionicBasis::generate_combinations(operators_type& current,
                                           size_t first_orbital, size_t depth) {
  if (depth == m_particles) {
    operators_type current_copy(current);
    std::sort(current_copy.begin(), current_copy.end());
    m_index_map.emplace_back(calculate_normalization(current_copy),
                             current_copy);
    return;
  }

  for (size_t i = first_orbital; i < m_orbitals; i++) {
    for (int spin_index = 0; spin_index < 2; ++spin_index) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      if (current.empty() || current.back().orbital() < i ||
          ((current.back().orbital() == i && spin > current.back().spin()))) {
        current.push_back(Operator::Fermion::creation(spin, i));
        generate_combinations(current, i, depth + 1);
        current.pop_back();
      }
    }
  }
}

void BosonicBasis::generate_combinations(operators_type& current,
                                         size_t first_orbital, size_t depth) {
  if (depth == m_particles) {
    operators_type current_copy(current);
    std::sort(current_copy.begin(), current_copy.end());
    m_index_map.emplace_back(calculate_normalization(current_copy),
                             current_copy);
    return;
  }

  for (size_t i = first_orbital; i < m_orbitals; i++) {
    for (int spin_index = 0; spin_index < 1; ++spin_index) {
      Operator::Spin spin = static_cast<Operator::Spin>(spin_index);
      if (current.empty() || current.back().orbital() <= i) {
        current.push_back(Operator::Boson::creation(spin, i));
        generate_combinations(current, i, depth + 1);
        current.pop_back();
      }
    }
  }
}

}  // namespace qmutils
