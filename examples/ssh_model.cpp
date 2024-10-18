#include <armadillo>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/operator.h"

using namespace qmutils;

template <size_t... Dims>
class IndexND {
 public:
  static constexpr size_t Dimensions = sizeof...(Dims);

  IndexND() {
    static_assert(total_size() <= 64,
                  "Lattice size exceeds maximum number of orbitals (64)");
  }

  template <typename... Coords>
  uint8_t to_orbital(Coords... coords) const {
    static_assert(sizeof...(coords) == Dimensions,
                  "Invalid number of coordinates");
    std::array<size_t, Dimensions> coordinates = {
        static_cast<size_t>(coords)...};

    for (size_t i = 0; i < Dimensions; ++i) {
      if (coordinates[i] >= m_dimensions[i]) {
        throw std::out_of_range("Coordinates out of lattice bounds");
      }
    }

    uint8_t orbital = 0;
    size_t multiplier = 1;
    for (size_t i = 0; i < Dimensions; ++i) {
      orbital += static_cast<uint8_t>(coordinates[i] * multiplier);
      multiplier *= m_dimensions[i];
    }
    return orbital;
  }

  std::array<size_t, Dimensions> from_orbital(uint8_t orbital) const {
    if (orbital >= total_size()) {
      throw std::out_of_range("Orbital index out of lattice bounds");
    }

    std::array<size_t, Dimensions> coordinates;
    for (size_t i = 0; i < Dimensions; ++i) {
      coordinates[i] = orbital % m_dimensions[i];
      orbital /= m_dimensions[i];
    }
    return coordinates;
  }

  constexpr const std::array<size_t, Dimensions>& dimensions() const {
    return m_dimensions;
  }

 private:
  static constexpr std::array<size_t, Dimensions> m_dimensions = {Dims...};

  static constexpr size_t total_size() {
    size_t size = 1;
    for (size_t dim : m_dimensions) {
      size *= dim;
    }
    return size;
  }
};

template <size_t L>
class SSHModel {
 public:
  SSHModel(float t1, float t2, float delta)
      : m_t1(t1), m_t2(t2), m_delta(delta), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }

  const IndexND<L, 2>& lattice() const { return m_index; }

 private:
  float m_t1, m_t2, m_delta;
  IndexND<L, 2> m_index;
  Expression m_hamiltonian;

  void construct_hamiltonian() {
    // Intra-cell hopping
    for (size_t i = 0; i < L; ++i) {
      add_hopping(i, 0, i, 1, m_t1 + m_delta);
    }

    // Inter-cell hopping
    for (size_t i = 0; i < L; ++i) {
      add_hopping(i, 1, (i + 1) % L, 0, m_t2 - m_delta);
    }
  }

  void add_hopping(size_t i1, size_t j1, size_t i2, size_t j2, float t) {
    uint8_t orbital1 = m_index.to_orbital(i1, j1);
    uint8_t orbital2 = m_index.to_orbital(i2, j2);
    m_hamiltonian +=
        t * Expression::hopping(orbital1, orbital2, Operator::Spin::Up);
    m_hamiltonian +=
        t * Expression::hopping(orbital1, orbital2, Operator::Spin::Down);
  }
};

template <size_t L>
Expression fourier_transform_operator(const Operator& op,
                                      const IndexND<L, 2>& lattice) {
  Expression result;
  const float type_sign =
      (op.type() == Operator::Type::Annihilation) ? -1.0f : 1.0f;

  constexpr float pi = std::numbers::pi_v<float>;
  const float factor = 2.0f * pi / static_cast<float>(L);

  auto [unit_cell, sublattice] = lattice.from_orbital(op.orbital());

  for (size_t k = 0; k < L; ++k) {
    std::complex<float> coefficient(
        0.0f, -type_sign * factor * static_cast<float>(k * unit_cell));
    coefficient = std::exp(coefficient) / std::sqrt(static_cast<float>(L));

    Operator transformed_op(op.type(), op.spin(),
                            lattice.to_orbital(k, sublattice));
    result += Term(coefficient, {transformed_op});
  }

  return result;
}

template <size_t L>
Expression fourier_transform_term(const Term& term,
                                  const IndexND<L, 2>& lattice) {
  Expression result(term.coefficient());
  for (const auto& op : term.operators()) {
    result *= fourier_transform_operator(op, lattice);
  }
  return result;
}

template <size_t L>
Expression fourier_transform(const Expression& expr,
                             const IndexND<L, 2>& lattice) {
  Expression result;
  for (const auto& [ops, coeff] : expr.terms()) {
    result += fourier_transform_term(Term(coeff, ops), lattice);
  }
  return result;
}

static void print_hamiltonian(const Expression& hamiltonian,
                              float epsilon = 3e-4f) {
  for (const auto& [ops, coeff] : hamiltonian.terms()) {
    if (std::abs(coeff) > epsilon) {
      std::cout << "  " << Term(coeff, ops).to_string() << std::endl;
    }
  }
  std::cout << std::endl;
}

static Expression transform_operator(const Operator& op,
                                     const arma::cx_fmat& eigenvectors,
                                     const Basis& basis) {
  Expression result;
  Operator needle(Operator::Type::Creation, op.spin(), op.orbital());
  size_t i = static_cast<size_t>(basis.index({needle}));
  for (size_t k = 0; k < basis.size(); ++k) {
    std::complex<float> coeff = (op.type() == Operator::Type::Annihilation)
                                    ? eigenvectors(i, k)
                                    : std::conj(eigenvectors(i, k));
    // Get spin and orbital indices from basis index k
    Operator::Spin spin = basis.at(k)[0].spin();
    size_t orbital = basis.at(k)[0].orbital();
    Term new_term(coeff, {Operator(op.type(), spin, orbital)});
    result += new_term;
  }
  return result;
}

static Expression transform_term(const Term& term,
                                 const arma::cx_fmat& eigenvectors,
                                 const Basis& basis) {
  Expression result(term.coefficient());
  for (const auto& op : term.operators()) {
    result *= transform_operator(op, eigenvectors, basis);
  }
  return result;
}

static Expression transform_to_diagonal(const Expression& expr,
                                        const arma::cx_fmat& eigenvectors,
                                        const Basis& basis) {
  Expression result;
  for (const auto& [ops, coeff] : expr.terms()) {
    result += transform_term(Term(coeff, ops), eigenvectors, basis);
  }
  return result;
}

int main() {
  const size_t L = 10;  // 10 unit cells
  float t1 = 1.0f, t2 = 0.5f, delta = 0.1f;

  SSHModel<L> ssh_model(t1, t2, delta);

  std::cout << "SSH Hamiltonian in real space:" << std::endl;
  print_hamiltonian(ssh_model.hamiltonian());

  Expression momentum_hamiltonian =
      fourier_transform(ssh_model.hamiltonian(), ssh_model.lattice());
  std::cout << "SSH Hamiltonian in momentum space:" << std::endl;
  print_hamiltonian(momentum_hamiltonian);

  Basis basis(2 * L, 1);  // Single-particle basis

  std::cout << "Basis:" << std::endl;
  for (const auto& elm : basis) {
    std::cout << Term(elm).to_string() << std::endl;
  }

  auto H_matrix =
      compute_matrix_elements<arma::cx_fmat>(basis, momentum_hamiltonian);
  arma::fvec eigenvalues;
  arma::cx_fmat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, H_matrix);

  std::cout << "Eigenvalues:" << std::endl;
  for (size_t i = 0; i < eigenvalues.n_elem; ++i) {
    std::cout << Term(basis.at(i)).to_string() << "  " << eigenvalues(i)
              << std::endl;
  }

  Expression diagonal_hamiltonian =
      transform_to_diagonal(momentum_hamiltonian, eigenvectors, basis);
  std::cout << "SSH Hamiltonian in diagonal form:" << std::endl;
  print_hamiltonian(diagonal_hamiltonian);

  return 0;
}
