#include <armadillo>

#include "qmutils/expression.h"
#include "qmutils/functional.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/sparse_matrix.h"
#include "qmutils/term.h"

using namespace qmutils;

template <size_t Lx, size_t Ly>
class SquareLatticeLandauLevel {
 public:
  using Index = StaticIndex<Lx, Ly>;

  SquareLatticeLandauLevel(float t, float phi) : m_t(t), m_phi(phi), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }

  const Index& index() const { return m_index; }

  float t() const { return m_t; }
  float phi() const { return m_phi; }

 private:
  float m_t;
  float m_phi;
  Index m_index;
  Expression m_hamiltonian;

  void construct_hamiltonian() {
    Operator::Spin spin = Operator::Spin::Up;  // Spin is not relevant here
    for (size_t i = 0; i < Lx; ++i) {
      std::complex<float> phase(0.0f, 2.0f * std::numbers::pi_v<float> * m_phi *
                                          static_cast<float>(i) /
                                          static_cast<float>(Lx));
      for (size_t j = 0; j < Ly; ++j) {
        size_t index_center = m_index.to_orbital(i, j);
        size_t index_right = m_index.to_orbital((i + 1) % Lx, j);
        size_t index_up = m_index.to_orbital(i, (j + 1) % Ly);
        m_hamiltonian +=
            m_t * Expression::Fermion::hopping(index_right, index_center, spin);
        m_hamiltonian +=
            m_t * Expression::Fermion::hopping(std::exp(phase), index_up,
                                               index_center, spin);
      }
    }
  }
};
static bool check_spin_conservation(const Term::container_type& ops) {
  int spin_up_count = 0;
  int spin_down_count = 0;

  for (const auto& op : ops) {
    if (op.spin() == Operator::Spin::Up) {
      spin_up_count += (op.type() == Operator::Type::Creation) ? 1 : -1;
    } else {
      spin_down_count += (op.type() == Operator::Type::Creation) ? 1 : -1;
    }
  }

  return spin_up_count == 0 && spin_down_count == 0;
}

static bool check_particle_conservation(const Term::container_type& ops) {
  int creation_count = 0;
  for (const auto& op : ops) {
    creation_count += (op.type() == Operator::Type::Creation) ? 1 : -1;
  }
  return creation_count == 0;
}

static std::vector<std::pair<float, float>> calculate_dos(
    const arma::fvec& eigenvalues, float sigma = 0.1f, size_t num_points = 1000,
    float padding_factor = 0.1f) {
  float E_min = eigenvalues.min();
  float E_max = eigenvalues.max();
  float padding = (E_max - E_min) * padding_factor;
  E_min -= padding;
  E_max += padding;

  float dE = (E_max - E_min) / static_cast<float>(num_points - 1);
  std::vector<std::pair<float, float>> dos(num_points);

  // Calculate DOS using Gaussian broadening
  const float normalization =
      1.0f / (sigma * std::sqrt(2.0f * std::numbers::pi_v<float>));

#pragma omp parallel for
  for (size_t i = 0; i < num_points; ++i) {
    float E = E_min + static_cast<float>(i) * dE;
    float rho = 0.0f;

    for (size_t j = 0; j < eigenvalues.n_elem; ++j) {
      float delta_E = (E - eigenvalues(j)) / sigma;
      rho += std::exp(-0.5f * delta_E * delta_E);
    }

    dos[i] = {E, rho * normalization};
  }

  return dos;
}

static std::vector<std::pair<float, float>> calculate_integrated_dos(
    const std::vector<std::pair<float, float>>& dos) {
  std::vector<std::pair<float, float>> integrated_dos(dos.size());
  float integral = 0.0f;
  float dE = dos[1].first - dos[0].first;

  for (size_t i = 0; i < dos.size(); ++i) {
    integral += dos[i].second * dE;
    integrated_dos[i] = {dos[i].first, integral};
  }

  return integrated_dos;
}

int main() {
  const size_t Lx = 5;
  const size_t Ly = 5;
  const size_t particles = 2;

  const float t = 1.0f;
  const float phi = 1.0f;

  SquareLatticeLandauLevel<Lx, Ly> model(t, phi);

  QMUTILS_ASSERT(check_expression_collective_predicate(check_spin_conservation,
                                                       model.hamiltonian()));
  QMUTILS_ASSERT(check_expression_collective_predicate(
      check_particle_conservation, model.hamiltonian()));

  for (const auto& [ops, coeff] : model.hamiltonian().terms()) {
    std::cout << term_to_string(Term(coeff, ops), model.index()) << "\n";
  }

  Basis basis(Lx * Ly, particles);
  auto H_matrix =
      compute_matrix_elements<arma::cx_fmat>(basis, model.hamiltonian());
  QMUTILS_ASSERT(arma::approx_equal(H_matrix, H_matrix.t(), "absdiff", 1e-4f));
  arma::fvec eigenvalues;
  arma::cx_fmat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, H_matrix);

  auto dos = calculate_dos(eigenvalues);
  auto integrated_dos = calculate_integrated_dos(dos);

  std::cout << "# DOS computed" << std::endl;

  std::ofstream dos_file("dos.dat");

  for (size_t i = 0; i < dos.size(); ++i) {
    dos_file << dos[i].first << " "
             << dos[i].second / static_cast<float>(basis.size()) << " "
             << integrated_dos[i].second / static_cast<float>(basis.size())
             << "\n";
  }

  return 0;
}
