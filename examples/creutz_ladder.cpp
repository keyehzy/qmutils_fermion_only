#include <armadillo>
#include <fstream>

#include "qmutils/assert.h"
#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"

using namespace qmutils;

template <size_t L>
class CreutzLadderModel {
 public:
  using Index = StaticIndex<L, 2>;  // L unit cells, 2 sites per cell

  CreutzLadderModel(float J, float theta, float U)
      : m_J(J), m_theta(theta), m_U(U), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }
  const Index& index() const { return m_index; }

 private:
  float m_J, m_theta, m_U;
  Index m_index;
  Expression m_hamiltonian;

  void construct_hamiltonian() {
    auto s = Operator::Spin::Up;

    for (size_t i = 0; i < L; ++i) {
      size_t A_site = m_index.to_orbital(i, 0);
      size_t B_site = m_index.to_orbital(i, 1);

      size_t next_A_site = m_index.to_orbital((i + 1) % L, 0);
      size_t next_B_site = m_index.to_orbital((i + 1) % L, 1);

      m_hamiltonian += m_J * Expression::Boson::hopping(A_site, next_B_site, s);
      m_hamiltonian += m_J * Expression::Boson::hopping(B_site, next_A_site, s);

      std::complex<float> phase(std::cos(m_theta), std::sin(m_theta));

      Term A_A = m_J * phase * Term::Boson::one_body(s, A_site, s, next_A_site);
      m_hamiltonian += A_A;
      m_hamiltonian += A_A.adjoint();

      Term B_B = m_J * std::conj(phase) *
                 Term::Boson::one_body(s, B_site, s, next_B_site);
      m_hamiltonian += B_B;
      m_hamiltonian += B_B.adjoint();
    }

    for (size_t i = 0; i < L; ++i) {
      size_t A_site = m_index.to_orbital(i, 0);
      size_t B_site = m_index.to_orbital(i, 1);

      m_hamiltonian += 0.5f * m_U *
                       Term({Operator::Boson::creation(s, A_site),
                             Operator::Boson::creation(s, A_site),
                             Operator::Boson::annihilation(s, A_site),
                             Operator::Boson::annihilation(s, A_site)});

      m_hamiltonian += 0.5f * m_U *
                       Term({Operator::Boson::creation(s, B_site),
                             Operator::Boson::creation(s, B_site),
                             Operator::Boson::annihilation(s, B_site),
                             Operator::Boson::annihilation(s, B_site)});
    }
  }
};

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
  const size_t L = 24;
  const size_t P = 2;
  const float J = 1.0f;
  const float U = 6.0f * J;

  CreutzLadderModel<L> model(J, 0.5f * std::numbers::pi_v<float>, U);
  BosonicBasis basis(2 * L, P);  // 2 sites per unit cell, 2 particles

  std::cout << "# Basis size: " << basis.size() << std::endl;

  QMUTILS_ASSERT(model.hamiltonian().is_purely(Operator::Statistics::Bosonic));

  auto H_matrix =
      compute_matrix_elements<arma::cx_fmat>(basis, model.hamiltonian());

  std::cout << "# Matrix elements computed" << std::endl;

  QMUTILS_ASSERT(arma::approx_equal(H_matrix, H_matrix.t(), "absdiff", 1e-4f));

  arma::fvec eigenvalues;
  arma::cx_fmat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, H_matrix);

  std::cout << "# Eigenvalues computed" << std::endl;

  auto dos = calculate_dos(eigenvalues, 0.1f * J);
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
