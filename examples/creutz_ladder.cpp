#include <armadillo>

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

    for (size_t i = 0; i < L - 1; ++i) {
      uint8_t A_site = m_index.to_orbital(i, 0);
      uint8_t B_site = m_index.to_orbital(i, 1);

      uint8_t next_A_site = m_index.to_orbital((i + 1) % L, 0);
      uint8_t next_B_site = m_index.to_orbital((i + 1) % L, 1);

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
      uint8_t A_site = m_index.to_orbital(i, 0);
      uint8_t B_site = m_index.to_orbital(i, 1);

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

int main() {
  for (float U = 0.0; U < 5.0f; U += 0.5f) {
    const size_t L = 16;
    CreutzLadderModel<L> model(1.0f, 0.5f * std::numbers::pi_v<float>, U);
    BosonicBasis basis(2 * L, 2);  // 2 sites per unit cell, 2 particles

    QMUTILS_ASSERT(
        model.hamiltonian().is_purely(Operator::Statistics::Bosonic));

    auto H_matrix =
        compute_matrix_elements<arma::sp_cx_fmat>(basis, model.hamiltonian());

    QMUTILS_ASSERT(
        arma::approx_equal(H_matrix, H_matrix.t(), "absdiff", 1e-4f));

    arma::fvec eigenvalues;
    arma::cx_fmat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, arma::cx_fmat(H_matrix));

    std::cout << U << eigenvalues.t();
  }

  return 0;
}