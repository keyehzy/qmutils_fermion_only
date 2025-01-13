#include <armadillo>
#include <iostream>

#include "dos_utils.h"
#include "qmutils/expression.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/operator.h"

using namespace qmutils;

// Kagome lattice model

template <size_t Lx, size_t Ly>
class KagomeModel {
 public:
  using Index = StaticIndex<Lx, Ly, 3>;

  KagomeModel(float t1, float U) : m_t1(t1), m_U(U), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }
  const Index& index() const { return m_index; }

 private:
  float m_t1, m_U;
  Index m_index;
  Expression m_hamiltonian;

  Term density_density_interaction(size_t site) const {
    auto s = Operator::Spin::Up;
    return 0.5f * m_U *
           Term({Operator::Boson::creation(s, site),
                 Operator::Boson::creation(s, site),
                 Operator::Boson::annihilation(s, site),
                 Operator::Boson::annihilation(s, site)});
  }

  void construct_hamiltonian() {
    auto s = Operator::Spin::Up;

    for (size_t i = 0; i < Lx; ++i) {
      for (size_t j = 0; j < Ly; ++j) {
        size_t A_site = m_index.to_orbital(i, j, 0);
        size_t B_site = m_index.to_orbital(i, j, 1);
        size_t C_site = m_index.to_orbital(i, j, 2);

        size_t next_A_y_site = m_index.to_orbital(i, (j + 1) % Ly, 0);
        size_t next_A_x_site = m_index.to_orbital((i + 1) % Lx, j, 0);
        size_t next_C_x_y_site =
            m_index.to_orbital((i + 1) % Lx, (j + Ly - 1) % Ly, 2);

        m_hamiltonian += m_t1 * Expression::Boson::hopping(A_site, B_site, s);
        m_hamiltonian += m_t1 * Expression::Boson::hopping(B_site, C_site, s);
        m_hamiltonian += m_t1 * Expression::Boson::hopping(C_site, A_site, s);

        m_hamiltonian +=
            m_t1 * Expression::Boson::hopping(B_site, next_A_x_site, s);
        m_hamiltonian +=
            m_t1 * Expression::Boson::hopping(C_site, next_A_y_site, s);
        m_hamiltonian +=
            m_t1 * Expression::Boson::hopping(B_site, next_C_x_y_site, s);

        m_hamiltonian += density_density_interaction(A_site);
        m_hamiltonian += density_density_interaction(B_site);
        m_hamiltonian += density_density_interaction(C_site);
      }
    }
  }
};

int main() {
  const size_t Lx = 9;
  const size_t Ly = 9;
  const size_t P = 2;
  KagomeModel<Lx, Ly> model(1.0f, 0.0f);
  BosonicBasis basis(Lx * Ly * 3, P);  // 3 sites per unit cell, P particles
  auto H_matrix =
      compute_matrix_elements<arma::sp_cx_fmat>(basis, model.hamiltonian());
  arma::fvec eigenvalues;
  arma::eig_sym(eigenvalues, arma::cx_fmat(H_matrix));
  std::cout << eigenvalues << std::endl;

  //   Dos_Utils dos_ctx(eigenvalues, 0.1f, 1000, 0.1f);

  //   std::ofstream dos_file("dos.dat");
  //   for (size_t i = 0; i < dos_ctx.dos().size(); ++i) {
  //     dos_file << dos_ctx.dos()[i].first << " " << dos_ctx.dos()[i].second <<
  //     " "
  //              << dos_ctx.integrated_dos()[i].second << std::endl;
  //   }

  return 0;
}
