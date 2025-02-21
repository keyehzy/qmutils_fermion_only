#include <armadillo>
#include <fstream>

#include "dos_utils.h"
#include "qmutils/assert.h"
#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/functional.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "timeintegrator.h"

using namespace qmutils;

// Reference: arXiv:1007.4640

template <size_t L>
class SawtoothModel {
 public:
  using Index = StaticIndex<L, 2>;  // L unit cells, 2 sites per cell

  SawtoothModel(float t1, float t2, float U)
      : m_t1(t1), m_t2(t2), m_U(U), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }
  const Index& index() const { return m_index; }

 private:
  float m_t1, m_t2, m_U;
  Index m_index;
  Expression m_hamiltonian;

  void construct_hamiltonian() {
    auto s = Operator::Spin::Up;

    for (size_t i = 0; i < L; ++i) {
      size_t A_site = m_index.to_orbital(i, 0);
      size_t B_site = m_index.to_orbital(i, 1);
      size_t next_A_site = m_index.to_orbital((i + 1) % L, 0);

      m_hamiltonian += m_t1 * Expression::Boson::hopping(A_site, B_site, s);
      m_hamiltonian +=
          m_t1 * Expression::Boson::hopping(B_site, next_A_site, s);
      m_hamiltonian +=
          m_t2 * Expression::Boson::hopping(A_site, next_A_site, s);
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

 public:
  Expression s(size_t i) {
    auto s = Operator::Spin::Up;
    Expression result;
    size_t A_site = m_index.to_orbital(i, 0);
    result += Term::Boson::creation(s, A_site);
    return result;
  }

  Expression cls(size_t i) {
    auto s = Operator::Spin::Up;
    Expression result;
    size_t A_site = m_index.to_orbital(i, 0);
    size_t B_site = m_index.to_orbital(i, 1);
    size_t previous_B_site = m_index.to_orbital((i + L - 1) % L, 1);
    result += 0.5f * std::sqrt(2.0f) * Term::Boson::creation(s, A_site);
    result -= 0.5f * Term::Boson::creation(s, B_site);
    result -= 0.5f * Term::Boson::creation(s, previous_B_site);
    return result;
  }

  Expression density_operator(size_t i) {
    auto s = Operator::Spin::Up;
    Expression result;
    size_t A_site = m_index.to_orbital(i, 0);
    result += Term::Boson::density(s, A_site);
    return result;
  }
};

int main() {
  const size_t L = 35;
  const float t2 = 1.0f;
  const float t1 = t2 * std::sqrt(2.0f);
  const float U = 1.0f;
  BosonicBasis basis(2 * L, 2);  // 2 sites per unit cell, P particles

  SawtoothModel<L> model(t1, t2, U);

  auto H_matrix =
      compute_matrix_elements<arma::sp_cx_fmat>(basis, model.hamiltonian());

  std::vector<arma::sp_cx_fmat> density_operator_matrix;
  for (size_t i = 0; i < L; i++) {
    density_operator_matrix.push_back(compute_matrix_elements<arma::sp_cx_fmat>(
        basis, model.density_operator(i)));
  }

  auto state_vector = compute_vector_elements_serial<arma::cx_fvec>(
      basis, model.cls(L / 2 + 1) * model.cls(L / 2 + 1));
  state_vector /= arma::norm(state_vector);

  float initial_time = 0.0f;
  float final_time = 100.0f / t2;
  size_t num_time_steps = 2500;
  TimeIntegrator integrator(initial_time, final_time, num_time_steps, H_matrix);

  while (integrator.step(state_vector)) {
    for (size_t i = 0; i < L; i++) {
      std::cout << arma::cdot(state_vector,
                              density_operator_matrix[i] * state_vector)
                       .real();
      if (i < L - 1) std::cout << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
