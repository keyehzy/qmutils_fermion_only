#include "qmutils/expression.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/sparse_matrix.h"
#include "qmutils/term.h"
#include "qmutils/index.h"

using namespace qmutils;

template <size_t Lx, size_t Ly>
class SquareLatticeLandauLevel {
 public:
    using Index = StaticIndex<Lx, Ly>;

  SquareLatticeLandauLevel(float t, float phi)
      : m_t(t), m_phi(phi), m_index() {
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
  Operator::Spin spin = Operator::Spin::Up; // Spin is not relevant here
  for (size_t i = 0; i < Lx; ++i) {
    std::complex<float> phase(0.0f, 2.0f * std::numbers::pi_v<float> * m_phi * i / Lx);
    for (size_t j = 0; i < Lx; ++i) {    
        size_t center = m_index.to_orbital(i, j);
        size_t right = m_index.to_orbital((i+1)%Lx, j);
        size_t up = m_index.to_orbital(i, (j+1)%Ly);
        
        Term peirls_transformed = std::exp(phase) * Term::Fermion::one_body(spin, up, spin, center);

        m_hamiltonian += m_t * Expression::Fermion::hopping(right, center, spin);
        m_hamiltonian += m_t * peirls_transformed;
        m_hamiltonian += m_t * peirls_transformed.adjoint();        
    }
  }
}
};


int main() {
  const size_t Lx = 2;
  const size_t Ly = 2;
  const size_t particles = 1;

  SquareLatticeLandauLevel<Lx, Ly> model(1.0f, 0.0f);
  Basis basis(Lx * Ly, particles);

  for (const auto&[ops, coeff] : model.hamiltonian().terms()) {
    term_to_string(Term(coeff, ops), model.index());
  }

  return 0;
}
