#include <armadillo>
#include <fstream>

#include "qmutils/assert.h"
#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/functional.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"

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
};

[[maybe_unused]] static std::vector<std::pair<float, float>> calculate_dos(
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

[[maybe_unused]] static std::vector<std::pair<float, float>>
calculate_integrated_dos(const std::vector<std::pair<float, float>>& dos) {
  std::vector<std::pair<float, float>> integrated_dos(dos.size());
  float integral = 0.0f;
  float dE = dos[1].first - dos[0].first;

  for (size_t i = 0; i < dos.size(); ++i) {
    integral += dos[i].second * dE;
    integrated_dos[i] = {dos[i].first, integral};
  }

  return integrated_dos;
}

static Expression commutator(const Expression& A, const Expression& B,
                             NormalOrderer& orderer) {
  return orderer.normal_order(A * B) - orderer.normal_order(B * A);
}

static constexpr uint64_t factorial(uint64_t n) {
  if (n == 0 || n == 1) return 1;
  return n * factorial(n - 1);
}

static Expression bakerCampbellHausdorff(const Expression& A,
                                         const Expression& B, float lambda,
                                         NormalOrderer& orderer,
                                         size_t order = 10) {
  std::vector<Expression> terms(order + 1);

  terms[0] = B;
  terms[1] = commutator(A, B, orderer);

  for (size_t n = 2; n <= order; ++n) {
    terms[n] = commutator(A, terms[n - 1], orderer);
  }

  Expression result;
  for (size_t n = 0; n <= order; ++n) {
    float coeff =
        powf(lambda, static_cast<float>(n)) / static_cast<float>(factorial(n));
    result += terms[n] * coeff;
  }

  return orderer.normal_order(result);
}

[[maybe_unused]] static float find_optimal_theta(
    const Expression& H, const Term& target, float min_theta, float max_theta,
    float tol, NormalOrderer& orderer, size_t level = 10) {
  auto calculate_coefficient = [&](float theta) -> std::complex<float> {
    Expression R;
    R -= target;
    R += target.adjoint();
    auto e = bakerCampbellHausdorff(R, H, theta, orderer, level);
    return e[target.operators()];
  };

  while (max_theta - min_theta > tol) {
    float mid1 = min_theta + (max_theta - min_theta) / 3.0f;
    float mid2 = max_theta - (max_theta - min_theta) / 3.0f;

    std::complex<float> coeff1 = calculate_coefficient(mid1);
    std::complex<float> coeff2 = calculate_coefficient(mid2);

    if (std::norm(coeff1) < std::norm(coeff2)) {
      max_theta = mid2;
    } else {
      min_theta = mid1;
    }
  }

  return (min_theta + max_theta) / 2.0f;
}

[[maybe_unused]] static float find_optimal_theta_quadratic(
    const Expression& H, const Term& target, NormalOrderer& orderer,
    float tol = 1e-3f) {
  Expression R;
  R -= target;
  R += target.adjoint();

  auto comm1 = commutator(R, H, orderer);
  auto comm2 = commutator(R, comm1, orderer);

  float alpha = comm1[target.operators()].real();
  float beta = comm2[target.operators()].real();
  float gamma = target.coefficient().real();

  float pi = std::numbers::pi_v<float>;

  if (std::abs(beta) < 1e-6f) {
    std::cout << "Oooops, we have to bail!" << std::endl;
    return find_optimal_theta(H, target, -pi, pi, tol, orderer, 4);
  }

  float alpha_beta = alpha / beta;
  float gamma_beta = gamma / beta;

  float Delta = alpha_beta * alpha_beta - 2.0f * gamma_beta;

  if (Delta < 0.0f) {
    std::cout << "Oooops, we found negative delta!" << std::endl;
    return find_optimal_theta(H, target, -pi, pi, tol, orderer, 4);
  }

  float r1 = -alpha_beta + std::sqrt(Delta);
  float r2 = -alpha_beta - std::sqrt(Delta);

  if (std::norm(r1) < std::norm(r2)) {
    return r1;
  } else {
    return r2;
  }
}

[[maybe_unused]] static Expression jacobi_diagonalization(
    const Expression& H, size_t nbody, bool& is_diagonal,
    std::vector<Expression>& rotate_ops, float tol = 1e-3f) {
  float off_diag_norm = 0.0f;
  Term largest_off_diag_term;
  float largest_off_diag = 0.0f;
  size_t biggest_term_so_far = 0;

  for (const auto& [ops, coeff] : H.terms()) {
    biggest_term_so_far = std::max(biggest_term_so_far, ops.size());
    if (ops.size() == 2 * nbody && !Term::is_diagonal(ops)) {
      off_diag_norm += std::norm(coeff);
      if (std::norm(coeff) > largest_off_diag) {
        largest_off_diag = std::norm(coeff);
        largest_off_diag_term = Term(coeff, ops);
      }
    }
  }

  std::cout << "Off-diagonal norm: " << off_diag_norm << "\n";
  std::cout << "Biggest term so far: " << biggest_term_so_far << "\n";

  if (std::norm(largest_off_diag_term.coefficient()) < tol * tol) {
    is_diagonal = true;
    return H;
  }
  std::cout << "term: " << largest_off_diag_term.to_string() << std::endl;

  // float pi = std::numbers::pi_v<float>;
  NormalOrderer orderer;
  // float optimal_theta =
  //   find_optimal_theta(H, largest_off_diag_term, -pi / 4.0f, pi / 4.0f, tol,
  //   orderer);
  float optimal_theta =
      find_optimal_theta_quadratic(H, largest_off_diag_term, orderer);

  std::cout << "optimal_theta: " << optimal_theta << "\n";

  Expression R;
  R -= largest_off_diag_term;
  R += largest_off_diag_term.adjoint();

  for (Expression& ops : rotate_ops) {
    ops = bakerCampbellHausdorff(R, ops, optimal_theta, orderer);
  }

  return bakerCampbellHausdorff(R, H, optimal_theta, orderer);
}

template <typename Index>
[[maybe_unused]] static void print_expression(const Expression& expr,
                                              const Index& index,
                                              float tol = 1e-3f) {
  std::vector<Term> terms;
  for (const auto& [ops, coeff] : expr.terms()) {
    terms.push_back(Term(coeff, ops));
  }
  std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
    if (a.operators().size() != b.operators().size()) {
      return a.operators().size() < b.operators().size();
    }
    return std::norm(a.coefficient()) > std::norm(b.coefficient());
  });
  for (const auto& term : terms) {
    if (std::norm(term.coefficient()) < tol * tol) continue;
    if (term.is_diagonal()) {
      std::cout << "Diagonal: " << index.to_string(term) << "\n";
    } else {
      std::cout << "Not Diagonal: " << index.to_string(term) << "\n";
    }
  }
}

// template <size_t L>
// static Expression transform_operator_to_AB(const Operator& op, const typename
// SawtoothModel<L>::Index& index) {
//   Expression result;

//   auto [cell, i] = index.from_orbital(op.orbital());

//   Operator A(op.type(), op.spin(), index.to_orbital(cell, 0),
//   op.statistics()); Operator B(op.type(), op.spin(), index.to_orbital(cell,
//   1), op.statistics()); Operator A_next(op.type(), op.spin(),
//   index.to_orbital((cell + 1) % L, 0), op.statistics()); Operator
//   B_next(op.type(), op.spin(), index.to_orbital((cell + 1) % L, 1),
//   op.statistics());

//   if (i == 0) {
//     result += Term(1.0f / std::sqrt(3.0f), {A});
//     result -= Term(1.0f / std::sqrt(2.0f), {A_next});
//     result -= Term(1.0f / std::sqrt(6.0f), {B_next});
//   } else {
//     result += Term(1.0f / std::sqrt(6.0f), {A});
//     result += Term(1.0f / std::sqrt(2.0f), {B});
//     result += Term(1.0f / std::sqrt(3.0f), {A_next});
//   }

//   return result;
// }

int main() {
  const size_t L = 2;
  const float t2 = 1.0f;
  const float t1 = t2 * std::sqrt(2.0f);
  BosonicBasis basis(2 * L, 2);  // 2 sites per unit cell, P particles
  std::cout << "# Basis size: " << basis.size() << std::endl;

#if 1
  {
    SawtoothModel<L> model(t1, t2, 2.0f);

    auto H_matrix =
        compute_matrix_elements<arma::cx_fmat>(basis, model.hamiltonian());

    std::cout << "# Matrix elements computed" << std::endl;

    QMUTILS_ASSERT(
        arma::approx_equal(H_matrix, H_matrix.t(), "absdiff", 1e-4f));

    arma::fvec eigenvalues;
    arma::cx_fmat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, H_matrix);

    std::cout << "# Eigenvalues computed" << std::endl;
    std::cout << eigenvalues << std::endl;

    // auto dos = calculate_dos(eigenvalues, 0.1f);
    // auto integrated_dos = calculate_integrated_dos(dos);

    // std::cout << "# DOS computed" << std::endl;

    // std::ofstream dos_file("dos.dat");

    // for (size_t i = 0; i < dos.size(); ++i) {
    //   dos_file << dos[i].first << " "
    //            << dos[i].second / static_cast<float>(basis.size()) << " "
    //            << integrated_dos[i].second / static_cast<float>(basis.size())
    //            << "\n";
    // }
  }
#endif

#if 0
  {
    SawtoothModel<L> foo(1.0f, 2.0f, 3.0f);

    std::vector<Term> cterms;

    for (size_t i = 0; i < foo.index().total_size(); i++) {
      for (size_t j = 0; j < foo.index().total_size(); j++) {
        for (size_t k = 0; k < foo.index().total_size(); k++) {
          cterms.emplace_back(1.0f, std::vector<Operator> {
              Operator::Boson::creation(Operator::Spin::Up, i),
              Operator::Boson::creation(Operator::Spin::Up, j),
              Operator::Boson::annihilation(Operator::Spin::Up, k)
            });
        }
      }
    }

    std::ofstream u_file("u.dat");
    std::ofstream v_file("v.dat");

    std::unordered_map<size_t, std::vector<float>> ddd;

    std::vector<float> Us;
    for (float U = 0.0f; U < 8.01f; U += 0.5f) {
      Us.push_back(U);
    }

    for (float U : Us) {
      SawtoothModel<L> model(t1, t2, U);
      QMUTILS_ASSERT(model.hamiltonian().is_purely(Operator::Statistics::Bosonic));
      Expression Hp = model.hamiltonian();

      std::vector<Expression> rotate_ops {
        Expression(Operator::Boson::creation(Operator::Spin::Up, 0))
      };

      size_t jacobi_iters = 0;
      for (size_t nbody = 1; nbody < 3; nbody++) {
        std::cout << "Nbody: " << nbody << std::endl;
        bool is_diagonal = false;

        while (!is_diagonal) {
          std::cout << "Jacobi Iterations: " << jacobi_iters++ << "\n";
          Hp = jacobi_diagonalization(Hp, nbody, is_diagonal, rotate_ops);
        }
        std::cout << "Diagonal at level: " << nbody << std::endl;
      }

      u_file << "\"U=" << U << "\"" << std::endl;
      for (size_t i = 0; i < cterms.size(); i++) {
        u_file << i << " " << rotate_ops[0][cterms[i].operators()].real() << std::endl;
        ddd[i].push_back(rotate_ops[0][cterms[i].operators()].real());
      }
      u_file << "\n\n";
    }

    for (const auto& [i, z] : ddd) {
      v_file << "\"term:" << foo.index().to_string(cterms[i]) << "\"" << std::endl;;
      for (size_t ui = 0; ui < z.size(); ui++) {
        v_file << Us[ui] << " " << z[ui] << std::endl;
      }
      v_file << "\n\n";
    }
  }
#endif

  return 0;
}
