#include <armadillo>
#include <cmath>
#include <random>
#include <unordered_set>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/functional.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/normal_order.h"

using namespace qmutils;

static NormalOrderer orderer;
// static float pi = std::numbers::pi_v<float>;

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
    float coeff = powf(lambda, n) / static_cast<float>(factorial(n));
    result += terms[n] * coeff;
  }

  return orderer.normal_order(result);
}

[[maybe_unused]] static float find_optimal_theta(const Expression& H,
                                                 const Term& target,
                                                 float min_theta,
                                                 float max_theta, float tol,
                                                 size_t level = 10) {
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

    if (std::abs(coeff1) < std::abs(coeff2)) {
      max_theta = mid2;
    } else {
      min_theta = mid1;
    }
  }

  return (min_theta + max_theta) / 2.0f;
}

[[maybe_unused]] static float find_optimal_theta_quadratic(const Expression& H,
                                                           const Term& target) {
  Expression R;
  R -= target;
  R += target.adjoint();

  auto comm1 = commutator(R, H, orderer);
  auto comm2 = commutator(R, comm1, orderer);

  float alpha = comm1[target.operators()].real();
  float beta = comm2[target.operators()].real();
  float gamma = target.coefficient().real();

  float alpha_beta = alpha / beta;
  float gamma_beta = gamma / beta;

  float r1 =
      -alpha_beta + std::sqrt(alpha_beta * alpha_beta - 2.0f * gamma_beta);
  float r2 =
      -alpha_beta - std::sqrt(alpha_beta * alpha_beta - 2.0f * gamma_beta);

  if (std::abs(r1) < std::abs(r2)) {
    return r1;
  } else {
    return r2;
  }
}

#if 0
static float find_optimal_theta_newton(Expression& H, Term& target,
                                float tol) {
  auto calculate_coefficient = [&](float theta) -> std::complex<float> {
    Expression R;
    R -= target;
    R += target.adjoint();
    auto e = bakerCampbellHausdorff(R, H, theta, orderer, /*level=*/4);
    return e[target.operators()];
  };

  auto calculate_coefficient_prime = [&](float theta) -> std::complex<float> {
    Expression R;
    R -= target;
    R += target.adjoint();
    Expression A = commutator(R, H, orderer);
    auto e =  bakerCampbellHausdorff(R, A, theta, orderer, /*level=*/4);
    return e[target.operators()];
  };

  // Perform a sweep to find a good first guess
  float result =
    find_optimal_theta(H, target, -pi / 4.0f, pi / 4.0f, 0.1f, /*level=*/1);
  float prev;
  float damp_factor = 0.1f;

  // Perform newton step once
  do {
    prev = result;
    std::complex<float> f = calculate_coefficient(prev);
    std::complex<float> df = calculate_coefficient_prime(prev);
    result = prev - (f / df).real();
    if (std::abs(result - prev) > damp_factor) {
      result = prev + damp_factor * std::copysign(1.0f, result - prev);
    }
    result = std::clamp(result, -pi / 4.0f, pi / 4.0f);
    std::cout << "initial guess: " << prev << " newton step: " << result << "\n";
  } while (std::abs(prev - result) > tol);

  return result;
}
#endif

#if 0
static float find_optimal_theta_gd(const Expression& H, const Term& target, size_t nbody, float tol,
                                   size_t level = 10, float eta = 0.01f, float beta = 0.9f) {
  auto compute_derivative = [&](float theta) -> Expression {
    Expression R;
    R -= target;
    R += target.adjoint();
    return bakerCampbellHausdorff(R, H, theta, orderer, level);
  };

  float theta = find_optimal_theta(H, target, -pi/4.0f, pi/4.0f, 0.1f, /*level=*/1);
  float prev_theta = theta;
  float velocity = 0.0f;
  float prev_off_diag_norm = std::numeric_limits<float>::max();

  size_t max_iters = 10;

  for (size_t iter = 0; iter < max_iters; iter++) {
    Expression df = compute_derivative(theta);
    float current_off_diag_norm = 0.0f;

    // Calculate current off-diagonal norm
    for (const auto& [ops, coeff] : df.terms()) {
      Term term(coeff, ops);
      if (term.size() == 2 * nbody) {
        if (!term.is_diagonal()) {
          current_off_diag_norm += std::norm(coeff);
        }
      }
    }

    // Check if off-diagonal norm is increasing
    if (current_off_diag_norm > prev_off_diag_norm) {
      // Return the previous theta since it gave a better result
      theta = prev_theta;
      break;
    }

    // Store current state before update
    prev_theta = theta;
    prev_off_diag_norm = current_off_diag_norm;

    // Update with momentum
    velocity = beta * velocity + eta * current_off_diag_norm;
    theta -= velocity;

    // Optional: Clamp theta to desired range
    // theta = std::clamp(theta, -pi / 4.0f, pi / 4.0f);

    std::cout << "Off-diagonal norm: " << current_off_diag_norm
              << " Theta: " << theta
              << " Previous theta: " << prev_theta
              << " Velocity: " << velocity << "\n";

    if (iter == max_iters-1) {
      std::cout << "Reached max iterations" << "\n";
    }
  }

  return theta;
}
#endif

static Expression jacobi_diagonalization(const Expression& H, size_t nbody,
                                         bool& is_diagonal,
                                         std::vector<Expression>& rotate_ops,
                                         float tol = 1e-3f) {
  float off_diag_norm = 0.0f;
  Term largest_off_diag_term;
  float largest_off_diag = 0.0f;

  for (const auto& [ops, coeff] : H.terms()) {
    if (ops.size() == 2 * nbody && !Term::is_diagonal(ops)) {
      off_diag_norm += std::norm(coeff);
      if (std::norm(coeff) > largest_off_diag) {
        largest_off_diag = std::norm(coeff);
        largest_off_diag_term = Term(coeff, ops);
      }
    }
  }

  std::cout << "Off-diagonal norm: " << off_diag_norm << "\n";

  if (std::norm(largest_off_diag_term.coefficient()) < tol) {
    is_diagonal = true;
    return H;
  }
  std::cout << "term: " << largest_off_diag_term.to_string() << std::endl;

  float optimal_theta = find_optimal_theta_quadratic(H, largest_off_diag_term);
  // float optimal_theta = find_optimal_theta(H, largest_off_diag_term,
  // -M_PI/4.0f, M_PI/4.0f, tol);

  std::cout << "optimal_theta: " << optimal_theta << "\n";

  Expression R;
  R -= largest_off_diag_term;
  R += largest_off_diag_term.adjoint();

  for (Expression& ops : rotate_ops) {
    ops = bakerCampbellHausdorff(R, ops, optimal_theta, orderer);
  }

  return bakerCampbellHausdorff(R, H, optimal_theta, orderer);
}

#if 1
static void print_expression(const Expression& expr) {
  std::vector<Term> terms;
  for (const auto& [ops, coeff] : expr.terms()) {
    terms.push_back(Term(coeff, ops));
  }
  std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
    if (a.operators().size() != b.operators().size()) {
      return a.operators().size() < b.operators().size();
    }
    return std::abs(a.coefficient()) > std::abs(b.coefficient());
  });
  for (const auto& term : terms) {
    if (term.is_diagonal()) {
      std::cout << "Diagonal: " << term.to_string() << "\n";
    } else {
      std::cout << "Not Diagonal: " << term.to_string() << "\n";
    }
  }
}
#endif

#if 0
template <typename Index>
static void print_expression(const Expression& expr, const Index& index) {
  std::vector<Term> terms;
  for (const auto& [ops, coeff] : expr.terms()) {
    terms.push_back(Term(coeff, ops));
  }
  std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
    if (a.operators().size() != b.operators().size()) {
      return a.operators().size() < b.operators().size();
    }
    return std::abs(a.coefficient()) > std::abs(b.coefficient());
  });
  for (const auto& term : terms) {
    if (term.is_diagonal()) {
      std::cout << "Diagonal: " << index.to_string(term) << "\n";
    } else {
      std::cout << "Not Diagonal: " << index.to_string(term) << "\n";
    }
  }
}
#endif

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

template <typename Index>
static Expression transform_operator_n_p_basis(const Operator& op,
                                               const Index& index) {
  Expression result;
  const float type_sign =
      (op.type() == Operator::Type::Annihilation) ? -1.0f : 1.0f;

  auto [cell, site] = index.from_orbital(op.orbital());

  Operator n(op.type(), op.spin(), index.to_orbital(cell, 0), op.statistics());
  Operator p(op.type(), op.spin(), index.to_orbital(cell, 1), op.statistics());
  if (site == 0) {
    result += Term(std::complex<float>(0, type_sign / std::sqrt(2.0f)), {n});
    result += Term(std::complex<float>(1.0f / std::sqrt(2.0f), 0), {p});
  } else {
    result += Term(std::complex<float>(1.0f / std::sqrt(2.0f), 0), {n});
    result += Term(std::complex<float>(0, type_sign / std::sqrt(2.0f)), {p});
  }

  return result;
}

#if 0
int main() {
  const size_t L = 3;
  std::vector<float> t = {1.0f, 1.0f, 1.0f};
  float U = 20.0f;

  Expression hamiltonian;
  
  for (size_t i = 0; i < L; i++) {
    hamiltonian += -t[i] * Expression::Fermion::hopping(i, (i+1)%L, Operator::Spin::Up);
  }
  hamiltonian += U * Term::Fermion::density(Operator::Spin::Up, L-1);

  // for (const auto& [ops, coeff] : hamiltonian.terms()) {
  //   std::cout << Term(coeff, ops).to_string() << std::endl;
  // }

  auto off_diag_number = [](const Expression& H)
      -> float {
    float off_diag = 0.0f;
    for (const auto& [ops, coeff] : H.terms()) {
      Term term(coeff, ops);
      if (!term.is_diagonal()) {
        off_diag += std::norm(coeff);
      }
    }
    return off_diag;
  };

  auto calculate_off_diag = [](const Expression& H, const Term&
  target, float theta, size_t level = 10) -> float {
    Expression R;
    R -= target;
    R += target.adjoint();
    auto e = bakerCampbellHausdorff(R, H, theta, orderer, level);

    float off_diag = 0.0f;
    for (const auto& [ops, coeff] : e.terms()) {
    Term term(coeff, ops);
    if (!term.is_diagonal()) {
      off_diag += std::norm(coeff);
    }
  }
  return off_diag;
  };

  auto compute_derivative = [](const Expression& H, const Term& target,
                               float theta, size_t level = 10) -> float {
    Expression R;
    R -= target;
    R += target.adjoint();
    Expression A = H * commutator(R, H, orderer) + commutator(R, H, orderer) * H;
    auto e = bakerCampbellHausdorff(R, A, theta, orderer, level);

    float off_diag = 0.0f;
    for (const auto& [ops, coeff] : e.terms()) {
      Term term(coeff, ops);
      if (!term.is_diagonal()) {
        off_diag += std::real(coeff);
      }
    }
    return off_diag;
  };

  Term target1 = Term(1.0f, {Operator::Fermion::creation(Operator::Spin::Up, 0),
                             Operator::Fermion::annihilation(Operator::Spin::Up, 1)});
  Term target2 = Term(1.0f, {Operator::Fermion::creation(Operator::Spin::Up, 1),
                             Operator::Fermion::annihilation(Operator::Spin::Up, 2)});
  Term target3 = Term(1.0f, {Operator::Fermion::creation(Operator::Spin::Up, 2),
                             Operator::Fermion::annihilation(Operator::Spin::Up, 0)});

  float theta1 = 0.1f;
  float theta2 = 0.2f;
  float theta3 = 0.3f;

  float off_diag = off_diag_number(hamiltonian);
  std::cout << "Off-diagonal norm: " << off_diag << "\n";

  for (size_t iter = 0; iter < 10; iter++) {
    float deriv1 = compute_derivative(hamiltonian, target1, theta1);

    theta1 -= 0.0001f * deriv1;

    Expression R1;
    R1 -= target1;
    R1 += target1.adjoint();

    hamiltonian = bakerCampbellHausdorff(R1, hamiltonian, theta1, orderer);
    
    float deriv2 = compute_derivative(hamiltonian, target2, theta2);
    theta2 -= 0.0001f * deriv2;
    Expression R2;
    R2 -= target2;
    R2 += target2.adjoint();
    hamiltonian = bakerCampbellHausdorff(R2, hamiltonian, theta2, orderer);
    
    float deriv3 = compute_derivative(hamiltonian, target3, theta3);
    theta3 -= 0.0001f * deriv3;
    Expression R3;
    R3 -= target3;
    R3 += target3.adjoint();
    hamiltonian = bakerCampbellHausdorff(R3, hamiltonian, theta3, orderer);

    float new_off_diag = off_diag_number(hamiltonian);
    std::cout << "Off-diagonal norm: " << new_off_diag << "\n";
    off_diag = new_off_diag;
  }

}
#endif

#if 1
int main() {
  // const size_t L = 2;
  // const float J = 1.0f;
  // const float U = 4.0f * J;

  // CreutzLadderModel<L> model(J, 0.5f * std::numbers::pi_v<float>, U);

  float t = 1.0f;
  float U = 20.0f;
  Expression hamiltonian;

  const size_t L = 3;

  for (size_t i = 0; i < L - 1; i++) {
    hamiltonian +=
        -t * Expression::Fermion::hopping(i, (i + 1) % L, Operator::Spin::Up);
    hamiltonian +=
        -t * Expression::Fermion::hopping(i, (i + 1) % L, Operator::Spin::Down);
  }

  for (size_t i = 0; i < L; i++) {
    hamiltonian += U * Term::Fermion::density_density(Operator::Spin::Up, i,
                                                      Operator::Spin::Down, i);
  }

  {
    // Diagonalize hamiltonian with two particles
    FermionicBasis basis(L, 2);
    auto H_matrix = compute_matrix_elements<arma::cx_fmat>(basis, hamiltonian);

    arma::fvec eigenvalues;
    arma::cx_fmat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, H_matrix);

    std::cout << "Eigenvalues:" << std::endl;
    for (size_t i = 0; i < eigenvalues.n_elem; ++i) {
      std::cout << basis.at(i).to_string() << "  " << eigenvalues(i)
                << std::endl;
    }
  }

  Expression Hp = hamiltonian;

  std::cout << "Hamiltonian in real space:" << std::endl;
  for (const auto& [ops, coeff] : hamiltonian.terms()) {
    Term term(coeff, ops);
    std::cout << term.to_string() << std::endl;
  }

  // std::cout << "Hamiltonian in real space:" << std::endl;
  // for (const auto& [ops, coeff] : model.hamiltonian().terms()) {
  //   Term term(coeff, ops);
  //   std::cout << model.index().to_string(term) << std::endl;
  // }

  // Expression Hp = transform_expression(
  //    transform_operator_n_p_basis<typename CreutzLadderModel<L>::Index>,
  //    model.hamiltonian(), model.index());

  // std::cout << "Hamiltonian in n-p basis:" << std::endl;
  // for (const auto& [ops, coeff] : Hp.terms()) {
  //   Term term(coeff, ops);
  //   std::cout << model.index().to_string(term) << std::endl;
  // }

  // std::vector<Expression> rotate_ops = {
  //   Expression(Operator::Boson::creation(Operator::Spin::Up,
  //   model.index().to_orbital(0, 0))), // n
  //   Expression(Operator::Boson::creation(Operator::Spin::Up,
  //   model.index().to_orbital(0, 1))), // p
  // };

  std::vector<Expression> rotate_ops = {
      Expression(Operator::Fermion::creation(Operator::Spin::Up, 0)),
      Expression(Term::Fermion::density(Operator::Spin::Up, 0)),
  };

  size_t jacobi_iters = 0;
  for (size_t nbody = 1; nbody < 10; nbody++) {
    std::cout << "Nbody: " << nbody << std::endl;
    bool is_diagonal = false;

    while (!is_diagonal) {
      std::cout << "Jacobi Iterations: " << jacobi_iters++ << "\n";
      Hp = jacobi_diagonalization(Hp, nbody, is_diagonal, rotate_ops);
    }
    std::cout << "Diagonal at level: " << nbody << std::endl;
  }

  std::cout << "Hamiltonian final form:" << std::endl;
  print_expression(Hp);
  std::cout << "\n\n";

  for (const auto& op : rotate_ops) {
    std::cout << "Rotate operator:" << std::endl;
    print_expression(op);
    std::cout << "\n\n";
  }

  std::cout << "Time derivative" << std::endl;

  Expression dt = commutator(
      Hp, Expression(Term::Fermion::density(Operator::Spin::Up, 0)), orderer);

  print_expression(dt);

  orderer.print_cache_stats();
  return 0;
}
#endif