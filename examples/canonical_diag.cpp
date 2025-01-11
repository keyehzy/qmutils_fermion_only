#include <armadillo>
#include <cmath>
#include <random>
#include <unordered_set>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
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
    float coeff =
        powf(lambda, static_cast<float>(n)) / static_cast<float>(factorial(n));
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
                                                           const Term& target,
                                                           float tol = 1e-6f) {
  Expression R;
  R -= target;
  R += target.adjoint();

  auto comm1 = commutator(R, H, orderer);
  auto comm2 = commutator(R, comm1, orderer);

  float alpha = comm1[target.operators()].real();
  float beta = comm2[target.operators()].real();
  float gamma = target.coefficient().real();

  if (std::norm(beta) < tol * tol) {
    return std::numbers::pi_v<float> / 4.0f;
  }

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

[[maybe_unused]] static float find_optimal_theta_newton(const Expression& H,
                                                        Term& target,
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
    auto e = bakerCampbellHausdorff(R, A, theta, orderer, /*level=*/4);
    return e[target.operators()];
  };

  float pi = std::numbers::pi_v<float>;

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
    std::cout << "initial guess: " << prev << " newton step: " << result
              << "\n";
  } while (std::abs(prev - result) > tol);

  return result;
}

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

  if (std::abs(largest_off_diag_term.coefficient()) < tol) {
    is_diagonal = true;
    return H;
  }
  std::cout << "term: " << largest_off_diag_term.to_string() << std::endl;

  float optimal_theta = find_optimal_theta_quadratic(H, largest_off_diag_term);
  // float pi = std::numbers::pi_v<float>;
  // float optimal_theta =
  //     find_optimal_theta(H, largest_off_diag_term, -pi / 4.0f, pi / 4.0f,
  //     tol);
  //  float optimal_theta =
  //      find_optimal_theta_newton(H, largest_off_diag_term, tol);

  std::cout << "optimal_theta: " << optimal_theta << "\n";

  Expression R;
  R -= largest_off_diag_term;
  R += largest_off_diag_term.adjoint();

  for (Expression& ops : rotate_ops) {
    ops = bakerCampbellHausdorff(R, ops, optimal_theta, orderer);
  }

  return bakerCampbellHausdorff(R, H, optimal_theta, orderer);
}

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
    if (std::abs(term.coefficient()) < 1e-4f) continue;
    if (term.is_diagonal()) {
      std::cout << "Diagonal: " << term.to_string() << "\n";
    } else {
      std::cout << "Not Diagonal: " << term.to_string() << "\n";
    }
  }
}

int main() {
  float t = 1.0f;
  float U = 4.0f;
  Expression hamiltonian;

  const size_t L = 2;

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
    FermionicBasis basis(L, 1);
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
      Hp, Expression(Operator::Fermion::creation(Operator::Spin::Up, 0)),
      orderer);

  print_expression(dt);

  orderer.print_cache_stats();
  return 0;
}
