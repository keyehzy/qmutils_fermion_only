#include <armadillo>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "qmutils/basis.h"
#include "qmutils/expression.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/operator.h"

using namespace qmutils;

class HeisenbergModel1D {
 public:
  explicit HeisenbergModel1D(size_t sites, float J) : m_sites(sites), m_J(J) {}

  size_t sites() const { return m_sites; }
  float J() const { return m_J; }

  Expression hamiltonian() const {
    Expression result;
    for (uint8_t i = 0; i < m_sites; ++i) {
      result += m_J * Expression::Spin::dot_product(i, (i + 1) % m_sites);
    }
    return result;
  }

 private:
  size_t m_sites;
  float m_J;
};

class VariationalWavefunction {
 public:
  VariationalWavefunction(size_t sites) : m_sites(sites), m_params(sites) {}

  float evaluate(const std::vector<int>& configuration) const {
    float result = 0.0f;
    for (size_t i = 0; i < m_sites; ++i) {
      result += m_params[i] * configuration[i];
    }
    return std::exp(result);
  }

  const arma::fvec& get_params() const { return m_params; }

  size_t num_params() const { return m_sites; }

  void update_params(const arma::fvec& new_params) { m_params = new_params; }

  void print_wavefunction() const {
    std::cout << "Variational Wavefunction Parameters:\n";
    for (size_t i = 0; i < m_sites; ++i) {
      std::cout << "Site " << i << ": " << m_params[i] << "\n";
    }
    std::cout << "\nWavefunction form: exp(";
    for (size_t i = 0; i < m_sites; ++i) {
      if (i > 0) std::cout << " + ";
      std::cout << m_params[i] << " * s_" << i;
    }
    std::cout << ")\n";
    std::cout << "Where s_i is the spin at site i (+1 or -1)\n";
  }

 private:
  size_t m_sites;
  arma::fvec m_params;
};

class VMCSolver {
 public:
  VMCSolver(const HeisenbergModel1D& model,
            const VariationalWavefunction& wavefunction, size_t num_samples)
      : m_model(model),
        m_wavefunction(wavefunction),
        m_num_samples(num_samples),
        m_rng(std::random_device{}()) {}

  float compute_energy() {
    float energy_sum = 0.0f;
    float weight_sum = 0.0f;

    for (size_t i = 0; i < m_num_samples; ++i) {
      auto config = generate_random_configuration();
      float weight = m_wavefunction.evaluate(config);
      weight_sum += weight * weight;

      float local_energy = compute_local_energy(config);
      energy_sum += local_energy * weight * weight;
    }

    return energy_sum / weight_sum;
  }

  arma::fvec compute_energy_gradient() {
    arma::fvec gradient(m_wavefunction.num_params(), arma::fill::zeros);
    float energy_sum = 0.0f;
    float weight_sum = 0.0f;

    for (size_t i = 0; i < m_num_samples; ++i) {
      auto config = generate_random_configuration();
      float weight = m_wavefunction.evaluate(config);
      weight_sum += weight * weight;

      float local_energy = compute_local_energy(config);
      energy_sum += local_energy * weight * weight;

      // Compute gradient terms
      for (size_t j = 0; j < m_wavefunction.num_params(); ++j) {
        gradient(j) += 2.0f * config[j] *
                       (local_energy - energy_sum / weight_sum) * weight *
                       weight;
      }
    }

    return gradient / weight_sum;
  }

 private:
  std::vector<int> generate_random_configuration() {
    std::vector<int> config(m_model.sites());
    std::uniform_int_distribution<> dist(0, 1);
    for (auto& spin : config) {
      spin = 2 * dist(m_rng) - 1;  // Map 0,1 to -1,1
    }
    return config;
  }

  float compute_local_energy(const std::vector<int>& config) {
    float energy = 0.0f;
    for (size_t i = 0; i < m_model.sites(); ++i) {
      size_t j = (i + 1) % m_model.sites();
      energy += 0.5f * m_model.J() * (config[i] * config[j]);
    }
    return energy;
  }

  const HeisenbergModel1D& m_model;
  const VariationalWavefunction& m_wavefunction;
  size_t m_num_samples;
  std::mt19937 m_rng;
};

static void print_state(const arma::cx_fvec& eigvec, const Basis& basis,
                        size_t count = 10) {
  std::vector<Term> terms;
  for (size_t i = 0; i < basis.size(); i++) {
    terms.emplace_back(std::norm(eigvec[i]), basis.at(i).operators());
  }
  std::sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
    return std::norm(a.coefficient()) > std::norm(b.coefficient());
  });
  for (size_t i = 0; i < std::min(count, terms.size()); i++) {
    std::cout << terms[i].to_string() << "\n";
  }
}

int main() {
  const size_t num_sites = 8;
  const float J = 1.0f;

  HeisenbergModel1D model(num_sites, J);

  // Exact diagonalization
  {
    // Construct the basis
    FermionicBasis basis(num_sites, num_sites, 0);
    std::cout << "Basis size: " << basis.size() << std::endl;

    // Compute Hamiltonian matrix
    auto H_matrix =
        compute_matrix_elements<arma::sp_cx_fmat>(basis, model.hamiltonian());

    // Find eigenvalues and eigenvectors
    arma::cx_fvec eigenvalues;
    arma::cx_fmat eigenvectors;
    arma::eigs_gen(eigenvalues, eigenvectors, H_matrix, 1, "sr");

    std::cout << "Exact ground state energy: " << eigenvalues[0] << std::endl;
    std::cout << "Exact ground state: " << std::endl;

    print_state(eigenvectors.col(0), basis);
  }

  // Variatonal Monte Carlo
  {
    const size_t num_samples = 10000;
    const size_t num_iterations = 1000;
    const float learning_rate = 0.01f;

    VariationalWavefunction wavefunction(num_sites);
    VMCSolver solver(model, wavefunction, num_samples);

    for (size_t iter = 0; iter < num_iterations; ++iter) {
      float energy = solver.compute_energy();
      arma::fvec gradient = solver.compute_energy_gradient();

      // Update wavefunction parameters
      arma::fvec new_params =
          wavefunction.get_params() - learning_rate * gradient;
      wavefunction.update_params(new_params);

      std::cout << "Iteration " << iter << ": Energy = " << energy << std::endl;

      if (iter % 10 == 0) {
        wavefunction.print_wavefunction();
        std::cout << std::endl;
      }
    }
  }

  return 0;
}
