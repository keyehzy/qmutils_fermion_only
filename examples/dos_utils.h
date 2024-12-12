#pragma once

#include <armadillo>

class Dos_Utils {
 public:
  Dos_Utils(const arma::fvec& eigenvalues, float sigma = 0.1f,
            size_t num_points = 1000, float padding_factor = 0.1f) {
    float E_min = eigenvalues.min();
    float E_max = eigenvalues.max();
    float padding = (E_max - E_min) * padding_factor;
    E_min -= padding;
    E_max += padding;

    float dE = (E_max - E_min) / static_cast<float>(num_points - 1);

    m_dos.resize(num_points);
    m_integrated_dos.resize(num_points);

    // Calculate DOS using Gaussian broadening
    const float normalization =
        1.0f / (sigma * std::sqrt(2.0f * std::numbers::pi_v<float>));

    float integral = 0.0f;

#pragma omp parallel for
    for (size_t i = 0; i < num_points; ++i) {
      float E = E_min + static_cast<float>(i) * dE;
      float rho = 0.0f;

      for (size_t j = 0; j < eigenvalues.n_elem; ++j) {
        float delta_E = (E - eigenvalues(j)) / sigma;
        rho += std::exp(-0.5f * delta_E * delta_E);
      }

      m_dos[i] = {E, rho * normalization};
      integral += rho * normalization * dE;
      m_integrated_dos[i] = {E, integral};
    }
  }

  const std::vector<std::pair<float, float>>& dos() const { return m_dos; }
  const std::vector<std::pair<float, float>>& integrated_dos() const {
    return m_integrated_dos;
  }

 private:
  std::vector<std::pair<float, float>> m_dos;
  std::vector<std::pair<float, float>> m_integrated_dos;
};
