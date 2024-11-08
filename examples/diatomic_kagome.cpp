#include <armadillo>
#include <fstream>
#include <vector>

#include "qmutils/assert.h"
#include "qmutils/basis.h"
#include "qmutils/eigsys.h"
#include "qmutils/expression.h"
#include "qmutils/functional.h"
#include "qmutils/index.h"
#include "qmutils/matrix_elements.h"
#include "qmutils/operator.h"

using namespace qmutils;

template <size_t Lx, size_t Ly>
class DiatomicKagome {
 public:
  using Index = StaticIndex<Lx, Ly, 6>;

  DiatomicKagome(float t1, float t2, float t3)
      : m_t1(t1), m_t2(t2), m_t3(t3), m_index() {
    construct_hamiltonian();
  }

  const Expression& hamiltonian() const { return m_hamiltonian; }
  const Index& index() const { return m_index; }

 private:
  float m_t1;
  float m_t2;
  float m_t3;
  Index m_index;
  Expression m_hamiltonian;

  void construct_hamiltonian() {
    auto intra_hopping = [&](size_t i, size_t j, size_t k1, size_t k2) {
      return Expression::Fermion::hopping(m_index.to_orbital(i, j, k1),
                                          m_index.to_orbital(i, j, k2),
                                          Operator::Spin::Up);
    };

    auto inter_hopping_x = [&](size_t i, size_t j, size_t k1, size_t k2) {
      return Expression::Fermion::hopping(
          m_index.to_orbital(i, j, k1), m_index.to_orbital((i + 1) % Lx, j, k2),
          Operator::Spin::Up);
    };

    auto inter_hopping_y = [&](size_t i, size_t j, size_t k1, size_t k2) {
      return Expression::Fermion::hopping(
          m_index.to_orbital(i, j, k1), m_index.to_orbital(i, (j + 1) % Ly, k2),
          Operator::Spin::Up);
    };

    for (size_t i = 0; i < Lx; ++i) {
      for (size_t j = 0; j < Ly; ++j) {
        // Intra-cell hopping
        // m_t1 hopping
        m_hamiltonian += m_t1 * intra_hopping(i, j, 2, 3);

        // m_t2 hopping
        m_hamiltonian += m_t2 * intra_hopping(i, j, 0, 1);
        m_hamiltonian += m_t2 * intra_hopping(i, j, 0, 2);
        m_hamiltonian += m_t2 * intra_hopping(i, j, 1, 2);

        m_hamiltonian += m_t2 * intra_hopping(i, j, 3, 4);
        m_hamiltonian += m_t2 * intra_hopping(i, j, 3, 5);
        m_hamiltonian += m_t2 * intra_hopping(i, j, 4, 5);

        // m_t3 hopping
        m_hamiltonian += m_t3 * intra_hopping(i, j, 2, 4);
        m_hamiltonian += m_t3 * intra_hopping(i, j, 2, 5);

        m_hamiltonian += m_t3 * intra_hopping(i, j, 3, 0);
        m_hamiltonian += m_t3 * intra_hopping(i, j, 3, 1);

        // Inter-cell hopping
        // m_t1 hopping
        m_hamiltonian += m_t1 * inter_hopping_y(i, j, 0, 5);
        m_hamiltonian += m_t1 * inter_hopping_x(i, j, 4, 1);

        // m_t2 hopping
        // Nothing

        // m_t3 hopping
        m_hamiltonian += m_t3 * inter_hopping_y(i, j, 0, 3);
        m_hamiltonian += m_t3 * inter_hopping_y(i, j, 0, 4);

        m_hamiltonian += m_t3 * inter_hopping_x(i, j, 4, 0);
        m_hamiltonian += m_t3 * inter_hopping_x(i, j, 4, 2);
      }
    }
  }
};

// struct Vec2 {
//   float x;
//   float y;

//   friend Vec2 operator+(const Vec2& a, const Vec2& b) {
//     return {a.x + b.x, a.y + b.y};
//   }

//   friend Vec2 operator-(const Vec2& a, const Vec2& b) {
//     return {a.x - b.x, a.y - b.y};
//   }
// };

enum class FourierTransformDirection : int { Forward = 1, Inverse = -1 };

template <size_t Lx, size_t Ly>
static Expression fourier_transform_operator_2d(
    const Operator& op, const typename DiatomicKagome<Lx, Ly>::Index& lattice,
    FourierTransformDirection direction) {
  Expression result;
  const float type_sign =
      (op.type() == Operator::Type::Annihilation) ? 1.0f : -1.0f;
  const float normalization = 1.0f / std::sqrt(static_cast<float>(Lx * Ly));
  constexpr float pi = std::numbers::pi_v<float>;
  const float factor_x = 2.0f * pi * type_sign / static_cast<float>(Lx);
  const float factor_y = 2.0f * pi * type_sign / static_cast<float>(Ly);

  auto [i, j, site] = lattice.from_orbital(op.orbital());

  //   Vec2 a1 = {1.5f + std::sqrt(3.0f), 1.0f + 0.5f * std::sqrt(3.0f)};
  //   Vec2 a2 = {-1.5f - std::sqrt(3.0f), 1.0f + std::sqrt(3.0f)};

  //   Vec2 d1 = {0.5f * std::sqrt(3.0f), 0.5f};
  //   Vec2 d2 = {0.5f * std::sqrt(3.0f), -0.5f};
  //   Vec2 d3 = {-1.0f, 0.0f};

  //   std::vector<Vec2> positions(6);
  //   positions[0] = Vec2{0.0f, 0.0f};
  //   positions[1] = positions[0] + (d2 - d1);
  //   positions[2] = positions[1] + d1;
  //   positions[3] = positions[2] - d3;
  //   positions[4] = positions[3] + d1;
  //   positions[5] = positions[4] + (d2 - d1);

  for (size_t kx = 0; kx < Lx; ++kx) {
    for (size_t ky = 0; ky < Ly; ++ky) {
      //   Vec2 R_ijs = {i * a1.x + j * a2.x + positions[site].x,
      //                 i * a1.y + j * a2.y + positions[site].y};
      std::complex<float> exponent(0.0f, factor_x * kx * i + factor_y * ky * j);
      std::complex<float> coefficient =
          std::exp(static_cast<float>(direction) * exponent) * normalization;

      Operator transformed_op(op.type(), op.spin(),
                              lattice.to_orbital(kx, ky, site),
                              op.statistics());
      result += Term(coefficient, {transformed_op});
    }
  }

  return result;
}

struct BandStructure {
  std::vector<std::array<float, 2>> k_points;
  std::vector<arma::fvec> eigenvalues;
  std::vector<arma::cx_fmat> eigenvectors;
};

template <size_t Lx, size_t Ly>
arma::cx_fmat construct_k_matrix(
    const Expression& H, const typename DiatomicKagome<Lx, Ly>::Index& indexer,
    size_t kx, size_t ky, size_t num_bands) {
  arma::cx_fmat matrix(num_bands, num_bands);
  // We don't have spin in this model
  Operator::Spin s = Operator::Spin::Up;
  for (size_t band1 = 0; band1 < num_bands; ++band1) {
    for (size_t band2 = 0; band2 < num_bands; ++band2) {
      const std::vector<Operator> needle{
          Operator::creation(s, indexer.to_orbital(kx, ky, band1)),
          Operator::annihilation(s, indexer.to_orbital(kx, ky, band2))};
      if (H.terms().contains(needle)) {
        matrix(band1, band2) = H.terms().at(needle);
      } else {
        matrix(band1, band2) = 0.0f;
      }
    }
  }
  return matrix;
}

template <size_t Lx, size_t Ly>
BandStructure diagonalize_multi_band_model(
    const Expression& H, const typename DiatomicKagome<Lx, Ly>::Index& indexer,
    size_t num_bands) {
  BandStructure result;

  for (size_t kx = 0; kx < Lx; ++kx) {
    for (size_t ky = 0; ky < Ly; ++ky) {
      arma::cx_fmat H_k_matrix =
          construct_k_matrix<Lx, Ly>(H, indexer, kx, ky, num_bands);
      arma::fvec eigenvalues;
      arma::cx_fmat eigenvectors;
      arma::eig_sym(eigenvalues, eigenvectors, H_k_matrix);

      float kx_value = 2.0f * std::numbers::pi_v<float> *
                       static_cast<float>(kx) / static_cast<float>(Lx);
      float ky_value = 2.0f * std::numbers::pi_v<float> *
                       static_cast<float>(ky) / static_cast<float>(Ly);
      result.k_points.push_back({kx_value, ky_value});
      result.eigenvalues.push_back(eigenvalues);
      result.eigenvectors.push_back(eigenvectors);
    }
  }
  return result;
}

template <size_t Lx, size_t Ly>
Expression transform_to_band_basis(
    const Operator& op, const BandStructure& band_structure,
    const typename DiatomicKagome<Lx, Ly>::Index& indexer, size_t num_bands) {
  Expression result;
  auto [kx, ky, alpha] = indexer.from_orbital(op.orbital());
  for (size_t n = 0; n < num_bands; ++n) {
    std::complex<float> coeff =
        band_structure.eigenvectors[kx * Ly + ky](alpha, n);
    Operator transformed_op(op.type(), op.spin(), indexer.to_orbital(kx, ky, n),
                            op.statistics());
    if (op.type() == Operator::Type::Annihilation) {
      result += Term(coeff, {transformed_op});
    } else {
      result += Term(std::conj(coeff), {transformed_op});
    }
  }
  return result;
}

template <size_t Lx, size_t Ly>
Expression inverse_transform_to_band_basis(
    const Operator& op, const BandStructure& band_structure,
    const typename DiatomicKagome<Lx, Ly>::Index& indexer, size_t num_bands) {
  Expression result;
  auto [kx, ky, n] = indexer.from_orbital(op.orbital());
  for (size_t alpha = 0; alpha < num_bands; ++alpha) {
    std::complex<float> coeff =
        band_structure.eigenvectors[kx * Ly + ky](alpha, n);
    Operator transformed_op(op.type(), op.spin(),
                            indexer.to_orbital(kx, ky, alpha), op.statistics());
    if (op.type() == Operator::Type::Annihilation) {
      result += Term(coeff, {transformed_op});
    } else {
      result += Term(std::conj(coeff), {transformed_op});
    }
  }
  return result;
}

int main() {
  const size_t Lx = 2;
  const size_t Ly = 2;
  const float t1 = 0.535f;
  const float t2 = 0.0258f;
  const float t3 = 0.0261f;

  DiatomicKagome<Lx, Ly> model(t1, t2, t3);

  std::cout << "Diatomic Kagome Hamiltonian in real space:" << "\n";

  for (const auto& [ops, coeff] : model.hamiltonian().terms()) {
    std::cout << "  " << term_to_string(Term(coeff, ops), model.index())
              << "\n";
  }
  std::cout << "\n";

  std::cout << "Diatomic Kagome Hamiltonian in k-space:" << "\n";

  Expression momentum_hamiltonian = transform_expression(
      fourier_transform_operator_2d<Lx, Ly>, model.hamiltonian(), model.index(),
      FourierTransformDirection::Forward);

  for (const auto& [ops, coeff] : momentum_hamiltonian.terms()) {
    std::cout << "  " << term_to_string(Term(coeff, ops), model.index())
              << "\n";
  }
  std::cout << "\n";

  std::cout << "Diatomic Kagome Hamiltonian in diagonal form:" << "\n";

  const size_t num_bands = 6;
  BandStructure band_structure = diagonalize_multi_band_model<Lx, Ly>(
      momentum_hamiltonian, model.index(), num_bands);

  std::ofstream out_file("band_structure.dat");
  for (size_t i = 0; i < band_structure.k_points.size(); ++i) {
    out_file << band_structure.k_points[i][0] << " "
             << band_structure.k_points[i][1];
    for (size_t j = 0; j < num_bands; ++j) {
      out_file << " " << band_structure.eigenvalues[i][j];
    }
    out_file << "\n";
  }

  Expression transformed_op = transform_expression(
      transform_to_band_basis<Lx, Ly>, momentum_hamiltonian, band_structure,
      model.index(), num_bands);

  for (const auto& [ops, coeff] : transformed_op.terms()) {
    std::cout << "  " << term_to_string(Term(coeff, ops), model.index())
              << "\n";
  }
  std::cout << "\n";

  std::cout << "Inverse transformation:" << "\n";

  Expression inverse_transformed_op = transform_expression(
      inverse_transform_to_band_basis<Lx, Ly>, transformed_op, band_structure,
      model.index(), num_bands);

  for (const auto& [ops, coeff] : inverse_transformed_op.terms()) {
    std::cout << "  " << term_to_string(Term(coeff, ops), model.index())
              << "\n";
  }
  std::cout << "\n";

  return 0;
}