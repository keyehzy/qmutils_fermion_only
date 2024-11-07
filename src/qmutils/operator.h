#pragma once

#include <bit>
#include <cstdint>
#include <string>
#include <type_traits>

namespace qmutils {

class Operator {
 public:
  using int_type = uint16_t;

  enum class Type : int_type { Creation = 0, Annihilation = 1 };
  enum class Spin : int_type { Up = 0, Down = 1 };
  enum class Statistics : int_type { Fermionic = 0, Bosonic = 1 };

  static constexpr size_t ORBITAL_BITFIELD_WIDTH = sizeof(int_type) * 8 - 3;

  static constexpr size_t max_orbital_size() {
    return 1 << ORBITAL_BITFIELD_WIDTH;
  }

  Operator() = default;

  constexpr Operator(Type type, Spin spin, size_t orbital,
                     Statistics stat = Statistics::Fermionic) noexcept
      : data_{.orbital = static_cast<int_type>(orbital),
              .spin = static_cast<int_type>(spin),
              .type = static_cast<int_type>(type),
              .statistics = static_cast<int_type>(stat)} {}

  [[nodiscard]] constexpr Type type() const noexcept {
    return static_cast<Type>(data_.type);
  }

  [[nodiscard]] constexpr size_t orbital() const noexcept {
    return data_.orbital;
  }

  [[nodiscard]] constexpr Spin spin() const noexcept {
    return static_cast<Spin>(data_.spin);
  }

  [[nodiscard]] constexpr Statistics statistics() const noexcept {
    return static_cast<Statistics>(data_.statistics);
  }

  constexpr bool operator<(const Operator &other) const noexcept {
    return data() < other.data();
  }

  constexpr bool operator>(const Operator &other) const noexcept {
    return data() > other.data();
  }

  constexpr bool operator==(const Operator &other) const noexcept {
    return data() == other.data();
  }

  constexpr bool operator>=(const Operator &other) const noexcept {
    return data() >= other.data();
  }

  constexpr bool operator<=(const Operator &other) const noexcept {
    return data() <= other.data();
  }

  [[nodiscard]] constexpr bool commutes_with(
      const Operator &other) const noexcept {
    return (data_.statistics != other.data_.statistics) ||
           ((data_.type == other.data_.type) ||
            ((data_.spin != other.data_.spin) ||
             (data_.orbital != other.data_.orbital)));
  }

  [[nodiscard]] constexpr Operator adjoint() const noexcept {
    return Operator(static_cast<Type>(!data_.type),
                    static_cast<Spin>(data_.spin), data_.orbital,
                    static_cast<Statistics>(data_.statistics));
  }

  [[nodiscard]] constexpr Operator flip_spin() const noexcept {
    return Operator(static_cast<Type>(data_.type),
                    static_cast<Spin>(!data_.spin), data_.orbital,
                    static_cast<Statistics>(data_.statistics));
  }

  [[nodiscard]] std::string to_string() const;

  [[nodiscard]] static constexpr Operator creation(Spin spin,
                                                   int_type orbital) noexcept {
    return Operator(Type::Creation, spin, orbital);
  }

  [[nodiscard]] static constexpr Operator annihilation(
      Spin spin, int_type orbital) noexcept {
    return Operator(Type::Annihilation, spin, orbital);
  }

  [[nodiscard]] constexpr int_type data() const noexcept {
    return std::bit_cast<int_type>(data_);
  }

  static constexpr bool is_fermion(Operator op) noexcept {
    return op.statistics() == Statistics::Fermionic;
  }

  static constexpr bool is_boson(Operator op) noexcept {
    return op.statistics() == Statistics::Bosonic;
  }

  struct Fermion;
  struct Boson;

 private:
  struct Data {
    int_type orbital : ORBITAL_BITFIELD_WIDTH;
    int_type spin : 1;
    int_type type : 1;
    int_type statistics : 1;
  };
  static_assert(std::has_unique_object_representations_v<Data>);

  Data data_;
};
static_assert(std::is_trivially_copyable_v<Operator>,
              "Operator must be trivially copyable");
static_assert(sizeof(Operator) == 2, "Operator must be 2 byte in size");

struct Operator::Fermion {
  static constexpr Operator creation(Spin spin, int_type orbital) noexcept {
    return Operator(Type::Creation, spin, orbital, Statistics::Fermionic);
  }

  static constexpr Operator annihilation(Spin spin, int_type orbital) noexcept {
    return Operator(Type::Annihilation, spin, orbital, Statistics::Fermionic);
  }
};

struct Operator::Boson {
  static constexpr Operator creation(Spin spin, int_type orbital) noexcept {
    return Operator(Type::Creation, spin, orbital, Statistics::Bosonic);
  }

  static constexpr Operator annihilation(Spin spin, int_type orbital) noexcept {
    return Operator(Type::Annihilation, spin, orbital, Statistics::Bosonic);
  }
};
}  // namespace qmutils

template <>
struct std::hash<qmutils::Operator> {
  size_t operator()(const qmutils::Operator &op) const noexcept {
    return static_cast<size_t>(op.data());
  }
};
