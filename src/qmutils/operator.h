#pragma once

#include <bit>
#include <cstdint>
#include <string>
#include <type_traits>

namespace qmutils {

class Operator {
 public:
  using int_type = uint8_t;

  enum class Type : int_type { Creation = 0, Annihilation = 1 };
  enum class Spin : int_type { Up = 0, Down = 1 };

  static constexpr size_t ORBITAL_BITFIELD_WIDTH = sizeof(int_type) * 8 - 2;

  static constexpr size_t max_orbital_size() {
    return 1 << ORBITAL_BITFIELD_WIDTH;
  }

  Operator() = default;

  constexpr Operator(Type type, Spin spin, size_t orbital) noexcept
      : data_{static_cast<int_type>(orbital), static_cast<int_type>(spin),
              static_cast<int_type>(type)} {}

  [[nodiscard]] constexpr Type type() const noexcept {
    return static_cast<Type>(data_.type);
  }

  [[nodiscard]] constexpr size_t orbital() const noexcept {
    return data_.orbital;
  }

  [[nodiscard]] constexpr Spin spin() const noexcept {
    return static_cast<Spin>(data_.spin);
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
    return (data_.type == other.data_.type) ||
           ((data_.spin != other.data_.spin) ||
            (data_.orbital != other.data_.orbital));
  }

  [[nodiscard]] constexpr Operator adjoint() const noexcept {
    return Operator(static_cast<Type>(!data_.type),
                    static_cast<Spin>(data_.spin), data_.orbital);
  }

  [[nodiscard]] constexpr Operator flip_spin() const noexcept {
    return Operator(type(), spin() == Spin::Up ? Spin::Down : Spin::Up,
                    orbital());
  }

  std::string to_string() const;

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

 private:
  struct Data {
    int_type orbital : ORBITAL_BITFIELD_WIDTH;
    int_type spin : 1;
    int_type type : 1;
  };
  static_assert(std::has_unique_object_representations_v<Data>);

  Data data_;
};
static_assert(std::is_trivially_copyable_v<Operator>,
              "Operator must be trivially copyable");
static_assert(sizeof(Operator) == 1, "Operator must be 1 byte in size");

}  // namespace qmutils

template <>
struct std::hash<qmutils::Operator> {
  size_t operator()(const qmutils::Operator &op) const noexcept {
    return static_cast<size_t>(op.data());
  }
};
