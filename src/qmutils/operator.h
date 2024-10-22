#pragma once

#include <bit>
#include <cstdint>
#include <string>
#include <type_traits>

namespace qmutils {

class Operator {
 public:
  enum class Type : uint8_t { Creation = 0, Annihilation = 1 };
  enum class Spin : uint8_t { Up = 0, Down = 1 };

  static constexpr uint8_t ORBITAL_BITFIELD_WIDTH = 6;

  static constexpr size_t max_orbital_size() {
    return 1 << ORBITAL_BITFIELD_WIDTH;
  }

  Operator() = default;

  constexpr Operator(Type type, Spin spin, uint8_t orbital) noexcept
      : data_{orbital, static_cast<uint8_t>(spin), static_cast<uint8_t>(type)} {
  }

  [[nodiscard]] constexpr Type type() const noexcept {
    return static_cast<Type>(data_.type);
  }

  [[nodiscard]] constexpr uint8_t orbital() const noexcept {
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
                                                   uint8_t orbital) noexcept {
    return Operator(Type::Creation, spin, orbital);
  }

  [[nodiscard]] static constexpr Operator annihilation(
      Spin spin, uint8_t orbital) noexcept {
    return Operator(Type::Annihilation, spin, orbital);
  }

  [[nodiscard]] constexpr uint8_t data() const noexcept {
    return std::bit_cast<uint8_t>(data_);
  }

 private:
  struct Data {
    uint8_t orbital : ORBITAL_BITFIELD_WIDTH;
    uint8_t spin : 1;
    uint8_t type : 1;
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
