#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace qmutils {

template <typename T>
class SparseRow {
 private:
  std::vector<uint8_t> mask_;
  std::vector<T> values_;

 public:
  explicit SparseRow(size_t size = 0) { resize(size); }

  void resize(size_t size) { mask_.resize((size + 7) / 8, 0); }

  void set(size_t pos, const T& value) {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;

    if (byte_pos >= mask_.size())
      throw std::out_of_range("Position exceeds row size");
    bool has_value = mask_[byte_pos] & (static_cast<uint8_t>(1U << bit_pos));

    if (value != T{}) {
      if (!has_value) {
        size_t insert_pos = 0;
        for (size_t i = 0; i < byte_pos; ++i) {
          insert_pos += static_cast<size_t>(
              __builtin_popcount(static_cast<unsigned int>(mask_[i])));
        }
        insert_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[byte_pos]) &
                               ((1U << bit_pos) - 1U)));

        values_.insert(values_.begin() + insert_pos, value);
        mask_[byte_pos] |= static_cast<uint8_t>(1U << bit_pos);
      } else {
        size_t value_pos = 0;
        for (size_t i = 0; i < byte_pos; ++i) {
          value_pos += static_cast<size_t>(
              __builtin_popcount(static_cast<unsigned int>(mask_[i])));
        }
        value_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[byte_pos]) &
                               ((1U << bit_pos) - 1U)));
        values_[value_pos] = value;
      }
    } else if (has_value) {
      size_t value_pos = 0;
      for (size_t i = 0; i < byte_pos; ++i) {
        value_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[i])));
      }
      value_pos += static_cast<size_t>(__builtin_popcount(
          static_cast<unsigned int>(mask_[byte_pos]) & ((1U << bit_pos) - 1U)));

      values_.erase(values_.begin() + value_pos);
      mask_[byte_pos] &= static_cast<uint8_t>(~(1U << bit_pos));
    }
  }

  T get(size_t pos) const {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;

    if (byte_pos >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    if (mask_[byte_pos] & (static_cast<uint8_t>(1U << bit_pos))) {
      size_t value_pos = 0;
      for (size_t i = 0; i < byte_pos; ++i) {
        value_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[i])));
      }
      value_pos += static_cast<size_t>(__builtin_popcount(
          static_cast<unsigned int>(mask_[byte_pos]) & ((1U << bit_pos) - 1U)));
      return values_[value_pos];
    }
    return T{};
  }

  T& get_or_create_ref(size_t pos) {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;

    if (byte_pos >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    bool has_value = mask_[byte_pos] & (static_cast<uint8_t>(1U << bit_pos));
    if (!has_value) {
      // Create new value position
      size_t insert_pos = 0;
      for (size_t i = 0; i < byte_pos; ++i) {
        insert_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[i])));
      }
      insert_pos += static_cast<size_t>(__builtin_popcount(
          static_cast<unsigned int>(mask_[byte_pos]) & ((1U << bit_pos) - 1U)));

      values_.insert(values_.begin() + static_cast<ptrdiff_t>(insert_pos), T{});
      mask_[byte_pos] |= static_cast<uint8_t>(1U << bit_pos);
      return values_[insert_pos];
    } else {
      // Return reference to existing value
      size_t value_pos = 0;
      for (size_t i = 0; i < byte_pos; ++i) {
        value_pos += static_cast<size_t>(
            __builtin_popcount(static_cast<unsigned int>(mask_[i])));
      }
      value_pos += static_cast<size_t>(__builtin_popcount(
          static_cast<unsigned int>(mask_[byte_pos]) & ((1U << bit_pos) - 1U)));
      return values_[value_pos];
    }
  }

  size_t size() const { return mask_.size() * 8; }
  size_t value_count() const { return values_.size(); }
};

template <typename T>
class SparseMatrix {
 private:
  std::vector<SparseRow<T>> rows_;
  size_t num_cols_;

 public:
  SparseMatrix(size_t rows = 0, size_t cols = 0) : num_cols_(cols) {
    resize(rows, cols);
  }

  void resize(size_t rows, size_t cols) {
    rows_.resize(rows);
    num_cols_ = cols;
    for (auto& row : rows_) {
      row.resize(cols);
    }
  }

  T operator()(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      throw std::out_of_range("Matrix index out of bounds");
    }
    return rows_[row].get(col);
  }

  T& operator()(size_t row, size_t col) {
    if (row >= rows_.size() || col >= num_cols_) {
      throw std::out_of_range("Matrix index out of bounds");
    }
    return rows_[row].get_or_create_ref(col);
  }

  size_t rows() const { return rows_.size(); }
  size_t cols() const { return num_cols_; }

  size_t non_zero_count() const {
    size_t count = 0;
    for (const auto& row : rows_) {
      count += row.value_count();
    }
    return count;
  }

  const SparseRow<T>& row(size_t index) const {
    if (index >= rows_.size()) {
      throw std::out_of_range("Row index out of bounds");
    }
    return rows_[index];
  }
};

}  // namespace qmutils
