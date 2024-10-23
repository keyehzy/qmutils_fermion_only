#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace qmutils {

template <typename T>
class SparseMatrix;

template <typename T>
class SparseMatrixIterator;

template <typename T>
class SparseMatrixConstIterator;

template <typename T>
struct SparseElement {
  size_t row;
  size_t col;
  T& value;

  SparseElement(size_t r, size_t c, T& v) : row(r), col(c), value(v) {}
};

template <typename T>
struct SparseConstElement {
  size_t row;
  size_t col;
  const T& value;

  SparseConstElement(size_t r, size_t c, const T& v)
      : row(r), col(c), value(v) {}
};

template <typename T>
class SparseRow {
 private:
  std::vector<uint8_t> mask_;
  std::vector<T> values_;

 public:
  explicit SparseRow(size_t size = 0) { resize(size); }

  void resize(size_t size) { mask_.resize((size + 7) / 8, 0); }

  size_t calculate_value_pos(size_t pos) const {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;
    int value_pos = 0;

    for (size_t i = 0; i < byte_pos; ++i) {
      value_pos += __builtin_popcount(static_cast<unsigned int>(mask_[i]));
    }
    value_pos += __builtin_popcount(static_cast<unsigned int>(mask_[byte_pos]) &
                                    ((1U << bit_pos) - 1U));

    return static_cast<size_t>(value_pos);
  }

  void set(size_t pos, const T& value) {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;

    if (byte_pos >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    bool has_value = mask_[byte_pos] & (1U << bit_pos);

    if (value != T{}) {
      if (!has_value) {
        size_t insert_pos = calculate_value_pos(pos);
        values_.insert(values_.begin() + static_cast<ptrdiff_t>(insert_pos),
                       value);
        mask_[byte_pos] |= (1U << bit_pos);
      } else {
        size_t value_pos = calculate_value_pos(pos);
        values_[value_pos] = value;
      }
    } else if (has_value) {
      size_t value_pos = calculate_value_pos(pos);
      values_.erase(values_.begin() + static_cast<ptrdiff_t>(value_pos));
      mask_[byte_pos] &= ~(1U << bit_pos);
    }
  }

  T get(size_t pos) const {
    size_t byte_pos = pos / 8;
    size_t bit_pos = pos % 8;

    if (byte_pos >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    if (mask_[byte_pos] & (1U << bit_pos)) {
      size_t value_pos = calculate_value_pos(pos);
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

    bool has_value = mask_[byte_pos] & (1U << bit_pos);
    if (!has_value) {
      size_t insert_pos = calculate_value_pos(pos);
      values_.insert(values_.begin() + static_cast<ptrdiff_t>(insert_pos), T{});
      mask_[byte_pos] |= (1U << bit_pos);
      return values_[insert_pos];
    } else {
      size_t value_pos = calculate_value_pos(pos);
      return values_[value_pos];
    }
  }

  size_t size() const { return mask_.size() * 8; }
  size_t value_count() const { return values_.size(); }

  inline uint8_t get_mask(size_t byte_pos) const {
    return byte_pos < mask_.size() ? mask_[byte_pos] : 0;
  }

  const std::vector<T>& values() const { return values_; }
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

  T get(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      throw std::out_of_range("Matrix index out of bounds");
    }
    return rows_[row].get(col);
  }

  void set(size_t row, size_t col, const T& value) {
    if (row >= rows_.size() || col >= num_cols_) {
      throw std::out_of_range("Matrix index out of bounds");
    }
    return rows_[row].set(col, value);
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

  void clear() {
    size_t num_rows = rows_.size();
    rows_.clear();
    rows_.resize(num_rows);
    for (auto& row : rows_) {
      row.resize(num_cols_);
    }
  }

  using iterator = SparseMatrixIterator<T>;
  using const_iterator = SparseMatrixConstIterator<T>;

  iterator begin() { return iterator(*this, 0, 0); }

  iterator end() { return iterator(*this, rows(), 0); }

  const_iterator begin() const { return const_iterator(*this, 0, 0); }

  const_iterator end() const { return const_iterator(*this, rows(), 0); }

  const_iterator cbegin() const { return const_iterator(*this, 0, 0); }

  const_iterator cend() const { return const_iterator(*this, rows(), 0); }

  iterator find(size_t row, size_t col) {
    if (row >= rows_.size() || col >= num_cols_) {
      return end();
    }

    const auto& sparse_row = rows_[row];
    size_t byte_pos = col / 8;
    size_t bit_pos = col % 8;

    if (!(sparse_row.get_mask(byte_pos) & (1U << bit_pos))) {
      return end();
    }

    size_t value_index = sparse_row.calculate_value_pos(col);
    return iterator(*this, row, value_index);
  }

  const_iterator find(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      return cend();
    }

    const auto& sparse_row = rows_[row];
    size_t byte_pos = col / 8;
    size_t bit_pos = col % 8;

    if (!(sparse_row.get_mask(byte_pos) & (1U << bit_pos))) {
      return cend();
    }

    size_t value_index = sparse_row.calculate_value_pos(col);
    return const_iterator(*this, row, value_index);
  }

  bool contains(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      return false;
    }

    const auto& sparse_row = rows_[row];
    size_t byte_pos = col / 8;
    size_t bit = col % 8;

    return (sparse_row.get_mask(byte_pos) & (1U << bit)) != 0;
  }
};

template <typename T>
class SparseMatrixIterator {
 private:
  SparseMatrix<T>& matrix_;
  size_t current_row_;
  size_t current_value_index_;

  void find_next_value() {
    while (current_row_ < matrix_.rows()) {
      const auto& row = matrix_.row(current_row_);
      if (current_value_index_ < row.value_count()) {
        return;
      }
      ++current_row_;
      current_value_index_ = 0;
    }
  }

 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = SparseElement<T>;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  SparseMatrixIterator(SparseMatrix<T>& matrix, size_t row, size_t value_index)
      : matrix_(matrix), current_row_(row), current_value_index_(value_index) {
    find_next_value();
  }

  value_type operator*() {
    const auto& row = matrix_.row(current_row_);
    size_t col = 0;
    size_t count = 0;

    for (size_t byte_pos = 0; byte_pos < (matrix_.cols() + 7) / 8; ++byte_pos) {
      uint8_t mask = row.get_mask(byte_pos);
      for (size_t bit = 0; bit < 8 && (byte_pos * 8 + bit) < matrix_.cols();
           ++bit) {
        if (mask & (1U << bit)) {
          if (count == current_value_index_) {
            col = byte_pos * 8 + bit;
            return value_type(current_row_, col, matrix_(current_row_, col));
          }
          ++count;
        }
      }
    }
    throw std::runtime_error("Invalid iterator state");
  }

  SparseMatrixIterator& operator++() {
    ++current_value_index_;
    find_next_value();
    return *this;
  }

  SparseMatrixIterator operator++(int) {
    SparseMatrixIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const SparseMatrixIterator& other) const {
    return &matrix_ == &other.matrix_ && current_row_ == other.current_row_ &&
           current_value_index_ == other.current_value_index_;
  }

  bool operator!=(const SparseMatrixIterator& other) const {
    return !(*this == other);
  }
};

template <typename T>
class SparseMatrixConstIterator {
 private:
  const SparseMatrix<T>& matrix_;
  size_t current_row_;
  size_t current_value_index_;

  void find_next_value() {
    while (current_row_ < matrix_.rows()) {
      const auto& row = matrix_.row(current_row_);
      if (current_value_index_ < row.value_count()) {
        return;
      }
      ++current_row_;
      current_value_index_ = 0;
    }
  }

 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = SparseConstElement<T>;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  SparseMatrixConstIterator(const SparseMatrix<T>& matrix, size_t row,
                            size_t value_index)
      : matrix_(matrix), current_row_(row), current_value_index_(value_index) {
    find_next_value();
  }

  value_type operator*() const {
    const auto& row = matrix_.row(current_row_);
    size_t col = 0;
    size_t count = 0;

    for (size_t byte_pos = 0; byte_pos < (matrix_.cols() + 7) / 8; ++byte_pos) {
      uint8_t mask = row.get_mask(byte_pos);
      for (size_t bit = 0; bit < 8 && (byte_pos * 8 + bit) < matrix_.cols();
           ++bit) {
        if (mask & (1U << bit)) {
          if (count == current_value_index_) {
            col = byte_pos * 8 + bit;
            return value_type(current_row_, col, matrix_(current_row_, col));
          }
          ++count;
        }
      }
    }
    throw std::runtime_error("Invalid iterator state");
  }

  SparseMatrixConstIterator& operator++() {
    ++current_value_index_;
    find_next_value();
    return *this;
  }

  SparseMatrixConstIterator operator++(int) {
    SparseMatrixConstIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const SparseMatrixConstIterator& other) const {
    return &matrix_ == &other.matrix_ && current_row_ == other.current_row_ &&
           current_value_index_ == other.current_value_index_;
  }

  bool operator!=(const SparseMatrixConstIterator& other) const {
    return !(*this == other);
  }
};

}  // namespace qmutils
