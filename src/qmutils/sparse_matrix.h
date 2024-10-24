#pragma once

#include <bitset>
#include <complex>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace qmutils {

template <typename T>
class SparseMatrix;

// Basic numeric sparse matrices
using SpMat_i = SparseMatrix<int>;
using SpMat_f = SparseMatrix<float>;
using SpMat_d = SparseMatrix<double>;

// Complex number sparse matrices
using SpMat_cf = SparseMatrix<std::complex<float>>;
using SpMat_cd = SparseMatrix<std::complex<double>>;

template <typename T, bool IsConst>
class SparseMatrixIteratorImpl;

template <typename T>
using SparseMatrixIterator = SparseMatrixIteratorImpl<T, false>;
template <typename T>
using SparseMatrixConstIterator = SparseMatrixIteratorImpl<T, true>;

template <typename T, bool IsConst>
struct SparseElementImpl {
  size_t row;
  size_t col;
  typename std::conditional_t<IsConst, const T&, T&> value;

  SparseElementImpl(size_t r, size_t c,
                    typename std::conditional_t<IsConst, const T&, T&> v)
      : row(r), col(c), value(v) {}
};

template <typename T>
class SparseRow {
 private:
  static constexpr size_t BITS_PER_BYTE = 8;
  std::vector<uint8_t> mask_;
  std::vector<T> values_;

  size_t byte_index(size_t pos) const { return pos / BITS_PER_BYTE; }
  size_t bit_index(size_t pos) const { return pos % BITS_PER_BYTE; }

  uint8_t bit_mask(size_t bit_pos) const {
    return static_cast<uint8_t>(1U << bit_pos);
  }

 public:
  explicit SparseRow(size_t size = 0) { resize(size); }

  void resize(size_t size) {
    mask_.resize((size + BITS_PER_BYTE - 1) / BITS_PER_BYTE, 0);
  }

  size_t calculate_value_index(size_t pos) const {
    size_t byte_pos = byte_index(pos);
    size_t bit_pos = bit_index(pos);
    size_t count = 0;

    // Count set bits in previous bytes
    for (size_t i = 0; i < byte_pos; ++i) {
      count += std::bitset<BITS_PER_BYTE>(mask_[i]).count();
    }

    // Count set bits in current byte up to target position
    uint8_t partial_mask = mask_[byte_pos] & ((1U << bit_pos) - 1U);
    count += std::bitset<BITS_PER_BYTE>(partial_mask).count();

    return count;
  }

  void set(size_t pos, const T& value) {
    if (byte_index(pos) >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    bool has_value = test(pos);
    if (value != T{}) {
      if (!has_value) {
        insert_value(pos, value);
      } else {
        update_value(pos, value);
      }
    } else if (has_value) {
      remove_value(pos);
    }
  }

  T get(size_t pos) const {
    if (byte_index(pos) >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    return test(pos) ? values_[calculate_value_index(pos)] : T{};
  }

  T& get_or_create_ref(size_t pos) {
    if (byte_index(pos) >= mask_.size()) {
      throw std::out_of_range("Position exceeds row size");
    }

    if (!test(pos)) {
      insert_value(pos, T{});
    }
    return values_[calculate_value_index(pos)];
  }

  bool test(size_t pos) const {
    size_t byte_pos = byte_index(pos);
    return (mask_[byte_pos] & bit_mask(bit_index(pos))) != 0;
  }

  size_t size() const { return mask_.size() * BITS_PER_BYTE; }
  size_t value_count() const { return values_.size(); }
  uint8_t get_mask(size_t byte_pos) const {
    return byte_pos < mask_.size() ? mask_[byte_pos] : 0;
  }
  const std::vector<T>& values() const { return values_; }

 private:
  void insert_value(size_t pos, const T& value) {
    size_t insert_pos = calculate_value_index(pos);
    values_.insert(values_.begin() + static_cast<ptrdiff_t>(insert_pos), value);
    mask_[byte_index(pos)] |= bit_mask(bit_index(pos));
  }

  void update_value(size_t pos, const T& value) {
    values_[calculate_value_index(pos)] = value;
  }

  void remove_value(size_t pos) {
    size_t value_pos = calculate_value_index(pos);
    values_.erase(values_.begin() + static_cast<ptrdiff_t>(value_pos));
    mask_[byte_index(pos)] &= ~bit_mask(bit_index(pos));
  }
};

template <typename T>
class SparseMatrix {
 private:
  std::vector<SparseRow<T>> rows_;
  size_t num_cols_;

 public:
  using iterator = SparseMatrixIterator<T>;
  using const_iterator = SparseMatrixConstIterator<T>;

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

  void validate_indices(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      throw std::out_of_range("Matrix index out of bounds");
    }
  }

  T get(size_t row, size_t col) const {
    validate_indices(row, col);
    return rows_[row].get(col);
  }

  void set(size_t row, size_t col, const T& value) {
    validate_indices(row, col);
    rows_[row].set(col, value);
  }

  T operator()(size_t row, size_t col) const { return get(row, col); }
  T& operator()(size_t row, size_t col) {
    validate_indices(row, col);
    return rows_[row].get_or_create_ref(col);
  }

  size_t rows() const { return rows_.size(); }
  size_t cols() const { return num_cols_; }

  size_t non_zero_count() const {
    return std::accumulate(
        rows_.begin(), rows_.end(), size_t{0},
        [](size_t sum, const auto& row) { return sum + row.value_count(); });
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

  iterator begin() { return iterator(*this, 0, 0); }
  iterator end() { return iterator(*this, rows(), 0); }
  const_iterator begin() const { return const_iterator(*this, 0, 0); }
  const_iterator end() const { return const_iterator(*this, rows(), 0); }
  const_iterator cbegin() const { return const_iterator(*this, 0, 0); }
  const_iterator cend() const { return const_iterator(*this, rows(), 0); }

  iterator find(size_t row, size_t col) {
    return find_impl<iterator>(row, col);
  }

  const_iterator find(size_t row, size_t col) const {
    return find_impl<const_iterator>(row, col);
  }

  bool contains(size_t row, size_t col) const {
    if (row >= rows_.size() || col >= num_cols_) {
      return false;
    }
    return rows_[row].test(col);
  }

 private:
  template <typename IteratorType>
  IteratorType find_impl(size_t row, size_t col) {
    if (!contains(row, col)) {
      return IteratorType(*this, rows(), 0);
    }
    return IteratorType(*this, row, rows_[row].calculate_value_index(col));
  }
};

template <typename T, bool IsConst>
class SparseMatrixIteratorImpl {
  using MatrixType =
      std::conditional_t<IsConst, const SparseMatrix<T>, SparseMatrix<T>>;
  using ElementType = SparseElementImpl<T, IsConst>;

 private:
  MatrixType& matrix_;
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
  using value_type = ElementType;
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;

  SparseMatrixIteratorImpl(MatrixType& matrix, size_t row, size_t value_index)
      : matrix_(matrix), current_row_(row), current_value_index_(value_index) {
    find_next_value();
  }

  ElementType operator*() const {
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
            return ElementType(
                current_row_, col,
                const_cast<MatrixType&>(matrix_)(current_row_, col));
          }
          ++count;
        }
      }
    }
    throw std::runtime_error("Invalid iterator state");
  }

  SparseMatrixIteratorImpl& operator++() {
    ++current_value_index_;
    find_next_value();
    return *this;
  }

  SparseMatrixIteratorImpl operator++(int) {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }

  bool operator==(const SparseMatrixIteratorImpl& other) const {
    return &matrix_ == &other.matrix_ && current_row_ == other.current_row_ &&
           current_value_index_ == other.current_value_index_;
  }

  bool operator!=(const SparseMatrixIteratorImpl& other) const {
    return !(*this == other);
  }
};

}  // namespace qmutils
