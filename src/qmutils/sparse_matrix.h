#pragma once

#include <stdexcept>
#include <unordered_map>

namespace qmutils {

template <typename T>
class SparseMatrix {
 public:
  using value_type = T;

  SparseMatrix() : m_rows(0), m_cols(0) {}

  SparseMatrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {}

  SparseMatrix(const SparseMatrix& other)
      : m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data) {}

  T& operator()(size_t row, size_t col) {
    check_bounds(row, col);
    return m_data[{row, col}];
  }

  const T& operator()(size_t row, size_t col) const {
    check_bounds(row, col);
    auto it = m_data.find({row, col});
    if (it == m_data.end()) {
      return T{};
    }
    return it->second;
  }

  bool contains(size_t row, size_t col) const {
    check_bounds(row, col);
    return m_data.find({row, col}) != m_data.end();
  }

  size_t rows() const { return m_rows; }
  size_t cols() const { return m_cols; }
  size_t count_non_zero() const { return m_data.size(); }
  void clear() { m_data.clear(); }

  auto begin() { return m_data.begin(); }
  auto end() { return m_data.end(); }
  auto begin() const { return m_data.begin(); }
  auto end() const { return m_data.end(); }

 private:
  void check_bounds(size_t row, size_t col) const {
    if (row >= m_rows || col >= m_cols) {
      throw std::out_of_range("Matrix indices out of range");
    }
  }

  struct Index {
    size_t row;
    size_t col;

    bool operator==(const Index& other) const {
      return row == other.row && col == other.col;
    }
  };

  struct IndexHash {
    size_t operator()(const Index& index) const {
      return std::hash<size_t>()(index.row) ^
             (std::hash<size_t>()(index.col) << 1);
    }
  };

  size_t m_rows;
  size_t m_cols;
  std::unordered_map<Index, T, IndexHash> m_data;
};

}  // namespace qmutils
