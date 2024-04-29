// Array2D.hpp
#ifndef ARRAY2D_HPP
#define ARRAY2D_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

template <typename T> class Array2D {
private:
  std::vector<T> data;
  size_t rows, cols;

public:
  Array2D(size_t rows, size_t cols);
  T &operator()(size_t row, size_t col);
  const T &operator()(size_t row, size_t col) const;
  Array2D<T> operator-() const;
  T *dataPtr();
  size_t getRows() const;
  size_t getCols() const;
  void fill(const T &value);
  void fillRange(size_t startRow, size_t startCol, size_t endRow, size_t endCol,
                 const T &value);
  void print() const;
  void write(std::string filename) const;
  void insert(const Array2D<T> &source, size_t sourceStartRow,
              size_t sourceStartCol, size_t subsetRows, size_t subsetCols,
              size_t targetStartRow, size_t targetStartCol);
  std::vector<T> getRowVec(size_t row) const;
  std::vector<T> getColVec(size_t col) const;
};

#include "Array2D.tpp" // Include template implementation

#endif // ARRAY2D_HPP
