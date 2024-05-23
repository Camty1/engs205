#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct Slice {
  size_t start, stop, step;
  Slice(size_t start, size_t stop, size_t step = 1) : start(start), stop(stop), step(step) {};
};

template <typename T> class Array {
private:
  std::vector<T> data;
  size_t rows, cols;

public:
  Array(size_t rows, size_t cols = 1);

  T &operator()(size_t row, size_t col = 0);
  const T &operator()(size_t row, size_t col = 0) const;
  T *dataPtr();

  const size_t get_rows() const;
  const size_t get_cols() const;

  const Array<T> get_slice(Slice slice);
  const Array<T> get_slice(Slice row_slice, Slice col_slice);
  const Array<T> get_slice(Slice row_slice, size_t col);
  const Array<T> get_slice(size_t row, Slice col_slice);

  void set_slice(const Array<T> &A, Slice slice);
  void set_slice(const Array<T> &A, Slice row_slice, Slice col_slice);
  void set_slice(const Array<T> &A, size_t row, Slice slice);
  void set_slice(const Array<T> &A, Slice slice, size_t col);

  // Unary operators
  void operator+=(const Array<T> &B);
  void operator-=(const Array<T> &B);
  void operator*=(const Array<T> &B);
  void operator/=(const Array<T> &B);

  const Array<T> operator+(const Array<T> &B);
  const Array<T> operator-(const Array<T> &B);
  const Array<T> operator*(const Array<T> &B);
  const Array<T> operator/(const Array<T> &B);

  void fill(const T &value);
  void print() const;
  void write(const std::string filename) const;

private:
  // Slicing stuff
  
  void check_slice(Slice slice);

};

#include "Array.tpp"

#endif
