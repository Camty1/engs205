#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

template <typename T> class Array {
private:
  std::vector<T> data;
  size_t rows, cols;

public:
  // Modified constructor to allow for 1D or  array creation
  Array(size_t rows, size_t cols = 1);

  // Accessor methods modified for 1D or  access
  T &operator()(size_t row, size_t col = 0);
  const T &operator()(size_t row, size_t col = 0) const;

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
};

#endif
