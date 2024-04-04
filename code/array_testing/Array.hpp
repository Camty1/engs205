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

template <typename T>
Array<T>::Array(size_t rows, size_t cols) : rows(rows), cols(cols) {
  data.resize(rows * cols);
}

template <typename T> T &Array<T>::operator()(size_t row, size_t col) {
  return data[col + row * cols];
}

template <typename T>
const T &Array<T>::operator()(size_t row, size_t col) const {
  return data[col + row * cols];
}

template <typename T> void Array<T>::operator+=(const Array<T> &B) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*this)(i, j) += B(i, j);
    }
  }
}

template <typename T> void Array<T>::operator-=(const Array<T> &B) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*this)(i, j) -= B(i, j);
    }
  }
}

template <typename T> void Array<T>::operator*=(const Array<T> &B) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*this)(i, j) *= B(i, j);
    }
  }
}

template <typename T> void Array<T>::operator/=(const Array<T> &B) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*this)(i, j) /= B(i, j);
    }
  }
}

template <typename T> const Array<T> Array<T>::operator+(const Array<T> &B) {
  Array<T> out(rows, cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      out(i, j) = (*this)(i, j) + B(i, j);
    }
  }

  return out;
}

template <typename T> const Array<T> Array<T>::operator-(const Array<T> &B) {
  Array<T> out(rows, cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      out(i, j) = (*this)(i, j) - B(i, j);
    }
  }

  return out;
}

template <typename T> const Array<T> Array<T>::operator*(const Array<T> &B) {
  Array<T> out(rows, cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      out(i, j) = (*this)(i, j) * B(i, j);
    }
  }

  return out;
}

template <typename T> const Array<T> Array<T>::operator/(const Array<T> &B) {
  Array<T> out(rows, cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      out(i, j) = (*this)(i, j) / B(i, j);
    }
  }

  return out;
}

template <typename T> void Array<T>::fill(const T &value) {
  std::fill(data.begin(), data.end(), value);
}

template <typename T> void Array<T>::print() const {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std::cout << std::setw(10) << std::setprecision(4) << (*this)(i, j)
                << " ";
    }
    std::cout << "\n";
  }
}

template <typename T> void Array<T>::write(std::string filename) const {
  std::ofstream file;
  file.open(filename);

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols - 1; j++) {
      file << (*this)(i, j) << ",";
    }
    file << (*this)(i, cols - 1) << std::endl;
  }

  file.close();
}
#endif
