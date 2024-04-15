#include "Array.hpp"

template <typename T>
Array<T>::Array(size_t rows, size_t cols)
    : data(rows * cols), rows(rows), cols(cols) {}

template <typename T> T &Array<T>::operator()(size_t row, size_t col) {
  return data[row + col * rows];
}

template <typename T>
const T &Array<T>::operator()(size_t row, size_t col) const {
  return data[row + col * rows];
}

template <typename T>
void Array<T>::check_slice(Slice slice) {
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  if (start < 0 || start > rows - 1) {
    throw std::range_error("Start out of range");
  }
  if (stop < 0 || stop > rows - 1) {
    throw std::range_error("Stop out of range");
  }
  if (start > stop && step > 0) {
    throw std::invalid_argument("Start is below stop and step isn't negative");
  }
  if (step == 0) {
    throw std::invalid_argument("Step can't be 0");
  }
  
}


template <typename T>
const Array<T> Array<T>::get_slice(Slice slice) {
  if (rows != 1 && cols != 1) {
    throw std::invalid_argument("1D slice given to 2D array");
  }

  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop - start + 1) / step;

  Array<T> out(size);

  for (size_t i = 0; i < size; i++) {
    out(i) = (*this)(start + step * i);
  }

  return out;
}

template <typename T>
const Array<T> Array<T>::get_slice(Slice row_slice, Slice col_slice) {

  check_slice(row_slice);
  size_t row_start = (row_slice.start == -1) ? 0 : row_slice.start;
  size_t row_stop = (row_slice.stop == -1) ? rows - 1 : row_slice.stop;
  size_t row_step = row_slice.step;

  size_t row_size = (row_stop >= row_start) ? (row_stop - row_start + 1) / row_step : (row_stop - row_start - 1) / row_step;

  check_slice(col_slice);
  size_t col_start = (col_slice.start == -1) ? 0 : col_slice.start;
  size_t col_stop = (col_slice.stop == -1) ? cols - 1 : col_slice.stop;
  size_t col_step = col_slice.step;

  size_t col_size = (col_stop >= col_start) ? (col_stop - col_start + 1) / col_step : (col_stop - col_start - 1) / col_step;

  Array<T> out(row_size, col_size);
  
  for (size_t i = 0; i < row_size; i++) {
    for (size_t j = 0; j < col_size; j++) {
      out(i, j) = (*this)(row_start + i * row_step, col_start + j * col_step);
    }
  }

  return out;
}

template <typename T>
const Array<T> Array<T>::get_slice(size_t row, Slice slice) {
  if (row > rows - 1 || rows < 0) {
    throw std::range_error("Invalid row given");
  }

  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? cols - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop >= start) ? (stop - start + 1) / step : (stop - start - 1) / step;

  Array<T> out(1, size);
  
  for (size_t j = 0; j < size; j++) {
    out(0, j) = (*this)(row, start + j * step);
  }

  return out;
}

template <typename T>
const Array<T> Array<T>::get_slice(Slice slice, size_t col) {
  if (col > cols - 1 || cols < 0) {
    throw std::range_error("Invalid col given");
  }

  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop >= start) ? (stop - start + 1) / step : (stop - start - 1) / step;

  Array<T> out(size);
  
  for (size_t i = 0; i < size; i++) {
    out(i) = (*this)(start + i * step, col);
  }

  return out;
}

template <typename T> void Array<T>::set_slice(const Array<T> &A, Slice slice) {
  if (cols != 1) {
    throw std::invalid_argument("1D slice passed to 2D Array");
  }
  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop >= start) ? (stop - start + 1) / step : (stop - start - 1) / step;

  if (size != A.get_rows()) {
    throw std::out_of_range("Size of slice does not match that of Array");
  }

  for (size_t i = 0; i < size; i++) {
    (*this)(start + i * step) = A(i);
  }
}

template <typename T> void Array<T>::set_slice(const Array<T> &A, Slice row_slice, Slice col_slice) {

  check_slice(row_slice);
  size_t row_start = (row_slice.start == -1) ? 0 : row_slice.start;
  size_t row_stop = (row_slice.stop == -1) ? rows - 1 : row_slice.stop;
  size_t row_step = row_slice.step;

  size_t row_size = (row_stop >= row_start) ? (row_stop - row_start + 1) / row_step : (row_stop - row_start - 1) / row_step;

  check_slice(col_slice);
  size_t col_start = (col_slice.start == -1) ? 0 : col_slice.start;
  size_t col_stop = (col_slice.stop == -1) ? cols - 1 : col_slice.stop;
  size_t col_step = col_slice.step;

  size_t col_size = (col_stop >= col_start) ? (col_stop - col_start + 1) / col_step : (col_stop - col_start - 1) / col_step;

  if (A.get_rows() != row_size || A.get_cols() != col_size) {
    throw std::out_of_range("Given shape does not match slice shape");
  }

  for (size_t i = 0; i < row_size; i++) {
    for (size_t j = 0; j < col_size; j++) {
      (*this)(row_start + i * row_step, col_start + j * col_step) = A(i,j);
    }
  }
}

template <typename T> void Array<T>::set_slice(const Array<T> &A, size_t row, Slice slice) {
  if (row < 0 || row >= rows) {
    throw std::out_of_range("Row given is invalid");
  }
  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop >= start) ? (stop - start + 1) / step : (stop - start - 1) / step;

  if (size != A.get_cols()) {
    throw std::invalid_argument("Size of slice does not match that of Array");
  }

  for (size_t i = 0; i < size; i++) {
    (*this)(row, start + i * step) = A(0, i);
  }
}

template <typename T> void Array<T>::set_slice(const Array<T> &A, Slice slice, size_t col) {
  if (col < 0 || col >= rows) {
    throw std::out_of_range("Col given is invalid");
  }
  check_slice(slice);
  size_t start = (slice.start == -1) ? 0 : slice.start;
  size_t stop = (slice.stop == -1) ? rows - 1 : slice.stop;
  size_t step = slice.step;

  size_t size = (stop >= start) ? (stop - start + 1) / step : (stop - start - 1) / step;

  if (size != A.get_rows()) {
    throw std::invalid_argument("Size of slice does not match that of Array");
  }

  for (size_t i = 0; i < size; i++) {
    (*this)(start + i * step, col) = A(i);
  }
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
