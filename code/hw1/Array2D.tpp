// Array2D.cpp

template <typename T>
Array2D<T>::Array2D(size_t rows, size_t cols)
    : data(rows * cols), rows(rows), cols(cols) {}

template <typename T> T &Array2D<T>::operator()(size_t row, size_t col) {
  return data[col * rows + row];
}

template <typename T>
const T &Array2D<T>::operator()(size_t row, size_t col) const {
  return data[col * rows + row];
}

template <typename T>
Array2D<T> Array2D<T>::operator-() const {
  Array2D<T> out(rows, cols);
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      out(i, j) = -(*this)(i, j);
    }
  }
  return out;
}

template <typename T> T *Array2D<T>::dataPtr() { return data.data(); }

template <typename T> size_t Array2D<T>::getRows() const { return rows; }

template <typename T> size_t Array2D<T>::getCols() const { return cols; }

template <typename T>
void Array2D<T>::fillRange(size_t startRow, size_t startCol, size_t endRow,
                           size_t endCol, const T &value) {
  for (size_t col = startCol; col <= endCol && col < cols; ++col) {
    for (size_t row = startRow; row <= endRow && row < rows; ++row) {
      (*this)(row, col) = value;
    }
  }
}

template <typename T>
void Array2D<T>::fill(const T &value) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*this)(i, j) = value;
    }
  }
}

template <typename T> void Array2D<T>::print() const {
  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      std::cout << std::setw(10) << std::setprecision(4) << (*this)(row, col)
                << " ";
    }
    std::cout << "\n";
  }
}

template <typename T> void Array2D<T>::write(std::string filename) const {
  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      file << (*this)(i, j) << ',';
    }
    file << std::endl;
  }
  file.close();
}

template<typename T>
void Array2D<T>::insert(const Array2D<T>& source,
                              size_t sourceStartRow, size_t sourceStartCol, 
                              size_t subsetRows, size_t subsetCols, 
                              size_t targetStartRow, size_t targetStartCol) {
    for (size_t r = 0; r < subsetRows; ++r) {
        for (size_t c = 0; c < subsetCols; ++c) {
            size_t sourceRow = sourceStartRow + r;
            size_t sourceCol = sourceStartCol + c;
            size_t targetRow = targetStartRow + r;
            size_t targetCol = targetStartCol + c;

            // Check bounds for both source and target arrays
            if (sourceRow < source.getRows() && sourceCol < source.getCols() &&
                targetRow < this->getRows() && targetCol < this->getCols()) {
                (*this)(targetRow, targetCol) = source(sourceRow, sourceCol);
            }
        }
    }
}

template<typename T>
std::vector<T> Array2D<T>::getRowVec(size_t row) const {
  std::vector<T> out(cols);
  for (int j = 0; j < cols; j++) {
    out[j] = (*this)(row, j);
  }
  return out;
}
template<typename T>
std::vector<T> Array2D<T>::getColVec(size_t col) const {
  std::vector<T> out(rows);
  for (int i = 0; i < rows; i++) {
    out[i] = (*this)(i, col);
  }
  return out;
}
