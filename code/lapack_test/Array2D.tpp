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

template <typename T> T *Array2D<T>::dataPtr() { return data.data(); }

template <typename T> size_t Array2D<T>::getRows() const { return rows; }

template <typename T> size_t Array2D<T>::getCols() const { return cols; }

// Implementation of fillRange
template <typename T>
void Array2D<T>::fill(size_t startRow, size_t startCol, size_t endRow,
                           size_t endCol, const T &value) {
  for (size_t col = startCol; col <= endCol && col < cols; ++col) {
    for (size_t row = startRow; row <= endRow && row < rows; ++row) {
      (*this)(row, col) = value;
    }
  }
}

// Implementation of the print method
template <typename T> void Array2D<T>::print() const {
  for (size_t row = 0; row < rows; ++row) {
    for (size_t col = 0; col < cols; ++col) {
      std::cout << std::setw(10) << std::setprecision(4) << (*this)(row, col)
                << " ";
    }
    std::cout << "\n";
  }
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
