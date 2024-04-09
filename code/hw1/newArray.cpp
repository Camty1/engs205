#include <vector>

struct Slice {
  std::size_t start, stop, step;
  Slice(std::size_t start, std::size_t stop, std::size_t step = 1) : start(start), stop(stop), step(step) {};
};

template <typename T> class Array {
  private:
    std::vector<T> data;
    std::size_t rows, cols;

  public:
    Array(std::size_t rows, std::size_t cols) : rows(rows), cols(cols) {
      data.resize(rows * cols);
    }

    const T &operator()(std::size_t row, std::size_t col) {
      return data[row + col * rows];
    }

    void operator()(std::size_t row, std::size_t col) {
      data
    }

};
