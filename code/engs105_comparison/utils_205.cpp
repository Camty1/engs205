#include "utils_205.hpp"

int general_matrix_index(char order, int num_rows, int num_cols, int row,
                         int col) {

  if (row < 0 || row >= num_rows) {
    std::cerr << "Invalid row: " << row;
    return EXIT_FAILURE;
  }
  if (col < 0 || col >= num_cols) {
    std::cerr << "Invalid column: " << col;
    return EXIT_FAILURE;
  }
  if (order == 'c') {
    return (col - 1) * num_rows + row;

  } else if (order == 'r') {
    return (row - 1) * num_cols + col;
  } else {
    std::cerr << "Invalid storage order ([r]ow or [c]olumn): " << order;
    return EXIT_FAILURE;
  }
}

int banded_matrix_index(char order, int num_rows, int num_cols, int bandwidth,
                        int row, int col) {
  std::cerr << "banded_matrix_index not implemented.";
  return EXIT_FAILURE;
}

int fill(char order, double *matrix, int num_rows, int num_cols, int row_start,
         int row_stop, int col_start, int col_stop, double value) {
  if (row_start < 0 || row_start >= num_rows) {
    std::cerr << "Invalid row start: " << row_start;
    return EXIT_FAILURE;
  }
  if (row_stop < 0 || row_stop >= num_rows || row_stop < row_start) {
    std::cerr << "Invalid row stop: " << row_stop;
    return EXIT_FAILURE;
  }
  if (col_start < 0 || col_start >= num_cols) {
    std::cerr << "Invalid col start: " << col_start;
    return EXIT_FAILURE;
  }
  if (col_stop < 0 || col_stop >= num_cols || col_stop < col_start) {
    std::cerr << "Invalid col stop: " << col_stop;
    return EXIT_FAILURE;
  }
  if (order == 'c' || order == 'r') {
    for (int r = row_start; r <= row_stop; r++) {
      for (int c = col_start; c <= col_stop; c++) {
        int idx = general_matrix_index(order, num_rows, num_cols, r, c);
        matrix[idx] = value;
      }
    }
    return EXIT_SUCCESS;
  } else {
    std::cerr << "Invalid storage order ([r]ow or [c]olumn): " << order;
    return EXIT_FAILURE;
  }
}

int fill(double *vector, int num_rows, int start, int stop, double value) {
  if (start < 0 || start >= num_rows) {
    std::cerr << "Invalid start: " << start << std::endl;
    return EXIT_FAILURE;
  }
  if (stop < 0 || stop >= num_rows || stop < start) {
    std::cerr << "Invalid stop: " << stop << std::endl;
    return EXIT_FAILURE;
  }

  for (int i = start; i <= stop; i++) {
    vector[i] = value;
  }
  return EXIT_SUCCESS;
}
void print_matrix(char order, double *matrix, int num_rows, int num_cols) {
  if (order == 'r' || order == 'c') {
    std::cout.precision(4);
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_cols; j++) {
        int idx = general_matrix_index(order, num_rows, num_cols, i, j);
        std::cout.width(10);
        std::cout << matrix[idx];
      }
      std::cout << "\n";
    }
  } else {
    std::cerr << "Invalid storage order([r]ow or [c]olumn): " << order;
  }
}

void print_matrix(char order, int *matrix, int num_rows, int num_cols) {
  if (order == 'r' || order == 'c') {
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_cols; j++) {
        int idx = general_matrix_index(order, num_rows, num_cols, i, j);
        std::cout.width(10);
        std::cout << matrix[idx];
      }
      std::cout << "\n";
    }
  } else {
    std::cerr << "Invalid storage order([r]ow or [c]olumn): " << order;
  }
}

void write_matrix(char order, double *matrix, int num_rows, int num_cols,
                  std::string filename) {

  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      int idx = general_matrix_index(order, num_rows, num_cols, i, j);
      file << matrix[idx] << ',';
    }
    file << std::endl;
  }

  file.close();
}

void write_matrix(char order, int *matrix, int num_rows, int num_cols,
                  std::string filename) {

  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      int idx = general_matrix_index(order, num_rows, num_cols, i, j);
      file << matrix[idx] << ',';
    }
    file << std::endl;
  }

  file.close();
}

void print_vector(double *vector, int num_rows) {
  std::cout.precision(4);
  for (int i = 0; i < num_rows; i++) {
    std::cout.width(10);
    std::cout << vector[i] << std::endl;
  }
}

void print_vector(int *vector, int num_rows) {
  for (int i = 0; i < num_rows; i++) {
    std::cout.width(10);
    std::cout << vector[i] << std::endl;
  }
}

void write_vector(double *vector, int num_rows, std::string filename) {
  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < num_rows; i++) {
    file << vector[i] << std::endl;
  }

  file.close();
}

void write_vector(int *vector, int num_rows, std::string filename) {
  std::ofstream file;
  file.open(filename);

  for (int i = 0; i < num_rows; i++) {
    file << vector[i] << std::endl;
  }

  file.close();
}
