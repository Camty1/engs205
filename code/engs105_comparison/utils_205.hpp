#ifndef utils_205
#define utils_205
#include <fstream>
#include <iostream>

int general_matrix_index(char order, int num_rows, int num_cols, int row,
                         int col);
int banded_matrix_index(char order, int num_rows, int num_cols, int bandwidth,
                        int row, int col);
int fill(char order, double *matrix, int num_rows, int num_cols, int row_start,
         int row_stop, int col_start, int col_stop, double value);
int fill(double *vector, int num_rows, int start, int stop, double value);
void print_matrix(char order, double *matrix, int num_rows, int num_cols);
void print_matrix(char order, int *matrix, int num_rows, int num_cols);
void write_matrix(char order, double *matrix, int num_rows, int num_cols,
                  std::string filename);
void write_matrix(char order, int *matrix, int num_rows, int num_cols,
                  std::string filename);
void print_vector(double *vector, int num_rows);
void print_vector(int *vector, int num_rows);
void write_vector(double *vector, int num_rows, std::string filename);
void write_vector(int *vector, int num_rows, std::string filename);
#endif
