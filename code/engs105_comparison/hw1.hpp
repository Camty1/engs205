#include "Array2D.hpp"
#include <cblas.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <lapacke.h>
#include <sstream>
#include <string>
#include <vector>

double calc_delta_s(double *x_elem, double *y_elem, int num_rows, int elem);

int calc_basis(double eta, double *phi);

template <typename T> void print_vector(std::vector<T> vec);
template <typename T>
void write_vector(std::vector<T> vec, std::string filename);
template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b);
template <typename T>
void operator+=(std::vector<T> &a, const std::vector<T> &b);
template <typename T>
void operator-=(std::vector<T> &a, const std::vector<T> &b);
template <typename T>
std::vector<T> operator*(const std::vector<T> &vec, const T &scalar);
template <typename T>
std::vector<T> operator*(const T &scalar, const std::vector<T> &vec);
