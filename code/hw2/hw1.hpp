#ifndef HW1_HPP
#define HW1_HPP
#include "Array2D.hpp"
#include <cblas.h>
#include <lapacke.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <lapacke.h>
#include <sstream>
#include <string>
#include <vector>

double calc_alpha(double x_behind, double y_behind, double x, double y, double x_ahead, double y_ahead);

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
#endif
