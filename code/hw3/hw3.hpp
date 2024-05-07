#ifndef HW3_HPP
#define HW3_HPP
#include "Array.hpp"
#include <cblas.h>
#include <lapacke.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

double calc_alpha(double x_behind, double y_behind, double x, double y, double x_ahead, double y_ahead);
std::complex<double> hankel1(double nu, double x);
int factorial(int n);
int factorial2(int n);
double struve(double nu, double x);
std::complex<double> int_hankel1(double x);

#endif
