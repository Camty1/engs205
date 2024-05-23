#include "Array.hpp"
#include <fftw3.h>
#include <cblas.h>
#include <lapacke.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

double f(double x);
double dfdx(double x);
void fem_elem_matrix(double delta_x, double f1, double f2, Array<double>* A, Array<double>* b);
void hanning_window(int N, Array<double>* w);
