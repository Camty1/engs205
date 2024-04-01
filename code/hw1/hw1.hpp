#include "utils_205.hpp"
#include <cblas.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

double calc_delta_s(double *x_elem, double *y_elem, int num_rows, int elem);

int calc_basis(double eta, double *phi);
