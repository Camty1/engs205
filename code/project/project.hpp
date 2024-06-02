#include <fftw3.h>
#include <cblas.h>
#include <lapacke.h>
#include <cmath>
#include <fstream>
#include <complex>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

void print_real_vector(int n, double* vec);
void write_real_vector(int n, double* vec, std::string filename);
void print_complex_vector(int n, std::complex<double>* vec);
