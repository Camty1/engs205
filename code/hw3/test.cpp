#include "Array.hpp"
#include <cblas.h>
#include <cmath>
#include <complex>
#include <iostream>
#include <lapacke.h>

int main(int argc, char **argv) {
  Array<std::complex<double>> A(3, 3);
  Array<std::complex<double>> b(3);
  A.fill(std::complex<double>(0.0, 0.0));
  A(0, 0) = std::complex<double>(1.0, 0.0);
  A(1, 1) = std::complex<double>(0.0, 1.0);
  A(2, 2) = std::complex<double>(1.0, 0.0);

  b.fill(std::complex<double>(0.0, 0.0));
  b(1) = std::complex<double>(1.0, 0.0);
  b(2) = std::complex<double>(1.0, 1.0);


  lapack_int info;
  lapack_int *ipiv = new lapack_int[A.get_rows()];
  info = LAPACKE_zgesv(LAPACK_COL_MAJOR, A.get_rows(), 1, reinterpret_cast <__complex__ double*> (A.dataPtr()), A.get_rows(), ipiv, reinterpret_cast <__complex__ double*> (b.dataPtr()), b.get_rows());

  b.print();
}
