#include "Array2D.hpp"
#include <cblas.h>
#include <iostream>
#include <lapacke.h>

int main(int argc, char **argv) {
  int n = 3;
  Array2D<double> A(n, n);
  std::vector<double> b(n);
  Array2D<double> C(2, 2);
  C.fill(0, 0, 1, 1, 1);
  A(0, 0) = 1;
  A.insert(C, 0, 0, 2, 2, 0, 1);
  A(2, 2) = 1;
  A.print();
  b[0] = 1;
  b[1] = 2;
  b[2] = 3;
  for (int i = 0; i < n; i++) {
    std::cout << b[i] << std::endl;
  }
  lapack_int info;
  lapack_int *ipiv = new lapack_int[A.getRows()];
  info = LAPACKE_dgesv(LAPACK_COL_MAJOR, A.getRows(), 1, A.dataPtr(),
                       A.getRows(), ipiv, b.data(), n);
  for (int i = 0; i < n; i++) {
    std::cout << b[i] << std::endl;
  }

  return 0;
}
