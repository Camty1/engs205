#include <cblas.h>
#include <iostream>

int main(int argc, char **argv) {
  int n = 3;
  double *A = new double[n];
  double *B = new double[n];

  for (int i = 0; i < n; i++) {
    A[i] = i + 1;
    B[i] = n - i;
  }

  double C = 0.0;

  C = cblas_ddot(n, A, 1, B, 1);

  std::cout << C << "\n";
  return 0;
}
