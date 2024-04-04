#include "Array.hpp"

int main(int argc, char **argv) {
  Array<double> A(3, 3);
  Array<double> b(3);
  A.fill(0.0);
  A(0, 0) = 1.0;
  A(1, 1) = 1.0;
  A(2, 2) = 1.0;

  b(0) = 1;
  b(1) = 2;
  b(2) = 3;

  A.print();
  b.print();
}
