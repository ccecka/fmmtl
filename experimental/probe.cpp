
#include <iostream>

#include "util/Probe.hpp"

#include "fmmtl/numeric/random.hpp"

int main() {
  using value_type = double;

  const unsigned N = 10;
  const unsigned M = 10;
  const unsigned r = 5;
  const double eps = 1e-10;

  flens::matrix<value_type> A(N,M);
#if 0
  A = 2.0;
#endif
#if 0
  fillRandom(A);
#endif
#if 1
  A(1,1) = 1;
  A(1,3) = 2;
  A(3,1) = 1;
  A(3,3) = 4;
#endif

  std::cout << "A(" << A.rows() << ", " << A.cols() << ") = " << A << std::endl;

  flens::matrix<double> U, V;

  std::tie(U, V) = probe_svd(A, r, eps);

  std::cout << "U(" << U.rows() << ", " << U.cols() << ") = " << U << std::endl;
  std::cout << "V(" << V.rows() << ", " << V.cols() << ") = " << V << std::endl;

  if (num_rows(U) != 0)
    std::cout << flens::matrix<double>(A - U*V) << std::endl;

  return 0;
}
