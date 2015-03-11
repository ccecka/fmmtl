
#include <iostream>

#include "util/ACA.hpp"

#include "fmmtl/numeric/random.hpp"

int main() {
  using value_type = double; //std::complex<double>;

  const unsigned N = 5;
  const unsigned M = 5;
  const unsigned r = 3;

  flens::matrix<value_type> A(N,M);
#if 1
  A = 2.0;
#endif
#if 0
  fillRandom(A);
#endif
#if 0
  A(1,1) = {1,2};
  A(1,3) = {2,1};
  A(3,1) = 1;
  A(3,3) = 4;
#endif

  std::cout << "A = \n" << A << std::endl;

  flens::matrix<value_type> U, V;

  std::tie(U, V) = adaptive_cross_approx(A, 1e-10, r);

  std::cout << "FINAL RANK = " << num_cols(U) << std::endl;
  std::cout << "U = \n" << U << std::endl;
  std::cout << "V = \n" << V << std::endl;
  std::cout << "UV = \n" << flens::matrix<value_type>(U*V) << std::endl;

  flens::matrix<value_type> Res;
  Res = A - U*V;

  std::cout << "Residual = \n" << Res << std::endl;
  std::cout << "norm_F = " << norm_f(Res) << std::endl;

  return 0;
}
