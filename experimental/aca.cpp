
#include <iostream>

#include "util/ACA.hpp"

#include "fmmtl/numeric/random.hpp"

int main() {
  using value_type = double;

  mtl::dense2D<value_type> A(5,5);
#if 0
  A = 2.0;
#endif
#if 1
  for (unsigned i = 0; i < num_rows(A); ++i)
    for (unsigned j = 0; j < num_cols(A); ++j)
      A(i,j) = fmmtl::random<value_type>::get();
#endif
#if 0
  A[1][1] = 1;
  A[1][3] = 2;
  A[3][1] = 1;
  A[3][3] = 4;
#endif

  std::cout << "A = \n" << A << std::endl;

  mtl::dense2D<double> U, V;

  std::tie(U, V) = adaptive_cross_approx(A, 1e-10, 1);

  std::cout << "FINAL RANK = " << num_cols(U) << std::endl;
  std::cout << "U = \n" << U << std::endl;
  std::cout << "V = \n" << V << std::endl;
  std::cout << "UV = \n" << U*V << std::endl;

  mtl::dense2D<value_type> R;
  R = A - U*V;

  std::cout << "Residual = \n" << R << std::endl;
  std::cout << "norm_F = " << frobenius_norm(R) << std::endl;

  return 0;
}
