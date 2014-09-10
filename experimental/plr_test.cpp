
#include "plr.hpp"

int main() {
  constexpr unsigned n = 6;  //< Rows
  constexpr unsigned m = 7;  //< Cols

  // Data in the matrix
  double M[n*m] = {1, 1, 1, 2, 2, 2, 2,
                   1, 1, 1, 2, 2, 2, 2,
                   1, 1, 1, 2, 2, 2, 2,
                   1, 1, 1, 2, 2, 2, 2,
                   3, 3, 3, 4, 4, 4, 4,
                   3, 3, 3, 4, 4, 4, 4};

  // 2D points corresponding to matrix cols
  double s[2*m] = {0, 0,
                   1, 0,
                   0, 1,
                   100, 100,
                   101, 100,
                   100, 101,
                   101, 101};

  // 3D points corresponding to matrix rows
  double t[3*n] = {0, 0, 0,
                   1, 0, 0,
                   0, 1, 0,
                   0, 0, 1,
                   100, 0, 0,
                   101, 0, 0};

  // Construct the compressed matrix
  auto plr_mat = plr_compression<3,2>(M, n, m,
                                      t, s,
                                      1, 1e-10);


  // Use the compressed matrix to perform a product
  double x[m] = {1,
                 2,
                 3,
                 4,
                 5,
                 6,
                 7};

  double y[n] = {0,
                 0,
                 0,
                 0,
                 0,
                 0};

  double y_exact[n] = {0,
                       0,
                       0,
                       0,
                       0,
                       0};

  // y += A * x
  prod_acc(plr_mat, x, y);

  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      y_exact[i] += M[i*m+j] * x[j];

  // Print
  for (unsigned i = 0; i < n; ++i)
    std::cout << i << ":\t" << y[i] << "\t" << y_exact[i] << std::endl;

  // y += A * x
  prod_acc(plr_mat, x, y);

  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < m; ++j)
      y_exact[i] += M[i*m+j] * x[j];

  // Print
  for (unsigned i = 0; i < n; ++i)
    std::cout << i << ":\t" << y[i] << "\t" << y_exact[i] << std::endl;

  return 0;
}
