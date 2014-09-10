#include "fmmtl/numeric/flens.hpp"


int main() {
  using T = double;

  // Define types from FLENS
  using namespace flens;
  using ZeroBased   = IndexOptions<int, 0>;
  //using matrix_type = GeMatrix<FullStorage<T, RowMajor, ZeroBased> >;
  using matrix_type = GeMatrix<FullStorage<T> >;
  //using vector_type = DenseVector<Array<T, ZeroBased> >;
  using vector_type = DenseVector<Array<T> >;
  const Underscore<typename matrix_type::IndexType>  _;

  unsigned n = 4;
  unsigned m = 5;
  unsigned r = std::min(n, m);

  matrix_type A(n, m);
  /*
  A = 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1,
      1, 1, 1, 1, 1,
      1, 1, 1, 1, 1;
  */
  /*
  A = 1, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 1, 0;
  */
  A = 1, 0, 0, 0, 2,
      0, 0, 3, 0, 0,
      0, 0, 0, 0, 0,
      0, 4, 0, 0, 0;
  std::cout << "A(" << A.rows() << ", " << A.cols() << ") = " << A << std::endl;

  matrix_type U(n,r), VT(r,m);
  vector_type D(r);
  flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::Save,
                     A, D, U, VT);

  const DiagMatrix<ConstArrayView<double> > DM(D);

  std::cout << "SVD Complete" << std::endl;

  std::cout << "A(" << A.rows() << ", " << A.cols() << ") = " << A << std::endl;

  std::cout << "U(" << U.rows() << ", " << U.cols() << ") = " << U << std::endl;

  std::cout << "D(..,..) = " << DM << std::endl;

  std::cout << "VT(" << VT.rows() << ", " << VT.cols() << ") = " << VT << std::endl;

  //std::cout << "U U^T = " << matrix_type(U * transpose(U)) << std::endl;
  //std::cout << "V V^T = " << matrix_type(transpose(VT) * VT) << std::endl;

  std::cout << matrix_type(U * matrix_type(DM * VT)) << std::endl;

  return 0;
}
