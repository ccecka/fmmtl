/** Simple C-library functions for PLR matrices
 */

#include "PLR.hpp"
#include <complex>

extern "C" {

#include "PLR.h"

// An opaque type that we'll use as a handle
struct ZPLR2D_Handle {
  using Tree = fmmtl::NDTree<2>;
  using PLR  = flens::PLR_Matrix<std::complex<double>, Tree, Tree>;
};


/** Complex general matrix plr over two-dimension domains */
ZPLR2D_Handle*
zplr2d(char order, double _Complex* data, int n, int m, int lda,
       const double* trgs, const double* srcs,
       unsigned max_rank, double eps_tol, unsigned init_depth)
{
  using PLR = typename ZPLR2D_Handle::PLR;

  std::complex<double>* raw_d = reinterpret_cast<std::complex<double>*>(data);
  PLR* plr = new PLR(geplr<2,2>(order, raw_d, n, m, lda, trgs, srcs,
                                max_rank, eps_tol, init_depth));

  return reinterpret_cast<ZPLR2D_Handle*>(plr);
}

/** y = beta*y + alpha*M*x */
void
zplr2dmv(char trans, const ZPLR2D_Handle* M, int n, int m,
         const double _Complex alpha, const double _Complex* x, int incx,
         const double _Complex beta, double _Complex* y, int incy)
{
  using PLR = typename ZPLR2D_Handle::PLR;

  // Wrap the pointers into vectors
  using namespace flens;
  using T   = std::complex<double>;
  using AV  = ArrayView<T>;
  using CAV = ConstArrayView<T>;

  DenseVector<CAV> xv = CAV(m, reinterpret_cast<const T*>(x), incx);
  DenseVector<AV>  yv =  AV(n, reinterpret_cast<      T*>(y), incy);

  blas::mv(cxxblas::getCxxBlasEnum<cxxblas::Transpose>(trans),
           reinterpret_cast<const T&>(alpha), *reinterpret_cast<const PLR*>(M),
           xv, reinterpret_cast<const T&>(beta), yv);
}

/** Destroy the PLR matrix pointed to by the handle */
void
free_zplr2d(ZPLR2D_Handle* plr)
{
  using PLR = typename ZPLR2D_Handle::PLR;

  if (plr)
    delete reinterpret_cast<PLR*>(plr);
}

} // extern C
