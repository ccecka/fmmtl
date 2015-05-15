#pragma once

#include <utility>
#include <cmath>

#include "fmmtl/numeric/flens.hpp"

/** Construct a low-rank approximation A ~= U * V^T
 * for maximum rank @a max_rank and error tolerance @a eps_tol at the rank.
 * @brief Uses a randomized "probing" SVD that approximates a true SVD with high
 * probability. The algorithm:
 *   - Input A (n x m)
 *   - Construct a random matrix R (m x rc) with rc < n.
 *   - Compute SVD of A^T R (n x rc) ~= U1 D1 V1T
 *   - Compute SVD of A  U1 (m x rc) ~= U2 D2 V2T
 *   - Then A ~= U2 D2 V2T U1T with D2 the approximate singular values of A.
 *
 * @param[in] A         Matrix of size n-by-m to be factorized
 * @param[in] max_rank  Maximum rank of the resulting factorizaton
 * @param[in] eps_tol   Maximum magnitude of the @a max_rank-th singular value.
 * @returns Tuple of <DenseMatrix, DenseMatrix>:
 *    DenseMatrix:  (U2 D2)   matrix of size n-by-r, where r <= max_rank.
 *    DenseMatrix:  (V2T U1T) matrix of size r-by-m, where r <= max_rank.
 *
 * @post num_rows(std::get<0>(result)) == 0 if no such factorization.
 *
 * @note Only uses matrix-matrix products of A and transpose(A).
 */
template <class MatrixIn,
          class MatrixOut = typename MatrixIn::NoView>
std::tuple<MatrixOut,MatrixOut>
probe_svd(const MatrixIn& A, const unsigned max_rank, const double eps_tol) {
  using T = typename MatrixIn::ElementType;

  using size_type  = typename MatrixOut::IndexType;
  flens::Underscore<size_type> _;

  using namespace flens;
  // XXX: Fix with Allocator
  using VectorType = DenseVector<Array<double> >;

  // A is an n x m matrix
  unsigned n = num_rows(A);
  unsigned m = num_cols(A);
  unsigned rc = std::min(max_rank + 10, std::min(n,m));

  // Construct a random matrix
  MatrixOut R1(n, rc);
  fillRandom(R1);

  // Factor A^T * R1 (which is m x rc)
  // TODO: Use an over-write?
  MatrixOut  U1(m,rc), VT(rc,rc);
  VectorType D(rc);
  flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::None,
                     MatrixOut(conjTrans(A) * R1),
                     D, U1, VT);

  // Factor A * U1 (which is n x rc)
  // Reuse VT and D
  MatrixOut U2(n,rc);
  flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::Save,
                     MatrixOut(A * U1),
                     D, U2, VT);

  // Find the eps-rank
  while (D(rc)/D(1) < eps_tol) --rc;

  if (rc <= max_rank) {
    auto r = _(1,rc);
    const DiagMatrix<ConstArrayView<double> > DM = D(r);

    return std::make_tuple(U2(_,r) * DM,
                           VT(r,r) * conjTrans(U1(_,r)));
  }

  return std::make_tuple(MatrixOut(), MatrixOut());
}
