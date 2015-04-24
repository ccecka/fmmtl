#pragma once

#include <complex>
#include <utility>
#include <vector>
#include <algorithm>

#include "fmmtl/numeric/random.hpp"
#include "fmmtl/numeric/flens.hpp"


/** Partial pivoted LU to construct low-rank.
 * @brief Obtains the low-rank decomposition of the matrix to a desired tolerance
 * using the partial pivoting LU algorithm. Given a matrix 'A' and tolerance
 * 'eps_tol', computes matrices 'U' and 'V' such that ||A-UV||_F < epsilon.
 * The norm is Frobenius norm.
 *
 * @param[in] ...
 *
 * @param[out] ...
 *
 * TODO: Could change the return type to yield a compressed permuted LU like sv
 */
template <class MatrixIn,
          class MatrixOut = flens::GeMatrix<flens::FullStorage<typename MatrixIn::ElementType> > >
std::tuple<MatrixOut,MatrixOut>
adaptive_cross_approx(const MatrixIn& A,
                      const double eps_tol,
                      typename MatrixIn::IndexType max_rank) {
  using value_type = typename MatrixIn::ElementType;
  using size_type  = typename MatrixIn::IndexType;
  flens::Underscore<size_type> _;
  using std::abs;
  using std::real;
  using flens::blas::dot;
  using flens::blas::iamax;

  const double eps_tol_sq = eps_tol * eps_tol;

  /*  INITIALIZATION  */

  // Matrices to construct
  MatrixOut U(num_rows(A), max_rank);
  MatrixOut V(max_rank, num_cols(A));

  // A stack of row indices to use
  std::vector<size_type> row_idx(num_rows(A));
  std::iota(row_idx.begin(), row_idx.end(), U.firstRow());

  size_type current_rank = U.firstRow();

  // Initialize the matrix and vector norms
  double matrix_normf_sq = 0;
  double row_norm_sq;
  double col_norm_sq;

  //std::cout << "NEW ACA " << num_rows(A) << "x" << num_cols(A) << std::endl;

  // Repeat until the desired tolerance is obtained
  while(true) {
    // Index range of the previous rows/cols
    auto prev = _(U.firstRow(), current_rank-1);

    // Construct an alias for the current row of V
    auto row = V(current_rank, _);
    // Construct an alias for the current col of U
    auto col = U(_, current_rank);

    // The pivot row and col
    size_type i, j;

    // Repeat until we find a good row
    do {
      if (row_idx.empty())
        goto return_statement;

      // Choose a random row from the remaining rows
      size_type k = fmmtl::random<size_type>::get(0,row_idx.size()-1);
      i = row_idx[k];
      // Remove this row from consideration in the future
      row_idx[k] = row_idx.back();
      row_idx.pop_back();

      // Set the row to the residiuum of the row of A
      row = A(i,_) - U(i,prev) * V(prev,_);

      // Find the largest element of the row
      j = iamax(row);

      // If this row is effectively zero, try again
    } while (abs(row(j)) < 1e-13);

    // Normalization
    row *= value_type(1) / row(j);

    // Set the col to the residiuum of the col of A
    col = A(_,j) - U(_,prev) * V(prev,j);

    // Increment and break if we can't continue
    ++current_rank;
    if (current_rank > U.lastCol())
      break;

    // Compute the vector norms
    row_norm_sq = real(dot(row,row));
    col_norm_sq = real(dot(col,col));

    // If this update will make a small contribution to the matrix, stop.
    if (row_norm_sq * col_norm_sq < eps_tol_sq * matrix_normf_sq)
      break;

    // Update the frobenius norm of U * V
    matrix_normf_sq += row_norm_sq * col_norm_sq;
    for (size_type k = U.firstRow(); k < current_rank-1; ++k)
      matrix_normf_sq += 2 * abs(dot(U(_,k),col)) * abs(dot(V(k,_),row));
  } // end while

return_statement:

  auto r = _(U.firstRow(), current_rank-1);
  return std::make_tuple(U(_,r), V(r,_));
};
