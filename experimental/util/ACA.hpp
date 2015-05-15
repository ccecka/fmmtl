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
          class MatrixOut = typename MatrixIn::NoView>
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
  using flens::blas::asum;
  using flens::blas::iamax;

  const double eps_tol_sq = eps_tol * eps_tol;

  /*  INITIALIZATION  */

  // Matrices to construct
  MatrixOut U, V;
  U.reserve(num_rows(A), max_rank);
  V.reserve(num_cols(A), max_rank);

  // Scratch
  typename MatrixOut::Vector s;
  s.reserve(max_rank);

  // A stack of row indices to use
  std::vector<size_type> row_idx(num_rows(A));
  std::iota(row_idx.begin(), row_idx.end(), U.firstRow());

  size_type current_rank = U.firstCol();

  // Initialize the matrix and vector norms
  double matrix_normf_sq = 0;

  //std::cout << "NEW ACA " << num_rows(A) << "x" << num_cols(A) << std::endl;

  // Repeat until the desired tolerance is obtained
  while (true) {
    // Index range of the previous rows/cols
    auto prev = _(U.firstCol(), current_rank-1);

    // Construct an alias for the current row of V
    auto row = V(_, current_rank);

    // The pivot row and col
    size_type i, j;

    // Repeat until we find a good row
    do {
      if (row_idx.empty())
        goto return_statement;

      // Choose a random row from the remaining rows
      size_type k = fmmtl::random<size_type>::get(0,row_idx.size()-1);
      i = row_idx[k];
      row_idx[k] = row_idx.back();
      row_idx.pop_back();

      // Set the row to the residiuum of the row of A
      row = A(i,_) - U(i,prev) * transpose(V(_,prev));

      // Find the largest element of the row
      j = iamax(row);

      // If this row is effectively zero, try again
    } while (abs(row(j)) < 1e-14);

    // Normalization
    row *= value_type(1) / row(j);

    // Construct an alias for the current col of U
    auto col = U(_, current_rank);
    // Set the col to the residiuum of the col of A
    col = A(_,j) - U(_,prev) * V(j,prev);

    // Increment and break if we can't continue
    ++current_rank;
    if (current_rank > U.lastCol())
      break;

    // Compute the vector norms
    double rowcol_norm_sq = real(dot(row,row)) * real(dot(col,col));

    // If this update will make a small contribution to the matrix, stop.
    if (rowcol_norm_sq < eps_tol_sq * matrix_normf_sq)
      break;

    // Update the frobenius norm of U * V
    matrix_normf_sq += rowcol_norm_sq;
    s(prev) = transpose(U(_,prev)) * col;
    matrix_normf_sq += 2. * asum(s(prev));
    s(prev) = transpose(V(_,prev)) * row;
    matrix_normf_sq += 2. * asum(s(prev));
  } // end while

return_statement:

  auto r = _(U.firstCol(), current_rank-1);
  return std::make_tuple(U(_,r), V(_,r));
};
