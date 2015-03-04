#pragma once

#include <complex>
#include <utility>
#include <vector>
#include <algorithm>
#include <random>

#include <fmmtl/numeric/flens.hpp>


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
  // A stack of col indices to use
  std::vector<size_type> col_idx(num_cols(A));
  std::iota(col_idx.begin(), col_idx.end(), U.firstRow());

  {
  std::random_device rd;
  std::mt19937 generator(rd());

  std::shuffle(row_idx.begin(), row_idx.end(), generator);
  std::shuffle(col_idx.begin(), col_idx.end(), generator);
  }

  size_type current_rank = U.firstRow();

  // Initialize the matrix and vector norms
  double matrix_norm = 0;
  double row_norm_sq;
  double col_norm_sq;

  // Repeat until the desired tolerance is obtained
  do {

    auto row = V(current_rank, _);

    // Repeat until we find a good row
    while (true) {
      // Initialize the row to the row of A
      row = A(row_idx.back(), _);
      // Row of the residuum and the pivot column
      for (size_type k = U.firstRow(); k < current_rank; ++k)
        row -= U(row_idx.back(),k) * V(k,_);
      // Done with this row
      row_idx.pop_back();

      // Find the largest element of the row
      size_type pivot_idx = iamax(row);

      // Good row
      if (abs(row(pivot_idx)) > eps_tol) {
        // Normalization
        row *= value_type(1) / row(pivot_idx);
        // Bring the pivot_idx to the back for processing
        auto cit = std::find(col_idx.begin(), col_idx.end(), pivot_idx);
        assert(cit != col_idx.end());
        //if (cit != col_idx.end())
        std::swap(*cit, col_idx.back());
        break;
      }

      // Bad row, try again if there's more
      if (row_idx.empty())
        goto return_statement;
    }


    auto col = U(_, current_rank);

    // Repeat until we find a good col
    while (true) {
      // Initialize the col to the col of A
      col = A(_, col_idx.back());
      // Column of the residuum and the pivot row
      for (size_type k = U.firstRow(); k < current_rank; ++k)
        col -= U(_,k) * V(k,col_idx.back());
      // Done with this column
      col_idx.pop_back();

      // Find the largest element of the col
      size_type pivot_idx = iamax(col);

      // Good col, bring the pivot to the front
      if (abs(col(pivot_idx)) > eps_tol) {
        // Bring the pivot to the back for processing if not already considered
        auto rit = std::find(row_idx.begin(), row_idx.end(), pivot_idx);
        if (rit != row_idx.end())
          std::swap(*rit, row_idx.back());
        break;
      }

      // Bad col, try again if there's more
      if (col_idx.empty())
        goto return_statement;
    }

    // New approximation of matrix norm
    row_norm_sq = real(dot(row,row));
    col_norm_sq = real(dot(col,col));

    matrix_norm += row_norm_sq * col_norm_sq;
    for (size_type k = U.firstRow(); k < current_rank; ++k)
      matrix_norm += 2.0 * abs(dot(U(_,k),col)) * abs(dot(V(k,_),row));

    ++current_rank;

  } while (row_norm_sq * col_norm_sq > eps_tol_sq * matrix_norm * matrix_norm &&
           current_rank <= U.lastCol() &&
           !row_idx.empty() && !col_idx.empty());

return_statement:

  auto r = _(U.firstRow(), current_rank-1);
  return std::make_tuple(U(_,r), V(r,_));
};
