#pragma once

#include <utility>
#include <vector>
#include <algorithm>

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
 */
template <class MatrixIn,
          class MatrixOut = flens::matrix<typename MatrixIn::ElementType> >
std::tuple<MatrixOut,MatrixOut>
adaptive_cross_approx(const MatrixIn& A,
                      const double eps_tol, const unsigned max_rank) {
  using value_type = typename MatrixIn::ElementType;
  using size_type  = typename MatrixIn::IndexType;
  flens::Underscore<size_type> _;
  using std::abs;
  using flens::blas::dot;

  const double eps_tol_sq = eps_tol * eps_tol;

  const size_type n_rows = num_rows(A);
  const size_type n_cols = num_cols(A);

  // A stack of row indices to use
  std::vector<size_type> row_idx(n_rows);
  std::iota(row_idx.begin(), row_idx.end(), 0);
  std::random_shuffle(row_idx.begin(), row_idx.end());
  // A stack of col indices to use
  std::vector<size_type> col_idx(n_cols);
  std::iota(col_idx.begin(), col_idx.end(), 0);
  std::random_shuffle(col_idx.begin(), col_idx.end());

  /*  INITIALIZATION  */

  // Initialize the matrix and vector norms
  double matrix_norm = 0;
  double row_norm_sq;
  double col_norm_sq;

  // Matrices to construct
  size_type current_rank = 0;

  MatrixOut U(n_rows, max_rank);
  MatrixOut V(max_rank, n_cols);

  // Repeat till the desired tolerance is obtained
  do {

    auto&& row = V(current_rank, _);

    // Repeat until we find a good row
    while (true) {
      row = A(row_idx.back(), _);
      /// Row of the residuum and the pivot column
      for (size_type k = 0; k < current_rank; ++k)
        row -= U(row_idx.back(),k) * V(k,_);
      // Find the largest element of the row (larger than eps_tol)
      size_type pivot_idx = n_cols;
      double abs_pivot_val = eps_tol;
      for (size_type k = 0; k < n_cols; ++k) {
        if (abs_pivot_val < abs(row(k))) {
          abs_pivot_val = abs(row(k));
          pivot_idx = k;
        }
      }

      // Done with this row
      row_idx.pop_back();

      // Good row
      if (pivot_idx != n_cols) {
        // Normalization
        row *= 1 / row(pivot_idx);
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

    auto&& col = U(_, current_rank);

    // Repeat until we find a good col
    while (true) {
      col = A(_, col_idx.back());
      // Column of the residuum and the pivot row
      for (size_type k = 0; k < current_rank; ++k)
        col -= U(_,k) * V(k,col_idx.back());
      // Find the largest element of the col (larger than eps_tol)
      size_type pivot_idx = n_rows;
      double abs_pivot_val = eps_tol;
      for (size_type k = 0; k < n_rows; ++k) {
        if (abs_pivot_val < abs(col(k))) {
          abs_pivot_val = abs(col(k));
          pivot_idx = k;
        }
      }

      // Done with this column
      col_idx.pop_back();

      // Good col, bring the pivot to the front
      if (pivot_idx != n_rows) {
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

    /// New approximation of matrix norm
    row_norm_sq = dot(row, row);
    col_norm_sq = dot(col, col);

    matrix_norm += row_norm_sq * col_norm_sq;
    for (size_type k = 0; k < current_rank; ++k)
      matrix_norm += 2.0 * abs(dot(U(_,k),col)) * abs(dot(V(k,_),row));

    ++current_rank;

  } while (row_norm_sq * col_norm_sq > eps_tol_sq * matrix_norm * matrix_norm &&
           current_rank < max_rank &&
           !row_idx.empty() && !col_idx.empty());

return_statement:

  U.resize(n_rows, current_rank);
  V.resize(current_rank, n_cols);

  return std::make_tuple(U, V);
};
