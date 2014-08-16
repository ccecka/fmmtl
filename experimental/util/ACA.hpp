#pragma once

#include <utility>
#include <vector>
#include <algorithm>

#define MTL_WITH_AUTO
#define MTL_WITH_DEFAULTIMPL
#define MTL_WITH_INITLIST
#define MTL_WITH_MOVE
#define MTL_WITH_RANGEDFOR
#define MTL_WITH_STATICASSERT
#define MTL_WITH_TEMPLATE_ALIAS
#define MTL_WITH_VARIADIC_TEMPLATE
#include <boost/numeric/mtl/mtl.hpp>


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
          class MatrixOut = mtl::dense2D<typename MatrixIn::value_type> >
std::tuple<MatrixOut,MatrixOut>
adaptive_cross_approx(const MatrixIn& A,
                      const double eps_tol, const unsigned max_rank) {
  using value_type = typename MatrixIn::value_type;
  using size_type  = typename MatrixIn::size_type;
  auto _ = mtl::iall;
  using std::abs;

  // A stack of row indices to use
  std::vector<size_type> row_idx(num_rows(A));
  std::iota(row_idx.begin(), row_idx.end(), 0);
  std::random_shuffle(row_idx.begin(), row_idx.end());
  // A stack of col indices to use
  std::vector<size_type> col_idx(num_cols(A));
  std::iota(col_idx.begin(), col_idx.end(), 0);
  std::random_shuffle(col_idx.begin(), col_idx.end());

  /*  INITIALIZATION  */

  /// Initialize the matrix norm and the the first row index
  double matrix_norm = 0;
  double row_norm;
  double col_norm;
  value_type pivot_val;

  size_type current_rank = 0;

  MatrixOut U(num_rows(A), max_rank);
  MatrixOut V(max_rank, num_cols(A));

  /// Repeat till the desired tolerance is obtained
  do {

    auto&& row = V[current_rank][_];

    // Repeat until we find a good row
    while (true) {
      row = A[row_idx.back()][_];
      /// Row of the residuum and the pivot column
      for (size_type k = 0; k < current_rank; ++k)
        row -= U[row_idx.back()][k] * V[k][_];
      // Find the largest element of the row
      double abs_pivot_val = eps_tol;
      size_type pivot_idx = num_cols(row);
      for (size_type k = 0; k < num_cols(row); ++k) {
        if (abs_pivot_val < abs(row(k))) {
          abs_pivot_val = abs(row(k));
          pivot_idx = k;
        }
      }

      // Done with this row
      row_idx.pop_back();

      // Good row
      if (pivot_idx != num_cols(row)) {
        // Bring the pivot_idx to the back for processing
        auto cit = std::find(col_idx.begin(), col_idx.end(), pivot_idx);
        if (cit != col_idx.end())
          std::swap(*cit, col_idx.back());
        break;
      }

      // Bad row, try again if there's more
      if (row_idx.empty())
        goto return_statement;
    }

    pivot_val = V[current_rank][col_idx.back()];
    value_type gamma = 1.0 / pivot_val;

    auto&& col = U[_][current_rank];

    // Repeat until we find a good col
    while (true) {
      col = A[_][col_idx.back()];
      /// Column of the residuum and the pivot row
      for (size_type k = 0; k < current_rank; ++k)
        col -= U[_][k] * V[k][col_idx.back()];
      /// Normalizing constant
      col *= gamma;
      // Find the largest element of the col
      double abs_pivot_val = eps_tol;
      size_type pivot_idx = num_rows(col);
      for (size_type k = 0; k < num_cols(row); ++k) {
        if (abs_pivot_val < abs(row(k))) {
          abs_pivot_val = abs(row(k));
          pivot_idx = k;
        }
      }

      // Done with this column
      col_idx.pop_back();

      // Good col, bring the pivot to the front
      if (pivot_idx != num_rows(col)) {
        // Bring the pivot to the back for processing
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
    row_norm = two_norm(row);
    col_norm = two_norm(col);

    matrix_norm += row_norm * row_norm * col_norm * col_norm;
    for (size_type k = 0; k < current_rank; ++k)
      matrix_norm += 2.0 * abs(dot(U[_][k],col)) * abs(dot(V[k][_],row));

    ++current_rank;

  } while (row_norm * col_norm > abs(pivot_val) * eps_tol * matrix_norm &&
           current_rank < max_rank &&
           !row_idx.empty() && !col_idx.empty());

return_statement:

  U.change_dim(num_rows(A), current_rank, /* keep_elements: */ true);
  V.change_dim(current_rank, num_cols(A), /* keep_elements: */ true);

  return std::make_tuple(U, V);
};
