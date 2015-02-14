#include "plr.hpp"

using namespace flens;

/** A type representing an implicit matrix
 * In this example, we represent the matrix
 * | 2 1 1 1 ..
 * | 1 2 1 1 ..
 * | 1 1 2 1 ..
 * . . . . .
 * of any size.
 */
struct MyMatrix {
  //! ElementType is the field this matrix is defined over
  using ElementType = double;

  /** Compute and return the explicit matrix that M(rows,cols) represents
   * @param[in] rows  A subset of the rows of this matrix
   * @param[in] cols  A subset of the cols of this matrix
   *
   * @result A dense matrix indexed from [1,rows.length()] x [1,cols.length()]
   */
  template <typename Rows, typename Cols>
  GeMatrix<FullStorage<ElementType> >
  operator()(const Rows& rows, const Cols& cols) const {
    unsigned I = rows.length();
    unsigned J = cols.length();
    GeMatrix<FullStorage<ElementType> > R(I,J);

    for (unsigned i = 1; i <= I; ++i)
      for (unsigned j = 1; j <= J; ++j)
        R(i,j) = (rows(i) == cols(j) ? 2 : 1);

    return R;
  }
};


/** Compute B += M(rows,cols) * A
 * @param[in] transM  Whether M should be treated as transposed or not.
 *                    Values: Trans, NoTrans
 * @param[in]      M  The rhs matrix of the product
 * @param[in]      A  The lhs matrix of the product
 * @param[in,out]  B  The matrix to accumulate the results into
 * @param[in]   rows  A subset of the rows of M involved in the product
 * @param[in]   cols  A subset of the cols of M involved in the product
 *
 * @note flens matrices are (by default) indexed from 1 rather than 0.
 */
template <typename MA, typename MB, typename RI, typename CI>
void
mm(Transpose transM, const MyMatrix& M,
   const GeMatrix<MA>& A, GeMatrix<MB>& B,
   const DenseVector<RI>& rows, const DenseVector<CI>& cols)
{
  (void) M; // unused -- quiet compiler

  const flens::Underscore<typename GeMatrix<MB>::IndexType> _;

  if (transM == flens::NoTrans) {
    unsigned I = rows.length();
    unsigned K = cols.length();
    for (unsigned i = 1; i <= I; ++i) {
      for (unsigned k = 1; k <= K; ++k) {
        auto Mik = (rows(i) == cols(k) ? 2 : 1);
        B(i,_) += Mik * A(k,_);
      }
    }
  } else {
    unsigned I = cols.length();
    unsigned K = rows.length();
    for (unsigned i = 1; i <= I; ++i) {
      for (unsigned k = 1; k <= K; ++k) {
        auto Mki = (cols(i) == rows(k) ? 2 : 1);
        B(i,_) += Mki * A(k,_);
      }
    }
  }
}



int main(int argc, char** argv) {
  unsigned N = 1000;
  unsigned M = 1000;
  bool checkErrors = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-M") == 0) {
      M = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // PLR Parameters for the submatrix blocks
  static constexpr std::size_t max_rank       = 10;
  static constexpr double      eps_tol        = 1e-10;
  static constexpr unsigned    start_at_level = 1;

  // Define an implicit matrix to work with
  MyMatrix Mat;

  // Define the dimension of the sources and targets
  constexpr unsigned SD = 2;
  constexpr unsigned TD = 2;

  // Define the points (sources) corresponding to the columns
  using source_type = Vec<SD,double>;
  std::vector<source_type> sources = fmmtl::random_n(M);
  // Define the points (targets) corresponding to the rows
  using target_type = Vec<TD,double>;
  std::vector<target_type> targets = fmmtl::random_n(N);

  // Call the PLR Compression
  auto plr_plan
      = plr_compression<TD,SD>(Mat, N, M,
                               reinterpret_cast<double*>(targets.data()),
                               //reinterpret_cast<double*>(sources.data()),
                               reinterpret_cast<double*>(targets.data()),
                               max_rank, eps_tol, start_at_level);

  // Construct a set of charges and results
  using charge_type = double;
  using result_type = double;
  std::vector<charge_type> charges = fmmtl::random_n(M);
  std::vector<result_type> results(N);

  // Perform the matvec
  { ScopeClock timer("PLR MatVec: ");

    prod_acc(plr_plan, charges, results);

  } // timer


  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M);
    { ScopeClock timer("Direct: ");
      for (auto& e : exact)
        for (const auto& c : charges)
          e += c;
      for (unsigned i = 0; i < N; ++i)
        exact[i] += charges[i];
    }

    double tot_error_sq = 0;
    double tot_norm_sq = 0;
    double tot_ind_rel_err = 0;
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < results.size(); ++k) {
      //std::cout << results[k] << "\t" << exact[k] << std::endl;

      // Individual relative error
      double rel_error = norm_2(exact[k] - results[k]) / norm_2(exact[k]);
      tot_ind_rel_err += rel_error;
      // Maximum relative error
      max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

      // Total relative error
      tot_error_sq += norm_2_sq(exact[k] - results[k]);
      tot_norm_sq  += norm_2_sq(exact[k]);
    }
    double tot_rel_err = std::sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / results.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }
}
