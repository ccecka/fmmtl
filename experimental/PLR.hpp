#pragma once

#include "PLR_Matrix.hpp"

#include "util/Probe.hpp"
#include "util/ACA.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/KDTree.hpp"

#include "fmmtl/numeric/Vec.hpp"

struct ProbeSVDOrIdentity {
  int max_rank;
  double eps_tol;

  template <typename Matrix>
  auto operator()(const Matrix& A)
      -> decltype(probe_svd(A, max_rank, eps_tol)) {
    if (std::min(num_rows(A), num_cols(A)) <= max_rank)
      return {{}, A};

    return probe_svd(A, max_rank, eps_tol);
  }
};


/** Construct an PLR Matrix from:
 * 1) a FLENS Matrix,
 * 2) partitioned according to the provided trees,
 * 3) approximating the off-diagonal blocks with the provided Interpolative Decomposition.
 *
 *
 */
template <typename MA, typename STree, typename TTree, typename ID>
flens::PLR_Matrix<typename MA::Impl::ElementType, STree, TTree>
plr(const flens::Matrix<MA>& A, STree&& stree, TTree&& ttree, ID id) {
  return {A.impl(), std::move(stree), std::move(ttree), id};
}

template <typename MA, typename STree, typename TTree>
flens::PLR_Matrix<typename MA::Impl::ElementType, STree, TTree>
plr(const flens::Matrix<MA>& A, STree&& stree, TTree&& ttree) {
  return {A.impl(), std::move(stree), std::move(ttree), ProbeSVDOrIdentity{20,1e-10}};
}

template <typename MA>
flens::PLR_Matrix<typename MA::Impl::ElementType,
                  fmmtl::NDTree<1>, fmmtl::NDTree<1> >
plr(const flens::Matrix<MA>& A) {
  // Generate an integer tree -- TODO: Make implicit integer tree
  using Tree = fmmtl::NDTree<1>;
  std::vector<Vec<1,double> > rows(A.impl().numRows());
  for (unsigned i = 0; i < rows.size(); ++i) rows[i] = Vec<1,double>(i);
  Tree ttree(rows, 32);

  std::vector<Vec<1,double> > cols(A.impl().numCols());
  for (unsigned i = 0; i < cols.size(); ++i) cols[i] = Vec<1,double>(i);
  Tree stree(cols, 32);
  // NDTree is stable, no need to permute data
  //assert(std::equal(ints.begin(), ints.end(), permute_begin(tree, ints.begin())));

  return plr(A.impl(), std::move(stree), std::move(ttree));
}


// Simple API

/** Construct an PLR Matrix from raw data representing a full matrix
 */
template <typename T>
flens::PLR_Matrix<T, fmmtl::NDTree<1>, fmmtl::NDTree<1> >
geplr(char order, T* data, int n, int m, int lda)
{
  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, m, data, lda);
    return plr(A);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, m, data, lda);
    return plr(A);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}

/** Construct an PLR Matrix from raw data representing a symmetric matrix
 */
template <typename T>
flens::PLR_Matrix<T, fmmtl::NDTree<1>, fmmtl::NDTree<1> >
syplr(char order, char uplo, T* data, int n, int lda)
{
  assert(uplo == 'U' || uplo == 'L');
  flens::StorageUpLo _uplo = (uplo == 'U' ? flens::Upper : flens::Lower);

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return plr(A);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return plr(A);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}

/** Construct an PLR Matrix from raw data representing a symmetric matrix
 */
template <typename T>
flens::PLR_Matrix<T, fmmtl::NDTree<1>, fmmtl::NDTree<1> >
heplr(char order, char uplo, T* data, int n, int lda)
{
  assert(uplo == 'U' || uplo == 'L');
  flens::StorageUpLo _uplo = (uplo == 'U' ? flens::Upper : flens::Lower);

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::HeMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return plr(A);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::HeMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return plr(A);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}


//
// Spacially enriched PLR
//

template <unsigned DT, unsigned DS,
          typename T, typename TT, typename TS>
flens::PLR_Matrix<T, fmmtl::NDTree<DT>, fmmtl::NDTree<DS> >
geplr(char order, T* data, int n, int m, int lda,
      const TT* trgs, const TS* srcs,
      unsigned max_rank, double eps_tol, unsigned init_depth = 0)
{
  ScopeClock plr_construction_timer("PLR Matrix Construction: ");

  const flens::Underscore<int> _;

  const Vec<DT,TT>* targets = reinterpret_cast<const Vec<DT,TT>*>(trgs);
  const Vec<DS,TS>* sources = reinterpret_cast<const Vec<DS,TS>*>(srcs);

  fmmtl::NDTree<DT> ttree(targets, targets+n, max_rank);
  fmmtl::NDTree<DS> stree(sources, sources+m, max_rank);

  // Permute matrix to tree order
  std::vector<int> ipiv(std::max(n,m));

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    flens::GeMatrix<Storage> A = Storage(n, m, data, lda);

    // IN PLACE version -- ugly, factor out
    auto tperm = ttree.permute_begin();
    for (int i = 0; i < n; ++i) {
      int k = tperm[i];
      while (k < i) k = ipiv[k];
      ipiv[i] = k;
      if (k != i)
        // Swap row i and k
        flens::blas::swap(A(i+A.firstRow(),_), A(k+A.firstRow(),_));
    }
    auto sperm = stree.permute_begin();
    for (int i = 0; i < m; ++i) {
      int k = sperm[i];
      while (k < i) k = ipiv[k];
      ipiv[i] = k;
      if (k != i)
        // Swap col i and k
        flens::blas::swap(A(_,i+A.firstCol()), A(_,k+A.firstCol()));
    }

    return plr(A, std::move(stree), std::move(ttree),
               ProbeSVDOrIdentity{max_rank,eps_tol});
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    flens::GeMatrix<Storage> A = Storage(n, m, data, lda);

    // IN PLACE version -- ugly, factor out
    auto tperm = ttree.permute_begin();
    for (int i = 0; i < n; ++i) {
      int k = tperm[i];
      while (k < i) k = ipiv[k];
      ipiv[i] = k;
      if (k != i)
        // Swap row i and k
        flens::blas::swap(A(i+A.firstRow(),_), A(k+A.firstRow(),_));
    }
    auto sperm = stree.permute_begin();
    for (int i = 0; i < m; ++i) {
      int k = sperm[i];
      while (k < i) k = ipiv[k];
      ipiv[i] = k;
      if (k != i)
        // Swap col i and k
        flens::blas::swap(A(_,i+A.firstCol()), A(_,k+A.firstCol()));
    }

    return plr(A, std::move(stree), std::move(ttree),
               ProbeSVDOrIdentity{max_rank,eps_tol});
  }
}


//
// Implicit kernel matrix PLR
//

#if 0
template <class Kernel, class T>
flens::PLR_Matrix
plr_matrix(const Kernel& k, const T* point, unsigned n) {

}
#endif
