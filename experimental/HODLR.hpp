#pragma once

#include "HODLR_Matrix.hpp"

#include "util/Probe.hpp"
#include "util/ACA.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/KDTree.hpp"

#include "fmmtl/numeric/Vec.hpp"

struct ProbeSVD {
  int max_rank;
  double eps_tol;

  template <typename Matrix>
  auto operator()(const Matrix& A) -> decltype(probe_svd(A, max_rank, eps_tol)) {
    return probe_svd(A, max_rank, eps_tol);
  }
};

struct ACA {
  int max_rank;
  double eps_tol;

  template <typename Matrix>
  auto operator()(const Matrix& A) -> decltype(adaptive_cross_approx(A, eps_tol, max_rank)) {
      return adaptive_cross_approx(A, eps_tol, max_rank);
    }
};



/** Construct an HODLR Matrix from:
 * 1) a FLENS Matrix,
 * 2) partitioned according to the provided tree,
 * 3) approximating the off-diagonal blocks with the provided Interpolative Decomposition.
 *
 *
 */
template <typename MA, typename Tree, typename ID>
flens::HODLR_Matrix<typename MA::Impl::ElementType, Tree>
hodlr(const flens::Matrix<MA>& A, Tree&& tree, ID id) {
  return {A.impl(), std::move(tree), id};
}

template <typename MA, typename Tree>
flens::HODLR_Matrix<typename MA::Impl::ElementType, Tree>
hodlr(const flens::Matrix<MA>& A, Tree&& tree) {
  return hodlr(A, std::move(tree), ACA{20,1e-10});
}

template <typename MA>
flens::HODLR_Matrix<typename MA::Impl::ElementType, fmmtl::NDTree<1> >
hodlr(const flens::Matrix<MA>& A, int leaf_size) {
  // Generate an integer tree -- TODO: Make implicit integer tree
  using Tree = fmmtl::NDTree<1>;
  std::vector<Vec<1,double> > ints(A.impl().numRows());
  for (int i = 0; i < ints.size(); ++i) ints[i] = Vec<1,double>(i);
  Tree tree(ints, leaf_size);
  // NDTree is stable, no need to permute data
  //assert(std::equal(ints.begin(), ints.end(), permute_begin(tree, ints.begin())));

  return hodlr(A, std::move(tree));
}

template <typename MA>
flens::HODLR_Matrix<typename MA::Impl::ElementType, fmmtl::NDTree<1> >
hodlr(const flens::Matrix<MA>& A) {
  return hodlr(A, 256);
}


// Simple API

/** Construct an HODLR Matrix from raw data representing a full matrix
 */
template <typename T>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
gehodlr(char order, T* data, int n, int lda, int leaf_size)
{
  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    return hodlr(A, leaf_size);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    return hodlr(A, leaf_size);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}

/** Construct an HODLR Matrix from raw data representing a symmetric matrix
 */
template <typename T>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
syhodlr(char order, char uplo, T* data, int n, int lda, int leaf_size)
{
  assert(uplo == 'U' || uplo == 'L');
  flens::StorageUpLo _uplo = (uplo == 'U' ? flens::Upper : flens::Lower);

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return hodlr(A, leaf_size);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return hodlr(A, leaf_size);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}

/** Construct an HODLR Matrix from raw data representing a symmetric matrix
 */
template <typename T>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
hehodlr(char order, char uplo, T* data, int n, int lda, int leaf_size)
{
  assert(uplo == 'U' || uplo == 'L');
  flens::StorageUpLo _uplo = (uplo == 'U' ? flens::Upper : flens::Lower);

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::HeMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return hodlr(A, leaf_size);
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::HeMatrix<Storage> A(Storage(n, n, data, lda), _uplo);
    return hodlr(A, leaf_size);
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}


//
// Spacially enriched HODLR
//


template <typename T, typename P>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
gehodlr(char order, T* data, int n, int lda, const P* point, int leaf_size) {
  // Generate a tree from the spacial hint -- TODO: Arbitrary dimension
  using Tree = fmmtl::NDTree<1>;
  auto p_begin = reinterpret_cast<const Vec<1,P>*>(point);
  Tree tree(p_begin, p_begin+n, leaf_size);

  /*
  // Construct the ipiv array to permute the matrix -- TODO: expose as output
  flens::DenseVector<flens::Array<int> > ipiv(tree.bodies());
  // Permute array to ipiv array
  int i = 1;
  for (auto it = tree.permute_begin(); it != tree.permute_end(); ++it, ++i) {
    int k = *it + 1;
    while (k < i)
      k = ipiv(k);
    ipiv(i) = k;
  }
  // Permute -- XXX: This is N^2, need lazy?
  */


  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    flens::GeMatrix<Storage> A = Storage(n, n, data, lda);

    // Permute

    return hodlr(A, std::move(tree));
  } else if (order == 'r' || order == 'R') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    flens::GeMatrix<Storage> A = Storage(n, n, data, lda);

    // Permute

    return hodlr(A, std::move(tree));
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}


//
// Implicit kernel matrix HODLR
//

#if 0
template <class Kernel, class T>
flens::HODLR_Matrix
hodlr_matrix(const Kernel& k, const T* point, unsigned n) {

}
#endif
