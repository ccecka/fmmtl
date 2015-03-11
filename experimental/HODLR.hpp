#pragma once

#include "HODLR_Matrix.hpp"

#include "util/Probe.hpp"
#include "fmmtl/tree/NDTree.hpp"

#include "fmmtl/numeric/Vec.hpp"

struct ProbeSVD {
  int max_rank;
  double eps_tol;

  template <typename Matrix>
  auto operator()(const Matrix& A) -> decltype(probe_svd(A, max_rank, eps_tol)) {
    return probe_svd(A, max_rank, eps_tol);
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
  return hodlr(A, std::move(tree), ProbeSVD{20,1e-10});
}


// Simple API

/** Construct an HODLR Matrix from raw data representing a full matrix
 */
template <typename T>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
gehodlr(char order, T* data, int n, int lda, int leaf_size) {
  // Generate an integer tree -- TODO: Make implicit tree
  using Tree = fmmtl::NDTree<1>;
  std::vector<Vec<1,double> > ints(n);
  for (int i = 0; i < n; ++i) ints[i] = Vec<1,double>(i);
  Tree tree(ints, leaf_size);
  // NDTree is stable, no need to permute data
  //assert(std::equal(ints.begin(), ints.end(), tree.body_permute(ints.begin())))

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    return hodlr(A, std::move(tree));
  } else {
    assert(order == 'r' || order == 'R');
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    return hodlr(A, std::move(tree));
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}

template <typename T, typename P>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
gehodlr(char order, T* data, int n, int lda, const P* point, int leaf_size) {
  // Generate a tree from the spacial hint
  using Tree = fmmtl::NDTree<1>;
  auto p_begin = reinterpret_cast<const Vec<1,P>*>(point);
  Tree tree(p_begin, p_begin+n, leaf_size);

  // PERMUTE

  if (order == 'n' || order == 'N') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    // Construct the interpolative decomposition -- TODO: Expose
    return hodlr(A, std::move(tree));
  } else if (order == 't' || order == 'T' || order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::GeMatrix<Storage> A = Storage(n, n, data, lda);
    // Construct the interpolative decomposition -- TODO: Expose
    return hodlr(A, std::move(tree));
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}


#if 0
/** Construct an HODLR Matrix from raw data representing a symmetric matrix
 */
template <typename T>
flens::HODLR_Matrix<T, fmmtl::NDTree<1> >
syhodlr(char order, char uplo, T* data, int n, int lda, int leaf_size) {
  // Generate an integer tree -- TODO: Make implicit tree
  using Tree = fmmtl::NDTree<1>;
  std::vector<Vec<1,double> > ints(n);
  for (int i = 0; i < n; ++i) ints[i] = Vec<1,double>(i);
  Tree tree(ints, leaf_size);
  // NDTree is stable, no need to permute data
  //assert(std::equal(ints.begin(), ints.end(), tree.body_permute(ints.begin())))

  assert(uplo == 'U' || uplo == 'L');

  if (order == 'c' || order == 'C') {
    using Storage = flens::FullStorageView<T, flens::ColMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), uplo);
    return hodlr(A, std::move(tree));
  } else {
    assert(order == 'r' || order == 'R');
    using Storage = flens::FullStorageView<T, flens::RowMajor>;
    const flens::SyMatrix<Storage> A(Storage(n, n, data, lda), uplo);
    return hodlr(A, std::move(tree));
  }
  fprintf(stderr, "Assertion failed: %s, %s(), %d at \'%s\'\n",
          __FILE__, __func__, __LINE__, "Invalid parameter: order");
  abort();
}
#endif






// Implicit Matrix API
#if 0
template <class Kernel, class T>
flens::HODLR_Matrix
hodlr_matrix(const Kernel& k, const T* point, unsigned n) {

}
#endif
