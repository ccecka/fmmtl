#pragma once

#include <type_traits>
#include <cmath>
#include <limits>

//#include <thrust/device_malloc_allocator.h>

#include "fmmtl/numeric/flens.hpp"

#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/traversal/DualTraversal.hpp"

namespace fmmtl {

// Helper functions for transforming tree boxes to flens index ranges
template <typename Box>
flens::Range<int> range(const Box& b) {
  return flens::Range<int>(b.body_begin().index()+1, b.body_end().index());
}

} // end namespace fmmtl



namespace flens {

// Representation of an ND Hierarchically Partitioned Low-Rank Matrix
template <typename T, typename STree, typename TTree>
class PLR_Matrix
    : public GeneralMatrix<PLR_Matrix<T,STree,TTree> > {
 public:
  typedef T     ElementType;
  typedef int   IndexType;

  using TargetTreeType = TTree;
  using SourceTreeType = STree;

  using TargetBoxType = typename TargetTreeType::box_type;
  using SourceBoxType = typename SourceTreeType::box_type;

  // CPU
  using MatrixType  = GeMatrix<FullStorage<T> >;
  using VectorType  = DenseVector<Array<T> >;
  using IndexVector = DenseVector<Array<IndexType> >;
  // GPU
  //using MatrixType  = GeMatrix<FullStorage<T, ColMajor, IndexOptions<>, thrust::device_malloc_allocator<T> > >;
  //using VectorType  = DenseVector<Array<T, IndexOptions<>, thrust::device_malloc_allocator<T> > >;
  //using IndexVector = DenseVector<Array<IndexType, IndexOptions<>, thrust::device_malloc_allocator<IndexType> > >;

  // The trees defining the partitioning of the matrix
  TargetTreeType ttree;
  SourceTreeType stree;

  struct DyadicTreeLeaf {
    SourceBoxType s;      //< Or charge_range for simplicity
    MatrixType U, V;      //< Low-rank approximation of this block
    template <class MatrixU, class MatrixV>
    DyadicTreeLeaf(const SourceBoxType& _s, MatrixU&& _u, MatrixV&& _v)
        : s(_s), U(std::forward<MatrixU>(_u)), V(std::forward<MatrixV>(_v)) {}
  };

  // List of dyadic blocks by target box index
  std::vector<std::vector<DyadicTreeLeaf>> leaf_nodes;
  // Count of blocks by level
  std::vector<char> leaf_count;

  /** Constructor
   */
  template <typename Matrix, typename ID>
  PLR_Matrix(const Matrix& A, TTree&& tt, STree&& st, ID&& interp_decomp,
             int init_depth = 0)
      : ttree(std::move(tt)), stree(std::move(st)),
        leaf_nodes(ttree.size()), leaf_count(ttree.size())
  {
    // Construct an evaluator for the traversal algorithm
    // See fmmtl::traverse_if documentation
    auto evaluator = [&] (const SourceBoxType& s, const TargetBoxType& t) {
      std::cout << "TargetBox: " << t << std::endl;
      std::cout << "SourceBox: " << s << std::endl;

      // Attempt to factor the block of the matrix
      MatrixType U, V;
      std::tie(U, V) = interp_decomp(A(range(t),range(s)));

      // Failure is signalled by a 0-sized V -- Other specification?
      if (num_rows(V) == 0) {
        int flag = ((!s.is_leaf()) << 1) | ((!t.is_leaf()) << 0);

        if (flag == 0) {
          std::cout << "LEAF ACCEPTED, " << range(t).length() << "-by-" << range(s).length() << std::endl;
          // No recursion, so save the matrix
          add_leaf(s, t, MatrixType(), A(range(t),range(s)));
        }

        // Recurse by splitting the source and/or target boxes
        return flag;
      } else {
        std::cout << "ACCEPTED BLOCK, Rank " << num_cols(U) << " from " << range(t).length() << "-by-" << range(s).length() << std::endl;
        // Accept the block and store low-rank approximation
        add_leaf(s, t, std::move(U), std::move(V));
        // Do not recurse further
        return 0;
      }
    };

    // Perform the traversal, starting at init_depth
    for (auto tbox : boxes(init_depth, ttree)) {
      for (auto sbox : boxes(init_depth, stree)) {
        fmmtl::traverse_if(sbox, tbox, evaluator);
      }
    }
  } // end constructor

  IndexType
  numRows() const { return ttree.bodies(); }

  IndexType
  numCols() const { return stree.bodies(); }

 private:

  template <class MatrixU, class MatrixV>
  void add_leaf(const SourceBoxType& s, const TargetBoxType& t,
                MatrixU&& u, MatrixV&& v) {
    leaf_nodes[t.index()].emplace_back(s,
                                       std::forward<MatrixU>(u),
                                       std::forward<MatrixV>(v));
    leaf_count[t.level()] = 1;
  }

};


// y = beta*y + alpha*H*x
template <typename T, typename ST, typename TT,
          typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans, const ALPHA &alpha,
   const PLR_Matrix<T,ST,TT> &H, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
  // Don't handle transpose for now
  (void) trans;
  ASSERT(trans==NoTrans);

  if (y.numRows() == 0) {
    ASSERT(beta == BETA(0));
    y.resize(H.numRows());
  } else if (beta != BETA(1)) {
    y *= beta;
  }

  using VectorType = typename PLR_Matrix<T,ST,TT>::VectorType;

  // Instead of transposing y in and out, create a temporary to accumulate into
  VectorType py(y.length());
  VectorType px(x.length());

  // In-permutation
  auto sperm = H.stree.permute_begin();
  for (int i = 0; i < H.numCols(); ++i)
    px(i+px.firstIndex()) = x(sperm[i]+x.firstIndex());


#pragma omp parallel default(shared)
  for (unsigned level = 0; level < H.ttree.levels(); ++level) {
    if (H.leaf_count[level] == 0)  // Nothing to be done on this level
      continue;

    // OMP requires iterator form of this loop
    auto tb_end = H.ttree.box_end(level);
    // In parallel over all of the target boxes of this level
#pragma omp for
    for (auto tb_i = H.ttree.box_begin(level); tb_i < tb_end; ++tb_i) {
      // The target box
      auto tbox = *tb_i;
      auto cols = range(tbox);

      // For all the source boxes in this tbox-level list
      for (auto& leaf : H.leaf_nodes[tbox.index()]) {
        auto rows = range(leaf.s);

        // Apply the (t,s) block
        if (num_rows(leaf.U) == 0) {
          py(cols) += alpha * leaf.V * px(rows);
        } else {
          // XXX revisit temp
          VectorType temp = alpha * leaf.V * px(rows);
          py(cols) += leaf.U * temp;
        }
      }
    }
  }

  // Out-permutation
  auto tperm = H.ttree.permute_begin();
  for (int i = 0; i < H.numRows(); ++i)
    y(tperm[i]+y.firstIndex()) += py(i+py.firstIndex());
}


// B = beta*B + alpha*H*A
template <typename T, typename ST, typename TT,
          typename ALPHA, typename MA, typename BETA, typename MB>
void
mm(Transpose transH, Transpose transA, const ALPHA &alpha,
   const PLR_Matrix<T,ST,TT> &H, const GeMatrix<MA> &A,
   const BETA &beta, GeMatrix<MB> &B)
{
  using MatrixType = typename PLR_Matrix<T,ST,TT>::MatrixType;
  using IndexType  = typename PLR_Matrix<T,ST,TT>::IndexType;
  const Underscore<IndexType> _;

  // Don't handle the transposes for now
  (void) transH; (void) transA;
  ASSERT(transH == NoTrans && transA == NoTrans);

  if (B.numRows() == 0 || B.numCols() == 0) {
    B.resize(A.numRows(), A.numCols());
  } else if (beta != BETA(1)) {
    B *= beta;
  }

  // Instead of transposing y in and out, create a temporary to accumulate into
  MatrixType pB(B.numRows(), B.numCols());
  MatrixType pA(A.numRows(), A.numCols());

  // In-permutation
  auto sperm = H.stree.permute_begin();
  for (int i = 0; i < H.numCols(); ++i)
    pA(i+pA.firstRow(),_) = A(sperm[i]+A.firstRow(),_);

#pragma omp parallel default(shared)
  for (unsigned level = 0; level < H.ttree.levels(); ++level) {
    if (H.leaf_count[level] == 0)  // Nothing to be done on this level
      continue;

    // OMP requires iterator form of this loop
    auto tb_end = H.ttree.box_end(level);
    // In parallel over all of the target boxes of this level
#pragma omp for
    for (auto tb_i = H.ttree.box_begin(level); tb_i < tb_end; ++tb_i) {
      // The target box
      auto tbox = *tb_i;
      auto cols = range(tbox);

      // For all the source boxes in this tbox-level list
      for (auto& leaf : H.leaf_nodes[tbox.index()]) {
        auto rows = range(leaf.s);

        // Apply the (t,s) block
        if (num_rows(leaf.U) == 0) {
          pB(cols,_) += alpha * leaf.V * pA(rows,_);
        } else {
          // XXX revisit temp
          MatrixType temp = alpha * leaf.V * pA(rows,_);
          pB(cols,_) += leaf.U * temp;
        }
      }
    }
  }

  // Out-permutation
  auto tperm = H.ttree.permute_begin();
  for (int i = 0; i < H.numRows(); ++i)
    B(tperm[i]+B.firstRow(),_) += pB(i+pB.firstRow(),_);
}

namespace blas {
using flens::mv;
using flens::mm;
}
namespace lapack {
}

} // end namespace flens
