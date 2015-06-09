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


// Representation of a 1D Hierarchically Off-Diagonal Low-Rank Matrix
template <typename T, typename Tree>
class HODLR_Matrix
    : public GeneralMatrix<HODLR_Matrix<T,Tree> > {
 public:
  typedef T     ElementType;
  typedef int   IndexType;

  // CPU
  using MatrixType  = GeMatrix<FullStorage<T> >;
  using IndexVector = DenseVector<Array<IndexType> >;
  // GPU
  //using MatrixType  = GeMatrix<FullStorage<T,ColMajor,IndexOptions<>,thrust::device_malloc_allocator<T> > >;
  //using IndexVector = DenseVector<Array<IndexType, IndexOptions<>,thrust::device_malloc_allocator<IndexType> > >;


  /** Constructor
   */
  template <typename Matrix, typename ID>
  HODLR_Matrix(const Matrix& A, Tree&& t, ID&& interp_decomp)
      : tree(std::move(t)), Aii(tree), U(tree), V(tree), max_r(0)
  {
    using box_type = typename Tree::box_type;

    //
    // Interpolative Decomposition of off-diagonal blocks
    //

    auto box_end = tree.box_end();
#pragma omp parallel for default(shared)
    for (auto it = tree.box_begin(); it < box_end; ++it) {
      auto box = *it;
      if (box.is_leaf()) {
        // Store the diagonal block
        auto rows = range(box);
        Aii[box] = A(rows,rows);
        continue;
      }

      // Binary trees only for now
      // TODO: Change data structure and logic for all non-diag pairs of boxes
      auto sbox = box.child_begin()[0];
      auto tbox = box.child_begin()[1];

      std::tie(U[tbox],V[sbox]) = interp_decomp(A(range(tbox),range(sbox)));
      ASSERT(num_cols(U[tbox]) == num_cols(V[sbox]));
      max_r = std::max(max_r, num_cols(U[tbox]));

      // Compile-time optimizations
      if (IsSymmetricMatrix<Matrix>::value) {
        U[sbox] = transpose(V[sbox]);
        V[tbox] = transpose(U[tbox]);
      } else if (IsHermitianMatrix<Matrix>::value) {
        U[sbox] = conjTrans(V[sbox]);
        V[tbox] = conjTrans(U[tbox]);
      } else {
        std::tie(U[sbox],V[tbox]) = interp_decomp(A(range(sbox),range(tbox)));
        ASSERT(num_cols(U[sbox]) == num_cols(V[tbox]));
        max_r = std::max(max_r, num_cols(U[sbox]));
      }
    }
  } // end constructor

  IndexType
  numRows() const { return tree.bodies(); }

  IndexType
  numCols() const { return tree.bodies(); }

  double
  compression() const {
    double result = 0;

    for (auto& a : Aii) {
      result += num_rows(a) * num_cols(a);
    }
    for (auto& a : U) {
      result += num_rows(a) * num_cols(a);
    }
    for (auto& a : V) {
      result += num_rows(a) * num_cols(a);
    }
    return result / (numRows()*numCols());
  }


  // private:  TODO: friends

  Tree tree;  // TODO: Copyable/shared_ptr?

  fmmtl::BoxBind<MatrixType,Tree>   Aii;   // Diagonal blocks (leaves only)
  fmmtl::BoxBind<MatrixType,Tree>   U;     // L2T blocks (idx by target box)
  fmmtl::BoxBind<MatrixType,Tree>   V;     // S2M blocks (idx by source box)
  int max_r;                               // Max rank of decompositions

  fmmtl::BoxBind<IndexVector,Tree>  ipiv;  // Pivoting info for Aii LU
};



// y = beta*y + alpha*H*x
template <typename T, typename TR,
          typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans, const ALPHA &alpha,
   const HODLR_Matrix<T,TR> &H, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
  using box_type   = typename TR::box_type;
  using MatrixType = typename HODLR_Matrix<T,TR>::MatrixType;
  using VectorType = typename MatrixType::Vector;
  using IndexType  = typename HODLR_Matrix<T,TR>::IndexType;
  const Underscore<IndexType> _;

  // Don't handle the transposes for now
  (void) trans;
  ASSERT(trans == NoTrans);

  if (y.numRows() == 0) {
    y.resize(H.numRows());
  } else if (beta != BETA(1)) {
    y *= beta;
  }

  // Define the dyadic matvec traversal operator -- TODO: custom traversal?
  auto offdiag = [&](const box_type& sbox, const box_type& tbox) {
    if (sbox == tbox) {           // This is a diagonal block
      if (tbox.is_leaf()) {       // If diagonal leaf, direct matvec
        auto rows = range(tbox);
        y(rows) += alpha * H.Aii[tbox] * x(rows);
        return 0;
      }                           // If diagonal non-leaf, partition both
      return 3;
    }

    // Compute the off-diag block matvec using the low-rank approximation
    // XXX: Revisit temp
    VectorType temp = alpha * transpose(H.V[sbox]) * x(range(sbox));
    y(range(tbox)) += H.U[tbox] * temp;
    return 0;                    // Done with this block
  };

  // Traverse the dyadic tree and perform the matvec
  fmmtl::traverse_if(H.tree.root(), H.tree.root(), offdiag);
}


// B = beta*B + alpha*H*A
template <typename T, typename TR,
          typename ALPHA, typename MA, typename BETA, typename MB>
void
mm(Transpose transH, Transpose transA, const ALPHA &alpha,
   const HODLR_Matrix<T,TR> &H, const GeMatrix<MA> &A,
   const BETA &beta, GeMatrix<MB> &B)
{
  using box_type   = typename TR::box_type;
  using MatrixType = typename HODLR_Matrix<T,TR>::MatrixType;
  using IndexType  = typename HODLR_Matrix<T,TR>::IndexType;
  const Underscore<IndexType> _;

  // Don't handle the transposes for now
  (void) transH; (void) transA;
  ASSERT(transH == NoTrans && transA == NoTrans);

  if (B.numRows() == 0 || B.numCols() == 0) {
    B.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
  } else if (beta != BETA(1)) {
    B *= beta;
  }

  // By level insures independent updates
  //#pragma omp parallel default(shared)
  for (int L = 0; L < H.tree.levels(); ++L) {

    auto box_end = H.tree.box_end(L);
    //#pragma omp for
    for (auto bit = H.tree.box_begin(L); bit < box_end; ++bit) {
      auto box = *bit;

      if (box.is_leaf()) {
        auto rows = range(box);
        B(rows,_) += alpha * H.Aii[box] * A(rows,_);
        continue;
      }

      auto sbox = box.child_begin()[0];
      auto tbox = box.child_begin()[1];
      auto rs = range(sbox);
      auto rt = range(tbox);

      MatrixType Tmp = alpha * transpose(H.V[sbox]) * A(rs,_);
      B(rt,_) += H.U[tbox] * Tmp;

      // Can avoid re-allocation in the common case
      Tmp.reserve(num_cols(H.V[tbox]), num_cols(Tmp));
      Tmp = alpha * transpose(H.V[tbox]) * A(rt,_);
      B(rs,_) += H.U[sbox] * Tmp;
    }
  }
}

//== (gehodlr)trf =======
template <typename T, typename TR, typename VPIV>
typename RestrictTo<IsIntegerDenseVector<VPIV>::value,
typename HODLR_Matrix<T,TR>::IndexType>::Type
trf(HODLR_Matrix<T,TR>& H, VPIV &&)
{
  using MatrixType  = typename HODLR_Matrix<T,TR>::MatrixType;
  using IndexType   = typename HODLR_Matrix<T,TR>::IndexType;
  using IndexVector = typename HODLR_Matrix<T,TR>::IndexVector;
  const Underscore<IndexType> _;

  // Initialize data structures
  H.ipiv    = fmmtl::make_box_binding<IndexVector>(H.tree);

  // For the diagonal boxes of the tree, from leaves to root
#pragma omp parallel default(shared)
  for (int L = H.tree.levels() - 1; L >= 0; --L) {

    auto box_end = H.tree.box_end(L);
#pragma omp for
    for (auto bit = H.tree.box_begin(L); bit < box_end; ++bit) {
      auto box = *bit;
      // Get the index range of this block
      auto rows = range(box);
      auto& AiiInv = H.Aii[box];
      auto& ipiv = H.ipiv[box];

      if (box.is_leaf()) {
        // This block is a diagonal leaf
        // -- Perform a dense solve: A(I,I) * X = B
        // -- Apply A(I,I)^{-1} to all U_L(I,_) up the tree

        // Compute LU factors into Aii
        flens::lapack::trf(AiiInv, ipiv, B(rows,_));

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == H.tree.root()); box = box.parent()) {
          auto& pU = H.U[box];
          auto urows = rows - box.body_begin().index();
          // Solve AX = U and store into U -- Uses previously computed LU
          flens::lapack::trs(NoTrans, AiiInv, ipiv, pU(urows,_));
        }
      } else {
        // This block is a diagonal non-leaf.
        // All descendent diagonal blocks have been processed.
        // This block is now of the form:
        // |     I      U_0 V_1^T |
        // | U_1 V_0^T       I    |
        // Which can be inverted via the Sherman-Morrison-Woodbury formula:
        // (I + U V^T)^{-1} = I - U (I + V^T U)^{-1} V^T
        // where U = [U_0, 0; 0, U_1] and V^T = [0, V_1^T; V_0^T, 0].
        // Thus, we invert an R x R matrix rather than the full N x N block.

        // Get the child boxes and factorizations
        auto cbox0 = *(box.child_begin()+0);
        auto cbox1 = *(box.child_begin()+1);
        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);
        auto& U0 = H.U[cbox0];
        auto& U1 = H.U[cbox1];
        auto& V0 = H.V[cbox0];
        auto& V1 = H.V[cbox1];

        unsigned r = num_cols(V1), R = r + num_cols(V0);
        unsigned c = num_cols(U0), C = c + num_cols(U1);
        ASSERT(R == C);

        // Construct the matrix K = I + [0,V1T;V0T,0] * [U0,0;0,U1]
        AiiInv = MatrixType(R,C);
        AiiInv.diag(0) = 1;
        AiiInv(_(r+1,R), _(  1,c)) = transpose(V0) * U0;   // Upper right block
        AiiInv(_(  1,r), _(c+1,C)) = transpose(V1) * U1;   // Lower left block

        // Compute LU factors into Aii
        flens::lapack::trf(AiiInv, ipiv);

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == H.tree.root()); box = box.parent()) {
          auto& pU = H.U[box];
          auto urows0 = rows0 - box.body_begin().index();
          auto urows1 = rows1 - box.body_begin().index();

          // Apply the SMW to pU: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
          MatrixType Tmp(R,num_cols(pU));
          Tmp(_(  1,r),_) = transpose(V1) * pU(urows1,_);
          Tmp(_(r+1,R),_) = transpose(V0) * pU(urows0,_);
          // Solve KX = T and store into T -- Uses previously computed LU
          flens::lapack::trs(NoTrans, AiiInv, ipiv, Tmp);
          // Apply U to complete the SMW
          pU(urows0,_) -= U0 * Tmp(_(  1,r),_);
          pU(urows1,_) -= U1 * Tmp(_(r+1,R),_);
        }
      }
    }
  }

  return 0;
}

//== (gehodlr)sv =======
template <typename T, typename TR, typename VPIV, typename MB>
typename RestrictTo<IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
typename RemoveRef<MB>::Type::IndexType>::Type
sv(HODLR_Matrix<T,TR>& H, VPIV &&, MB &&B)
{
  using MatrixType  = typename HODLR_Matrix<T,TR>::MatrixType;
  using IndexType   = typename HODLR_Matrix<T,TR>::IndexType;
  using IndexVector = typename HODLR_Matrix<T,TR>::IndexVector;
  const Underscore<IndexType> _;

  // Initialize data structures
  H.ipiv    = fmmtl::make_box_binding<IndexVector>(H.tree);

  // For the diagonal boxes of the tree, from leaves to root
#pragma omp parallel default(shared)
  for (int L = H.tree.levels() - 1; L >= 0; --L) {

    auto box_end = H.tree.box_end(L);
#pragma omp for
    for (auto bit = H.tree.box_begin(L); bit < box_end; ++bit) {
      auto box = *bit;
      // Get the index range of this block
      auto rows = range(box);
      auto& AiiInv = H.Aii[box];
      auto& ipiv = H.ipiv[box];

      if (box.is_leaf()) {
        // This block is a diagonal leaf
        // -- Perform a dense solve: A(I,I) * X = B
        // -- Apply A(I,I)^{-1} to all U_L(I,_) up the tree

        // Solve AX = B and store into X -- Also stores LU factors into Aii
        flens::lapack::sv(AiiInv, ipiv, B(rows,_));

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == H.tree.root()); box = box.parent()) {
          auto& pU = H.U[box];
          auto urows = rows - box.body_begin().index();
          // Solve AX = U and store into U -- Uses previously computed LU
          flens::lapack::trs(NoTrans, AiiInv, ipiv, pU(urows,_));
        }
      } else {
        // This block is a diagonal non-leaf.
        // All descendent diagonal blocks have been processed.
        // This block is now of the form:
        // |     I      U_0 V_1^T |
        // | U_1 V_0^T       I    |
        // Which can be inverted via the Sherman-Morrison-Woodbury formula:
        // (I + U V^T)^{-1} = I - U (I + V^T U)^{-1} V^T
        // where U = [U_0, 0; 0, U_1] and V^T = [0, V_1^T; V_0^T, 0].
        // Thus, we invert an R x R matrix rather than the full N x N block.

        // Get the child boxes and factorizations
        auto cbox0 = *(box.child_begin()+0);
        auto cbox1 = *(box.child_begin()+1);
        auto& U0 = H.U[cbox0];
        auto& U1 = H.U[cbox1];
        auto& V0 = H.V[cbox0];
        auto& V1 = H.V[cbox1];

        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);
        unsigned r = num_cols(V1), R = r + num_cols(V0);
        unsigned c = num_cols(U0), C = c + num_cols(U1);
        ASSERT(R == C);

        // Construct the matrix K = I + [0,V1T;V0T,0] * [U0,0;0,U1]
        AiiInv = MatrixType(R,C);
        AiiInv.diag(0) = 1;
        AiiInv(_(r+1,R), _(  1,c)) = transpose(V0) * U0;   // Upper right block
        AiiInv(_(  1,r), _(c+1,C)) = transpose(V1) * U1;   // Lower left block

        // Apply the SMW to X: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
        MatrixType Tmp(R,num_cols(B));
        Tmp(_(  1,r),_) = transpose(V1) * B(rows1,_);
        Tmp(_(r+1,R),_) = transpose(V0) * B(rows0,_);
        // Solve AiiInv X = Tmp and store into Tmp -- Also stores LU factors
        flens::lapack::sv(AiiInv, ipiv, Tmp);
        // Apply U to complete the SMW
        B(rows0,_) -= U0 * Tmp(_(  1,r),_);
        B(rows1,_) -= U1 * Tmp(_(r+1,R),_);

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == H.tree.root()); box = box.parent()) {
          auto& pU = H.U[box];
          auto urows0 = rows0 - box.body_begin().index();
          auto urows1 = rows1 - box.body_begin().index();

          // Apply the SMW to pU: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
          MatrixType Tmp(R,num_cols(pU));
          Tmp(_(  1,r),_) = transpose(V1) * pU(urows1,_);
          Tmp(_(r+1,R),_) = transpose(V0) * pU(urows0,_);
          // Solve KX = T and store into T -- Uses previously computed LU
          flens::lapack::trs(NoTrans, AiiInv, ipiv, Tmp);
          // Apply U to complete the SMW
          pU(urows0,_) -= U0 * Tmp(_(  1,r),_);
          pU(urows1,_) -= U1 * Tmp(_(r+1,R),_);
        }
      }
    }
  }

  return 0;
}


//== (gehodlr)trs =======
template <typename T, typename TR, typename VPIV, typename MB>
typename RestrictTo<IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
typename RemoveRef<MB>::Type::IndexType>::Type
trs(Transpose trans, HODLR_Matrix<T,TR>& H, VPIV &&, MB &&B)
{
  using MatrixType  = typename HODLR_Matrix<T,TR>::MatrixType;
  using IndexType   = typename HODLR_Matrix<T,TR>::IndexType;
  const Underscore<IndexType> _;

  // Don't mess with transpose right now
  (void) trans;
  ASSERT(trans == NoTrans);

  // Reuse computed inverses

  // For the diagonal boxes of the tree, from leaves to root
  for (int L = H.tree.levels() - 1; L >= 0; --L) {

    for (auto box : boxes(L, H.tree)) {
      // Get the index range of this block
      auto rows = range(box);
      auto& AiiInv = H.Aii[box];
      auto& ipiv = H.ipiv[box];

      if (box.is_leaf()) {
        // Solve AX = B and store into B -- Uses previously computed LU
        flens::lapack::trs(NoTrans, AiiInv, ipiv, B(rows, _));
      } else {
        // Get the child boxes and factorizations
        auto cbox0 = *(box.child_begin()+0);
        auto cbox1 = *(box.child_begin()+1);
        auto& U0 = H.U[cbox0];
        auto& U1 = H.U[cbox1];
        auto& V0 = H.V[cbox0];
        auto& V1 = H.V[cbox1];

        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);
        unsigned r = num_cols(V1), R = r + num_cols(V0);

        // Apply the SMW to X: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
        MatrixType Tmp(R,num_cols(B));
        Tmp(_(  1,r),_) = transpose(V1) * B(rows1,_);
        Tmp(_(r+1,R),_) = transpose(V0) * B(rows0,_);
        // Solve KX = T and store into T -- Use previously computed LU
        flens::lapack::trs(NoTrans, AiiInv, ipiv, Tmp);
        // Apply U to complete the SMW
        B(rows0,_) -= U0 * Tmp(_(  1,r),_);
        B(rows1,_) -= U1 * Tmp(_(r+1,R),_);
      }
    }
  }

  return 0;
}

//== (gehodlr)det =======
template <typename T, typename TR>
std::pair<T, int>
det(HODLR_Matrix<T,TR>& H)
{
  using std::abs;

  // Compute the determinant in the form a*10^b to prevent over/underflow
  // and avoid computing log10() each iteration.
  T a = T(1);
  int b = 0;

  // For each factored diagonal block
  for (auto box : boxes(H.tree)) {
    const auto& ipiv = H.ipiv[box];
    const auto Aii = H.Aii[box].diag(0);

    for (auto i = Aii.firstIndex(); i <= Aii.lastIndex(); ++i) {
      a *= Aii(i);

      if (ipiv(i) != i)
        a = -a;

      // TODO: Consider case a == 0
      while (abs(a) <  1e-15) { a *= 1e+16; b -= 16; }
      while (abs(a) >= 1e+16) { a *= 1e-16; b += 16; }
    }
  }

  // Force 1 < abs(a) < 10
  while (abs(a) <  1e0) { a *= 1e+1; b -= 1; }
  while (abs(a) >= 1e1) { a *= 1e-1; b += 1; }

  return {a,b};
}

namespace blas {
using flens::mv;
using flens::mm;
}
namespace lapack {
using flens::trf;
using flens::trs;
using flens::sv;
using flens::det;
}

} // end namespace flens
