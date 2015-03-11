#pragma once

#include <type_traits>

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

template <typename T, typename Tree, char Class>
class HODLR_Matrix;


template <typename M>
struct MatrixChar;

template <typename S>
struct MatrixChar<GeneralMatrix<S> > {
  static constexpr char value = 'G';
};

template <typename S>
struct MatrixChar<SymmetricMatrix<S> > {
  static constexpr char value = 'S';
};

template <typename S>
struct MatrixChar<HermitianMatrix<S> > {
  static constexpr char value = 'H';
};


template <typename H>
struct HODLRMatrixBase;

template <typename T, typename TR>
struct HODLRMatrixBase<HODLR_Matrix<T,TR,'G'> > {
  using type = GeneralMatrix<HODLR_Matrix<T,TR,'G'> >;
  using matrix = GeMatrix<FullStorage<T> >;
};

template <typename T, typename TR>
struct HODLRMatrixBase<HODLR_Matrix<T,TR,'S'> > {
  using type = SymmetricMatrix<HODLR_Matrix<T,TR,'S'> >;
  using matrix = SyMatrix<FullStorage<T> >;
};

template <typename T, typename TR>
struct HODLRMatrixBase<HODLR_Matrix<T,TR,'H'> > {
  using type = HermitianMatrix<HODLR_Matrix<T,TR,'H'> >;
  using matrix = HeMatrix<FullStorage<T> >;
};


// Representation of a 1D Hierarchically Off-Diagonal Low-Rank Matrix
template <typename T, typename Tree, char Class = 'G'>
class HODLR_Matrix
    : public HODLRMatrixBase<HODLR_Matrix<T,Tree,Class> >::type {
 public:
  typedef T     ElementType;
  typedef int   IndexType;

  using ThisType       = HODLR_Matrix<T,Tree,Class>;
  using DiagMatrixType = typename HODLRMatrixBase<ThisType>::matrix;
  using GeMatrixType   = GeMatrix<FullStorage<T> >;
  using IndexVector    = DenseVector<Array<IndexType> >;

  /** Constructor
   */
  template <typename Matrix, typename ID>
  HODLR_Matrix(const Matrix& A, Tree&& t, ID&& interp_decomp)
      : tree(std::move(t)), Aii(tree), U(tree), VT(tree)
  {
    using box_type = typename Tree::box_type;

    //
    // Interpolative Decomposition of off-diagonal blocks
    //

    // Define the dyadic decomp traversal operator -- TODO: Use custom traversal?
    auto decomp_offdiag = [&](const box_type& sbox, const box_type& tbox) {
      if (sbox == tbox) {  // This is a diagonal block
        if (tbox.is_leaf()) {         // If diagonal leaf, store diagonal block and done
          auto rows = range(tbox);
          Aii[tbox] = A(rows,rows);
          return 0;
        }
        return 3;                     // Else, split both and recurse
      }

      // If this is a symmetric HODLR, skip lower tri blocks
      if ((Class == 'S' || Class == 'H') && (sbox < tbox))
        return 0;

      // Apply the interpolative decomposition to this off-diag block
      std::tie(U[tbox],VT[sbox]) = interp_decomp(A(range(tbox),range(sbox)));
      //std::cout << "Level " << tbox.level() << ": " << num_rows(tU) << "," << num_cols(tU) << " -- " << num_rows(sVT) << "," << num_cols(sVT) << std::endl;
      assert(num_cols(U[tbox]) != 0 && num_rows(VT[sbox]) != 0);
      return 0;                      // Done with this block
    };

    // Traverse the dyadic tree and decompose the off-diagonal blocks
    fmmtl::traverse_if(tree.root(), tree.root(), decomp_offdiag);
  } // end constructor

  IndexType
  numRows() const { return tree.bodies(); }

  IndexType
  numCols() const { return tree.bodies(); }



  // private:  TODO: friends

  Tree tree;  // TODO: Copyable/shared_ptr?

  fmmtl::BoxBind<DiagMatrixType,Tree>  Aii;   // Diagonal blocks (leaves only)
  fmmtl::BoxBind<GeMatrixType,Tree>    U;     // L2T blocks (idx by target box)
  fmmtl::BoxBind<GeMatrixType,Tree>    VT;    // S2M blocks (idx by source box)

  fmmtl::BoxBind<IndexVector,Tree>     ipiv;  // Pivoting info for Aii LU
};


template <typename T, typename Tree>
using HODLR_GeMatrix = HODLR_Matrix<T,Tree,'G'>;

template <typename T, typename Tree>
using HODLR_SyMatrix = HODLR_Matrix<T,Tree,'S'>;

template <typename T, typename Tree>
using HODLR_HeMatrix = HODLR_Matrix<T,Tree,'H'>;




/** HODLR matvec forwards to HODLR matmat */
template <typename T, typename TR, char HC,
          typename ALPHA, typename VX, typename BETA, typename VY>
void
mv(Transpose trans, const ALPHA &alpha,
   const HODLR_Matrix<T,TR,HC> &H, const DenseVector<VX> &x,
   const BETA &beta, DenseVector<VY> &y)
{
  typedef typename DenseVector<VX>::ElementType  XElementType;
  const auto nx = x.length();
  const StorageOrder Xorder = DenseVector<VX>::Engine::order;
  GeMatrix<FullStorageView<XElementType, Xorder> >  X(nx, 1, x, nx);

  typedef typename DenseVector<VY>::ElementType  YElementType;
  const auto ny = y.length();
  const StorageOrder Yorder = DenseVector<VY>::Engine::order;
  GeMatrix<FullStorageView<YElementType, Yorder> >  Y(ny, 1, y, ny);

  return mm(trans, NoTrans, alpha, H, X, beta, Y);
}


// B = beta*B + alpha*H*A
template <typename T, typename TR, char HC,
          typename ALPHA, typename MA, typename BETA, typename MB>
void
mm(Transpose transH, Transpose transA, const ALPHA &alpha,
   const HODLR_Matrix<T,TR,HC> &H, const GeMatrix<MA> &A,
   const BETA &beta, GeMatrix<MB> &B)
{
  using box_type   = typename TR::box_type;
  using MatrixType = typename HODLR_Matrix<T,TR,HC>::GeMatrixType;
  using IndexType  = typename HODLR_Matrix<T,TR,HC>::IndexType;
  const Underscore<IndexType> _;

  // Don't handle the transposes for now
  ASSERT(transH == NoTrans && transA == NoTrans);

  if (B.numRows() == 0 || B.numCols() == 0) {
    B.resize(A.numRows(), A.numCols(), A.firstRow(), A.firstCol());
  } else if (beta != BETA(1)) {
    B *= beta;
  }

  // Define the dyadic matvec traversal operator -- TODO: custom traversal?
  auto offdiag = [&](const box_type& sbox, const box_type& tbox) {
    if (sbox == tbox) {           // This is a diagonal block
      if (tbox.is_leaf()) {       // If diagonal leaf, direct matvec
        auto rows = range(tbox);
        B(rows,_) += alpha * H.Aii[tbox] * A(rows,_);
        return 0;
      }                          // If diagonal non-leaf, partition both
      return 3;
    }

    // Compute the off-diag block matvec using the interpolative decomposition.

    // If symmetric and lower tri, use the tranpose of the upper tri decomp
    if (HC == 'S' && sbox < tbox) {
      MatrixType temp = alpha * transpose(H.U[sbox]) * A(range(tbox),_);
      B(range(sbox),_) += transpose(H.VT[tbox]) * temp;
      return 0;
    }
    // If hermetian and lower tri, use the conjtran of the upper tri decomp
    if (HC == 'H' && sbox < tbox) {
      MatrixType temp = alpha * conjTrans(H.U[sbox]) * A(range(tbox),_);
      B(range(sbox),_) += conjTrans(H.VT[tbox]) * temp;
      return 0;
    }
    // Else, use the decomposition of this block
    MatrixType temp = alpha * H.VT[sbox] * A(range(sbox),_);
    B(range(tbox),_) += H.U[tbox] * temp;
    return 0;                    // Done with this block
  };

  // Traverse the dyadic tree and perform the matvec
  fmmtl::traverse_if(H.tree.root(), H.tree.root(), offdiag);
}

//== (gehodlr)sv =======
template <typename T, typename TR, typename VPIV, typename MB>
typename RestrictTo<IsIntegerDenseVector<VPIV>::value
                 && IsGeMatrix<MB>::value,
typename RemoveRef<MB>::Type::IndexType>::Type
sv(HODLR_GeMatrix<T,TR>& H, VPIV &&, MB &&B)
{
  using DiagMatrixType = typename HODLR_GeMatrix<T,TR>::DiagMatrixType;
  using MatrixType  = typename HODLR_GeMatrix<T,TR>::GeMatrixType;
  using IndexType   = typename HODLR_GeMatrix<T,TR>::IndexType;
  using IndexVector = typename HODLR_GeMatrix<T,TR>::IndexVector;
  const Underscore<IndexType> _;

  // Initialize data structures
  H.ipiv    = fmmtl::make_box_binding<IndexVector>(H.tree);

  // For the diagonal boxes of the tree, from leaves to root
  for (int L = H.tree.levels() - 1; L >= 0; --L) {

    for (auto box : boxes(L, H.tree)) {
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
        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);
        auto& U0  = H.U[cbox0];      // box_u_inv[target child 0]
        auto& U1  = H.U[cbox1];      // box_u_inv[target child 1]
        auto& V0T = H.VT[cbox0];     // box_v_inv[source child 0]
        auto& V1T = H.VT[cbox1];     // box_v_inv[source child 1]

        unsigned r = num_rows(V1T), R = r + num_rows(V0T);
        unsigned c = num_cols(U0),  C = c + num_cols(U1);

        // Construct the matrix K = I + [0,V1T;V0T,0] * [U0,0;0,U1]
        AiiInv = MatrixType(R,C);
        AiiInv.diag(0) = 1;
        AiiInv(_(r+1,R), _(  1,c)) = V0T * U0;   // Upper right block
        AiiInv(_(  1,r), _(c+1,C)) = V1T * U1;   // Lower left block

        // Apply the SMW to X: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
        MatrixType Tmp(R,num_cols(B));
        Tmp(_(  1,r),_) = V1T * B(rows1,_);
        Tmp(_(r+1,R),_) = V0T * B(rows0,_);
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
          Tmp(_(  1,r),_) = V1T * pU(urows1,_);
          Tmp(_(r+1,R),_) = V0T * pU(urows0,_);
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
trs(Transpose trans, HODLR_GeMatrix<T,TR>& H, VPIV &&, MB &&B)
{
  using MatrixType  = typename HODLR_GeMatrix<T,TR>::GeMatrixType;
  using IndexType   = typename HODLR_GeMatrix<T,TR>::IndexType;
  const Underscore<IndexType> _;

  // Don't mess with transpose right now
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
        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);
        auto& U0  = H.U[cbox0];      // box_u_inv[target child 0]
        auto& U1  = H.U[cbox1];      // box_u_inv[target child 1]
        auto& V0T = H.VT[cbox0];     // box_v_inv[source child 0]
        auto& V1T = H.VT[cbox1];     // box_v_inv[source child 1]

        unsigned r = num_rows(V1T), R = r + num_rows(V0T);

        // Apply the SMW to X: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
        MatrixType Tmp(R,num_cols(B));
        Tmp(_(  1,r),_) = V1T * B(rows1,_);
        Tmp(_(r+1,R),_) = V0T * B(rows0,_);
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




namespace blas {
using flens::mv;
using flens::mm;
}
namespace lapack {
using flens::sv;
using flens::trs;
}

} // end namespace flens
