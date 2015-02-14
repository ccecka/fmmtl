#include <memory>
#include <iterator>

#include "fmmtl/numeric/flens.hpp"

#include "fmmtl/numeric/Vec.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/traversal/DualTraversal.hpp"
#include "util/Probe.hpp"


int main() {
  //
  // Parameters
  //

  unsigned N = 1 << 14;     // rows
  unsigned M = N;           // cols -- square for inversion
  unsigned leaf_size = 64;  // maximum size of the tree leaves
  // TODO: Make interpolative decomposition closure
  double tolerance = 1e-10; // tolerance of the interpolative decomposition
  unsigned max_rank = 30;   // maximum rank of the interpolative decomposition

  // Define types from FLENS
  using namespace flens;
  //using ZeroBased      = IndexOptions<int, 0>;
  using MatrixType     = GeMatrix<FullStorage<double, ColMajor> >;
  using VectorType     = DenseVector<Array<double> >;
  using IndexType      = typename MatrixType::IndexType;
  using IndexVector    = DenseVector<Array<IndexType> >;
  const Underscore<IndexType> _;

  //
  // Setup
  //

  // Initialize the sources/targets as random values
  using source_type = Vec<1,double>;
  using target_type = Vec<1,double>;
  std::vector<source_type> sources = fmmtl::random_n(M);
  std::vector<target_type>& targets = sources;

  // Construct the trees
  // TODO: Implicit trees if over the integers, though this could be any 1D Tree
  using source_tree_type = fmmtl::NDTree<1>;
  using target_tree_type = fmmtl::NDTree<1>;
  using target_box_type  = typename target_tree_type::box_type;
  using source_box_type  = typename source_tree_type::box_type;
  source_tree_type  source_tree(sources, leaf_size);
  target_tree_type& target_tree = source_tree; // Only one -- force square

  // TODO: Null-op if implicit trees
  // Permute the sources/targets to match the body order in the tree
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_targets = make_body_binding(target_tree, targets);

  // Create a test matrix
  MatrixType A(N,M);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= M; ++j)
      //A(i,j) = 1;
      //A(i,j) = std::exp(-norm_2_sq(p_targets[i-1] - p_sources[j-1]));
      A(i,j) = std::exp(-2*norm_2_sq(std::sin(2*3.141592653*(p_targets[i-1][0]-p_sources[i-1][0])))/2);
  A.diag(0) = 2;

  // Initialize a random RHS,   A*X = B
  unsigned Nb = 1;
  MatrixType B(N,Nb);
  for (unsigned i = 1; i <= N; ++i)
    for (unsigned j = 1; j <= Nb; ++j)
      B(i,j) = fmmtl::random<double>::get();

  // Helper functions for transforming tree boxes to flens index ranges
  auto range = [](const source_box_type& b) {
    using flens::Range;
    return Range<IndexType>(b.body_begin().index()+1, b.body_end().index());
  };
  auto shift = [](const flens::Range<IndexType>& r, const source_box_type& b) {
    IndexType i = b.body_begin().index();
    return flens::Range<IndexType>(r.firstIndex()-i, r.lastIndex()-i);
  };

  //
  // Interpolative Decomposition of off-diagonal blocks
  //

  auto box_u = make_box_binding<MatrixType>(target_tree);
  auto box_v = make_box_binding<MatrixType>(source_tree);

  // Define the dyadic decomp traversal operator -- TODO: Use custom traversal?
  auto decomp_offdiag = [&](const source_box_type& sbox,
                            const target_box_type& tbox) {
    if (sbox == tbox) {  // This is a diagonal block
      if (sbox.is_leaf())
        return 0;        // If diagonal leaf, don't do anything
      else
        return 3;        // If diagonal non-leaf, partition both boxes
    }

    // Apply the interpolative decomposition to this off-diag block
    MatrixType& U  = box_u[tbox];
    MatrixType& VT = box_v[sbox];
    std::tie(U,VT) = probe_svd(A(range(tbox),range(sbox)), max_rank, tolerance);
    //std::cout << "Level " << tbox.level() << ": " << num_rows(U) << "," << num_cols(U) << " -- " << num_rows(VT) << "," << num_cols(VT) << std::endl;
    return 0; // Done with this block
  };

  std::cout << "Computing interpolative decompositions..." << std::endl;
  { ScopeClock timer;
  // Traverse the dyadic tree and decompose the off-diagonal blocks
  fmmtl::traverse_if(source_tree.root(), target_tree.root(), decomp_offdiag);
  } // end scopeclock

#if 1
  //
  // HODLR matvec using the interpolative decompositions
  //

  MatrixType testR(N,Nb), exactR(N,Nb);

  // Define the dyadic matvec traversal operator -- TODO: Use custom traversal?
  auto matvec_offdiag = [&](const source_box_type& sbox,
                            const target_box_type& tbox) {
    if (sbox == tbox) {  // This is a diagonal block
      if (sbox.is_leaf()) {       // If diagonal leaf, direct matvec
        auto rows = range(tbox);
        testR(rows,_) += A(rows,rows) * B(rows,_);
        return 0;
      } else {                    // If diagonal non-leaf, partition both boxes
        return 3;
      }
    }

    // Compute the off-diag block matvec using the interpolative decomposition.
    auto rows = range(tbox);
    auto cols = range(sbox);
    testR(rows,_) += box_u[tbox] * MatrixType(box_v[sbox] * B(cols,_));
    return 0; // Done with this block
  };

  // Traverse the dyadic tree and perform the matvec
  { ScopeClock timer("HODLR  Matvec: ");
  fmmtl::traverse_if(source_tree.root(), target_tree.root(), matvec_offdiag);
  }

  { ScopeClock timer(std::string("Direct Matvec: "));
  exactR = A * B;
  }
  MatrixType ResR = exactR - testR;
  std::cout << "Matvec residual norm_F = " << frobenius_norm(ResR) << std::endl;
#endif

  //
  // HODLR inverse using interpolative decompositions
  //

  // Copy the low-rank block approximations to the inverted factor approximations
  auto  box_u_inv = box_u;
  auto& box_v_inv = box_v;   // No need to copy

  // Solution matrix
  MatrixType testX = B, exactX = B;

  std::cout << "Computing Inverse..." << std::endl;
  { ScopeClock whole_timer("HODLR  Solve: ");

  // For the diagonal boxes of the tree, from leaves to root
  for (int L = source_tree.levels() - 1; L >= 0; --L) {
    ScopeClock timer(std::string("Level ") + std::to_string(L) + ": ");

    for (auto box : boxes(L, source_tree)) {
      // Get the index range of this block
      auto rows = range(box);

      if (box.is_leaf()) {
        // This block is a diagonal leaf
        // -- Perform a dense solve: A(I,I) * X = B
        // -- Apply A(I,I)^{-1} to all U_L(I,_) up the tree

        // Solve AX = B and store into X -- Also stores LU factors into Aii
        MatrixType Aii = A(rows,rows);
        IndexVector ipiv(num_rows(Aii));
        flens::lapack::sv(Aii, ipiv, testX(rows,_));

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == source_tree.root()); box = box.parent()) {
          auto& pU = box_u_inv[box];
          auto urows = shift(rows, box);
          // Solve AX = U and store into U -- Uses previously computed LU
          flens::lapack::trs(NoTrans, Aii, ipiv, pU(urows, _));
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
        auto& U0  = box_u_inv[cbox0];  // box_u_inv[target child 0]
        auto& U1  = box_u_inv[cbox1];  // box_u_inv[target child 1]
        auto& V0T = box_v_inv[cbox0];  // box_v_inv[source child 0]
        auto& V1T = box_v_inv[cbox1];  // box_v_inv[source child 1]

        unsigned r = num_rows(V1T), R = r + num_rows(V0T);
        unsigned c = num_cols(U0),  C = c + num_cols(U1);
        auto rows0 = range(cbox0);
        auto rows1 = range(cbox1);

        // Construct the matrix K = I + [0,V1T;V0T,0] * [U0,0;0,U1]
        MatrixType K = MatrixType(R,C);
        K.diag(0) = 1;
        K(_(  1,r), _(c+1,C)) = V1T * U1;
        K(_(r+1,R), _(  1,c)) = V0T * U0;

        // Apply the SMW to X: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
        MatrixType T(R,num_cols(testX));
        T(_(  1,r),_) = V1T * testX(rows1,_);
        T(_(r+1,R),_) = V0T * testX(rows0,_);
        // Solve KX = T and store into T -- Also stores LU factors into K
        IndexVector ipiv(num_rows(K));
        flens::lapack::sv(K, ipiv, T);
        // Apply U to complete the SMW
        testX(rows0,_) -= U0 * T(_(  1,r),_);
        testX(rows1,_) -= U1 * T(_(r+1,R),_);

        // Apply A^{-1} to the appropriate rows of all of the parent U's
        for ( ; !(box == source_tree.root()); box = box.parent()) {
          auto& pU = box_u_inv[box];
          auto urows0 = shift(rows0, box);
          auto urows1 = shift(rows1, box);

          // Apply the SMW to pU: (I + U*VT)^{-1} = I - U * (I+VT*U)^{-1} * VT
          MatrixType T(R,num_cols(pU));
          T(_(  1,r),_) = V1T * pU(urows1,_);
          T(_(r+1,R),_) = V0T * pU(urows0,_);
          // Solve KX = T and store into T -- Uses previously computed LU
          flens::lapack::trs(NoTrans, K, ipiv, T);
          // Apply U to complete the SMW
          pU(urows0,_) -= U0 * T(_(  1,r),_);
          pU(urows1,_) -= U1 * T(_(r+1,R),_);
        }
      }
    }
  }

  } // end scopeclock

  { ScopeClock timer(std::string("Direct Solve: "));
  flens::lapack::sv(A, IndexVector(), exactX);
  }

  MatrixType ResX = exactX - testX;
  //for (unsigned i = 1; i <= N; ++i) std::cout << i << ":\t" << exactX(i,1) << "\t" << testX(i,1) << "\t" << ResX(i,1) << std::endl;
  std::cout << "Solved residual norm_F = " << frobenius_norm(ResX) << std::endl;
}
