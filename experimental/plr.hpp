
#include <memory>
#include <iterator>

#include "fmmtl/numeric/flens.hpp"
#include "util/Probe.hpp"

#include "fmmtl/numeric/Vec.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/traversal/DualTraversal.hpp"


// A block-compressed representation of a matrix
template <typename T, unsigned DT, unsigned DS>
struct PLR_Matrix {
  using target_tree_type = fmmtl::NDTree<DT>;
  using source_tree_type = fmmtl::NDTree<DS>;

  using target_box_type = typename target_tree_type::box_type;
  using source_box_type = typename source_tree_type::box_type;

  using matrix_type     = flens::GeMatrix<flens::FullStorage<T> >;

  struct dyadic_tree_leaf {
    source_box_type s;    //< Or charge_range for simplicity
    matrix_type U, V;     //< Low-rank approximation of this block
    template <class MatrixU, class MatrixV>
    dyadic_tree_leaf(const source_box_type& _s, MatrixU&& _u, MatrixV&& _v)
        : s(_s), U(std::forward<MatrixU>(_u)), V(std::forward<MatrixV>(_v)) {}
  };

  // Trees representing the dyadic block structure
  // RI Breaking: Required only because source_box/target_box are proxies...
  std::unique_ptr<target_tree_type> target_tree;
  std::unique_ptr<source_tree_type> source_tree;

  // List of dyadic blocks by target box index
  std::vector<std::vector<dyadic_tree_leaf>> leaf_nodes;
  // Count of blocks by level
  std::vector<char> leaf_count;

  template <class MatrixU, class MatrixV>
  void add_leaf(const source_box_type& s, const target_box_type& t,
                MatrixU&& u, MatrixV&& v) {
    leaf_nodes[t.index()].emplace_back(s,
                                       std::forward<MatrixU>(u),
                                       std::forward<MatrixV>(v));
    leaf_count[t.level()] = 1;
  }

  // Y += A*X
  template <typename YA, typename XA>
  void prod_acc(const flens::GeMatrix<XA>& X, flens::GeMatrix<YA>& Y) {
    using MatrixType = flens::GeMatrix<flens::FullStorage<T> >;
    const flens::Underscore<typename MatrixType::IndexType> _;

    // Permute X
    MatrixType pX(X.numRows(), X.numCols());
    for (auto body : bodies(*source_tree))
      pX(body.index()+1,_) = X(body.number()+1,_);
    // Temporary Y
    MatrixType pY(Y.numRows(), Y.numCols());

#pragma omp parallel default(shared)
    for (unsigned level = 0; level < target_tree->levels(); ++level) {
      if (leaf_count[level] == 0)  // Nothing to be done on this level
        continue;
      auto tb_end = target_tree->box_end(level);
      // In parallel over all of the target boxes of this level
#pragma omp for
      for (auto tb_i = target_tree->box_begin(level); tb_i < tb_end; ++tb_i) {
        // The target box
        target_box_type t = *tb_i;
        auto t_range = _(t.body_begin().index()+1, t.body_end().index());
        // The range of associated results
        auto y = pY(t_range, _);
        // For all the source boxes in this tbox-level list
        for (auto& leaf : leaf_nodes[t.index()]) {
          // The range of charges
          auto s = leaf.s;
          auto s_range = _(s.body_begin().index()+1, s.body_end().index());
          auto x = pX(s_range, _);

          // Apply the (t,s) block
          if (num_rows(leaf.U) == 0) {
            y += leaf.V * x;
          } else {
            MatrixType Vx = leaf.V * x;
            y += leaf.U * Vx;
          }
        }
      }
    }

    // Permute and accumulate into Y
    for (auto body : bodies(*target_tree))
      Y(body.number()+1,_) += pY(body.index()+1,_);
  }


  // y += A*x
  template <typename CIterator, typename RIterator>
  void prod_acc(CIterator x_first, CIterator x_last,
                RIterator y_first, RIterator y_last) {
    // TODO: size check
    (void) x_last; (void) y_last;
    using charge_type = typename std::iterator_traits<CIterator>::value_type;
    using result_type = typename std::iterator_traits<RIterator>::value_type;
    // Wrap the range in flens vectors for linear algebra
    // Requires contiguous data
    using charge_ref = flens::ArrayView<charge_type>;
    using result_ref = flens::ArrayView<result_type>;

    // Permute the charges to match the body order in the tree
    auto p_charges = make_body_binding(*source_tree, x_first);
    // Create permuted results, but don't bother initializing to existing
    auto p_results = make_body_binding<result_type>(*target_tree);

#pragma omp parallel default(shared)
    for (unsigned level = 0; level < target_tree->levels(); ++level) {
      if (leaf_count[level] == 0)  // Nothing to be done on this level
        continue;
      auto tb_end = target_tree->box_end(level);
      // In parallel over all of the target boxes of this level
#pragma omp for
      for (auto tb_i = target_tree->box_begin(level); tb_i < tb_end; ++tb_i) {
        // The target box
        target_box_type t = *tb_i;
        // The range of associated results
        auto r = p_results[t];
        flens::DenseVector<result_ref> y = result_ref(r.size(), &*std::begin(r));
        // For all the source boxes in this tbox-level list
        for (auto& leaf : leaf_nodes[t.index()]) {
          // The range of charges
          auto c = p_charges[leaf.s];
          flens::DenseVector<charge_ref> x = charge_ref(c.size(), &*std::begin(c));
          // Apply the (t,s) block
          if (num_rows(leaf.U) == 0) {
            y += leaf.V * x;
          } else {
            y += leaf.U * (leaf.V * x);
          }
        }
      }
    }

    // Copy back permuted results
    auto pri = target_tree->body_permute(y_first);
    for (const auto& ri : p_results) {
      *pri += ri;
      ++pri;
    }
  }
};


/* y += A*x
 */
template <typename PLR_M, typename RangeC, typename RangeR>
void prod_acc(PLR_M& plr, RangeC& x, RangeR& y) {
  return prod_acc(plr, std::begin(x), std::end(x), std::begin(y), std::end(y));
}

/* y += A*x
 */
template <typename PLR_M, typename CIterator, typename RIterator>
void prod_acc(PLR_M& plr,
              CIterator x_first, CIterator x_last,
              RIterator y_first, RIterator y_last) {
  return plr.prod_acc(x_first, x_last, y_first, y_last);
}

/* y += A*x   for DenseVector
 */
template <typename PLR_M, typename XM, typename YM>
void prod_acc(PLR_M& plr, flens::GeMatrix<XM>& X, flens::GeMatrix<YM>& Y) {
  return plr.prod_acc(X, Y);
}


/** Simpler C-like interface to the PLR matrix decomposition
 * @param[in] data    Row-major matrix to compress
 *                      data[i*m + j] represents the i-jth matrix entry.
 * @param[in] n       Number of rows of the matrix.
 * @param[in] m       Number of columns of the matrix
 * @param[in] trgs    Coordinate-major points corresponding to rows of the matrix
 * @param[in] srcs    Coordinate-major points corresponding to cols of the matrix
 * @param[in] max_rank The maximum rank of a block in the PLR structure
 * @param[in] eps_tol  The maximum err tolerance of a block in the PLR structure
 *
 *
 * @tparam DS The number of dimensions in the source data
 * @tparam DT The number of dimensions in the target data
 * @pre size(data) == n*m
 * @pre size(targets) == DT*m
 * @pre size(sources) == DS*n
 *
 * Usage example:
 *    mat = ...;
 *    t = ...; s = ...;
 *    auto plr_matrix = plr_compression<3,2>(mat, 5, 7, t, s);
 *
 *
 *
 */
template <unsigned DT, unsigned DS,
          typename T, typename TT, typename TS>
PLR_Matrix<T,DT,DS>
plr_compression(T* data, unsigned n, unsigned m,
                const TT* trgs, const TS* srcs,
                unsigned max_rank, double eps_tol,
                unsigned init_depth = 0) {
  ScopeClock plr_construction_timer("PLR Matrix Construction: ");

  const Vec<DT,TT>* targets = reinterpret_cast<const Vec<DT,TT>*>(trgs);
  const Vec<DS,TS>* sources = reinterpret_cast<const Vec<DS,TS>*>(srcs);

  // Construct compressed matrix and trees
  using PLRType = PLR_Matrix<T,DT,DS>;
  using target_tree_type = typename PLRType::target_tree_type;
  using source_tree_type = typename PLRType::source_tree_type;

  PLRType plr_m;

  { ScopeClock timer("Trees Complete: ");

  // C++11 overlooked make_unique
  // TODO: Constructor
  plr_m.target_tree =
      std::unique_ptr<target_tree_type>(
          new target_tree_type(targets, targets + n, max_rank));
  plr_m.source_tree =
      std::unique_ptr<source_tree_type>(
          new source_tree_type(sources, sources + m, max_rank));
  plr_m.leaf_nodes.resize(plr_m.target_tree->boxes());
  plr_m.leaf_count.resize(plr_m.target_tree->levels());

  }

  auto& source_tree = *(plr_m.source_tree);
  auto& target_tree = *(plr_m.target_tree);

  // Get the tree types
  using target_box_type  = typename target_tree_type::box_type;
  using source_box_type  = typename source_tree_type::box_type;

  // Permute the sources to match the body order in the tree
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_targets = make_body_binding(target_tree, targets);

  // Define types from FLENS
  using namespace flens;
  using MatrixType     = GeMatrix<FullStorage<T> >;
  using MatrixRefType  = GeMatrix<FullStorageView<T, RowMajor> >;
  const Underscore<typename MatrixType::IndexType>  _;

  // Wrap the data into a FLENS matrix for easy access
  MatrixRefType pA = data_ref_type(n, m, data, m);

  // Permute the matrix data to match the source/target order
  // SLOW WAY
  MatrixType A(n,m);
  unsigned i = 0;
  for (auto tb : bodies(target_tree)) {
    auto Ai = A(++i, _);
    auto pAi = pA(tb.number()+1, _);
    unsigned j = 0;
    for (auto sb : bodies(source_tree)) {
      Ai(++j) = pAi(sb.number()+1);
    }
  }

  std::cout << "Entering traversal" << std::endl;

  // Construct an evaluator for the traversal algorithm
  // See fmmtl::traverse_if documentation
  auto evaluator = [&] (const source_box_type& s, const target_box_type& t) {
    //std::cout << "TargetBox: " << t << std::endl;
    //std::cout << "SourceBox: " << s << std::endl;

    // Compute an approximate SVD of the submatrix defined by the boxes
    // TODO: Replace by fast rank-revealing alg
    unsigned n = t.num_bodies();
    unsigned m = s.num_bodies();
    unsigned r = std::min(n, m);

    auto rows = _(t.body_begin().index()+1, t.body_end().index());
    auto cols = _(s.body_begin().index()+1, s.body_end().index());

    if (max_rank >= r) {
      // Accept the block
      //std::cout << "AUTO ACCEPTED BLOCK, Rank " << r << std::endl;

      // Auto-accept A
      plr_m.add_leaf(s, t, MatrixType(), A(rows,cols));
      // Do not recurse
      return 0;
    }

    // Attempt to factor the block to max_rank
    MatrixType U, VT;
    std::tie(U, VT) = probe_svd(A(rows,cols), max_rank, eps_tol);

    // Accept block and store low-rank approximation or recurse
    if (num_rows(U) != 0) {
      // Accept the block and store low-rank approx
      //std::cout << "ACCEPTED BLOCK, Rank " << num_cols(U) <<
      //    " from " << n << "-by-" << m << std::endl;

      // Store the low rank decomposition of this node
      plr_m.add_leaf(s, t, std::move(U), std::move(VT));
      // Do not recurse further
      return 0;
    } else {
      // Recurse by splitting the source and target boxes (dyadically)
      //std::cout << "REJECTED BLOCK" << std::endl;
      return 3;
    }
  };

  // Perform the traversal, starting at init_depth
  for (auto tbox : boxes(init_depth, target_tree)) {
    for (auto sbox : boxes(init_depth, source_tree)) {
      fmmtl::traverse_if(sbox, tbox, evaluator);
    }
  }

  return plr_m;
}



/** Black Box PLR
 *
 *
 */
template <unsigned DT, unsigned DS,
          typename MatrixFn,
          typename TT, typename TS>
PLR_Matrix<typename MatrixFn::ElementType,DT,DS>
plr_compression(const MatrixFn& M, unsigned n, unsigned m,
                const TT* trgs, const TS* srcs,
                unsigned max_rank, double eps_tol,
                unsigned init_depth = 0) {
  ScopeClock plr_construction_timer("PLR Matrix Construction: ");
  using value_type = typename MatrixFn::ElementType;

  const Vec<DT,TT>* targets = reinterpret_cast<const Vec<DT,TT>*>(trgs);
  const Vec<DS,TS>* sources = reinterpret_cast<const Vec<DS,TS>*>(srcs);

  // Construct compressed matrix and trees
  using PLRType = PLR_Matrix<value_type,DT,DS>;
  using target_tree_type = typename PLRType::target_tree_type;
  using source_tree_type = typename PLRType::source_tree_type;

  PLRType plr_m;

  { ScopeClock timer("Trees Complete: ");

  // C++11 overlooked make_unique
  // TODO: Constructor
  plr_m.target_tree =
      std::unique_ptr<target_tree_type>(
          new target_tree_type(targets, targets + n, max_rank));
  plr_m.source_tree =
      std::unique_ptr<source_tree_type>(
          new source_tree_type(sources, sources + m, max_rank));
  plr_m.leaf_nodes.resize(plr_m.target_tree->boxes());
  plr_m.leaf_count.resize(plr_m.target_tree->levels());

  }

  auto& source_tree = *(plr_m.source_tree);
  auto& target_tree = *(plr_m.target_tree);

  // Get the tree types
  using target_box_type  = typename target_tree_type::box_type;
  using source_box_type  = typename source_tree_type::box_type;

  // Permute the sources/targets to match the body order in the tree
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_targets = make_body_binding(target_tree, targets);

  // Define types from FLENS
  using namespace flens;
  using MatrixType   = GeMatrix<FullStorage<value_type> >;
  using VectorType   = DenseVector<Array<value_type> >;
  using IndexType    = typename MatrixType::IndexType;
  const Underscore<IndexType>  _;

  // Inverse permutation vector of source indices
  DenseVector<Array<IndexType> > ips_idx(m);
  for (auto s : bodies(source_tree))
    ips_idx(s.index()+1) = s.number();

  // Inverse permutation vector of target indices
  DenseVector<Array<IndexType> > ipt_idx(n);
  for (auto t : bodies(target_tree))
    ipt_idx(t.index()+1) = t.number();

  std::cout << "Entering traversal" << std::endl;

  // Construct an evaluator for the traversal algorithm
  // See fmmtl::traverse_if documentation
  auto evaluator = [&] (const source_box_type& s, const target_box_type& t) {
    //std::cout << "TargetBox: " << t << std::endl;
    //std::cout << "SourceBox: " << s << std::endl;

    // Compute an approximate SVD of the submatrix defined by the boxes
    // TODO: Replace by fast rank-revealing alg
    unsigned n = t.num_bodies();
    unsigned m = s.num_bodies();
    unsigned r = std::min(n, m);

    auto rows = _(t.body_begin().index()+1, t.body_end().index());
    auto cols = _(s.body_begin().index()+1, s.body_end().index());
    auto prows = ipt_idx(rows);
    auto pcols = ips_idx(cols);

    if (max_rank >= r) {
      // Accept the block
      //std::cout << "AUTO ACCEPTED BLOCK, Rank " << r << std::endl;

      // Auto-accept the block
      plr_m.add_leaf(s, t, MatrixType(), M(prows,pcols));
      // Do not recurse
      return 0;
    }

    // Attempt to factor the block to max_rank
    // INLINE PROBE SVD WITH CALLBACK FOR MATVEC
    unsigned rc = std::min(max_rank + 10, r);

    // Construct a random matrix
    MatrixType R1(n, rc);
    fillRandom(R1);

    // Factor A^T * R1 (which is m x rc)   OVERWRITE MTR1?
    MatrixType MTR1(m,rc);
    mm(Trans, M, R1, MTR1, prows, pcols);   // Dispatch to mm callback
    MatrixType U1(m,rc), VT2(rc,rc);
    VectorType D(rc);
    flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::None,
                       MTR1, D, U1, VT2);

    // Factor A * U1 (which is n x rc)
    // Reuse VT and D
    MatrixType MU1(n,rc);
    mm(NoTrans, M, U1, MU1, prows, pcols);    // Dispatch to mm callback
    MatrixType U2(n,rc);
    flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::Save,
                       MU1, D, U2, VT2);

    // Find the eps-rank
    while (D(rc)/D(1) < eps_tol) --rc;

    MatrixType U, VT;
    if (rc <= max_rank) {
      auto r = _(1,rc);
      const DiagMatrix<ConstArrayView<double> > DM = D(r);
      U  = U2(_,r) * DM;
      VT = VT2(r,r) * conjTrans(U1(_,r));
    } else {
      U = VT = MatrixType();
    }

    // end inline probe SVD

    // Accept block and store low-rank approximation or recurse
    if (num_rows(U) != 0) {
      // Accept the block and store low-rank approx
      //std::cout << "ACCEPTED BLOCK, Rank " << num_cols(U) <<
      //    " from " << n << "-by-" << m << std::endl;

      // Store the low rank decomposition of this node
      plr_m.add_leaf(s, t, std::move(U), std::move(VT));
      // Do not recurse further
      return 0;
    } else {
      // Recurse by splitting the source and target boxes (dyadically)
      //std::cout << "REJECTED BLOCK" << std::endl;
      return 3;
    }
  };

  // Perform the traversal, starting at init_depth
  for (auto tbox : boxes(init_depth, target_tree)) {
    for (auto sbox : boxes(init_depth, source_tree)) {
      fmmtl::traverse_if(sbox, tbox, evaluator);
    }
  }

  return plr_m;
}
