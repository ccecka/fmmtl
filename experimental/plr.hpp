
#include <memory>

#include "fmmtl/numeric/flens.hpp"

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/executor/Traversal.hpp"

#include "fmmtl/Direct.hpp"

// A block-compressed representation of a matrix
template <typename T, unsigned DT, unsigned DS>
struct PLR_Matrix {
  // Avoid copy/move of a tree for now...
  fmmtl::NDTree<DT>* target_tree;
  fmmtl::NDTree<DS>* source_tree;

  ~PLR_Matrix() {
    delete target_tree; target_tree = nullptr;
    delete source_tree; source_tree = nullptr;
  }

  using target_box_type = typename fmmtl::NDTree<DT>::box_type;
  using source_box_type = typename fmmtl::NDTree<DS>::box_type;

  using matrix_type     = flens::GeMatrix<flens::FullStorage<T> >;

  // TODO: Still need a compressed (+parallel) storage type for dyadic boxes...
  struct dyadic_tree_leaf {
    source_box_type s;    //< Or charge_range for simplicity
    target_box_type t;    //< Or result_range for simplicity
    matrix_type U, V;     //< Low-rank approximation of this block
    template <class MatrixU, class MatrixVT>
    dyadic_tree_leaf(const source_box_type& _s, const target_box_type& _t,
                     MatrixU&& _u, MatrixVT&& _v)
        : s(_s), t(_t),
          U(std::forward<MatrixU>(_u)), V(std::forward<MatrixVT>(_v)) {}
  };

  // List of blocks, could use a better data structure...
  std::vector<dyadic_tree_leaf> leaf_node;
};


/* y += A*x
 */
template <typename PLR_M, typename RangeC, typename RangeR>
void prod_acc(PLR_M& plr, RangeC&& x, RangeR& y) {
  return prod_acc(plr, std::begin(x), std::end(x), std::begin(y), std::end(y));
}

/* y += A*x
 */
template <typename PLR_M, typename CIterator, typename RIterator>
void prod_acc(PLR_M& plr,
              CIterator x_first, CIterator x_last,
              RIterator y_first, RIterator y_last) {
  // TODO: size check
  (void) x_last;
  (void) y_last;
  using charge_type = typename std::iterator_traits<CIterator>::value_type;
  using result_type = typename std::iterator_traits<RIterator>::value_type;
  // Wrap the range in flens vectors for linear algebra
  // Requires contiguous data
  using charge_ref = flens::ArrayView<charge_type>;
  using result_ref = flens::ArrayView<result_type>;

  // Permute the charges to match the body order in the tree
  auto p_charges = make_body_binding(*(plr.source_tree), x_first);
  // Permute the results to match the body order in the tree
  auto p_results = make_body_binding(*(plr.target_tree), y_first);

  // Perform the matvec
  // Evaluate the leaf low-rank decompositions
  for (auto&& leaf : plr.leaf_node) {
    // Get the range of charges and results associated with the boxes
    auto c = p_charges[leaf.s];
    auto r = p_results[leaf.t];

    flens::DenseVector<charge_ref> x = charge_ref(c.size(), &*std::begin(c));
    flens::DenseVector<result_ref> y = result_ref(r.size(), &*std::begin(r));

    // Compute the product
    y += leaf.U * (leaf.V * x);
  }

  // Copy back permuted results
  auto pri = plr.target_tree->body_permute(y_first, plr.target_tree->body_begin());
  for (auto&& ri : p_results) {
    *pri = ri;
    ++pri;
  }

}


/** Simlpler C-like interface to the PLR matrix decomposition
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
 *    auto plr_matrix = plr_compression(mat, 5, 7, t, s);
 *
 *
 *
 */
template <unsigned DT, unsigned DS,
          typename T, typename TT, typename TS>
PLR_Matrix<T,DT,DS>
plr_compression(T* data, unsigned n, unsigned m,
                const TT* trgs, const TS* srcs,
                unsigned max_rank, double eps_tol) {

  const Vec<DT,TT>* targets = reinterpret_cast<const Vec<DT,TT>*>(trgs);
  const Vec<DS,TS>* sources = reinterpret_cast<const Vec<DS,TS>*>(srcs);

  // Construct compressed matrix and trees
  PLR_Matrix<T,DT,DS> plr_m;
  plr_m.source_tree = new fmmtl::NDTree<DS>(sources, sources + m, max_rank);
  plr_m.target_tree = new fmmtl::NDTree<DT>(targets, targets + n, max_rank);

  std::cout << "Trees complete" << std::endl;

  // Get the tree types
  typedef typename fmmtl::NDTree<DT>::box_type target_box_type;
  typedef typename fmmtl::NDTree<DS>::box_type source_box_type;
  typedef typename fmmtl::NDTree<DT>::body_type target_body_type;
  typedef typename fmmtl::NDTree<DS>::body_type source_body_type;

  // Permute the sources to match the body order in the tree
  auto p_sources = make_body_binding(*(plr_m.source_tree), sources);

  // Permute the targets to match the body order in the tree
  auto p_targets = make_body_binding(*(plr_m.target_tree), targets);

  // Define types from FLENS
  using namespace flens;
  using matrix_type     = GeMatrix<FullStorage<T> >;
  using vector_type     = DenseVector<Array<double> >;
  using data_ref_type   = FullStorageView<T, RowMajor>;
  using matrix_ref_type = GeMatrix<data_ref_type>;
  const Underscore<typename matrix_type::IndexType>  _;

  // Wrap the data into a FLENS matrix
  matrix_ref_type pA = data_ref_type(n, m, data, m);

  // Permute the matrix data to match the source/target order
  // SLOW WAY
  matrix_type A(n,m);
  unsigned i = 0;
  for (target_body_type tb : bodies(*(plr_m.target_tree))) {
    auto Ai = A(++i, _);
    auto pAi = pA(tb.number()+1, _);
    unsigned j = 0;
    for (source_body_type sb : bodies(*(plr_m.source_tree))) {
      Ai(++j) = pAi(sb.number()+1);
    }
  }

  std::cout << "Entering traversal" << std::endl;

  // Construct an evaluator for the traversal algorithm
  // See fmmtl::traverse_if documentation
  auto evaluator = [&] (const source_box_type& s, const target_box_type& t) {
    std::cout << "TargetBox: " << t << std::endl;
    std::cout << "SourceBox: " << s << std::endl;

    // Compute an approximate SVD of the submatrix defined by the boxes
    // TODO: Replace by fast rank-revealing alg
    unsigned n = t.num_bodies();
    unsigned m = s.num_bodies();
    unsigned r = std::min(n, m);

    matrix_type U(n,r), VT(r,m);
    vector_type D(r);
    auto rows = _(t.body_begin().index()+1, t.body_end().index());
    auto cols = _(s.body_begin().index()+1, s.body_end().index());
    // A is mutated by SVD so we need a copy
    matrix_type A_ST = A(rows,cols);
    flens::lapack::svd(flens::lapack::SVD::Save, flens::lapack::SVD::Save,
                       A_ST, D, U, VT);

    // Find the rank-cutoff, if exists
    unsigned rank = 1;
    while (rank <= r && D(rank)/D(1) > eps_tol) ++rank;
    --rank;

    // If either of the boxes are leaves, force-accept
    if (s.is_leaf() || t.is_leaf())
      rank = std::min(rank, std::min(max_rank, r));

    // Accept block and store low-rank approximation or recurse
    if (rank <= max_rank) {
      // Accept the block and store low-rank approx
      std::cout << "ACCEPTED BLOCK, Rank " << rank << std::endl;

      auto _rank = _(1,rank);

      const DiagMatrix<ConstArrayView<double> > DM = D(_rank);

      // Store the low rank decomposition of this node
      plr_m.leaf_node.emplace_back(s, t,
                                   U(_(1,n), _rank),
                                   DM * VT(_rank, _(1,m)));
      // Do not recurse further
      return 0;
    } else {
      // Recurse by splitting the source and target boxes (dyadically)
      std::cout << "REJECTED BLOCK, Rank " << rank << std::endl;
      return 3;
    }

  };

  // Perform the traversal
  fmmtl::traverse_if(plr_m.source_tree->root(),
                     plr_m.target_tree->root(),
                     evaluator);

  std::cout << "PLR_Matrix RETURNING" << std::endl;

  return plr_m;
}
