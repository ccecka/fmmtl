
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/executor/Traversal.hpp"

#include "fmmtl/Direct.hpp"

#include "fmmtl/util/Clock.hpp"

#include "fmmtl/numeric/mtl.hpp"



struct MyKernel {
  typedef double value_type;

  typedef Vec<2,value_type> source_type;
  typedef Vec<2,value_type> target_type;
  typedef value_type        charge_type;
  typedef value_type        result_type;

  typedef value_type        kernel_value_type;

  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    return std::exp(-norm_2(s - t));
  }
};



template <typename Kernel>
class kernel_matrix_functor {
  typedef kernel_matrix_functor                              self_type;
 public:
  typedef std::size_t                                        size_type;
  typedef typename KernelTraits<Kernel>::kernel_value_type   result_type;

  typedef typename KernelTraits<Kernel>::source_type         source_type;
  typedef typename KernelTraits<Kernel>::target_type         target_type;

  typedef Kernel                                             kernel_type;

  // TODO: Abstract range types
  typedef std::vector<source_type> source_container;
  typedef std::vector<target_type> target_container;

  kernel_matrix_functor()
      : K_(), s_(), t_() {}
  kernel_matrix_functor(const Kernel& K)
      : K_(K), s_(), t_() {}
  template <typename TargetRange, typename SourceRange>
  kernel_matrix_functor(const Kernel& K,
                        const TargetRange& t,
                        const SourceRange& s)
      : K_(K), t_(std::begin(t), std::end(t)), s_(std::begin(s), std::end(s)) {}

  friend size_type inline num_rows(const self_type& A) { return A.t_.size(); }
  friend size_type inline num_cols(const self_type& A) { return A.s_.size(); }

  result_type operator()(size_type r, size_type c) const {
    return K_(t_[r], s_[c]);
  }


  kernel_type&            kernel()        { return K_; }
  kernel_type const&      kernel()  const { return K_; }
  source_container&       sources()       { return s_; }
  source_container const& sources() const { return s_; }
  target_container&       targets()       { return t_; }
  target_container const& targets() const { return t_; }

 private:
  Kernel K_;
  target_container t_;
  source_container s_;
};





int main(int argc, char** argv) {
  int N = 200;
  int M = 200;
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
  static constexpr std::size_t max_rank = 5;
  static constexpr double      eps_tol  = 1e-5;

  // Define an (mtl) implicit matrix to work with
  using kernel_type = MyKernel;
  using kernel_matrix_f_type = kernel_matrix_functor<kernel_type>;
  using kernel_matrix_type = mtl::mat::implicit_dense<kernel_matrix_f_type>;

  using source_type = typename kernel_type::source_type;
  using target_type = typename kernel_type::target_type;
  using charge_type = typename kernel_type::charge_type;
  using result_type = typename kernel_type::result_type;

  using kernel_value_type = typename kernel_type::kernel_value_type;

  // Construct a set of sources and targets
  std::vector<source_type> sources = fmmtl::random_n(M);
  std::vector<target_type> targets = fmmtl::random_n(N);

  // Construct a set of sources and targets
  std::vector<charge_type> charges = fmmtl::random_n(M);
  std::vector<result_type> results(M);

  // Construct the source and target trees
  constexpr unsigned SD = fmmtl::dimension<source_type>::value;
  constexpr unsigned TD = fmmtl::dimension<target_type>::value;
  fmmtl::NDTree<SD> source_tree(sources, max_rank);
  fmmtl::NDTree<TD> target_tree(targets, max_rank);

  std::cout << "Trees complete" << std::endl;

  typedef typename fmmtl::NDTree<TD>::box_type target_box_type;
  typedef typename fmmtl::NDTree<SD>::box_type source_box_type;
  typedef typename fmmtl::NDTree<TD>::body_type target_body_type;
  typedef typename fmmtl::NDTree<SD>::body_type source_body_type;

  // Permute the sources and charges to match the body order in the tree
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_charges = make_body_binding(source_tree, charges);

  // Permute the targets and results to match the body order in the tree
  auto p_targets = make_body_binding(target_tree, targets);
  auto p_results = make_body_binding(target_tree, results);

  // Construct a kernel and matrix wrapper
  kernel_type kernel;
  kernel_matrix_f_type kernel_matrix_f(kernel, p_targets, p_sources);
  kernel_matrix_type kernel_matrix(kernel_matrix_f);

  // Evaluate the entire kernel matrix (H-Matrix algebraic system)
  using matrix_type = mtl::matrix<kernel_value_type>;
  matrix_type A(kernel_matrix);


  std::cout << "Entering traversal" << std::endl;


  // TODO: Still need a compressed (+parallel) storage type for dyadic boxes...
  struct dyadic_tree_leaf {
    source_box_type s;    //< Or charge_range for simplicity
    target_box_type t;    //< Or result_range for simplicity
    matrix_type U, V;     //< Low-rank approximation of this block
    dyadic_tree_leaf(const source_box_type& _s, const target_box_type& _t,
                     matrix_type&& _u, matrix_type&& _v)
        : s(_s), t(_t),
          U(std::forward<matrix_type>(_u)), V(std::forward<matrix_type>(_v)) {}
  };
  std::vector<dyadic_tree_leaf> leaf_node;


  // Construct an evaluator for the traversal algorithm
  // See fmmtl::traverse_if documentation
  auto evaluator = [&] (const source_box_type& s, const target_box_type& t) {
    std::cout << "TargetBox: " << t << std::endl;
    std::cout << "SourceBox: " << s << std::endl;

    // Construct the submatrix defined by the source and target box
    mtl::irange rows(t.body_begin().index(), t.body_end().index());
    mtl::irange cols(s.body_begin().index(), s.body_end().index());
    // Clone shouldn't be required... stupid MTL/SVD
    matrix_type Ats = mtl::clone(A[rows][cols]);

    // Compute an approximate SVD   TODO: Replace by fast rank-revealing alg
    matrix_type U, D, VT;
    boost::tie(U, D, VT) = svd(Ats, eps_tol);

    if (num_rows(D) <= max_rank || num_cols(D) <= max_rank ||
        D[max_rank][max_rank] < eps_tol) {
      std::cout << "ACCEPTED BLOCK" << std::endl;

      // Store the low rank decomposition of this node
      const unsigned r = std::min(max_rank, std::min(num_rows(D), num_cols(D)));
      U.change_dim(num_rows(D), r, true);
      VT.change_dim(num_cols(D), r, true);
      D.change_dim(r, r, true);

      leaf_node.emplace_back(s,t,std::move(U),matrix_type(D*trans(VT)));
      // Do not recurse further
      return 0;
    } else {
      // Recurse by splitting the source and target boxes (dyadically)
      return 3;
    }
  };

  // Perform the traversal
  fmmtl::traverse_if(source_tree.root(), target_tree.root(), evaluator);


  // Perform the matvec
  { ScopeClock timer("PLR MatVec: ");

  // Evaluate the leaf low-rank decompositions
  for (auto&& leaf : leaf_node) {
    // Get the range of charges and results associated with the boxes
    auto c = p_charges[leaf.s];
    auto r = p_results[leaf.t];

    // Wrap the range in mtl::vectors for linear algebra
    // Requires contiguous data
    mtl::vector<charge_type> cv(c.size(), &(*(c.begin())));
    mtl::vector<result_type> rv(r.size(), &(*(r.begin())));

    // Compute the product
    rv += leaf.U * mtl::vector<kernel_value_type>(leaf.V * cv);
  }

  } // timer


  // Copy back permuted results
  auto pri = target_tree.body_permute(results.begin(), target_tree.body_begin());
  for (auto&& ri : p_results) {
    *pri = ri;
    ++pri;
  }


  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M);
    { ScopeClock timer("Direct: ");
      fmmtl::direct(kernel, sources, charges, targets, exact);
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
    double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / results.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }
}
