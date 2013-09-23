#pragma once
/** @file Context
 *
 * A class which stores the Kernel, tree(s), and data.
 */

/** TODO: Write Context and Tree concepts */

#include "fmmtl/KernelTraits.hpp"
#include "fmmtl/TreeTraits.hpp"
#include "fmmtl/FMMOptions.hpp"

#include <boost/iterator/transform_iterator.hpp>

#include "fmmtl/meta/dimension.hpp"

// General TreeContext declarations
template <typename TreeType>
class SingleTreeContext;

template <typename SourceTreeType,
          typename TargetTreeType>
class DualTreeContext;


/** @struct SingleTreeContext
 * Single tree context specialized for an NDTree
 */
template <unsigned DIM>
class SingleTreeContext<NDTree<DIM> > {
 public:
  typedef NDTree<DIM> source_tree_type;
  typedef NDTree<DIM> target_tree_type;
  FMMTL_IMPORT_TREEPAIR_TRAITS(source_tree_type, target_tree_type);

 protected:
  template <typename Iter>
  using stree_permute_iterator =
      typename source_tree_type::template body_permuted_iterator<Iter>::type;
  template <typename Iter>
  using ttree_permute_iterator =
      typename target_tree_type::template body_permuted_iterator<Iter>::type;

  //! The tree of sources and targets
  source_tree_type source_tree_;

  /** Permute iterators */
  template <typename Iterator>
  inline stree_permute_iterator<Iterator>
  source_tree_permute(Iterator it, const source_body_iterator& sbi) const {
    return source_tree().body_permute(it, sbi);
  }
  template <typename Iterator>
  inline ttree_permute_iterator<Iterator>
  target_tree_permute(Iterator it, const target_body_iterator& tbi) const {
    return target_tree().body_permute(it, tbi);
  }

 public:
  //! Constructor
  template <typename KernelMatrix, typename Options>
  SingleTreeContext(const KernelMatrix& mat, Options& opts)
      : source_tree_(mat.sources().begin(), mat.sources().end(), opts.ncrit) {
  }

  // Tree accessors
  inline source_tree_type& source_tree() {
    return source_tree_;
  }
  inline const source_tree_type& source_tree() const {
    return source_tree_;
  }
  inline target_tree_type& target_tree() {
    return source_tree_;
  }
  inline const target_tree_type& target_tree() const {
    return source_tree_;
  }
};


/** @struct DualTreeContext
 * Dual tree context specialized for two NDTree trees
 */
template <unsigned SOURCEDIM, unsigned TARGETDIM>
class DualTreeContext<NDTree<SOURCEDIM>,
                      NDTree<TARGETDIM> > {
 public:
  typedef NDTree<SOURCEDIM> source_tree_type;
  typedef NDTree<TARGETDIM> target_tree_type;
  FMMTL_IMPORT_TREEPAIR_TRAITS(source_tree_type, target_tree_type);

 protected:
  template <typename Iter>
  using stree_permute_iterator =
      typename source_tree_type::template body_permuted_iterator<Iter>::type;
  template <typename Iter>
  using ttree_permute_iterator =
      typename target_tree_type::template body_permuted_iterator<Iter>::type;

  //! The tree of sources
  source_tree_type source_tree_;
  //! The tree of targets
  target_tree_type target_tree_;

  /** Permuted iterators */
  template <typename Iterator>
  inline stree_permute_iterator<Iterator>
  source_tree_permute(Iterator it, const source_body_iterator& sbi) const {
    return source_tree().body_permute(it, sbi);
  }
  template <typename Iterator>
  inline ttree_permute_iterator<Iterator>
  target_tree_permute(Iterator it, const target_body_iterator& tbi) const {
    return target_tree().body_permute(it, tbi);
  }

 public:
  //! Constructor
  template <typename KernelMatrix, typename Options>
  DualTreeContext(const KernelMatrix& mat, Options& opts)
      : source_tree_(mat.sources().begin(), mat.sources().end(), opts.ncrit),
        target_tree_(mat.targets().begin(), mat.targets().end(), opts.ncrit) {
  }

  inline source_tree_type& source_tree() {
    return source_tree_;
  }
  inline const source_tree_type& source_tree() const {
    return source_tree_;
  }
  inline target_tree_type& target_tree() {
    return target_tree_;
  }
  inline const target_tree_type& target_tree() const {
    return target_tree_;
  }
};



template <typename KernelMatrix,
          typename TreeContext>
class DataContext
    : public TreeContext {
 public:
  typedef KernelMatrix kernel_matrix_type;
  FMMTL_IMPORT_EXPANSION_TRAITS(typename kernel_matrix_type::expansion_type);

  FMMTL_IMPORT_TREEPAIR_TRAITS(typename TreeContext::source_tree_type,
                               typename TreeContext::target_tree_type);
 private:
  //! The kernel matrix this context is built for
  const kernel_matrix_type& mat_;

  //! Source and target iterator types in the kernel_matrix
  typedef typename kernel_matrix_type::source_array::const_iterator source_container_iterator;
  typedef typename kernel_matrix_type::target_array::const_iterator target_container_iterator;

  //! The "multipole acceptance criteria" to decide which boxes to interact
  std::function<bool(const source_box_type&, const target_box_type&)> mac_;

  //! Iterator to the start of the charge vector
  typedef std::vector<charge_type> charge_container;
  typedef typename charge_container::const_iterator charge_container_iterator;
  charge_container_iterator c_;
  //! Iterator to the start of the result vector
  typedef std::vector<result_type> result_container;
  typedef typename result_container::iterator result_container_iterator;
  result_container_iterator r_;

  //! Multipole expansions corresponding to Box indices in Tree
  typedef std::vector<multipole_type> multipole_container;
  multipole_container M_;
  //! Local expansions corresponding to Box indices in Tree
  typedef std::vector<local_type> local_container;
  local_container L_;

 public:
  template <class Options>
  DataContext(const kernel_matrix_type& mat, Options& opts)
      : TreeContext(mat, opts),
        mat_(mat),
        mac_(opts.MAC()),
        // TODO: only allocate if used...
        M_(this->source_tree().boxes()),
        L_(this->target_tree().boxes()) {
  }

  template <typename Executor>
  inline void execute(const std::vector<charge_type>& charges,
                      std::vector<result_type>& results,
                      Executor* exec) {
    c_ = charges.begin();
    r_ = results.begin();

    exec->execute(*this);
  }

  const expansion_type& expansion() const {
    return mat_.expansion();
  }
  const kernel_type& kernel() const {
    return mat_.kernel();
  }

  // Accessors to the multipole expansion of a source box
  inline multipole_type& multipole(const source_box_type& box) {
    return M_[box.index()];
  }
  inline const multipole_type& multipole(const source_box_type& box) const {
    return M_[box.index()];
  }
  // Accessors to the local expansion of a target box
  inline local_type& local(const target_box_type& box) {
    return L_[box.index()];
  }
  inline const local_type& local(const target_box_type& box) const {
    return L_[box.index()];
  }

  // Accept or reject the interaction of this source-target box pair
  inline bool mac(const source_box_type& sbox,
                  const target_box_type& tbox) const {
    return mac_(sbox, tbox);
  }

  // Define the body data iterators
  typedef typename TreeContext::template stree_permute_iterator<source_container_iterator> source_iterator;
  typedef typename TreeContext::template stree_permute_iterator<charge_container_iterator> charge_iterator;
  typedef typename TreeContext::template ttree_permute_iterator<target_container_iterator> target_iterator;
  typedef typename TreeContext::template ttree_permute_iterator<result_container_iterator> result_iterator;

  // Accessor to the source data of a source body
  inline source_iterator source(const source_body_iterator& sbi) const {
    return this->source_tree_permute(mat_.sources().begin(), sbi);
  }
  // Convenience methods for the sources
  inline source_iterator source_begin(const source_box_type& b) const {
    return this->source(b.body_begin());
  }
  inline source_iterator source_end(const source_box_type& b) const {
    return this->source(b.body_end());
  }
  inline source_iterator source_begin() const {
    return this->source(this->source_tree().body_begin());
  }
  inline source_iterator source_end() const {
    return this->source(this->source_tree().body_end());
  }

  // Accessor to the charge of a source body
  inline charge_iterator charge(const source_body_iterator& sbi) const {
    return this->source_tree_permute(c_, sbi);
  }
  // Convenience methods for charges
  inline charge_iterator charge_begin(const source_box_type& b) const {
    return this->charge(b.body_begin());
  }
  inline charge_iterator charge_end(const source_box_type& b) const {
    return this->charge(b.body_end());
  }
  inline charge_iterator charge_begin() const {
    return this->charge(this->source_tree().body_begin());
  }
  inline charge_iterator charge_end() const {
    return this->charge(this->source_tree().body_end());
  }

  // Accessor to the target data of a target body
  inline target_iterator target(const target_body_iterator& sbi) const {
    return this->target_tree_permute(mat_.targets().begin(), sbi);
  }
  // Convenience methods for the targets
  inline target_iterator target_begin(const target_box_type& b) const {
    return this->target(b.body_begin());
  }
  inline target_iterator target_end(const target_box_type& b) const {
    return this->target(b.body_end());
  }
  inline target_iterator target_begin() const {
    return this->target(this->target_tree().body_begin());
  }
  inline target_iterator target_end() const {
    return this->target(this->target_tree().body_end());
  }

  // Accessor to the result of a target body
  inline result_iterator result(const target_body_iterator& tbi) const {
    return this->target_tree_permute(r_, tbi);
  }
  // Convenience methods for results
  inline result_iterator result_begin(const target_box_type& b) const {
    return this->result(b.body_begin());
  }
  inline result_iterator result_end(const target_box_type& b) const {
    return this->result(b.body_end());
  }
  inline result_iterator result_begin() const {
    return this->result(this->target_tree().body_begin());
  }
  inline result_iterator result_end() const {
    return this->result(this->target_tree().body_end());
  }
};
