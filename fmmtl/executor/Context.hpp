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

template <class Expansion>
class ExpansionContext {
  Expansion E_;

 public:
  FMMTL_IMPORT_EXPANSION_TRAITS(Expansion);

  template <typename Options>
  ExpansionContext(const Expansion& E,
                   Options&)
      : E_(E) {
  }

  Expansion& expansion() {
    return E_;
  }
  const Expansion& expansion() const {
    return E_;
  }
  typename Expansion::kernel_type& kernel() {
    return E_;
  }
  const typename Expansion::kernel_type& kernel() const {
    return E_;
  }
};


// General TreeContxt declarations
template <typename Expansion,
          typename SourceTree,
          typename TargetTree = SourceTree>
class DualTreeContext;

template <typename Expansion,
          typename Tree>
class SingleTreeContext;


/** @struct SingleTreeContext
 * Single tree context specialized for an NDTree
 */
template <typename Expansion>
class SingleTreeContext<Expansion,
                        NDTree<ExpansionTraits<Expansion>::dimension> > {
 public:
  FMMTL_IMPORT_EXPANSION_TRAITS(Expansion);

  typedef NDTree<ExpansionTraits<Expansion>::dimension> tree_type;
  FMMTL_IMPORT_TREEPAIR_TRAITS(tree_type, tree_type);

 protected:
  template <typename Iter>
  using stree_permute_iterator =
      typename source_tree_type::template body_permuted_iterator<Iter>::type;
  template <typename Iter>
  using ttree_permute_iterator =
      typename target_tree_type::template body_permuted_iterator<Iter>::type;

  //! The tree of sources and targets
  source_tree_type source_tree_;

  // TODO: If std::is_same<source_type,point_type>,
  // then we don't need to store sources, the source tree is doing it for us

  //! The sources associated with bodies in the source_tree_
  typedef const std::vector<source_type> source_container;
  typedef typename source_container::const_iterator source_container_iterator;
  source_container sources_;
  //! The targets associated with bodies in the target_tree_
  typedef source_container target_container;
  typedef source_container_iterator target_container_iterator;

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
  template <typename SourceIter, typename Options>
  SingleTreeContext(SourceIter sfirst, SourceIter slast,
                    Options& opts)
      : source_tree_(sfirst, slast, opts.ncrit),
        sources_(sfirst, slast) {
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

  // Accessor to the source of a source body
  typedef stree_permute_iterator<source_container_iterator> source_iterator;
  inline source_iterator source(const source_body_iterator& sbi) const {
    return source_tree_permute(sources_.cbegin(), sbi);
  }
  // Accessor to the target of a target body
  typedef ttree_permute_iterator<target_container_iterator> target_iterator;
  inline target_iterator target(const target_body_iterator& tbi) const {
    return target_tree_permute(sources_.cbegin(), tbi);
  }
};


/** @struct DualTreeContext
 * Dual tree context specialized for NDTree trees
 */
template <typename Expansion>
class DualTreeContext<Expansion,
                      NDTree<ExpansionTraits<Expansion>::dimension> > {
 public:
  FMMTL_IMPORT_EXPANSION_TRAITS(Expansion);

  typedef NDTree<ExpansionTraits<Expansion>::dimension> tree_type;
  FMMTL_IMPORT_TREEPAIR_TRAITS(tree_type, tree_type);

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

  // TODO: If std::is_same<source_type,point_type>,
  // then we don't need to store sources, the source tree is doing it for us

  //! The sources associated with bodies in the source_tree_
  typedef const std::vector<source_type> source_container;
  typedef typename source_container::const_iterator source_container_iterator;
  source_container sources_;
  //! The targets associated with bodies in the target_tree_
  typedef const std::vector<target_type> target_container;
  typedef typename target_container::const_iterator target_container_iterator;
  target_container targets_;

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
  template <typename SourceIter, typename TargetIter, typename Options>
  DualTreeContext(SourceIter sfirst, SourceIter slast,
                  TargetIter tfirst, TargetIter tlast,
                  Options& opts)
      : source_tree_(sfirst, slast, opts.ncrit),
        target_tree_(tfirst, tlast, opts.ncrit),
        sources_(sfirst, slast),
        targets_(tfirst, tlast) {
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

  // Accessor to the source of a source body
  typedef stree_permute_iterator<source_container_iterator> source_iterator;
  inline source_iterator source(const source_body_iterator& sbi) const {
    return source_tree_permute(sources_.cbegin(), sbi);
  }
  // Accessor to the target of a target body
  typedef ttree_permute_iterator<target_container_iterator> target_iterator;
  inline target_iterator target(const target_body_iterator& tbi) const {
    return target_tree_permute(targets_.cbegin(), tbi);
  }
};



template <typename ExpansionContext,
          typename TreeContext>
class DataContext
    : public ExpansionContext,
      public TreeContext {
 public:
  FMMTL_IMPORT_EXPANSION_TRAITS(typename ExpansionContext::expansion_type);
  FMMTL_IMPORT_TREEPAIR_TRAITS(typename TreeContext::source_tree_type,
                               typename TreeContext::target_tree_type);
 private:
  //! Multipole expansions corresponding to Box indices in Tree
  typedef std::vector<multipole_type> multipole_container;
  multipole_container M_;
  //! Local expansions corresponding to Box indices in Tree
  typedef std::vector<local_type> local_container;
  local_container L_;

  //! Iterator to the start of the charge vector
  typedef std::vector<charge_type> charge_container;
  typedef typename charge_container::const_iterator charge_container_iterator;
  charge_container_iterator c_;
  //! Iterator to the start of the result vector
  typedef std::vector<result_type> result_container;
  typedef typename result_container::iterator result_container_iterator;
  result_container_iterator r_;

 public:
  template <typename SourceIter, typename TargetIter,
            typename Options>
  DataContext(const expansion_type& E,
              SourceIter sfirst, SourceIter slast,
              TargetIter tfirst, TargetIter tlast,
              Options& opts)
      : ExpansionContext(E, opts),
        TreeContext(sfirst, slast, tfirst, tlast, opts),
        M_(this->source_tree().boxes()),
        L_(this->target_tree().boxes()) {
  }

  template <typename SourceIter,
            typename Options>
  DataContext(const expansion_type& E,
              SourceIter sfirst, SourceIter slast,
              Options& opts)
      : ExpansionContext(E, opts),
        TreeContext(sfirst, slast, opts),
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

  // Define the body data iterators
  typedef typename TreeContext::source_iterator source_iterator;
  typedef typename TreeContext::template stree_permute_iterator<charge_container_iterator> charge_iterator;
  typedef typename TreeContext::target_iterator target_iterator;
  typedef typename TreeContext::template ttree_permute_iterator<result_container_iterator> result_iterator;

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
