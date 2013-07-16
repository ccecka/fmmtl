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
  ExpansionContext(const Expansion& E)
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

// TODO: Better factorization:
// ExpansionContext, TreeContext, MACContext, DataContext, etc

template <typename Expansion,
          typename SourceTree,
          typename TargetTree = SourceTree>
class DualTreeContext;

template <typename Expansion,
          typename Tree>
class SingleTreeContext;

/** @class DualTreeContext
 * @brief A very general TreeContext class holds a tree and endows it with
 * expansion-specific information.
 * This provides a context to any tree that provides the following interface:
 * TODO
 */
template <typename Expansion>
class DualTreeContext<Expansion,
                      Octree<typename ExpansionTraits<Expansion>::point_type> >
    : public ExpansionContext<Expansion>
{
 public:
  typedef ExpansionTraits<Expansion> expansion_traits;
  FMMTL_IMPORT_EXPANSION_TRAITS(expansion_traits);

  typedef TreePairTraits<Octree<point_type>, Octree<point_type>> tree_traits;
  FMMTL_IMPORT_TREEPAIR_TRAITS(tree_traits);

 protected:
  // Transform a body to a value
  template <typename BodyIter, typename Indexable>
  struct Transformer {
    typedef typename Indexable::reference result_type;
    typedef typename BodyIter::value_type argument_type;
    Transformer(const Indexable& value) : value_(value) {}
    result_type operator()(const argument_type& body) const {
      return value_[body.number()]; // TODO: TEMP to avoid permutation
    }
   private:
    Indexable value_;
  };

  template <typename BodyIter, typename Indexable>
  using body_transform_iterator =
      boost::transform_iterator<Transformer<BodyIter,Indexable>,BodyIter>;

  template <typename BodyIter, typename Indexable>
  inline body_transform_iterator<BodyIter, Indexable>
  transform_iterator(BodyIter bi, Indexable v) const {
    return boost::make_transform_iterator(bi,
                                          Transformer<BodyIter,Indexable>(v));
  }

  //! The tree of sources
  source_tree_type source_tree_;
  //! The tree of targets
  target_tree_type target_tree_;

  //! Multipole expansions corresponding to Box indices in Tree
  typedef std::vector<multipole_type> multipole_container;
  multipole_container M_;
  //! Local expansions corresponding to Box indices in Tree
  typedef std::vector<local_type> local_container;
  local_container L_;

  // TODO: If std::is_same<source_type,point_type>,
  // then we don't need to store sources, the source tree is doing it for us

  //! The sources associated with bodies in the source_tree_
  typedef const std::vector<source_type> source_container;
  source_container sources_;
  //! The targets associated with bodies in the target_tree_
  typedef const std::vector<target_type> target_container;
  target_container targets_;

  //! Iterator to the start of the source vector
  typedef typename source_container::const_iterator source_iterator;
  source_iterator s_;
  //! Iterator to the start of the charge vector
  typedef std::vector<charge_type> charge_container;
  typedef typename charge_container::const_iterator charge_iterator;
  charge_iterator c_;
  //! Iterator to the start of the target vector
  typedef typename target_container::const_iterator target_iterator;
  target_iterator t_;
  //! Iterator to the start of the result vector
  typedef std::vector<result_type> result_container;
  typedef typename result_container::iterator result_iterator;
  result_iterator r_;

 public:
  //! Constructor
  template <typename SourceIter, typename TargetIter, typename Options>
  DualTreeContext(const Expansion& K,
                  SourceIter sfirst, SourceIter slast,
                  TargetIter tfirst, TargetIter tlast,
                  Options& opts)
      : ExpansionContext<Expansion>(K),
        source_tree_(sfirst, slast, opts.ncrit),
        target_tree_(tfirst, tlast, opts.ncrit),
        M_(source_tree_.boxes()),
        L_((opts.evaluator == FMMOptions::TREECODE ? 0 : target_tree_.boxes())),
        sources_(sfirst, slast),
        targets_(tfirst, tlast) {
    s_ = sources_.begin();
    t_ = targets_.begin();
  }

  template <typename Executor>
  inline void execute(const std::vector<charge_type>& charges,
                      std::vector<result_type>& results,
                      Executor* exec) {
    c_ = charges.begin();
    r_ = results.begin();

    exec->execute(*this);
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

  // Accessors to make this Executor into a BoxContext
  inline multipole_type& multipole(const source_box_type& box) {
    return M_[box.index()];
  }
  inline const multipole_type& multipole(const source_box_type& box) const {
    return M_[box.index()];
  }
  inline local_type& local(const target_box_type& box) {
    return L_[box.index()];
  }
  inline const local_type& local(const target_box_type& box) const {
    return L_[box.index()];
  }

  typedef body_transform_iterator<source_body_iterator,
                                  source_iterator>      body_source_iterator;
  inline body_source_iterator source_begin(const source_box_type& b) const {
    return transform_iterator(b.body_begin(), s_);
  }
  inline body_source_iterator source_end(const source_box_type& b) const {
    return transform_iterator(b.body_end(), s_);
  }
  inline body_source_iterator source_begin() const {
    return transform_iterator(source_tree().body_begin(), s_);
  }
  inline body_source_iterator source_end() const {
    return transform_iterator(source_tree().body_end(), s_);
  }

  typedef body_transform_iterator<source_body_iterator,
                                  charge_iterator>      body_charge_iterator;
  inline body_charge_iterator charge_begin(const source_box_type& b) const {
    return transform_iterator(b.body_begin(), c_);
  }
  inline body_charge_iterator charge_end(const source_box_type& b) const {
    return transform_iterator(b.body_end(), c_);
  }
  inline body_charge_iterator charge_begin() const {
    return transform_iterator(source_tree().body_begin(), c_);
  }
  inline body_charge_iterator charge_end() const {
    return transform_iterator(source_tree().body_end(), c_);
  }

  typedef body_transform_iterator<target_body_iterator,
                                  target_iterator>      body_target_iterator;
  inline body_target_iterator target_begin(const target_box_type& b) const {
    return transform_iterator(b.body_begin(), t_);
  }
  inline body_target_iterator target_end(const target_box_type& b) const {
    return transform_iterator(b.body_end(), t_);
  }
  inline body_target_iterator target_begin() const {
    return transform_iterator(target_tree().body_begin(), t_);
  }
  inline body_target_iterator target_end() const {
    return transform_iterator(target_tree().body_end(), t_);
  }

  typedef body_transform_iterator<target_body_iterator,
                                  result_iterator>      body_result_iterator;
  inline body_result_iterator result_begin(const target_box_type& b) {
    return transform_iterator(b.body_begin(), r_);
  }
  inline body_result_iterator result_end(const target_box_type& b) {
    return transform_iterator(b.body_end(), r_);
  }
  inline body_result_iterator result_begin() {
    return transform_iterator(target_tree().body_begin(), r_);
  }
  inline body_result_iterator result_end() {
    return transform_iterator(target_tree().body_end(), r_);
  }
};



template <typename Expansion>
class SingleTreeContext<Expansion,
                        Octree<typename ExpansionTraits<Expansion>::point_type> >
    : public ExpansionContext<Expansion>
{
 public:
  typedef ExpansionTraits<Expansion> expansion_traits;
  FMMTL_IMPORT_EXPANSION_TRAITS(expansion_traits);

  typedef TreePairTraits<Octree<point_type>, Octree<point_type>> tree_traits;
  FMMTL_IMPORT_TREEPAIR_TRAITS(tree_traits);

 protected:
  // Transform a body to a value
  template <typename BodyIter, typename Indexable>
  struct Transformer {
    typedef typename Indexable::reference result_type;
    typedef typename BodyIter::value_type argument_type;
    Transformer(const Indexable& value) : value_(value) {}
    result_type operator()(const argument_type& body) const {
      return value_[body.number()]; // TODO: TEMP to avoid permutation
    }
   private:
    Indexable value_;
  };

  template <typename BodyIter, typename Indexable>
  using body_transform_iterator =
      boost::transform_iterator<Transformer<BodyIter,Indexable>,BodyIter>;

  template <typename BodyIter, typename Indexable>
  inline body_transform_iterator<BodyIter, Indexable>
  transform_iterator(BodyIter bi, Indexable v) const {
    return boost::make_transform_iterator(bi,
                                          Transformer<BodyIter,Indexable>(v));
  }

  //! The tree of sources and targets
  source_tree_type source_tree_;

  //! Multipole expansions corresponding to Box indices in Tree
  typedef std::vector<multipole_type> multipole_container;
  multipole_container M_;
  //! Local expansions corresponding to Box indices in Tree
  typedef std::vector<local_type> local_container;
  local_container L_;

  //! The sources (and targets) associated with bodies in the source_tree_
  typedef const std::vector<source_type> source_container;
  source_container sources_;

  //! Iterator to the start of the source vector
  typedef typename source_container::const_iterator source_iterator;
  source_iterator s_;
  typedef source_container target_container;
  typedef source_iterator  target_iterator;
  //! Iterator to the start of the charge vector
  typedef std::vector<charge_type> charge_container;
  typedef typename charge_container::const_iterator charge_iterator;
  charge_iterator c_;
  //! Iterator to the start of the result vector
  typedef std::vector<result_type> result_container;
  typedef typename result_container::iterator result_iterator;
  result_iterator r_;

 public:
  //! Constructor
  template <typename SourceIter, typename Options>
  SingleTreeContext(const Expansion& K,
                    SourceIter sfirst, SourceIter slast,
                    Options& opts)
      : ExpansionContext<Expansion>(K),
        source_tree_(sfirst, slast, opts.ncrit),
        M_(source_tree_.boxes()),
        L_((opts.evaluator == FMMOptions::TREECODE ? 0 : source_tree_.boxes())),
        sources_(sfirst, slast) {
    s_ = sources_.begin();
  }

  template <typename Executor>
  inline void execute(const std::vector<charge_type>& charges,
                      std::vector<result_type>& results,
                      Executor* exec) {
    c_ = charges.begin();
    r_ = results.begin();

    exec->execute(*this);
  }

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

  // Accessors to make this Executor into a BoxContext
  inline multipole_type& multipole(const source_box_type& box) {
    return M_[box.index()];
  }
  inline const multipole_type& multipole(const source_box_type& box) const {
    return M_[box.index()];
  }
  inline local_type& local(const target_box_type& box) {
    return L_[box.index()];
  }
  inline const local_type& local(const target_box_type& box) const {
    return L_[box.index()];
  }

  typedef body_transform_iterator<source_body_iterator,
                                  source_iterator>      body_source_iterator;
  inline body_source_iterator source_begin(const source_box_type& b) const {
    return transform_iterator(b.body_begin(), s_);
  }
  inline body_source_iterator source_end(const source_box_type& b) const {
    return transform_iterator(b.body_end(), s_);
  }
  inline body_source_iterator source_begin() const {
    return transform_iterator(source_tree().body_begin(), s_);
  }
  inline body_source_iterator source_end() const {
    return transform_iterator(source_tree().body_end(), s_);
  }

  typedef body_transform_iterator<source_body_iterator,
                                  charge_iterator>      body_charge_iterator;
  inline body_charge_iterator charge_begin(const source_box_type& b) const {
    return transform_iterator(b.body_begin(), c_);
  }
  inline body_charge_iterator charge_end(const source_box_type& b) const {
    return transform_iterator(b.body_end(), c_);
  }
  inline body_charge_iterator charge_begin() const {
    return transform_iterator(source_tree().body_begin(), c_);
  }
  inline body_charge_iterator charge_end() const {
    return transform_iterator(source_tree().body_end(), c_);
  }

  typedef body_transform_iterator<target_body_iterator,
                                  target_iterator>      body_target_iterator;
  inline body_target_iterator target_begin(const target_box_type& b) const {
    return transform_iterator(b.body_begin(), s_);
  }
  inline body_target_iterator target_end(const target_box_type& b) const {
    return transform_iterator(b.body_end(), s_);
  }
  inline body_target_iterator target_begin() const {
    return transform_iterator(target_tree().body_begin(), s_);
  }
  inline body_target_iterator target_end() const {
    return transform_iterator(target_tree().body_end(), s_);
  }

  typedef body_transform_iterator<target_body_iterator,
                                  result_iterator>      body_result_iterator;
  inline body_result_iterator result_begin(const target_box_type& b) {
    return transform_iterator(b.body_begin(), r_);
  }
  inline body_result_iterator result_end(const target_box_type& b) {
    return transform_iterator(b.body_end(), r_);
  }
  inline body_result_iterator result_begin() {
    return transform_iterator(target_tree().body_begin(), r_);
  }
  inline body_result_iterator result_end() {
    return transform_iterator(target_tree().body_end(), r_);
  }
};
