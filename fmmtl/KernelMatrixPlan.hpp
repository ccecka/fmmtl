#pragma once

#include "KernelTraits.hpp"
#include "executor/make_executor.hpp"
#include "tree/NDTree.hpp"
#include "executor/Context.hpp"

/** Abstract PlanBase class */
template <class Expansion>
class PlanBase {
 public:
  typedef typename ExpansionTraits<Expansion>::expansion_type expansion_type;
	typedef typename ExpansionTraits<Expansion>::source_type    source_type;
	typedef typename ExpansionTraits<Expansion>::target_type    target_type;
	typedef typename ExpansionTraits<Expansion>::charge_type    charge_type;
	typedef typename ExpansionTraits<Expansion>::result_type    result_type;

  /** Virtual destructor */
  virtual ~PlanBase() {};

  /** Execute this plan */
  virtual void execute(const std::vector<charge_type>& charges,
                       std::vector<result_type>& results) = 0;

  /** Accessors */

  /** The potentially reordered targets for this plan */
  virtual std::vector<target_type> targets() const = 0;

  /** The potentially reordered sources for this plan */
  virtual std::vector<source_type> sources() const = 0;
};


template <typename Expansion, typename Context>
struct Plan
    : public PlanBase<Expansion> {
  typedef typename ExpansionTraits<Expansion>::expansion_type expansion_type;
	typedef typename ExpansionTraits<Expansion>::source_type    source_type;
	typedef typename ExpansionTraits<Expansion>::target_type    target_type;
	typedef typename ExpansionTraits<Expansion>::charge_type    charge_type;
	typedef typename ExpansionTraits<Expansion>::result_type    result_type;

  template <typename KernelMatrix, typename Options>
  Plan(const KernelMatrix& mat, Options& opts)
      : context(mat, opts),
        executor(make_evaluator(context, opts)) {
    if (opts.print_tree) {
      std::cout << "Source Tree:\n" << context.source_tree() << std::endl;
      std::cout << "Target Tree:\n" << context.target_tree() << std::endl;
    }
  }

  virtual ~Plan() {
    delete executor;
  }

  virtual void execute(const std::vector<charge_type>& charges,
                       std::vector<result_type>& results) {
    return context.execute(charges, results, executor);
  }

  virtual std::vector<target_type> targets() const {
    return std::vector<target_type>(context.target_begin(),
                                    context.target_end());
  }
  virtual std::vector<source_type> sources() const {
    return std::vector<source_type>(context.source_begin(),
                                    context.source_end());
  }

 private:
  Context context;
  EvaluatorBase<Context>* executor;
};





template <class KernelMatrix, class Options>
PlanBase<typename KernelMatrix::expansion_type>*
make_kernel_matrix_plan(const KernelMatrix& mat, const Options& opts) {
  // Statically compute the plan type, potentially from Option types

  typedef typename KernelMatrix::expansion_type expansion_type;

  // The source and target tree types
  typedef typename expansion_type::point_type point_type;
  typedef NDTree<fmmtl::dimension<point_type>::value> source_tree_type;
  typedef NDTree<fmmtl::dimension<point_type>::value> target_tree_type;

  typedef typename expansion_type::source_type source_type;
  typedef typename expansion_type::target_type target_type;

  // Check if source and target sets are the same
  if (std::is_same<source_type, target_type>::value) {
    if (mat.sources() == mat.targets()) {    // TODO: fix O(N) with aliased test
#if defined(FMMTL_DEBUG)
      std::cout << "Using single tree context." << std::endl;
#endif
      typedef SingleTreeContext<source_tree_type> tree_context_type;
      typedef DataContext<KernelMatrix, tree_context_type> context_type;

      typedef Plan<expansion_type, context_type> plan_type;
      return new plan_type(mat, opts);
    }
  }

  // Source and target sets are unique
#if defined(FMMTL_DEBUG)
      std::cout << "Using dual tree context." << std::endl;
#endif
  typedef DualTreeContext<source_tree_type, target_tree_type> tree_context_type;
  typedef DataContext<KernelMatrix, tree_context_type> context_type;

  typedef Plan<expansion_type, context_type> plan_type;
  return new plan_type(mat, opts);
}
