#pragma once

#include "config.hpp"

//! Global logging
#include "Logger.hpp"
static Logger fmm_global_log;

// FMM includes
#include "KernelTraits.hpp"
#include "executor/make_executor.hpp"

#include <vector>

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

  /** The number of rows (number of targets) of the abstract matrix */
  virtual unsigned rows() const = 0;
  /** The (potentially reordered) targets generating the abstract matrix */
  virtual std::vector<target_type> targets() const = 0;

  /** The number of cols (number of sources) of the abstract matrix */
  virtual unsigned cols() const = 0;
  /** The (potentially reordered) sources generating the abstract matrix */
  virtual std::vector<source_type> sources() const = 0;

  /** Access to the kernel/expansion that generates this matrix */
  virtual expansion_type& expansion() = 0;
  virtual const expansion_type& expansion() const = 0;
};


template <typename Expansion, typename Context>
struct Plan
    : public PlanBase<Expansion> {
  typedef typename ExpansionTraits<Expansion>::expansion_type expansion_type;
	typedef typename ExpansionTraits<Expansion>::source_type    source_type;
	typedef typename ExpansionTraits<Expansion>::target_type    target_type;
	typedef typename ExpansionTraits<Expansion>::charge_type    charge_type;
	typedef typename ExpansionTraits<Expansion>::result_type    result_type;

  template <typename SourceIter, typename Options>
  Plan(const Expansion& K,
       SourceIter sbegin, SourceIter send,
       Options& opts)
      : context(K, sbegin, send, opts),
        executor(make_evaluator(context, opts)) {
    if (opts.print_tree) {
      std::cout << "Source/Target Tree:\n" << context.source_tree() << std::endl;
    }
  }

  template <typename SourceIter, typename TargetIter, typename Options>
  Plan(const Expansion& K,
       SourceIter sbegin, SourceIter send,
       TargetIter tbegin, TargetIter tend,
       Options& opts)
      : context(K, sbegin, send, tbegin, tend, opts),
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
    FMMTL_ASSERT(charges.size() == cols());
    FMMTL_ASSERT(results.size() == rows());
    return context.execute(charges, results, executor);
  }

  virtual unsigned rows() const {
    return context.target_tree().bodies();
  }
  virtual std::vector<target_type> targets() const {
    return std::vector<target_type>(context.target_begin(),
                                    context.target_end());
  }

  virtual unsigned cols() const {
    return context.source_tree().bodies();
  }
  virtual std::vector<source_type> sources() const {
    return std::vector<source_type>(context.source_begin(),
                                    context.source_end());
  }

  virtual expansion_type& expansion() {
    return context.expansion();
  }
  virtual const expansion_type& expansion() const {
    return context.expansion();
  }

 private:
  Context context;
  EvaluatorBase<Context>* executor;
};

template <class Expansion>
class fmm_matrix {  // Inherit ExpasionTraits?
 public:
  typedef typename ExpansionTraits<Expansion>::expansion_type expansion_type;
  typedef typename ExpansionTraits<Expansion>::kernel_type    kernel_type;

	typedef typename ExpansionTraits<Expansion>::source_type    source_type;
	typedef typename ExpansionTraits<Expansion>::target_type    target_type;
	typedef typename ExpansionTraits<Expansion>::charge_type    charge_type;
	typedef typename ExpansionTraits<Expansion>::result_type    result_type;

  fmm_matrix(PlanBase<Expansion>* _plan)
      : plan(_plan) {
  }
  ~fmm_matrix() {
    delete plan;
    FMMTL_PRINT_LOG(std::cout);
    fmmtl_logger.clear();
  };

  inline unsigned rows() const {
    return plan->rows();
  }
  inline unsigned cols() const {
    return plan->cols();
  }
  inline unsigned size() const {
    return rows() * cols();
  }

  inline void execute(const std::vector<charge_type>& charges,
                      std::vector<result_type>& results) {
    return plan->execute(charges, results);
  }
  inline std::vector<result_type> execute(const std::vector<charge_type>& charges) {
    std::vector<result_type> results(rows());
    this->execute(charges, results);
    return results;
  }
  // Syntactic sugar, return a proxy that can be
  // manipulated or passed to fmm::direct, fmm::treecode, etc?
  inline std::vector<result_type> operator*(const std::vector<charge_type>& c) {
    return this->execute(c);
  }

  inline expansion_type& expansion() {
    return plan->expansion();
  }
  const expansion_type& expansion() const {
    return plan->expansion();
  }
  inline std::vector<target_type> targets() const {
    return plan->targets();
  }
  inline std::vector<source_type> sources() const {
    return plan->sources();
  }

 private:
  PlanBase<Expansion>* plan;
};



#include "tree/NDTree.hpp"
#include "executor/Context.hpp"

template <class Expansion,
          class Options = FMMOptions>
fmm_matrix<Expansion>
make_fmm_matrix(const Expansion& E,
                const std::vector<typename ExpansionTraits<Expansion>::source_type>& sources,
                Options opts = FMMOptions()) {
  // Statically compute the context type, potentially from Options
  typedef ExpansionContext<Expansion> expansion_context_type;

  typedef NDTree<ExpansionTraits<Expansion>::dimension> tree_type;
  typedef SingleTreeContext<Expansion, tree_type> tree_context_type;

  typedef DataContext<expansion_context_type, tree_context_type> context_type;

  // Construct the plan
  typedef Plan<Expansion, context_type> plan_type;
  plan_type* plan = new plan_type(E,
                                  sources.begin(), sources.end(),
                                  opts);
  return fmm_matrix<Expansion>(plan);
}

template <class Expansion,
          class Options = FMMOptions>
fmm_matrix<Expansion>
make_fmm_matrix(const Expansion& E,
                const std::vector<typename ExpansionTraits<Expansion>::source_type>& sources,
                const std::vector<typename ExpansionTraits<Expansion>::target_type>& targets,
                Options opts = FMMOptions()) {
  // Statically compute the context type, potentially from Option types
  typedef ExpansionContext<Expansion> expansion_context_type;

  typedef NDTree<ExpansionTraits<Expansion>::dimension> tree_type;
  typedef DualTreeContext<Expansion, tree_type> tree_context_type;

  typedef DataContext<expansion_context_type, tree_context_type> context_type;

  // Construct the plan
  typedef Plan<Expansion, context_type> plan_type;
  plan_type* plan = new plan_type(E,
                                  sources.begin(), sources.end(),
                                  targets.begin(), targets.end(),
                                  opts);
  return fmm_matrix<Expansion>(plan);
}


/*
struct FMM_NEAR_FIELD_ONLY {};
struct FMM_FAR_FIELD_ONLY {};

struct FMM_DYNAMIC_TRAVERSAL {};
struct FMM_INTERACTION_LIST {};
etc etc


template <class Expansion, typename... CompileOptions>
struct FMM_Builder {
  // Compute the FMM_Plan type from the CompileOptions

  typedef typename ExpansionTraits<Expansion>::expansion_type expansion_type;
  typedef typename ExpansionTraits<Expansion>::kernel_type    kernel_type;

	typedef typename ExpansionTraits<Expansion>::source_type    source_type;
	typedef typename ExpansionTraits<Expansion>::target_type    target_type;
	typedef typename ExpansionTraits<Expansion>::charge_type    charge_type;
	typedef typename ExpansionTraits<Expansion>::result_type    result_type;

  typedef typename ExpansionTraits<Expansion>::point_type     point_type;

  typedef Octree<point_type>                                  tree_type;



  template <class RunOptions>
  FMM_Plan<Expansion,
           Context,Executor>
  make(const Expansion& K,
       const std::vector<source_type>& s,
       const std::vector<target_type>& t,
       const RunOptions& opts) {
    // Construct the Tree, Context, Evaluators, etc from the Runtime options

    // Pass the runtime options on to the plan_type
  }

  template <class RunOptions>
  SOMETHING make(const Expansion& K,
                 const std::vector<typename Expansion::source_type>& s,
                 const RunOptions& opts) {
    // Construct the Tree, Context, Evaluators, etc from the Runtime options

    // Pass the runtime options on to the plan_type
  }
};
*/
