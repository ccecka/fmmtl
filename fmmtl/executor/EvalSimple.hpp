#pragma once

#include "Evaluator.hpp"

#include "Upward.hpp"
#include "Traversal.hpp"
#include "Downward.hpp"

#include "Dispatchers.hpp"

template <class Context>
class EvalSimple
    : public EvaluatorBase<Context>
{
  typedef typename Context::source_box_type source_box;
  typedef typename Context::target_box_type target_box;

  Context& c_;  // TODO: Fix this design
  std::function<bool(const source_box&, const target_box&)> mac_;

 public:

  template <class MAC>
  EvalSimple(Context& c, const MAC& mac)
      : c_(c), mac_(mac) {
    // Do precomputation (precompute interaction lists, etc)
  }

  void execute(Context&) {
    // Initialize all the multipoles and locals (not all may be needed)
    auto s_end = c_.source_tree().box_end();
    for (auto bi = c_.source_tree().box_begin(); bi != s_end; ++bi)
      INITM::eval(c_, *bi);
    auto t_end = c_.target_tree().box_end();
    for (auto bi = c_.target_tree().box_begin(); bi != t_end; ++bi)
      INITL::eval(c_, *bi);

    // Perform the upward pass (not all may be needed)
    UpwardPass::eval(c_.source_tree(), *this);

    // Perform the box interactions
    Traverse::eval(c_.source_tree().root(), c_.target_tree().root(), *this);

    // Perform the box interactions (not all may be needed)
    DownwardPass::eval(c_.target_tree(), *this);
  }

  /*******************/
  /** Functions called by the UpwardPass algorithm */
  /*******************/

  void up_process(const source_box& box) {
    if (box.is_leaf()) {
      // If leaf, make P2M calls
      P2M::eval(c_, box);
    } else {
      // If not leaf, then for all the children M2M
      auto c_end = box.child_end();
      for (auto cit = box.child_begin(); cit != c_end; ++cit)
        M2M::eval(c_, *cit, box);
    }
  }

  /*******************/
  /** Functions called by the Traverse algorithm */
  /*******************/

  void near_field(const source_box& s, const target_box& t) {
    P2P::eval(c_, s, t, P2P::ONE_SIDED());
  }

  bool far_field(const source_box& s, const target_box& t) {
    if (mac_(s,t)) {
      M2L::eval(c_, s, t);
      return true;
    }
    return false;
  }

  /*******************/
  /** Functions called by the DownPass algorithm */
  /*******************/

  void down_process(const target_box& box) {
    if (box.is_leaf()) {
      // If leaf, make L2P calls
      L2P::eval(c_, box);
    } else {
      // If not leaf, then for all children L2L
      auto c_end = box.child_end();
      for (auto cit = box.child_begin(); cit != c_end; ++cit)
        L2L::eval(c_, box, *cit);
    }
  }
};


template <class Context, class Options>
EvaluatorBase<Context>* make_eval_simple(Context& c, Options& opts) {
  return new EvalSimple<Context>(c, opts.MAC());
}
