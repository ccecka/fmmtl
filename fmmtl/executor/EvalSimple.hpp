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

  P2P_Batch<Context> p2p_batch;
  M2L_Batch<Context> m2l_batch;

  std::function<bool(source_box, target_box)> mac_;

  // TEMP: experimental
  //P2P_Compressed<typename Context::kernel_type> p2p_comp;

  struct Dispatch {
    Context& c_;
    Dispatch(Context& c) : c_(c) {}

    void up_process(typename Context::source_box_type& box) {
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
    void down_process(typename Context::target_box_type& box) {
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

 public:

  template <class MAC>
  EvalSimple(Context& c, const MAC& mac)
      : mac_(mac) {
    // Perform the box interactions
    Traverse::eval(c.source_tree().root(), c.target_tree().root(), *this);
    // TEMP: experimental
    //p2p_comp = p2p_batch.compressed(c);
  }

  void execute(Context& c) {
    // Initialize all the multipoles and locals (not all may be needed)
    auto s_end = c.source_tree().box_end();
    for (auto bi = c.source_tree().box_begin(); bi != s_end; ++bi)
      INITM::eval(c, *bi);
    auto t_end = c.target_tree().box_end();
    for (auto bi = c.target_tree().box_begin(); bi != t_end; ++bi)
      INITL::eval(c, *bi);

    // Perform the upward pass (not all may be needed)
    Dispatch updown(c);
    UpwardPass::eval(c.source_tree(), updown);

    m2l_batch.execute(c);

    // Perform the box interactions (not all may be needed)
    DownwardPass::eval(c.target_tree(), updown);

    p2p_batch.execute(c);

    // TEMP: experimental
    //FMMTL_LOG("P2P Compressed");
    //p2p_comp.execute(c);
  }

  /*******************/
  /** Functions called by the Traverse algorithm */
  /*******************/

  void near_field(const source_box& s, const target_box& t) {
    p2p_batch.insert(s,t);
  }

  bool far_field(const source_box& s, const target_box& t) {
    if (mac_(s,t)) {
      m2l_batch.insert(s,t);
      return true;
    }
    return false;
  }
};


template <class Context, class Options>
EvaluatorBase<Context>* make_eval_simple(Context& c, Options& opts) {
  return new EvalSimple<Context>(c, opts.MAC());
}
