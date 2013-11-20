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

  P2P_Batch<Context> p2p_;
  M2L_Batch<Context> m2l_;

  struct UpDispatch {
    Context& c_;
    UpDispatch(Context& c) : c_(c) {}

    inline void operator()(const source_box& box) {
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
  };
  struct DownDispatch {
    Context& c_;
    DownDispatch(Context& c) : c_(c) {}

    inline void operator()(const target_box& box) {
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

  EvalSimple(Context& c) {
    // Construct functors for dispatched near and far operators
    auto far_batcher = [&c,this](const source_box& s, const target_box& t) {
      if (MAC::eval(c,s,t)) {
        m2l_.insert(s,t);
        return true;
      }
      return false;
    };
    auto near_batcher = [this](const source_box& s, const target_box& t) {
      p2p_.insert(s,t);
    };
    // Determine the box interactions
    Traverse::eval(c.source_tree().root(), c.target_tree().root(),
                   near_batcher, far_batcher);
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
    UpDispatch up(c);
    UpwardPass::eval(c.source_tree(), up);

    // Perform the source-target box interactions
    m2l_.execute(c);
    p2p_.execute(c);

    // Perform the downward pass (not all may be needed)
    DownDispatch down(c);
    DownwardPass::eval(c.target_tree(), down);
  }
};


template <class Context, class Options>
EvaluatorBase<Context>* make_eval_simple(Context& c, Options&) {
  return new EvalSimple<Context>(c);
}
