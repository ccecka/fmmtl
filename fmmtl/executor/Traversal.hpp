#pragma once
/** @file Traversal
 * @brief A simple queue-based tree traversal that passes box pairs to
 * near and far field evaluators.
 */

/**
 * concept box {
 *   int index() const;
 *   bool is_leaf() const;
 *   double volume() const;
 *   box_iterator child_begin() const;
 *   box_iterator child_end() const;
 * }
 *
 * concept NearEvaluator {
 *   // Handles a box pair with near field evaluation
 *   void operator()(source_box&, target_box&);
 * }
 * concept FarEvaluator {
 *   // @returns whether box pair will be handled as a far field or not
 *   bool operator()(source_box&, target_box&);
 * }
 */

#include <utility>
#include <queue>
#include <tuple>

template <typename SourceBox, typename TargetBox>
struct Traversal {
  typedef SourceBox source_box;
  typedef TargetBox target_box;

  typedef std::pair<source_box, target_box> box_pair;

  typedef std::queue<box_pair> box_queue;

  template <class NearEvaluator, class FarEvaluator>
  inline static void eval(source_box sbox, target_box tbox,
                          NearEvaluator& near_eval, FarEvaluator& far_eval) {
    // Queue based tree traversal
    box_queue pairQ;
    interact(sbox, tbox, far_eval, pairQ);

    while (!pairQ.empty()) {
      std::tie(sbox, tbox) = pairQ.front();
      pairQ.pop();

      char code = (sbox.is_leaf() << 1) | (tbox.is_leaf() << 0);
      switch (code) {
        case 0: {             // sbox and tbox are not leaves
          // Split the larger of the two into children and interact
          if (sbox.volume() > tbox.volume()) {
            case 1:           // tbox is a leaf, sbox is not a leaf
              split_source(sbox, tbox, far_eval, pairQ);
          } else {
            case 2:           // sbox is a leaf, tbox is not a leaf
              split_target(sbox, tbox, far_eval, pairQ);
          }
        } continue;

        case 3: {             // sbox and tbox are leaves
          near_eval(sbox, tbox);
        } continue;
      } // end switch
    } // end while
  }

 private:

  template <class FarEvaluator>
  inline static void split_source(source_box& sbox, target_box& tbox,
                                  FarEvaluator& far_eval, box_queue& pairQ) {
    auto c_end = sbox.child_end();
    for (auto cit = sbox.child_begin(); cit != c_end; ++cit) {
      source_box cbox = *cit;
      interact(cbox, tbox, far_eval, pairQ);
    }
  }

  template <class FarEvaluator>
  inline static void split_target(source_box& sbox, target_box& tbox,
                                  FarEvaluator& far_eval, box_queue& pairQ) {
    auto c_end = tbox.child_end();
    for (auto cit = tbox.child_begin(); cit != c_end; ++cit) {
      target_box cbox = *cit;
      interact(sbox, cbox, far_eval, pairQ);
    }
  }

  template <class FarEvaluator>
  inline static void interact(source_box& sbox, target_box& tbox,
                              FarEvaluator& far_eval, box_queue& pairQ) {
    if (!far_eval(sbox, tbox))
      pairQ.emplace(sbox, tbox);
  }
};


struct Traverse {
  template <class SBox, class TBox, class NearEvaluator, class FarEvaluator>
  static void eval(SBox s, TBox t, NearEvaluator& n, FarEvaluator& f) {
    Traversal<SBox, TBox>().eval(s, t, n, f);
  }
  template <class Box, class NearEvaluator, class FarEvaluator>
  static void eval(Box b, NearEvaluator& n, FarEvaluator& f) {
    Traverse::eval(b, b, n, f);
  }
};
