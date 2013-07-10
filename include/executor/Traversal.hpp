#pragma once
/** @file Traversal
 * @brief A simple queue-based tree traversal that passes box pairs to an
 * @a Evaluators near_field and far_field operators.
 */

/**
 * concept box {
 *   int index() const;
 *   bool is_leaf() const;
 *   double size_length() const;  // Maybe level()?
 *   box_iterator child_begin() const;
 *   box_iterator child_end() const;
 * }
 *
 * concept Evaluator {
 *   void near_field(source_box&, target_box&);
 *   bool far_field(source_box&, target_box&);
 * }
 */

#include <queue>

template <typename SourceBox, typename TargetBox = SourceBox>
struct Traversal {
  typedef SourceBox source_box;
  typedef TargetBox target_box;

  typedef typename source_box::box_iterator source_box_iterator;
  typedef typename target_box::box_iterator target_box_iterator;

  typedef std::pair<source_box, target_box> box_pair;

  typedef std::queue<box_pair> box_queue;

  template <typename Evaluator>
  void eval(source_box& sbox, Evaluator& eval) {
    Traversal<source_box,target_box>::eval(sbox, sbox, eval);
  }

  template <typename Evaluator>
  void eval(source_box& sbox, target_box& tbox, Evaluator& eval) {
    // Queue based tree traversal
    box_queue pairQ;
    interact(sbox, tbox, eval, pairQ);

    while (!pairQ.empty()) {
      source_box sbox = pairQ.front().first;
      target_box tbox = pairQ.front().second;
      pairQ.pop();

      char code = (sbox.is_leaf() << 1) | (tbox.is_leaf() << 0);
      switch (code) {
        case 3: {    // sbox and tbox are leaves
          eval.near_field(sbox, tbox);
        } break;

        case 2: {    // sbox is a leaf, tbox is not leaf
          // Split tbox into children and interact
          split_target(sbox, tbox, eval, pairQ);
        } break;

        case 1: {    // sbox is not leaf, tbox is leaf
          // Split sbox into children and interact
          split_source(sbox, tbox, eval, pairQ);
        } break;

        case 0: {    // sbox and tbox are not leaves
          // Split the larger of the two into children and interact
          if (sbox.side_length() > tbox.side_length())
            split_source(sbox, tbox, eval, pairQ);
          else
            split_target(sbox, tbox, eval, pairQ);
        } break;
      } // end switch
    } // end while
  }

 private:

  template <typename Evaluator>
  inline static void split_source(source_box& sbox, target_box& tbox,
                                  Evaluator& eval, box_queue& pairQ) {
    source_box_iterator c_end = sbox.child_end();
    for (source_box_iterator cit = sbox.child_begin(); cit != c_end; ++cit) {
      source_box cbox = *cit;
      interact(cbox, tbox, eval, pairQ);
    }
  }

  template <typename Evaluator>
  inline static void split_target(source_box& sbox, target_box& tbox,
                                  Evaluator& eval, box_queue& pairQ) {
    target_box_iterator c_end = tbox.child_end();
    for (target_box_iterator cit = tbox.child_begin(); cit != c_end; ++cit) {
      target_box cbox = *cit;
      interact(sbox, cbox, eval, pairQ);
    }
  }

  template <typename Evaluator>
  inline static void interact(source_box& sbox, target_box& tbox,
                              Evaluator& eval, box_queue& pairQ) {
    if (!eval.far_field(sbox, tbox)) {
      pairQ.push(box_pair(sbox, tbox));
    }
  }
};


struct Traverse {
  template <typename Box, typename Evaluator>
  static void eval(Box b, Evaluator& e) {
    Traversal<Box>().eval(b, e);
  }
  template <typename SourceBox, typename TargetBox, typename Evaluator>
  static void eval(SourceBox s, TargetBox t, Evaluator& e) {
    Traversal<SourceBox, TargetBox>().eval(s, t, e);
  }
};
