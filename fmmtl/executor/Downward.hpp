#pragma once

/**
 * concept Evaluator {
 *   void up_process(box_type& b)
 * }
 */
struct DownwardPass {
  template <class Tree, class Evaluator>
  inline static void eval(Tree& tree, Evaluator& e) {
    // For the highest level down to the lowest level
    for (unsigned l = 0; l < tree.levels(); ++l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        e.down_process(box);
      }
    }
  }
};
