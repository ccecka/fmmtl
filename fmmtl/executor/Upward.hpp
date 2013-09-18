#pragma once


/** Helper for performing a full upward pass
 */
/**
 * concept Evaluator {
 *   void up_process(box_type& b)
 * }
 */
struct UpwardPass {
  template <typename Tree, class Evaluator>
	inline static void eval(Tree& tree, Evaluator e) {
		// For the lowest level up to the highest level
		for (int l = tree.levels()-1; l >= 0; --l) {
			// For all boxes at this level
			auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        e.up_process(box);
			}
		}
	}
};


#include "INITM.hpp"
#include "P2M.hpp"
#include "M2M.hpp"

/** Helper for computing the multipole for a box and all sub-boxes
 */
struct ComputeM {
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox) {
    INITM::eval(c, sbox);

    if (sbox.is_leaf()) {
      // Compute the multipole from the box's sources
      P2M::eval(c, sbox);
    } else {
      auto c_end = sbox.child_end();
      for (auto cit = sbox.child_begin(); cit != c_end; ++cit) {
        typename Context::source_box_type child = *cit;

        // Recursively initialize the multipole
        ComputeM::eval(c, child);
        // Accumulate the child's multipole into sbox
        M2M::eval(c, child, sbox);
      }
    }
  }
};
