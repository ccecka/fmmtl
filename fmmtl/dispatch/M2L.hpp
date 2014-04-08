#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

/** Default behavior gives a warning -- using non-existent method */
template <bool has_l2p>
struct M2L_Helper {
  inline static void apply(...) {
    std::cerr << "WARNING: Expansion does not have a correct M2L!\n";
  }
  inline static void eval(...) {
    std::cerr << "WARNING: Expansion does not have a correct M2L!\n";
  }
};

/** Expansion has an M2L method to dispatch to */
template <>
struct M2L_Helper<true> {
  /** The Expansion provides an M2L accumulator */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& translation) {
    K.M2L(M, L, translation);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.multipole(sbox),
          c.local(tbox),
          tbox.center() - sbox.center());
  }
};

// Public M2L dispatcher
struct M2L {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           typename Expansion::local_type& L,
                           const typename Expansion::point_type& translation) {
    typedef M2L_Helper<ExpansionTraits<Expansion>::has_M2L> M2L_H;
    M2L_H::apply(K, M, L, translation);
  }

  /** Forward to M2L_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2L");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef M2L_Helper<expansion_traits::has_M2L> M2L_H;
    M2L_H::eval(c, sbox, tbox);
  }
};


/** A lazy M2L evaluator which saves a list of pairs of boxes
 * That are sent to the M2L dispatcher on demand.
 */
template <typename Context>
class M2L_Batch {
  //! Type of box
  typedef typename Context::source_box_type source_box_type;
  typedef typename Context::target_box_type target_box_type;

  std::vector<target_box_type> target_box_list;
  std::vector<std::vector<source_box_type>> source_boxes;

  //! Box list for M2L interactions    TODO: could further compress these...
  //typedef std::pair<source_box_type, target_box_type> box_pair;
  //std::vector<box_pair> m2l_list;

 public:

  /** Insert a source-target box interaction to the interaction list */
  void insert(const source_box_type& s, const target_box_type& t) {
    if (source_boxes.size() <= t.index()) {
      source_boxes.resize(t.index() + 1);
      target_box_list.push_back(t);
    } else if (source_boxes[t.index()].empty()) {
      target_box_list.push_back(t);
    }

    source_boxes[t.index()].push_back(s);

    //m2l_list.push_back(std::make_pair(s,t));
  }

  /** Compute all interations in the interaction list */
  void execute(Context& c) {
    FMMTL_LOG("M2L Batch");
#pragma omp parallel for
    for (auto ti = target_box_list.begin(); ti < target_box_list.end(); ++ti) {
      target_box_type& tb = *ti;
      auto s_end = source_boxes[tb.index()].end();
      for (auto si = source_boxes[tb.index()].begin(); si != s_end; ++si) {
        M2L::eval(c, *si, tb);
      }
    }

    /*
    auto b_end = m2l_list.end();
    //#pragma omp parallel for   TODO: Make thread safe!
    for (auto bi = m2l_list.begin(); bi < b_end; ++bi)
      M2L::eval(c, bi->first, bi->second);
    */
  }
};
