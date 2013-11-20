#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/KernelTraits.hpp"
#include <type_traits>

class M2L
{
  /** If no other M2L dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct M2L!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_M2L>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& source,
       typename Expansion::local_type& target,
       const typename Expansion::point_type& translation) {
    K.M2L(source, target, translation);
  }

 public:

  /** Unwrap data from Context and dispatch to the M2L evaluator
   */
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

    M2L::eval(c.expansion(),
              c.multipole(sbox),
              c.local(tbox),
              tbox.center() - sbox.center());
  }
};


#include "Evaluator.hpp"

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
