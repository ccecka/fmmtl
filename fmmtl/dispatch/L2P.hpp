#pragma once
/** @file L2P.hpp
 * @brief Dispatch methods for the L2P stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

/** Default behavior gives a warning -- using non-existent method */
template <bool has_l2p>
struct L2P_Helper {
  inline static void apply(...) {
    std::cerr << "WARNING: Expansion does not have a correct L2P!\n";
  }
  inline static void eval(...) {
    std::cerr << "WARNING: Expansion does not have a correct L2P!\n";
  }
};

/** Expansion has an L2P method (scalar/vector) to dispatch to */
template <>
struct L2P_Helper<true> {
  /** The Expansion provides a vector L2P accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_L2P>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.L2P(L, center, t_begin, t_end, r_begin);
  }

  /** The Expansion provides a scalar L2P accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_L2P &
                          !ExpansionTraits<Expansion>::has_vector_L2P>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        const typename Expansion::target_type& target,
              typename Expansion::result_type& result) {
    K.L2P(L, center, target, result);
  }

  /** The Expansion provides a scalar L2P accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_L2P &
                          !ExpansionTraits<Expansion>::has_vector_L2P>::type
  apply(const Expansion& K,
        const typename Expansion::local_type& L,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      apply(K, L, center, *t_begin, *r_begin);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.local(tbox),
          tbox.center(),
          c.target_begin(tbox), c.target_end(tbox),
          c.result_begin(tbox));
  }
};

// Public L2P dispatcher
struct L2P {
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::local_type& L,
                           const typename Expansion::point_type& center,
                           const typename Expansion::target_type& target,
                           typename Expansion::result_type& result) {
    typedef L2P_Helper<ExpansionTraits<Expansion>::has_L2P> L2P_H;
    L2P_H::apply(K, L, center, target, result);
  }

  /** Forward to L2P_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "L2P:"
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("L2P");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef L2P_Helper<expansion_traits::has_L2P> L2P_H;
    L2P_H::eval(c, tbox);
  }
};
