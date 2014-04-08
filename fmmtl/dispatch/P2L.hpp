#pragma once
/** @file P2L.hpp
 * @brief Dispatch methods for the P2L stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

/** Default behavior gives a warning -- using non-existent method */
template <bool has_p2l>
struct P2L_Helper {
  inline static void apply(...) {
    std::cerr << "WARNING: Expansion does not have a correct P2L!\n";
  }
  inline static void eval(...) {
    std::cerr << "WARNING: Expansion does not have a correct P2L!\n";
  }
};

/** Expansion has an P2L method (scalar/vector) to dispatch to */
template <>
struct P2L_Helper<true> {
  /** The Expansion provides a vector P2L accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_P2L>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    K.P2L(s_begin, s_end, c_begin, center, L);
  }

  /** The Expansion provides a scalar P2L accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2L &
                          !ExpansionTraits<Expansion>::has_vector_P2L>::type
  apply(const Expansion& K,
        const typename Expansion::source_type& source,
        const typename Expansion::charge_type& charge,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    K.P2L(source, charge, center, L);
  }

  /** The Expansion provides a scalar P2L accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2L &
                          !ExpansionTraits<Expansion>::has_vector_P2L>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::local_type& L) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      apply(K, *s_begin, *c_begin, center, L);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.source_begin(sbox), c.source_end(sbox), c.charge_begin(sbox),
          tbox.center(),
          c.local(tbox));
  }
};

/** Public P2L dispatcher */
class P2L {
  /** Forward to P2L_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::source_type& source,
                           const typename Expansion::charge_type& charge,
                           const typename Expansion::point_type& center,
                           typename Expansion::local_type& L) {
    typedef P2L_Helper<ExpansionTraits<Expansion>::has_P2L> P2L_H;
    P2L_H::apply(K, source, charge, center, L);
  }

  /** Forward to P2L_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "P2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("P2L");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef P2L_Helper<expansion_traits::has_P2L> P2L_H;
    P2L_H::eval(c, sbox, tbox);
  }
};
