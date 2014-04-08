#pragma once
/** @file P2M.hpp
 * @brief Dispatch methods for the P2M stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

/** Default behavior gives a warning -- using non-existent method */
template <bool has_p2m>
struct P2M_Helper {
  inline static void apply(...) {
    std::cerr << "WARNING: Expansion does not have a correct P2M!\n";
  }
  inline static void eval(...) {
    std::cerr << "WARNING: Expansion does not have a correct P2M!\n";
  }
};

/** Expansion has an P2M method (scalar/vector) to dispatch to */
template <>
struct P2M_Helper<true> {
  /** The Expansion provides a vector P2M accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_P2M>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    K.P2M(s_begin, s_end, c_begin, center, M);
  }

  /** The Expansion provides a scalar P2M accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2M &
                          !ExpansionTraits<Expansion>::has_vector_P2M>::type
  apply(const Expansion& K,
        const typename Expansion::source_type& source,
        const typename Expansion::charge_type& charge,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    K.P2M(source, charge, center, M);
  }

  /** The Expansion provides a scalar P2M accumulator. */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2M &
                          !ExpansionTraits<Expansion>::has_vector_P2M>::type
  apply(const Expansion& K,
        SourceIter s_begin, SourceIter s_end, ChargeIter c_begin,
        const typename Expansion::point_type& center,
        typename Expansion::multipole_type& M) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      apply(K, *s_begin, *c_begin, center, M);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox) {
    apply(c.expansion(),
          c.source_begin(sbox), c.source_end(sbox), c.charge_begin(sbox),
          sbox.center(),
          c.multipole(sbox));
  }
};

/** Public P2M dispatcher */
struct P2M {
  /** Forward to P2M_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::source_type& source,
                           const typename Expansion::charge_type& charge,
                           const typename Expansion::point_type& center,
                           typename Expansion::multipole_type& M) {
    typedef P2M_Helper<ExpansionTraits<Expansion>::has_P2M> P2M_H;
    P2M_H::apply(K, source, charge, center, M);
  }

  /** Forward to P2M_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "P2M:"
              << "\n  " << sbox << std::endl;
#endif
    FMMTL_LOG("P2M");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef P2M_Helper<expansion_traits::has_P2M> P2M_H;
    P2M_H::eval(c, sbox);
  }
};
